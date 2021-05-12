namespace BurgersMlogN {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arithmetic;


    /// # Summary
    /// The first register has size M*log(N),
    /// corresponding to the M binary encodings of the N possible velocities at each position.
    /// The second register has size log(2M+2),
    /// corresponding to the 2M neighbors of a node
    /// + an extra state representing the node as a neighbor to itself (for self transitions). By convention it is state |2M+1>. 
    /// + an extra state representing the node (all other states represent directed edges). By convention it is state |0>.
    /// This last state flags what is also called the flat subspace.
    newtype WalkSpace = (node : Node,
    neighbor : LittleEndian,
    nPositions : Int,
    nVelocities : Int); 

    newtype Node = LittleEndian[]; // Array of sength N. Node[i] is a LittleEndian register encoding M velocities.

    /// # Summary
    /// returns One if the measurement projects to the WalkFlatSubspace.
    /// returns Zero if the measurement projects to a subspace orthogonal to the WalkFlatSubspace.
    operation ProjectToWalkFlatSubspace (walkState : WalkSpace) : Result {
        let neighbor = walkState::neighbor;
        use auxQubitFlagFlatSubspace = Qubit();      // use an auxiliary qubit that will be flipped iff neighbor is in the flat subspace.
        ApplyControlledOnInt (0, X, neighbor!, auxQubitFlagFlatSubspace);    // neighbor is in the flat subspace iff it is in state 0 in Little Endian representation.
        return MResetZ(auxQubitFlagFlatSubspace);
    }

    /// # Summary
    /// reflects about the walk flat subspace.
    operation ReflectAboutWalkFlatSubspace (walkState : WalkSpace) : Unit is Ctl {
        let neighbor = walkState::neighbor;
        ReflectAboutInteger (0 , neighbor);
    }

    operation ShiftBurgers (walkState : WalkSpace) : Unit is Ctl + Adj {
        // Note that ShiftBurgers is self adjoint. TODO: Give this information to the compiler.
        let (node, neighbor, nPositions, nVelocities) = walkState!;
        for i in 1 .. nPositions { // Index 0 corresponds to the flat subspace. 
            ApplyControlledOnInt (i, IncrementByInteger (1, _), neighbor!, node![i]);
            ApplyControlledOnInt (i, IncrementByInteger (-1, _), neighbor!, node![i+1]);
            X (neighbor![Length(neighbor!)-1]); // exchanges the +1 -1 transition and the -1 +1 transition.
        }
        for i in nPositions+1 .. 2*nPositions { // Index 2*nPositions+1 corresponds to a selftransition. 
            ApplyControlledOnInt (i, IncrementByInteger (-1, _), neighbor!, node![i]);
            ApplyControlledOnInt (i, IncrementByInteger (1, _), neighbor!, node![i+1]);
            X (neighbor![Length(neighbor!)-1]); // exchanges the +1 -1 transition and the -1 +1 transition.
        }
    }

    operation CoinTossBurgers (walkState : WalkSpace) : Unit is Ctl + Adj {
        let (node, neighbor, nPositions, nVelocities) = walkState!;
        let eta = 5; // eta is the number of bits used to encode transition probabilities.
        let nNeighbors = 2*nPositions+1;
        mutable transitionProbabilitiesArray = new Qubit[][nNeighbors];

        for neighborIndex in 1 .. nNeighbors {
            use auxQubitsTransitionProbability = Qubit[eta];
            set transitionProbabilitiesArray w/= neighborIndex <- auxQubitsTransitionProbability;
            // TODO: make this operation adjointable.
            ApplyTransitionProbabilitiesBurgers (node, neighborIndex, auxQubitsTransitionProbability, nPositions);
            // TODO: check that it works well with Q# to apply the transition probabilities after assigning a qubit register to an array.
        }

        let nLevels = 0; // TODO: ceil(log_2 (nNeighbors);)
        mutable transitionProbabilitiesTree = new Qubit[][][nLevels];
        for level in nLevels-1.. -1 .. 0 { // decreasing order, last level is zero.
            let nRegisters = 2^level;
            mutable transitionProbabilitiesLevel = new Qubit[][nRegisters];
            for i in 0 .. nRegisters {
                use register = Qubit[eta];
                AddI (transitionProbabilitiesTree[2*i][level+1], register);
                AddI (transitionProbabilitiesTree[2*i+1][level+1], register);
                set transitionProbabilitiesLevel w/= i <- register;
            }
            set transitionProbabilitiesTree w/= level <- transitionProbabilitiesLevel;
        }

        // Compute the transitions in amplitude following https://arxiv.org/pdf/0903.3465.pdf
        // we use the register neighbor as the log(d) qubits of the paper.
        // theta_0 = arccos(\sqrt(transitionProbabilitiesTree[0][1]));
        // R (theta_0)(neighbor[0]);
        for level in 0 .. nLevels-1 {
            for i in 0 .. 2^level {
                use flagRegister = Qubit[2];
                let cRegister = transitionProbabilitiesTree[2*i][level+1];
                let bRegister = transitionProbabilitiesTree[i][level];
                use thetaRegister = Qubit[eta]; // we may want something else than eta for the number of bits of precision of theta.
                DetermineAngleCircuit (flagRegister, cRegister, bRegister, thetaRegister);
                R (thetaRegister) (neighbor![level]);
                // TODO : is it possible to control a rotation by a qubit register ?
                // -> find a reversible circuit for this.
            }
        }

    }

    operation DetermineAngleCircuit (
        flagRegister : Qubit[],
        cRegister : Qubit[],
        bRegister : Qubit[],
        thetaRegister : Qubit[]
    )
    : Unit is Ctl + Adj {
        SpecialCaseCircuit (thetaRegister, flagRegister[0], cRegister, flagRegister[1], bRegister);
        // Controlled ThetaCircuit (flagRegister, (cRegister, bRegister, thetaRegister));
        // wrong : we want to implement the control defined below :
        // 00 the circuit θ computes normally ,
        // 01, 11 the circuit θ does nothing (keeps angle = 0, as b = c) ,
        // 10 the circuit θ outputs θ = π/2, as c = 0.
    }

    operation SpecialCaseCircuit (
        thetaRegister : Qubit[],
        flagRegister0 : Qubit,
        cRegister : Qubit[],
        flagRegister1 : Qubit,
        bRegister : Qubit[]
    )
    : Unit is Ctl + Adj {
        Equality (thetaRegister, cRegister, flagRegister0);
        Equality (cRegister, bRegister, flagRegister1);
    }

    operation Equality (
        firstRegister : Qubit[],
        secondRegister : Qubit[],
        outputQubit : Qubit
    )
    : Unit is Ctl + Adj {
        // Assumes that firstRegister and secondRegister have the same length.
        let n = Length(firstRegister);
        use auxQubits = Qubit[n];
        within {
            for i in 0 .. n-1 {
                CNOT (firstRegister[i], auxQubits[i]);
                CNOT (secondRegister[i], auxQubits[i]);
                X (auxQubits[i]);
            }
        } apply {
            Controlled X (auxQubits, outputQubit);
        }
    }

    operation ThetaCircuit (
        cRegister : Qubit[],
        bRegister : Qubit[],
        thetaRegister : Qubit[]
    )
    : Unit is Ctl + Adj {
        // TODO : reversible cicruit for arccos(\sqrt{cRegister/bRegister})
    }



    operation ApplyTransitionProbabilitiesBurgers (
        node : Node,
        neighborIndex : Int,
        auxQubitsTransitionProbabilities : Qubit[],
        nPositions : Int
    )
    : Unit is Ctl + Adj {
        if 1 <= neighborIndex and neighborIndex <= nPositions {
            // +1 -1 TransitionProbability
        }
        elif nPositions+1 <= neighborIndex and neighborIndex <= 2*nPositions {
            // -1 +1 TransitionPobability
        }
        elif nPositions==2*nPositions+1 {
            // SelfTransitionProbability
        }
        else {
            Message("Unexpected neighborIndex");
        }
    }
    


    operation WalkOperatorBurgers (walkState : WalkSpace) : Unit is Ctl {
        within {
            CoinTossBurgers(walkState);
        } apply {
            ShiftBurgers(walkState);
        }
        ReflectAboutWalkFlatSubspace(walkState);
    }


}