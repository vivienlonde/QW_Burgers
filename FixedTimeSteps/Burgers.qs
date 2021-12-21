namespace Burgers {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;

    open AmplitudeTransduction;
    open ArithmeticOperations;


    /// # Summary
    /// nPositions is M.
    /// nVelocities is N.
    /// The first register "node" has size M*log(N),
    /// corresponding to the M binary encodings of the N possible velocities at each position.
    /// The second register "neighbor" has size log(2M+2),
    /// corresponding to the 2M neighbors of a node
    /// (i.e the node with velocity +1 and the node with velocity -1 for each position i \in {1, .., M}),
    /// + an extra state representing the node as a neighbor to itself (for self transitions). By convention it is the state |2M+1>. 
    /// + an extra state representing the node (all other states represent directed edges). By convention it is hte state |0>.
    /// This last state flags what is also called the flat subspace.
    newtype WalkSpace = (
        node : Node,
        neighbor : Neighbor,
        nPositions : Int,
        nVelocities : Int
    ); 

    ////////////////// Node:
    // Array of length M.
    // For a position lambda, Node[Lambda] is a LittleEndian register on log(N) qubits encoding a velocity.
    newtype Node = LittleEndian[]; // change to SignedLittleEndian ?

    // deprecated:
    // A value i \in {1, .., M} for the neighbor register
    // corresponds to the node with velocity +1 at position i and -1 at position i+1.
    // A value i \in {M+1, .., 2M} for the neighbor register
    // corresponds to the node with velocity -1 at position i-M and +1 at position i+1-M.

    ///////////////// Neighbor:
    // the first element (of type LittleEndian) indexes a position lambda ($0 \leq \lambda \leq M-1$).
    // the second element (of type Qubit) discriminates between the two following neghbors:
    //      |0> corresponds to the LambdaPlus neighbor:
    //          (Nlambda, NlambdaPlusOne) -> (Nlambda + 1, NlambdaPlusOne - 1)
    //      |1> corresponds to the LambdaMinus neighbor:
    //          (Nlambda, NlambdaPlusOne) -> (Nlambda - 1, NlambdaPlusOne + 1)
    newtype Neighbor = (LittleEndian, Qubit);

    /// # Summary
    /// returns One if the measurement projects to the WalkFlatSubspace.
    /// returns Zero if the measurement projects to a subspace orthogonal to the WalkFlatSubspace.
    operation ProjectToWalkFlatSubspace (
        walkState : WalkSpace
    ) : Result {
        let neighbor = walkState::neighbor;
        use auxQubitFlagFlatSubspace = Qubit();      // use an auxiliary qubit that will be flipped iff neighbor is in the flat subspace.
        ApplyControlledOnInt (0, X, neighbor!, auxQubitFlagFlatSubspace);    // a walkstate is in the flat subspace iff its neighbor register is in state 0 in Little Endian representation.
        return MResetZ(auxQubitFlagFlatSubspace);
    }

    /// # Summary
    /// reflects about the walk flat subspace.
    operation ReflectAboutWalkFlatSubspace (
        walkState : WalkSpace
    ) : Unit is Ctl + Adj {
        body (...) {
            let neighbor = walkState::neighbor;
            ReflectAboutInteger (0 , neighbor);
        }
        adjoint self;
    }

    operation WalkOperatorBurgers (
        walkState : WalkSpace
    ) : Unit is Ctl + Adj {
        within {
            CoinTossBurgers(walkState);
        } apply {
            ShiftBurgers(walkState);
        }
        ReflectAboutWalkFlatSubspace(walkState);
    }
    // TODO: write a function "ConstructWalkOperator" that takes as input a CoinToss operator and a Shift operator
    // and outputs a Walk operator.
    // Then specialize it to define WalkOperatorBurgers.

    operation ShiftBurgers (
        walkState : WalkSpace
    ) : Unit is Ctl + Adj {
        body (...) {
            let (node, neighbor, nPositions, nVelocities) = walkState!;

            let PlusOne = IncrementByInteger (1, _);
            let MinusOne = IncrementByInteger (-1, _);

            for i in 1 .. nPositions { // Index 0 corresponds to the flat subspace.
                // modify the node register:
                ApplyControlledOnInt (i, PlusOne, neighbor!, node![i]);
                ApplyControlledOnInt (i, MinusOne, neighbor!, node![i+1]);
                // modify the neighbor register:
                X (neighbor![Length(neighbor!)-1]); // exchanges the +1 -1 transition and the -1 +1 transition.
            }
            for i in nPositions+1 .. 2*nPositions { // Index 2*nPositions+1 corresponds to a selftransition.
                // modify the node register:
                ApplyControlledOnInt (i, MinusOne, neighbor!, node![i-nPositions]);
                ApplyControlledOnInt (i, PlusOne, neighbor!, node![i+1-nPositions]);
                // modify the neighbor register:
                X (neighbor![Length(neighbor!)-1]); // exchanges the +1 -1 transition and the -1 +1 transition.
            }
            // The self-transition edge is left unchanged.
        }
        adjoint self;
    }

    operation CoinTossBurgers (
        walkState : WalkSpace
    ) : Unit is Ctl + Adj {

        let (node, neighbor, nPositions, nVelocities) = walkState!;
        let eta = 5; // eta is the number of bits used to encode transition probabilities.

        let (mu, epsilon) = neighbor!;

        let n = Length(node![0]!);
        use transitionProbabilityQubitRegister = Qubit[n];
        let transitionProbability = FixedPoint(0, transitionProbabilityQubitRegister);

        within {
            // Compute transition probabilities.
            for lambda in 0 .. (nPositions-1) {
                ControlledOnInt (lambda, ApplyTransitionProbabilitiesBurgers) (mu!, (epsilon, node, lambda, transitionProbability));
            }
        }
        
        apply {
            // Do the amplitude transduction following https://arxiv.org/pdf/1807.03206v2.pdf
            AmplitudeTransduction (transitionProbability);
        }
    }

    operation ApplyTransitionProbabilitiesBurgers (
        epsilon : Qubit,
        node : Node,
        lambda : Int,
        transitionProbability : FixedPoint
    ): Unit is Ctl + Adj {
        let nu = 1.; // numerical value to be changed.
        let deltaL = 1.; // numerical value to be changed.
        let deltaU = 1.; // numerical value to be changed.
        let diffusiveConstant = nu/(deltaL*deltaL);
        let convectiveConstant = deltaU/(4.*deltaL);

        use signNlambda = Qubit();
        let NlambdaMinusOne = node![lambda-1];
        let Nlambda = node![lambda];
        let NlambdaPlusOne = node![lambda+1];
        ComputeSign(Nlambda, signNlambda);

        // case: Nlambda >= 0
        Controlled IncrementDiffusiveTransition ([signNlambda], (Nlambda, diffusiveConstant, transitionProbability)); // minus transition
        X (epsilon); // minus transitions are controlled by epsilon = |O> and plus transitions by epsilon = |1>.
        Controlled IncrementDiffusiveTransition ([signNlambda, epsilon], (NlambdaMinusOne, diffusiveConstant, transitionProbability)); // plus transition
        // case: Nlambda < 0
        X (signNlambda);
        Controlled IncrementDiffusiveTransition ([signNlambda], (NlambdaMinusOne, diffusiveConstant, transitionProbability)); // plus transition
        X (epsilon);
        Controlled IncrementDiffusiveTransition ([signNlambda], (Nlambda, diffusiveConstant, transitionProbability)); // minus transition

        IncrementConvectiveTransition (Nlambda, NlambdaPlusOne, convectiveConstant, transitionProbability);
    }

    operation ComputeSign(
        littleEndianVariable : LittleEndian,
        resultQubit : Qubit
    ): Unit is Ctl + Adj {
        let qubitRegister = littleEndianVariable!;
        CNOT(qubitRegister[0], resultQubit);
    }

    operation IncrementDiffusiveTransition (
        Nlambda : LittleEndian,
        multiplicativeConstant : Double,
        transitionProbability : FixedPoint
    ) : Unit is Ctl + Adj {
        let NlambdaFxP = FixedPoint(0, Nlambda!);
        within {
            MultiplyConstantFxPInPlace(multiplicativeConstant, NlambdaFxP);
        } apply {
            AddFxP(NlambdaFxP, transitionProbability);
        }
    }

    operation IncrementConvectiveTransition (
        Nlambda : LittleEndian,
        NlambdaPlusOne : LittleEndian,
        multiplicativeConstant : Double,
        transitionProbability : FixedPoint
    ) : Unit is Ctl + Adj {
        let (pointPosition , transitionProbabilityQubitRegister) = transitionProbability!;
        let nQubits = Length(transitionProbabilityQubitRegister);
        use additionalTransitionProbabilityQubitRegister = Qubit[nQubits];
        let additionalTransitionProbability = FixedPoint(pointPosition, additionalTransitionProbabilityQubitRegister);
        within {
            let NlambdaFxP = FixedPoint(0, Nlambda!);
            let NlambdaPlusOneFxP = FixedPoint(0, NlambdaPlusOne!);
            SquareFxP(NlambdaFxP, additionalTransitionProbability);
            SquareFxP(NlambdaPlusOneFxP, additionalTransitionProbability);
            MultiplyConstantFxPInPlace(multiplicativeConstant, additionalTransitionProbability);
        } apply {
            AddFxP(additionalTransitionProbability, transitionProbability);
        }
    }

    


}