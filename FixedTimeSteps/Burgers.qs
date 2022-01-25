namespace Burgers {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;

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
        neighbor : Neighbor
    ); 

    ////////////////// Node ////////////////////
    // Array of length M.
    // For a position lambda, Node[Lambda] is a SignedLittleEndian register on log(N) qubits encoding a velocity.
    newtype Node = SignedLittleEndian[];
    // the sign bit is the last one.

    ///////////////// Neighbor /////////////////
    // the first element (of type LittleEndian) indexes a position lambda ($0 \leq \lambda \leq M-1$).
    // the second element (of type Qubit) discriminates between the two following neghbors:
    //      |0> corresponds to the (lambda, +) neighbor:
    //          (Nlambda, NlambdaPlusOne) -> (Nlambda + 1, NlambdaPlusOne - 1)
    //      |1> corresponds to the (lambda, -) neighbor:
    //          (Nlambda, NlambdaPlusOne) -> (Nlambda - 1, NlambdaPlusOne + 1)
    newtype Neighbor = (LittleEndian, Qubit);

    /// # Summary
    /// returns One if the measurement projects to the WalkFlatSubspace.
    /// returns Zero if the measurement projects to a subspace orthogonal to the WalkFlatSubspace.
    ///
    operation ProjectToWalkFlatSubspace (
        walkState : WalkSpace
    ) : Result {
        let neighbor = walkState::neighbor;
        let (mu, _) = neighbor!;
        use auxQubitFlagFlatSubspace = Qubit();      // use an auxiliary qubit that will be flipped iff neighbor is in the flat subspace.
        ApplyControlledOnInt (0, X, mu!, auxQubitFlagFlatSubspace);    // a walkstate is in the flat subspace iff its neighbor register is in state 0 in Little Endian representation.
        return MResetZ(auxQubitFlagFlatSubspace);
    }

    /// # Summary
    /// reflects about the walk flat subspace.
    ///
    operation ReflectAboutWalkFlatSubspace (
        walkState : WalkSpace
    ) : Unit is Ctl + Adj {
        body (...) {
            let neighbor = walkState::neighbor;
            let (mu, _) = neighbor!;
            ReflectAboutInteger (0 , mu);
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
            let (node, neighbor) = walkState!;
            let (mu, epsilon) = neighbor!;
            let nPositions = Length(node!)-1; 

            for lambda in 1 .. nPositions { // Index 0 corresponds to the flat subspace.
                let Nlambda = node![lambda];
                let NlambdaPlusOne = node![lambda+1];
            ///// Operations controlled on epsilon = |1> and mu = |lambda>:
                // modify the node register:
                Controlled ApplyControlledOnInt ([epsilon], (lambda, PlusOne, mu!, Nlambda));
                Controlled ApplyControlledOnInt ([epsilon], (lambda, MinusOne, mu!, NlambdaPlusOne));
                // modify the neighbor register:
                Controlled ApplyControlledOnInt ([epsilon], (lambda, X, mu!, epsilon)); // exchanges the +1 -1 transition and the -1 +1 transition by flipping epsilon.
            ///// Operations controlled on epsilon = |0> and mu = |lambda>:
                within {X (epsilon);
                } apply {    
                    // modify the node register:
                    Controlled ApplyControlledOnInt ([epsilon], (lambda, MinusOne, mu!, Nlambda));
                    Controlled ApplyControlledOnInt ([epsilon], (lambda, PlusOne, mu!, NlambdaPlusOne));
                    // modify the neighbor register:
                    Controlled ApplyControlledOnInt ([epsilon], (lambda, X, mu!, epsilon)); // exchanges the +1 -1 transition and the -1 +1 transition by flipping epsilon.
                }
            }
            ///// The self-transition edge is left unchanged.
        }
        adjoint self;
    }

    operation PlusOne(
        z : SignedLittleEndian
    ) : Unit is Ctl + Adj {
        IncrementByInteger (1, z!);
    }

    operation MinusOne(
        z : SignedLittleEndian
    ) : Unit is Ctl + Adj {
        IncrementByInteger (-1, z!);
    }

    operation CoinTossBurgers (
        walkState : WalkSpace
    ) : Unit is Ctl + Adj {

        let (node, neighbor) = walkState!;
        let (mu, epsilon) = neighbor!;

        let outRegister = mu! + [epsilon]; // same information as "neighbor" but with type Qubit[].

        let n = Length(node![0]!!); // n \approx log(N): number of bits used to encode velocities.
        let lengthDataRegister = n; // number of bits used to encode transition probabilities.
        let DigitalOracle = WrappedApplyAllTransitions(node, _, _); // DigitalOracle is called $\mathfrak{amp}$ in https://arxiv.org/pdf/1807.03206v2.pdf

        AmplitudeTransduction (outRegister, lengthDataRegister, DigitalOracle);

        // within {
        //     // Create a superposition of all neighbors (assuming that the qubits are all in state |0>):
        //     ApplyToEachCA(H, mu!);
        //     H (epsilon);
        //     // Compute transition probabilities.
        //     for lambda in 0 .. (nPositions-1) {
        //         ControlledOnInt (lambda, ApplyOneTransition) (mu!, (epsilon, node, lambda, transitionProbability));
        //     }
        // }
        
        // apply {
        //     // Do the amplitude transduction following https://arxiv.org/pdf/1807.03206v2.pdf
        //     AmplitudeTransduction (outRegister, lengthDataRegister, DigitalOracle);
        // }
    }

    operation WrappedApplyAllTransitions (
        // additional input
        node : Node,
        // inputs used to call AmplitudeTransduction.
        outRegister : Qubit[], // corresponds to the neigborRegister, i.e. to mu and epsilon.
        dataRegister : Qubit[] // corresponds to the transitionProbability register.
    ): Unit is Ctl + Adj {

        let epsilon = Tail(outRegister);
        let mu = Most(outRegister);
        let neighbor = Neighbor(LittleEndian(mu), epsilon);

        let transitionProbability = FixedPoint(0, dataRegister);

        ApplyAllTransitions (node, neighbor, transitionProbability);
    }


    operation ApplyAllTransitions (
        node : Node,
        neighbor : Neighbor,
        transitionProbability : FixedPoint
    ): Unit is Ctl + Adj {
        let (mu, epsilon) = neighbor!;
        let nPositions = Length(node!)-1;

        for lambda in 0 .. nPositions {
            // lambda (Int) controls mu (Qubit[]).
            let LambdaTransition = ControlledOnInt (lambda, ApplyOneTransition) (mu!, (lambda, _, _, _)); 
            for eta in [false, true] {
                within {
                    if not eta {X(epsilon);} // eta (Bool) controls epsilon (Qubit). 
                } apply {
                    let LambdaEtaTransition = Controlled LambdaTransition ([epsilon], (eta, _, _));
                    LambdaEtaTransition (node, transitionProbability);
                }
            }
        }
    }

    operation ApplyOneTransition (
        lambda : Int,
        eta : Bool,
        node : Node,
        transitionProbability : FixedPoint
    ): Unit is Ctl + Adj {
        let nu = 1.; // numerical value to be changed.
        let deltaL = 1.; // numerical value to be changed.
        let deltaU = 1.; // numerical value to be changed.
        let diffusiveConstant = nu/(deltaL*deltaL);
        let convectiveConstant = deltaU/(4.*deltaL);

        let NlambdaMinusOne = node![lambda-1];
        let Nlambda = node![lambda];
        let NlambdaPlusOne = node![lambda+1];
        let signNlambda = Tail(Nlambda!!);

        let CtrldOnSignNlambdaDiffusiveTransition = Controlled DiffusiveTransition ([signNlambda], (diffusiveConstant, _, _));

        // case: Nlambda >= 0 
        if eta { // (lambda-1, +) transition:
            CtrldOnSignNlambdaDiffusiveTransition (NlambdaMinusOne, transitionProbability); 
        } else { // (lambda, -) transition:
            CtrldOnSignNlambdaDiffusiveTransition (Nlambda, transitionProbability);
        }

        // case: Nlambda < 0
        within {
            X (signNlambda);
        } apply {
            if eta { // (lambda, +) transition:
                CtrldOnSignNlambdaDiffusiveTransition (Nlambda, transitionProbability);
            } else { // (lambda-1, -) transition:
                CtrldOnSignNlambdaDiffusiveTransition (NlambdaMinusOne, transitionProbability);
            }
        }
        
        ConvectiveTransition (convectiveConstant, Nlambda, NlambdaPlusOne, transitionProbability);
    }

    operation DiffusiveTransition (
        diffusiveConstant : Double,
        Nlambda : SignedLittleEndian,
        transitionProbability : FixedPoint
    ) : Unit is Ctl + Adj {
        let AbsNlambdaFxP = FixedPoint(0, Most(Nlambda!!));
        within {
            MultiplyConstantFxPInPlace(diffusiveConstant, AbsNlambdaFxP);
        } apply {
            AddFxP(AbsNlambdaFxP, transitionProbability);
        }
    }

    operation ConvectiveTransition (
        convectiveConstant : Double,
        Nlambda : SignedLittleEndian,
        NlambdaPlusOne : SignedLittleEndian,
        transitionProbability : FixedPoint
    ) : Unit is Ctl + Adj {
        let (pointPosition , transitionProbabilityRegister) = transitionProbability!;
        let nQubits = Length(transitionProbabilityRegister);
        use additionalTransitionProbabilityRegister = Qubit[nQubits];
        let additionalTransitionProbability = FixedPoint(pointPosition, additionalTransitionProbabilityRegister);
        within {
            let AbsNlambdaFxP = FixedPoint(0, Most(Nlambda!!));
            let AbsNlambdaPlusOneFxP = FixedPoint(0, Most(NlambdaPlusOne!!));
            SquareFxP(AbsNlambdaFxP, additionalTransitionProbability);
            SquareFxP(AbsNlambdaPlusOneFxP, additionalTransitionProbability);
            MultiplyConstantFxPInPlace(convectiveConstant, additionalTransitionProbability);
        } apply {
            AddFxP(additionalTransitionProbability, transitionProbability);
        }
    }

    


}