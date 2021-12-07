namespace Burgers {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;

    open WalkOperations;
    // open ArithmeticOperations;


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
        neighbor : LittleEndian,
        nPositions : Int,
        nVelocities : Int
    ); 

    // Array of length M.
    // For a position i, Node[i] is a LittleEndian register on log(N) qubits encoding a velocity.
    newtype Node = LittleEndian[]; // change to SignedLittleEndian.

    // A value i \in {1, .., M} for the neighbor register
    // corresponds to the node with velocity +1 at position i and -1 at position i+1.
    // A value i \in {M+1, .., 2M} for the neighbor register
    // corresponds to the node with velocity -1 at position i-M and +1 at position i+1-M.
    
    /// # Summary
    /// returns One if the measurement projects to the WalkFlatSubspace.
    /// returns Zero if the measurement projects to a subspace orthogonal to the WalkFlatSubspace.
    operation ProjectToWalkFlatSubspace (
        walkState : WalkSpace
    ) : Result {
        let neighbor = walkState::neighbor;
        use auxQubitFlagFlatSubspace = Qubit();      // use an auxiliary qubit that will be flipped iff neighbor is in the flat subspace.
        ApplyControlledOnInt (0, X, neighbor!, auxQubitFlagFlatSubspace);    // neighbor is in the flat subspace iff it is in state 0 in Little Endian representation.
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
        let nNeighbors = 2*nPositions+1;
        use transitionProbabilitiesQubits = Qubit[eta*nNeighbors*2]; // the factor 2 is here to build the tree of sums.
        let nLevels = Ceiling(Lg(IntAsDouble(nNeighbors)));
        let pointPosition = 0;

        within {
            // Compute transition probabilities.
            for lambda in 0 .. (nPositions-1) {
                // LambdaPlusTransition is: (Nlambda, NlambdaPlusOne) -> (Nlambda + 1, NlambdaPlusOne - 1)
                let transitionProbabilityLambdaPlus = FixedPoint(pointPosition, transitionProbabilitiesQubits[TreeIndicesToArrayRange(nLevels-1, lambda, eta)]);
                // LambdaMinusTransition is: (Nlambda, NlambdaPlusOne) -> (Nlambda - 1, NlambdaPlusOne + 1)
                let transitionProbabilityLambdaMinus = FixedPoint(pointPosition, transitionProbabilitiesQubits[TreeIndicesToArrayRange(nLevels-1, 2*nPositions - lambda , eta)]);
                ApplyTransitionProbabilitiesBurgers (node, lambda, transitionProbabilityLambdaPlus, transitionProbabilityLambdaMinus, nPositions);
            }

            // Compute binary tree made of sums of transition probabilities.
            ComputeBinaryTree(nLevels, eta, pointPosition, transitionProbabilitiesQubits);
        }

        apply {
            // Compute the transitions in amplitude following https://arxiv.org/pdf/0903.3465.pdf
            // we use the register neighbor as the log(d) qubits of the paper.
            // theta_0 = arccos(\sqrt(transitionProbabilitiesTree[0][1]));
            // R (theta_0)(neighbor[0]);
            for level in 1..(nLevels-1) { 
                for i in 0..(2^(level-1)) {
                    let cRegister = transitionProbabilitiesQubits[TreeIndicesToArrayRange(level, 2*i, eta)];
                    let c = FixedPoint(pointPosition, cRegister);
                    let bRegister = transitionProbabilitiesQubits[TreeIndicesToArrayRange(level-1, i, eta)];
                    let b = FixedPoint(pointPosition, bRegister);
                    use thetaOverPiRegister = Qubit[eta]; // we may want something different from eta for the number of bits of precision of theta.
                    let thetaOverPi = FixedPoint(pointPosition, thetaOverPiRegister);
                    within {
                        DetermineAngleCircuit (c, b, thetaOverPi);
                    } apply {
                        for ithBitTheta in 0..(Length(thetaOverPiRegister)-1) {
                            let controlledRotation = ControlledOnInt (i, RFrac);
                            let firstQubits = neighbor![0..(level-1)]; // level also indexes the bits of the register neighbor.
                            let currentQubit = neighbor![level];
                            controlledRotation (firstQubits, (PauliX, 1, 2^(ithBitTheta+1), currentQubit)); 
                        }
                    }
                }
            }
        }

    }

    operation ComputeBinaryTree(
        nLevels : Int,
        eta: Int,
        pointPosition : Int,
        transitionProbabilitiesQubits : Qubit[]
    ): Unit is Adj + Ctl { 
        for level in (nLevels-2)..(-1)..0 {
            let nRegisters = 2^level;
            for i in 0..(nRegisters-1) {

                let firstSummandRange = TreeIndicesToArrayRange(level+1, 2*i, eta);
                let secondSummandRange = TreeIndicesToArrayRange(level+1, 2*i+1, eta);
                let outputRange = TreeIndicesToArrayRange(level, i, eta);

                let firstSummand = FixedPoint(pointPosition, transitionProbabilitiesQubits[firstSummandRange]);
                let secondSummand = FixedPoint(pointPosition, transitionProbabilitiesQubits[secondSummandRange]);
                let output = FixedPoint(pointPosition, transitionProbabilitiesQubits[outputRange]);

                AddFxOutOfplace (firstSummand, secondSummand, output);
            }
        }
    }

    /// # Summary
    /// first: duplicates b with CNOTS.
    /// then: uses an Inplace adder.
    operation AddFxOutOfplace (
        a : FixedPoint,
        b : FixedPoint,
        result : FixedPoint
    ): Unit is Adj + Ctl {
        // assumes that a and b have the same pointPosition.
        // assumes that result is in the AllZero state.
        let (pointPosition, aRegister) = a!;
        let (_, bRegister) = b!;
        let (_, resultRegister) = result!;
        for i in 0..Length(bRegister){
            CNOT(bRegister[i], resultRegister[i]);
        }
        AddFxP(a, result);
    }

    operation ApplyTransitionProbabilitiesBurgers (
        node : Node,
        lambda : Int,
        transitionProbabilityLambdaPlus : FixedPoint,
        transitionProbabilityLambdaMinus : FixedPoint,
        nPositions : Int
    ): Unit is Ctl + Adj {
        let nu = 1.; // numerical value to be changed.
        let deltaL = 1.; // numerical value to be changed.
        let deltaU = 1.; // numerical value to be changed.
        let diffusiveConstant = nu/(deltaL*deltaL);
        let convectiveConstant = deltaU/(4.*deltaL);

        use signAuxiliaryQubit = Qubit();
        // let NlambdaMinusOne = node![lambda-1];
        let Nlambda = node![lambda];
        let NlambdaPlusOne = node![lambda+1];
        ComputeSign(Nlambda, signAuxiliaryQubit);

        // case: Nlambda >= 0
        Controlled IncrementDiffusiveTransition ([signAuxiliaryQubit], (Nlambda, diffusiveConstant, transitionProbabilityLambdaMinus));
        Controlled IncrementDiffusiveTransition ([signAuxiliaryQubit], (Nlambda, diffusiveConstant, transitionProbabilityLambdaMinus));
        X(signAuxiliaryQubit);
        // case: Nlambda < 0
        Controlled IncrementDiffusiveTransition ([signAuxiliaryQubit], (Nlambda, diffusiveConstant, transitionProbabilityLambdaPlus));
        Controlled IncrementDiffusiveTransition ([signAuxiliaryQubit], (Nlambda, diffusiveConstant, transitionProbabilityLambdaPlus));

        IncrementConvectiveTransition (Nlambda, NlambdaPlusOne, convectiveConstant, transitionProbabilityLambdaMinus);
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

    operation MultiplyConstantFxPInPlace(
        a : Double,
        x : FixedPoint
    ) : Unit is Ctl + Adj {
        let (pointPosition, xQubitRegister) = x!;
        let n = Length(xQubitRegister);
        use constantFxPQubitRegister = Qubit[n];
        let constantFxP = FixedPoint(pointPosition, constantFxPQubitRegister);
        within {
            PrepareFxP(a, constantFxP);
        }
        apply {
            MultiplyFxPInPlace(constantFxP, x);
        }
    }


}