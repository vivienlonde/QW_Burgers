namespace Walkoperations {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Math;

    open ArithmeticOperations;


    function TreeIndicesToArrayRange (
        level : Int,
        horizontalIndex : Int,
        eta : Int
    ): Range {
        let leftValue = 2^level - 1 + horizontalIndex;
        let rightValue = leftValue + eta - 1;
        return leftValue..rightValue;
    }

    // theta = arccos(sqrt(c/b))
    operation DetermineAngleCircuit (
        c : FixedPoint,
        b : FixedPoint,
        thetaOverPi : FixedPoint
    ): Unit is Ctl + Adj {
        let (_, cRegister) = c!;
        let (_, bRegister) = b!;
        use flagRegister = Qubit[2];
        SpecialCaseCircuit (cRegister, bRegister, flagRegister);
        // flagRegister[1] = |1> iff c = 0.
        // flagRegister[0] = |1> iff c = b.
        within {
            X (flagRegister[1]); // if flagRegister was |?1> after SCcircuit, do nothing.
        } apply {
            Controlled ApplyPiOverTwo (flagRegister, (thetaOverPi)); // if flagRegister was |10> after SCcircuit, θ = π/2, as c = 0.
            within {
                X (flagRegister[0]); // if flagRegister was |00> after SCcircuit, compute θ normally.
            } apply {
                Controlled ThetaCircuit (flagRegister, (c, b, thetaOverPi));
            }
        }
    }

    operation ApplyPiOverTwo (
        thetaOverPi : FixedPoint
    ) : Unit is Ctl + Adj {
        // Assumes that theta contains the decimal value θ/π.
        // Assumes that thetaOverPiRegister is initially in state |0..0>.
        let (pointPosition, thetaOverPiRegister) = thetaOverPi!;
        if pointPosition==0 {
            X (thetaOverPiRegister[0]);
        }
        else {
            fail "angle theta is larger than pi";
        }
    }

    operation SpecialCaseCircuit (
        cRegister : Qubit[],
        bRegister : Qubit[],
        flagRegister : Qubit[]
    ): Unit is Ctl + Adj {
        let n = Length(cRegister);
        use zeroRegister = Qubit[n];
        Equality (cRegister, zeroRegister, flagRegister[0]);
        Equality (cRegister, bRegister, flagRegister[1]);
    }

    operation Equality (
        firstRegister : Qubit[],
        secondRegister : Qubit[],
        outputQubit : Qubit
    ): Unit is Ctl + Adj {
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

    // theta = arccos(sqrt(c/b))
    operation ThetaCircuit (
        c : FixedPoint,
        b : FixedPoint,
        thetaOverPi : FixedPoint
    ): Unit is Ctl + Adj {

        let (pointPosition, cRegister) = c!;
        let (_, bRegister) = b!;
        let n = Length(cRegister);

        use temp1 = Qubit[n];
        use temp2 = Qubit[n];
        use thetaRegister = Qubit[n];
        let theta = FixedPoint(pointPosition, thetaRegister);

        within {
            DivideI (LittleEndian(cRegister), LittleEndian(bRegister), LittleEndian(temp1)); 
            SqrtFxP (FixedPoint(pointPosition, temp1), FixedPoint(pointPosition, temp2));
            ApplyArccos (FixedPoint(pointPosition, temp2), theta);
        } apply {
            MultiplyConstantFxP(1./PI(), theta, thetaOverPi);
        }
    }

    

}