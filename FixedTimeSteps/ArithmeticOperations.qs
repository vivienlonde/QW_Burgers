namespace ArithmeticOperations {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    // open Microsoft.Quantum.Bitwise;
    // open Microsoft.Quantum.Convert;
    // open Microsoft.Quantum.Math;
    // open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;

    ///////////////////////////// Multiplications /////////////////////////

    /// # Summary
    /// a is left unchanged.
    /// b becomes a*b if a!=0.
    /// if a==0, refer to DivideI to understand what happens.
    operation MultiplyFxPInPlace(
        a : FixedPoint,
        b : FixedPoint
    ): Unit is Ctl + Adj {

        let (pointPosition, aQubitRegister) = a!;
        let n = Length(aQubitRegister);
        let (_, bQubitRegister) = b!;

        use auxQubitRegister = Qubit[n];
        let auxFxP = FixedPoint(pointPosition, auxQubitRegister);

        MultiplyFxP(a, b, auxFxP); // auxFxP now contains a*b.
        Adjoint DivideI(LittleEndian(auxQubitRegister), LittleEndian(aQubitRegister), LittleEndian(bQubitRegister)); // reinitializes b to zero.
        SwapFxP(b, auxFxP);
    }

    operation SwapFxP(
        x : FixedPoint,
        y : FixedPoint
    ): Unit is Ctl + Adj {
        let (_, xQubitRegister) = x!;
        let (_, yQubitRegister) = y!;
        for (xQubit, yQubit) in Zipped(xQubitRegister, yQubitRegister) {
            SWAP(xQubit, yQubit);
        }
    }

    operation MultiplyConstantFxP(
        a : Double,
        x : FixedPoint,
        result : FixedPoint
    ): Unit is Ctl + Adj {
        let (pointPosition, xQubitRegister) = x!;
        let n = Length(xQubitRegister);
        use constantFxPQubitRegister = Qubit[n];
        let constantFxP = FixedPoint(pointPosition, constantFxPQubitRegister);
        within {
            PrepareFxP(a, constantFxP);
        }
        apply {
            MultiplyFxP(constantFxP, x, result);
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