namespace ArithmeticOperations {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Bitwise;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;
    
    ///////////////////////////// Additions ///////////////////////////////

    /// # Summary
    /// first: duplicates b with CNOTS.
    /// then: uses an Inplace adder.
    operation AddFxPOutOfplace (
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

    ///////////////////////////// Numeric Functions /////////////////////////


    // Following https://arxiv.org/pdf/1805.12445.pdf :
    // If x \in [-0.5, 0.5], use a polynomial approximation.
    // If x \in [-1, -0.5] or [0.5, 1], use the identity arccos(x) = pi - 2arccos(sqrt((1-x)/2)).
    operation ApplyArccos (
        x : FixedPoint,
        theta : FixedPoint
    ): Unit is Ctl + Adj {
        // Assumes x \in [-1, 1].
        // Assumes that theta is initialysed to 0.
        // Outputs theta = arccos(x) \in [0, pi].
        use AbsxIsGreaterThanOneHalfFlag = Qubit(); 
        within {
            AbsxGreaterThanOneHalf(x, AbsxIsGreaterThanOneHalfFlag);
        } apply {
            Controlled ApplyArccosByPolynomialEvaluation ([AbsxIsGreaterThanOneHalfFlag], (x, theta));
            within {
                X (AbsxIsGreaterThanOneHalfFlag);
            } apply {
                Controlled ApplyArccosByTransformation ([AbsxIsGreaterThanOneHalfFlag], (x, theta));
            }
        }
    }

    operation AbsxGreaterThanOneHalf (
        x : FixedPoint,
        flag : Qubit
    ): Unit is Adj + Ctl {
        // depends on the two first bits of x:
        // if 00 -> absolute value is less than one half.
        // if 01 -> absolute value is more than one half.
        // if 10 -> absolute value is more than one half.
        // if 11 -> absolute value is less than one half.
        // therefore flag is the XOR of these two bits.
        let (_, xQubitRegister) = x!;
        CNOT(xQubitRegister[0], flag);
        CNOT(xQubitRegister[1], flag);
    }

    // For instance, use the command lolremez -d 5 -r -0.5:0.5 "acos(x)" --progress
    // to obtain the degree 5 approximation:
    // 1.5707963267948966 -1.0001136325718906*x -1.6323064163532001e-1*x**3 -1.0015336213855335e-1*x**5.
    operation ApplyArccosByPolynomialEvaluation (
        x : FixedPoint,
        theta : FixedPoint
    ): Unit is Ctl + Adj {
        let coefficients = [1.5707963267948966, -1.0001136325718906, 0., -0.16323064163532001, 0., -0.10015336213855335]; // from Remez algorithm
        EvaluatePolynomialFxP(coefficients, x, theta);
    }

    // if x \in [-1, -0.5] or [0.5, 1], use the identity arccos(x) = pi - 2arccos(sqrt((1-x)/2)).
    operation ApplyArccosByTransformation (
        x : FixedPoint,
        theta : FixedPoint
    ): Unit is Ctl + Adj {

        let pi = PI();
        let (pointPosition, xQubitRegister) = x!;
        let n = Length(xQubitRegister);

        use sqrtxQubitRegister = Qubit[n];
        let sqrtx = FixedPoint(pointPosition, sqrtxQubitRegister);
        use tmp1QubitRegister = Qubit[n];
        let tmp1 = FixedPoint(pointPosition, tmp1QubitRegister);
        use tmp2QubitRegister = Qubit[n];
        let tmp2 = FixedPoint(pointPosition, tmp2QubitRegister);

        within {
            MultiplyConstantFxP(-0.5, x, tmp1);
            AddConstantFxP(0.5, tmp1);
            SqrtFxP(tmp1, sqrtx);
            ApplyArccosByPolynomialEvaluation(sqrtx, tmp2);
        } apply {
            MultiplyConstantFxP(-2., tmp2, theta);
        }

        AddConstantFxP(pi, theta);
    }

    // x -> (1-x)/2 maps [0.5, 1] to [0, 0.5].
    // and x -> sqrt((1-x)/2) maps [0.5, 1] to [0, sqrt(0.5)].
    // therefore we need to implement the function sqrt on [0, 0.5].
    operation SqrtFxP (
        x : FixedPoint,
        result : FixedPoint
    ): Unit is Ctl + Adj {
        InverseSqrtFxP(x, result);
        MultiplyFxPInPlace(x, result);
    }

    // Following https://arxiv.org/pdf/1805.12445.pdf to implement Newton's method.
    // Since sqrt(x) = x/sqrt(x), it is sufficient to implement the function 1/sqrt on [0, 0.5].
    operation InverseSqrtFxP (
        a : FixedPoint,
        result : FixedPoint
    ): Unit is Ctl + Adj {
        // Assumes a \in [0, 0.5].
        // Assumes that result is initialised to 0 (FixedPoint with the same nb of qubits and point position as a).
        // Assumes that the next operation will be to multiply y by x,
        // therefore, a low-precision output is acceptable when x is close to 0.

        let nbIterations = 3;    // arbitrary value. Defines a precision vs time tradoff.
        let (pointPosition, aQubitRegister) = a!;
        let n = Length(aQubitRegister);

        use xQubitRegister = Qubit[n]; 
        let x = FixedPoint(pointPosition, xQubitRegister);
        use minusOneOverTwoQubitRegister = Qubit[n];
        let minusOneOverTwo = FixedPoint(pointPosition, minusOneOverTwoQubitRegister);
        use minusaOverTwoQubitRegister = Qubit[n];
        let minusaOverTwo = FixedPoint(pointPosition, minusaOverTwoQubitRegister);
        
        within {
            InitialGuessInverseSqrt(a, x);
            PrepareFxP(-0.5, minusOneOverTwo);
            MultiplyFxP(a, minusOneOverTwo, minusaOverTwo);
            // TODO : Prepare minusaOverTwo more efficiently.
        }

        apply {
            use tmpQubitRegister = Qubit[n];
            let tmp = FixedPoint(pointPosition, tmpQubitRegister);

            for _ in 1 .. nbIterations {
                within {
                    SquareFxP(x, tmp);
                    MultiplyFxPInPlace(minusaOverTwo, tmp);
                    AddConstantFxP(1.5, tmp);
                } apply {
                    MultiplyFxPInPlace(tmp, x); 
                }    
            }
        }
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

    operation OppositeFxP(
        x : FixedPoint
    ): Unit is Ctl + Adj {
        let (pointPosition, xQubitRegister) = x!;
        ApplyToEachCA(X, xQubitRegister);
        IncrementByInteger(1, LittleEndian(xQubitRegister));
    }

    operation InitialGuessInverseSqrt(
        a : FixedPoint,
        x : FixedPoint
    ): Unit is Ctl + Adj {
        let (pointPosition, aQubitRegister) = a!;
        let n = Length(aQubitRegister);

        let initialValue = 1.5;
        for i in n-1 .. -1 .. 0 {
            let k = DoubleAsInt( (IntAsDouble(pointPosition) - 2.*IntAsDouble(i)) / 2.); // there must be a better way..
            let shiftedInitialValue = initialValue * IntAsDouble( LeftShiftedI(2, k ) ); 
            PrepareFxP(initialValue, x);
            let b = LeftShiftFxP(a, 3*k-1);
            within {
                OppositeFxP(b);
            } apply {
                AddFxP(b, x);
            }
        }
    }

    function LeftShiftFxP (
        a : FixedPoint,
        shift : Int 
    ): FixedPoint {
        let (pointPosition, aQubitRegister) = a!;
        let newPointPosition = pointPosition + shift;
        return FixedPoint(newPointPosition, aQubitRegister);
    }

    // not used
    operation SignAsFxP (
        x : SignedLittleEndian,
        y : FixedPoint
    ) : Unit is Ctl + Adj {
        let (_, yQubitRegister) = y!;
        X(yQubitRegister[0]);
        let signBit = x!![-1];
        Controlled Z ([signBit], yQubitRegister[0]);
    }


}