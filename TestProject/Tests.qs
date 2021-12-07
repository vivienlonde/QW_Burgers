namespace WalkOperations {
    
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;

    // open ArithmeticOperations;


    @Test("QuantumSimulator") // passes.
    operation TestEquality() : Unit {

        let n=2;
        use firstRegister = Qubit[n];
        use secondRegister = Qubit[n];
        use outputQubit = Qubit();

        Equality (firstRegister, secondRegister, outputQubit);

        AssertMeasurement([PauliZ], [outputQubit], One, "Newly allocated registers must equal.");

        Reset(outputQubit);
        Message("TestEquality passed.");
    }

    // @Test("QuantumSimulator") // fails.
    operation TestDetermineAngleCircuit() : Unit {
        
        let n=2;
        use cQubitRegister = Qubit[n];
        let c = FixedPoint(0, cQubitRegister);
        use bQubitRegister = Qubit[n];
        let b = FixedPoint(0, bQubitRegister);
        use thetaOverPiQubitRegister = Qubit[n];
        let thetaOverPi = FixedPoint(0, thetaOverPiQubitRegister);

        within {
            X(cQubitRegister[1]);
            X(bQubitRegister[0]);
            // c/b = 1/2
            // \sqrt(c/b) = \sqrt(2)/2
            // arccos(\sqrt(c/b)) = pi/4
        } apply {
            DetermineAngleCircuit (c, b, thetaOverPi);
        }

        AssertMeasurement([PauliZ], [thetaOverPiQubitRegister[0]], Zero, "Theta is less than Pi/2");
        AssertMeasurement([PauliZ], [thetaOverPiQubitRegister[1]], One, "Theta is Pi/4");

        ResetAll(thetaOverPiQubitRegister);
        Message("TestDetermineAngleCircuit passed.");
    }

    operation TestSqrtFxP() : Unit {
        
        let n = 5;
        use xQubitRegister = Qubit[n];
        let x = FixedPoint(n, xQubitRegister);
        use resultQubitRegister = Qubit[n];
        let result = FixedPoint(n, resultQubitRegister);

        within {
            X (xQubitRegister[0]);
            X (xQubitRegister[3]);
            X (xQubitRegister[4]);
            // x = 25.
        } apply {
            SqrtFxP (x, result);
            // result = 5.
        }

        AssertMeasurement([PauliZ], [resultQubitRegister[0]], One, "result[0] is 1");
        AssertMeasurement([PauliZ], [resultQubitRegister[1]], Zero, "result[1] is 0");
        AssertMeasurement([PauliZ], [resultQubitRegister[2]], One, "result[2] is 1");
        AssertMeasurement([PauliZ], [resultQubitRegister[3]], Zero, "result[3] is 0");
        AssertMeasurement([PauliZ], [resultQubitRegister[4]], Zero, "result[4] is 0");

        ResetAll(resultQubitRegister);
        Message("TestSqrtFxP passed.");
    }




}
