namespace Estimate {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    open AmplitudeTransduction;

    operation MyDigitalOracle (outRegister: Qubit[], dataRegister: Qubit[]) : Unit is Adj + Ctl {
        // data encoded on 3 bits -> Length(dataRegister) = 3.
        // distribution on 4 = 2^2 elements -> Length(outRegister) = 2.
        // alpha_0 = 0.125   -> LE: 100
        // alpha_1 = 0.5     -> LE: 001
        // alpha_2 = 0.      -> LE: 000
        // alpha_3 = 0.375   -> LE: 110
        // We consider that outRegister is LE encoded as well:
        // 00 : alpha_0
        // 10 : alpha_1
        // 01 : alpha_2
        // 11 : alpha_3

        // let probabilityDistribution = [0.125, 0.5, 0., 0.375];
        let probabilityDistributionAsIntegers = [1, 4, 0, 3]; // x = p * 2^(lengthDataRegister)
        for x in probabilityDistributionAsIntegers {
            ControlledOnInt(x, IncrementByInteger(x, _)) (outRegister, LittleEndian(dataRegister));
        }
    }

    operation TestAmplitudeTransduction () : Unit {
        let nOutQubits = 1;
        use outRegister = Qubit[nOutQubits];
        let lengthDataRegister = 1;
        AmplitudeTransduction (MyDigitalOracle, outRegister, lengthDataRegister);
        ResetAll(outRegister);
    }

    @EntryPoint()
    operation TestUnifPrime() : Unit {
        // allocate qubits.
        let nPrecisionBits = 3;
        use dataRegister = Qubit[nPrecisionBits];
        use referenceRegister = Qubit[nPrecisionBits];
        use flagQubit = Qubit();

        // set data value.
        let classicalValues = [false, false, true];
        for (classicalValue, dataQubit) in Zipped(classicalValues, dataRegister){
            if classicalValue {
                X(dataQubit);
            }
        }

        // apply UnifPrime and Assert referenceRegister state.
        UnifPrime([flagQubit] + referenceRegister, dataRegister);

        // // do the measurement tests.
        // mutable HgateApplied = false;
        // for i in nPrecisionBits-1..-1..0 {
        //     if classicalValues[i] {
        //         set HgateApplied = true;
        //     }
        //     if HgateApplied {
        //         AssertMeasurement([PauliX], [referenceRegister[i]], Zero, $"Expected |+> state for referenceRegister[{i}].");
        //     } else {
        //         AssertMeasurement([PauliZ], [referenceRegister[i]], Zero, $"Expected |0> state for referenceRegister[{i}].");
        //     }
        // }
        // DumpMachine("dump.txt");
        
        // reset every qubit.
        ResetAll(dataRegister); ResetAll(referenceRegister); Reset(flagQubit);
    }

    operation WrappedUnifPrime (register : Qubit[]) : Unit is Adj + Ctl {
        let n = Length(register);
        Message($"n = {n}");
        let d = n / 2;
        Message($"d = {d}");
        let auxiliaryRegister = register[0..d];
        let systemRegister = register[d+1..n-1];
        UnifPrime(auxiliaryRegister, systemRegister);
    } 

    operation WrappedUnifWithFlagQubit (register : Qubit[]) : Unit is Adj + Ctl {
        let n = Length(register);
        let d = n / 2;
        let dataRegister = register[0..d-1];
        let referenceRegister = register[d..2*d-1];
        let auxiliaryQubit = register[n-1];
        UnifWithFlagQubit(dataRegister, referenceRegister, auxiliaryQubit);
    }

    operation WrappedUnif (register : Qubit[]) : Unit is Adj + Ctl {
        let n = Length(register);
        let d = n / 2;
        let dataRegister = register[0..d-1];
        let referenceRegister = register[d..n-1];
        Unif(dataRegister, referenceRegister);
    }

    operation WrappedATWithAuxQubits (register : Qubit[]) : Unit is Adj + Ctl {
        let lengthOutRegister = 2;
        let lengthDataRegister = 3;
        let lengthReferenceRegister = 3;
        let outRegister = register[0..lengthOutRegister-1];
        let dataRegister = register[lengthOutRegister..lengthOutRegister+lengthDataRegister-1];
        let referenceRegister = register[lengthOutRegister+lengthDataRegister..lengthOutRegister+lengthDataRegister+lengthReferenceRegister-1];
        let flagQubit = register[lengthOutRegister+lengthDataRegister+lengthReferenceRegister];
        
        AmplitudeTransductionWithAuxQubits (
            MyDigitalOracle,
            outRegister,
            dataRegister,
            referenceRegister,
            flagQubit
        );
    }

    operation WrappedAT (outRegister : Qubit[]) : Unit is Adj + Ctl {
        let lengthDataRegister = 3;
        AmplitudeTransduction(MyDigitalOracle, outRegister, lengthDataRegister);
    }

}
