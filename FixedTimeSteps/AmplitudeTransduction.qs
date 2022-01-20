namespace AmplitudeTransduction {
    // follows https://arxiv.org/pdf/1807.03206v2.pdf

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.AmplitudeAmplification;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;


    operation AmplitudeTransduction (
        outRegister : Qubit[],
        lengthDataRegister : Int,
        DigitalOracle : (Qubit[], Qubit[]) => Unit is Ctl + Adj
    ) : Unit is Ctl + Adj {
        let lengthOutRegister = Length(outRegister);
        use dataRegister = Qubit[lengthDataRegister];
        use referenceRegister = Qubit[lengthDataRegister]; // referenceregister and dataregister have the same length.
        use flagQubit = Qubit();

        ApplyToEachCA (H, referenceRegister);

        let MyOperation = PartialStatePreparation(_, _,
            lengthOutRegister, lengthDataRegister, DigitalOracle);
        let MyStateOracle = StateOracle(MyOperation);

        let flagIndex = lengthOutRegister + 2*lengthDataRegister;
        let nIterations = Floor(Sqrt(IntAsDouble(lengthDataRegister)));
        let AmplifiedOperation = StandardAmplitudeAmplification(nIterations, MyStateOracle, flagIndex);
        let register = outRegister + dataRegister + referenceRegister + [flagQubit];
        AmplifiedOperation(register);
        // following Q# terminology, the prepared state is marked by |1>_{flagQubit}.
        // (whereas in https://arxiv.org/pdf/1807.03206v2.pdf, it is marked by |0>_{flagQubit}.)  

        Adjoint Unif(dataRegister, referenceRegister);

        // Summary:
        // outRegister is approximately processed according to |0> -> \sum_i \sqrt(p_i) |i>.
        // dataRegister is approximately in state |0>.
        // referenceRegister is approximately in state |0>.
        // flagQubit is approximately in state |1>.
        // -> what happens when auxiliary qubits are not exactly in state |0> ?
    }

    operation PartialStatePreparation (
        flagIndex : Int,
        register : Qubit[],
        lengthOutRegister : Int,
        lengthDataRegister : Int,
        DigitalOracle : (Qubit[], Qubit[]) => Unit is Ctl + Adj
    ) : Unit is Ctl + Adj {
        let outRegister = register[0 .. lengthOutRegister-1];
        let dataRegister = register[lengthOutRegister .. lengthOutRegister+lengthDataRegister-1];
        let referenceRegister = register[lengthOutRegister+lengthDataRegister .. lengthOutRegister+2*lengthDataRegister-1];
        let flagQubit = register[flagIndex];
        within {
            DigitalOracle (outRegister, dataRegister);
        } apply {
            // we want to flip flagQubit if dataRegister < referenceRegister.
            // GreaterThan flips flagQubit if firstRegister > secondRegister.
            GreaterThan(LittleEndian(referenceRegister), LittleEndian(dataRegister), flagQubit);
        }
    }

    operation Unif(
        //its inverse sends the amplitude of values of referenceRegister between |0> and |value_{dataRegister} - 1> on |0>_{referenceRegister}.
        dataRegister : Qubit[],
        referenceRegister : Qubit[]
    ) : Unit is Ctl + Adj {

        use auxiliaryQubit = Qubit();
        let n = Length(dataRegister); // n = Length(referenceRegister) as well
        let auxiliaryIndex = 2*n;
        let register = dataRegister + referenceRegister + [auxiliaryQubit];
        
        let MyOperation = UnifPrime(_, _, n);
        let MyStateOracle = StateOracle(MyOperation);

        let nQueries = n; // up to a multiplicative constant.
        Message($"number of queries: {nQueries}");
        let epsilon = PowD(2., IntAsDouble(-n)); // up to a multiplicative constant.
        let successMin = 1. - epsilon;
        Message($"Minimum success after amplification: {successMin}");
        let phases = FixedPointReflectionPhases(nQueries, successMin);
        Message($"{phases}");

        let AmplifiedOperation = AmplitudeAmplificationFromStatePreparation(phases, MyStateOracle, auxiliaryIndex);
        AmplifiedOperation(register);

        // auxiliaryQubit is approximately in state |1>.
    }

    operation UnifPrime(
        auxiliaryIndex : Int,
        register : Qubit[],
        lengthDataRegister : Int
    ) : Unit is Adj + Ctl {

        let dataRegister = register[0 .. lengthDataRegister-1];
        let referenceRegister = register[lengthDataRegister .. 2*lengthDataRegister-1];
        let auxiliaryQubit = register[auxiliaryIndex];

        // apply an H gate to the $\lceil \log_2(L) \rceil$ least significant qubits of referenceRegister.
        // L is the integer encoded in dataRegister.
        let n = Length(dataRegister);
        // the most significant qubit doesn't need auxiliary control qubits.
        Controlled H ([dataRegister[0]], referenceRegister[0]);
        use controlQubits = Qubit[n-1];
        within {
            // set control qubits:
            // controlQubits[i] = dataRegister[0] or dataRegister[1] or ... or dataRegister[i] or dataRegister[i+1].
            within {
                X(dataRegister[0]); X(dataRegister[1]);
            } apply {
                Controlled X ([dataRegister[0], dataRegister[1]], controlQubits[0]);
            }
            for i in 2..n-1 {
                within {
                    X(dataRegister[i]);
                } apply {
                    Controlled X ([controlQubits[i-2], dataRegister[i]], controlQubits[i-1]);
                }
            }
        } apply {
            // appply Hadamard gates.
            for i in 1..n-1 {
                Controlled H ([controlQubits[i-1]], referenceRegister[i]);
            }
        }
        
        // we want to flip auxiliaryQubit if dataRegister < referenceRegister.
        // GreaterThan flips auxiliaryQubit if firstRegister > secondRegister.
        GreaterThan(LittleEndian(referenceRegister), LittleEndian(dataRegister), auxiliaryQubit);
    }  
    
}