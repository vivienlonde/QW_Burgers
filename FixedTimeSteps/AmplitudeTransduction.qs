namespace AmplitudeTransduction {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.AmplitudeAmplification;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Amplitude transduction from a digital oracle.
    ///
    /// # Input
    /// ## outRegister
    /// The digital oracle computes $2^{Length(outRegister)}$ real positive numbers $\alpha_i$ in superposition and writes them in the dataRegister. 
    /// ## lengthDataRegister
    /// Precision of amplitudes is $2^{-\text{lengthDataRegister}}$.
    /// ## DigitalOracle
    /// For state |i> of outRegister, adds $\alpha_i$ to the dataRegister.
    ///
    /// # References
    /// See [ *Y.R. Sanders, G.H. Low, A. Scherer, D.W. Berry* ](https://arxiv.org/pdf/1807.03206v2.pdf)
    operation AmplitudeTransduction (
        outRegister : Qubit[],
        lengthDataRegister : Int,
        DigitalOracle : (Qubit[], Qubit[]) => Unit is Ctl + Adj
    ) : Unit is Ctl + Adj {
        let lengthOutRegister = Length(outRegister);
        use dataRegister = Qubit[lengthDataRegister];
        use referenceRegister = Qubit[lengthDataRegister]; // referenceregister and dataregister have the same length.
        use flagQubit = Qubit();

        ApplyToEachCA (H, outRegister);

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
        // -> how can we release auxiliary qubits that are not exactly in state |0>
        //    and still have an is Ctl + Adj operation ?
    }

    /// # Summary
    /// Partially applies amplitude transduction in the subspace flagged by $\ket{1}_{\text{flag}}$.
    ///
    internal operation PartialStatePreparation (
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

    /// # Summary
    /// When dataRegister is in state $\ket{L}$,
    /// $\text{Unif}^{\dagger}$ sends the amplitude of $\ket{i}_{\text{referenceRegister}}$ for $i \in \{0, \cdots, L-1$ on $\ket{0}_{\text{referenceRegister}}$.
    ///
    internal operation Unif(
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

    /// # Summary
    /// Partially applies $\text{Unif}$ in the subspace flagged by $\ket{1}_{\text{auxiliary}}$.
    ///
    internal operation UnifPrime(
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