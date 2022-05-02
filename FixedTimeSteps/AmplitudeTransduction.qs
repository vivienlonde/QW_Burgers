namespace AmplitudeTransduction {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.AmplitudeAmplification;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Amplitude transduction from a digital oracle.
    /// This version of amplitude transduction outputs amplitudes equal to the square root of the values that the digital oracle computes.
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
        DigitalOracle : (Qubit[], Qubit[]) => Unit is Adj + Ctl,
        outRegister : Qubit[],
        lengthDataRegister : Int
    ) : Unit is Adj + Ctl {

        body (...) {
            use dataRegister = Qubit[lengthDataRegister];
            use referenceRegister = Qubit[lengthDataRegister]; // referenceregister and dataregister have the same length.
            use flagQubit = Qubit();
            AmplitudeTransductionWithAuxQubits (
                DigitalOracle, outRegister, dataRegister, referenceRegister, flagQubit
            );

            let register = dataRegister + referenceRegister + [flagQubit];
            ResetAllWithFailOnZeroResult (register);
        }

        adjoint (...) {
            use dataRegister = Qubit[lengthDataRegister];
            use referenceRegister = Qubit[lengthDataRegister]; // referenceregister and dataregister have the same length.
            use flagQubit = Qubit();
            Adjoint AmplitudeTransductionWithAuxQubits (
                DigitalOracle, outRegister, dataRegister, referenceRegister, flagQubit
            );

            let register = dataRegister + referenceRegister + [flagQubit];
            ResetAllWithFailOnZeroResult (register);
        }

        controlled (cs, ...) {
            use dataRegister = Qubit[lengthDataRegister];
            use referenceRegister = Qubit[lengthDataRegister]; // referenceregister and dataregister have the same length.
            use flagQubit = Qubit();
            Controlled AmplitudeTransductionWithAuxQubits (cs,
                (DigitalOracle, outRegister, dataRegister, referenceRegister, flagQubit)
            );

            let register = dataRegister + referenceRegister + [flagQubit];
            ResetAllWithFailOnZeroResult (register);
        }

        controlled adjoint (cs, ...) {
            use dataRegister = Qubit[lengthDataRegister];
            use referenceRegister = Qubit[lengthDataRegister]; // referenceregister and dataregister have the same length.
            use flagQubit = Qubit();
            Controlled Adjoint AmplitudeTransductionWithAuxQubits (cs,
                (DigitalOracle, outRegister, dataRegister, referenceRegister, flagQubit)
            );

            let register = dataRegister + referenceRegister + [flagQubit];
            ResetAllWithFailOnZeroResult (register);
        }  
    }

    operation AmplitudeTransductionWithAuxQubits (
        DigitalOracle : (Qubit[], Qubit[]) => Unit is Adj + Ctl,
        outRegister : Qubit[],
        dataRegister : Qubit[],
        referenceRegister : Qubit[],
        flagQubit : Qubit
    ) : Unit is Adj + Ctl {

        let lengthOutRegister = Length(outRegister);
        let lengthDataRegister = Length(dataRegister);
        let lengthReferenceRegister = Length(referenceRegister);

        ApplyToEachCA (H, outRegister);

        ApplyToEachCA (H, referenceRegister);

        let MyOperationToAmplify = PartiallySuccessfulScatteredAmplitudeTransduction(_, _,
            DigitalOracle, lengthOutRegister, lengthDataRegister);
        let PartiallySuccessfulScatteredAmplitudeTransductionStateOracle = StateOracle(MyOperationToAmplify);

        // We amplify amplitudes by $\sqrt(2^{lengthOutRegister})$ to counteract the effect of ApplyToEachCA(H, outRegister).
        let flagIndex = lengthOutRegister + lengthDataRegister + lengthReferenceRegister;
        let nIterations = Floor(Sqrt(IntAsDouble(lengthDataRegister)));
        let MyAmplifiedOperation = StandardAmplitudeAmplification(
            nIterations,
            PartiallySuccessfulScatteredAmplitudeTransductionStateOracle,
            flagIndex
        );
        let register = outRegister + dataRegister + referenceRegister + [flagQubit];
        MyAmplifiedOperation(register);
        // following Q# terminology, the prepared state is marked by |1>_{flagQubit}.
        // (whereas in https://arxiv.org/pdf/1807.03206v2.pdf, it is marked by |0>_{flagQubit}.)  

        // The operation Unif groups amplitudes that were scattered over different values of referenceRegister.
        Adjoint Unif(dataRegister, referenceRegister);

        // outRegister is approximately processed according to: |0> -> \sum_i \sqrt(p_i) |i>.
        // dataRegister is approximately in state |0>.
        // referenceRegister is approximately in state |0>.
        // flagQubit is approximately in state |1>.
    }

    internal operation ResetAllWithFailOnOneResult (
        register : Qubit[]
    ) : Unit {
        for q in register {
            if M(q) == One {
                X(q);
                fail "auxiliary qubit measurement output is One: failure.";
            }
        }
    }

    internal operation ResetAllWithFailOnZeroResult (
        register : Qubit[]
    ) : Unit {
        for q in register {
            if M(q) == Zero {
                fail "auxiliary qubit measurement output is Zero: failure.";
            }
            else{
                X(q);
            }
        }
    }

    /// # Summary
    /// Applies amplitude transduction in the subspace flagged by $\ket{1}_{\text{flag}}$ with amplitudes divided by $\sqrt{2^{Length{outRegister}}}$.
    /// Note that amplitudes are scattered over different values of the referenceRegister.
    ///
    internal operation PartiallySuccessfulScatteredAmplitudeTransduction (
        flagIndex : Int,
        globalRegister : Qubit[],
        DigitalOracle : (Qubit[], Qubit[]) => Unit is Adj + Ctl,
        lengthOutRegister : Int,
        lengthDataRegister : Int
    ) : Unit is Adj + Ctl {

        // We grouped all quantum registers in a single globalRegister
        // in order to have the correct signature to call the StandardAmplitudeAmplification function
        // from the Microsoft.Quantum.AmplitudeAmplification namespace.
        // Inside the operation, we split globalRegister back into four registers: 
        let outRegister = globalRegister[0 .. lengthOutRegister-1];
        let dataRegister = globalRegister[lengthOutRegister .. lengthOutRegister+lengthDataRegister-1];
        let referenceRegister = globalRegister[lengthOutRegister+lengthDataRegister .. lengthOutRegister+2*lengthDataRegister-1];
        let flagQubit = globalRegister[flagIndex];

        within {
            DigitalOracle (outRegister, dataRegister);
        } apply {
            // We want to flip flagQubit if referenceRegister \leq dataRegister.
            // In other words if not referenceRegister > dataRegister
            // GreaterThan flips flagQubit if firstRegister > secondRegister.
            GreaterThan(LittleEndian(referenceRegister), LittleEndian(dataRegister), flagQubit);
            X(flagQubit);
        }
    }

    /// # Summary
    /// When dataRegister is in state $\ket{L}$,
    /// $\text{Unif}^{\dagger}$ sends the amplitude of $\ket{i}_{\text{referenceRegister}}$ for $i \in \{0, \cdots, L-1$ on $\ket{0}_{\text{referenceRegister}}$.
    ///
    operation Unif(
        dataRegister : Qubit[],
        referenceRegister : Qubit[]
    ) : Unit is Adj + Ctl {

        body (...) {
            use flagQubit = Qubit();
            UnifWithFlagQubit(dataRegister, referenceRegister, flagQubit);
            ResetAllWithFailOnOneResult ([flagQubit]);
        }

        adjoint (...) {
            use flagQubit = Qubit();
            Adjoint UnifWithFlagQubit(dataRegister, referenceRegister, flagQubit);
            ResetAllWithFailOnOneResult ([flagQubit]);
        }

        controlled (cs, ...) {
            use flagQubit = Qubit();
            Controlled UnifWithFlagQubit (cs, (dataRegister, referenceRegister, flagQubit));
            ResetAllWithFailOnOneResult ([flagQubit]);
        }

        controlled adjoint (cs, ...) {
            use flagQubit = Qubit();
            Controlled Adjoint UnifWithFlagQubit (cs, (dataRegister, referenceRegister, flagQubit));
            ResetAllWithFailOnOneResult ([flagQubit]);
        }
        
    }

    operation UnifWithFlagQubit(
        dataRegister : Qubit[],
        referenceRegister : Qubit[],
        flagQubit : Qubit
    ) : Unit is Adj + Ctl {

        let d = Length(dataRegister); // d = Length(referenceRegister) too.

        // Group registers to comply to the signature of ObliviousAmplitudeAmplificationFromPartialReflections.
        let systemRegister = dataRegister;
        let auxiliaryRegister = [flagQubit] + referenceRegister;
        
        // let UnifPrimeAsStateOracle = StateOracle(UnifPrime(_, _, d));

        // Compute phases.
        let epsilon = 0.0001; // user-defined: may be modified
        let successMin = 1. - epsilon;
        let nQueries = 9; // minimal nQueries depends only on epsilon since Unif needs to be amplified on a range independent of the sizes of qubit registers.
        let phases = FixedPointReflectionPhases(nQueries, successMin);

        let startStateReflection = ReflectionStart(); // reflection about |0,0>_{flagQubit, referenceRegister}.
        let targetStateReflection = TargetStateReflectionOracle(0); // reflection about |1>_{flagQubit}. 0 is the index of flagQubit in auxiliary register
        let signalOracle = ObliviousOracle(UnifPrime);

        let MyAmplifiedOperation = ObliviousAmplitudeAmplificationFromPartialReflections(
            phases, startStateReflection, targetStateReflection, signalOracle
        );
        MyAmplifiedOperation(auxiliaryRegister, systemRegister);

        // flagQubit is approximately in state |1>.
    }

    /// # Summary
    /// Partially applies $\text{Unif}$ in the subspace flagged by $\ket{1}_{\text{auxiliary}}$.
    ///
    operation UnifPrime(
        auxiliaryRegister : Qubit[],
        systemRegister : Qubit[]
    ) : Unit is Adj + Ctl {
        // applies an H gate to the $\lceil \log_2(L) \rceil$ least significant qubits of referenceRegister.
        // L is the integer encoded in dataRegister.

        let d = Length(systemRegister);    
        let flagQubit = auxiliaryRegister[0];
        let referenceRegister = auxiliaryRegister[1..d];
        let dataRegister = systemRegister;

        // the most significant qubit doesn't need auxiliary control qubits.
        Controlled H ([dataRegister[d-1]], referenceRegister[d-1]);

        use controlQubits = Qubit[d-1];
        within {
            ApplyToEachCA(X, controlQubits);
            // all controlQubits are initially on.
            // let's switch off controlQubits when the most significant bits of dataRegister ar all set to zero:
            // controlQubits[i] = dataRegister[d-1] ∨ dataRegister[d-2] ∨ ... ∨ dataRegister[i].
            // i.e. controlQubits[i] = controlQubits[i+1] ∨ dataRegister[i]. 
            if d >= 2 {
                let firstControl = dataRegister[d-1];
                let secondControl = dataRegister[d-2];
                let target = controlQubits[d-2];
                ControlledOnZeroNot ([firstControl, secondControl], target);
            }
            for i in d-3..-1..0 {
                let firstControl = controlQubits[i+1];
                let secondControl = dataRegister[i];
                let target = controlQubits[i];
                ControlledOnZeroNot ([firstControl, secondControl], target);
            }
        } apply {
            // appply Hadamard gates.
            for (control, target) in Zipped(controlQubits, referenceRegister[0..d-2]) {
                Controlled H ([control], target);
            }
            // for i in d-2..-1..0 {
            //     Controlled H ([controlQubits[i]], referenceRegister[i]);
            // }
        }
        
        // we want to flip flagQubit if referenceRegister <= dataRegister.
        // in other words if not referenceRegister > dataRegister
        // GreaterThan flips flagQubit if firstRegister > secondRegister.
        GreaterThan(LittleEndian(referenceRegister), LittleEndian(dataRegister), flagQubit);
        X(flagQubit); // to mark by |1> and not by |0>.
    }

    operation ControlledOnZeroNot (
        controls : Qubit[],
        target : Qubit
    ) : Unit is Adj + Ctl {
        within {
            ApplyToEachCA(X, controls);
        } apply {
            Controlled X (controls, target);
        }
    }

}