namespace AmplitudeTransduction {
    // follows https://arxiv.org/pdf/1807.03206v2.pdf

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.AmplitudeAmplification;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Convert;

    operation AmplitudeTransduction (
        transitionProbability : FixedPoint  // transitionProbability contains |p_i>
    ) : Unit is Ctl + Adj {
        let (_ , transitionProbabilityQubitRegister) = transitionProbability!;
        let d = Length(transitionProbabilityQubitRegister);
        use referenceQubitRegister = Qubit[d];
        let reference = LittleEndian(referenceQubitRegister);
        use flag = Qubit();

        ApplyToEachCA (H, referenceQubitRegister);

        let GreaterThanStateOracle = StateOracle(GreaterThanOracle);
        
        let nIterations = d; // up to a multiplicative constant. TODO: give the exact value.
        let flagIndex = 2*d;
        let amplifiedComparison = StandardAmplitudeAmplification(nIterations, GreaterThanStateOracle, flagIndex);
        
        let allQubitRegister = transitionProbabilityQubitRegister + referenceQubitRegister + [flag];
        amplifiedComparison(allQubitRegister);

        Adjoint Unif (transitionProbabilityQubitRegister, referenceQubitRegister);

        // At this point referenceQubitRegister is in state |0> (approximately).
        // transitionProbabilityQubitRegister has acquired the amplitude \sqrt(p_i).
        // flag is in state |0> (approximately).
        // -> what happens when auxiliary qubits are only approximately in state |0> ?
    } 

    operation GreaterThanOracle(
        flagIndex : Int, // not used because we assume that the flag qubit is the last qubit of qubitArray.
        qubitArray : Qubit[]
    ) : Unit is Ctl + Adj {
        let flagQubit = Tail(qubitArray);
        let otherQubits = Most(qubitArray);
        let n = Length(otherQubits);
        let transitionProbabilityQubitRegister = otherQubits[0..n/2-1];
        let referenceQubitRegister = otherQubits[n/2..n-1];
        GreaterThan(LittleEndian(transitionProbabilityQubitRegister), LittleEndian(referenceQubitRegister), flagQubit);
    }

    operation Unif(
        transitionProbabilityQubitRegister : Qubit[],
        referenceQubitRegister : Qubit[]
    ) : Unit is Ctl + Adj { 
        //its inverse sends the amplitude on referenceQubitRegisters between indices 0 and 2^n*transitionProbabilityQubitRegister on |0>_{referenceQubitRegister}
        use flagQubit = Qubit();
        let myFlagIndex = Length(transitionProbabilityQubitRegister) + Length(referenceQubitRegister);
        let allQubitRegister = [flagQubit] + transitionProbabilityQubitRegister + referenceQubitRegister;
        
        let statePrepOracle = StateOracle(UniformStatePreparationOracleWithSomeAmplitudeCorrectType);
        let n = Length(referenceQubitRegister);
        let nQueries = n; // up to a multiplicative constant.
        let successMin = IntAsDouble(n); // up to a multiplicative constant.
        let phases = FixedPointReflectionPhases(nQueries, successMin);
        let UnifOnSingleRegister = AmplitudeAmplificationFromStatePreparation(phases, statePrepOracle, myFlagIndex);
        UnifOnSingleRegister(allQubitRegister);
    }

    operation UniformStatePreparationOracleWithSomeAmplitude(
        L : Qubit[],
        referenceQubitRegister : Qubit[],
        flagQubit : Qubit
    ) : Unit is Adj + Ctl {
        ApplyToEachCA (H, referenceQubitRegister);
        GreaterThan(LittleEndian(L), LittleEndian(referenceQubitRegister), flagQubit);
    }
    
    operation UniformStatePreparationOracleWithSomeAmplitudeCorrectType(
        flagIndex : Int, // not used because we assume that the flag qubit is the last qubit of qubitArray.
        qubitRegister : Qubit[]
    ) : Unit is Adj + Ctl {
        let flagQubit = Tail(qubitRegister);
        let otherQubits = Most(qubitRegister);
        let n = Length(otherQubits);
        let L = otherQubits[0..n/2-1];
        let referenceQubitRegister = otherQubits[n/2..n-1];
        UniformStatePreparationOracleWithSomeAmplitude(L, referenceQubitRegister, flagQubit);
    }
    
   

    
}