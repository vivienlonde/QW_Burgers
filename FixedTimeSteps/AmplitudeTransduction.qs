namespace AmplitudeTransduction {
    // follows https://arxiv.org/pdf/1807.03206v2.pdf

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.AmplitudeAmplification;
    open Microsoft.Quantum.Oracles;

    operation AmplitudeTransduction (
        transitionProbability : FixedPoint  // transitionProbability contains p_i
    ) : Unit is Ctl + Adj {
        let (_ , transitionProbabilityQubitRegister) = transitionProbability!;
        let d = Length(transitionProbabilityQubitRegister);
        use referenceQubitRegister = Qubit[d];
        let reference = LittleEndian(referenceQubitRegister);
        use flag = Qubit();


        ApplyToEachCA (H, referenceQubitRegister);

        let myStateOracle = StateOracle(myOracle);
        
        let nIterations = d; // up to a multiplicative constant. TODO: give the exact value.
        let flagIndex = 2*d+1; // TODO: verifythathtis is the correct value.
        let amplifiedComparison = StandardAmplitudeAmplification (nIterations, myStateOracle, flagIndex);
        
        let allQubitsArray = transitionProbabilityQubitRegister + referenceQubitRegister + [flag];
        amplifiedComparison (allQubitsArray);

        Adjoint Unif (transitionProbabilityQubitRegister, referenceQubitRegister);


        // At this point referenceQubitRegister is in state |0> (approximately).
        // transitionProbabilityQubitRegister has acquired the amplitude \sqrt(p_i).
        // flag is in state |0> (approximately).
    } 

    operation myOracle(
        flagIndex : Int, // not used because we assume that the flag qubit is the last qubit of qubitArray.
        qubitArray : Qubit[]
    ) : Unit is Ctl + Adj {
        let flag = Tail(qubitArray);
        let otherQubits = Most(qubitArray);
        let n = Length(otherQubits);
        let transitionProbabilityQubitRegister = otherQubits[0..n/2-1];
        let referenceQubitRegister = otherQubits[n/2-1..n/2];
        GreaterThan(LittleEndian(transitionProbabilityQubitRegister), LittleEndian(referenceQubitRegister), flag);
    }



    operation Unif(
        transitionProbabilityQubitRegister : Qubit[],
        referenceQubitRegister : Qubit[]
    ) : Unit is Adj + Ctl {
        //sends the amplitude on referenceQubitRegisters between indices 0 and 2^n*transitionProbabilityQubitRegister on |0>_{referenceQubitRegister}
    }
   

    
}