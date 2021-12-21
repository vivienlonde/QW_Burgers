namespace AmplitudeTransduction {
    // follows https://arxiv.org/pdf/1807.03206v2.pdf

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;

    operation AmplitudeTransduction (
        transitionProbability : FixedPoint  // transitionProbability contains p_i
    ) : Unit is Ctl + Adj {
        let (_ , transitionProbabilityQubitRegister) = transitionProbability!;
        use referenceQubitRegister = Qubit[Length(transitionProbabilityQubitRegister)];
        let reference = LittleEndian(referenceQubitRegister);
        use flag = Qubit();


        ApplyToEachCA (H, referenceQubitRegister);
        
        AmplitudeAmplification (order of sqrt(2M) iterations, COMP) (transitionProbabilityQubitRegister, referenceQubitRegister, flag);

        Adjoint Unif (transitionProbabilityQubitRegister, referenceQubitRegister);


        // At this point referenceQubitRegister is in state |0> (approximately).
        // transitionProbabilityQubitRegister has acquired the amplitude \sqrt(p_i).
        // flag is in state |0> (approximately).
    } 

    operation Unif(
        transitionProbabilityQubitRegister : Qubit[],
        referenceQubitRegister : Qubit[]
    ) : Unit is Adj + Ctl {
        //sends the amplitude on referenceQubitRegisters between indices 0 and 2^n*transitionProbabilityQubitRegister on |0>_{referenceQubitRegister}
    }
   

    
}