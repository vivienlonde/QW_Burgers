namespace Estimate {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;

    open AmplitudeTransduction;
    
    operation Resources() : Unit {
        let n = 5;
        use L = Qubit[n];
        use referenceQubitRegister = Qubit[n];
        use flagQubit = Qubit();
        UniformStatePreparationOracleWithSomeAmplitude(L, referenceQubitRegister, flagQubit);
    }

    @EntryPoint()
    operation TestAmplitudeTransduction () : Unit {
        let nbQubits = 5;
        use transitionProbabilityQubitRegister = Qubit[nbQubits];
        let transitionProbability = FixedPoint(0, transitionProbabilityQubitRegister);
        AmplitudeTransduction (transitionProbability);
    }
}
