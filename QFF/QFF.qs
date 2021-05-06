namespace QFF {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Simulation;

    /// # Summary
    /// The first register has size N*log(M),
    /// corresponding to the N positions and the binary encoding of the M velocities.
    /// The second register has size 2M+1,
    /// corresponding to the 2M neighbors of a node + an extra neighbor representing the node itself.
    newtype WalkSpace = (node : Qubit[], neighborIndex : Qubit[]);

    /// # Summary
    /// Inputs : WalkOperator U, a controlRegister |\phi> and a walkState |\psi>.
    /// Output : Applies the operator \sum_j |j><j| \otimes U^j to |\phi> \otimes |\psi>.
    operation Select (WalkOperator : (WalkSpace => Unit is Ctl),
        controlRegister : Qubit[],
        walkState : WalkSpace
    )
    : Unit {
        let nQubits = Length(controlRegister);

        for idxControlQubit in 0 .. nQubits - 1 {
            let control = controlRegister[idxControlQubit];
            let power = 2 ^ ((nQubits - idxControlQubit) - 1);
            let IterateOfWalkOperator = OperationPowC (WalkOperator, power);
            Controlled IterateOfWalkOperator ([control], (walkState));
        }
    }
    
    function Factorial (n : Int) : Int {
        if n==0 {return 1;}
        else {return TimesI(n, Factorial(n-1));}
    } 

    function  BinomialCoefficient (n : Int, k : Int) : Int {
        let numerator = Factorial(n);
        let denominator = TimesI(Factorial(k), Factorial(n-k));
        return DividedByI(numerator, denominator);
    }

    /// # Summary
    /// Computes the first \tau coefficients of the Chebyshev expansion of X^t.
    /// See https://arxiv.org/pdf/1804.02321.pdf (Quantum Fast Forwarding)
    function ComputeChebyshevCoefficients (t : Int, tau : Int) : Double[] {
        mutable chebyshevCoefficients = new Double[tau];
        
        for i in IndexRange(chebyshevCoefficients) {
            if (i > 0 and t == ModI(i,2)) {
                let numerator = IntAsDouble(BinomialCoefficient(t, DividedByI(t-i,2)));
                let denominator = IntAsDouble(PowI(2,t-1));
                set chebyshevCoefficients w/= i <- DividedByD(numerator, denominator);
            }
            elif (i == 0 and t == ModI(0,2)) {
                let numerator = IntAsDouble(BinomialCoefficient(t, DividedByI(t,2)));
                let denominator = IntAsDouble(PowI(2,t));
                set chebyshevCoefficients w/= i <- DividedByD(numerator, denominator);
            }
            else {
                set chebyshevCoefficients w/= i <- 0.; 
            }
        }
        return chebyshevCoefficients;
    }

    /// # Summary
    /// Approximates the evolution of walkstate by WalkOperator^t.
    /// Uses a polynomial of degree tau in WalkOperator to compute the approximation.
    operation QuantumFastForwarding(WalkOperator : (WalkSpace => Unit is Ctl), walkState: WalkSpace, t : Int, tau : Int) : Unit {
        let chebyshevCoefficients = ComputeChebyshevCoefficients(t, tau);
        let nAuxQubitsLCU = DoubleAsInt(Log(IntAsDouble(tau))) + 1; // nAuxQubits = ceil(log_2(tau))
        use auxQubitsLCU = Qubit[nAuxQubitsLCU];
        let auxQubitsLCUlittleEndian = LittleEndian(auxQubitsLCU);
        within {
            PrepareArbitraryStateD (chebyshevCoefficients, auxQubitsLCUlittleEndian);
        } apply {
            Select (WalkOperator, auxQubitsLCU, walkState); // check that it's ok to act on the qubit register as well as on its LittleEndian wrapped version. 
        }
    }

}
