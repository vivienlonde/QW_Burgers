namespace QFF {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Simulation;
    open Burgers; // define in this file a generic type instead  of the Burgers specific WalkSpace type defined in file Burgers.qs


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

    operation ProjectToLCUFlatSubspace (walkState : WalkSpace) : Result {

        return One; // return Zero if the measurement projects to a subspace orthogonal to the LCUFlatSubspace. 
    }
    
    function Factorial (n : Int) : Int {
        if n==0 {return 1;}
        else {return n*Factorial(n-1);}
    } 

    function  BinomialCoefficient (n : Int, k : Int) : Int {
        return Factorial(n)/(Factorial(k)*Factorial(n-k));
    }

    /// # Summary
    /// Computes the first \tau coefficients of the Chebyshev expansion of X^t.
    /// See https://arxiv.org/pdf/1804.02321.pdf (Quantum Fast Forwarding)
    function ComputeChebyshevCoefficients (t : Int, tau : Int) : Double[] {
        mutable chebyshevCoefficients = new Double[tau];
        
        for i in IndexRange(chebyshevCoefficients) {
            if (i > 0 and t == ModI(i,2)) {
                let numerator = IntAsDouble(BinomialCoefficient(t, t-i/2));
                let denominator = IntAsDouble(2^(t-1));
                set chebyshevCoefficients w/= i <- numerator/denominator;
            }
            elif (i == 0 and t == ModI(0,2)) {
                let numerator = IntAsDouble(BinomialCoefficient(t, t/2));
                let denominator = IntAsDouble(2^t);
                set chebyshevCoefficients w/= i <- numerator/denominator;
            }
            else {
                set chebyshevCoefficients w/= i <- 0.; 
            }
        }
        return chebyshevCoefficients;
    }

    /// # Summary
    /// Approximates the evolution of walkstate by WalkOperator^t.
    /// Uses a polynomial of degree \tau in WalkOperator to compute the approximation.
    /// \tau is in \Theta(\sqrt(t)).
    operation QuantumFastForwarding (WalkOperator : (WalkSpace => Unit is Ctl),
        walkState: WalkSpace,
        t : Int,
        tau : Int
    )
    : String {
        let chebyshevCoefficients = ComputeChebyshevCoefficients(t, tau);
        let nAuxQubitsLCU = DoubleAsInt(Log(IntAsDouble(tau+1)))+1;    // nAuxQubits = ceil(log_2(tau+1))
        // we need an extra state (tau+1) to represent the LCUFlatSubspace.
        use auxQubitsLCU = Qubit[nAuxQubitsLCU];
        let auxQubitsLCUlittleEndian = LittleEndian(auxQubitsLCU);

        within {
            PrepareArbitraryStateD (chebyshevCoefficients, auxQubitsLCUlittleEndian);
        } apply {
            Select (WalkOperator, auxQubitsLCU, walkState);
            // TODO: check that it works well to act on the qubit register with one operation 
            // and on its LittleEndian wrapped version with another operation. 
        }

        // Test whether we are in the flat subspaces :
        // - for the LCU (i.e. auxQubitsLCU is in state |0>)
        // - and for the Walkoperator (i.e. walkstate has its neighborIndex register in state |0>).
        // If both conditions are satisfied, output "Success".
        // Otherwise, output "Failure".
        let inLCUFlatSubspace = ProjectToLCUFlatSubspace (walkState);
        let inWalkFlatSubspace = ProjectToWalkFlatSubspace (walkState);
        if  (inLCUFlatSubspace==Zero or inWalkFlatSubspace==Zero) {
            return "Failure";
        } else {
            return "Success";
        } 
    }

}
