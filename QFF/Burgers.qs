namespace Burgers {

    /// # Summary
    /// The first register has size N*log(M),
    /// corresponding to the N positions and the binary encoding of the M velocities.
    /// The second register has size 2M+2,
    /// corresponding to the 2M neighbors of a node
    /// + an extra state representing the node as a neighbor to itself (for self transitions). By convention it is state |2M+1>. 
    /// + an extra state representing the node (all other states represent directed edges). By convention it is state |0>.
    /// This last state flags what is also called the flat subspace.
    newtype WalkSpace = (node : Qubit[], neighborIndex : Qubit[]);

    operation ProjectToWalkFlatSubspace (walkState : WalkSpace) : Result {

        return One; // return Zero if the measurement projects to a subspace orthogonal to the WalkFlatSubspace. 
    }

    operation ReflectAboutWalkFlatSubspace (walkState : WalkSpace) : Unit is Ctl {
        // https://docs.microsoft.com/en-us/qsharp/api/qsharp/microsoft.quantum.arithmetic.reflectaboutinteger
    }

    operation ShiftBurgers (walkState : WalkSpace) : Unit is Ctl {
        // Note that ShiftBurgers is self adjoint.

    }

    operation CoinTossBurgers (walkState : WalkSpace) : Unit is Ctl + Adj {

    }

    operation WalkOperatorBurgers (walkState : WalkSpace) : Unit is Ctl {
        within {
            CoinTossBurgers(walkState);
        } apply {
            ShiftBurgers(walkState);
        }
        ReflectAboutWalkFlatSubspace(walkState);
    }


}