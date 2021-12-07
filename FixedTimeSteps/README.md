
# Questions

- How do we choose tau, the truncation order in Chebyshev expansion? In the QFF paper, it depends on ||D^t|v>||. In this implementation, it is a user-defined parameter for now.

# Useful Quantum Walk operations

- In WalkOperations.qs, the operation DetermineAngleCircuit is useful for any quantum walk. It relies on many operations in WalkOperations.qs and in ArithmeticOperations.qs.  Would it make sense to add them to Q# as a library?

# Useful QFF operations

- The operations in QFF.qs are useful for any implementation of QFF. Would it make sense to add them to Q# as a library?

-> Add unit tests.