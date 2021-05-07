
# Questions

- How do we choose tau, the truncation order in Chebyshev extension? In the QFF paper, it depends on ||D^t|v>||. In this implementation, it is a user-defined parameter for now.

- QFF outputs a quantum state whose amplitudes are a norm-2 normalised vector of probabilities. In other words, they are not square roots of the original probability vector but are directly proportional to the original probabilities.
-> How do we use this quantum state?
