# Fidelity vs Heavy output
Code used in the paper Fidelity decay and error accumulation in random quantum circuits [arXiv:2404.11444](https://arxiv.org/abs/2404.11444).

This repository serves to give an example of how to compute the average fidelity and heavy output frequency of a random quantum circuit, which consists on the action of (faulty) permutations and then (faulty) two-qubit gates.

The file libs.jl contains all the needed functions, such as functions that evolve the wave vector after a random permutation/two-qubit gate and generate random matrices. The file example.jl contains an example of how to compute the average fidelity/heavy output over the number of layers for fixed parameters $\alpha$ and $p$ in a 1D architecture. The output is a file with the desired quantity and respective statistical uncertainty.

## Further reading

The obtained results and a more detailed description can be found in [Fidelity decay and error accumulation in quantum volume circuits](https://arxiv.org/).
