# Ising2D
This rust code consists of two modules.
In Ising module, the Wolff sampling algorithm for 2-dimensional classical Ising model in a square latticce is implemented.
RBM module represents general restricted Boltzmann machine. Especially, Ising to RBM mapping is implemented.

In the main function, we use both modules to calculate the Fisher information matrix for coherent Gibbs state for 2D Ising models. The result can be found in [http://arxiv.org/abs/1910.11163].
