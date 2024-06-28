# MC_coherent_SU2_BF
An importance sampling Monte Carlo algorithm for computing the SU(2) BF coherent vertex amplitude

This is a simple Monte Carlo algorithm for sampling coherent intertwiners, here applied to computing / approximating the coherent SU(2) BF vertex amplitude.
So, this is a concrete example for sampling 5 intertwiners and use these samples to compute the amplitude. It provides several different sets of boundary data to try,
as well as an algorithm to perform the full calculation (costly in terms of memory and computational time).

The paper for this code can be found here: https://arxiv.org/abs/2403.04836 (accepted for publication in Physical Review D)

The algorithm is fairly simple to run:

File "su2vertex.jl": Contains all necessary definitions of SU(2) representation theory, including the computation of the vertex amplitude, definition of coherent states etc.

File "sampling_coherent_amplitude.jl": Contains the definitions to perform the sampling of intertwiners, including thermalization, sampling and approximation of the vertex amplitude.

File "run_file.jl": Contains a few interesting examples of boundary data for the amplitude to simulate. Boundary states are defined, boundary spins can be chosen.
