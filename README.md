# Computer simulation of Uld Dieckmann's model of stationary biological community.

### Implementation
Simulations are written using C++11. Data processing is implemented on Python 3.

#### Algorithm description
The algorithm is based on idea of coordinate discretization. There might be arbitrary number of individuals in every cell of the grid.
One step of the algorithm consists of computing birth probability matrix of new individuals and possible death of each one. Birth and death kernels are considered to be Gaussian.

### Unimplemented
 * Arbitrary birth and death kernels option
 * Multispecies models
 * Small issues (watch code comments)
