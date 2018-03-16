# Computer simulation of Ulf Dieckmann's model of stationary biological community.

### Algorithm description
Simulations are written using C++11. Data processing is implemented on Python 3.

There are two versions of simulations.

The first one implements the algorithm based on the plane discretization idea. There might be arbitrary number of individuals in every cell of the grid. One step of the algorithm consists of computing birth and death probabilities matrices of individuals. Now birth and death kernels are considered to be Gaussian with zero expectation, arbitrary variance (can be set in the code parameters) and normalizing factor.

The second implementation based on the Poisson process idea. We simulate continiuos process via Poisson process.

### Usage
* To run the first implementation for one dimention case run
```
g++ sim_one_dim.cpp -o sim -std=c++11 -O3
./sim [field size] [discretization] [iterations] [initial population] [birth rate] [death rate] [birth variance] [death variance]"
```
or you can use make file
```
./s1d.x [field size] [discretization] [iterations] [initial population] [birth rate] [death rate] [birth variance] [death variance]"
```

* To run the second implementation (one and two dimensions cases are in one file) download [Boost](http://www.boost.org). 
For windows, download visual studio c++ build tools. Then run
```
cl *.cpp -bigobj -EHsc -I D:\boost_1_65_1 -o simulator.exe
```
For linux, not tested.

* First implementation for two dimentions case doesn't work correctly at the moment (sim_two_dim.cpp).

### Unimplemented
 * Arbitrary birth and death kernels option (done in poisson sim)
 * Multispecies models
