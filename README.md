# PurcellCpp

PurcellCpp is a C++ code to compute both radiative and non radiative decays rates of dipoles place above an infinitely spread 2D materials or dielectric slabs or combination of both stacked over each other. 



## Mandatory Requirments

1) A C++ 11 Compiler
2) Flint                (Required for compiling Symengine)
3) Symengine            (Compile with flint integer class for performance)
5) Eigen 3              (Required for finding roots in the code. Can be downloaded from https://gitlab.com/libeigen/eigen )
6) GMP                  (Required for compiling Symengine)
7) GSL                  (Required for compiling Symengine)
8) MPFR                 (Required for compiling Symengine)
9) MPC                  (Required for compiling Symengine)
10) Quadpack++          (Already included. Downloaded from https://github.com/drjerry/quadpackpp )
## Optional Requirments

1) OpenMP
2) TCMALLOC 

The perfomance of the code heavly depends on how well Symengine is build.

## Installation
```
g++ -std=c++11 -O3 -I/path/to/include/dirs -L/path/to/lib/ -lsymengine -lgmp -lflint -lmpfr -lmpc purcell.cpp -o purcell.x
```


## Authors and acknowledgment
Please consider citing the following paper if you use the code

Muralidhar Nalabothula, Pankaj K. Jha, Tony Low, and Anshuman Kumar Phys. Rev. B 102, 045416 (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.045416)




