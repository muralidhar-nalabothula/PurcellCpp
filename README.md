# PurcellCpp

PurcellCpp is a C++ code to compute both radiative and non radiative decays rates of dipoles place above an infinitely spread 2D materials or isotropic dielectric slabs or combination of both stacked over each other. 



## Mandatory Requirments

1) A C++ 11 Compiler
2) Flint                (Required for compiling Symengine)
3) Symengine            (Compile with flint integer class for performance)
5) Eigen 3              (Required for finding roots in the code. Already included. Downloaded from https://gitlab.com/libeigen/eigen )
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
cd dir
unzip eigen-master.zip -d eigen
g++ -std=c++11 -O3 -I/path/to/include/dirs -L/path/to/lib/ -lgmp -lflint -lmpfr -lmpc -lgomp -ltcmalloc purcell.cpp -o purcell.x
```

## Usage

The function "Purcell" computes the total, radiative and non-radiative decay (and their numerical integration errors).

To use it, you must include the header "purcell.h" and link the library "libpurcell.so"

******** Setting up the system: *********

"System" is a datatype which stores the system configuration. It accepts two kinds of elements:
1) "Dielectric_slab" . Dielectric slab takes the following inputs in same order: (dielectric constant, thinkness of dielectric in (meter-1))
2) "Interface" . Interface takes the following input: (dielectric constant of incoming medium, dielectric constant of outgoing medium, sigma_xx, sigma_xy, sigma_yx, sigma_yy)

where sigma is optical conductuvy tensor(units Siemens) at the interface.

 Note that the order of the elements is very importrant. The top most elements (closer to dipole ) are arranged to left and bottom elements (far from dipole) are arranged at left.

For Example, ff you have the following system: dipole ---> dielectric ---->2D material (with 2D surface optical conductivity) ---> Substrate, then you system is 

``` System elements = {Interface(dielectric slab, 2d material and substrate),  Dielectric_slab(), Interface(air dielectric slab) };```

To understand more about setting up the system, refer: https://iopscience.iop.org/article/10.1088/0953-8984/28/37/375802 

******** Using Purcell function *********

Function Purcell takes the inputs in the following order:

(wavevector in (meter^{-1}); System configuation; distance of dipole from the top interface; compontent of the Greens tensor given in string {xx, xy, xz, ....}

Check the main.cpp for an example on how to calculate

Note: All the units are must be given is SI units. 

## Authors and acknowledgment
Please consider citing the following paper if you use the code

Muralidhar Nalabothula, Pankaj K. Jha, Tony Low, and Anshuman Kumar Phys. Rev. B 102, 045416 (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.045416)




