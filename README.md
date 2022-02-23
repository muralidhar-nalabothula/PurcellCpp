# PurcellCpp

PurcellCpp is a C++ code to compute both radiative and non radiative decays rates of dipoles place above an infinitely spread 2D materials or isotropic dielectric slabs or combination of both stacked over each other. 



## Mandatory Requirments

1) A C++ 11 Compiler
2) Flint                
3) Symengine           
5) Eigen 3              ( https://gitlab.com/libeigen/eigen )
6) GMP                  
8) MPFR                 
9) MPC                  
10) Quadpack++          (Download from https://github.com/drjerry/quadpackpp )
11) TCMALLOC 

The perfomance of the code heavly depends on how well Symengine is build.

## Installation

Incase if the above libraries are not present,

**Mac:**

```
Using homebrew: brew install flint, gmp, mpfr, libmpc, gperftools
Using Macports: sudo port install flint, gmp, mpfr, libmpc, gperftools, 
```
for macports:

/path/to/lib is generally "/opt/local/lib"

/path/to/include is generally "/opt/local/include"

**Linux:**

```
sudo apt-get update
sudo apt-get install libflint-dev
sudo apt-get install libgmp-dev
sudo apt-get install libmpfr-dev
sudo apt install libmpc-dev
sudo apt install google-perftools
```



```
export PATH="/path/to/lib:$PATH"
export PATH="/path/to/include:$PATH"
cd directory
unzip symengine-0.8.1.zip -d symengine
unzip eigen-master.zip -d eigen
cd symengine/symengine-0.8.1/
mkdir build && cd build
cmake DWITH_MPFR=off -DWITH_MPC=on -DINTEGER_CLASS=flint CXXFLAGS="-O3 -std=c++11" -DWITH_FLINT:BOOL=on -DCMAKE_INSTALL_PREFIX:PATH="./../../../" -DWITH_COTIRE=off -DWITH_SYMENGINE_THREAD_SAFE=no  ..
make
make install
g++ -O3 -std=c++11 -I/path/to/include -L/path/to/lib -lgmp -lflint -lmpfr -lmpc -ltcmalloc -I./include -L./lib -lsymengine -fPIC -shared purcell.cpp -o libpurcell.so


```
Now include the "purcell.h" in your code and link the compiled libpurcell.so library:
for example: to compile the main program: 
```
g++ -O3 main.cpp -L. -lpurcell -o purcell.exe
./purcell.exe
```
## Usage

The function "Purcell" computes the total, radiative and non-radiative decay (and their numerical integration errors).

To use it, you must include the header "purcell.h" and link the library "libpurcell.so"

**Setting up the system:**

"System" is a datatype which stores the system configuration. It accepts two kinds of elements:
1) "Dielectric_slab" . Dielectric slab takes the following inputs in same order: (dielectric constant, thinkness of dielectric in (meter-1))
2) "Interface" . Interface takes the following input: (dielectric constant of incoming medium, dielectric constant of outgoing medium, sigma_xx, sigma_xy, sigma_yx, sigma_yy)

where sigma is optical conductivity tensor(units Siemens) at the interface.

 Note that the order of the elements is very importrant. The top most elements (closer to dipole ) are arranged to right and bottom elements (far from dipole) are arranged towards left.

For Example, ff you have the following system: dipole ---> dielectric ---->2D material (with 2D surface optical conductivity) ---> Substrate, then you system is 

``` System elements = {Interface(dielectric slab, 2d material and substrate),  Dielectric_slab(), Interface(air dielectric slab) };```

To understand more about setting up the system, refer: https://iopscience.iop.org/article/10.1088/0953-8984/28/37/375802 

**Using Purcell function**

Function Purcell takes the inputs in the following order:

(wavevector (in meter^{-1}); System configuation; distance of dipole from the top interface; compontent of the Greens tensor given in string {xx, xy, xz, ....}

Check the main.cpp for an example on how to calculate

Note: All the units must be given is SI units. 

## Authors and acknowledgment
Please consider citing the following paper if you use the code

Muralidhar Nalabothula, Pankaj K. Jha, Tony Low, and Anshuman Kumar Phys. Rev. B 102, 045416 (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.045416)




