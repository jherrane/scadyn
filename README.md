# Scattering Dynamics

`scadyn` is a code for scattering dynamics calculations, which utilizes a volume integral equation solution for the T-matrices of non-spherical scatterers. The main motivations for the development of this code are the study of
* grain alignment dynamics in interstellar environments,
* optical tweezers and traps.

A more complete description of the algorithm for scattering dynamics can be found in [Herranen et. al.](https://tuhat.helsinki.fi/portal/services/downloadRegister/90932227/RS2017.pdf) ([DOI:10.1002/2017RS006333](https://dx.doi.org/10.1002/2017RS006333)) and references therein. 

## Installation
**Linux/OSX**
Run `install` if you dare. The installation is dependent on the LAPACK, BLAS, HDF5, and FFTW libraries. If the libraries are found, simply invoking `make` will do. An example run with minimal input can also be found from the installation script.

**Windows**
Library locations and compiler probably need changing. Never tried, so good luck!

## Run
Run install if you dare. All errors should indicate what dependencies I forgot about. The basic command is:
```
./scadyn --mesh shape.h5 -T T.h5 --paramsfile params.in
```
The geometry meshes are to be `tetgen`-compatible. Geometry generation routines are available, and can be run e.g. in `tetgen` or `quartet` mode. Former tends to generate unevenly sized tetrahedra inside the geometry while the latter is much more optimal on the inside, though the surface can be poor. Optimal geometry generation thus depends on the choice of initial surface refinement level, tetrahedralization refinement level and the tetrahedralization engine.
