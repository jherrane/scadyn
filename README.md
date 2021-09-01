# Scattering Dynamics

`scadyn` is a code for scattering dynamics calculations, which utilizes a volume integral equation solution for the _T_-matrices of non-spherical scatterers. The main motivations for the development of this code are the study of

* grain alignment dynamics in interstellar environments,
* optical tweezers and traps.

A more complete description of the algorithm for scattering dynamics can be found in [Herranen et. al.](https://tuhat.helsinki.fi/portal/services/downloadRegister/90932227/RS2017.pdf) ([DOI:10.1002/2017RS006333](https://dx.doi.org/10.1002/2017RS006333)) and references therein. 

## Installation
**Linux/OSX**
Run `install` if you dare. The installation is dependent on the LAPACK, BLAS, HDF5, and FFTW libraries. If the libraries are found, simply invoking `make` will do. An example run with minimal input can also be found from the installation script.

**Windows**
Library locations and compiler probably need changing. Never tried, so good luck!

### PyMesh

The farce also known as maintaining of PyMesh (for further confusion, the actual package is `pymesh2` in package maintaining services) makes it really annoying to use on contemporary systems with `python>v3.6`. I got it working using conda.

With conda installed, create an Python 3.6 environment via

`conda create -n py36 python=3.6 anaconda`

then activate the new environment with

`conda activate py36`

Finally, we need to install two necessary libraries, PyMesh and MeshIO, with

`conda install -c conda-forge pymesh2 meshio`

Nowhere do we have official instructions that work, so just trust me and use this! Given that we can still create an 3.6 environment, this approach should work till the end of time.

## Run
Run install if you dare. All errors should indicate what dependencies I forgot about. The basic command is:
```
./scadyn --mesh shape.h5 -T T.h5 --paramsfile params.in
```
## Description of the pipeline

### Geometry generation

The geometry meshes are to be `tetgen`-compatible volumetric tetrahedral meshes. Geometry generation routines are available, and can be run e.g. in `tetgen` or `quartet` mode. Former tends to generate unevenly sized tetrahedra inside the geometry while the latter is much more optimal on the inside, though the surface can be poor. Optimal geometry generation thus depends on the choice of initial surface refinement level, tetrahedralization refinement level and the tetrahedralization engine.

All in all, any `tetgen`-compatible software should work. The package installation may be a drag, `PyMesh` is actually `PyMesh2` in the `pip`-repository. Further, things tend to break between releases of `PyMesh`, making included routines readily obsolete. 

The `.h5` files are expected to contain the `tetgen` vertices as dataset called `coord` and the edge data as `etopol` (edge topology). Further, the mesh should contain complex permittivity for each tetrahedron, as datasets `param_r` and `param_i` (keep in mind that the tetrahedra are defined by `etopol`).

### _T_-matrix calculation

Usually the main bulk of computational efforts go into calculating the _T_-matrix given a mesh, and a `params.in` file describing the radiative environment. The code was initially serialized to a dumb degree, and most of the inconveniences and bugs relate to workarounds to the initial ~~problems~~ features.

In `params.in`, any intrinsic property of the scatterer must be regarded with caution when using a precalculated _T_-matrix. As of now, only the density of the scatterer can be altered after precalculation. The same goes for the available wavelengths in the radiative environment. 

There are multiple ways to construct a radiative environment. In simplest of cases, where only a single wavelength is considered, `bars` = 1 and `lambda1` = `lambda2`, where the lambdas define the desired wavelengths, the code will happily produce a _T_-matrix with no problem. After there is a _T_--matrix, and the wavelength setup is altered, all hell will probably break loose.

To this day there still is no parallelized computational setup, neither when calculating columns of a single _T_-matrix nor when calculating separate _T_-matrices in a `bars`-discretized wavelength range. To combat the latter problem, there is a possibility to set a `whichbar` flag to the desired wavelength and perform the calculations in separate processes. The first calculation will create an empty _T_-matrix, to which each process with a certain `whichbar` knows to replace the correct elements. Now, any programmer worth their salt immediately recognizes the possible issue of file locks being on, which in the worst case promptly raise a fatal error in the code. Luckily, different wavelengths are very likely to take a different amount of time, so when the stars align (not always, but *most* of the time), nothing should fail.

### “What to do with a mesh and a precalculated _T_-matrix?”

**Mainly, scattering dynamics calculations**

For this, one is strongly recommended to look at the `main.f90`file, particularly the `integrate()` function, and compare the contents with leftover parameters in `params.in`. 

The default function `integrate()`will perform a dynamical integration with details defined by flag `int_mode`. Funnily enough, the flag values 0 and 1 correspond to ideas toyed around with an unpublished article, and the flag value 2 corresponds to an explicit dynamical integration described in [Herranen (2017)](https://dx.doi.org/10.1002/2017RS006333).  As hinted in `params.in`, flag value 2 is useful when considering optical tweezers ([Herranen, 2019](https://doi.org/10.1371/journal.pone.0225773)), for which there are Laguerre-Gaussian and Bessel beams available. With `run_test`=3​, even the actual field expansions used within the code can be plotted using `out/field_print.py`. The log files can be plotted nicely enough with `out/evo.py` and `out/plot.py`. Nice animations can be produced from the log files with `out/animate.py`.

**“Ok, what else?”**

In `tests()` lies the answer. The most prized property of the code is in `run_test`=2​, which provides radiative torques (RATs) from a _T_-matrix. The RATs are output either (when `call torque_efficiency()` is uncommented) in so-called scattering coordinates, in a file `Q`, modified by the command-line argument `-o`. The format is exactly as in [Lazarian & Hoang (2007)](https://dx.doi.org/10.1111/j.1365-2966.2007.11817.x).  On the other hand, RATs in alignment coordinates (in format of [Draine & Weingartner (1997)](https://doi.org/10.1086/304008)) are output in a file `F` (when `call RAT_efficiency()` is uncommented). One should note that alignment frame RATs are trivially reproducible from scattering frame RATs, provided that interpolation is used.

The final possibility is to uncomment `RAT_alignment()`, which produces a rather unflatteringly dated attempt to analyze the alignment frame RAT behavior, akin to results in early papers of Lazarian & Hoang.  Again, more elegant analysis can be done using a torque efficiency `Q` sufficiently dense in its domain to provide reliable interpolation.

### Summary

The code is by my own definition a version 1.0. Still, it is plagued with issues related to the quality of code and its structure, an artifact of the fact that during the early life of the code the original author was a typical physicist coding: a rather horrible one.

The code is, apart from the parts mentioned here explicitly, rather well tested and everything *should* work correctly. If not, a needed amount of comments or pull requests can be sent to the original repository to provide some miserable soul days of backtracking where things went wrong from the times of original publications (results were cross-checked at those times).

### Future considerations

In order of decreasing tediousness:

+ Converting the code from `Fortran` to any modern language. Not really relevant, if the code remains usable to even a single person from whom things can be asked and who is willing to dedicate time helping instead of writing a solid application with proper documentation for future developers.
+ Including a surface integral equation solution for a _T_-matrix to significantly speed up calculations for scatterers of uniform composition.
+ Parallelizing and/or restructuring the logic of _T_-matrix calculations, or even go as far as to tweak the `JVIE` solution.
+ Removing any and all redundant code and provide post-processing tools that make up for any loss of features.

