**DRAFT VERSION**

# Colony Software

## Authors

## Build

Compile C++ project with QtCreator. The different libraries needed for compilation are listed in COLONY.pro project file at the line beginning with `LIBS +=` and also some headers at line `INCLUDEPATH +=`. QtCreator is needed in order to compile, with the following libraries:

+ [Blitz++](https://sourceforge.net/projects/blitz/) version 0.10, a fast implementation of arrays in C++. This one is used throughout the code to define arrays. Note that there is an important difference in the array reading method (as reading from the input file) between version 0.9 and 0.10.
+ [Boost](http://www.boost.org) system and Boost filesystem, version 1.48. Used for handling files for output.
+ [GSL library](https://www.gnu.org/software/gsl/), version 1.15, the GNU scientific library. Used in the random number generator.
+ [Open Dynamics Engine (ODE) library](http://www.ode.org/), version 0.11. Used in the computation of the spatial movement of the cells in the growing colony. Note that even if the GUI option is not set, the spatial dynamics of the cells can still be computed.
+ libstdc++, additional runtime library for C++ compiled with the GNU compiler, version 4.6.3.
+ For the Graphical User Interface (GUI):
  + [Qt libraries](https://www.qt.io/), >=4.8.6
+ For the Chemical Langevin integration algorithm:
  + [ATLAS library](http://math-atlas.sourceforge.net/atlas_install/), an optimized BLAS library from ATLAS (Automatically Tuned Linear Algebra Software) which is architecture-specific. Used in the matrix-vector multiplication in the Langevin algorithm. Speed up 1000x this computation step. This library has to be compiled for each architecture/machine, and result in optimized ATLAS, CBLAS and LAPACK libraries. GFORTRAN and f77BLAS seem to be needed also.
  + Alternatively, OpenBLAS. Easier to install (no need to compile for specific architecture). There is an ubuntu package libopenblas-base and libopenblas-dev.

Important options can be set in the `compilation_options.h` file that will result in different code sections being compiled and thus in different executables. We use this hard-coded options mainly for maintaining a high level of computational performance, especially to avoid conditionals that would appear in the most-inner loop of the simulation. The most important option is to choose between fixed cell volume simulation and growing/dividing cell simulation, reflected in the option `TIME_DEPENDENT_PROPENSITIES`. Note that the `compilation_options.h` file is written by the Mathematica notebook (see section "Write C++ compilation options") such that options should be set in the Mathematica notebook, otherwise changes to the header file will be overwritten when running the notebook.

------------------------

## Licensing

## Remarks

- A lot of documentation and comments are already included in the mathematica notebook.

## List of future improvements

- Implement XML input file support (much more flexible than hard-coded by line text file reading).