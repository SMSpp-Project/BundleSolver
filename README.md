# BundleSolver

BunldeSolver class, which implements the CDASolver interface within the SMS++
framework. The class implements a NonDifferentiable Optimization Solver using
a "Generalized Bundle" algorithm; cf. e.g.

- [A. Frangioni "Generalized Bundle Methods" SIAM Journal on Optimization
  13(1), p. 117 - 156,
  2002](http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02)

- [A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
  with "Easy" Components: Applications to Multicommodity Network Design"
  Mathematical Programming 145(1), 133–161,
  2014](http://pages.di.unipi.it/frangio/abstracts.html#MP11c)

- [W. van Ackooij, A. Frangioni "Incremental Bundle Methods Using Upper
  Models" SIAM Journal on Optimization 28(1), 379–410,
  2018](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT16)

- [A. Frangioni "Standard Bundle Methods: Untrusted Models and Duality" in
  Numerical Nonsmooth Optimization: State of the Art Algorithms, A.M. Bagirov,
  M. Gaudioso, N. Karmitsa, M. Mäkelä, S. Taheri (Eds.), 61—116, Springer,
  2020](http://pages.di.unipi.it/frangio/abstracts.html#NDOB18)


## Getting started

These instructions will let you build MCFBlock and MCFSolver on your system.


### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.

- The [MILPSolver](https://gitlab.com/smspp/milpsolver) SMS++ module.

- [The NDOSolver/FiOracle
  project](https://gitlab.com/frangio68/ndosolver_fioracle_project) and its
  requirements (depending on the actual MPSolver built); note that this
  dependency is supposed to be removed down the line.


### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
sudo make install
```


### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(BundleSolver)
target_link_libraries(<my_target> SMS++::BundleSolver)
```


### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. General instructions are:

- The arrangements of folders must be that envisioned by the
  [Umbrella SMS++ Project](https://gitlab.com/smspp/smspp-project)

- The main step is to edit the makefiles into ../extlib/. There is one for
  each of the external libraries that any module requires, starting with

  = [Boost](https://www.boost.org)

  = [Eigen](http://eigen.tuxfamily.org)

  = [netCDF-C++](https://www.unidata.ucar.edu/software/netcdf)

  that are required by the "core" SMS++ library and therefore by everyone.
  Setting the

```make
lib*INC = -I<paths to include files directories>
lib*LIB = -L<paths to lib files directories> -l<libs>
```

  in each allows one to set any non-standard path if the library is not
  installed in the system (or leave them empty if they are).

- The "core" SMS++ classes have a makefile for building the corresponding
  library in

```sh
SMS++/lib/makefile-lib
```

  The makefile allow to choose the compiler name and the optimization/debug.
  This builds the lib/libSMS++.a that can be linked upon. Also, the

```sh
SMS++/lib/makefile-inc
```

  file is provided for allowing external makefiles to ensure that the library
  is up-to-date (useful in case one is actually developing it). The simplest
  way to learn how to use it is to check the makefiles of the "main" file

```sh
Main/makefile
```

  Note that the "basic" makefile macros

```make
CC =
SW =
```

  for setting the c++ compiler and its options are "automatically forwarded"
  from the makefile to these of the other SMS++ components, and therefore
  (possibly at the cost of a make clean) ensure consistency during the
  building process.

- The [The NDOSolver/FiOracle
  project](https://gitlab.com/frangio68/ndosolver_fioracle_project)
  has a similar arrangement with its own extlib/ folder that must be
  independently edited in an analogous way.

## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/bundlesolver/-/issues/new).

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Gorgone**  
  Dipartimento di Matematica ed Informatica  
  Università di Cagliari

### Contributors

## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
