# BundleSolver

This module provides `BundleSolver`, an SMS++ `CDASolver` for the minimization
(or maximization) of a sum of convex (respectively, concave) `C05Function`,
possibly over the polyhedral set described by the `Constraint` of the `Block`.
These are known through an oracle only, i.e., one that returns their value and
a linearization at a given point. It implements the "generalized bundle"
method of [1], extended to sum-functions with "easy" components [2] and to
inexact oracles with upper models [3], in the framework of [4]. In SMS++ it is
typically used to solve the Lagrangian dual of a `Block`, whose components are
the `LagBFunction` built by `LagrangianDualSolver`, and as the master of a
Benders' decomposition, whose components are `BendersBFunction`; however, any
`C05Function` can be minimized.

Its master problem is a `MasterProblemBlock` (in the SMS++ core), which is
solved by the `Solver` registered to it, usually one of the `:MILPSolver` of
the [MILPSolver](https://gitlab.com/smspp/milpsolver) module. This allows for
four stabilizations: (i) the proximal one, whose master is a convex quadratic
program, (ii) the level one, (iii) the doubly stabilized one of [5], which
combines the two, and (iv) the trust region in the infinity norm, whose master
is a linear program. The master can be written in primal or in dual form. Each
component can be kept in it as a cutting-plane model (a "hard" component) or,
when a compact convex description of it is known (say, a `LagBFunction` whose
sub-`Block` is a linear program), it can be copied into the master as it is
(an "easy" component); since the master then has the exact description of that
component, this may reduce the number of iterations by orders of magnitude
[2].

Several devices adapt the method to the application, mostly through the
parameters of the `ComputeConfig`. They include the long-term and the
heuristic strategies that manage the proximal parameter t, and the incremental
evaluation of the components, which takes a null step as soon as a part of
them has shown it (unless the remaining evaluations are cheap with respect to
a master problem). One can also aggregate the components into fewer ones
before they enter the master problem, pass targets and accuracy to inexact
oracles, and share the global pools of linearizations with the other `Solver`
of the same functions. Furthermore, the integer variables can be kept integer,
which turns the method into a stabilized Benders' approach along the lines of
[6]. Both the bundle and the stability center are kept across calls of
`compute()`, and the `Modification` of the `Block` (e.g., new costs or
constraints of the sub-`Block`, or variables added and removed) are reflected
in them; hence, a sequence of related problems is solved by re-optimization.

Also, the lower bound returned by `get_lb()` (`get_ub()` for a maximization)
is certified. When the variables with a nonzero aggregated subgradient have
finite bounds, the master problem of each iteration gives one through the
Hölder inequality, and it also gives a stopping test that does not depend on
the threshold below which that subgradient is taken as null. Otherwise a bound
is only available when the subgradient vanishes, and it is thus only as
reliable as the threshold is small (the level and the doubly stabilized
masters provide their bounds differently). Note that `ParallelBundleSolver` is
a `BundleSolver` that evaluates the components in parallel, with a simple
master/slave scheme, and that it is otherwise the same.

The module also provides `BundleSolverML`, a `BundleSolver` whose proximal
parameter t is predicted by a small neural network (implemented with the Torch
C++ API). The network can be trained online across solves, its weights can be
shared among several instances, and they can be saved to and loaded from
TorchScript archives that are compatible with Python. `BundleSolverML` is only
built if Torch is available [see Requirements].

### References

1. [A. Frangioni "Generalized Bundle Methods" SIAM Journal on Optimization
   13(1), 117–156, 2002](http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02)

2. [A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
   with "Easy" Components: Applications to Multicommodity Network Design"
   Mathematical Programming 145(1), 133–161,
   2014](http://pages.di.unipi.it/frangio/abstracts.html#MP11c)

3. [W. van Ackooij, A. Frangioni "Incremental Bundle Methods Using Upper
   Models" SIAM Journal on Optimization 28(1), 379–410,
   2018](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT16)

4. [A. Frangioni "Standard Bundle Methods: Untrusted Models and Duality" in
   Numerical Nonsmooth Optimization: State of the Art Algorithms, A.M.
   Bagirov, M. Gaudioso, N. Karmitsa, M. Mäkelä, S. Taheri (Eds.), 61–116,
   Springer, 2020](http://pages.di.unipi.it/frangio/abstracts.html#NDOB18)

5. W. de Oliveira, M. Solodov "A Doubly Stabilized Bundle Method for
   Nonsmooth Convex Optimization" Mathematical Programming 156(1),
   125–159, 2016

6. W. van Ackooij, A. Frangioni, W. de Oliveira "Inexact Stabilized
   Benders' Decomposition Approaches with Application to Chance-Constrained
   Problems with Finite Support" Computational Optimization and Applications
   65(3), 637–669, 2016


## Getting started

These instructions will let you build the `BundleSolver` module on your system.

The module also comes ready-made, in any of

```sh
sudo add-apt-repository ppa:smspp-project/smspp   # Ubuntu
sudo apt install libsmspp-bundle-dev

conda install -c conda-forge smspp-project        # Linux, macOS, Windows

brew tap SMSpp-Project/smspp                      # macOS, Linux
brew install smspp

vcpkg install "smspp[core,bundle]"                # from the sources
```

where apt and the port give the module alone, while conda and the tap carry
the whole framework. What follows is about building it yourself.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its requirements.

- The [MILPSolver](https://gitlab.com/smspp/milpsolver) SMS++ module.

- Optionally, [Torch](https://pytorch.org/get-started/locally/) (the
  PyTorch C++ API), which is required by (and only by) BundleSolverML: with
  the makefiles it is searched for at `$(Torch_ROOT)` (see
  `extlib/makefile-default-paths-*`), with CMake via `find_package(Torch)`;
  if it is not found, BundleSolverML is simply not built.

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```
The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(BundleSolver)
target_link_libraries(<my_target> SMS++::BundleSolver)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `BundleSolver` or `ParallelBundleSolver`,
such as the [tester for `PolyhedralFunction` comparing `BundleSolver` with
`MILPSolver`](https://gitlab.com/smspp/tests/-/tree/develop/PolyhedralFunction?ref_type=heads),
has to include a "main makefile" of the module, which typically is either
[makefile-c](makefile-c) including all necessary libraries comprised the
"core SMS++" one, or [makefile-s](makefile-s) including all necessary
libraries but not the "core SMS++" one (for the common case in which this is
used together with other modules that already include them). These in turn
recursively include all the required other makefiles, hence one should only
need to edit the "main makefile" for compilation type (C++ compiler and its
options) and it all should be good to go. In case some of the external
libraries are not at their default location, it should only be necessary to
create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.

## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/bundlesolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Calandrini**  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Gorgone**  
  Dipartimento di Matematica ed Informatica  
  Università di Cagliari

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa

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
