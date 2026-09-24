# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `test/osbdo`, the multicommodity instances of OSBDO solved in the same
  resource-directive form, with a `BendersBFunction` per commodity and the
  coupling as an easy `LagBFunction`, against the whole problem as one LP;
  the generator is the one of OSBDO and a script runs OSBDO on the same
  instance, so that the two are compared on what they both solve

### Changed

- the noise reduction has a global memory: it sets t to dblmxIncr times the
  largest among t and the values it has set since the last serious step,
  and it stops the solver when that value is already dbltMaior, so that the
  decreases of t in a sequence of null steps can no longer undo it and let
  it happen forever; the two noise reductions, before the stopping test and
  after an evaluation that gives neither a serious nor a null step, are one
  method

- the makefile asks for `-O3 -DNDEBUG` and nothing else, the macro of the
  patch for `boost::any` on macOS having no reason to be there since there is
  no `boost::any` left in the core

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

- the master of the tests asks Gurobi for its least numerical care and not
  for none of it, and declares the residual zero on the scale of the model:
  the extra care is paid at every one of the thousands of solves of a run

- the tests of the ML variant link `SMS++::BundleSolverML`, which is where
  `BundleSolverML` now lives

- the configuration that is installed no longer looks for NDOSolver/FiOracle,
  which the library does not link any more

- the three places where the parallel loop names a failing component write
  the log at the same verbosity, and the log says which status the component
  that stopped the loop returned, a stop with no reason having left whoever
  read it to guess among the components

- `ParallelBundleSolver` keeps its threads alive and wakes up when an
  evaluation ends, rather than starting and joining a thread per component at
  every iteration, and never lets a thread wait on one that has already
  finished

- the master problem is a `MasterProblemBlock`, i.e., a Block of the model
  solved by whichever Solver is attached to it, in place of the `MPSolver`
  hierarchy of the previous versions: what the bundle asks of it is said
  through the abstract representation and the Modification, so that the
  master is built, solved and changed as any other Block, and the parameters
  that only the old hierarchy understood are refused where they no longer
  mean anything. With the MPSolver hierarchy go the NDOSolver/FiOracle
  submodule and the Osi and Clp requirements: a BundleSolver needs the core
  library and MILPSolver, and the Solver of the master is the one the
  BlockSolverConfig that `strMPBSolverCfg` points at attaches

- the Modification that `MasterProblemBlock` issues in a loop travel in one
  channel: the shift of the constant of every cut at a move of the reference
  [see `set_reference()`], the refresh of the box, the linear part of the
  0-th component and the quadratic term of every z at a change of t. A
  Solver able to write a whole set of coefficients, or of sides, in one
  operation then does that instead of one call per cut, per variable or per
  coupling row

- the groups of components take the threads they spend on their members from
  the pool of the parallel solver driving them, through
  `C05SumFunction::set_submitter()`, instead of starting one per member: a
  thread costs of the order of 100 microseconds to start and the members of
  the instances of interest cost less than that, which is why the hand-down
  was worth it only above a threshold. The two levels share one pool, sized
  for both, and a group evaluates one of its members in the thread that is
  waiting for the others anyway, so that no thread of the pool is ever held
  doing nothing and the members cannot be starved by the components

- the fifth loop of `process_outstanding_Modification()` asks whether any of
  the Variable a Modification speaks of is still "active" in the component it
  comes from, that component being reloaded from its global pool whatever has
  changed. It used to actualise the names against the first component, "which
  is fairly taken as a representative since all the C05Function have the same
  active Variable", into a range that nothing read: an assumption that no
  longer holds when the Lambda is sparse, in sixty lines whose only effect
  was to decide whether the component is reloaded at all

- BundleSolverML is a library of its own, SMS++::BundleSolverML: Torch is
  some hundreds of megabytes of shared objects, and a program linking
  BundleSolver paid the loading of every one of them at each start, 0.2 s per
  process on our machines, whether or not the ML variant was ever used.
  Whoever wants that variant links the new library, which brings BundleSolver
  along with it

### Fixed

- the combination of linearizations given back to a component when its
  bundle is full is divided by the mass its diagonal rows carry, lambda
  minus the share of its individual lower bound, rather than by 1 - r, so
  that it is a convex combination also with a level row or an individual
  lower bound; with all the mass on the lower bounds the item is replaced
  without any aggregation
- a full bundle frees an item in base that is a vertical row only when fewer
  than two diagonal ones are in base, the aggregate being made of the
  diagonal rows alone and the multiplier of the vertical one being otherwise
  lost
- `is_subgradient_global()` asks the master by the hard component and the
  global name of the item, as the other helpers do, rather than by the
  position in the pool, which with easy components made every row vertical

- the duals of an easy component, and the reduced costs of its columns, are
  the ones of the last solve of the master: the sub-Block of an easy
  component is a Block of the model, which any Solver may write into between
  that solve and the question, so what it held when asked was not what the
  master had left there, as was already the case for its primal. They are
  saved only when `intDoEasy` says that they will be asked for

- a master that the Solver could not solve, and whose bundle held no item
  to remove, was the end of the run: it is now solved once more with t as
  it was before the empty bundle brought it down to its minimum, where the
  quadratic term vanishes against the data of the easy components and
  leaves the master a badly scaled problem

- the part of an easy component in the solution of a LagrangianDualSolver
  was whatever its inner Block held when the solution was asked for: the
  Solver of the master writes there after each solve, but the inner Block
  is a Block of the model, which any other Solver may write into after it,
  so that the solution could mix the multipliers of the master for the hard
  components with the solution of another Solver for the easy ones. The
  master now saves what it writes there, and `get_dual_solution()` writes
  it back [see `MasterProblemBlock::restore_easy_primal()`]

- a variable change of one component that is not wrapped in a
  `GroupModification` (a "naked" `FunctionModVars`) threw when there were
  more components, and the throw only hid what went wrong behind it. It is
  now processed per component, and three things it needs have been put
  right: (i) with the sparse representation of the Variable the
  linearizations reloaded into the master were read as dense, which put
  their coefficients on the wrong global Variable; (ii) a component that
  stops depending on some of its Variable, while the others still do, has
  its linearizations reloaded from its global pool (a reset of that
  component only, the others being no longer reset), and its value in the
  stability centre is marked unknown if one of those Variable is not 0
  there; (iii) an easy component has its Lagrangian terms dropped from the
  coupling rows of the Variable it no longer has, through the new
  `MasterProblemBlock::drop_easy_coupling()`, and the master is given the
  new local-to-global maps of the easy components

- the component index passed to `MasterProblemBlock`: it keeps a bundle only
  for the "hard" components and indexes them by their position among those,
  while BundleSolver counts the components of Fi globally, the "easy" ones
  included. With no easy component the two numbers coincide, so the
  difference showed up only under `intDoEasy`, where it addressed the wrong
  `PolyhedralFunctionBlock` or ran past the last one; the translation is now
  made explicit by `hard_k()` and, for the easy components, by `easy_k()`

- a Modification changing the linearizations of a component without
  changing its Variable only asked the Solver of the master to read the
  subgradients again, but the master holds copies of them and so gave back
  the stale ones; the component is now reset, which reloads its
  linearizations from its global pool as they are

- when a component had used up its share of the bundle with constraints
  (vertical linearizations, which are never removed), a new one was still
  put in a free spot elsewhere in the bundle, beyond the global pool of the
  component, and the C05Function then refused to store it with "invalid
  linearization name"; this happens when a Lagrangian subproblem is unbounded
  over and over, and now the solve ends with kError and says that the bundle
  of that component is full of constraints

- a component answering kLowPrecision stopped the whole solve with an error:
  that code sorts after kError among the return codes, and the two places
  that read the status of a component only asked whether it was at least
  kError. A component saying kLowPrecision has found a solution and says it
  could not prove it optimal, which is inexact information and not a
  failure, so it is now let through and used as such

## [0.5.0] - 2026-09-12

### Added

- BundleSolverML, a BundleSolver variant whose step-size t is predicted by
  a neural network (Torch) trainable online across solves, with shared
  weights among multiple instances and TorchScript model save / load; it is
  only built if Torch is available.

### Changed

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Fixed

- the documentation of the int parameters, which did not say that
  intMaxThread is ignored, the implementation being sequential

- BPar5 and RstAlgPrm were never initialized in the BundleSolver
  constructor, leading to nondeterministic behavior (up to complete stalls
  of the algorithm) whenever dblBPar5 / intRstAlg were not explicitly set
  by the ComputeConfig.

- the package configuration file finds the libraries the module links, so that
  a project using the installed module needs nothing more than find_package()

## [0.4.5] - 2025-12-12

### Added

- check for when all components are easy

- temporary solution to management of putting an
  inconsistent BundleSolverState (see comments to
  intFrcLstSS)

- useful option to get\_var\_solution() for duals/reduced
  costs of all easy components

- support for reading primal and dual solutions out of
  easy components

- CHECK_BAD_F (warnings to show possible oracle issues)

- algo information getters

### Changed

- Heuristic() are protected and virtual (so that we can
  develop new ML-based ones)

- support to un\_any\_thing\_count\_*

- adapted to new standard organization of makefiles

- avoided const vectors prone to static initialization fiasco

- deactivated useless warnings

### Fixed

- huge logical flaw whereby a problem that was empty due to
  easy components making it so was not properly recognised
  as being so

- allow to terminate in exactly one iteration

## [0.4.4] - 2024-02-27

### Changed

- makefiles and Cmake files updated to new global SMS++ scheme

- documentation updated accordingly

- default OSISolver is Gurobi

- improved BundleSolverState by keeping all upper and lower values for each component

### Fixed

- fixed sloppy setting of OSI parameters

- fixed `has_*_solution()` and `has_*_feasible()`

## [0.4.3] - 2023-05-23

### Added

- handling the constant term Modification for hard
  components

- handling constant term for "easy" components

- support for OsiGrbSolverInterface

### Changed

- better printing in case of errors in a Fi

### Fixed

- missing check of max number of iterations in one case

- issues in put_State()

- issue in one case of process_outstanding_Modification (adding
  Linearization was not properly managed)

- management of numerical errors in the MPSolver

- typos and spaces

## [0.4.2] - 2022-06-28

### Added

- added new parameter intTrgtMng controlling different aspects of how models
  are used to pass targets and accuracy to the C05Function and to forecast
  the value (using the Lipschitz_constant)

### Changed

- generate the abstract representation in set_Block()

- minor improvements in logging and comments

- quite significant rehaul of t heuristics: reshuffled the bits in tSPar1
  (baiscally keeping compatibility with develop branch) and added the "reversal
  form of the poorman's quasi-Newton update". This requires the computation of
  the norm of the "aggregated representative subgradient", doing which called
  for some nontrivial changes in the computation of Alfa1 and ScPr1

### Fixed

- ParallelBundleSolver now properly counts the computed components; previously
  those of the ramp-down phase were not, which may lead to getting an error if
  the last step was not a SS and intFrcLstSS was true

- fixed get\_dual\_solution() when there is no easy component

## [0.4.1] - 2021-12-07

### Changed

- improved Modification handling (no over-reacting to easy ones, important bugfix)

- better log, printing times for each component in verbosity 4 and higher

- fixed several flaws

- improved namespace handling and similar stuff

## [0.4.0] - 2021-05-02

### Changed

Several major improvements:

- restructured stopping conditions, with a new one on the norm (INF, 1 or 2)
  of zStar, possibly relative to subgradient ; tStar anyway remains the driver
  of t-strategies.

- BundleSolver now detects when "zStar == 0" (see previous point) and uses
  this to algorithmically produce globally valid lower bounds.

- BundleSolver now checks if the objectives are (all) convex or concave and
  adapts accordingly (by sneakily converting the concave problem internally
  into a convex one and adapting the results on the fly when they go out).

- Box constraints are now handled, either by means of NNConstraint/BoxConstraint
  or as NN constraints "inherent" in the ColVariable.

- Changes in the costs or RHS/bounds in easy components now are allowed and managed.

- Changes in the 0-th component now allowed.

- Added hooks for perspective handling of changes in easy components.

- Complete rehaul of Modification in "easy" components.

- Better bad termination reporting.

- Allowing to switch away Cplex.

- Added ComputeConfig for hard and easy and NoEasy.

- Added individual Configurations.

### Fixed

- Many fixes throughout the code, too many to list.

## [0.3.0] - 2020-09-16

### Changed

- Support for better modification framework.

### Fixed

- Minor fixes.

## [0.2.0] - 2020-03-06

### Fixed

- Minor fixes.

## [0.1.2] - 2020-03-02

### Fixed

- Minor fixes.

## [0.1.1] - 2020-02-17

### Fixed

- Minor bug in CMake configuration file.

## [0.1.0] - 2020-02-10

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/bundlesolver/-/compare/0.5.0...develop
[0.5.0]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.5...0.5.0
[0.4.5]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.4...0.4.5
[0.4.4]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.3...0.4.4
[0.4.3]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.2...0.4.3
[0.4.2]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.1...0.4.2
[0.4.1]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/bundlesolver/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/bundlesolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/bundlesolver/-/compare/0.1.2...0.2.0
[0.1.2]: https://gitlab.com/smspp/bundlesolver/-/compare/0.1.1...0.1.2
[0.1.1]: https://gitlab.com/smspp/bundlesolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/bundlesolver/-/tags/0.1.0
