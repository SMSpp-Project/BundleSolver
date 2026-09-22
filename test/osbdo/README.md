# osbdo

The examples of

> T. Parshakova, F. Zhang, S. Boyd, "Implementation of an Oracle-Structured
> Bundle Method for Distributed Optimization", *Optimization and
> Engineering*, 2023

and of its code OSBDO (<https://github.com/cvxgrp/OSBDO>), solved by
`BundleSolver` in the same form, so that the two methods can be compared on
the same problems. Each example has a generator that is a verbatim copy of the
one of OSBDO, with the same seed, so that the instances are the same, and a
tester that solves the instance with `BundleSolver` and checks the value
against a reference.

## The examples

- **Multicommodity flow** (`test_osbdo_mcf.cpp`, `gen_mcf.py`). The master
  variables are the capacities `x_i` assigned to each commodity `i` on the
  edges, with `0 <= x_ie <= cap_e`. Each commodity is an agent whose function
  is the value of its single-commodity flow problem,

      f_i( x_i ) = min { - b_i t : A z + t F_i = 0 , 0 <= z <= x_i , t >= 0 } ,

  a `BendersBFunction` whose rows `z_e <= x_ie` take their right-hand side
  from the master. The coupling `sum_i x_i <= cap` is the indicator

      g( x ) = max { mu' ( sum_i x_i - cap ) : mu >= 0 } ,

  a `LagBFunction` whose inner Block has only the `mu` and a linear objective,
  which `BundleSolver` treats as an "easy" component, i.e., exactly in the
  master problem, as the paper does with its structured function. The
  reference is the whole problem as one LP. With the optional argument
  `scale` equal to 1 the master variables are `x_ie = cap_e y_ie` with `y_ie`
  in `[ 0 , 1 ]`, i.e., the diagonal preconditioning `D = diag( u - l )` of
  the paper done in the model. With `nseq > 0` the instance is then re-solved
  `nseq` more times, each time with the capacities of the coupling drawn in
  `[ ( 1 - spread ) cap , cap ]`: the changes reach `BundleSolver` as
  Modification of the `LagBFunction`, and with the `intRstAlg` of the
  configuration one chooses between re-optimization and solving every problem
  from scratch [see `intRstAlg` in `BundleSolver.h`].

- **Supply chain** (`test_osbdo_sc.cpp`, `gen_sc.py`). A chain of agents, each
  the value of a quadratic transportation problem in its public variables
  (the flows at its input and output nodes), a `BendersBFunction` whose inner
  problem is a QP; the coupling is the linear cost of the sale and retail
  prices plus the equalities between consecutive agents, the latter an easy
  `LagBFunction` with free multipliers. The reference is the whole problem as
  one QP.

- **Intersection of convex sets** (`test_osbdo_ics.cpp`, `gen_ics.py`). Each
  agent is the squared distance from a polyhedron, a `BendersBFunction` whose
  inner problem is the QP `min { || d ||^2 : d + y = x , A_i y <= b_i }`; all
  the agents share the same master variables, which is the consensus coupling
  of the paper written without copies. The reference is the whole problem as
  one QP, whose value is 0.

- **Federated learning** (`test_osbdo_fl.cpp`, `gen_fl.py`). Each agent is the
  logistic loss of its samples, a `LogisticFunction` (a `C05Function` with the
  closed form of value and gradient and a global pool) on shared master
  variables; the regularization `lambda || x ||_1` is an easy `LagBFunction`
  with multipliers in `[ -1 , 1 ]`. There is no LP to compare with:
  `ref_fl.py` computes the reference with CVXPY, and the tester takes it as an
  optional argument.

The fourth example of the paper, the allocation of resources, is not here:
its agents are geometric means, which no `:MILPSolver` represents.

## Usage

    osbdo_mcf_test instance BSPar LPPar LPPar_inner [ scale [ tol [ nseq [ spread [ seed [ cold ] ] ] ] ] ]
    osbdo_sc_test  instance BSPar QPPar QPPar_inner [ tol ]
    osbdo_ics_test instance BSPar QPPar QPPar_inner [ tol ]
    osbdo_fl_test  instance BSPar [ ref [ tol ] ]

Every tester fails if the status of `BundleSolver` is not optimal or if its
value differs from the reference by more than `tol` (default `1e-6`),
relative to `max( 1 , | ref | )`.

- `gen_*.py` write the instances; with the default arguments they write the
  instances of the paper (of the notebooks of OSBDO for the intersection of
  convex sets), and the format is described at the top of each of them.
- `run_osbdo_mcf.py` and `run_osbdo.py` (supply chain, intersection of convex
  sets, federated learning) run OSBDO on an instance as in the notebooks of
  its examples, after checking that the instance its own generator builds is
  the one in the file, and print the bounds at every iteration with the time
  elapsed.
- `mcf-20-100-3.txt`, `sc-5.txt`, `ics-20-30-4.txt` and `fl-400-20-4.txt`
  are the instances `ctest` runs; `mcf-100-1000-10.txt` is the multicommodity
  instance of the paper, and the others are written by the generators.
- `BSPar.txt` (`BundleSolver`, default parameters), `MPBCfg.txt` (the Solver
  of its master problem), `LPPar_inner.txt` (the subproblems) and `LPPar.txt`
  (the reference LP or QP).

The current `main` of OSBDO does not import: `osbdo/__init__.py` still
imports `osbdo/quasi_newton.py`, which a later commit removed. Restoring it
from the commit that added it, `git show 007b65c:osbdo/quasi_newton.py >
osbdo/quasi_newton.py`, is enough, since the bundle method does not use it.
With numpy 2 OSBDO also calls `np.product`, which no longer exists:
`run_osbdo.py` aliases it to `np.prod`.
