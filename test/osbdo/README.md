# osbdo

The multicommodity flow instances of

> T. Parshakova, F. Zhang, S. Boyd, "Implementation of an Oracle-Structured
> Bundle Method for Distributed Optimization", *Optimization and
> Engineering*, 2023

and of its code OSBDO (<https://github.com/cvxgrp/OSBDO>), solved by
`BundleSolver` in the same resource-directive form, so that the two methods
can be compared on the same problem.

## The model

The master variables are the capacities `x_i` assigned to each commodity `i`
on the edges, with `0 <= x_ie <= cap_e`. Each commodity is an agent whose
function is the value of its single-commodity flow problem,

    f_i( x_i ) = min { - b_i t : A z + t F_i = 0 , 0 <= z <= x_i , t >= 0 } ,

a `BendersBFunction` whose rows `z_e <= x_ie` take their right-hand side from
the master. The coupling `sum_i x_i <= cap` is the indicator

    g( x ) = max { mu' ( sum_i x_i - cap ) : mu >= 0 } ,

a `LagBFunction` whose inner Block has only the `mu` and a linear objective,
which `BundleSolver` treats as an "easy" component, i.e., exactly in the
master problem, as the paper does with its structured function. The same
instance is also solved as one LP, whose value is the reference.

With the optional argument `scale` equal to 1 the master variables are
`x_ie = cap_e y_ie` with `y_ie` in `[ 0 , 1 ]`: this is the diagonal
preconditioning `D = diag( u - l )` of the paper, done in the model (the
mapping of each `BendersBFunction` and the dual pairs of the `LagBFunction`
get the coefficient `cap_e` instead of 1).

## Files

- `test_osbdo_mcf.cpp`, the tester:

      osbdo_mcf_test instance BSPar LPPar LPPar_inner [ scale [ tol ] ]

  it prints the value and time of the LP and of `BundleSolver`, and fails if
  the two values differ by more than `tol` relative (default `1e-6`).
- `gen_mcf.py`, which writes an instance; its generator is the one of OSBDO,
  so with the default arguments (100 nodes, 1000 edges, 10 commodities, seed
  1001) it writes the instance of the paper, `mcf-100-1000-10.txt`. The format
  is `V E M`, then one line `tail head cap` per edge and one line
  `src dst b` per commodity.
- `mcf-20-100-3.txt`, a small instance from the same generator, which is the
  one `ctest` runs.
- `run_osbdo_mcf.py`, which runs OSBDO on an instance, after checking that the
  instance its own generator builds is the one in the file:

      python run_osbdo_mcf.py instance [ max_iter [ rel_gap [ abs_gap ] ] ]

- `BSPar.txt` (`BundleSolver`, default parameters), `MPBCfg.txt` (the Solver
  of its master problem), `LPPar_inner.txt` (the subproblems) and `LPPar.txt`
  (the reference LP).

The current `main` of OSBDO does not import: `osbdo/__init__.py` still
imports `osbdo/quasi_newton.py`, which a later commit removed. Restoring it
from the commit that added it, `git show 007b65c:osbdo/quasi_newton.py >
osbdo/quasi_newton.py`, is enough, since the bundle method does not use it.
