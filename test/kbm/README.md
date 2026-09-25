# kbm

The large-scale nonsmooth test problems of

> N. Karmitsa, A. Bagirov, M.M. Makela, "Comparing different nonsmooth
> minimization methods and software", *Optimization Methods and Software*
> 27(1):131-153, 2012

solved by `BundleSolver`. The functions are those of `tnsunc.f`, the file of
test problems distributed with the limited memory bundle method LMBM of
Karmitsa (<https://napsu.karmitsa.fi/lmbm/>); the tester is linked with that
very file and takes the starting point from its `STARTX` and the value and a
subgradient from its `FUNC`, so that `BundleSolver` solves the same functions
from the same points as the Fortran codes. The file is not part of this
repository: the tester is built if the CMake variable `KBM_TNSUNC` gives its
path, which requires a Fortran compiler.

The function is a single component of an unconstrained Block, an
`OracleFunction`, i.e., a `C05Function` whose oracle is a callback, with a
global pool.

    kbm_test problem n BSPar [ fstar [ tol ] ]

prints the value, the iterations, the evaluations and the time, and fails if
the status is not optimal or, if `fstar` is given, if the error
`( f - fstar ) / ( 1 + | fstar | )` is larger than `tol` (default `1e-4`, the
criterion of the paper). Problems 1 to 5 are the convex ones, with optimal
values `0`, `0`, `- sqrt( 2 ) ( n - 1 )`, `2 ( n - 1 )` and `2 ( n - 1 )`.

## Reproducing the comparison

The one parameter that changes the outcome on these functions is the scaling
of the master problem, so the comparison is run for each of the four values of
`intMPHScaling` (`0` none, `1` local row scaling, which is the default, `2`
global epigraph scaling, `3` both) on each of problems 1 to 5, with the
relative accuracy `1e-6` and `intMaxIter 100000`, and no other parameter
touched. Since the oracle is an analytic routine while each iteration solves a
quadratic program of the same size, the seconds are dominated by the master
problem; `intLogVerb 1` adds one summary line per call whose last field is the
time spent in the oracle, which on these runs is between 1 and 12 % of the
total.

The times are only comparable if the machine is free, so the load is checked
first and every arm is run twice, the spread of the two rounds being the noise
to quote. The two Fortran codes are driven differently: the LMBM tester loops
over its own problems 1 to 10 and writes the times to its standard output,
while MPBNGC takes the problem number as its only argument.
