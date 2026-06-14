# test

A self-contained internal test which provides validation for
`BundleSolverML`, the `BundleSolver` variant whose step-size t is predicted
by a neural network (implemented with the Torch C++ API) instead of the
classical rule-based heuristics, trained online across successive solves.

This executable constructs a "random" convex `PolyhedralFunction` and puts
it as the only `Objective` of an otherwise "empty" `AbstractBlock`, save for
the `ColVariable` active in the `PolyhedralFunction`.

The `AbstractBlock` is first solved by a standard `BundleSolver`, whose
optimal value is taken as the reference. Then a `BundleSolverML` is attached
and the `AbstractBlock` is repeatedly solved with it; configured with
`intMLTrainOnline`, the solver trains its network online by itself at the
end of each solve, so no explicit `Backward()` is needed. At each epoch the
optimal value is compared with the reference one and the network parameters
are checked to change (which also certifies that the network actually drove
the step-size during the solve). The tester also checks that the
`SaveModel()` / `LoadModel()` round-trip exactly restores the weights, and
that the shared-network mechanism (`set_shared_net()`, `get_shared_net()`,
`clear_shared_net()`) correctly redirects the active network among multiple
`BundleSolverML` objects.

The usage of the executable is the following:

       ./BundleSolverML_test -S <plain-BSC> -L <ML-BSC> [options]
        -S <file>   plain BundleSolver BlockSolverConfig (mandatory)
        -L <file>   BundleSolverML BlockSolverConfig (mandatory)
        -e <n>      pseudo-random generator seed [0]
        -N <n>      number of variables [10]
        -d <x>      rows / variables [4]
        -E <n>      training epochs [5]

A batch file is provided that runs a small set of tests with different sizes
and seeds of the random generator; all these passing is a good sign that no
regressions have been introduced in `BundleSolverML`.

A makefile is also provided that builds the executable including the
`BundleSolver` module and all its dependencies, in particular `MILPSolver`
and the core SMS++ library. Note that `BundleSolverML` is only compiled into
`BundleSolver` if Torch is available, see the BundleSolver README, hence this
test requires Torch to be installed at `$(Torch_ROOT)` (makefile builds) or
findable by `find_package(Torch)` (CMake builds).


## Authors

- **Francesco Demelas**  
  Laboratoire d'Informatique de Paris Nord  
  Université Sorbonne Paris Nord

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa

## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](LICENSE) file for details.
