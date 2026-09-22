# The reference value of a federated learning instance written by gen_fl.py:
# the whole problem min sum_s log( 1 + exp( - y_s a_s x ) ) + lambda ||x||_1
# solved by CVXPY with an interior point conic solver at high accuracy.
# Usage: python ref_fl.py instance.txt [ solver ]
import sys, time
import numpy as np
import cvxpy as cp

fn = sys.argv[1]
solver = sys.argv[2] if len(sys.argv) > 2 else "CLARABEL"
with open(fn) as f:
    n, p, M, lam = f.readline().split()
    n, p, M, lam = int(n), int(p), int(M), float(lam)
    D = np.loadtxt(f)
A, y = D[:, :p], D[:, p]
x = cp.Variable(p)
obj = cp.sum(cp.logistic(cp.multiply(-y, A @ x))) + lam * cp.norm(x, 1)
prob = cp.Problem(cp.Minimize(obj))
t0 = time.perf_counter()
if solver == "CLARABEL":
    prob.solve(solver=solver, tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)
else:
    prob.solve(solver=solver)
print("reference:", solver, prob.status, repr(prob.value), "time", time.perf_counter() - t0, "s")
