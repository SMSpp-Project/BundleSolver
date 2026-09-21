# Runs OSBDO on its multicommodity instance exactly as in
# examples/multicommodity_flow/multicommodity_flow.ipynb, after checking that
# the instance is the one in the given text file (written by gen_mcf.py).
# Usage: python run_osbdo_mcf.py instance.txt [max_iter] [rel_gap] [abs_gap]
import sys, time, random
import numpy as np
import osbdo as ob

fn = sys.argv[1]
max_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 40
rel_gap = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5
abs_gap = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-5

with open(fn) as f:
    V, E, M = map(int, f.readline().split())
    edges = [f.readline().split() for _ in range(E)]
    comms = [f.readline().split() for _ in range(M)]

np.random.seed(1001); random.seed(1001)
params = ob.mcf_params(num_vertices=V, num_edges=E, M=M)
A = params[0]['incidence']
for e, (t, h, c) in enumerate(edges):
    assert A[int(t), e] == 1 and A[int(h), e] == -1
    assert abs(params[0]['upper_bound'][e] - float(c)) <= 1e-15 * max(1, abs(float(c)))
for i, (s, d, b) in enumerate(comms):
    assert params[i]['f'][int(s)] == 1 and params[i]['f'][int(d)] == -1
    assert abs(params[i]['b'] - float(b)) <= 1e-15
print("instance check: identical to", fn, flush=True)

agents = ob.mcf_agents(params)
g = ob.mcf_coupling(agents, params)
prob = ob.Problem(agents=agents, g=g)
t0 = time.perf_counter()
prob.solve(rel_gap=rel_gap, abs_gap=abs_gap, max_iter=max_iter, print_freq=1)
tt = time.perf_counter() - t0
L, U = prob.lower_bnd[-1], prob.upper_bnd[-1]
print("OSBDO: iterations", len(prob.upper_bnd), "time", tt, "s")
print("OSBDO: lower", L, "upper", U)
for k, (l, u) in enumerate(zip(prob.lower_bnd, prob.upper_bnd)):
    print("it", k, "L", l, "U", u)
