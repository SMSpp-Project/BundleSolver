# Runs OSBDO on its supply chain ("sc"), intersection-of-convex-sets
# ("ics") or federated learning ("fl") instance exactly as in the timed cell
# of the notebook in examples/, after checking that the instance its own
# generator builds is the one in the given text file (written by gen_sc.py,
# gen_ics.py or gen_fl.py).
# Usage: python run_osbdo.py sc|ics|fl instance.txt [max_iter] [rel_gap] [abs_gap]
import sys, time, random
import numpy as np
import osbdo as ob

# numpy 2 has no np.product, which OSBDO still calls when it discovers rho
if not hasattr(np, "product"):
    np.product = np.prod

# every line printed from here on carries the seconds since the start
class Stamped:
    def __init__(self, out): self.out, self.t0, self.bol = out, time.perf_counter(), True
    def write(self, s):
        for part in s.splitlines(keepends=True):
            if self.bol: self.out.write("[%10.3f] " % (time.perf_counter() - self.t0))
            self.out.write(part); self.bol = part.endswith("\n")
    def flush(self): self.out.flush()
sys.stdout = Stamped(sys.stdout)

kind, fn = sys.argv[1], sys.argv[2]
max_iter = int(sys.argv[3]) if len(sys.argv) > 3 else {"sc": 150, "ics": 30, "fl": 100}[kind]
rel_gap = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-5
abs_gap = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-5

np.random.seed(1001); random.seed(1001)
nums = np.array([float(a) for a in open(fn).read().split()])
if kind == "sc":
    params = ob.sc_params([20, 30, 40, 25, 35], [30, 40, 25, 35, 20])
    agents = ob.sc_agents(params)
    mine = [params[0]["sale"].ravel(), params[-1]["retail"].ravel()]
    for p in params[1:-1]:
        mine += [np.array([p["m"], p["n"], p["mu"][0]]), p["upper_bound"],
                 p["lin"].ravel(), p["quad"].ravel(), p["cap"].ravel()]
    mine = np.concatenate([[len(params) - 2]] + mine)
    g = ob.sc_coupling(params, agents)
elif kind == "fl":
    n, p_, M = map(int, nums[:3])
    params = ob.fl_params(num_samples=n, num_agents=M, size=p_)
    X = np.concatenate([q["x"] for q in params]); Y = np.concatenate([q["y"] for q in params])
    mine = np.concatenate([[n, p_, M, 5.0], np.hstack([X, Y[:, None]]).ravel()])
    agents = ob.fl_agents(params)
    g = ob.fl_coupling(agents, params)
else:
    R, C, M = map(int, nums[:3])
    params = ob.ics_params(num_row=R, num_col=C, num_agents=M)
    A = np.concatenate([p["A"] for p in params]); b = np.concatenate([p["b"] for p in params])
    mine = np.concatenate([[R, C, M], np.hstack([A, b[:, None]]).ravel()])
    agents = ob.ics_agents(params)
    g = ob.ics_coupling(agents, params)
assert mine.shape == nums.shape and np.array_equal(mine, nums)
print("instance check: identical to", fn, flush=True)

prob = ob.Problem(agents=agents, g=g)
t0 = time.perf_counter()
prob.solve(rel_gap=rel_gap, abs_gap=abs_gap, max_iter=max_iter, print_freq=1)
tt = time.perf_counter() - t0
print("OSBDO: iterations", len(prob.upper_bnd), "time", tt, "s")
print("OSBDO: lower", prob.lower_bnd[-1], "upper", prob.upper_bnd[-1])
for k, (l, u) in enumerate(zip(prob.lower_bnd, prob.upper_bnd)):
    print("it", k, "L", l, "U", u, flush=True)
