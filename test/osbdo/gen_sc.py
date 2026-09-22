# Writes the supply chain instance of OSBDO (Parshakova, Zhang, Boyd 2023),
# generated exactly as in examples/supply_chain/supply_chain.ipynb (seed
# 1001, ms = [20, 30, 40, 25, 35], ns = [30, 40, 25, 35, 20]), in a plain
# text file:
#   N                               the number of agents (the middle ones)
#   s_1 ... s_m0                    sale prices of the inputs of agent 0
#   r_1 ... r_nN                    retail prices of the outputs of agent N-1
#   then for each agent:
#   m n mu                          inputs, outputs, weight of the penalty
#   m + n upper bounds of the public variable (the lower ones are 0)
#   n lines of m entries: lin, then n lines: quad, then n lines: cap
# The generator below is a verbatim copy of osbdo.create_params.sc_params
# and of the mu of osbdo.create_params.sc_agents (osbdo.utils.get_mus).
import sys, random, copy
import numpy as np

def sc_params(ms, ns):
    def generate_costs(n, m):
        lin  = np.exp(np.random.normal(loc = 0.07, scale = 0.7, size = (n, m)))
        cap  = np.exp(np.random.normal(loc = 0,   scale = 1, size = (n, m)))
        quad = 0.5 * np.divide(lin, cap)
        return [lin, quad, cap]

    N = len(ms) + 2
    params = [0] * N
    sale   = np.random.uniform(low = 8,  high = 10, size = (ms[0], 1))
    retail = np.random.uniform(low = 10, high = 12, size = (ns[-1], 1))
    params[0]   = {"m":0,      "n":ms[0],  "sale":sale, "dimension" : ms[0] }
    params[N-1] = {"m":ns[-1], "n":0, "retail":retail, "dimension" : ns[-1]}

    ranges = [0] * N
    for i, (m, n) in enumerate(zip(ms, ns)):
        lin, quad, cap = generate_costs(n, m)
        idx = i + 1
        params[idx] = {"cap":cap, "lin":lin, "quad":quad, "dimension": m + n, "m":m, "n":n, "norm":"l1"}
    for i, (m, n) in enumerate(zip(ms, ns)):
        idx = i + 1
        if N-1 > idx > 0:
            a = params[idx]["cap"].sum(axis=0).T
            b = params[idx]["cap"].sum(axis=1)
            if idx != N-2:
                b = np.maximum(b, params[idx+1]["cap"].sum(axis=0).T)
            if idx != 1:
                a = np.maximum(a, params[idx-1]["cap"].sum(axis=1))
            assert a.shape[0] == m and b.shape[0] == n
            ranges[idx] = np.concatenate([a, b], axis=0)
    ranges[0] = params[1]["cap"].sum(axis=0).T
    ranges[N-1] = params[N-2]["cap"].sum(axis=1)
    for i in range(N):
        assert np.zeros((params[i]["dimension"], )).shape == ranges[i].shape
        params[i]["lower_bound"] = np.zeros((params[i]["dimension"], ))
        params[i]["upper_bound"] = copy.deepcopy(ranges[i])
    return params

if __name__ == "__main__":
    ms, ns, seed = [20, 30, 40, 25, 35], [30, 40, 25, 35, 20], 1001
    np.random.seed(seed); random.seed(seed)
    P = sc_params(ms, ns)
    f = lambda v: " ".join(repr(float(a)) for a in np.ravel(v))
    print(len(ms))
    print(f(P[0]["sale"]))
    print(f(P[-1]["retail"]))
    for p in P[1:-1]:
        print(p["m"], p["n"], 50.0)   # get_mus( 50 , ... )
        print(f(p["upper_bound"]))
        for key in ("lin", "quad", "cap"):
            for row in p[key]:
                print(f(row))
