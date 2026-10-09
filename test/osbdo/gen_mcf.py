# Writes the multicommodity instance of OSBDO (Parshakova, Zhang, Boyd 2023),
# generated exactly as in examples/multicommodity_flow/multicommodity_flow.ipynb
# (seed 1001, 100 vertices, 1000 edges, 10 commodities), in a plain text file:
#   V E M
#   E lines: tail head cap       (A[tail,e] = +1, A[head,e] = -1)
#   M lines: src dst b           (F[i,src] = +1, F[i,dst] = -1)
# The generator below is a verbatim copy of osbdo.create_params.mcf_params.
import sys, random
import numpy as np

def mcf_params(*, num_vertices, num_edges, M):
    n = num_vertices
    A = np.zeros((num_vertices, num_edges))
    unsampled_edges = list(range(n*(n-1)))
    cur_edge = 0; global_idx = -1
    for i in range(num_vertices):
        for j in range(num_vertices):
            if i == j: continue
            global_idx += 1
            idx = (i+1) % (num_vertices)
            if j != idx: continue
            A[i, cur_edge] =  1
            A[j, cur_edge] = -1
            unsampled_edges.remove(global_idx)
            cur_edge += 1
    assert (len(unsampled_edges) == n*(n-1) - num_vertices) and (cur_edge == num_vertices)
    chosen_edges = random.sample(unsampled_edges, num_edges - num_vertices)
    chosen_pairs = random.sample(list(range(n*(n-1))), M)
    idx = 0; cur_edge = num_vertices
    F = np.zeros((M, num_vertices))
    cur_pair = 0
    for i in range(num_vertices):
        for j in range(num_vertices):
            if i == j: continue
            if idx in chosen_edges:
                A[i, cur_edge] =  1
                A[j, cur_edge] = -1
                cur_edge += 1
            if idx in chosen_pairs:
                F[cur_pair, i] =  1
                F[cur_pair, j] = -1
                cur_pair += 1
            idx += 1
    assert idx == (n*(n-1)) and (cur_edge  == len(chosen_edges) + num_vertices)
    assert cur_pair == M and (A.sum(axis=0) == 0).all()
    cap_min = 0.2; cap_max = 2
    cap = np.random.uniform(low=cap_min, high=cap_max, size=(num_edges, 1))
    b_min =  0.5; b_max = 1.5
    b = [np.random.uniform(low=b_min, high=b_max) for i in range(M)]
    return A, F, cap.squeeze(), b

if __name__ == "__main__":
    V, E, M, seed = 100, 1000, 10, 1001
    if len(sys.argv) > 1: V, E, M, seed = map(int, sys.argv[1:5])
    np.random.seed(seed); random.seed(seed)
    A, F, cap, b = mcf_params(num_vertices=V, num_edges=E, M=M)
    out = sys.stdout
    print(V, E, M, file=out)
    for e in range(E):
        t = int(np.where(A[:, e] == 1)[0][0]); h = int(np.where(A[:, e] == -1)[0][0])
        print(t, h, repr(float(cap[e])), file=out)
    for i in range(M):
        s = int(np.where(F[i] == 1)[0][0]); d = int(np.where(F[i] == -1)[0][0])
        print(s, d, repr(float(b[i])), file=out)
