# Writes the intersection-of-convex-sets instance of OSBDO (Parshakova,
# Zhang, Boyd 2023), generated exactly as in
# examples/intersection_cvx_sets/intersection_cvx_sets.ipynb (seed 1001, 200
# rows, 300 columns, 20 agents), in a plain text file:
#   R C M
#   R lines: the C entries of a row of A, then its b
# agent i owns the rows i * ( R / M ) , ... , ( i + 1 ) * ( R / M ) - 1.
# The generator below is a verbatim copy of osbdo.create_params.ics_params.
import sys, random
import numpy as np

def ics_params(num_row, num_col, num_agents):
    size = num_col
    A = np.random.random((num_row,num_col))-0.5
    x_feas = np.random.random(num_col)-0.5
    b = A@x_feas + 0.1*np.random.random(num_row)
    return A, b

if __name__ == "__main__":
    R, C, M, seed = 200, 300, 20, 1001
    if len(sys.argv) > 1: R, C, M, seed = map(int, sys.argv[1:5])
    np.random.seed(seed); random.seed(seed)
    A, b = ics_params(num_row=R, num_col=C, num_agents=M)
    print(R, C, M)
    for r in range(R):
        print(" ".join(repr(float(a)) for a in A[r]), repr(float(b[r])))
