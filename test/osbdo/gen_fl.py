# Writes the federated learning instance of OSBDO (Parshakova, Zhang, Boyd
# 2023), generated exactly as in
# examples/federated_learning/federated_learning.ipynb (seed 1001, 10000
# samples, 10 agents, 500 features, lambda = 5 as in fl_coupling), in a plain
# text file:
#   n p M lambda
#   n lines: the p features of a sample, then its label in { -1 , 1 }
# agent i owns the samples i * ( n / M ) , ... , ( i + 1 ) * ( n / M ) - 1.
# The generator below is a verbatim copy of osbdo.create_params.fl_params
# (and of osbdo.utils.sign).
import sys, random
import numpy as np

def sign(z):
    z[z>=0] = 1
    z[z<0] = -1
    return z

def fl_params(num_samples, num_agents, size):
    assert(num_samples >= num_agents)
    cols = random.sample(range(size), int(size/10))
    theta_true = np.zeros((size,))
    theta_true[cols] = np.random.normal()
    X = np.random.normal(size=(num_samples, size))
    Y = sign(X @ theta_true + 0.1*np.random.normal())
    return X, Y

if __name__ == "__main__":
    n, p, M, seed = 10000, 500, 10, 1001
    if len(sys.argv) > 1: n, p, M, seed = map(int, sys.argv[1:5])
    np.random.seed(seed); random.seed(seed)
    X, Y = fl_params(num_samples=n, num_agents=M, size=p)
    out = sys.stdout
    print(n, p, M, 5.0, file=out)
    for s in range(n):
        print(" ".join(repr(float(a)) for a in X[s]), repr(float(Y[s])), file=out)
