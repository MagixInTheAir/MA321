import numpy as np
import matplotlib.pyplot as plt

def Lagrange(X, Y, t):
    P = 0
    n = len(X)
    for k in range(n):
        Lk = 1
        for i in range(n):
            if i != k:
                Lk *= (t - X[i]) / (X[k] - X[i])
        P = P + (Y[k] * Lk)
    return P


n = [5, 10, 25, 50, 100, 250, 500, 1000, 5000, 10000]
y = [0.03182225, 0.05200005, 0.13144445, 0.5536226500000001, 1.12726755, 5.93487195, 25.3301558, 156.5459837, 7248.128887299999, 46940.1731024]

t = np.linspace(0, 10000, 10)
P1 = np.interp(10000, n, y)

P = []
for i in t:
    P.append(np.interp(i, n, y))

print(P)
plt.plot(t, P)
plt.show()
