import numpy as np
import scipy.sparse.linalg as sc


def GMRES(A, b, x0):
    r0 = b - A @ x0
    v = list()
    v.append(r0 / np.linalg.norm(r0))
    h = np.zeros((np.shape(A)[0] + 1, np.shape(A)[0]))
    for k in range(np.shape(A)[0]):
        w = A @ v[k]
        for i in range(k + 1):
            h[i, k] = w.T @ v[i]
            w = w - h[i, k] * v[i]
        h[k+1, k] = np.linalg.norm(w)
        if h[k+1, k] == 0:
            print(True)
            m = k
            break
        if k != np.shape(A)[0]-1:
            v.append(w / h[k+1][k])
    be1 = np.zeros(np.shape(A)[0] + 1)
    be1[0] = np.linalg.norm(r0)
    ym = np.linalg.lstsq(h, be1, rcond=None)[0]
    return x0 + ym @ np.asarray(v)


A = np.array([[1., 5, 4, 8, 9, 12], [2, 3, 5, 8, 4, 3], [0, 5, 3, 8, 9, 1], [2, 4, 6, 8, 1, 3], [5, 8, 6, 4, 7, 10], [0, 0, 5, 6, 8, 6]])
b = np.array([[10], [12], [8], [51], [42], [27]])
x0 = np.zeros((6, 1))
x = GMRES(A, b, x0)
result = sc.gmres(A, b, x0)
print(result[0])
print(x)
print(np.linalg.norm(A @ result[0].T - b))
print(np.linalg.norm(A @ x - b))
