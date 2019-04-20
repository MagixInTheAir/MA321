import numpy as np
import scipy.sparse.linalg as sc
import matplotlib.pyplot as plt


def GMRES(A, b, espilon):
    max_iter = A.shape[0] #Number of maxiter
    mat_q = np.zeros((max_iter, max_iter + 1))
    mat_h = np.zeros((max_iter + 1, max_iter))
    norm_b = np.linalg.norm(b)
    be1 = np.zeros(max_iter + 1)
    be1[0] = norm_b
    mat_q[:, 0] = 1 / np.linalg.norm(b) * b.T # On définit ici que l'on a forcément x0 = 0
    for j in range(max_iter):
        mat_q[:, j+1] = A @ mat_q[:, j]
        for i in range(j+1):
            mat_h[i, j] = mat_q[:, i] @ mat_q[:, j + 1]
            mat_q[:, j+1] -= mat_h[i, j] * mat_q[:, i]
        mat_h[j + 1, j] = np.linalg.norm(mat_q[:, j + 1])
        mat_q[:, j + 1] /= mat_h[j + 1, j]
        y = np.linalg.lstsq(mat_h, be1, rcond=None)[0]
        print(np.linalg.norm(mat_h@y-be1)
        residue = np.linalg.norm(y) / norm_b
        if residue < espilon:
            return mat_q[:max_iter, :max_iter] @ y, residue
    return mat_q[:max_iter, :max_iter] @ y


n = 300
"""
A = np.array([[1., 5, 4, 8, 9, 12], [2, 3, 5, 8, 4, 3], [0, 5, 3, 8, 9, 1], [2, 4, 6, 8, 1, 3], [5, 8, 6, 4, 7, 10], [0, 0, 5, 6, 8, 6]])
b = np.array([[172], [93], [102], [83], [190], [115]])
x0 = np.zeros((6, 1))"""
l_r = list()
l_n = list()
for n in np.arange(5, 1000, 100):
    A = 100 * np.random.rand(n, n)
    b = 100 * np.random.rand(n, 1)
    x0 = np.zeros((n, 1))

    x = GMRES(A, b, 1e-8)
    mat_exp = A @ x
    r = np.linalg.norm(mat_exp.reshape((n, 1)) - b)
    l_r.append(r)
    l_n.append(n)
    """result = sc.gmres(A, b, x0)
    mat_the = np.array(A @ result[0])
    print(np.linalg.norm(mat_the.reshape((n, 1)) - b))"""

plt.plot(l_n, l_r)
plt.yscale('log')
plt.xlabel('Taille de la matrice')
plt.ylabel('Erreur sur le calcul de Ax - b')
plt.show()
plt.savefig('fig.png')