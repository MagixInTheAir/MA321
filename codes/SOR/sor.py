from numpy import tril, dot, diag, array, allclose, sqrt, pi
from numpy.linalg import norm, inv

def SOR(A, B, omega): 
    M = (1/omega)*diag(diag(A))
    N = M-A
    J = dot(inv(M), N)
    K = dot(inv(M), B)
    return J,K

 

def test_SOR():
    A = array([[1, 2, -2], [1, 1, 1], [2, 2, 1]])
    B = array([[-1], [6], [9]])
    X0 = array([[0], [0], [0]])
    X = array([[1], [2], [3]])
    X_SOR, iters = res_iter(SOR, A, B, X0, 1e-6, 1e6, 1/sqrt(pi))
    assert allclose(X_SOR, X)
    return

def res_iter(decomp, A, B, X0, epsilon, itermax, omega=1):
    (J,K) = decomp(A,B,omega)
    Xprec = 0 
    X = X0
    iters = 0

    test1 = True
    test2 = True

    while test1 and iters < itermax:
        Xprec = X
        X = J@Xprec+K
        iters = iters+1

        test1 = ( norm((A@X)-B) > epsilon )
        test2 = ( norm(X-Xprec) > epsilon )
    return X,iters


test_SOR()