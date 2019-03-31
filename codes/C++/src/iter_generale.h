#include "matrix.h"
#include <algorithm>
#include <cmath>

template<typename T>
void iter_generale(Matrix<T> const& m, Matrix<T> const& iter, Matrix<T> const& b, 
                    T x0, T epsilon, long long maxiter) {

    rho = iter.highest_eigen();
    if(rho >= 1) {
        throw std::logic_error("La méthode ne converge pas");
    }

	Matrix<T> m_inv = m.inv();
    auto xs = x0;
    auto xp = iter.dot(xs) + m_inv.dot(b);
    long long iter = 1;
    T err = (xp - xs).norm();
    
    while(err > epsilon && iter < maxiter) {
        xs = xp;
        xo = iter.dot(xs) + m_inv.dot(b);
        iter++;
        err = (xp - xs).norm();
    }

    return std::make_tuple(xp, iter, err);
}