#include "matrix.h"
#include <algorithm>
#include <cmath>

template<typename T>
std::tuple<Matrix<T>, long long, T> iter_generale(Matrix<T> const& m_inv, Matrix<T> const& iter, Matrix<T> const& b, 
                    Matrix<T> const& x0, T const& epsilon, long long maxiter) {

    Matrix<T> xs(x0);
    Matrix<T> xp(iter.dot(xs) + m_inv.dot(b));
    long long nb_iter = 1;
	T errp(0);
	Matrix<T> err_m(xp - xs);
    T err(err_m.norm());
    
    while(err > epsilon && nb_iter < maxiter) {
        xs = xp;
        xp = iter.dot(xs) + m_inv.dot(b);
        nb_iter++;
		errp = err;
		Matrix<T> err_m(xp - xs);
        err = err_m.norm();

		// std::cout << err << std::endl;

		if (err > errp) {
			throw std::logic_error("Method is not convergent");
		}
    }


    return std::make_tuple(xp, nb_iter, err);
}