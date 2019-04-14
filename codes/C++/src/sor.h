#include "matrix.h"

template<class T>
std::tuple<Matrix<T>, Matrix<T>> sor(Matrix<T> const& A, Matrix<T> const& b, long double omega, long double precision) {
	Matrix<T> M = (1 / omega) * Matrix<T>::gen_diag(A.diag());
	Matrix<T> N = M - A;
	Matrix<T> J = M.inv_LU().dot(N);
	Matrix<T> K = M.inv_LU().dot(B);

	return J, K;
}