#include <tuple>
#include "matrix.h"

template<typename T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> gaussSeidel(Matrix<T> const& A) {
	if (A.lines() != A.cols()) {
		throw std::logic_error("Matrix must be square for Gauss-Seidel");
	}

	Matrix<T> m(A.tri_lo(true));
	Matrix<T> n(m - A);
	Matrix<T> m_inv(Matrix<T>::solve_descent(m, Matrix<T>::gen_diag(A.lines(), A.cols(), T(1))));


	Matrix<T> g(m_inv.dot(n)); // iteration matrix
	return std::make_tuple(m, m_inv, g);
}