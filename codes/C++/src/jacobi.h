
#include <tuple>
#include "matrix.h"

template<typename T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> jacobi(Matrix<T> const& A) {
	if (A.lines() != A.cols()) {
		throw std::logic_error("Matrix must be square");
	}

	Matrix<T> m(Matrix<T>::gen_diag(A.diag()));
	Matrix<T> n(m - A);

	Matrix<T> m_inv(A.lines(), A.cols());
	for (unsigned i = 0; i < A.lines(); i++) {
		m_inv[i][i] = (1 / m[i][i]);
	}

    Matrix<T> j(m_inv.dot(n));

	return std::make_tuple(m, m_inv, j);
}