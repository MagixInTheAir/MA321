#include "matrix.h"

template<class T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> sor(Matrix<T> const& A,long double omega) {
	if (A.lines() != A.cols()) {
		throw std::logic_error("Matrix must be square for SOR resolution");
	}

	Matrix<T> m(Matrix<T>::gen_diag(A.diag()) * (1/omega));
	Matrix<T> n(m - A);
	Matrix<T> m_inv(m.lines(), m.cols());
	for (unsigned i = 0; i < m.lines(); i++) {
		m_inv[i][i] = (1 / m[i][i]);
	}
	
	Matrix<T> j(m_inv.dot(n));

	return std::make_tuple(m, m_inv, j);
}