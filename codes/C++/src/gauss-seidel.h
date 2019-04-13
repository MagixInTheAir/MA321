#include <tuple>
#include "matrix.h"

template<typename T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> gaussSeidel(Matrix<T> const& A) {
	Matrix<T> m(A.tri_lo(true));
	Matrix<T> n(m - A);
	Matrix<T> m_inv(m.inv());
	Matrix<T> g(m_inv.dot(n)); // iteration matrix
	return std::make_tuple(m, m_inv, g);
}