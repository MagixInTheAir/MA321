#include <tuple>
#include "matrix.h"

template<typename T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> gaussSeidel(Matrix<T> const& A, Matrix<T> const& b, long double precision) {
	Matrix<T> m(A.tri_lo());
	auto n(m - A);
	auto m_inv(m.inv());
	auto g(m_inv.dot(n));
	auto k(m_inv.dot(b));
	return std::make_tuple(m, g, k);

}