
#include <tuple>
#include "matrix.h"

template<typename T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> jacobi(Matrix<T> const& A, Matrix<T> const& b, long double precision) {
	Matrix<T> m(Matrix<T>::gen_diag(A.diag()));
	auto n(m - A);
	auto m_inv(m.inv());
    auto j(m_inv.dot(n));
    auto k(m_inv.dot(b));

	return std::make_tuple(m, j, k);
}