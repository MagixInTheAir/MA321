#include "matrix.h"

#include <iostream>

template<class T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> gmres(Matrix<T> const& A, Matrix<T> const& b, long double x0) {
	Matrix<T> r0(b - A.dot(x0));
	std::vector<Matrix<T>> v;
	v.push_back(r0 / r0.norm());
	Matrix<T> h(Matrix<T>::gen_full(A.lines() + 1, A.lines()));
	for (unsigned k = 0; k < A.lines(); k++) {
		Matrix<T> w(A.dot(v[k]));
		for (unsigned i = 0; i < k + 1; k++) {
			h[i][k] = w.transp().dot(v[i]);
			w -= h[i][k] * v[i];
		}
		h[k + 1][k] = w.norm();
		if (h[k + 1][k] == 0) {
			std::cout << true << std::endl;
			auto m(k);
			break;
		}
		if (k != A.lines() - 1) {
			v.push_back(w / h[k + 1][k]);
		}
	}
	auto be1(Matrix<T>::gen_full(A.lines() + 1));
	be1[0] = r0.norm();
	auto ym = Matrix<T>::leastSquares(h, be1)[0];
	return x0 + ym.dot(v);
}