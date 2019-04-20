#include "matrix.h"
#include "NumCpp/NumCpp.hpp"

#include <iostream>
#include <tuple>

template<class T>
std::tuple<Matrix<T>, long long, T> gmres(Matrix<T> const& A, Matrix<T> const& b, T const& epsilon) {
	
	size_t max_iter = A.lines();
	Matrix<T> mat_q(Matrix<T>::gen_full(max_iter, max_iter + 1, T(0)));
	Matrix<T> mat_h(Matrix<T>::gen_full(max_iter + 1, max_iter, T(0)));
	Matrix<T> be1(Matrix<T>::gen_full(max_iter + 1, 1, T(0)));
	be1[0][0] = b.norm();
	mat_q.setCol(0, b * (1 / b.norm()));
	long long nbiter(0);
	T residue(0);
	Matrix<T> y;
	for (size_t j = 0; j < max_iter; j++) {
		nbiter++;
		mat_q.setCol(j + 1, A.dot(mat_q.col(j)));
		

		for (size_t i = 0; i < j + 1; i++) {
			// equivalent : mat_h[i][j] = mat_q.col(i).dot(mat_q.col(j + 1))[0][0];
			T sum(0);
			for (size_t k = 0; k < mat_q.lines(); k++) {
				sum += mat_q[i][k] * mat_q[k][j+1];
			}
			mat_h[i][j] = sum;

			// equivalent : mat_q.setCol(j+1, mat_q.col(j+1) - mat_h[i][j]*mat_q.col(i))
			for(size_t k = 0; k < mat_q.lines(); k++) {
				mat_q[k][j + 1] -= mat_h[i][j] * mat_q[k][i];
			}
		}


		mat_h[j + 1][j] = mat_q.col(j + 1).norm();
		mat_q.setCol(j + 1, mat_q.col(j + 1) / mat_h[j + 1][j]);

		nc::NdArray<T> mat_h_2(mat_h.lines(), mat_h.cols());
		for (size_t i = 0; i < mat_h.lines(); i++) {
			for (size_t j = 0; j < mat_h.cols(); j++) {
				mat_h_2(i, j) = mat_h[i][j];
			}
		}

		nc::NdArray<T> be1_2(be1.lines(), be1.cols());
		for (size_t i = 0; i < be1.lines(); i++) {
			for (size_t j = 0; j < be1.cols(); j++) {
				be1_2(i, j) = be1[i][j];
			}
		}

		auto rep_2(nc::linalg::lstsq(mat_h_2, be1_2));
		
		y = Matrix<T>(rep_2.shape().rows, rep_2.shape().cols);
		for (size_t i = 0; i < rep_2.shape().rows; i++) {
			for (size_t j = 0; j < rep_2.shape().cols; j++) {
				y[i][j] = rep_2(i, j);
			}
		}
		y = y.transp();

		residue = y.norm() / b.norm();
		if (residue < epsilon) {
			break;
		}
	}

	return std::make_tuple(mat_q.truncate(max_iter, max_iter).dot(y), nbiter, residue);
}