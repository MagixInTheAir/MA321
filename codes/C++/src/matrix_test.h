#ifndef MATRIX_TEST_H
#define MATRIX_TEST_H

#include <stdexcept>
#include <exception>
#include "matrix_def.h"
#include "matrix_test_def.h"

template<typename T>
void Matrix_test<T>::run_all() {
	try {
		Matrix_test<T>::test_constructors();
		Matrix_test<T>::test_operators();
		Matrix_test<T>::test_utilities();
		Matrix_test<T>::test_mathematics();
		Matrix_test<T>::test_generators();
		Matrix_test<T>::test_comparators();
		std::cout << "Tests successful" << std::endl;
	}
	catch (std::exception const& err) {
		std::cout << err.what() << std::endl;
	}
}


template<typename T>
void Matrix_test<T>::test_constructors() {

}

template<typename T>
void Matrix_test<T>::test_operators() {

}

template<typename T>
void Matrix_test<T>::test_utilities() {

}

template<typename T>
void Matrix_test<T>::test_mathematics() {
	{
		Matrix<T> A({ {1,2,3}, {4,5,6}, {7,8,9} });
		Matrix<T> B({ {1,2,3}, {4,5,6},{7,8,9} });
		Matrix<T> res = A.dot(B);
		Matrix<T> conf = Matrix<T>({ {10,12,14}, {22,27,32}, {34,42,50} }) * 3;
		
		if (!res.allclose(conf, 0.001, 0.001)) {
			throw std::logic_error("Test of Matrix::dot failed");
		}
	}

	{
		Matrix<T> A = Matrix<T>::gen_diag({ 3,4,5 });
		T res = A.det();
		T conf = 3*4*5;

		if (!Matrix<T>::close(res, conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::det failed");
		}
	}

	{
		Matrix<T> A ({ {1,2}, {3,4} });
		Matrix<T> res = A.transp();
		Matrix<T> conf({ {1,3}, {2, 4} });

		if (!res.allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::transp failed");
		}
	}

	{
		Matrix<T> A({ {1,2}, {3,4} });
		std::vector<T> res = A.diag();
		std::vector<T> conf({ 1,4 });

		if (!Matrix<T>::allclose(res, conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::diag failed");
		}
	}

	{
		Matrix<T> A(Matrix<T>::gen_diag({ 1,2,3,4 }));
		Matrix<T> res = A.inv();
		Matrix<T> conf(Matrix<T>::gen_diag({ 1, 1./2., 1./3., 1./4. }));

		if (!res.allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::inv failed");
		}
	}
}

template<typename T>
void Matrix_test<T>::test_generators() {

}

template<typename T>
void Matrix_test<T>::test_comparators() {
	{
		{
			Matrix<T> A({ {1,2,3}, {4,5,6}, {7,8,9} });
			Matrix<T> B({ {1,2,3}, {4,5,6},{7,8,9} });

			if (!A.allclose(B, 0.001, 0.001)) {
				throw std::logic_error("Test of Matrix::allclose failed");
			}
		}
	}
}
#endif