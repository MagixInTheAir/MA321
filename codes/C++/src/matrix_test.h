#ifndef MATRIX_TEST_H
#define MATRIX_TEST_H

#include <stdexcept>
#include "matrix_def.h"
#include "matrix_test_def.h"

template<typename T>
void Matrix_test<T>::run_all() {
	Matrix_test<T>::test_constructors();
	Matrix_test<T>::test_operators();
	Matrix_test<T>::test_utilities();
	Matrix_test<T>::test_mathematics();
	Matrix_test<T>::test_generators();
	Matrix_test<T>::test_comparators();
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

}

template<typename T>
void Matrix_test<T>::test_generators() {

}

template<typename T>
void Matrix_test<T>::test_comparators() {
	{
		Matrix<T> A({ {1,2,3}, {4,5,6}, {7,8,9} });
		Matrix<T> B({ {1,2,3}, {4,5,6},{7,8,9} });

		if (!A.allclose(B, 0.001, 0.001)) {
			throw std::logic_error("Test of Matrix::allclose failed");
		}
	}
}
#endif