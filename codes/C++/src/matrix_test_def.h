#ifndef MATRIX_TEST_DEF_H
#define MATRIX_TEST_DEF_H

#include "matrix_def.h"

template<typename T>
class Matrix_test {
public:
	static void run_all();

	static void test_constructors();
	static void test_operators();
	static void test_utilities();
	static void test_mathematics();
	static void test_generators();
	static void test_comparators();
};

#endif