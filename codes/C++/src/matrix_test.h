#ifndef MATRIX_TEST_H
#define MATRIX_TEST_H

template<typename T>
class Matrix_test {
public:
	void run_all();
	
	void test_constructors();
	void test_operators();
	void test_utilities();
	void test_mathematics();
	void test_generators();
	void test_comparators();
};

template<typename T>
void Matrix_test<T>::run_all() {
	test_constructors();
	test_operators();
	test_utilities();
	test_mathematics();
	test_generators();
	test_comparators();
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

}
#endif