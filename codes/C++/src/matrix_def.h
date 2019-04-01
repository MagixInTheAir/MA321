#ifndef MATRIX_DEF_H
#define MATRIX_DEF_H

#include "matrix_test_def.h"

template<class T>
class Matrix {
	friend class matrix_test;

protected:
	std::vector<std::vector<T>> data;

public:
	// CONSTRUCTORS & DESTRUCTORS
	Matrix() : data{ {T{}} } {}; // Implemented
	Matrix(std::vector<std::vector<T>> _data) : data{ _data } {}; // Implemented
	Matrix(unsigned int const lines, unsigned int const cols) { // Implemented
		for (unsigned int i = 0; i < lines; i++) {
			this->data.push_back(std::vector<T>(cols, T()));
		};
	}
	template<class T2> Matrix(Matrix<T2> const& other) : data{ other.data } {}; // Implemented

	// UTILITIES
	unsigned int cols() const; // Implemented
	unsigned int lines() const; // Implemented

	// OPERATORS
	std::vector<T> operator[](unsigned int line) const; // Implemented
	template<class T2> auto operator+(Matrix<T2> const& other) const; // Implemented
	template<class T2> auto operator-(Matrix<T2> const& other) const; // Implemented
	Matrix<T> operator-() const; // Implemented
	template<class T2> bool operator==(Matrix<T2> const& other) const; // Implemented
	template<class T2> Matrix<T> operator=(Matrix<T2> const& other); // Implemented

	// MATHEMATICS
	T det() const; // Implemented
	template<class T2> auto dot(Matrix<T2> const& other) const;  // Implemented
	Matrix<T> transp() const; // Implemented
	std::vector<T> diag() const; // Implemented
	Matrix<T> inv() const; // Implemented
	Matrix<T> tri_lo(bool include_diag = false) const; // Implemented
	Matrix<T> tri_up(bool include_diag = false) const; // Implemented
	T highest_eigenval_iteratedPower(T x0, T precision, unsigned long long maxiter) const; // Implemented
	T lowest_eigenval_invIteratedPower(T x0, T precision, unsigned long long maxiter) const; // Implemented
	Matrix<T> abs() const; // Implemented
	T norm() const; // Implemented
	Matrix<T> comatrix() const; // Implemented

	// GENERATORS
	static Matrix<T> gen_random(unsigned int size, T min, T max); // Implemented
	static Matrix<T> gen_random(unsigned int lines, unsigned int cols, T min, T max); // Implemented
	static Matrix<T> gen_diag(unsigned int size, T value = T()); // Implemented
	static Matrix<T> gen_diag(unsigned int lines, unsigned int cols, T value = T()); // Implemented
	static Matrix<T> gen_diag(std::vector<T> values); // Implemented
	static Matrix<T> gen_full(unsigned int size, T value = T()); // Implemented
	static Matrix<T> gen_full(unsigned int lines, unsigned int cols, T value = T()); // Implemented

	// COMPARATORS
	bool allclose(Matrix<T> other, T abs_precision, T rel_precision) const; // Implemented

	// MISC
	static void run_tests(); // Implemented
};

#endif