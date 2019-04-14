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
	Matrix(std::vector<std::vector<T>> _data) : data(_data) {}; // Implemented
	Matrix(unsigned int const lines, unsigned int const cols) : data(lines, std::vector<T>(cols, T())) {}; // Implemented
	template<class T2> Matrix(Matrix<T2> const& other) : data(other.data) {}; // Implemented

	// UTILITIES
	unsigned int cols() const; // Implemented
	unsigned int lines() const; // Implemented
	std::string str() const; // Implemented

	// OPERATORS
	template<class T2> auto operator+(Matrix<T2> const& other) const; // Implemented
	template<class T2> auto operator-(Matrix<T2> const& other) const; // Implemented
	Matrix<T> operator-() const; // Implemented
	template<class T2> bool operator==(Matrix<T2> const& other) const; // Implemented
	template<class T2> Matrix<T> operator=(Matrix<T2> const& other); // Implemented
	template<class T2> auto operator*(T2 const& other) const; // Implemented
	auto operator*(T const& other) const; // Implemented
	std::vector<T>& operator[](unsigned pos); // Implemented

	// MATHEMATICS
	T det() const; // Implemented, tested
	template<class T2> auto dot(Matrix<T2> const& other) const;  // Implemented, tested
	Matrix<T> transp() const; // Implemented, tested
	std::vector<T> diag() const; // Implemented, tested
	Matrix<T> inv_LU() const; // Implemented, tested
	Matrix<T> tri_lo(bool include_diag = false) const; // Implemented, tested
	Matrix<T> tri_up(bool include_diag = false) const; // Implemented, tested
	T highest_eigenval_iteratedPower(std::vector<T> x0, T precision, unsigned long long maxiter) const; // Implemented, tested
	T lowest_eigenval_invIteratedPower(std::vector<T> x0, T precision, unsigned long long maxiter) const; // Implemented, tested
	Matrix<T> abs() const; // Implemented, tested
	T norm() const; // Implemented
	Matrix<T> adj() const; // Implemented, tested
	T cofactor(unsigned int const line, unsigned int const col) const; // Implemented, tested
	std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> decomp_PLU() const; // Implemented, tested
	Matrix<T> decomp_cholesky() const; // Implemented
	std::tuple<Matrix<T>, Matrix<T>> decomp_QR() const;
	// TODO : ?? decomp_SVD() const;
	bool isDiagonal() const; // Implemented
	std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> diagonalize() const;
	Matrix<T> eigenvects() const;
	std::vector<T> eigenvals() const;	
	T rank() const;
	T trace() const; // Implemented
	Matrix<T> pivot() const; // Implemented

	// GENERATORS
	static Matrix<T> gen_random(unsigned int size, T min, T max); // Implemented
	static Matrix<T> gen_random(unsigned int lines, unsigned int cols, T min, T max); // Implemented
	static Matrix<T> gen_diag(unsigned int size, T value = T()); // Implemented
	static Matrix<T> gen_diag(unsigned int lines, unsigned int cols, T value = T()); // Implemented
	static Matrix<T> gen_diag(std::vector<T> values); // Implemented
	static Matrix<T> gen_full(unsigned int size, T value = T()); // Implemented
	static Matrix<T> gen_full(unsigned int lines, unsigned int cols, T value = T()); // Implemented
	static Matrix<T> gen_col(std::vector<T> values); // Implemented
	static Matrix<T> gen_line(std::vector<T> values); // Implemented

	// COMPARATORS
	bool allclose(Matrix<T> other, T abs_precision, T rel_precision) const; // Implemented
	static bool close(T lhs, T rhs, T abs_precision, T rel_precision); // Implemented
	static bool allclose(std::vector<T> lhs, std::vector<T> rhs, T abs_precision, T rel_precision); // Implemented

	// SOLVERS
	static std::vector<T> solve_descent_col(Matrix<T> A, std::vector<T> B); // Implemented, NEEDS OPTIMIZATION (MOVE INTO SOLVE_DESCENT ?)
	static std::vector<T> solve_climb_col(Matrix<T> A, std::vector<T> B); // Implemented, NEEDS OPTIMIZATION (MOVE INTO SOLVE_CLIMB ?)
	static Matrix<T> solve_descent(Matrix<T> A, Matrix<T> B); // Implemented
	static Matrix<T> solve_climb(Matrix<T> A, Matrix<T> B); // Implemented
	static Matrix<T> solve_LU(Matrix<T> A, Matrix<T> B); // Implemented

	// MISC
	static void run_tests(); // Implemented
};

#endif