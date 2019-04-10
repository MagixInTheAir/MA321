/**********
 GENERIC MATRIX IMPLEMENTATION
 CHRIS DE CLAVERIE <c.de-claverie@pm.me>
 VERSION 1.0
**********/

#ifndef MATRIX_HEADER_INCLUDED
#define MATRIX_HEADER_INCLUDED

#include <vector>
#include <cassert>
#include <utility>
#include <random>
#include <functional>
#include <stdexcept>
#include <complex>
#include <tuple>

#include "matrix_def.h"
#include "matrix_test_def.h"


template<class T>
unsigned int Matrix<T>::cols() const { return this->data[0].size(); };

template<class T>
unsigned int Matrix<T>::lines() const { return this->data.size(); };


template<class T>
template<class T2>
auto Matrix<T>::dot(Matrix<T2> const& other) const {
	if (this->cols() != other.lines()) {
		throw std::logic_error("Dimensions must match for dot product");
	}

    using T3 = decltype(std::declval<T>() * std::declval<T2>());
	Matrix<T3> result(this->lines(), other.cols());

    for(unsigned i = 0; i < result.lines(); i++) {
        for(unsigned j = 0; j < result.cols(); j++) {

            T3 total = 0;
            for(unsigned pos = 0; pos < this->cols(); pos++) {
                total += this->data[i][pos] * other.data[pos][j];
            }
            result[i][j] = total;
        }
    }

    return result;
};

template<class T>
Matrix<T> Matrix<T>::transp() const {
	Matrix<T> result(*this);

    for(unsigned int i = 0; i < this->lines(); i++) {
        for(unsigned int j = 0; j < this->cols(); j++) {
            result[i][j] = this->data[j][i];
        }
    }

    return result;
};

template<class T>
Matrix<T> Matrix<T>::gen_random(unsigned int lines, unsigned int cols, T min, T max) {
	std::default_random_engine generator;
	std::uniform_real_distribution<T> distribution(min, max);
	auto generate = std::bind(distribution, generator);

	Matrix<T> result(lines, cols);

	for (unsigned int i = 0; i < lines; i++) {
		for (unsigned int j = 0; j < cols; j++) {
			result[i][j] = generate();
		}
	}

	return result;
}

template<class T>
Matrix<T> Matrix<T>::gen_random(unsigned int size, T min, T max) { return Matrix<T>::gen_random(size, size, min, max); };

template<class T>
Matrix<T> Matrix<T>::gen_diag(unsigned int lines, unsigned int cols, T value) {
    Matrix<T> result(lines, cols);
    
    for(unsigned int i = 0; i < lines; i++) {
        for(unsigned int j = 0; j < cols; j++) {
			long l(0);
            if (i == j) {
                result[i][j] = value;
            } else {
                result[i][j] = T();
            }
        }
    }

    return result;
};

template<class T>
Matrix<T> Matrix<T>::gen_diag(unsigned int size, T value) {
    return Matrix<T>::gen_diag(size, size, value);
}

template<class T>
Matrix<T> Matrix<T>::gen_diag(std::vector<T> values) {
	auto size = values.size();
    Matrix<T> result(size, size);
    unsigned int pos = 0;

    for(unsigned int i = 0; i < size; i++) {
        for(unsigned int j = 0; j < size; j++) {
            if (i == j) {
                result[i][j] = values[pos];
				pos++;
                i++;
            } else {
                result[i][j] = T(0);
            }
        }
    }

    return result;
}

template<class T>
Matrix<T> Matrix<T>::gen_full(unsigned int lines, unsigned int cols, T value){
	std::vector<std::vector<T>> res_data;

    for(unsigned int i = 0; i < lines; i++) {
        std::vector<T> line;
        for(unsigned int j = 0; j < cols; j++) {
            line.push_back(value);
        }
        res_data.push_back(line);
    }

	return Matrix<T>(res_data);
};

template<class T>
Matrix<T> Matrix<T>::gen_full(unsigned int size, T value) {
    return Matrix<T>::gen_full(size, size, value);
}

template<class T>
T Matrix<T>::det() const {
	if (this->lines() != this->cols()) {
		throw std::logic_error("Matrix must be square");
	}

	// HOT PATH : SMALL MATRICES
	if (this->lines() == 1) {
		return this->data[0][0];
	} else if (this->lines() == 2) {
		return ((this->data[0][0] * this->data[1][1]) - (this->data[0][1] * this->data[1][0]));
	}

	// HOT PATH : DIAGONAL MATRICES
	if (this->isDiagonal()) {
		T det(1);
		for (unsigned p = 0; p < this->lines(); p++) {
			det *= this->data[p][p];
		}
		return det;
	}

	// OTHER CASES
	
	throw std::exception("Path not handled for det computation");

}

template<class T>
T Matrix<T>::cofactor(unsigned int const line, unsigned int const col) const {
	std::vector<std::vector<T>> tempdata;
	for (unsigned int k = 0; k < this->lines(); k++) {
		if (k == line) { continue; }
		std::vector<T> thisline;
		for (unsigned int l = 0; l < this->cols(); l++) {
			if (l == col) { continue; }

			thisline.push_back(this->data[k][l]);
		}
		tempdata.push_back(thisline);
	}

	Matrix<T> temp(tempdata);
	return temp.det();
}

template<class T>
std::vector<T> Matrix<T>::diag() const {
	std::vector<T> diag;
	for (unsigned int i = 0; i < std::min(this->lines(), this->cols()); i++) {
		diag.push_back(this->data[i][i]);
	}
	return diag;
}

template<class T>
Matrix<T> Matrix<T>::adj() const {
	Matrix<T> res(this->lines(), this->cols());
	for (unsigned i = 0; i < this->lines(); i++) {
		for (unsigned j = 0; j < this->cols(); j++) {
			res[i][j] = std::pow(-1, i+j) * this->cofactor(i, j);
		}
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::inv() const {
	std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> PLU(this->decomp_PLU());
	Matrix<T> P(std::get<0>(PLU));
	Matrix<T> L(std::get<1>(PLU));
	Matrix<T> U(std::get<2>(PLU));

	Matrix<T> Y(Matrix<T>::solve_climb(L, P));
	Matrix<T> X(Matrix<T>::solve_descent(U, Y));

	return X;
}

template<class T>
Matrix<T> Matrix<T>::tri_lo(bool include_diag) const {
	Matrix<T> res(this->lines(), this->cols());

	for (unsigned int i = 0; i < this->lines(); i++) {
		for (unsigned int j = 0; j < this->cols(); j++) {
			if (i >= j) {
				if (i == j && !include_diag) {
					continue;
				}
				res[i][j] = this->data[i][j];
			}
			else {
				res[i][j] = T();
			}
		}
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::tri_up(bool include_diag) const {
	Matrix<T> res(this->lines(), this->cols());
	for (unsigned i = 0; i < this->lines(); i++) {
		for (unsigned j = 0; j < this->cols(); j++) {
			if (i <= j) {
				if (i == j && !include_diag) {
					continue;
				}
				res[i][j] = this->data[i][j];
			}
			else {
				res[i][j] = T();
			}
		}
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::abs() const {
	Matrix<T> res(*this);

	for (unsigned i = 0; i < res.lines(); i++) {
		for (unsigned j = 0; j < res.cols(); j++) {
			res[i][j] = std::abs(res[i][j]);
		}
	}

	return res;
}

template<class T>
T Matrix<T>::norm() const {
	if (this->cols() == this->lines()) {
		//auto eigs = this->eigenvals();
		//return std::max(eigs) / std::min(eigs);
		return T();
	}
	else if (this->cols() == 1) {
		T total(0);
		for (unsigned i = 0; i < this->lines(); i++) {
			total += std::pow(this->data[i][0], 2);
		}
		return std::sqrt(total);
	}

	throw std::exception("Path not handled");
};

template<class T>
template<class T2>
auto Matrix<T>::operator+(Matrix<T2> const& other) const {
	Matrix<T> res(this->lines(), this->cols());

	for (unsigned int i = 0; i < this->lines(); i++) {
		for (unsigned int j = 0; j < this->cols(); j++) {
			res[i][j] = this->data[i][j] + other.data[i][j];
		}
	}
	return res;
}

template<class T>
template<class T2>
auto Matrix<T>::operator-(Matrix<T2> const& other) const {
	if (this->lines() != other.lines() || this->cols() != other.cols()) {
		throw std::logic_error("The two matrices are not the same size");
	}

	using T3 = decltype(std::declval<T>() - std::declval<T2>());

	Matrix<T3> res(*this);

	for (unsigned int i = 0; i < this->lines(); i++) {
		for (unsigned int j = 0; j < this->cols(); j++) {
			res[i][j] -= other.data[i][j];
		}
	}

	return res;
};

template<class T>
Matrix<T> Matrix<T>::operator-() const {
	Matrix<T> res(*this);
	for (unsigned int const i = 0; i < this->lines(); i++) {
		for (unsigned int const j = 0; j < this->cols(); j++) {
			res[i][j] *= -1;
		}
	}
	return res;
}

template<class T>
template<class T2>
bool Matrix<T>::operator==(Matrix<T2> const& other) const {
	return this->data == other;
}

template<class T>
template<class T2>
Matrix<T> Matrix<T>::operator=(Matrix<T2> const& other) {
	this->data = other;
	return *this;
};

template<typename T>
T Matrix<T>::highest_eigenval_iteratedPower(std::vector<T> x0, T precision, unsigned long long maxiter) const {
	if (x0.size() != this->cols()) {
		throw std::logic_error("x0 is not a valid vector");
	}

	Matrix<T> Y;
	T norm_m(1);
	T norm_p(1);
	Matrix<T> X(Matrix<T>::gen_col(x0));
	Matrix<T> X_prec(X);
	for (unsigned long long i = 0; i < maxiter; i++) {
		if (norm_m <= precision || norm_p <= precision) {
			break;
		}

		Y = this->dot(X);
		X_prec = X;
		X = Y*(1 / Y.norm());
		norm_m = (X - X_prec).norm();
		norm_p = (X + X_prec).norm();
	}
	
	return this->dot(X).norm();

};

template<typename T>
T Matrix<T>::lowest_eigenval_invIteratedPower(std::vector<T> x0, T precision, unsigned long long maxiter) const {
	return this->inv().highest_eigenval_iteratedPower(x0, precision, maxiter);
};

template<typename T>
void Matrix<T>::run_tests() {
	Matrix_test<T>::run_all();
}

template<typename T>
bool Matrix<T>::allclose(Matrix<T> other, T abs_precision, T rel_precision) const {
	if (this->lines() != other.lines() || this->cols() != other.cols()) {
		throw std::length_error("Matrix must be the same size");
	}

	bool close = true;
	for (unsigned i = 0; i < other.lines(); i++) {
		for (unsigned j = 0; j < other.cols(); j++) {
			if (!Matrix<T>::close(this->data[i][j], other[i][j], abs_precision, rel_precision)) {
				close = false;
			}
		}
	}
	return close;
}

template<class T>
template<class T2>
auto Matrix<T>::operator*(T2 const& other) const {
	using T3 = decltype(std::declval<T>() * std::declval<T2>());
	Matrix<T3> res(this->lines(), this->cols());
	for (unsigned i = 0; i < this->lines(); i++) {
		for (unsigned j = 0; j < this->cols(); j++) {
			res[i][j] = this->data[i][j] * other;
		}
	}

	return res;
}

template<class T>
auto Matrix<T>::operator*(T const& other) const {
	Matrix<T> res(this->lines(), this->cols());
	for (unsigned i = 0; i < this->lines(); i++) {
		for (unsigned j = 0; j < this->cols(); j++) {
			res[i][j] = this->data[i][j] * other;
		}
	}

	return res;
}

template<class T>
bool Matrix<T>::close(T lhs, T rhs, T abs_precision, T rel_precision) {
	if (std::abs(lhs - rhs) <= abs_precision + (rel_precision * std::abs(rhs))) {
		return true;
	}

	return false;
}

template<class T>
bool Matrix<T>::allclose(std::vector<T> lhs, std::vector<T> rhs, T abs_precision, T rel_precision) {
	if (lhs.size() != rhs.size()) {
		return false;
	}

	bool close = true;
	for (unsigned i = 0; i < lhs.size(); i++) {
		if (!Matrix<T>::close(lhs[i], rhs[i], abs_precision, rel_precision)) {
			close = false;
			break;
		}
	}

	return close;
}

template<class T>
Matrix<T> Matrix<T>::gen_col(std::vector<T> values) {
	Matrix<T> res(values.size(), 1);
	for (unsigned i = 0; i < values.size(); i++) {
		res[i][0] = values[i];
	}
	return res;
}

template<class T>
Matrix<T> Matrix<T>::gen_line(std::vector<T> values) {
	Matrix<T> res(1, values.size());
	for (unsigned i = 0; i < values.size(); i++) {
		res[0][i] = values[i];
	}
	return res;
}

template<class T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> Matrix<T>::decomp_PLU() const {
	// https://www.quantstart.com/articles/LU-Decomposition-in-Python-and-NumPy

	if (this->lines() != this->cols()) {
		throw std::logic_error("LUP decomposition impossible on non-square matrix");
	}

	unsigned n = this->lines();

	Matrix<T> L(Matrix<T>::gen_full(n, T(0)));
	Matrix<T> U(Matrix<T>::gen_full(n, T(0)));
	Matrix<T> P = this->pivot();
	auto PA = P.dot(*this);


	for (unsigned j = 0; j < n; j++) {
		L[j][j] = T(1);

		for (unsigned i = 0; i < j + 1; i++) {
			T sum(0);
			for (unsigned k = 0; k < i; k++) {
				sum += (U[k][j] * L[i][k]);
			}

			U[i][j] = PA[i][j] - sum;
		}

		for (unsigned i = j;  i < n; i++) {
			T sum(0);
			for (unsigned k = 0; k < j; k++) {
				sum += U[k][j] * L[i][k];
			}

			L[i][j] = (PA[i][j] - sum) / U[j][j];
		}
	}



	return std::make_tuple(P, L, U);
}

template<class T>
bool Matrix<T>::isDiagonal() const {
	bool diagonal = true;
	for (unsigned i = 0; i < this->lines(); i++) {
		for (unsigned j = 0; j < this->cols(); j++) {
			if (i == j) { continue; }
			if (this->data[i][j] != 0) {
				diagonal = false;
				break;
			}
		}
	}
	return diagonal;
}

template<class T>
T Matrix<T>::trace() const {
	T res(0);
	for (unsigned i = 0; i < std::min(this->lines(), this->cols()); i++) {
		res += this->data[i][i];
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::decomp_cholesky() const {
	Matrix<T> L(this->lines(), this->cols());
	L[0][0] = std::sqrt(this->data[0][0]);
	for (unsigned i = 0; i < this->lines(); i++) {
		for (unsigned j = i; j < this->lines(); j++) {
			if (i == j) {
				T S(0);
				for (unsigned k = 0; k < i; k++) {
					S += std::pow(L[i][k], 2);
				}
				L[i][j] = std::sqrt(this->data[i][j] - S);
			}
			else {
				T S(0);
				for (unsigned k = 0; k < i; k++) {
					S += L[i][k] * L[j][k];
				}
				L[j][i] = (this->data[i][j] - S) / L[i][i];
			}
		}
	}

	return L;
};

template<class T>
Matrix<T> Matrix<T>::solve_descent(Matrix<T> A, Matrix<T> B) {
	if (A.lines() != A.cols() || B.lines() != A.lines() || B.cols() != 1) {
		throw std::logic_error("Matrix dimensions failed for solve_descent");
	}

	// Creating Tiaug 
	Matrix<T> Taug(A.lines(), A.cols() + B.cols());
	for (unsigned i = 0; i < A.lines(); i++) {
		for (unsigned j = 0; j < A.cols(); j++) {
			Taug[i][j] = A[i][j];
		}
	}
	for (unsigned k = 0; k < B.lines(); k++) {
		for (unsigned l = 0; l < B.cols(); l++) {
			Taug[k][A.cols() + l] = B[k][l];
		}
	}

	// Solving
	std::vector<T> X(Taug.lines(), 0);
	for (unsigned i = 0; i < Taug.lines(); i++) {
		T somme(0);
		for (unsigned k = 0; k < i; k++) {
			somme += X[k] * Taug[i][k];
		}
		X[i] = (Taug[i][Taug.cols() - 1] - somme) / Taug[i][i];
	}

	return Matrix<T>::gen_col(X);
}

template<class T>
Matrix<T> Matrix<T>::solve_climb(Matrix<T> A, Matrix<T> B) {
	if (A.lines() != A.cols() || B.lines() != A.lines()) {
		throw std::logic_error("Matrix dimensions failed for solve_climb");
	}

	// Creating Taug 
	Matrix<T> Taug(A.lines(), A.cols() + B.cols());
	for (unsigned i = 0; i < A.lines(); i++) {
		for (unsigned j = 0; j < A.cols(); j++) {
			Taug[i][j] = A[i][j];
		}
	}
	for (unsigned k = 0; k < B.lines(); k++) {
		for (unsigned l = 0; l < B.cols(); l++) {
			Taug[k][A.cols() + l] = B[k][l];
		}
	}

	// Solving
	std::vector<T> X(Taug.lines(), 0);
	for (unsigned ip = Taug.lines(); ip > 0; ip--) { // VA BOUCLER A L'INFINI
		unsigned i = ip - 1;
		T somme(0);
		for (unsigned k = i; k < Taug.lines(); k++) {
			somme += X[k] * Taug[i][k];
		}
		X[i] = (Taug[i][Taug.cols()-1] - somme) / Taug[i][i];
	}

	return Matrix<T>::gen_col(X);
}

template<class T>
std::string Matrix<T>::str() const {
	std::string str;
	str += "{";
	for (unsigned i = 0; i < this->lines(); i++) {
		str += "{";
		for (unsigned j = 0; j < this->cols(); j++) {
			str += std::to_string(this->data[i][j]);
			if (j != this->cols() - 1) {
				str += ",";
			}
		}
		str += "}";
		if (i != this->lines() - 1) {
			str += ",";
		}
	}
	str += "}";
	return str;
}

template<class T>
std::vector<T>& Matrix<T>::operator[](unsigned pos) {
	return this->data[pos];
}

template<class T>
Matrix<T> Matrix<T>::pivot() const {
	if (this->lines() != this->cols()) {
		throw std::logic_error("Matrix must be square");
	}
	unsigned m = this->lines();
	Matrix<T> Id_mat(Matrix<T>::gen_diag(m, T(1)));

	for (unsigned j = 0; j < m; j++) {
		T maxval(this->data[0][j]);
		unsigned maxline(0);
		for (unsigned i = 0; i < m; i++) {
			if (std::abs(this->data[i][j]) > maxval) {
				maxval = this->data[i][j];
				maxline = i;
			}
		}

		if (j != maxline) {
			std::swap(Id_mat[j], Id_mat[maxline]);
		}
	}

	return Id_mat;
}

#endif