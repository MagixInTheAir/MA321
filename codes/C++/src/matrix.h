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


template<class T>
class Matrix {
protected:
	std::vector<std::vector<T>> data;

public:
	// CONSTRUCTORS & DESTRUCTORS
	Matrix() : data{ {T{}} } {}; // Implemented
	Matrix(std::vector<std::vector<T>> _data) : data{ _data } {}; // Implemented
	Matrix(unsigned int const lines, unsigned int const cols) { // Implemented
		for (unsigned int i = 0; i < lines; i++) { this->data.push_back(std::vector<T>(cols, T())); }
	};
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
	T highest_eigenval_iteratedPower(T precision) const;
	T lowest_eigenval_invIteratedPower(T precision) const;
    Matrix<T> abs() const; // Implemented
    T norm() const; // Implemented
	Matrix<T> comatrix() const; // Implemented

    // GENERATORS
    static Matrix<T> gen_random(unsigned int size, T min, T max); // Implemented
    static Matrix<T> gen_random(unsigned int lines, unsigned int cols, T min, T max); 
    static Matrix<T> gen_diag(unsigned int size, T value = T()); // Implemented
    static Matrix<T> gen_diag(unsigned int lines, unsigned int cols, T value = T()); // Implemented
    static Matrix<T> gen_diag(std::vector<T> values); // Implemented
    static Matrix<T> gen_full(unsigned int size, T value = T()); // Implemented
    static Matrix<T> gen_full(unsigned int lines, unsigned int cols, T value = T()); // Implemented

	// COMPARATORS
    static bool allclose(std::vector<Matrix<T>> matrices, long double abs_precision, long double rel_precision);
};

template<class T>
unsigned int Matrix<T>::cols() const { return this->data.size(); };

template<class T>
unsigned int Matrix<T>::lines() const { return this->data[0].size(); };

template<class T>
std::vector<T> Matrix<T>::operator[](unsigned int line) const { return this->data[line]; };

template<class T>
template<class T2>
auto Matrix<T>::dot(Matrix<T2> const& other) const { // Here, auto is decltype(std::declval<T>() * std::declval<T2>())
    //static_assert(this->cols() == other.lines());
	//static_assert(this->lines() == other.cols());

    using T3 = decltype(std::declval<T>() * std::declval<T2>());
    Matrix<T3> result = Matrix<T3>::gen_full(this->lines(), other.cols());

    for(unsigned int i = 0; i < result.lines(); i++) {
        for(unsigned int j = 0; j < result.cols(); i++) {

            T3 total = 0;
            for(unsigned int pos = 0; pos < this->cols(); pos++) {
                total += this->data[i][pos] * other[pos][j];
            }
            result.data[i][j] = total;
        }
    }

    return result;
};

template<class T>
Matrix<T> Matrix<T>::transp() const {
	Matrix<T> result(*this);

    for(unsigned int i = 0; i < this->lines(); i++) {
        for(unsigned int j = 0; j < this->cols(); j++) {
            result.data[i][j] = this->data[j][i];
        }
    }

    return result;
};

template<class T>
Matrix<T> Matrix<T>::gen_random(unsigned int lines, unsigned int cols, T min, T max) {
	std::default_random_engine generator;
	std::uniform_real_distribution<T> distribution(min, max);
	auto generate = std::bind(distribution, generator);

	Matrix<T> result;
	for (unsigned int i = 0; i < lines; i++) {
		std::vector<T> line;
		for (unsigned int j = 0; j < cols; j++) {
			line.push_back(generate());
		}
		result.data.push_back(line);
	}

	return result;
}

template<class T>
Matrix<T> Matrix<T>::gen_random(unsigned int size, T min, T max) { return Matrix<T>::gen_random(size, size, min, max); };

template<class T>
Matrix<T> Matrix<T>::gen_diag(unsigned int lines, unsigned int cols, T value) {
    Matrix<T> result;
    
    for(unsigned int i = 0; i < lines; i++) {
        for(unsigned int j = 0; j < cols; j++) {
			long l;
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
    Matrix<T> result;
    unsigned int pos = 0;

    for(unsigned int i = 0; i < values.size(); i++) {
        for(unsigned int j = 0; j < values.size(); j++) {
            if (i == j) {
                result[i][j] = values[pos];
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
	T det = T();
	if (this->data[0].size() == 1) {
		det = this->data[0][0];
	}
	else {
		for (unsigned int p = 0; p < this->cols(); p++) {

			std::vector<std::vector<T>> sub_data(this->data.size() - 1);
			for (unsigned int i = 0; i < this->lines(); i++) {
				std::vector<T> line;
				for (unsigned int j = 0; j < this->cols(); j++) {
					if (i != j) {
						line.push_back(this->data[i][j]);
					}
				}
				sub_data.push_back(line);
			}

			Matrix<T> sub(sub_data);
			det += std::pow(-1, p+1) * data[0][p] * sub.det();
		}
	}
	return det;
}

template<class T>
std::vector<T> Matrix<T>::diag() const {
	std::vector<T> diag;
	for (unsigned int i = 0; i < std::min(this->lines(), this->cols()); i++) {
		diag.push_back(this->data[i][i]);
	}
}

template<class T>
Matrix<T> Matrix<T>::comatrix() const {
	Matrix<T> res(this->lines(), this->cols());
	for (unsigned int i = 0; i < this->lines(); i++) {
		for (unsigned int j = 0; j < this->cols(); j++) {

			std::vector<std::vector<T>> tempdata;
			for (unsigned int k = 0; k < this->lines(); k++) {
				std::vector<T> line;
				for (unsigned int l = 0; l < this->cols(); l++) {
					if (k != l) {
						line.push_back(this->data[k][l]);
					}
				}
				tempdata.push_back(line);
			}

			Matrix<T> temp(tempdata);
			res.data[i][j] = temp.det();
		}
	}
	return res;
}

template<class T>
Matrix<T> Matrix<T>::inv() const {
	return this->comatrix().transp();
}

template<class T>
Matrix<T> Matrix<T>::tri_lo(bool include_diag) const {
	Matrix<T> res(this->lines(), this->cols());
	for (unsigned int i = 0; i < this->lines(); i++) {
		for (unsigned int j = 0; j < this->cols(); j++) {
			if (i <= j) {
				if (i == j && !include_diag) {
					continue;
				}
				res.data[i][j] = this->data[i][j];
			}
			else {
				res.data[i][j] = T();
			}
		}
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::tri_up(bool include_diag) const {
	Matrix<T> res(this->lines(), this->cols());
	for (int i = 0; i < this->lines(); i++) {
		for (int j = 0; j < this->cols(); j++) {
			if (i >= j) {
				if (i == j && !include_diag) {
					continue;
				}
				res.data[i][j] = this->data[i][j];
			}
			else {
				res.data[i][j] = T();
			}
		}
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::abs() const {
	Matrix<T> res(*this);

	for (int i = 0; i < res.lines(); i++) {
		for (int j = 0; j < res.cols(); j++) {
			res[i][j] = std::abs(res[i][j]);
		}
	}

	return res;
}

template<class T>
T Matrix<T>::norm() const {
	auto eigs = this->eigenvals();
	return std::max(eigs) / std::min(eigs);
};

template<class T>
template<class T2>
auto Matrix<T>::operator+(Matrix<T2> const& other) const {
	Matrix<T> res(this->lines(), this->cols());

	for (unsigned int i = 0; i < this->lines(); i++) {
		for (unsigned int j = 0; j < this->cols(); j++) {
			res.data[i][j] = this->data[i][j] + other.data[i][j];
		}
	}
	return res;
}

template<class T>
template<class T2>
auto Matrix<T>::operator-(Matrix<T2> const& other) const {
	// static_assert(this->lines() == other.lines(), "The two matrices are not the same size");
	// static_assert(this->cols() == other.cols(), "The two matrices are not the same size");

	using T3 = decltype(std::declval<T>() - std::declval<T2>());

	Matrix<T3> res(*this);

	for (unsigned int i = 0; i < this->lines(); i++) {
		for (unsigned int j = 0; j < this->cols(); j++) {
			res[i][j] -= other[i][j];
		}
	}

	return res;
};

template<class T>
Matrix<T> Matrix<T>::operator-() const {
	Matrix<T> res(*this);
	for (unsigned int const i = 0; i < this->lines(); i++) {
		for (unsigned int const j = 0; j < this->cols(); j++) {
			res.data[i][j] *= -1;
		}
	}
	return res;
}

template<class T>
template<class T2>
bool Matrix<T>::operator==(Matrix<T2> const& other) const {
	return this->data == other.data;
}

template<class T>
template<class T2>
Matrix<T> Matrix<T>::operator=(Matrix<T2> const& other) {
	this->data = other.data;
};

template<typename T>
T Matrix<T>::highest_eigenval_iteratedPower(T precision) const {

};

template<typename T>
T Matrix<T>::lowest_eigenval_invIteratedPower(T precision) const {
	
};



#endif