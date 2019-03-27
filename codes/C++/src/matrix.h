#include <vector>
#include <cassert>

template<typename T>
class Matrix {
protected :
    std::vector<std::vector<T>> data;

public :
    Matrix(); // Implemented
    template<typename T2> Matrix(Matrix<T2> const& other); // Implemented
    long long cols(); // Implemented
    long long lines(); // Implemented

    std::vector<T> operator[](long long pos); // Implemented
    Matrix operator+(Matrix const& other);
    Matrix operator-(Matrix const& other);
    Matrix operator-();
    bool operator==();
    Matrix const& operator=(Matrix const& other);

    long long det();
    template<typename T2> Matrix<decltype(std::declval<T>() * std::declval<T2>())> dot(Matrix<T2> const& other); // Implemented
    Matrix<T> transp(); // Implemented
    std::vector<T> diag();
    
    static Matrix<T> gen_random(long long size); // Implemented
    static Matrix<T> gen_random(long long lines, long long cols); 
    static Matrix<T> gen_diag(long long size, T value = T()); // Implemented
    static Matrix<T> gen_diag(long long lines, long long cols, T value = T()); // Implemented
    static Matrix<T> gen_diag(std::vector<T> values); // Implemented
    static Matrix<T> gen_full(long long size); // Implemented
    static Matrix<T> gen_full(long long lines, long long cols, T value = T()); // Implemented
    static Matrix<T> gen_uninitialized(long long size); // Implemented
    static Matrix<T> gen_uninitialized(long long lines, long long cols); // Implemented

    static bool allclose(std::vector<Matrix<T>> matrices, long double abs_precision, long_double rel_precision);
};

template<typename T>
Matrix<T>::Matrix() : data{{T()}} {};

template<typename T, typename T2>
Matrix<T>::Matrix(Matrix<T2> const& other) : data(other.data) {};

template<typename T>
long long Matrix<T>::cols() { return data.size(); };

template<typename T>
long long Matrix<T>::lines() { return data[0].size(); }

template<typename T>
std::vector<T> Matrix<T>::operator[](long long pos) { return this->data[pos]; }

template<typename T, typename T2>
Matrix<decltype (std::declval<T>() * std::declval<T2>())> Matrix<T>::dot(Matrix<T2> const& other) {
    assert(this->cols() == other.lines());

    decltype (std::declval<T>()*std::declval<T2>() T3;
    Matrix<T3> result = Matrix<T3>::gen_uninitialized(this->lines() * other.cols());

    for(long long i = 0; i < result.lines(); i++) {
        for(long long j = 0; j < result.cols(); i++) {

            long long total = 0;
            for(long long pos = 0; pos < this->cols(); pos++) {
                total += ((*this)[i][pos] * other[pos][j]);
            }
            result.data[i][j] = total;
        }
    }

    return result;
};

template<typename T>
Matrix<T> Matrix<T>::transp() {
    Matrix<T> result = Matrix<T>::gen_uninitialized(this->cols(), this->lines());
    for(long long i = 0; i < this->lines(); i++) {
        for(long long j = 0; j < this->cols(); j++) {
            result[i][j] = this[j][i];
        }
    }

    return result;
};



template<typename T>
Matrix<T>::gen_uninitialized(long long lines, long long cols) {
    Matrix<T> result;
    for(int i = 0; i < lines; i++) {
        std::vector<T> line;
        for(int j = 0; j < lines; j++) {
            line.push_back(T());
        }
        result.data.push_back(line);
    }
}

template<typename T>
Matrix<T>::gen_uninitialized(long long size) { 
    return gen_uninitialized(size, size);
}

template<typename T>
Matrix<T>::gen_random(long long lines, long long cols) {

}

template<typename T>
Matrix<T> Matrix<T>::gen_random(long long size) { return gen_random(size, size); };

template<typename T>
Matrix<T> Matrix<T>::gen_diag(long long lines, long long cols, T value = T()) {
    Matrix<T> result;
    
    for(long long i = 0; i < lines; i++) {
        for(long long j = 0; j < cols; j++) {
            long l
            if (i == j) {
                result[i][j] = value;
            } else {
                result[i][j] = T();
            }
        }
    }

    return result;
};

template<typename T>
Matrix<T> Matrix<T>::gen_diag(long long size, long long value = T()) {
    return Matrix<T>::gen_diag(size, size, value);
}

template<typename T>
Matrix<T> Matrix<T>::gen_diag(std::vector<T> values) {
    Matrix<T> result;
    let pos = 0;

    for(long long i = 0; i < values.size(); i++) {
        for(long long j = 0; j < values.size(); j++) {
            long l
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

template<typename T>
Matrix<T> Matrix<T>::gen_full(long long lines, long long cols, long long value = T()){
    for(long long i = 0; i < lines; i++) {
        std::vector<T> line;
        for(long long j = 0; j < cols; j++) {
            line.push_back(value);
        }
        this->data.push_back(line);
    }
};

template<typename T>
Matrix<T> Matrix<T>::gen_full(long long size, long long value = 1) {
    return Matrix<T>::gen_full(size, size, value);
}