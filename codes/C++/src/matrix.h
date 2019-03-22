
template<typename T>
struct Matrix {
    Matrix();
    Matrix(Matrix const& other);
    
    static Matrix<T> gen_random(long long size);
    static bool allclose(std::vector<Matrix<T>> matrices, long double abs_precision, long_double rel_precision);


};