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
		std::cout << "Matrix tests successful" << std::endl;
	}
	catch (std::exception const& err) {
		std::cout << "ERROR : " << err.what() << std::endl;
	}
};


template<typename T>
void Matrix_test<T>::test_constructors() {

};

template<typename T>
void Matrix_test<T>::test_operators() {

};

template<typename T>
void Matrix_test<T>::test_utilities() {

};

template<typename T>
void Matrix_test<T>::test_mathematics() {
	{
		Matrix<T> A({ {1,2,3}, {4,5,6}, {7,8,9} });
		Matrix<T> B({ {1,2,3}, {4,5,6},{7,8,9} });
		Matrix<T> conf(Matrix<T>({ {10,12,14}, {22,27,32}, {34,42,50} }) * 3);
		
		if (!A.dot(B).allclose(conf, 1e-9, 1e-9)) {
			throw std::logic_error("Test of Matrix::dot failed");
		}
	}

	{
		Matrix<T> A = Matrix<T>::gen_diag({ 3,4,5 });
		T res(A.det());
		T conf(3*4*5);

		if (!Matrix<T>::close(res, conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::det failed");
		}

		Matrix<T> B({ {1,2}, {3,4} });
		T res2(B.det());
		T conf2(-2);

		if (!Matrix<T>::close(res2, conf2, 1e-3, 1e-3)) {
			throw std::logic_error("Test 2 of Matrix::det failed");
		}
	}

	{
		Matrix<T> A ({ {1,2}, {3,4} });
		Matrix<T> conf({ {1,3}, {2, 4} });

		if (!A.transp().allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::transp failed");
		}
	}

	{
		Matrix<T> A({ {1,2,3,4} });
		Matrix<T> conf({ {1},{2},{3},{4} });

		if (!A.transp().allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::transp failed");
		}
	}

	{
		Matrix<T> A({ {1,2}, {3,4} });
		std::vector<T> conf({ 1,4 });

		if (!Matrix<T>::allclose(A.diag(), conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::diag failed");
		}
	}

	{
		Matrix<T> A({ {1,2,3},{4,5,6},{5,4,3} });
		T conf(-9);

		if (!Matrix<T>::close(A.cofactor(0, 0), conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::cofactor failed");
		}
	}

	{
		Matrix<T> A({ {1,2,3},{4,5,6},{5,4,3} });
		Matrix<T> conf({ {-9, 18, -9}, {6, -12, 6}, {-3, 6, -3} });

		if (!A.adj().allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::adj failed");
		}
	}

	{
		Matrix<T> A({ {3,4,3}, {4,8,6}, {3,6,9} });
		Matrix<T> conf({ {1.732051, 0, 0}, {2.309401, 1.632993, 0}, {1.7320508, 1.224745, 2.12132} });


		if (!A.decomp_cholesky().allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::decomp_cholesky failed");
		}
	}

	{
		Matrix<T> A(Matrix<T>::gen_random(10, 10, -100, 100).tri_lo(true));
		Matrix<T> B(Matrix<T>::gen_random(10, 1, -100, 100));
		Matrix<T> X(Matrix<T>::solve_descent(A, B));
		Matrix<T> conf(A.dot(X));


		if (!conf.allclose(B, 1e-9, 1e-9)) {
			throw std::logic_error("Test of Matrix::solve_descent failed");
		}
	}

	{
		Matrix<T> A({ {1,0,0}, {2,3,0}, {4,5,6} });
		Matrix<T> B(Matrix<T>::gen_col({ 7,8,9 }));
		Matrix<T> conf({ {7.}, {-2.}, {-1.5} });
		Matrix<T> X(Matrix<T>::solve_descent(A, B));


		if (!X.allclose(conf, 1e-9, 1e-9)) {
			throw std::logic_error("Test of Matrix::solve_descent failed");
		}
	}

	{

		Matrix<T> A(Matrix<T>::gen_random(10, 10, -100, 100).tri_up(true));
		Matrix<T> B(Matrix<T>::gen_random(10, 1, -100, 100));
		Matrix<T> X(Matrix<T>::solve_climb(A, B));
		Matrix<T> conf(A.dot(X));


		if (!conf.allclose(B, 1e-9, 1e-9)) {
			throw std::logic_error("Test of Matrix::solve_climb failed");
		}
	}

	{
		Matrix<T> A(Matrix<T>::gen_random(4, 4, -100, 100));
		std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> res(A.decomp_PLU());
		Matrix<T> P(std::get<0>(res));
		Matrix<T> L(std::get<1>(res));
		Matrix<T> U(std::get<2>(res));
		Matrix<T> PA(P.dot(A));
		Matrix<T> LU(L.dot(U));

		if (!PA.allclose(LU, 1e-3, 1e-3)) {
			std::string msg("Matrices :\n PA : " + P.str() + ", " + A.str() + "\nLU : " + L.str() + "," + U.str());
			throw std::logic_error("Test of Matrix::decomp_LUP failed\n" + msg);
		}
	}

	{
		Matrix<T> A(Matrix<T>::gen_random(3, 3, -100, 100));
		Matrix<T> I(Matrix<T>::gen_diag(3, T(1)));
		Matrix<T> res(A.dot(A.inv_LU()));
		
		if (!res.allclose(I, 1e-9, 1e-9)) {
			throw std::logic_error("Test of Matrix::inv_LU failed");
		}
	}

	{
		Matrix<T> A({ {1,2,3},{3,2,1},{2,1,3} });
		Matrix<T> conf({ {-5./12., 3./12., 4./12.}, {7./12., 3./12., -8./12.}, {1./12., -3./12., 4./12.} });


		if (!A.inv_LU().allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::inv_LU failed");
		}
	}

	{
		Matrix<T> A({ {1,2},{3,4} });
		Matrix<T> res(A.tri_lo(true));
		Matrix<T> conf({ {1,0},{3,4} });

		if (!res.allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::tri_lo failed");
		}
	}

	{
		Matrix<T> A({ {1,2},{3,4} });
		Matrix<T> res(A.tri_lo(false));
		Matrix<T> conf({ {0,0},{3,0} });

		if (!res.allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::tri_lo failed");
		}
	}

	{
		Matrix<T> A({ {1,2},{3,4} });
		Matrix<T> res(A.tri_up(true));
		Matrix<T> conf({ {1,2},{0,4} });

		if (!res.allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::tri_up failed");
		}
	}

	{
		Matrix<T> A({ {1,2},{3,4} });
		Matrix<T> res = A.tri_up(false);
		Matrix<T> conf({ {0,2},{0,0} });

		if (!res.allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::tri_up failed");
		}
	}

	{
		Matrix<T> A({ {1,-2},{-3,4} });
		Matrix<T> res = A.abs();
		Matrix<T> conf({ {1,2},{3,4} });

		if (!res.allclose(conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::abs failed");
		}
	}

	{
		Matrix<T> A(Matrix<T>::gen_diag({ 1,2,3 }));
		T res(A.highest_eigenval_iteratedPower({ {100},{100}, {100} }, 1e-6, 100000));
		T conf(3);

		if (!Matrix<T>::close(res, conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::highest_eigenval_iteratedPower failed");
		}
	}

	{
		Matrix<T> A(Matrix<T>::gen_diag({ 1,2,3 }));
		T res(A.lowest_eigenval_invIteratedPower({ {100},{100}, {100} }, 1e-6, 100000));
		T conf(1);

		if (!Matrix<T>::close(res, conf, 1e-3, 1e-3)) {
			throw std::logic_error("Test of Matrix::lowest_eigenval_invIteratedPower failed");
		}
	}


};

template<typename T>
void Matrix_test<T>::test_generators() {

};

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
};


#endif