#pragma once

#include <stdexcept>

#include "iter_generale.h"
#include "gauss-seidel.h"
#include "jacobi.h"
#include "sor.h"

template<typename T>
void test_gaussSeidel() {
	{
		Matrix<T> A({ {2,-1,1}, {2,2,2}, {-1,-1,2} });
		Matrix<T> B({ {3}, {12}, {3} });
		Matrix<T> conf({ {1}, {2}, {3} });

		std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> decomp(gaussSeidel(A));
		Matrix<T> m(std::get<0>(decomp));
		Matrix<T> m_inv(std::get<1>(decomp));
		Matrix<T> iter(std::get<2>(decomp));

		Matrix<T> x0({ {10}, {10}, {10} });
		T epsilon(1e-9);
		long long maxiter(1000000);

		std::tuple<Matrix<T>, long long, T> res(iter_generale(m_inv, iter, B, x0, epsilon, maxiter));
		Matrix<T> x(std::get<0>(res));

		if (!x.allclose(conf, 1e-9, 1e-9)) {
			throw std::logic_error("Test of gauss_seidel failed");
		}
	}
}

template<typename T>
void test_algos() {
	try {
		test_gaussSeidel<T>();
		std::cout << "Algorithm tests successful" << std::endl;
	}
	catch (std::exception const& err) {
		std::cout << "ERROR : " << err.what() << std::endl;
	}

}