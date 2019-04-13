#include <vector>
#include <cassert>
#include <iostream>
#include <chrono>
#include <tuple>
#include <fstream>

#include "matrix.h"
#include "matrix_test.h"
#include "algos.h"

unsigned tries = 100;
std::vector<unsigned int> sizes {5, 10, 25, 50, 100};
std::vector<long double> precisions {1e-1, 1e-2, 1e-3, 1e-4, 1e-5};
long double min_gen = -1e9;
long double max_gen = 1e9;

using T = long double;

void run_gaussSeidel_bench() {

	std::ofstream file;
	file.open("gauss-seidel_bench.txt", std::ios::out | std::ios::trunc);

	file << "index, matrix_size, precision, error, time, A, b, x, iterations, min_val, max_val" << std::endl;

	for (unsigned s = 0; s < sizes.size(); s++) {
		unsigned int cursize(sizes[s]);
		std::cout << "TESTING SIZE " << cursize << std::endl;
		for (unsigned p = 0; p < precisions.size(); p++) {
			T curprec(precisions[p]);
			std::cout << "TESTING PRECISION " << curprec << std::endl;
			for (unsigned i = 0; i < tries; i++) {
				auto m = Matrix<T>::gen_random(cursize, min_gen, max_gen);
				Matrix<T> A = Matrix<T>::gen_random(cursize, min_gen, max_gen);
				Matrix<T> b = Matrix<T>::gen_random(cursize, min_gen, max_gen);
				Matrix<T> x0 = Matrix<T>::gen_random(cursize, min_gen, max_gen);


				// START MEASURE TIME
				auto begin(std::chrono::high_resolution_clock::now());

				// CALL FUNCTION
				std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> res0(gaussSeidel(A));
				std::tuple<Matrix<T>, long long, T> result(iter_generale(std::get<1>(res0), std::get<2>(res0), b, x0, curprec, 1e9));

				// END MEASURE TIME
				auto end(std::chrono::high_resolution_clock::now());

				auto duration(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());

				// index, matrix_size, precision, error, time, A, b, x, iterations, min_val, max_val
				std::cout << i << std::endl;
				file << i << ", "
					<< cursize << ", "
					<< curprec << ", "
					<< duration << ", "
					<< A.str() << ", "
					<< b.str() << ", "
					<< std::get<0>(result).str() << ", "
					<< std::get<1>(result) << ", "
					<< std::get<2>(result) << ", "
					<< min_gen << ", "
					<< max_gen << ", "
					<< std::endl;
			}
		}
	}

	file.close();
}


int main() {

	Matrix<T>::run_tests();
	test_algos<T>();
	
	
	run_gaussSeidel_bench();

	
	system("pause");
}
