#include <vector>
#include <cassert>
#include <iostream>

#include "matrix.h"
#include "gauss-seidel.h"
#include "jacobi.h"
#include "richardson.h"
#include "sor.h"
#include "gmres.h"


std::vector<unsigned int> sizes {10, 50, 100, 500, 1000, 5000};
std::vector<long double> precisions {1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
long double min_gen = -1e9;
long double max_gen = 1e9;

// IMPORT ALL FUNCTIONS IN A std::vector<auto*>
template<typename T>
std::vector<std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> (*) (Matrix<T> const&, Matrix<T> const&, T)> functs {gaussSeidel<T>}; // , jacobi<T>, richardson<T>, sor<T>, gmres<T>
void output_result(long long size, long double precision, long long cycles) {};
long long start_cycles_measure() { return 0;  };
long long end_cycles_measure() { return 0;  };

using T = long double;


int main() {

    // FOR EACH MATRIX SIZE :
    for(unsigned int size : sizes) {

        // GENERATE MATRIX
        auto m = Matrix<T>::gen_random(size, min_gen, max_gen);
		std::cout << m.lines() << std::endl;

        // FOR EACH PRECISION :
        for(long double precision : precisions) {

            std::vector<std::tuple<Matrix<T>, Matrix<T>, Matrix<T>>> results;
            // FOR EACH ALGORITHM :
			for (unsigned int i = 0; i < functs<T>.size(); i++) {
				auto funct = functs<T>[i];
				Matrix<T> A = Matrix<T>::gen_random(size, min_gen, max_gen);
				Matrix<T> b = Matrix<T>::gen_random(size, min_gen, max_gen);
				
                // START MEASURE TIME
                //auto start = start_cycles_measure();

                // CALL FUNCTION (MATRIX_A, MATRIX_B, PRECISION)
				auto result = funct(A, b, precision);

                // END MEASURE TIME
                //auto end = end_cycles_measure();
				results.push_back(result);

                // WRITE (SIZE, PRECISION, TIME) TO FILE
                //auto time = end - start;
                //output_result(size, precision, time);
				
            }
        }
    }

	system("pause");
}
