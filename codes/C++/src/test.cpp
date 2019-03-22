#include <vector>
#include <cassert>

#include "matrix.h"
#include "gauss-seidel.h"
#include "jacobi.h"
#include "richardson.h"
#include "sor.h"
#include "gmres.h"

std::vector<long long> sizes = {10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000};
std::vector<long double> precisions = {1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};

// IMPORT ALL FUNCTIONS IN A std::vector<auto*>
std::vector<void (*) (Matrix, long double)> functs = {gaussSeidel, jacobi, richardson, sor, gmres};

void output_result(long long size, long double precision, long long cycles) {};
long long start_cycles_measure() {};
long long end_cycles_measure() {};

void main() {

    // FOR EACH MATRIX SIZE :
    for(long long size : sizes) {

        // GENERATE MATRIX
        Matrix m(Matrix::gen_random(size));

        // FOR EACH PRECISION :
        for(long double precision : precisions) {

            std::vector<Matrix m> results;

            // FOR EACH ALGORITHM :
            for(auto* funct : functs) {

                // START MEASURE TIME
                auto start = start_cycles_measure();

                // CALL FUNCTION (MATRIX, PRECISION, MAXITER)
                results.push_back(funct(m, precision));

                // END MEASURE TIME
                auto end = end_cycles_measure();

                // WRITE (SIZE, PRECISION, TIME) TO FILE
                auto time = end - start;
                output_result(size, precision, time);
            }

            assert(Matrix::allclose(results));
        }
    }
}