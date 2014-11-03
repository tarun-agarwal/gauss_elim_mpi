/** Ryan Ordille, 260399372, ECSE420lab2 **/

/**
 * Serial implementation of Gaussian Elimination
 */

#include "GEserial.h"
#include "tools.h"
#include "stdio.h"

double** gauss_elim_serial(double** A, int n) {
    // convert input matrix A to upper triangular
    // assumes A is square

    int pivot, row, col;
    double factor, denominator;

    for (pivot = 0; pivot < n; pivot++) {
        for (row = pivot+1; row < n; row++) {
            denominator = A[pivot][pivot];

            if (denominator) {
                factor = A[row][pivot] / denominator;
            } else {
                // preventing divide-by-zero problems
                factor = 0.0;
            }

            for (col = pivot; col < n; col++) {
                A[row][col] -= factor * A[pivot][col];
            }
        }
    }

    return A;
}

void time_serial(double** A, int n) {
    double start, end;

    start = timer();

    gauss_elim_serial(A, n);

    end = timer();

    printf("[Serial] (%dx%d) %e s\n", n, n, end - start);

}

void time_serial_all() {
    // time serial implementation for 1024,2048,3072,4096 matrices
    int i, n;
    int sizes[3] = {1024, 2048, 4096};

    for (i = 0; i < 3; i++) {
        n = sizes[i];
        time_serial(make_matrix(n), n);
    }
}

void test_serial() {
    int n = 6;

    double** A = make_matrix(n);

    print_matrix(A, n);

    A = gauss_elim_serial(A, n);

    printf("\n");
    print_matrix(A, n);
}
