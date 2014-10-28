/** Ryan Ordille, 260399372, ECSE420lab2 **/

/**
 * Serial implementation of Gaussian Elimination
 */

#include <stdio.h>
#include <stdlib.h>

#include "tools.h"

double** gauss_elim_serial(double** input, int n, bool make_copy) {
    int k, i, j;
    double factor;
    double** A;

    if (make_copy){
        A = deepcopy_matrix(input);
    } else {
        A = input;
    }


    for (k = 0; k < n -1; k++) {
        for (i = k + 1; i < n; i++) {
            factor = A[i][k] / A[k][k];

            for (j = k + 1; j < n; j++) {
                A[i][j] = A[i][j] - factor * A[k][j];
            }
        }
    }

    return A;
}
