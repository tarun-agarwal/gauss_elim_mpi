/** Ryan Ordille, 260399372, ECSE420lab2 **/

/**
 * Serial implementation of Gaussian Elimination
 */

double** gauss_elim_serial(double** A, int n) {
    // convert input matrix A to upper triangular
    // assumes A is square

    int k, i, j;
    double factor, denominator;

    for (k = 0; k < n; k++) {
        for (i = k+1; i < n; i++) {
            denominator = A[k][k];

            if (denominator) {
                factor = A[i][k] / denominator;
            } else {
                // preventing divide-by-zero problems
                factor = 0.0;
            }

            for (j = k; j < n; j++) {
                A[i][j] = A[i][j] - factor * A[k][j];
            }
        }
    }

    return A;
}
