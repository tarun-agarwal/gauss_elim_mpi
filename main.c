/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "tools.h"
#include "GESerial.h"

#include <stdio.h>
#include <stdlib.h>

int main() {

    int matrix_size =8;

    double** A = make_matrix(matrix_size);

    print_matrix(A, matrix_size);
    A = gauss_elim_serial(A, matrix_size);

    printf("\n");

    print_matrix(A, matrix_size);

    return 0;
}
