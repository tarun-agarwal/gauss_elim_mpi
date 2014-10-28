/** Ryan Ordille, 260399372, ECSE420lab2 **/

/**
 * Common functions and settings (in tools.h)
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#include "tools.h"

double** make_matrix(int n){
  int i, j;

  double** matrix;

  srand(time(NULL));

  matrix = (double**) malloc(n * sizeof(double*));

  for (i = 0; i < n; i++){
    matrix[i] = (double*) malloc(n * sizeof(double*));

    for (j = 0; j < n; j++){
      matrix[i][j] = rand();
    }
  }
  return matrix;
}

double** allocate_matrix(int n) {
  int i;

  double** matrix;

  matrix = (double**) malloc(n * sizeof(double));

  for (i = 0; i < n; i++){
    matrix[i] = (double*) malloc(n * sizeof(double));
  }
  return matrix;
}

void print_matrix(double** matrix, int n){
  print_submatrix(matrix, 0, n, 0, n);
}

void print_submatrix(double** matrix, int row_start, int row_end, int col_start, int col_end){
  int i,j;

  for (i = row_start; i < row_end; i++){
    printf("| ");
    for (j = col_start; j < col_end; j++){
      printf("%e, ", matrix[i][j]);
    }
    printf(" |\n");
  }
}

double timer(void){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (((double) tv.tv_usec)/1e6);
}

bool compare_matrix(double** A, double** B, int n) {
  int i, j;

  bool result = TRUE;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (A[i][j] != B[i][j]) {
        result = FALSE;
        break;
      }
    }
  }
  return result;
}

double** deepcopy_matrix(double** A, int n) {
  int i, j;
  double** B;

  B = allocate_matrix(n);

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      B[i][j] = A[i][j];
    }
  }

  return B;
}

