/** Ryan Ordille, 260399372, ECSE420lab2 **/

/**
 * Common functions and settings (in tools.h)
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#include "tools.h"

double** make_matrix(int n){
  // allocate and fill nxn matrix with random numbers
  int i, j;

  double** matrix;

  srand(time(NULL)); // set random seed

  matrix = (double**) malloc(n * sizeof(double*));

  for (i = 0; i < n; i++){
    matrix[i] = (double*) malloc(n * sizeof(double*));

    for (j = 0; j < n; j++){
      // generate random non-zero double
      matrix[i][j] = (double) (rand() % 99) + 1;
    }
  }
  return matrix;
}

double** allocate_matrix(int n) {
  return allocate_submatrix(n, n);
}

double** allocate_submatrix(int width, int height) {
  int i;

  double** submatrix;

  submatrix = (double**) malloc(height * sizeof(double));

  for (i = 0; i < height; i++) {
    submatrix[i] = (double*) malloc(width * sizeof(double));
  }

  return submatrix;
}

double** copy_submatrix(double** input, int width, int y_start, int y_end) {
  int i, j;

  double** submatrix = allocate_submatrix(width, y_end - y_start);

  for (i = 0; i < (y_end - y_start); i++) {
    for (j = 0; j < width; j++) {
      submatrix[i][j] = input[y_start + i][j];
    }
  }

  return submatrix;
}

void print_matrix(double** matrix, int n){
  print_submatrix(matrix, 0, n, 0, n);
}

void print_submatrix(double** matrix, int row_start, int row_end, int col_start, int col_end){
  int i,j;
  double to_print;

  for (i = row_start; i < row_end; i++){
    printf("| ");
    for (j = col_start; j < col_end; j++){
      to_print = matrix[i][j];
      if (i == j)
        printf(ANSI_COLOR_RED "%+e " ANSI_COLOR_RESET, to_print);
      else if (to_print == 0.0)
        printf(ANSI_COLOR_BLACK "%+e " ANSI_COLOR_RESET, to_print);
      else
        printf("%+e ", to_print);
    }
    printf("|\n");
  }
}

double timer(void){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (((double) tv.tv_usec)/1e6);
}

int compare_matrix(double** A, double** B, int n) {
  int i, j;

  int result = TRUE;

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

double** copy_matrix(double** A, int n) {
  return copy_submatrix(A, n, 0, n);
}

void free_matrix(double** A, int n) {
  int i;
  for (i = 0; i < n; i++) {
    free(A[i]);
  }
  free(A);
}

double** replace_submatrix(double** big_matrix, double** submatrix, int n, int y_start, int y_end){
  int i, j;
  for (i = y_start; i < y_end; i++) {
    for (j = 0; j < n; j++) {
      big_matrix[i][j] = submatrix[i - y_start][j];
    }
  }
  return big_matrix;
}
