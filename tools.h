/** Ryan Ordille, 260399372, ECSE420lab2 **/

#define TRUE 1
#define FALSE 0

double** make_matrix(int n);
double** allocate_matrix(int n);

void print_matrix(double **matrix, int n);
void print_submatrix(double **matrix, int row_start, int row_end, int col_start, int col_end);

double timer();

int compare_matrix(double** A, double** B, int n);
double** deepcopy_matrix(double** A, int n);

typedef int bool;
