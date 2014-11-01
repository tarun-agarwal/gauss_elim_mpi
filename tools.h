/** Ryan Ordille, 260399372, ECSE420lab2 **/

#define TRUE 1
#define FALSE 0

#define ANSI_COLOR_BLACK   "\x1b[30m"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

double** make_matrix(int n);
double** allocate_matrix(int n);

void print_matrix(double **matrix, int n);
void print_submatrix(double **matrix, int row_start, int row_end, int col_start, int col_end);

double timer();

int compare_matrix(double** A, double** B, int n);
double** deepcopy_matrix(double** A, int n);
