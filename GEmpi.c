/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "GEmpi.h"

#include <stdio.h> // for printf
#include <stdlib.h> // for malloc/free

#define BUFFER_TAG 0
#define SUBMATRIX_INIT_TAG 1
#define SUBMATRIX_DONE_TAG 2

#define ROOT_RANK 0

#define CONTINUOUS 0
#define CICRULAR 1

double** gauss_elim_parallel_p2p_continuous(double** A, int n){
    int num_processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign this process's rank
    MPI_Status status; // for MPI_Recv() calls

    int blocking = n / num_processes; // blocking factor

    int process; // loop iterator for sending stuff to processes
    int pivot, row, column, denominator, factor; // for gauss elim
    int i; // generic loop iterator

    int submatrix_size = n * (blocking) * sizeof(double); // for submatrix decomp
    int y_start, y_end; // for submatrix decomposition

    double* buffer = malloc(n * sizeof(double));
    double** submatrix;

    for (pivot = 0; pivot < n; pivot++) {
        denominator = A[pivot][pivot];

        if (rank == ROOT_RANK) {
            // parent, make buffer
            for (column = 0; column < n; i++) {
                buffer[column] = A[pivot][column];
            }

            // now send buffer to every other process
            for (process = ROOT_RANK; process < num_processes; process++) {
                MPI_Send(buffer, n, MPI_DOUBLE, process, BUFFER_TAG, MPI_COMM_WORLD);
            }

            // split matrix into local matrices using continuous decomp
            for (process = ROOT_RANK; process < num_processes; process++) {
                // make local copy of submatrix
                y_start = process * blocking;
                y_end = (process * blocking) + blocking;
                submatrix = copy_submatrix(A, n, y_start, y_end);

                // send copy to mapped process
                MPI_Send(submatrix, submatrix_size, MPI_DOUBLE, process, SUBMATRIX_INIT_TAG, MPI_COMM_WORLD);
            }

        }

        // now do elimination for each processor
        for (i = ROOT_RANK; i < num_processes; i++) {
            if (rank == i) {
                // get buffer
                MPI_Recv(buffer, n, MPI_DOUBLE, ROOT_RANK, BUFFER_TAG, MPI_COMM_WORLD, &status);
                // get submatrix
                MPI_Recv(submatrix, submatrix_size, MPI_DOUBLE, ROOT_RANK, BUFFER_TAG, MPI_COMM_WORLD, &status);

                // do elimination for process i
                for (row = 0; row < blocking; row++) {
                    if (denominator) {
                        factor = buffer[pivot] / denominator;
                    } else {
                        // prevent div by zero problems
                        factor = 0.0;
                    }

                    for (column = pivot; column < n; column++) {
                        submatrix[row][column] -= factor * buffer[column];
                    }
                }

                // return eliminated submatrix to root
                MPI_Send(submatrix, submatrix_size, MPI_DOUBLE, ROOT_RANK, SUBMATRIX_DONE_TAG, MPI_COMM_WORLD);
            }
        }

        if (rank == 0) {
            // combine all results
            for (process = ROOT_RANK; i < num_processes; i++) {
                MPI_Recv(submatrix, submatrix_size, MPI_DOUBLE, i, SUBMATRIX_DONE_TAG, MPI_COMM_WORLD, &status);

                // replace the relevant lines of A with this submatrix
                y_start = process*(n / blocking);
                y_end = (process + 1)*(n / blocking);
                replace_submatrix(A, submatrix, n, y_start, y_end);
            }
        }

    }

    free(buffer);
    free_matrix(submatrix, n);
    return A;
}

double** gauss_elim_parallel_broadcast(double** A, int n, int mode){
    return A;
}

void test_parallel(int argc, char** argv) {
    int n = 4;
    int rank;

    double** A = make_matrix(n);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0){
        print_matrix(A, n);
        printf("\n\n");
    }

    A = gauss_elim_parallel_p2p_continuous(A, n);

    if (rank == 0){
        print_matrix(A, n);
        printf("\n");
    }

    MPI_Finalize();
}

void time_parallel_all() {
    return;
}


