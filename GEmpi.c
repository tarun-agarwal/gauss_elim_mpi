/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "GEmpi.h"

#include <stdio.h> // for printf
#include <stdlib.h> // for malloc/free

#define TAG 0

#define CONTINUOUS 0
#define CICRULAR 1

double** gauss_elim_parallel_p2p(double** A, int n, int mode){
    int num_processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);


    int blocking = n / num_processes; // blocking factor

    int pivot, row, col, i, j, owner, start_row;
    int ok_send[1];
    int ok_recv[1];
    double factor, denominator;
    MPI_Status status;

    int decomposition;

    double *buffer = malloc(n * sizeof(double));

    for (pivot = 0; pivot < n; pivot++) {
        denominator = A[pivot][pivot];

        MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get rank

        if (rank == 0) {
            // parent, make buffer
            for (i = 0; i < n; i++) {
                buffer[i] = A[pivot][i];
            }

            // now send buffer to every other process
            for (i = 1; i < num_processes; i++) {
                MPI_Send(buffer, n, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            }

        } else {
            // get buffer
            MPI_Recv(buffer, n, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, &status);
        }

        if (buffer[pivot])
            denominator = buffer[pivot];
        else
            denominator = 1.0;

        if ((rank == 0) && (n == num_processes))
            // in the case where n == np, process 0 does nothing
            // so send an OK anyway
            MPI_Send(ok_send, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);

        // split matrix into local matrices

        // now do elimination
        for (start_row = pivot+1; start_row < n; start_row += blocking) {
            if (rank == (start_row / blocking)) {
                for (row = start_row; row < (start_row + blocking); row++) {
                    factor = A[row][pivot] / denominator;

                    for (col = pivot; col < n; col++) {
                        A[row][col] -= factor * A[pivot][col];
                    }

                    printf("Process %d operated on row %d\n", rank, row);

                    MPI_Send(ok_send, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
                }
            }
        }
    }

    if (rank == 0) {
        for (i = 0; i < num_processes; i++) {
            MPI_Recv(ok_recv, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
            printf("Parent received OK from %d\n", i);
        }
    }
    free(buffer);
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

    A = gauss_elim_parallel_p2p(A, n, CONTINUOUS);

    if (rank == 0){
        print_matrix(A, n);
        printf("\n");
    }

    MPI_Finalize();


}

void time_parallel_all() {
    return;
}


