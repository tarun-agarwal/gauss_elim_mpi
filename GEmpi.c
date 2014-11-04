/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "GEmpi.h"

#include <stdio.h> // for printf
#include <stdlib.h> // for malloc/free

#define TAG 0

double** gauss_elim_parallel_p2p_continuous(double** A, int n){
    int num_processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);


    int blocking = n / num_processes; // blocking factor

    int pivot, row, col, i, owner;
    int ok_send[1];
    int ok_recv[1];
    double factor, denominator;
    MPI_Status status;

    double *buffer = malloc(n * sizeof(double));

    for (pivot = 0; pivot < n; pivot++) {
        denominator = A[pivot][pivot];

        MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get rank

        if (rank == 0) {
            // parent, make buffer

            // pad with zeros
            for (i = 0; i < pivot; i++) {
                buffer[i] = 0.0;
            }

            for (i = pivot; i < n; i++) {
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

        denominator = buffer[pivot];

        // now do elimination

        for (row = pivot+1; row < n; row++) {
            //while ((rank >= row / blocking) && (rank < (row + blocking) / blocking))
                factor = A[row][pivot] / denominator;

                for (col = pivot; col < n; col++) {
                    A[row][col] -= factor * A[pivot][col];
                }

                MPI_Send(ok_send, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
        }

        if (rank == 0) {
            for (i = 0; i < num_processes; i++) {
                MPI_Recv(ok_recv, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
            }
        }
    }

    free(buffer);

    return A;
}

double** gauss_elim_parallel_p2p_circular(double** A, int n){
    return A;
}

double** gauss_elim_parallel_broadcast_continuous(double** A, int n){
    return A;
}

double** gauss_elim_parallel_broadcast_circular(double** A, int n){
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
        A = gauss_elim_parallel_p2p_continuous(A, n);

        print_matrix(A, n);
        printf("\n");
    }

    MPI_Finalize();


}

void time_parallel_all() {
    return;
}


