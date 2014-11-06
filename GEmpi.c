/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "GEmpi.h"

#include <stdio.h> // for printf
#include <stdlib.h> // for malloc/free

#define BUFFER_TAG 0
#define SUBMATRIX_INIT_TAG 1
#define SUBMATRIX_DONE_TAG 2
#define OK_TAG 3

#define ROOT_RANK 0

#define CONTINUOUS 0
#define CICRULAR 1

double** gauss_elim_parallel_p2p_continuous(double** matrix, int n) {
    // DECLARATIONS & INITALIZATIONS

    // for MPI calls:
    int np, rank; // from MPI_Init()
    MPI_Status status; // for MPI_Recv()

    MPI_Comm_size(MPI_COMM_WORLD, &np); // number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // this process' rank

    int target_row, sender, receiver; // target rank for send/recv

    int row_owner;

    // for matrix allocation
    int bf = n / np; // blocking factor, i.e. number of rows per local matrix
    int local_size = n * bf; // # of elements per MPI_Send of submatrix
    int buffer_size;

    double* buffer; // allocate row buffer
    double** local; // local submatrix for each process

    // for timer:
    double comp_start, comp_end, comm_start, comm_end;
    double total_comp, total_comm;

    // loop iterators:
    int process; // looping over processes
    int pivot, row, column; // looping through matrices

    // for gauss elimination:
    double denominator, factor;

    comp_start = timer(); // start total timer after init of global matrix;

    // ALLOCATE LOCAL SUBMATRICES
    for (process = 0; process < np; process++) {
        if (rank == process) {
            // local submatrix should be of size n*bf*(8 bytes)
            local = allocate_submatrix(n, bf);
        }
    }

    // INIT LOCAL SUBMATRICES
    for (row = 0; row < n; row++) {
        if (row >= bf) {
            // not owned by root, so send to relevant process
            if (rank == ROOT_RANK) {
                // send to applicable process
                receiver = row / bf;
                MPI_Send(matrix[row], n, MPI_DOUBLE, receiver, SUBMATRIX_INIT_TAG, MPI_COMM_WORLD);
            } else {
                // recv submatrix from root
                target_row = row % bf;

                comm_start = timer();
                MPI_Recv(local[target_row], n, MPI_DOUBLE, ROOT_RANK, SUBMATRIX_INIT_TAG, MPI_COMM_WORLD, &status);
                comm_end = timer();
                total_comm += comm_end - comm_start;
            }
        } else {
            // first bf rows so allocate root's portion directly
            if (rank == ROOT_RANK) {
                target_row = row % bf;
                for (column = 0; column < n; column++) {
                    local[target_row][column] = matrix[row][column];
                }
            }
            // else break?
        }
    }

    // PERFORM GAUSS ELIMINATION
    for (pivot = 0; pivot < n; pivot++) {
        // allocate buffer
        buffer_size = n - pivot;
        buffer = malloc(buffer_size * sizeof(double));

        // figure out who owns this pivot row
        row_owner = row / bf; // TODO

        // owner sends pivot row to everyone via buffer
        if (rank == row_owner) {
            // need to send this row to all other processes
            for (column = pivot; column < n; column++) {
                buffer[column - pivot] = local[pivot % bf][column];
            }

            for (process = 0; process < np; process++) {
                if (process != rank) {
                    // don't send to self but to everyone else
                    MPI_Send(buffer, buffer_size, MPI_DOUBLE, process, BUFFER_TAG, MPI_COMM_WORLD);
                }
            }
        } else {
            // get buffer from owner row
            comm_start = timer();

            MPI_Recv(buffer, buffer_size, MPI_DOUBLE, row_owner, BUFFER_TAG, MPI_COMM_WORLD, &status);

            comm_end = timer();
            total_comm += comm_end - comm_start;
        }

        // now do elimination on the rows below only using the relevant processes
        for (row = pivot+1; row < n; row++) {
            row_owner = row / bf;
            if (rank == row_owner) {
                denominator = buffer[0];

                target_row = row % bf;

                if (denominator) {
                    factor = local[target_row][pivot] / denominator;
                } else {
                    factor = 0.0; // prevent div by zero problems
                }

                for (column = pivot; column < n; column++) {
                    // pivot or pivot+1?
                    local[target_row][column] -= factor*buffer[column];
                }
            }
        }
    }

    // COMBINE RESULTS
    for (row = 0; row < n; row++) {
        if (row >= bf) {
            // not owned by root
            row_owner = row / bf;

            if (rank == ROOT_RANK) {
                // recieve results from row owner
                comm_start = timer();

                MPI_Recv(matrix[row], n, MPI_DOUBLE, row_owner, SUBMATRIX_DONE_TAG, MPI_COMM_WORLD, &status);

                comm_end = timer();
                total_comm += comm_end;
            } else if (rank == row_owner) {
                // send results to root
                target_row = row % bf;
                MPI_Send(local[target_row], n, MPI_DOUBLE, ROOT_RANK, SUBMATRIX_DONE_TAG, MPI_COMM_WORLD);
            }
        } else {
            if (rank == ROOT_RANK) {
                // root/P0 merges its own submatrix into the original
                target_row = row % bf;
                for (column = 0; column < n; column++) {
                    matrix[row][column] = local[target_row][column];
                }
            }
        }
    }


    // REPORT RESULTS
    comp_end = timer();
    total_comp = comp_end - comp_start;

    free(buffer);
    free_matrix(local, n);

    return matrix;
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


