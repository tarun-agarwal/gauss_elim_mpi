/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "GEmpi.h"

#include <stdio.h> // for printf
#include <stdlib.h> // for malloc/free

#define BUFFER_TAG 0
#define SUBMATRIX_INIT_TAG 1
#define SUBMATRIX_DONE_TAG 2
#define OK_TAG 3
#define COMM_TAG 4
#define COMP_TAG 5

#define ROOT_RANK 0

double** gauss_elim_parallel_p2p(double** matrix, int n, int mode) {
    // DECLARATIONS & INITALIZATIONS

    // for MPI calls:
    int np, rank; // from MPI_Init()
    MPI_Status status; // for MPI_Recv()

    MPI_Comm_size(MPI_COMM_WORLD, &np); // number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // this process' rank

    int target_row; // target rank for send/recv
    int row_owner, pivot_owner;

    // for matrix allocation
    int bf = n / np; // blocking factor, i.e. number of rows per local matrix
    int local_size = n * bf; // # of elements per MPI_Send of submatrix
    int buffer_size;

    double* buffer; // allocate row buffer
    double** local; // local submatrix for each process

    // for timer:
    double comm_start, comm_end;
    double comp_start, comp_end;
    double total_comm = 0.0;
    double total_comm_all = 0.0; // for all processes
    double total_comp = 0.0;
    double total_comp_all = 0.0;

    // loop iterators:
    int process; // looping over processes
    int pivot, row, column; // looping through matrices

    // for gauss elimination:
    double denominator, factor;

    comp_start = timer();

    // ALLOCATE LOCAL SUBMATRICES
    for (process = 0; process < np; process++)
        if (rank == process)
            // local submatrix should be of size n*bf*(8 bytes)
            local = allocate_submatrix(n, bf);

    // INIT LOCAL SUBMATRICES
    for (row = 0; row < n; row++) {
        if (mode == CONTINUOUS) {
            row_owner = row / bf;
            target_row = row % bf;
        }
        else {
            row_owner = row % np;
            target_row = row / np ;
        }

        if (row >= bf) {
            // not owned by root, so send to relevant process
            if (rank == ROOT_RANK) {
                // send to applicable process
                MPI_Send(matrix[row], n, MPI_DOUBLE, row_owner, SUBMATRIX_INIT_TAG, MPI_COMM_WORLD);
            } else if (rank == row_owner) {
                // recv submatrix from root

                comm_start = timer();

                MPI_Recv(local[target_row], n, MPI_DOUBLE, ROOT_RANK, SUBMATRIX_INIT_TAG, MPI_COMM_WORLD, &status);

                comm_end = timer();
                total_comm += comm_end - comm_start;
            }
        } else {
            // first bf rows so allocate root's portion directly
            if (rank == ROOT_RANK) {
                for (column = 0; column < n; column++) {
                    local[target_row][column] = matrix[row][column];
                }
            }
        }
    }

    // PERFORM GAUSS ELIMINATION
    for (pivot = 0; pivot < n; pivot++) {
        // allocate buffer
        buffer_size = n - pivot;
        buffer = malloc(buffer_size * sizeof(double));

        if (mode == CONTINUOUS) {
            pivot_owner = pivot / bf;
            target_row = pivot % bf;
        }
        else {
            pivot_owner = pivot % np;
            target_row = pivot / np;
        }

        // owner sends pivot row to everyone via buffer
        if (rank == pivot_owner) {
            // need to send this row to all other processes
            for (column = pivot; column < n; column++) {
                buffer[column - pivot] = local[target_row][column];
            }

            // printf("%d will now send buffer.\n", pivot_owner);
            for (process = 0; process < np; process++) {
                if (process != rank) { // don't send to self but to everyone else
                    MPI_Send(buffer, buffer_size, MPI_DOUBLE, process, BUFFER_TAG, MPI_COMM_WORLD);
                }
            }

        } else {
            // get buffer from owner row
            comm_start = timer();

            MPI_Recv(buffer, buffer_size, MPI_DOUBLE, pivot_owner, BUFFER_TAG, MPI_COMM_WORLD, &status);

            comm_end = timer();
            total_comm += comm_end - comm_start;
        }

        // now do elimination on the rows below only using the relevant processes
        for (row = pivot+1; row < n; row++) {
            if (mode == CONTINUOUS) {
                row_owner = row / bf;
                target_row = row % bf;
            }
            else {
                row_owner = row % np;
                target_row = row / np;
            }

            if (rank == row_owner) {
                denominator = buffer[0];

                factor = local[target_row][pivot] / denominator;

                for (column = pivot; column < n; column++)
                    local[target_row][column] -= factor*buffer[column - pivot];
            }
        }

        free(buffer);
    }

    // COMBINE RESULTS
    for (row = 0; row < n; row++) {
        if (mode == CONTINUOUS) {
            row_owner = row / bf;
            target_row = row % bf;
        } else {
            row_owner = row % np;
            target_row = row / np;
        }

        if (row_owner != ROOT_RANK) {
            // not owned by root
            if (rank == ROOT_RANK) {
                // recieve results from row owner
                comm_start = timer();

                MPI_Recv(matrix[row], n, MPI_DOUBLE, row_owner, SUBMATRIX_DONE_TAG, MPI_COMM_WORLD, &status);

                comm_end = timer();
                total_comm += comm_end - comm_start;

            } else if (rank == row_owner) {
                // send results to root
                MPI_Send(local[target_row], n, MPI_DOUBLE, ROOT_RANK, SUBMATRIX_DONE_TAG, MPI_COMM_WORLD);
            }
        } else {
            if (rank == ROOT_RANK)
                // root/P0 merges its own submatrix into the original
                for (column = 0; column < n; column++)
                    matrix[row][column] = local[target_row][column];
        }
    }

    comp_end = timer();
    total_comp = comp_end - comp_start;

    // REPORT RESULTS
    // send comp results to the root
    MPI_Reduce(&total_comm, &total_comm_all, 1, MPI_DOUBLE, MPI_SUM, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Reduce(&total_comp, &total_comp_all, 1, MPI_DOUBLE, MPI_SUM, ROOT_RANK, MPI_COMM_WORLD);

    if (rank == ROOT_RANK) {
        printf(ANSI_COLOR_RED "Point-to-Point, ");
        if (mode == CONTINUOUS)
            printf("continuous decomposition.");
        else
            printf("circular decomposition.");
        printf("\n" ANSI_COLOR_RESET);
        printf(ANSI_COLOR_BLUE "Matrix size (%d,%d), %d processes, blocking factor = %d\n" ANSI_COLOR_RESET, n, n, np, bf);
        printf("Communication time: %+e s\n", total_comm_all);
        printf("Computation time: %+e s\n", total_comp_all - total_comm_all);
        printf("Total computation + communication time: %+e s\n\n", total_comp_all);
    }

    free(local);

    return matrix;
}

double** gauss_elim_parallel_broadcast(double** matrix, int n, int mode){
    // DECLARATIONS & INITALIZATIONS

    // for MPI calls:
    int np, rank; // from MPI_Init()

    MPI_Comm_size(MPI_COMM_WORLD, &np); // number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // this process' rank

    int target_row; // target rank for send/recv
    int row_owner, pivot_owner;

    // for matrix allocation
    int bf = n / np; // blocking factor, i.e. number of rows per local matrix
    int local_size = n * bf; // # of elements per MPI_Send of submatrix
    int buffer_size;

    double* buffer; // allocate row buffer
    double** local; // local submatrix for each process
    double* local_temp = malloc(n * sizeof(double));

    // for timer:
    double comm_start, comm_end;
    double comp_start, comp_end;
    double total_comm = 0.0;
    double total_comm_all = 0.0; // for all processes
    double total_comp = 0.0;
    double total_comp_all = 0.0;

    // loop iterators:
    int process; // looping over processes
    int pivot, row, column; // looping through matrices

    // for gauss elimination:
    double denominator, factor;

    comp_start = timer();

    // ALLOCATE LOCAL SUBMATRICES
    for (process = 0; process < np; process++)
        if (rank == process)
            // local submatrix should be of size n*bf*(8 bytes)
            local = allocate_submatrix(n, bf);

    MPI_Barrier(MPI_COMM_WORLD);

    // INIT LOCAL SUBMATRICES
    for (row = 0; row < n; row++) {
        if (mode == CONTINUOUS) {
            row_owner = row / bf;
            target_row = row % bf;
        }
        else {
            row_owner = row % np;
            target_row = row / np ;
        }

        if (row >= bf) {
            // not owned by root, so send to relevant process
            // send to applicable process

            local_temp = matrix[row];

            comm_start = timer();
            MPI_Bcast(local_temp, n, MPI_DOUBLE, row_owner, MPI_COMM_WORLD);
            comm_end = timer();
            total_comm += comm_end - comm_start;

            if (rank == row_owner) {
                // recv submatrix from root
                local[target_row] = local_temp;
            }
        } else {
            // first bf rows so allocate root's portion directly
            if (rank == ROOT_RANK) {
                for (column = 0; column < n; column++) {
                    local[target_row][column] = matrix[row][column];
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // PERFORM GAUSS ELIMINATION
    for (pivot = 0; pivot < n; pivot++) {
        // allocate buffer
        buffer_size = n - pivot;
        buffer = malloc(buffer_size * sizeof(double));

        if (mode == CONTINUOUS) {
            pivot_owner = pivot / bf;
            target_row = pivot % bf;
        } else {
            pivot_owner = pivot % np;
            target_row = pivot / np;
        }

        // owner sends pivot row to everyone via buffer
        if (rank == pivot_owner) {
            // make buffer
            for (column = pivot; column < n; column++)
                buffer[column - pivot] = local[target_row][column];
        }

        comm_start = timer();

        // send this row to all other processes
        MPI_Bcast(buffer, buffer_size, MPI_DOUBLE, pivot_owner, MPI_COMM_WORLD);

        comm_end = timer();
        total_comm += comm_end - comm_start;

        // now do elimination on the rows below only using the relevant processes
        for (row = pivot+1; row < n; row++) {
            if (mode == CONTINUOUS) {
                row_owner = row / bf;
                target_row = row % bf;
            }
            else {
                row_owner = row % np;
                target_row = row / np;
            }

            if (rank == row_owner) {
                denominator = buffer[0];

                factor = local[target_row][pivot] / denominator;

                for (column = pivot; column < n; column++)
                    local[target_row][column] -= factor*buffer[column - pivot];
            }
        }

        free(buffer);

        // don't let other processes get too far ahead, i.e. avoid race conditions
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // COMBINE RESULTS
    for (row = 0; row < n; row++) {
        if (mode == CONTINUOUS) {
            row_owner = row / bf;
            target_row = row % bf;
        } else {
            row_owner = row % np;
            target_row = row / np;
        }

        if (row_owner != ROOT_RANK) {
            local_temp = matrix[row];

            comm_start = timer();

            // send this row to everyone
            MPI_Bcast(local_temp, n, MPI_DOUBLE, row_owner, MPI_COMM_WORLD);

            comm_end = timer();
            total_comm += comm_end - comm_start;

            if (rank == ROOT_RANK) {
                // root merges results from row owner
                matrix[row] = local_temp;
            }

        } else {
            if (rank == ROOT_RANK)
                // root/P0 merges its own submatrix into the original
                for (column = 0; column < n; column++)
                    matrix[row][column] = local[target_row][column];
        }
    }

    comp_end = timer();
    total_comp = comp_end - comp_start;

    // REPORT RESULTS
    // send comp results to the root
    MPI_Reduce(&total_comm, &total_comm_all, 1, MPI_DOUBLE, MPI_SUM, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Reduce(&total_comp, &total_comp_all, 1, MPI_DOUBLE, MPI_SUM, ROOT_RANK, MPI_COMM_WORLD);

    if (rank == ROOT_RANK) {
        printf(ANSI_COLOR_RED "Collective, ");
        if (mode == CONTINUOUS)
            printf("continuous decomposition.");
        else
            printf("circular decomposition.");
        printf("\n" ANSI_COLOR_RESET);
        printf(ANSI_COLOR_BLUE "Matrix size (%d,%d), %d processes, blocking factor = %d\n" ANSI_COLOR_RESET, n, n, np, bf);
        printf("Communication time: %+e s\n", total_comm_all);
        printf("Computation time: %+e s\n", total_comp_all - total_comm_all);
        printf("Total computation + communication time: %+e s\n\n", total_comp_all);
    }

    free(local);

    return matrix;
}

void test_parallel(int n, int mechanism, int decomp_mode) {
    int rank;

    double start, end;

    double** A = make_matrix(n);


    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == ROOT_RANK){
        // print_matrix(A, n);
        // printf("\n\n");
    }

    if (mechanism == P2P)
        A = gauss_elim_parallel_p2p(A, n, decomp_mode);
    else
        A = gauss_elim_parallel_broadcast(A, n, decomp_mode);

    if (rank == ROOT_RANK){
        // print_matrix(A, n);
        // printf("\n");
    }


}


