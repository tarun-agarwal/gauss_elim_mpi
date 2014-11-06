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

double** gauss_elim_parallel_p2p_continuous_old(double** A, int n){
    int num_processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign this process's rank
    MPI_Status status; // for MPI_Recv() calls

    int blocking = n / num_processes; // blocking factor

    int process; // loop iterator for sending stuff to processes
    int pivot, row, column; // for gauss elim
    double denominator, factor;

    int submatrix_size; // for submatrix decomp
    int y_start, y_end;
    int rows_per_process;

    double* buffer = malloc(n * sizeof(double));
    double** submatrix;

    int ok[1]; // tells children when it's ok to move to next pivot

    printf("Num processes=%d, blocking=%d, submatrix_size=%d, n=%d, rank=%d\n", num_processes, blocking, submatrix_size, n, rank);

    for (pivot = 0; pivot < n; pivot++) {
        denominator = A[pivot][pivot];

        submatrix_size = (n - pivot) / blocking;

        if (rank == ROOT_RANK) {
            printf("ROOT creates buffer: ");
            // parent, make buffer
            for (column = 0; column < n; column++) {
                buffer[column] = A[pivot][column];
                printf("%+e ", buffer[column]);
            }
            printf("\n");

            // now send buffer to every other process
            for (process = ROOT_RANK; process < num_processes; process++) {
                // printf("ROOT sends buffer to process %d\n", process);
                MPI_Send(buffer, n, MPI_DOUBLE, process, BUFFER_TAG, MPI_COMM_WORLD);
            }

            // split matrix into local matrices using continuous decomp
            for (process = ROOT_RANK; process < num_processes; process++) {
                // make local copy of submatrix

                y_start = process + pivot;
                y_end = y_start + blocking;
                submatrix = copy_submatrix(A, n, y_start, y_end);

                // send copy to mapped process
                MPI_Send(submatrix, submatrix_size, MPI_DOUBLE, process, SUBMATRIX_INIT_TAG, MPI_COMM_WORLD);
                printf("ROOT sends copy from %d to %d (size %d) to process %d\n", y_start, y_end, submatrix_size, process);

                for (row = 0; row < blocking; row++){
                    printf("\tSubmatrix %d row %d: ", process, row);
                    for (column = 0; column < n; column++) {
                        printf("%+e ", submatrix[row][column]);
                    }
                    printf("\n");
                }
            }

        }


        // now do elimination for each processor
        for (process = ROOT_RANK; process < num_processes; process++) {
            if (rank == process) {
                // printf("\tPROCESS %d: ", process);
                // get buffer
                MPI_Recv(buffer, n, MPI_DOUBLE, ROOT_RANK, BUFFER_TAG, MPI_COMM_WORLD, &status);
                // printf("gets buffer: ");

                // for (column = 0; column < n; column++) {
                //     printf("%+e ", buffer[column]);
                // }
                // printf("\n\t");

                // get submatrix
                MPI_Recv(submatrix, submatrix_size, MPI_DOUBLE, ROOT_RANK, SUBMATRIX_INIT_TAG, MPI_COMM_WORLD, &status);
                // printf("PROCESS %d gets submatrix of size %d:\n", process, submatrix_size);
                // for (row = 0; row < blocking; row++) {
                //     printf("\tSubmatrix %d row %d: ", process, row);
                //     for (column = 0; column < n; column++) {
                //         printf("%+e ", submatrix[row][column]);
                //     }
                //     printf("\n");
                // }
                printf("Doing PROCESS %d:\n", process);

                // do elimination for process i
                if (denominator) {
                    // printf("\tBuffer[pivot=%d]=%+e\n\tdenom=%+e\n", pivot, buffer[pivot], denominator);
                    factor = buffer[pivot] / denominator;
                } else {
                    // prevent div by zero problems
                    factor = 0.0;
                }
                // printf("\tFactor=%+e\n", factor);

                for (row = 0; row < blocking; row++) {
                    printf("\tNew sm%d row %d: ", process, row);
                    for (column = 0; column < n; column++) {
                        submatrix[row][column] -= factor * buffer[column];
                        printf("%+e ", submatrix[row][column]);
                    }
                    printf("\n");
                }
                printf("Ok.");

                // return eliminated submatrix to root
                MPI_Send(submatrix, submatrix_size, MPI_DOUBLE, ROOT_RANK, SUBMATRIX_DONE_TAG, MPI_COMM_WORLD);
                printf(" Sent %d to root.\n", process);
            }
        }

        // printf("Almost done pivot %d...", pivot);
        if (rank == 0) {
            // combine all results
            for (process = ROOT_RANK; process < num_processes; process++) {
                printf("Waiting for process %d...", process);
                MPI_Recv(submatrix, submatrix_size, MPI_DOUBLE, process, SUBMATRIX_DONE_TAG, MPI_COMM_WORLD, &status);

                // replace the relevant lines of A with this submatrix
                printf("done. Replacing lines %d to %d...", y_start, y_end);
                replace_submatrix(A, submatrix, n, y_start, y_end);
                printf("done.\n");
                printf("Intermediate A:\n\n");
                print_matrix(A, n);

                MPI_Send(&ok, 1, MPI_INT, process, OK_TAG, MPI_COMM_WORLD);
            }
        }

        // wait from OK from parent
        MPI_Recv(&ok, 1, MPI_INT, ROOT_RANK, OK_TAG, MPI_COMM_WORLD, &status);


        printf("Done pivot %d for all process.\n", rank);
    }

    printf("Freeing mallocs... ");
    free(buffer);
    free_matrix(submatrix, n);
    printf("Freed. Exiting.\n");
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


