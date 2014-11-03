/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include <mpi.h>

int gauss_elim_parallel_p2p_continuous(double** A, int n);

int gauss_elim_parallel_p2p_circular(double** A, int n);

int gauss_elim_parallel_broadcast_continuous(double** A, int n);

int gauss_elim_parallel_broadcast_circular(double** A, int n);
