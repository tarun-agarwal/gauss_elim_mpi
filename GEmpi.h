/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "tools.h"
#include <mpi.h>

#define P2P 0
#define BCAST 1

#define CONTINUOUS 0 // Process i gets rows blocks of rows
#define CICRULAR 1 // Process i gets every bf row

void test_parallel(int, int, int, int, char**);

