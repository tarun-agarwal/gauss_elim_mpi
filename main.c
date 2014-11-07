/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "tools.h"
#include "GESerial.h"
#include "GEmpi.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    // test_serial();
    // time_serial_all();

    int n = 1024;

    MPI_Init(&argc, &argv);

    test_parallel(n, P2P, CONTINUOUS);
    test_parallel(n, P2P, CIRCULAR);
    test_parallel(n, BCAST, CONTINUOUS);
    test_parallel(n, BCAST, CIRCULAR);

    MPI_Finalize();

    return 0;
}
