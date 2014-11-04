/** Ryan Ordille, 260399372, ECSE420lab2 **/

#include "tools.h"
#include "GESerial.h"
#include "GEmpi.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    // test_serial();
    // time_serial_all();

    test_parallel(argc, argv);

    return 0;
}
