/*

=======================================================================
  Parallel Sorting by Regular Sampling

  Quicksort Algorithm

  Quinn, Michael J. Parallel Programming in C with MPI and OpenMP. 2003
=======================================================================

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include "tools.h"
#include "quick_sort.h"

void Phase1(int* Elements, const int NumberOfElements, const int NumberOfProcs);

int main(int argc, char* argv[]) {
    // Checking input arguments
    if (argc < 3) {
        printf("Use: %s <Number of Elements> <Number of Processes>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int NumberOfElements = atoi(argv[1]);
    int NumberOfProcs = atoi(argv[2]);
    const int MAX_ELEMENT = 1000;

    int* Elements = create_vector(NumberOfElements);
    initialize_vector(Elements, NumberOfElements, MAX_ELEMENT);

    // PHASE 1
    Phase1(Elements, NumberOfElements, NumberOfProcs);

    // PHASE 2

    // PHASE 3

    // PHASE 4

    print_vector(Elements, NumberOfElements);
    destroy_vector(&Elements);
}

void Phase1(int* Elements, const int NumberOfElements, const int NumberOfProcs) {
    #pragma omp parallel num_threads(NumberOfProcs)
    {
        const int ProcessNumber = omp_get_thread_num();
        const int Start = (ProcessNumber) * NumberOfElements / NumberOfProcs;
        const int End = (ProcessNumber + 1) * NumberOfElements / NumberOfProcs;
        quick_sort(Elements, Start, End);
    }
}