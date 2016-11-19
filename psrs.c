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

int main(int argc, char *argv[])
{

  // Checking input arguments...
  if (argc < 3)
  {
    printf("Use: ./psrs <Numbers> <Processes/Threads>");
    return EXIT_FAILURE;
  }

  int n, p, max;
  int *v;

  /*
    Parsing arguments...

    n >> Number of elements to be sorted;

    p >> Number of processes and threads
  */

  n = atoi(argv[1]);
  p = atoi(argv[2]);

  v = create_vector(n);
  inicialize_vector(v, n);

  max = ceil(n / p);

  // PHASE 1

  // PHASE 2

  // PHASE 3

  // PHASE 4

  print_vector(v, n);
}
