/*

=======================================================================
  Parallel Sorting by Regular Sampling

  Quicksort Algorithm

  Quinn, Michael J. Parallel Programming in C with MPI and OpenMP. 2003
=======================================================================

*/

#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "quick_sort.h"
#include "tools.h"

int *FindSamples(int *Elements, const int NumberOfElements,
                 const int NumberOfProcs, int *NumberOfSamples);
int *FindPivots(int **Samples, int NumberOfProcs, int NumberOfSamples, int* NumberOfPivots);

int main(int argc, char *argv[]) {
  // Checking input arguments
  if (argc < 3) {
    printf("Use: %s <Number of Elements> <Number of Processes>\n", argv[0]);
    return EXIT_FAILURE;
  }

  int NumberOfElements = atoi(argv[1]);
  int NumberOfProcs = atoi(argv[2]);
  const int MAX_ELEMENT = 1000;

  //
  // IMPORTANTE:
  //
  int *Elements = create_vector(NumberOfElements);
  initialize_vector(Elements, NumberOfElements, MAX_ELEMENT);

  // Phases 1
  int NumberOfSamples;
  int *Samples =
      FindSamples(Elements, NumberOfElements, NumberOfProcs, &NumberOfSamples);

  // Phase 2
  int NumberOfPivots;
  int* Pivots = FindPivots(&Samples, NumberOfProcs, NumberOfSamples, &NumberOfPivots);

  // Phase 3

  // PHASE 4

  print_vector(Elements, NumberOfElements);
  destroy_vector(&Elements);
}

int *FindSamples(int *Elements, const int NumberOfElements,
                 const int NumberOfProcs, int *NumberOfSamples) {
  const int ProcsSqr = (NumberOfProcs * NumberOfProcs);

  *NumberOfSamples = ProcsSqr;
  int *Samples = (int *)malloc(*NumberOfSamples * sizeof(*Samples));

#pragma omp parallel num_threads(NumberOfProcs)
  {
    const int ProcessNumber = omp_get_thread_num();
    const int Start = (ProcessNumber)*NumberOfElements / NumberOfProcs;
    const int End = (ProcessNumber+1) * NumberOfElements / NumberOfProcs;
    quick_sort(Elements, Start, End - 1);

    for (int iLocal = 0; iLocal < NumberOfProcs; iLocal++) {
      const int iElement =
          Start + (int)(iLocal * (NumberOfElements / (float)ProcsSqr));
      const int iSample = (ProcessNumber * NumberOfProcs) + iLocal;
      Samples[iSample] = Elements[iElement];
    }
  }

  return Samples;
}

int *FindPivots(int **Samples, int NumberOfProcs, int NumberOfSamples, int* NumberOfPivots) {

  quick_sort(*Samples, 0, NumberOfSamples - 1);

  *NumberOfPivots = NumberOfProcs - 1;
  int *Pivots = (int *)malloc(*NumberOfPivots * sizeof(*Pivots));
  for (int iPivot = 0; iPivot < *NumberOfPivots; ++iPivot) {
    const int iSample = (iPivot+1) * NumberOfProcs + (NumberOfProcs / 2) - 1;
    Pivots[iPivot] = *Samples[iSample]; 
  }

  free(*Samples);
  *Samples = NULL;

  return Pivots;
}