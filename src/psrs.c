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
#include <limits.h>

#include "quick_sort.h"
#include "tools.h"

int* FindSamples(int* Elements,
                 const int NumberOfElements,
                 const int NumberOfProcs,
                 int* NumberOfSamples);

int* FindPivots(int* Samples,
                int NumberOfProcs,
                int NumberOfSamples,
                int* NumberOfPivots);

void MakeSublists(int* Pivots,
                  int NumberOfPivots,
                  int* Elements,
                  int NumberOfElements,
                  int NumberOfProcs,
                  int** SublistSizes,
                  int*** Sublists);

int* MergeSublists(int* SublistSizes,
                   int** Sublists,
                   int NumberOfProcs,
                   int NumberOfElements);

int main(int argc, char* argv[]) {
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
    int* Elements = create_vector(NumberOfElements);
    initialize_vector(Elements, NumberOfElements, MAX_ELEMENT);

    // Phases 1
    int NumberOfSamples;
    int* Samples = FindSamples(Elements, NumberOfElements, NumberOfProcs,
                               &NumberOfSamples);

    // Phase 2
    int NumberOfPivots;
    int* Pivots =
        FindPivots(Samples, NumberOfProcs, NumberOfSamples, &NumberOfPivots);

    // Phase 3
    int** Sublists;
    int* SublistSizes;
    MakeSublists(Pivots, NumberOfPivots, Elements, NumberOfElements,
                 NumberOfProcs, &SublistSizes, &Sublists);

    // PHASE 4
    int* Final = MergeSublists(SublistSizes, Sublists, NumberOfProcs, NumberOfElements);

    print_vector(Final, NumberOfElements);
    destroy_vector(&Elements);
}

int* FindSamples(int* Elements,
                 const int NumberOfElements,
                 const int NumberOfProcs,
                 int* NumberOfSamples) {
    const int ProcsSqr = (NumberOfProcs * NumberOfProcs);

    *NumberOfSamples = ProcsSqr;
    int* Samples = (int*)malloc(*NumberOfSamples * sizeof(*Samples));

#pragma omp parallel num_threads(NumberOfProcs)
    {
        const int ProcessNumber = omp_get_thread_num();
        const int Start = (ProcessNumber)*NumberOfElements / NumberOfProcs;
        const int End = (ProcessNumber + 1) * NumberOfElements / NumberOfProcs;
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

int* FindPivots(int* Samples,
                int NumberOfProcs,
                int NumberOfSamples,
                int* NumberOfPivots) {
    quick_sort(Samples, 0, NumberOfSamples - 1);

    *NumberOfPivots = NumberOfProcs - 1;
    int* Pivots = (int*)malloc(*NumberOfPivots * sizeof(*Pivots));
    for (int iPivot = 0; iPivot < *NumberOfPivots; ++iPivot) {
        const int iSample =
            (iPivot + 1) * NumberOfProcs + (NumberOfProcs / 2) - 1;
        Pivots[iPivot] = Samples[iSample];
    }

    return Pivots;
}

void MakeSublists(int* Pivots,
                  int NumberOfPivots,
                  int* Elements,
                  int NumberOfElements,
                  int NumberOfProcs,
                  int** SublistSizes,
                  int*** Sublists) {
    MPI_Init(NULL, NULL);
    MPI_Comm Communicator;

    char* Args[2];
    char buf[20];
    Args[0] = buf;
    sprintf(Args[0], "%d", NumberOfElements);
    Args[1] = (char*)0;
    // printf("%s : %s", Args[0], Args[1]);
    MPI_Comm_spawn("./worker", Args, NumberOfProcs, MPI_INFO_NULL, 0,
                   MPI_COMM_SELF, &Communicator, MPI_ERRCODES_IGNORE);

    for (int iProc = 0; iProc < NumberOfProcs; iProc++) {
        const int Start = (iProc)*NumberOfElements / NumberOfProcs;
        const int End = (iProc + 1) * NumberOfElements / NumberOfProcs;
        const int Size = End - Start;
        MPI_Send(Elements + Start, Size, MPI_INT, iProc, 0, Communicator);
        MPI_Send(Pivots, NumberOfPivots, MPI_INT, iProc, 0, Communicator);
    }

    const int NumberOfSublists = NumberOfProcs * NumberOfProcs;
    *Sublists = (int**)malloc(NumberOfSublists * sizeof(**Sublists));
    *SublistSizes = (int*)malloc(NumberOfSublists * sizeof(**SublistSizes));

    for (int iProc = 0; iProc < NumberOfProcs; iProc++) {
        for (int iLocal = 0; iLocal < NumberOfProcs; ++iLocal) {
            const int iSublist = iProc * NumberOfProcs + iLocal;

            int* SublistSize = (*SublistSizes + iSublist);
            MPI_Recv(SublistSize, 1, MPI_INT, iProc, MPI_ANY_TAG, Communicator,
                     NULL);

            (*Sublists)[iSublist] =
                (int*)malloc(*SublistSize * sizeof(***Sublists));
            MPI_Recv((*Sublists)[iSublist], *SublistSize, MPI_INT, iProc,
                     MPI_ANY_TAG, Communicator, NULL);
        }
    }

    MPI_Finalize();
}

int* MergeSublists(int* SublistSizes,
                   int** Sublists,
                   int NumberOfProcs,
                   int NumberOfElements) {
    const int BlockSize = NumberOfProcs;
    const int NumberOfSublists = NumberOfProcs * NumberOfProcs;

    int* Elements = (int*)malloc(NumberOfElements * sizeof(*Elements));
    int* Starts = (int*)malloc(NumberOfProcs * sizeof(*Starts));
    Starts[0] = 0;
    for (int iProc = 1; iProc < NumberOfProcs; ++iProc) {
        int Sum = 0;
        for (int iSublist = 0; iSublist < BlockSize; ++iSublist) {
            Sum += SublistSizes[(iProc-1) * BlockSize + iSublist];
        }
        Starts[iProc] = Starts[iProc-1] + Sum;
    }

    int* Iterators = calloc(NumberOfSublists, sizeof(*Iterators));
#pragma omp parallel num_threads(NumberOfProcs)
    {
        const int ProcNum = omp_get_thread_num();
        const int StartSublist = ProcNum * BlockSize;
        int Size;
        if (ProcNum == NumberOfProcs-1) {
            Size =  NumberOfElements - Starts[ProcNum];
        } else {
            Size = Starts[ProcNum+1] - Starts[ProcNum];
        }
        
        for(int iElement = 0; iElement < Size; ++iElement) {
            int MinimumSublist = StartSublist;
            while(Iterators[MinimumSublist] == SublistSizes[MinimumSublist]) {
                MinimumSublist++;
            }
            for(int iLocal = 0; iLocal < BlockSize; ++iLocal) {
                const int iSublist = StartSublist + iLocal;
                if(Iterators[iSublist] < SublistSizes[iSublist]) {
                    const int ThisNumber = Sublists[iSublist][Iterators[iSublist]];
                    const int CurrentMinimum = Sublists[MinimumSublist][Iterators[MinimumSublist]];
                    if(ThisNumber < CurrentMinimum) {
                        MinimumSublist = iSublist;
                    }
                }
            }

            const int Value = Sublists[MinimumSublist][Iterators[MinimumSublist]];
            Elements[Starts[ProcNum]+iElement] = Value;
            Iterators[MinimumSublist]++;
        }
    }

    return Elements;
}

// void Merge(int* A, int ASize, int *B, int BSize, int *Out) {
//     const int FullSize = ASize + BSize;
//     int iA = 0, iB = 0;
//     for(int iElement = 0; iElement < FullSize; ++iElement) {
//         if(iA == ASize) {
//             Out[iElement] = B[iB];
//             iB++;
//         } else if(iB == ASize) {
//             Out[iElement] = A[iA];
//             iA++;
//         } else {
//             if (A[iA] < B[iB]) {
//                 Out[iElement] = A[iA];
//                 iA++;
//             } else {
//                 Out[iElement] = B[iB];
//                 iB++;
//             }
//         }
//     }
// }