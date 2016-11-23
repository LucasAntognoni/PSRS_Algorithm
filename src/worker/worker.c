#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int NumberOfProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfProcs);
    int TotalNumberOfElements = atoi(argv[1]);

    int Rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    int NumberOfElements;
    {
        const int Start = (Rank)*TotalNumberOfElements / NumberOfProcs;
        const int End = (Rank + 1) * TotalNumberOfElements / NumberOfProcs;
        NumberOfElements = End - Start;
    }
    int* Elements = (int*)malloc(NumberOfElements * sizeof(*Elements));

    const int NumberOfPivots = NumberOfProcs - 1;
    int* Pivots = (int*)malloc(NumberOfPivots * sizeof(*Pivots));

    MPI_Comm Parent;
    MPI_Comm_get_parent(&Parent);
    MPI_Recv(Elements, NumberOfElements, MPI_INT, 0, MPI_ANY_TAG, Parent, NULL);
    MPI_Recv(Pivots, NumberOfPivots, MPI_INT, 0, MPI_ANY_TAG, Parent, NULL);

    int MyStart = 0, MyEnd = 0;
    int Start = 0;
    for (int iProc = 0; iProc < NumberOfProcs; ++iProc) {
        int iElement = Start;
        int End;
        if (iProc != NumberOfPivots) {
            while (Elements[iElement] < Pivots[iProc] &&
                   iElement < NumberOfElements) {
                iElement++;
            }
            End = iElement;
        } else {
            End = NumberOfElements;
        }

        const int NumberOfElementsToSend = End - Start;
        if (iProc != Rank) {
            MPI_Send(&NumberOfElementsToSend, 1, MPI_INT, iProc, 0,
                     MPI_COMM_WORLD);
            MPI_Send(Elements + Start, NumberOfElementsToSend, MPI_INT, iProc,
                     0, MPI_COMM_WORLD);
        } else {
            MyStart = Start;
            MyEnd = End;
        }
        Start = End;
    }

    for (int iProc = 0; iProc < NumberOfProcs; ++iProc) {
        int NumberOfElementsInSublist;
        int* Sublist;
        if (iProc != Rank) {
            MPI_Recv(&NumberOfElementsInSublist, 1, MPI_INT, iProc, MPI_ANY_TAG,
                     MPI_COMM_WORLD, NULL);

            Sublist =
                (int*)malloc(NumberOfElementsInSublist * sizeof(*Sublist));
            MPI_Recv(Sublist, NumberOfElementsInSublist, MPI_INT, iProc,
                     MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
        } else {
            NumberOfElementsInSublist = MyEnd - MyStart;
            Sublist = Elements + MyStart;
        }

        MPI_Send(&NumberOfElementsInSublist, 1, MPI_INT, 0, 0, Parent);
        MPI_Send(Sublist, NumberOfElementsInSublist, MPI_INT, 0, 0, Parent);
        if (iProc != Rank) {
            free(Sublist);
        }
    }

    free(Pivots);
    free(Elements);


    MPI_Finalize();
    return EXIT_SUCCESS;
}