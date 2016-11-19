#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "tools.h"

int *create_vector(int size)
{
  int *v = (int*) malloc (sizeof(int) * size);

  return v;
}

void inicialize_vector(int *vector, int size)
{
  int i = 0;

  srand((unsigned)time(NULL));

  for (i = 0; i < size; i++)
  {
    vector[i] = rand() % 100;
  }
}

void print_vector(int *vector, int size)
{
  int i = 0;

  for (i = 0; i < size; i++)
  {
    printf("%d ", vector[i]);
  }
  printf("\n");
}

void destroy_vector(int *vector)
{
  free(vector);
  vector = NULL;
}
