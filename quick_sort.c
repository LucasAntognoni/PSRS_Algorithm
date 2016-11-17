#include <stdio.h>
#include <stdlib.h>

#include "quick_sort.h"

void swap(int *array, int i, int j)
{
  int aux;

  aux = array[i];
  array[i] = array[j];
  array[j] = aux;
}

int partition(int *array, int p, int r)
{
  int x, i, j;

  x = array[p];
  i = p;

  for(j = p + 1; j <= r; j++)
  {
    if (array[j] <= x)
    {
      i = i + 1;
      swap(array, i, j);
    }
  }
  swap(array, i, p);

  return i;
}

void quick_sort(int *array, int p, int r)
{
  if (p < r)
  {
    int q;
    
    q = partition(array, p, r);
    quick_sort(array, p, q - 1);
    quick_sort(array, q + 1, r);
  }
}
