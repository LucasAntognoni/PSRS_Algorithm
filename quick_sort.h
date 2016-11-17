#ifndef QUICK_SORT_H
#define QUICK_SORT_H

void swap(int *array, int i, int j);

int random_pivot(int *array, int p, int r);

int median(int array[]);

int median_pivot(int *array, int p, int r);

int middle_pivot(int *array, int p, int r);

int last_pivot_partition(int *array, int p, int r);

int first_pivot_partition(int *array, int p, int r);

void quick_sort(int *array, int p, int r);

#endif
