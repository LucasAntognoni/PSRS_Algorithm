#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "tools.h"

int* create_vector(int size) {
    int* v = (int*)malloc(sizeof(int) * size);

    return v;
}

void initialize_vector(int* vector, int size, int max) {
    int i = 0;

    srand((unsigned)time(NULL));

    for (i = 0; i < size; i++) {
        vector[i] = 1 + (rand() % max);
    }
}

void print_vector(int* vector, int size) {
    int i = 0;

    char* separator = "";
    for (i = 0; i < size; i++) {
        printf("%s%d", separator, vector[i]);
        separator = ", ";
    }
    printf(".");
    fflush(stdout);
}

void destroy_vector(int** vector) {
    if (vector && *vector) {
        free(*vector);
        vector = NULL;
    }
}
