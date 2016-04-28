#ifndef LANCZOS_H_INCLUDED
#define LANCZOS_H_INCLUDED
/*include guards */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

typedef struct _coord{
     int row;
     int col;
}coord;

void matvec(coord* A, int* diagonal, int* scanned, int* nodeIndex,double* x, double* result,size_t n, size_t dim);
double dot(double* x, double* y, size_t n);
double norm(double* x, size_t n);
int eig(double* d, double* e, double* z,size_t n);
void scan(void* base, size_t n,size_t l, void(*oper)(void* x1, void* x2));
void rec_scan(void* base, size_t n,size_t l, void(*oper)(void* x1, void* x2));
void addInt(void* x1, void* x2);
void addDouble(void* x1, void* x2);
void addVec3(void* x1, void* x2);
void coord2Dot(coord* connections, int n, int process);
void node2Dot(int* myNodes, int n, int process);
void quicksort (void* base, size_t num, size_t size,
                    int (*comp)(const void*, const void*));
int coordCompare (const void* a, const void* b);

#endif
