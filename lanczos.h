#ifndef LANCZOS_H_INCLUDED
#define LANCZOS_H_INCLUDED
/*include guards */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

void matvec(char* matrix, double* x, double* result,size_t dim);
double dot(double* x, double* y, size_t n);
double norm(double* x, size_t n);
int eig(double* d, double* e, double* z,size_t n);

#endif
