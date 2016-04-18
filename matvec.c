#include "lanczos.h"

void matvec(coord* A, int* diagonal, double* x, double* result, size_t n, size_t dim) {

unsigned int i;
int j;

int* scanned = (int*)malloc(sizeof(int)*dim);

#pragma omp parallel for
for(i=0;i<dim;i++){
  result[i]=diagonal[i]*x[i];
  scanned[i] = diagonal[i];
}

scan(scanned,dim,sizeof(int),addInt);

// multiply x by laplacian of A
#pragma omp parallel for private(j)
  for(i=0;i<dim;i++) {
    if(i > 0)
      for(j=0;j<diagonal[i];j++){
        result[i]-=x[A[scanned[i-1] + j].col-1];
      }
    else
      for(j=0;j<diagonal[i];j++) {
        result[i]-=x[A[j].col-1];
      }
    }

    free(scanned);
    }


double dot(double* x, double* y, size_t n){
  unsigned int i;
  double a=0.0;
#pragma omp parallel for reduction(+: a)
  for(i=0;i<n;i++)
    a += x[i]*y[i];
  return a;
}

double norm(double* x, size_t n){
  unsigned int i;
  double nor = 0.0;
#pragma omp parallel for reduction(+: nor)
  for(i=0;i<n;i++)
    nor+=x[i]*x[i];
  nor = sqrt(nor);
  return nor;
}
