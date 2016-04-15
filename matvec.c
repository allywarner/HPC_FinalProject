#include "lanczos.h"

void matvec(int** A, double* x, double* result, size_t n, size_t dim) {

  int row, col;
  unsigned int i;
  double* diagonal = (double*)malloc(sizeof(double)*dim);

#pragma omp parallel for
  for(i=0;i<dim;i++) {
    result[i] = 0;
    diagonal[i] = 0;
  }

//#pragma omp parallel for -- fix reading
  for(i=0;i<n;i++) {
    row = A[0][i];
    col = A[1][i];
    if (row != col){
      result[row-1]-=x[col-1];
      result[col-1]-=x[row-1];
      diagonal[col-1]+=1;
      diagonal[row-1]+=1;
    }
  }

#pragma omp parallel for
  for(i=0;i<dim;i++)
    result[i]+=diagonal[i]*x[i];

  free(diagonal);
}

double dot(double* x, double* y, size_t n){
  unsigned int i;
  double a=0.0;
#pragma omp parallel for
  for(i=0;i<n;i++)
    a += x[i]*y[i];
  return a;
}

double norm(double* x, size_t n){
  unsigned int i;
  double nor = 0.0;
#pragma omp parallel for
  for(i=0;i<n;i++)
    nor+=x[i]*x[i];
  nor = sqrt(nor);
  return nor;
}
