#include "lanczos.h"

void matvec(char* matrix, double* x, double* result, size_t dim) {

  unsigned int row, col, i, n;
  double* diagonal = (double*)malloc(sizeof(double)*dim);

#pragma omp parallel for
  for(i=0;i<dim;i++) {
    result[i] = 0;
    diagonal[i] = 0;
  }

  FILE* A = fopen(matrix,"r");

  fscanf(A,"%d %d %d\n",&row,&col,&n);

  if(col != dim){
    fprintf(stderr, "error: cannot multiply %ux%u matrix by %lux1 vector\n"
                                                      , row,col,dim);
    exit(EXIT_FAILURE);
  }

//#pragma omp parallel for -- fix reading
  for(i=0;i<n;i++) {

    fscanf(A,"%d %d\n",&row,&col);

    if (row != col){
      result[row-1]-=x[col-1];
      result[col-1]-=x[row-1];
      diagonal[col-1]+=1;
      diagonal[row-1]+=1;
    }


    // for(j=0;j<dim;j++)
    //   printf("%lf\n", result[j]);

  }

#pragma omp parallel for
  for(i=0;i<dim;i++)
    result[i]+=diagonal[i]*x[i];


  free(diagonal);
  fclose(A);
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
