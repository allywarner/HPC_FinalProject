#include "lanczos.h"
#define ITER 100

int main(int argc,char* argv[]) {

  int i, j, row, col, N, dim;
  srand(time(NULL));

  char matrix[100];

  strcpy(matrix,"Matrices/finance256.dat");

  FILE* A = fopen(matrix,"r");

  fscanf(A,"%d %d %d\n",&row,&col,&N);
  fclose(A);

  if (row != col){
    printf("error: this is not a sqare matrix\n");
    exit(EXIT_FAILURE);
  }

  dim = row;

  double* q[ITER];
  double* z = (double*)malloc(sizeof(double)*dim);
  double a[ITER];
  double b[ITER];

  for(i=0;i<ITER;i++)
  q[i] = (double*)malloc(sizeof(double)*dim);

  for(i=0;i<dim;i++)
  q[0][i] = (double)rand()/(double)RAND_MAX * 100;

  b[0] = norm(q[0],dim);

  for(i=0;i<dim;i++)
  q[0][i] = q[0][i]/b[0];

  //lanczos iterations
  for(i=0;i<ITER;i++) {
    matvec(matrix,q[i],z,dim);
    a[i] = dot(q[i],z,dim);
    if(i > 0) {
      for(j=0;j<dim;j++)
        z[j] = z[j] - a[i]*q[i][j] - b[i-1]*q[i-1][j]; //no reorthogonalization
    } else {
        for(j=0;j<dim;j++)
          z[j] = z[j] - a[i]*q[i][j];
    }

    b[i] = norm(z,dim);

    if (b[i] == 0){
      printf("stopped short of %d iterations\n", ITER);
      break;
    }

    if(i < ITER-1)
      for(j=0;j<dim;j++)
        q[i+1][j] = z[j]/b[i];
  }

  free(z);
  for(i=0;i<ITER;i++)
    free(q[i]);

  return 0;
}
