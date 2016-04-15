#include "lanczos.h"
#define ITER 100

int main(int argc,char* argv[]) {

  int i, j, row, col, N, dim;
  srand(time(NULL));

  char matrix[100];

  if (argc > 1)
    strcpy(matrix,argv[1]);
  else {
    fprintf(stderr, "error: no matrix path provided\n");
    exit(EXIT_FAILURE);
  }

  FILE* A = fopen(matrix,"r");

  if (A == NULL){
    fprintf(stderr, "error: matrix path invalid\n");
    exit(EXIT_FAILURE);
  }

  fscanf(A,"%d %d %d\n",&row,&col,&N);
  fclose(A);

  if (row != col){
    fprintf(stderr,"error: this is not a sqare matrix\n");
    exit(EXIT_FAILURE);
  }

  dim = row;

  double* V = (double*)malloc(sizeof(double)*ITER*ITER);
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
      exit(EXIT_FAILURE);
    }

    if(i < ITER-1)
      for(j=0;j<dim;j++)
        q[i+1][j] = z[j]/b[i];
  }

  eig(a,b,V,ITER);

  printf("eigenvalues: \n\t");
  for(i=0;i<ITER;i++)
    printf("%lf\n\t", a[i]);
  printf("\n");

  free(z);
  free(V);
  for(i=0;i<ITER;i++)
    free(q[i]);

  return 0;
}
