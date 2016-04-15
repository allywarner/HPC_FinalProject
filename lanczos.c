#include "lanczos.h"
#define ITER 100

int main(int argc,char* argv[]) {

  int i, j, row, col, N, dim;
  srand(time(NULL));  // seed random number generator
  char mstr[100];

  // read name of matrix from arg list
  if (argc > 1)
    strcpy(mstr,argv[1]);
  else {
    fprintf(stderr, "error: no matrix path provided\n");
    exit(EXIT_FAILURE);
  }

  // open matrix file
  FILE* matrix = fopen(mstr,"r");
  if (matrix == NULL){
    fprintf(stderr, "error: matrix path invalid\n");
    exit(EXIT_FAILURE);
  }

  //scan first row of matrix for metadata
  fscanf(matrix,"%d %d %d\n",&row,&col,&N);
  if (row != col){
    fprintf(stderr,"error: this is not a sqare matrix\n");
    exit(EXIT_FAILURE);
  }
  dim = row;

  // define variables
  int* A[2]; //graph matrix coordinates
  A[0] = (int*)malloc(sizeof(int)*N);
  A[1] = (int*)malloc(sizeof(int)*N);
  double* V = (double*)malloc(sizeof(double)*ITER*ITER);// blank workspace for eig function
  double* q[ITER];          //Krylov space basis matrix
  double* z = (double*)malloc(sizeof(double)*dim); //new vector in each lanczos iteration
  double a[ITER]; // diagonal of lanczos matrix
  double b[ITER]; // subdiagonal of lanczos matrix

  //scan matrix coordinates and store in memory
  for(i=0;i<N;i++){
    fscanf(matrix,"%d %d\n",&A[0][i],&A[1][i]);
  }

  // allocate space for each Krylov space vector
  for(i=0;i<ITER;i++)
    q[i] = (double*)malloc(sizeof(double)*dim);

  // initial guess
  for(i=0;i<dim;i++)
    q[0][i] = (double)rand()/(double)RAND_MAX * 100;

  // normalize initial guess
  b[0] = norm(q[0],dim);
  for(i=0;i<dim;i++)
    q[0][i] = q[0][i]/b[0];

  //lanczos iterations
  for(i=0;i<ITER;i++) {
    matvec(A,q[i],z,N,dim);
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

  fclose(matrix);

  // compute eigenvalues and eigenvectors of lanczos matrix
  eig(a,b,V,ITER);


  // print eigenvalues
  printf("eigenvalues: \n\t");
  for(i=0;i<ITER;i++)
    printf("%lf\n\t", a[i]);
  printf("\n");

  free(A[1]);
  free(A[2]);
  free(z);
  free(V);
  for(i=0;i<ITER;i++)
    free(q[i]);

  return 0;
}
