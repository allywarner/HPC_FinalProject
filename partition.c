#include "lanczos.h"
#include <mpi.h>
#define ITER 100

void partition(coord* A, size_t dim, size_t N, int comm_size, MPI_Comm comm);

int main(int argc,char* argv[]) {

  MPI_Init(&argc,&argv);
  MPI_Comm world_comm = MPI_COMM_WORLD;
  int world_rank, world_size;
  MPI_Comm_rank(world_comm, &world_rank);
  MPI_Comm_size(world_comm, &world_size);

  int i, row, col, N, dim;
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

  coord* A = (coord*)malloc(sizeof(coord)*2*N); //ajacency matrix coordinates
  //scan matrix coordinates and store in memory
  for(i=0;i<N;i++)
    fscanf(matrix,"%d %d\n",&A[i].row,&A[i].col);

  fclose(matrix);

  partition(A, dim, N, world_size, world_comm);

}

void partition(coord* A, size_t dim, size_t N, int comm_size, MPI_Comm comm ) {

  unsigned int i, j, l, T_size=ITER;
  int* diagonal = (int*)malloc(sizeof(int)*dim); // #nonzeros specified per row
  int* scanned = (int*)malloc(sizeof(int)*dim);
  double* V = (double*)malloc(sizeof(double)*ITER*ITER);// blank workspace for eig function
  double* q[ITER];          //Krylov space basis matrix
  double* z = (double*)malloc(sizeof(double)*dim); //new vector in each lanczos iteration
  double a[ITER]; // diagonal of lanczos matrix
  double b[ITER]; // subdiagonal of lanczos matrix

  #pragma omp parallel for
  for(i=0;i<dim;i++)
    diagonal[i] = 0;

  for(i=0;i<N;i++)
    diagonal[A[i].row-1]+=1;

  #pragma omp parallel for
  for(i=0;i<dim;i++)
    scanned[i] = diagonal[i];

  scan(scanned,dim,sizeof(int),addInt);

  // allocate space for each Krylov space vector
  #pragma omp parallel for
  for(i=0;i<ITER;i++)
    q[i] = (double*)malloc(sizeof(double)*dim);

  // initial guess
  #pragma omp parallel for
  for(i=0;i<dim;i++)
    q[0][i] = (double)rand()/(double)RAND_MAX * 100;

  // normalize initial guess
  b[0] = norm(q[0],dim);
  #pragma omp parallel for
  for(i=0;i<dim;i++)
    q[0][i] = q[0][i]/b[0];

  //lanczos iterations
  double ortho[ITER];
  double begin = omp_get_wtime();
  for(j=0;j<ITER;j++) {
    matvec(A,diagonal,scanned,q[j],z,N,dim);
    a[j] = dot(q[j],z,dim);
    #pragma omp parallel for
    for(i=0;i<=j;i++)
      ortho[i] = dot(q[i],z,dim);

    if(j > 0) {
      for(l=0;l<=j;l++) {
        #pragma omp parallel for
        for(i=0;i<dim;i++)
          z[i] -= ortho[l]*q[l][i];
      }
      #pragma omp parallel for
      for(i=0;i<=j;i++)
        ortho[i] = dot(q[i],z,dim);
      for(l=0;l<=j;l++){
        #pragma omp parallel for
        for(i=0;i<dim;i++)
          z[i] -= ortho[l]*q[l][i];
      }

        //no reorthogonalization
          // z[i] = z[i] - a[j]*q[j][i] - b[j-1]*q[j-1][i];


      } else {
          #pragma omp for
          for(i=0;i<dim;i++)
            z[i] = z[i] - a[j]*q[j][i];
      }

    b[j] = norm(z,dim);
    if (b[j] < 1e-13){
      printf("stopped short of %d iterations\n", ITER);
      T_size = j+1;
      break;
    }
    if(j < ITER-1)
      // #pragma omp parallel for
      for(i=0;i<dim;i++)
        q[j+1][i] = z[i]/b[j];
  }
  double end = omp_get_wtime();

  printf("time: %lf\n", end-begin);

  // compute eigenvalues and eigenvectors of lanczos matrix
  eig(a,b,V,T_size);

  double* splitter = (double*)malloc(sizeof(double)*dim);

  #pragma omp parallel for
  for(i=0;i<dim;i++)
    splitter[i] = 0;

  #pragma omp parallel for private(j)
  for(i=0;i<dim;i++)
    for(j=0;j<T_size;j++)
      splitter[i] += q[j][i]*V[T_size + j];


  printf("splitter:\n");
  for(i=0;i<dim;i++)
    printf("%lf\n", splitter[i]);

  // print eigenvalues
  printf("eigenvalues: \n\t");
  for(i=0;i<T_size;i++)
    printf("%lf\n\t", a[i]);
  printf("\n");

  free(A);
  free(diagonal);
  free(scanned);
  free(V);
  free(z);

  #pragma omp parallel for
  for(i=0;i<ITER;i++)
    free(q[i]);

  MPI_Finalize();
}
