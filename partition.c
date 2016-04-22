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
  char mstr[100];

  // read name of matrix from arg list
  if (argc > 1)
  strcpy(mstr,argv[1]);
  else if (world_rank == 0){
    fprintf(stderr, "error: no matrix path provided\n");
    exit(EXIT_FAILURE);
  }

  FILE* matrix = fopen(mstr,"r");

  // open matrix file
  if (matrix == NULL && world_rank == 0){
    fprintf(stderr, "error: matrix path invalid\n");
    exit(EXIT_FAILURE);
  }

  //scan first row of matrix for metadata
  fscanf(matrix,"%d %d %d\n",&row,&col,&N);

  if (row != col && world_rank == 0){
    fprintf(stderr,"error: this is not a sqare matrix\n");
    exit(EXIT_FAILURE);
  }
  dim = row;

  coord* A = (coord*)malloc(sizeof(coord)*N); //ajacency matrix coordinates

  //scan matrix coordinates and store in memory
  for(i=0;i<N;i++)
    fscanf(matrix,"%d %d\n",&A[i].row,&A[i].col);

  fclose(matrix);

  partition(A, dim, N, world_size, world_comm);

}

void partition(coord* A, size_t dim, size_t N, int comm_size, MPI_Comm comm ) {

  int rank, size, color;
  MPI_Comm new_comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  color = rank % 2;

  unsigned int i, j, l, T_size=ITER;
  int* diagonal = (int*)malloc(sizeof(int)*dim); // #nonzeros specified per row
  int* scanned = (int*)malloc(sizeof(int)*dim);
  double* q[ITER]; //Krylov space basis matrix
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
  if (rank == 0)
  {
  #pragma omp parallel for
  for(i=0;i<dim;i++)
    q[0][i] = (double)rand()/(double)RAND_MAX * 100;
  }

  MPI_Bcast(&q[0][0],dim, MPI_DOUBLE, 0, comm);

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

  free(diagonal);
  free(scanned);
  free(z);

  double* V = (double*)malloc(sizeof(double)*ITER*ITER);// blank workspace for eig function
  double* splitter = (double*)malloc(sizeof(double)*dim);

  // compute eigenvalues and eigenvectors of lanczos matrix
  eig(a,b,V,T_size);

  #pragma omp parallel for
  for(i=0;i<dim;i++)
    splitter[i] = 0;

  #pragma omp parallel for private(j)
  for(i=0;i<dim;i++)
    for(j=0;j<T_size;j++)
      splitter[i] += q[j][i]*V[T_size + j];


  free(V);

  #pragma omp parallel for
  for(i=0;i<ITER;i++)
    free(q[i]);


  printf("splitter:\n");
  for(i=0;i<dim;i++)
    printf("%lf\n", splitter[i]);

  // print eigenvalues
  printf("eigenvalues: \n\t");
  for(i=0;i<T_size;i++)
    printf("%lf\n\t", a[i]);
  printf("\n");

  unsigned int AnewCount=0, newDim=0;

  #pragma omp parallel for reduction(+:AnewCount)
  for(i=0;i<N;i++)
    if (color == 0 && splitter[A[i].row] >= 0 && splitter[A[i].col] >= 0)
      AnewCount++;
    else if (color == 1 && splitter[A[i].row] < 0 && splitter[A[i].col] < 0)
      AnewCount++;

  #pragma omp parallel for reduction(+:newDim)
  for(i=0;i<dim;i++)
    if (color == 0 && splitter[i] >= 0)
      newDim++;
    else if (color == 1 && splitter[i] < 0)
      newDim++;

  coord* Anew = (coord*)malloc(sizeof(coord)*AnewCount);

  omp_lock_t writelock;
  omp_init_lock(&writelock);

  AnewCount=0;




  #pragma omp parallel for
  for(j=0;j<N;j++)
    if(color == 0 && splitter[A[j].row] >= 0 && splitter[A[j].col] >= 0){
      omp_set_lock(&writelock);
      Anew[AnewCount].row = A[j].row;
      Anew[AnewCount].col = A[j].col;
      AnewCount++;
      omp_unset_lock(&writelock);
    }
    else if(color == 1 && splitter[A[j].row] < 0 && splitter[A[j].col] < 0){
      omp_set_lock(&writelock);
      Anew[AnewCount].row = A[j].row;
      Anew[AnewCount].col = A[j].col;
      AnewCount++;
      omp_unset_lock(&writelock);
    }

  omp_destroy_lock(&writelock);

  FILE* fp;

  printf("%d %d\n", color, newDim);

  if(color == 0)
    fp = fopen("Matrices/Aplus.dat","w");
  else
    fp = fopen("Matrices/Aminus.dat","w");

  for(i=0;i<AnewCount;i++)
    fprintf(fp,"%d %d\n",Anew[i].row,Anew[i].col);



  // MPI_Comm_split(comm, color, rank, &new_comm);
  // partition(Anew,newDim,AnewCount,size,new_comm);


  free(A);
  free(Anew);
  free(splitter);

  MPI_Finalize();
}
