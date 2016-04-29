#include "lanczos.h"
#include <mpi.h>
#define ITER 120

void partition(coord* A, int* myNodes, int* nodeIndex,size_t dim, size_t N, MPI_Comm comm);

int main(int argc,char* argv[]) {

  //initialize MPI environment
  MPI_Init(&argc,&argv);
  MPI_Comm world_comm = MPI_COMM_WORLD;
  int world_rank;
  MPI_Comm_rank(world_comm, &world_rank);

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

  //ajacency matrix coordinates
  coord* A = (coord*)malloc(sizeof(coord)*N);

  //scan matrix coordinates and store in memory
  for(i=0;i<N;i++)
    fscanf(matrix,"%d %d\n",&A[i].row,&A[i].col);

  fclose(matrix);


  //arrays needed to reference nodes
  int* myNodes = (int*)malloc(sizeof(int)*dim);
  int* nodeIndex = (int*)malloc(sizeof(int)*(dim+1));

  #pragma omp parallel for
  for(i=0;i<dim;i++) {
    myNodes[i] = i+1;
    nodeIndex[myNodes[i]] = i;
  }

  //Initialize file
  FILE *dotFile, *discarded;

  //Opens new file to read
  dotFile = fopen("dotFile.gc","r");
  discarded = fopen("Matrices/discarded.dat","r");

  //if the files exist delete them
  if(dotFile)
    remove("dotFile.gc");
  if(discarded)
    remove("Matrices/discarded.dat");

  //creates files to write
  dotFile = fopen("dotFile.gc","a");
  discarded = fopen("Matrices/discarded.dat","a");


  if(world_rank == 0) {
    //Writes the first two lines
    fprintf(dotFile,"graph {\n");
    fprintf(dotFile,"node [shape = point]\n");
  }

  //Closes the file
  fclose(dotFile);

  // partition the graph that we just read in
  partition(A,myNodes,nodeIndex,dim,N,world_comm);

  // barrier so that nobody tries to finish of the dot file while other threads are still working
  MPI_Barrier(MPI_COMM_WORLD);


  // open and read list of discarded edges
  discarded = fopen("Matrices/discarded.dat","r");
  coord* disc = (coord*)malloc(sizeof(coord)*N);
  int disclen=0;
  while(fscanf(discarded,"%d %d\n",&disc[disclen].row,&disc[disclen].col) != EOF){
    disclen++;
  }

  // write those edges to the dotFile
  if (world_rank == 0)
    coord2Dot(disc,disclen,0);

  // finish off dotFile
  dotFile = fopen("dotFile.gc","a");
  if(world_rank == 0)
    fprintf(dotFile,"}");


  fclose(dotFile);

  //free memory
  free(myNodes);
  free(nodeIndex);
  free(A);
  MPI_Finalize();
}

void partition(coord* A,int* myNodes,int* nodeIndex,size_t dim, size_t N, MPI_Comm comm ) {

  //figure out comm info
  int rank, color, size, new_size;
  MPI_Comm new_comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  color = rank % 2;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

  //declare memory locations that we'll need
  unsigned int i, j, l, T_size=ITER;
  int* diagonal = (int*)malloc(sizeof(int)*dim); // #n onzeros specified per row
  int* scanned = (int*)malloc(sizeof(int)*dim);  // prefix sum of the diagonal
  double* q[ITER]; //Krylov space basis matrix
  double* z = (double*)malloc(sizeof(double)*dim); //new vector in each lanczos iteration
  double a[ITER]; // diagonal of lanczos matrix
  double b[ITER]; // subdiagonal of lanczos matrix

  //initialize diagonal to zero
  #pragma omp parallel for
  for(i=0;i<dim;i++)
    diagonal[i] = 0;

  //count the number of nonzeros in each row and assign to diagonal
  //not independent... will take some work to parallelize
  for(i=0;i<N;i++)
    diagonal[nodeIndex[A[i].row]]+=1;

  // scan the diagonal so we can use it to parallelize matvec
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

  //broadcast initial guess to all other nodes
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
    matvec(A,diagonal,scanned,nodeIndex,q[j],z,N,dim);
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
      // printf("stopped after %d iterations\n", j+1);
      T_size = j+1;
      break;
    }
    if(j < ITER-1)
      #pragma omp parallel for
      for(i=0;i<dim;i++)
        q[j+1][i] = z[i]/b[j];
  }



  double end = omp_get_wtime();
  printf("time: %lf\n", end-begin);

  free(diagonal);
  free(scanned);
  free(z);

  double* V = (double*)malloc(sizeof(double)*ITER*ITER);// blank workspace for eig function
  double* splitter = (double*)malloc(sizeof(double)*dim);//second largest eigenvector

  // compute eigenvalues and eigenvectors of lanczos matrix
  eig(a,b,V,T_size);

  //compute second largest eigenvector of the laplacian matrix
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


  //use splitter to split the graph among nodes
  unsigned int AnewCount=0, newDim=0;

  //count the numer of edges in the new matrix
  #pragma omp parallel for reduction(+:AnewCount)
  for(i=0;i<N;i++)
    if (color == 0 && splitter[nodeIndex[A[i].row]] >= 0 && splitter[nodeIndex[A[i].col]] >= 0)
      AnewCount++;
    else if (color == 1 && splitter[nodeIndex[A[i].row]] < 0 && splitter[nodeIndex[A[i].col]] < 0)
      AnewCount++;

  //allocate space for the new matrix
  coord* Anew = (coord*)malloc(sizeof(coord)*AnewCount);

  // open up file to store discarded edges
  FILE* discarded;
  discarded = fopen("Matrices/discarded.dat","a");
  AnewCount=0;

  // sort edges between odd and even nodes and discarded edges
  // don't parallelize because it causes more work in the end cause you have to sort
  for(j=0;j<N;j++)
    if(color == 0 && splitter[nodeIndex[A[j].row]] >= 0 && splitter[nodeIndex[A[j].col]] >= 0){
      Anew[AnewCount].row = A[j].row;
      Anew[AnewCount].col = A[j].col;
      AnewCount++;
    }
    else if(color == 1 && splitter[nodeIndex[A[j].row]] < 0 && splitter[nodeIndex[A[j].col]] < 0){
      Anew[AnewCount].row = A[j].row;
      Anew[AnewCount].col = A[j].col;
      AnewCount++;
    }
    else if((splitter[nodeIndex[A[j].row]] < 0 && splitter[nodeIndex[A[j].col]] >= 0) ||
                            (splitter[nodeIndex[A[j].row]] >= 0 && splitter[nodeIndex[A[j].col]] < 0) )
    {
        if (rank == 0)
          fprintf(discarded,"%d %d\n", A[j].row,A[j].col);
    }

    fclose(discarded);

    // count the new dimension of the resulting matrix
    for(i=0;i<dim;i++)
      if (color == 0 && splitter[i] >= 0) {
        myNodes[newDim] = myNodes[i];
        newDim++;
      }
      else if (color == 1 && splitter[i] < 0) {
        myNodes[newDim] = myNodes[i];
        newDim++;
      }

    // assign new node index based on the nodes you get
    #pragma omp parallel for
    for(i=0;i<newDim;i++)
      nodeIndex[myNodes[i]] = i;

    //split the comm
    MPI_Comm_split(comm, color, rank, &new_comm);
    MPI_Comm_size(new_comm,&new_size);

    //if the comm size is bigger than 1 recurse
    if(new_size > 1)
      partition(Anew,myNodes,nodeIndex,newDim,AnewCount,new_comm);
    else{ //otherwise print your output to the dotFile
      int k;
      for(k=0;k<size;k++){
        if (rank == k){
          coord2Dot(Anew,AnewCount,world_rank+1);
          node2Dot(myNodes,newDim,world_rank+1);
        }
        MPI_Barrier(comm);
      }
    }
  free(Anew);
  free(splitter);
}
