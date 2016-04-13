#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N 167936

double* matvec(FILE* A, double* x) {

  int status, row, col;

  do {
    status = fscanf(A,"%d %d\n",&row,&col);




  } while (status != EOF);

  return x;
}

int main(int argc,char* argv[]) {

  srand(time(NULL));
  double* x =  (double*)malloc(sizeof(double)*N);
  int i;

  for(i=0;i<N;i++)
    x[i] = ((double)rand()/(double)RAND_MAX * 100);


  FILE* A = fopen("finance256.dat","r");

  matvec(A,x);

  return 0;
}
