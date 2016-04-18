#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "quicksort.h"

int coordCompare (const void* a, const void* b);

int main(int argc,char* argv[]) {

  int i, j, row, col, N, status, num_row, num_col;
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
  fscanf(matrix,"%d %d %d\n",&num_row,&num_col,&N);
  if (num_row != num_col){
    fprintf(stderr,"error: this is not a sqare matrix\n");
    exit(EXIT_FAILURE);
  }

  typedef struct _coord{
       int row;
       int col;
  }coord;

  coord* A = (coord*)malloc(sizeof(coord)*2*N); //ajacency matrix coordinates
  char str[100];

  strcpy(str,"p_");

  FILE* processed = fopen(strcat(str,mstr),"w");

  //copy lower and upper triangular coordinates to memory
  j=0;
  status = fscanf(matrix,"%d %d\n",&row,&col);
  while(status == 2){
    if(row != col){
      A[j].row = row;
      A[j].col = col;

      A[j+1].row = col;
      A[j+1].col = row;

      j+=2;
      status = fscanf(matrix,"%d %d\n",&row,&col);
    }
  }

  fprintf(processed, "%d %d %d\n", num_row, num_col, j);

  // sort coordinates by row
  quicksort(A,j,sizeof(coord),coordCompare);

  // write sorted coordinates to a new file
  for(i=0;i<j;i++)
    fprintf(processed, "%d %d\n", A[i].row, A[i].col);

  return 0;
}

int coordCompare (const void* a, const void* b) {
  struct point {
    int x;
    int y;
  };

  if ((*(struct point*)a).x > (*(struct point*)b).x)
    return 1;
  if ((*(struct point*)a).x == (*(struct point*)b).x) {
    if ((*(struct point*)a).y > (*(struct point*)b).y)
      return 1;
    else if ((*(struct point*)a).y == (*(struct point*)b).y)
      return 0;
    else
      return -1;
  }
  else /*(a < b)*/
    return -1;
}
