#include "lanczos.h"

void scan(void* base, size_t n,size_t l, void(*oper)(void* x1, void* x2)) {

  /*use binary tree scan if array is small enough for available processors*/
  unsigned int num_chunks = omp_get_max_threads();
  if (n <= num_chunks)
    rec_scan(base,n,l,oper);

  /*sequentially scan the tails first (in parallel)
    and then use binary tree*/
  else {
    char* array = (char*)base;
    unsigned int i, j, chunk_size;
    chunk_size = n/num_chunks;

    /*sequentially scan tails*/
    #pragma omp parallel for private(i)
    for (j = 0; j < num_chunks; j++)
      for (i = j*chunk_size; i < (j+1)*chunk_size - 1; i++)
        oper(array + i*l,array + (i+1)*l);

     /*recursively scan the top of the tree*/
     rec_scan(array + l*(chunk_size - 1), num_chunks, l*chunk_size, oper);

     /*add the top of each tail to the subsequent tail*/
     for (j = 1; j < num_chunks; j++)
     #pragma omp parallel for
      for (i = j*chunk_size; i < (j+1)*chunk_size - 1; i++)
        oper(array + (j*chunk_size-1)*l,array + i*l);

     /*scan the remainder*/
     scan(array + (num_chunks*chunk_size - 1)*l, (n % num_chunks) + 1, l, oper);
  }

}

void rec_scan(void *base, size_t n,size_t l, void (*oper)(void* x1, void* x2)) {
  if (n > 1) {
    char* array = (char*)base;
    unsigned int i;

    #pragma omp parallel for //upsweep
    for (i = 0; i < n - 1; i+=2)
      oper(array + i*l,array + (i+1)*l);

    rec_scan(array + l, n/2, l*2, oper);

    #pragma omp parallel for //downsweep
    for (i = 1; i < n - 1; i+=2)
      oper(array + i*l, array + (i+1)*l);
  }
}

void addInt(void* x1, void* x2) {
  *(int*)x2 = *(int*)x1 + *(int*)x2;
}

void addDouble(void* x1, void* x2) {
    *(double*)x2 = *(double*)x1 + *(double*)x2;
}

void addVec3(void* x1, void* x2) {
  struct vec3 {
    double x;
    double y;
    double z;
  };

  (*(struct vec3*)x2).x = (*(struct vec3*)x1).x + (*(struct vec3*)x2).x;
  (*(struct vec3*)x2).y = (*(struct vec3*)x1).y + (*(struct vec3*)x2).y;
  (*(struct vec3*)x2).z = (*(struct vec3*)x1).z + (*(struct vec3*)x2).z;
}
