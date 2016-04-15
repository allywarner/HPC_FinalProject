#include <stdio.h>
#include <stdlib.h>

long dsteqr(char COMPZ, long N, double *D, double *E, double *Z, long LDZ,
      double *WORK);

int eig(double* d, double* e, double* z,size_t n)
{
  unsigned int i, j, info;
  double* work = (double*)malloc(sizeof(double)*n*n);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if (i == j)
        z[n*i+j] = 1;
      else
        z[n*i+j] = 0;

  info = dsteqr('V', n, d, e, z, n, work);
  if (info != 0)
    fprintf(stderr, "failure with error %d\n", info);

  free(work);

  return info;
}

long dsteqr(char COMPZ, long N, double *D, double *E, double *Z, long LDZ,
      double *WORK)
{
  extern void dsteqr_(const char *COMPZp, const long *Np, double *D,
		     double *E, double *Z, const long *LDZp, double *WORK,
		     long *INFOp);
  long info;
  dsteqr_(&COMPZ, &N, D, E, Z, &LDZ, WORK, &info);
  return info;
}
