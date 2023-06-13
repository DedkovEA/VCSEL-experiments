/* 
   Demonstrates creating a DLL with an exported function in a flexible and
   elegant way.
*/

#include "testclib.h"
#include <stdlib.h>
#include <stdio.h>

double ADDCALL Add(double a, double b)
{
  return (a + b);
};

double * ADDCALL AddArrays(double* a, double* b, int n)
{
  double* ans = (double*) malloc(n * sizeof(double));
  for(int i = 0; i < n; i++) {
    printf("%f %f\n", a[i], b[i]);
    ans[i] = a[i] + b[i];
  };
  return ans;
};

void ADDCALL freeArray(double* a)
{
  printf("Freing array...\n");
  free(a);
};