#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>


int main(int argc, char *argv[])
{
  int n=10;
  if (argc>1)
    n=atoi(argv[1]);
  float complex *vec=(float complex *)malloc(sizeof(float complex)*n);
  for (int i=0;i<n;i++)
    vec[i]=i;
  float complex fac=1-2*I;
  for (int i=0;i<n;i++)
    vec[i]*=fac;
  for (int i=0;i<n;i++)
    printf("%3d %4.0f %4.0f\n",n,creal(vec[i]),cimag(vec[i]));
  
  float complex crud=cexpf(fac);
  printf("crud is %14.5f %14.5f\n",creal(crud),cimag(crud));
}
