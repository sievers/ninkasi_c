#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  double i,j;
  double t1,t2;
  double sum=0;
  double n1=900*400;
  double n2=800;
  double params[4]={0.3, 0.5 ,0.7, -0.4};
  int k;
  for (i=0;i<n2;i++)
    for (j=0;j<n1;j++) {
      //sum+=cos(i+j);
      double x=i+j;
      double cur=params[3];
      for (k=2;k>=0;k--)
	cur=cur*x+params[k];
      if (i+j>1000)
	sum += cur;
    }
  printf("sum is %14.4e\n",sum);
}
