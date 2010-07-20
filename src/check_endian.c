#include <stdio.h>

void assign_byte(void *ptr)
{
  char *c=ptr;
  c[0]=1;
}


/*--------------------------------------------------------------------------------*/
void swap_4bytes(void *ptr)
{
  char c,*cc;
  cc=ptr;
  c=cc[0];
  cc[0]=cc[3];
  cc[3]=c;
  c=cc[1];
  cc[1]=cc[2];
  cc[2]=c;

}

/*--------------------------------------------------------------------------------*/
void swap_8bytes(void *ptr)
{
  char c,*cc;
  cc=ptr;
  c=cc[0];
  cc[0]=cc[7];
  cc[7]=c;

  c=cc[1];
  cc[1]=cc[6];
  cc[6]=c;

  c=cc[2];
  cc[2]=cc[5];
  cc[5]=c;


  c=cc[3];
  cc[3]=cc[4];
  cc[4]=c;

}


/*--------------------------------------------------------------------------------*/
int main(int argc, char *argv[])

{
  int i=0;
  float x=0;
  double dx=0;

  assign_byte(&i);
  assign_byte(&x);
  assign_byte(&dx);
  printf("values are %d %12.4g %12.4lg\n",i,x,dx);

  swap_4bytes(&i);
  swap_4bytes(&x);
  //float *dd=(float *)(&dx);
  //swap_4bytes((void *)dd);
  //swap_4bytes(&dd[1]);
  swap_8bytes(&dx);
  printf("post-swap, i and x are %d %12.4g %12.4lg\n",i,x,dx);
}
