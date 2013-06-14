#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <nk_clapack.h>

void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc, int uplolen, int translen);
void ssyrk_(char *uplo, char *trans, int *n, int *k, float *alpha, float *a, int *lda, float *beta, float *c, int *ldc, int uplolen, int translen);



/*--------------------------------------------------------------------------------*/

void clapack_dsyrk(char uplo, char trans, int n, int k, double alpha, double *a, int lda, double beta, double *c, int ldc)
{
  dsyrk_(&uplo,&trans,&n,&k,&alpha,a,&lda,&beta,c,&ldc,1,1);
}


/*--------------------------------------------------------------------------------*/

void clapack_ssyrk(char uplo, char trans, int n, int k, float alpha, float *a, int lda, float beta, float *c, int ldc)
{
  ssyrk_(&uplo,&trans,&n,&k,&alpha,a,&lda,&beta,c,&ldc,1,1);
}




/*--------------------------------------------------------------------------------*/

void clapack_sgesv(int n, int nrhs, float *a, int lda, int *ipiv, float *b, int ldb, int *info)
{
  sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info);
}

/*--------------------------------------------------------------------------------*/

void clapack_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb, int *info)
{
  dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info);
}

/*--------------------------------------------------------------------------------*/

void clapack_spotri(char uplo, int n, float *a, int lda, int *info)
{
  spotri_(&uplo,&n,a,&lda,info);
}
/*--------------------------------------------------------------------------------*/

void clapack_dpotri(char uplo, int n, double *a, int lda, int *info)
{
  dpotri_(&uplo,&n,a,&lda,info);
}

/*--------------------------------------------------------------------------------*/

void clapack_spotrf(char uplo, int n, float *a, int lda, int *info)
{
  spotrf_(&uplo,&n,a,&lda,info);
}

/*--------------------------------------------------------------------------------*/

void clapack_dpotrf(char uplo, int n, double *a, int lda, int *info)
{
  dpotrf_(&uplo,&n,a,&lda,info);
}


/*--------------------------------------------------------------------------------*/

void clapack_ssyev(char jobz, char uplo, int n, float *a, int lda, float *w, float *work, int lwork, int *info)
{
  ssyev_(&jobz,&uplo,&n,a,&lda,w,work,&lwork,info,1,1);
}

/*--------------------------------------------------------------------------------*/

void clapack_ssyevd(char jobz, char uplo, int n, float *a, int lda, float *w, float *work, int lwork,int *iwork, int liwork, int *info)
{
  ssyevd_(&jobz,&uplo,&n, a,&lda, w,work,&lwork,iwork,&liwork,info,1,1);
}


/*--------------------------------------------------------------------------------*/

void clapack_dsyevd(char jobz, char uplo, int n, double *a, int lda, double *w, double *work, int lwork,int *iwork, int liwork, int *info)
{
  dsyevd_(&jobz,&uplo,&n, a,&lda, w,work,&lwork,iwork,&liwork,info,1,1);
}


/*--------------------------------------------------------------------------------*/

void clapack_ssyevd_simple(char jobz, char uplo, int n, float *a, int lda, float *w, int *info)
/*allocate scratch space for the user & make all the calls*/
{
  float fwork;
  int liwork,lwork;

  clapack_ssyevd(jobz,uplo,n,a,lda,w,&fwork,-1,&liwork,-1,info);

  lwork=fwork+1;  //extra number just in case..
  float *work=(float *)malloc(sizeof(float)*lwork);
  int *iwork=(int *)malloc(sizeof(int)*liwork);

  clapack_ssyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info);

  free(work);
  free(iwork);
  
}

/*--------------------------------------------------------------------------------*/

void clapack_dsyevd_simple(char jobz, char uplo, int n, double *a, int lda, double *w, int *info)
/*allocate scratch space for the user & make all the calls*/
{
  double fwork;
  int liwork,lwork;

  clapack_dsyevd(jobz,uplo,n,a,lda,w,&fwork,-1,&liwork,-1,info);

  lwork=fwork+1;  //extra number just in case..
  double *work=(double *)malloc(sizeof(double)*lwork);
  int *iwork=(int *)malloc(sizeof(int)*liwork);

  clapack_dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info);

  free(work);
  free(iwork);
  
}
/*--------------------------------------------------------------------------------*/
int dinvsafe(double *a, int n, double thresh)
//invert a potentially singular positive definite matrix. Any eigenvalues below
//thresh*max(eigenvalues) will be replaced by zero in the inverse.
{
  double *myeigs=(double *)malloc(sizeof(double)*n);
  double *mytmp=(double *)malloc(sizeof(double)*n*n);
  //memcpy(mytmp,a,n*n*sizeof(double));
  for (int i=0;i<n*n;i++)
    mytmp[i]=a[i];
  int info;
  clapack_dsyevd_simple('V','U',n,mytmp,n,myeigs,&info);
  double max_eig=myeigs[0];
  for (int i=1;i<n;i++)
    if (myeigs[i]>max_eig)
      max_eig=myeigs[i];

  int i_am_sketchy=0;
  for (int i=0;i<n;i++) {
    if (myeigs[i]<thresh*max_eig)  {
      i_am_sketchy=1;
      myeigs[i]=0;      
  }
    else
      myeigs[i]=1.0/sqrt(myeigs[i]);
  }
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++)
      mytmp[i*n+j]*=myeigs[i];
  }
  clapack_dsyrk('U','N',n,n,1.0,mytmp,n,0.0,a,n);
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++) {
      a[i*n+j]=a[j*n+i];
    }
  free(myeigs);
  free(mytmp);
  return i_am_sketchy;
}
/*--------------------------------------------------------------------------------*/

void clapack_ssyev_simple(char jobz, char uplo, int n, float *a, int lda, float *w, int *info)
/*allocate scratch space for the user & make all the calls*/
{
  float fwork;

  clapack_ssyev(jobz,uplo,n,a,lda,w,&fwork,-1,info);
  int lwork=fwork+1;  //extra number just in case..
  float *work=(float *)malloc(sizeof(float)*lwork);
  clapack_ssyev(jobz,uplo,n,a,lda,w,work,lwork,info);
  free(work); 
}
/*--------------------------------------------------------------------------------*/
#if 0
//shouldn't need these, but they appear to be missing
void cblas_sgemm();

void cblas_sdot(int n, float *x, int incx, float *y, int incy)
{
  sdot_(&n,x,&incx,y,&incy);
}
#endif





