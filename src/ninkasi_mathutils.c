#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "ninkasi.h"
#include "ninkasi_mathutils.h"
#include "nk_clapack.h"

//#define MPI_DEBUG
#ifdef MPI_DEBUG
#include <mpi.h>
#endif

void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA, double *B, int *LDB, double  *beta, double *C, int *LDC, int transalen, int transblen );
void sgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA, float *A, int *LDA, float *B, int *LDB, float *beta, float *C, int *LDC, int transalen, int transblen );

/*--------------------------------------------------------------------------------*/

#define SSWAP(a,b) temp=(a);(a)=(b);(b)=temp;

/*!
 * Returns the kth smallest value in the input array arr[1...n].
 * Notice the twisted deal where the first valid array value is arr[1], not arr[0].  The caller has
 * to make sure not to screw that up.
 * The input array will be rearranged such that the return value lives in arr[k], and all values at
 * lower index are <=arr[k] while all values at higher index are >=arr[k].
 * If this routine looks kind of like it came from the Second Edition of a certain book of
 * programming recipes, section 8.5, then you are imagining things.  That would be unethical.
 * \param k   The index of the value that we want (counting into a sorted array).
 * \param n   Length of the data array.
 * \param arr The data array to select from.  WILL BE MODIFIED (partially sorted) by this operation.
 */

actData  sselect(unsigned long k, unsigned long n, actData *arr)
{
  unsigned long i,ir,j,l,mid;
  actData  a,temp;

  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SSWAP(arr[l],arr[ir]);
      }
      return arr[k];
    } else {
      mid=(l+ir) >> 1;
      SSWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
        SSWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SSWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SSWAP(arr[l],arr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SSWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}
#undef SSWAP


/*---------------------------------------------------------------------------------------------------------*/

/// Return the median value of an array.
/// Wraps the Num Rec sselect() function to find the Nth highest value in a data set of length M
/// (where M>=N).  Note that it accounts for the stupid 1-offset nature of arrays in Num Rec code.
/// \param ndata  Length of the data vector.
/// \param data   The data vector (0-offset, as is natural in C).
/// \return The median.

//inline static float compute_median(int ndata, float *data) 
actData compute_median(int ndata, actData *data) 
{
  return sselect(ndata/2, ndata, data-1);
}

/*--------------------------------------------------------------------------------*/
actData compute_median_inplace(int ndata, actData *data)
{
  actData *vec=vector(ndata);
  memcpy(vec,data,sizeof(actData)*ndata);
  actData val=compute_median(ndata,vec);
  free(vec);
  return val;
}

/*--------------------------------------------------------------------------------*/
actData compute_median_scat(int ndata, actData *data)
{
  actData median_val=compute_median_inplace(ndata,data);
  actData *vec=vector(ndata);
  for (int i=0;i<ndata;i++) 
    vec[i]=fabs(vec[i]-median_val);
  actData median_scat=compute_median_inplace(ndata,vec);
  free(vec);
  return median_scat;
}

/*--------------------------------------------------------------------------------*/
void add_outer_product(actData *col, int ncol, actData *row, int nrow, actData **mat)
{
#pragma omp parallel for shared(col,ncol,row,nrow,mat)
  for (int i=0;i<ncol;i++)
    for (int j=0;j<nrow;j++)
      mat[i][j]+=col[i]*row[j];
}

/*--------------------------------------------------------------------------------*/
actData *linfit(actData *d, double **vecs, actData *errs,  int n, int np)
//return best-fit parameters of d from vecs, with errors errs.  length of d is n, 
//and # of parameters is np.
{
  static int nfit=0;
  double *dcopy=dvector(n);
  for (int i=0;i<n;i++)
    dcopy[i]=d[i];
  double *atd=dvector(np);  
  memset(atd,0,np*sizeof(double));
  double **at_ni_a=dmatrix(np,np);
  memset(at_ni_a[0],0,sizeof(double)*np*np);
  if (errs)
    for (int i=0;i<n;i++)
      dcopy[i]/=(errs[i]*errs[i]);
  
  
  for (int i=0;i<n;i++)
    for (int j=0;j<np;j++)
      atd[j]+=vecs[j][i]*dcopy[i];
  
  if (errs) {
    for (int i=0;i<np;i++)
      for (int j=0;j<np;j++) 
	for (int k=0;k<n;k++) 
	  at_ni_a[i][j]+=vecs[i][k]*vecs[j][k]/(errs[k]*errs[k]);
  }
  else
    for (int i=0;i<np;i++)
      for (int j=0;j<np;j++) 
	for (int k=0;k<n;k++) 
	  at_ni_a[i][j]+=vecs[i][k]*vecs[j][k];
  
  if (invert_posdef_mat_double(at_ni_a,np))  {
    free(at_ni_a[0]);
    free(at_ni_a);
    free(dcopy);
    free(atd);
    return NULL;
  }
  
  double *pp=dvector(np);
  memset(pp,0,sizeof(double)*np);
  
  for (int i=0;i<np;i++)
    for (int j=0;j<np;j++)
      pp[i]+=at_ni_a[i][j]*atd[j];
  
  actData *params=vector(np);  
  for (int i=0;i<np;i++)
    params[i]=pp[i];
  free(pp);
  free(at_ni_a[0]);
  free(at_ni_a);  
  free(atd);
  free(dcopy);
  
  return params;
}

/*--------------------------------------------------------------------------------*/
void eval_2d_poly_inplace(actData *x, actData *y, int n, PolyParams2d *fit, actData *d)
{

  actData *xvals=vector(fit->nx);
  actData *yvals=vector(fit->ny);
  xvals[0]=1;
  yvals[0]=1;
  
  for (int i=0;i<n;i++) {
    actData xx=(x[i]-fit->xcent)/fit->xwidth;
    actData yy=(y[i]-fit->ycent)/fit->ywidth;
    d[i]=0;
    
    for (int j=1;j<fit->nx;j++)
      xvals[j]=xvals[j-1]*xx;       
    for (int j=1;j<fit->ny;j++)
      yvals[j]=yvals[j-1]*yy;
    
    for (int ix=0;ix<fit->nx;ix++)
      for (int iy=0;iy<fit->ny;iy++)
	d[i]+=xvals[ix]*yvals[iy]*fit->params[ix][iy];
  }
    
  free(xvals);
  free(yvals);
}

/*--------------------------------------------------------------------------------*/
void eval_2d_poly_pair_inplace(actData *x, actData *y, int n,const  PolyParams2d *fit, actData *d,  const PolyParams2d *fit2, actData *d2)
//do a pair of 2d poly evaluations.  Hopefully save time when evaluating pointing fits.
{
  
  actData *xvals=vector(fit->nx);
  actData *yvals=vector(fit->ny);
  xvals[0]=1;
  yvals[0]=1;
  
  for (int i=0;i<n;i++) {
    actData xx=(x[i]-fit->xcent)/fit->xwidth;
    actData yy=(y[i]-fit->ycent)/fit->ywidth;
    d[i]=0;
    d2[i]=0;
    
    for (int j=1;j<fit->nx;j++)
      xvals[j]=xvals[j-1]*xx;       
    for (int j=1;j<fit->ny;j++)
      yvals[j]=yvals[j-1]*yy;
    
    for (int ix=0;ix<fit->nx;ix++)
      for (int iy=0;iy<fit->ny;iy++) {
	d[i]+=xvals[ix]*yvals[iy]*fit->params[ix][iy];
	d2[i]+=xvals[ix]*yvals[iy]*fit2->params[ix][iy];
      }
  }
  free(xvals);
  free(yvals);
  
}

/*--------------------------------------------------------------------------------*/
actData *eval_2d_poly(actData *x, actData *y, int n, PolyParams2d *fit)
{
  actData *d=vector(n);
  eval_2d_poly_inplace(x,y,n,fit,d);
  return d;
}
/*--------------------------------------------------------------------------------*/
void destroy_2d_poly(PolyParams2d *fit)
{
  free_matrix(fit->params);
  free(fit);
}
/*--------------------------------------------------------------------------------*/
PolyParams2d *copy_2d_polyfit(PolyParams2d *fit)
{
  PolyParams2d *fit2=(PolyParams2d *)malloc(sizeof(PolyParams2d));
  memcpy(fit2,fit,sizeof(PolyParams2d));
  fit2->params=matrix(fit->nx,fit->ny);
  memcpy(fit2->params[0],fit->params[0],sizeof(actData)*fit->nx*fit->ny);
  return fit2;
}

/*--------------------------------------------------------------------------------*/
actData vec_mean(actData *vec, int n)
{
  actData mv=0;
  for (int i=0;i<n;i++)
    mv+=vec[i];
  mv /=(actData)n;
  return mv;
}
/*--------------------------------------------------------------------------------*/
actData vec_abs_err(actData *vec, int n)
{
  actData mv=vec_mean(vec,n);
  actData scat=0;
  for (int i=0;i<n;i++)
    scat +=fabs(vec[i]-mv);

  scat /=(actData)n;
  return scat;
}
/*--------------------------------------------------------------------------------*/
PolyParams2d *fit_2d_poly(actData *x, actData *y, actData *d, int n, actData *errs, int nx_param, int ny_param)
//fit a 2d polynomial to data.  x and y are the 2-d positions, d is the data, n is the length
// of x,y, and d.  nx_param and ny_param are the polynomial orders in x and y, respectively
//errs is the errors.  If they come in null, they'll get ignored.
{

  PolyParams2d *fit=(PolyParams2d *)malloc(sizeof(PolyParams2d));assert(fit);
  int np=nx_param*ny_param;
#if 1  //new way, below seems to be causing problems.
  fit->xcent=vec_mean(x,n);
  fit->ycent=vec_mean(y,n);
  fit->xwidth=vec_abs_err(x,n);
  fit->ywidth=vec_abs_err(y,n);
#else
  fit->xcent=compute_median_inplace(n,x);
  fit->ycent=compute_median_inplace(n,y);
  fit->xwidth=compute_median_scat(n,x);
  fit->ywidth=compute_median_scat(n,y);
#endif


#ifdef MPI_DEBUG
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  if (myid==0) {
    printf("cents and scats are %14.8e %14.8g %14.8g %14.8g\n",fit->xcent,fit->xwidth,fit->ycent,fit->ywidth);
    FILE *outfile=fopen("inside_stuff.txt","w");
    for (int i=0;i<n;i++) 
      fprintf(outfile,"%16.8e %16.8e %16.8e\n",x[i],y[i],d[i]);
    fclose(outfile);
  }
#endif

  fit->nx=nx_param;
  fit->ny=ny_param;
  double **vecs=dmatrix(np,n);
  for (int i=0;i<n;i++) {
    double xx= (x[i]-fit->xcent)/fit->xwidth;
    double yy= (y[i]-fit->ycent)/fit->ywidth;
    vecs[0][i]=1;
    for (int ix=0;ix<fit->nx;ix++) {
      if (ix>0)
	vecs[ix][i]=xx*vecs[ix-1][i];
      for (int iy=1;iy<fit->ny;iy++) {
	vecs[ix+iy*fit->nx][i]= vecs[ix+(iy-1)*fit->nx][i]*yy;
      }
    }    
  }
#if 1
  {
    bool am_i_broken=false;
    for (int i=0;i<n;i++) {
      if (!isfinite(x[i]))
	am_i_broken=true;
    }
    if (am_i_broken)
    printf("broken in fit on x from not finite.\n");
    
    
    am_i_broken=false;
    for (int i=0;i<n;i++) {
      if (!isfinite(y[i]))
	am_i_broken=true;
    }
    if (am_i_broken)
      printf("broken in fit on y from not finite.\n");
  }
#endif

  actData *params_1d=linfit(d,vecs,errs,n,np);


#if 1
  {
    bool am_i_broken=false;
    for (int i=0;i<n;i++) {
      if (!isfinite(x[i]))
	am_i_broken=true;
    }
    if (am_i_broken)
    printf("broken in fit on x from not finite after.\n");
    
    
    am_i_broken=false;
    for (int i=0;i<n;i++) {
      if (!isfinite(y[i]))
	am_i_broken=true;
    }
    if (am_i_broken)
      printf("broken in fit on y from not finite after.\n");
  }
#endif

  if (!params_1d) {
    fprintf(stderr,"Failure in fit_2d_poly.\n");    
    free(fit);
    fit=NULL;
  }
  else {
    fit->params=matrix(fit->nx,fit->ny);
    for (int i=0;i<fit->nx;i++)
      for (int j=0;j<fit->ny;j++)
	fit->params[i][j]=params_1d[i+fit->nx*j];
    free(params_1d);
  }

  free(vecs[0]);
  free(vecs);

  return fit;
  
}

/*--------------------------------------------------------------------------------*/
void act_gemm(char transa, char transb, int m, int n, int k, actData alpha, actData *a, int lda, actData *b, int ldb, actData beta, actData *c, int ldc)

{
#ifdef ACTDATA_DOUBLE
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
#else
  sgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
#endif

}

/*--------------------------------------------------------------------------------*/


int  invert_posdef_mat(actData **mat, int n)
{
  
  int info;
#ifdef ACTDATA_DOUBLE
  clapack_dpotrf('u', n,mat[0], n, &info);
#else
  clapack_spotrf('u', n,mat[0], n, &info);
#endif
  if (info)
    return info;
#ifdef ACTDATA_DOUBLE
  clapack_dpotri('u', n, mat[0], n, &info);
#else
  clapack_spotri('u', n, mat[0], n, &info);
#endif
  if (info)
    return info;
  for (int i = 0; i < n; i++)
    for (int j = i+1; j < n; j++)
      mat[i][j] = mat[j][i];
  return info;
}

/*--------------------------------------------------------------------------------*/


int  invert_posdef_mat_double(double **mat, int n)
{
  
  int info;
  clapack_dpotrf('u', n, mat[0], n, &info);
  if (info)
    return info;
  clapack_dpotri('u', n,mat[0], n, &info);
  if (info)
    return info;
  for (int i = 0; i < n; i++)
    for (int j = i+1; j < n; j++)
      mat[i][j] = mat[j][i];
  return info;
}


/*--------------------------------------------------------------------------------*/
actData vecmin(actData *vec, int n)
{
  actData minval;
  minval=vec[0];
  for (int i=1;i<n;i++)
    if (vec[i]<minval)
      minval=vec[i];
  return minval;
}
/*--------------------------------------------------------------------------------*/
float svecmin(float *vec, int n)
{
  float minval;
  minval=vec[0];
  for (int i=1;i<n;i++)
    if (vec[i]<minval)
      minval=vec[i];
  return minval;
}
/*--------------------------------------------------------------------------------*/
int ivecmin(int *vec, int n)
{
  int minval;
  minval=vec[0];
  for (int i=1;i<n;i++)
    if (vec[i]<minval)
      minval=vec[i];
  return minval;
}
/*--------------------------------------------------------------------------------*/
int ivecmax(int *vec, int n)
{
  int maxval;
  maxval=vec[0];
  for (int i=1;i<n;i++)
    if (vec[i]>maxval)
      maxval=vec[i];
  return maxval;
}
/*--------------------------------------------------------------------------------*/
actData vecmax(actData *vec, int n)
{
  actData maxval;
  maxval=vec[0];
  for (int i=1;i<n;i++)
    if (vec[i]>maxval)
      maxval=vec[i];
  return maxval;
}

/*--------------------------------------------------------------------------------*/
float svecmax(float *vec, int n)
{
  float maxval;
  maxval=vec[0];
  for (int i=1;i<n;i++)
    if (vec[i]>maxval)
      maxval=vec[i];
  return maxval;
}




/*--------------------------------------------------------------------------------*/

actData myrand(unsigned *seed)
{
  int ival=rand_r(seed);
  actData val=(actData)ival/(RAND_MAX);
  return val;
}
/*--------------------------------------------------------------------------------*/
actData mygasdev(unsigned *seed)
{
  actData r1=myrand(seed);
  actData r2=myrand(seed);
  actData r=r1*r1+r2*r2;
  while ((r>=1)||(r<=0)) {
    r1=myrand(seed);
    r2=myrand(seed);
    r=r1*r1+r2*r2;
  }
  actData r3=myrand(seed);
  if (r3>0.5)
    return (actData)(r1*sqrt(-2*log(r)/r));
  else
    return (actData)(-r1*sqrt(-2*log(r)/r));
  
}



