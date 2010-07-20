
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "ninkasi.h"
#include "noise.h"
#include "mbCommon.h"
#ifndef _SKIP_MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif
#include <nk_clapack.h>
#include "ninkasi_mathutils.h"

//#include <mkl.h>
#define NOISE_FIT_WIDTH 10  //yes, need to put this in a function somewhere...

/*--------------------------------------------------------------------------------*/
void SetNoiseType(NoiseParams1Pix *noise, mbNoiseType noise_type)
{
  noise->noise_type=noise_type;  
}
/*--------------------------------------------------------------------------------*/
actData get_tod_fft_delta(mbTOD *tod)
{
  return 1.0/(tod->deltat*tod->ndata);
}
/*--------------------------------------------------------------------------------*/
void SetMaxFreq(NoiseParams1Pix *noise, actData maxFreq)
{
  noise->maxfreq=maxFreq;  
}
/*--------------------------------------------------------------------------------*/
void SetMinFreq(NoiseParams1Pix *noise, actData minFreq)
{
  noise->minfreq=minFreq;  
}
/*--------------------------------------------------------------------------------*/
void SetPowlaw(NoiseParams1Pix *noise, actData powlaw)
{
  noise->powlaw=powlaw;
}
/*--------------------------------------------------------------------------------*/
void SubtractMedian(actData *vec, int n)
{
  actData val=compute_median_inplace(n,vec);
  for (int i=0;i<n;i++)
    vec[i]-=val;
}

/*--------------------------------------------------------------------------------*/
int fft_real2complex_nelem(int n)
{
  return (n/2+1);
}

/************************************************************************************************************/
/*!
 * Assign the knee to a noise structure
 */

void nkAssignNoiseKnee( mbTOD *tod, actData knee) 
{
  assert(tod->noise);
  assert(tod->cuts);
  for (int i=0;i<tod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
      NoiseParams1Pix *noise=&(tod->noise->noises[i]);
      if (noise->params[0]<=0)  //If we don't have any estimate of white, then put something non-zero in there.
	noise->params[0]=1;  
      noise->params[1]=noise->params[0]/pow(knee,noise->powlaw);
      noise->old_knee=noise->knee;
      noise->knee=knee;
  
  noise->converged=2;  //let the world know we think we have a noise model, but not one we found ourselves.

    }
  }
  //assert(noise->noise_type==MBNOISETYPE1PIX_1F);
  
}



/*--------------------------------------------------------------------------------*/

void set_tod_noise(mbTOD *tod, actData white, actData knee, actData powlaw)
//set the noise in a TOD.
{
  assert(tod);
  if (!tod->noise) {
    printf("allocating tod noise space.\n");
    tod->noise=(  mbNoiseVectorStruct *)calloc(1,sizeof(mbNoiseVectorStruct));
    tod->noise->ndet=tod->ndet;
    tod->noise->noises=(  NoiseParams1Pix *)calloc(tod->ndet,sizeof(NoiseParams1Pix));
  }
  for (int i=0;i<tod->ndet;i++) {
    NoiseParams1Pix *mynoise=&(tod->noise->noises[i]);
    SetNoiseType(mynoise,MBNOISE_LINEAR_POWLAW);
    SetPowlaw(mynoise,powlaw);
    mynoise->params[0]=white*white/(tod->deltat);  //the deltat is there 'cause white is noise per second, not per sample
    
  }
  nkAssignNoiseKnee(tod,knee);
  
}

/*--------------------------------------------------------------------------------*/


mbNoiseVectorStruct *nkFitTODNoise(mbTOD *tod, mbNoiseType noise_type, actData minFreq, actData maxFreq, actData powlaw)
{
  assert (tod != NULL);
  assert (tod->ndata > 1);
  assert(tod->have_data);
  assert(maxFreq>minFreq);
  
  createFFTWplans1TOD(tod);  //make sure the fft plans are there.
  mbNoiseVectorStruct *noises=(  mbNoiseVectorStruct *)calloc(1,sizeof(mbNoiseVectorStruct));
  noises->ndet=tod->ndet;
  noises->noises=(  NoiseParams1Pix *)calloc(noises->ndet,sizeof(NoiseParams1Pix));
  for (int i=0;i<noises->ndet;i++) {
    SetNoiseType(&(noises->noises[i]),noise_type);
    SetMaxFreq(&(noises->noises[i]),maxFreq);
    SetMinFreq(&(noises->noises[i]),minFreq);
    SetPowlaw(&(noises->noises[i]),powlaw);    
  }
#pragma omp parallel for shared(noises,tod) default(none) schedule(dynamic,1)
  for (int i=0;i<noises->ndet;i++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) 
      nkFitNoise1Det(&(noises->noises[i]),tod,i);
  }
  return noises;

}


/*--------------------------------------------------------------------------------*/

int nkFitNoise1Det(NoiseParams1Pix *noise,mbTOD *tod,int mydet) 
{
  switch(noise->noise_type) {
  case MBNOISE_LINEAR_POWLAW:
    return nkFitNoise_LinearPowlaw(noise,tod,mydet);
    break;
  }

  return -1;
}
/*--------------------------------------------------------------------------------*/

int nkFitNoise_LinearPowlaw(NoiseParams1Pix *noise,mbTOD *tod,int mydet) 
{
  assert(tod->have_data);
  assert(mydet<tod->ndet);

  int n=tod->ndata;
  //int nn=n/2+1;
  int nn=fft_real2complex_nelem(n);

  
  actData *dat=CopyData1Det(tod,mydet);
  SubtractMedian(dat,n);
  act_fftw_complex *fdata_raw=act_fftw_malloc(nn*sizeof(act_fftw_complex));
  act_fftw_execute_dft_r2c(tod->p_forward,dat,fdata_raw);
  int nn_start, nn_stop;
  actData delta;
  nkSetNoiseFreqRange(tod,noise,&nn_start,&nn_stop,&delta);
  act_fftw_complex *fdata;
  if (nn_start>0)
    fdata=fdata_raw+nn_start;
  else
    fdata=fdata_raw;
  actData nd=nn; 
  actData **vecs;
  int nparam;
  nkCalculateOneoverFVecs(tod,noise,&vecs,&nparam);
  nn=nn_stop-nn_start;
  actData  *fit_params=vector(nparam);
  if (!noise->initialized) {
    if (nkFitStartingSpecUncorrData(noise->params,(actData *)fdata,2*nn,vecs,nparam)) {  //we failed
      noise->converged=0;
      act_fftw_free(fdata_raw);
      psFree(dat);
      psFree(fit_params);
      free_matrix(vecs);
      return  PS_ERR_BAD_PARAMETER_VALUE;
    }
    for (int i=0;i<nparam;i++)
      fit_params[i]=noise->params[i];
  }

  if (nkFitSpecUncorrDataQuadratic(fit_params,(actData*)fdata,2*nn,vecs,nparam)==PS_ERR_NONE){
    //fprintf(stderr,"Quadratic fit did not puke.\n");
    for (int i=0;i<nparam;i++)
      noise->params[i]=fit_params[i];
    
  }
#if 1
  actData like;
  if ((nkFitNoiseUncorrData(fit_params,(actData *)fdata,NULL,vecs,2*nn,nparam,&like)==PS_ERR_NONE) &&(isfinite(fit_params[0]))) { //we converged
    for (int i=0;i<nparam;i++)
      noise->params[i]=fit_params[i]; 
    noise->converged=1;
    if (noise->params[0]<0)
      noise->converged=0;  //we probably didn't do the right thing in this case.
    else
      noise->knee=pow(noise->params[0]/noise->params[1],1.0/noise->powlaw);
  }
  else
    noise->converged=0;
#if 0
  if (noise->converged)
    printf("Noise fit succeeded on %d with knee %14.4e and white %14.4e\n",mydet,noise->knee,sqrt(noise->params[1]));
  else
    printf("Noise fit failed on %d\n",mydet);
#endif

#endif
  free_matrix(vecs);
  act_fftw_free(fdata_raw);
  psFree(dat);
  psFree(fit_params);
  return PS_ERR_NONE;
}



/************************************************************************************************************/
/*!
 *  Calculate the vectors for a 1/f fit.
 */

psErrorCode nkCalculateOneoverFVecs(mbTOD *tod, NoiseParams1Pix *noise_params,actData ***vecs_out,int *nparam_out)
{
  actData  **vecs;
  actData isq;
  
  assert (tod != NULL);
  assert (tod->ndata > 1);
  
  int n=tod->ndata;
  //int nn=n/2+1;  //how many element we need in the fft  
  int nn=fft_real2complex_nelem(n);
  actData nd=nn;

  actData delta;
  int nn_start,nn_max;
  nkSetNoiseFreqRange(tod,noise_params,&nn_start,&nn_max,&delta);
  //printf("nnstart and nn_max are %d %d, with delta %14.5e\n",nn_start,nn_max,delta);
  
  int nparam=2;  /*since we're fitting 1/f + constant*/
  
  vecs=matrix(nparam,2*(nn_max-nn_start));
  for (int i=0;i<nn_max-nn_start;i++) {
    vecs[0][2*i]=1.0*nd;
    vecs[0][2*i+1]=1.0*nd;
    isq=i+nn_start;
    if (isq>0)	{
      vecs[1][2*i]=pow(isq,noise_params->powlaw)*nd;
      vecs[1][2*i+1]=pow(isq,noise_params->powlaw)*nd;
    }
    else {
      vecs[1][2*i]=0;
      vecs[1][2*i+1]=pow(isq+0.5,noise_params->powlaw)*nd;
    }
    
    vecs[1][2*i]*=pow(delta,noise_params->powlaw);
    vecs[1][2*i+1]*=pow(delta,noise_params->powlaw);
  }
  
  *vecs_out=vecs;
  *nparam_out=nparam;
  return PS_ERR_NONE;
}



/************************************************************************************************************/
/*!
 *  Figure out starting/stopping pixel values for filtering across selected frequencies
 */
psErrorCode nkSetNoiseFreqRange(mbTOD *tod, NoiseParams1Pix *noise_params, int *nn_start_out, int *nn_max_out, actData *delta_out)
{
  int n,nn,nn_start,nn_max;
  

  assert(tod->deltat>0);
  //actData delta=1.0/(tod->deltat*tod->ndata);
  actData delta=get_tod_fft_delta(tod);
  
  n=tod->ndata;
  //nn=n/2+1;  //how many element we need in the fft
  nn=fft_real2complex_nelem(n);
  if (noise_params->minfreq>0)  {
      nn_start=noise_params->minfreq/delta;
      if (noise_params->minfreq-delta*(actData)nn_start>0)
	nn_start++;
    }
  else     /*no low-end to cut out here.*/      
    nn_start=0;
  
  if (noise_params->maxfreq>0) {  //Put in a maximum frequency to use, if so desired   
    nn_max=noise_params->maxfreq/delta;
    if (nn_max<nn)
      nn=nn_max;
  }
  else
    nn_max=nn;

  *nn_start_out=nn_start;
  *nn_max_out=nn_max;
  if (delta_out!=NULL)
    *delta_out=delta;

  return PS_ERR_NONE;

}


/************************************************************************************************************/
/*!
 *  Fit a quick-and-dirty set of parameters describing the variance of a set of uncorrelated data.  
 *  Returns the answer assuming that the variance on each point is the same, then does a linear fit to
 *  the square of the data.
 */

psErrorCode nkFitStartingSpecUncorrData(actData *params,actData *data,int  n,actData **vecs,int nparam)
{
  
  actData **mat=matrix(nparam,nparam);
  actData *ax=vector(nparam);
  
  for (int i=0;i<nparam;i++) {
    ax[i]=0;
    for (int j=0;j<n;j++)
      ax[i]+=vecs[i][j]*data[j]*data[j];
  }
  for (int i=0;i<nparam;i++)
    for (int j=0;j<nparam;j++)
      mat[i][j]=mbDot(n,vecs[i],1,vecs[j],1);

    
  if (mbInvertPosdefMat(mat,nparam))
    {
      fprintf(stderr,"Failure in inverting A^T A in mbFitStartingSpecUncorrData.\n");
      return PS_ERR_BAD_PARAMETER_VALUE;
    }
  
  for (int i=0;i<nparam;i++) {
    params[i]=0;
    for (int j=0;j<nparam;j++)
      params[i]+=mat[i][j]*ax[j];      
  }
  
  psFree(mat[0]);
  psFree(mat);
  psFree(ax);
  
  return PS_ERR_NONE;

}

/************************************************************************************************************/
/*!
 *  call ?dot
 */

actData mbDot(long n,actData *x,long incx,actData *y,long incy)
{
#ifndef ACTDATA_DOUBLE
  return cblas_sdot(n,x,incx,y,incy);
#else
  return cblas_ddot(n,x,incx,y,incy);
#endif
}

/************************************************************************************************************/
/*!
 *  call ?axpy
 */

void mbDaxpy(long n,actData a,actData *x,long incx,actData *y,long incy)
{
#ifndef ACTDATA_DOUBLE
  cblas_saxpy(n,a,x,incx,y,incy);
#else
  cblas_daxpy(n,a,x,incx,y,incy);
#endif
}


/************************************************************************************************************/
/*!
 *  Use a quadratic estimator to calculate the maximum-likelihood solution.
 *  Turns out that the standard linear least-squares fit to the data^2 produces the ML
 *  solution if the weights used are the noise^-4 instead of noise^-2.  Seems to work pretty well
 *  and is more robust than using Newton's method (at least in my simple tests).
 */

psErrorCode nkFitSpecUncorrDataQuadratic(actData *params,actData *data,long n,actData **vecs,int nparam)
{
  actData *cov,**mat,*ax,*last_params;
  int i,j,k,converged,niter,initialized; //cov_failed;

  
  /*first, if all parameters are equal to zero, assume we're starting off
    fresh.  Need to get an initial guess.*/
  initialized=0;
  for (i=0;i<nparam;i++)
    if (params[i]!=0)
      initialized=1;
  
  niter=0;
  if (initialized==0) {
    nkFitStartingSpecUncorrData(params,data,n,vecs,nparam);  //get a starting guess in.
    niter++;
  }
  

  cov=vector(n);
  mat=matrix(nparam,nparam);
  ax=vector(nparam);
  last_params=vector(nparam);
  converged=0;
  while (converged==0) {
    niter++;    
    for (i=0;i<nparam;i++)
      last_params[i]=params[i];
    
    /*first calculate covariance, since we use its inverse square as the weights.*/
    /*put the C^-2 in here since it gets used multiple times down the road.*/
    for (i=0;i<n;i++)  {
      cov[i]=0;
      for (j=0;j<nparam;j++)
	cov[i]+=params[j]*vecs[j][i];
      if (cov[i]<0)  {
	psTrace("moby.pcg",6,"Hit negative covariance in mbFitSpecUncorrDataQuadratic in iteration %d.\n",niter);
	free(last_params);
	free(cov);
	free(ax);
	free_matrix(mat);
      
	return PS_ERR_BAD_PARAMETER_VALUE;	      
      }
      cov[i]=1/cov[i]/cov[i];
    }
    
    /*Fill up A^T C^-2 A .*/
    for (i=0;i<nparam;i++)
      for (j=i;j<nparam;j++) {
	mat[i][j]=0;
	for (k=0;k<n;k++)
	  mat[i][j]+=vecs[i][k]*vecs[j][k]*cov[k];
	mat[j][i]=mat[i][j];
      }
    /*now do A^T C^-2 x^2*/
    for (i=0;i<nparam;i++)  {
      ax[i]=0;
      for (j=0;j<n;j++)
	ax[i]+=vecs[i][j]*cov[j]*data[j]*data[j];
    }
    if (mbInvertPosdefMat(mat,nparam)) {
      fprintf(stderr,"Failure in inverting A^T A in mbFitStartingSpecUncorrData.\n");
      free(last_params);
      free(cov);
      free(ax);
      free_matrix(mat);
      return PS_ERR_BAD_PARAMETER_VALUE;
    }
    
    for (i=0;i<nparam;i++) {
      params[i]=0;
      for (j=0;j<nparam;j++)
	params[i]+=mat[i][j]*ax[j];	  
      //fprintf(stdout,"On iteration %2d parameter %2d=%14.5e\n",niter,i,params[i]);
    }
    
#if 0
    /*lets just see how things are going simply for now.*/
    printf("Iter %2d: ",niter);
    for (int k=0;k<nparam;k++)
      printf(" %14.4e ",params[k]);
    printf("\n");
#endif

    if (niter>5)
      converged=1;
    
  }
  psFree(last_params);
  psFree(cov);
  psFree(ax);
  psFree(mat[0]);
  psFree(mat);
  return PS_ERR_NONE;
}




/************************************************************************************************************/
/*!
 *  Fit a set of parameters describing the variance of a set of uncorrelated data.
 */
psErrorCode nkFitNoiseUncorrData(actData *params, actData *data, actData *noise, actData **vecs, long n,  int nparam,actData *like )
{
  int i,iter,finished,converged;
  actData oldlike=1.0,lambda=0,*old_params,*shifts,*deriv, **curve,max_shift;



  //set the first two data points to zero for now, in case the zero points are left in
  data[0]=0;
  data[1]=0;
  
  actData tol=1e-2;
  if (tol<=0) /*haven't set tolerance, so do something sensible*/
    tol=1e-2; /*default to accurate to better than 1e-2 of a sigma.
		convergence should be quadratic, so this ought to be fast.*/
  int max_iter=MB_DEFAULT_MAX_ITER_NOISEFIT;
  
  finished=0;
  converged=0;
  iter=0;
  old_params=vector(nparam);
  shifts=vector(nparam);
  deriv=vector(nparam);
  curve=matrix(nparam,nparam);

  while (finished==0)
    {
      iter++;
      if ((mbGetCurveDerivUncorrData(params,data,noise,vecs,like,deriv,curve,n,nparam))&&(iter==1))
	{
	  //fprintf(stderr,"Error - Bad starting parameters in mbFitNoiseUncorrData.\n");
	  free(old_params);
	  free(shifts);
	  free(deriv);
	  free_matrix(curve);  
	  return PS_ERR_BAD_PARAMETER_VALUE;
	}
      //for (i=0;i<nparam;i++)
      //fprintf(stdout,"deriv[%d]=%14.4e\n",i,deriv[i]);
      if (iter==1)
	oldlike=*like;
      mbGetShiftsUncorrData(nparam, shifts,curve,deriv,*like,oldlike,&lambda,&max_shift);
      for (i=0;i<nparam;i++)
	{
	  if (iter>10)
	    params[i]-=shifts[i];
	  else
	    params[i]-=0.5*shifts[i];  //take it easy on the first few steps in case we would have gone haywire.
	}
      if ((max_shift<tol))
	{
	  converged=1;
	  finished=1;
	}
      if (iter>=max_iter)
	finished=1;     
      oldlike=*like;

      //for (i=0;i<nparam;i++)
      //fprintf(stderr,"parameter %d is %14.4e\n",i,params[i]);
 
    }
  //fprintf(stdout,"Took %d iterations to fit.\n",iter);
  free(shifts);
  free(old_params);
  free(deriv);
  free(curve[0]);
  free(curve);

  if (converged)
    return PS_ERR_NONE;
  else
    return PS_ERR_BAD_PARAMETER_VALUE;
}



/************************************************************************************************************/
/*!
 *  Get the gradient and curvature of the likelihood for a set of data with uncorrelated noise.
 *  We are modelling the variance, assuming that the mean of the data is expected to be zero.
 *  If there is no noise (e.g. we're trying to model it), send in NULL for the noise.
 */
psErrorCode mbGetCurveDerivUncorrData(actData *params, actData *data, actData *noise, actData **vecs, actData *like, actData *deriv, actData **curve, long n, int nparam)
{
  actData *cov, *cix, *trace_deriv, *data_deriv, **trace_curve, **data_curve,**cbci, *tempvec;
  int i,j;  

  cov=vector(n);
  if (noise)
    memcpy(cov,noise,sizeof(actData)*n);
  else
    memset(cov,0,sizeof(actData)*n);
  for (i=0;i<nparam;i++)
    mbDaxpy(n,params[i],vecs[i],1,cov,1);

  /*now take inverse of covariance, failing if something went negative.*/
  for (i=0;i<n;i++) {
    if (cov[i]<=0) {  /*true if we have a bad set of input parameters*/      
      *like=MB_BAD_LIKE;
      free(cov);
      return i;

    }
    else
      cov[i]=1.0/cov[i];      
  }
  
  
  /*these don't strictly need to be allocated as they are used in-place,
    but in case anything acts up, they can be useful to have around for 
    de-bugging purposes.*/
  trace_deriv=vector(nparam);
  data_deriv=vector(nparam);
  trace_curve=matrix(nparam,nparam);
  data_curve=matrix(nparam,nparam);

  tempvec=vector(n);
  cix=vector(n);
  cbci=matrix(nparam,n);

  mbVecMult(n,data,cov,cix);  /*calculate data*cov*/  

  /*calculate the likelihood*/
  *like=0;  
  for (i=0;i<n;i++)
    *like+=log(cov[i]); /*covariance is inverted, so the sign is flipped here*/
  *like=-0.5*mbDot(n,cix,1,data,1)+0.5*(*like);
  //fprintf(stdout,"Likelihood is %14.4f\n",*like);

  

  for (i=0;i<nparam;i++)
    mbVecMult(n,cov,vecs[i],cbci[i]);
  mbVecMult(n,cix,cix,cix);  /*square the data vector*/
  
  for (i=0;i<nparam;i++)
    {
      mbVecMult(n,cix,vecs[i],tempvec);  /*calculate the components of the data deriv*/
      data_deriv[i]=mbVecSum(n,tempvec);
      trace_deriv[i]=mbVecSum(n,cbci[i]);
      deriv[i]=0.5*data_deriv[i]-0.5*trace_deriv[i];
      for (j=i;j<nparam;j++) {
	trace_curve[i][j]=mbDot(n,cbci[i],1,cbci[j],1);
	trace_curve[j][i]=trace_curve[i][j];
	data_curve[i][j]=mbDot(n,tempvec,1,cbci[j],1);
	data_curve[j][i]=data_curve[i][j];
	curve[i][j]=0.5*trace_curve[i][j]-data_curve[i][j];
	curve[j][i]=curve[i][j];
      }
    }
  
  free(cov);
  free_matrix(trace_curve);
  free_matrix(data_curve);
  free(data_deriv);
  free(trace_deriv);
  free(cix);
  free_matrix(cbci);
  free(tempvec);

  return PS_ERR_NONE;
} 


/************************************************************************************************************/
/*!
 *  Do an element-wise multiplication of two vectors and store it in the third.
 *  z=x.*y (in matlab-speak).
 */
void mbVecMult(long n,actData *x,actData *y,actData *z)
{
  long i;
  for (i=0;i<n;i++)
    z[i]=x[i]*y[i];
}

/************************************************************************************************************/
/*!
 *  Sum the elements of a vector
 */
actData mbVecSum(long n,actData *x)
{
  long i;
  actData sum=0;
  for (i=0;i<n;i++)
    {
      sum+=x[i];
      //fprintf(stdout,"Adding %14.5e\n",x[i]);
    }
  return sum;
}


/************************************************************************************************************/
/*!
 *  Figure out where to go for the next step while fitting noise parameters.
 *  Lambda and oldlike will be used in case the fitting hits problems.
 *  In the first instance, they won't be used, but may be down the road.
 */


psErrorCode mbGetShiftsUncorrData(int nparam, actData *shifts, actData **curve,actData *deriv,actData like,actData oldlike,actData *lambda,actData *max_shift)
{
  int i;
  actData frac_shift;  


  memcpy(shifts,deriv,nparam*sizeof(actData));

#ifndef MB_HAVE_LAPACK
  if (nparam==2)
    {
      //fprintf(stdout,"Curvature matrix is:\n%14.4e %14.4e\n%14.4e %14.4e\n",curve[0][0],curve[0][1],curve[1][0],curve[1][1]);
      double a,det;
      det=curve[0][0]*curve[1][1]-curve[0][1]*curve[1][0];
      a=curve[0][0];
      curve[0][0]=curve[1][1]/det;
      curve[1][1]=a/det;
      a=curve[0][1];
      curve[0][1]=-1*curve[1][0]/det;
      curve[1][0]=-1*a/det;
      shifts[0]=(curve[0][0]*deriv[0]+curve[0][1]*deriv[1]);
      shifts[1]=(curve[1][0]*deriv[0]+curve[1][1]*deriv[1]);
      //fprintf(stdout,"Curvature matrix and derivs are:\n%14.4e %14.4e %14.4e\n%14.4e %14.4e %14.4e\n",curve[0][0],curve[0][1],deriv[0],curve[1][0],curve[1][1],deriv[1]);
    }
  else
    {
      fprintf(stderr,"Do not have lapack, so can't do much unless nparam==2.\n");
      return PS_ERR_BAD_PARAMETER_VALUE;
    }
#else  
  char uplo='U';
  int info;
  
  mbPotrf(uplo,nparam,curve[0],nparam,&info);
  if (info)
    {
      fprintf(stderr,"Error in mbGetShiftsUncorrData - curvature matrix is non-positive definite.\n");
      return PS_ERR_BAD_PARAMETER_VALUE;
    }
  mbPotrs(uplo,nparam,1,curve[0],nparam,shifts,nparam,&info);
  if (info)
    {
      fprintf(stderr,"Error in mbGetShiftsUncorrData - potrs failed.\n");
      return PS_ERR_BAD_PARAMETER_VALUE;
    }
  mbPotri(uplo,nparam,curve[0],nparam,&info);
  if (info)
    {
      fprintf(stderr,"Error in mbGetShiftsUncorrData - potri failed.\n");
      return PS_ERR_BAD_PARAMETER_VALUE;
    }
#endif  
  *max_shift=0;
  for (i=0;i<nparam;i++)
    {
      frac_shift=fabs(shifts[i]/sqrt(fabs(curve[i][i])));
      
      if (frac_shift>*max_shift)
	*max_shift=frac_shift;
      //fprintf(stderr,"Max shift is %14.5e\n",*max_shift);
    }
  
  //fprintf(stdout,"Shifts are:\n ");
  //for (i=0;i<nparam;i++)
  //  fprintf(stdout," %14.4e %14.4e\n",shifts[i],sqrt(fabs(curve[i][i])));
  //fprintf(stdout,"\n");

  return PS_ERR_NONE;
}

/*--------------------------------------------------------------------------------*/
void cut_badnoise_dets(mbTOD *tod) 
{
  assert(tod->noise);
  int ncut=0;
  for (int i=0;i<tod->ndet;i++) {
    if (tod->noise->noises[i].converged==0) {
      mbCutsSetAlwaysCut(tod->cuts, tod->rows[i], tod->cols[i]);
      ncut++;
	}
    if (tod->noise->noises[i].params[0]<0)
      mbCutsSetAlwaysCut(tod->cuts, tod->rows[i], tod->cols[i]);    
    if (tod->noise->noises[i].params[1]<0)
      mbCutsSetAlwaysCut(tod->cuts, tod->rows[i], tod->cols[i]);    
  }
  int ngood=0;
  actData *whites=vector(tod->ndet);
  FILE *outfile=fopen("noise_whites.txt","w");
  for (int i=0;i<tod->ndet;i++) 
    if (is_det_used(tod,i)) {
      whites[ngood]=tod->noise->noises[i].params[0];
      fprintf(outfile,"%4d %4d %4d %14.4e\n",i,tod->rows[i],tod->cols[i],whites[ngood]);
      ngood++;
    }
  fclose(outfile);
  actData medval=compute_median(ngood,whites);
  actData medscat=compute_median_scat(ngood,whites);
  actData thresh=2.5;
  int ncut_white=0;

  for (int i=0;i<tod->ndet;i++) {
    if (is_det_used(tod,i)) {
      actData white=tod->noise->noises[i].params[0];
      if ((white<medval-thresh*medscat)||(white>medval+thresh*medscat)) {
	mbCutsSetAlwaysCut(tod->cuts, tod->rows[i], tod->cols[i]);
	ncut_white++;
	
      }
    }
  }

 
  mprintf(stdout,"cut a total of %d detectors due to noise, and %d due to out-of-bounds white level on %s\n",ncut,ncut_white,tod->dirfile);
  free(whites);
}

/*--------------------------------------------------------------------------------*/
void write_ifilter(actComplex *vec, int n, char *fname)
{
  FILE *outfile=fopen(fname,"w");
  if (!outfile)
    return;
  for (int i=0;i<n;i++) 
    fprintf(outfile,"%14.5e %14.5e\n",creal(vec[i]),cimag(vec[i]));
  fclose(outfile);
}
/*--------------------------------------------------------------------------------*/
void CalculateIfilter(mbTOD *tod, int mydet, actComplex *ifilter)
{
  actData **vecs;
  int nparam;
  //int ndata=tod->ndata/2+1;  //since we may want to change this.
  int ndata=fft_real2complex_nelem(tod->ndata);  //since we may want to change this.
  memset(ifilter,0,ndata*sizeof(actComplex));

#if 1
  actData white=tod->ndata*tod->noise->noises[mydet].params[0];
  actData knee_coeff=tod->ndata*tod->noise->noises[mydet].params[1];
  actData delta=get_tod_fft_delta(tod);
  actData powlaw=tod->noise->noises[mydet].powlaw;
  actData fac=pow(delta,powlaw);

  for (int i=1;i<ndata;i++) {  //i starts at one since zero is a special case.
    ifilter[i]=white+pow(i,powlaw)*fac*knee_coeff;
    
  }
  ifilter[0]=2*ifilter[1];


#else  
  nkCalculateOneoverFVecs(tod,&(tod->noise->noises[mydet]),&vecs,&nparam);
  printf("calculated vecs with %d parameters.\n",nparam);
  
  for (int i=0;i<ndata;i++)
    for (int j=0;j<nparam;j++)
      ifilter[i]+=tod->noise->noises[mydet].params[j]*vecs[j][i];
  printf("summed into ifilter.\n");
  ifilter[0]=2*ifilter[1];  //put in something here.
  printf("set zero element.\n");
  free(vecs[0]);
  free(vecs);
#endif
}



/*--------------------------------------------------------------------------------*/
void filter_data_wnoise(mbTOD *tod)
{
  assert(tod->have_data);
  assert(tod->noise);
#pragma omp parallel shared(tod) default(none)
  {    

    actComplex *ifilter=cvector(tod->ndata);
    actComplex *vec=cvector(tod->ndata);
    int nn=fft_real2complex_nelem(tod->ndata);
#pragma omp for schedule(dynamic,1)    
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	CalculateIfilter(tod,i,ifilter);
	act_fftw_execute_dft_r2c(tod->p_forward,tod->data[i],vec);
	for (int j=0;j<nn;j++)
	  vec[j]/=(ifilter[j]*tod->ndata);
	act_fftw_execute_dft_c2r(tod->p_back,vec,tod->data[i]);
      }
    }
    free(vec);
    free(ifilter);
  }

}


/*--------------------------------------------------------------------------------*/
void add_noise_to_tod_new(mbTOD *tod)
{
  assert(tod);
  assert(tod->noise);
  assert(tod->data);
  
  
#pragma omp parallel shared(tod) default(none)
  {
    actComplex *vec=cvector(tod->ndata);
    actComplex *ifilter=(actComplex *)malloc(tod->ndata*sizeof(actComplex));
    int nn=fft_real2complex_nelem(tod->ndata);
#pragma omp for schedule(dynamic,1)    
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	unsigned idum=tod->seed+i+1;
	CalculateIfilter(tod,i,ifilter);
	act_fftw_execute_dft_r2c(tod->p_forward,tod->data[i],vec);
	for (int j=0;j<nn;j++) {
	  actData amp=sqrt(cabs(ifilter[j]));///((actData)tod->ndata);
	  actData re=amp*mygasdev(&idum);
	  actData im=amp*mygasdev(&idum);
	  vec[j]+=re+I*im;
	}
	act_fftw_execute_dft_c2r(tod->p_back,vec,tod->data[i]);
	for (int j=0;j<tod->ndata;j++)
	  tod->data[i][j]/=(actData)tod->ndata;
      }
    }
    free(vec);
    free(ifilter);
  }  
}


/*--------------------------------------------------------------------------------*/
void add_noise_to_tod_gaussian(mbTOD *tod)
{
  assert(tod);
  //assert(tod->noise);
  assert(tod->data);
  
  
#pragma omp parallel shared(tod) default(none)
  {
#pragma omp single
    printf("Using %d threads to add gaussian noise.\n",omp_get_num_threads());
    //schedule(dynamic,1)    
#pragma omp for 
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	unsigned idum=tod->seed+i+1;
	for (int j=0;j<tod->ndata;j++) {
	  tod->data[i][j]=mygasdev(&idum);
	}  
      }
    }
  }
}
/*--------------------------------------------------------------------------------*/

actComplex *nkMCEButterworth(mbTOD *tod)
{
  actComplex b_1_1=(-2.0*32092.0/pow(2,15));
  actComplex b_1_2=( 2.0*15750.0/pow(2,15));
  actComplex b_2_1=(-2.0*31238.0/pow(2,15));
  actComplex b_2_2=( 2.0*14895.0/pow(2,15));
  actData f_samp=1./(0.00000002*100.*33.);

  actData mynorm=4/(1+b_1_1+b_1_2)*4/(1+b_2_1+b_2_2);
  actComplex *filt=cvector(tod->ndata);


  for (int i=0;i<tod->ndata;i++) {
    actComplex omega=((actData)i)/(tod->deltat*((actData)tod->ndata))/f_samp*2.0*M_PI;
    actComplex h1_omega=(1+2*cexp(-I*omega)+cexp(-2*I*omega))/(1+b_1_1*cexp(-I*omega)+b_1_2*cexp(-2*I*omega));
    actComplex h2_omega=(1+2*cexp(-I*omega)+cexp(-2*I*omega))/(1+b_2_1*cexp(-I*omega)+b_2_2*cexp(-2*I*omega));
    filt[i]=h1_omega*h2_omega/mynorm;

  }

  return filt;
    
}

/*--------------------------------------------------------------------------------*/
void nkDeButterworth(mbTOD *tod)
{
  assert(tod);
  assert(tod->data);


#if 0

#pragma omp parallel shared(tod) default(none)
  {
    actComplex *vec=(actComplex *)malloc(tod->ndata*sizeof(actComplex));
    actComplex *filt=nkMCEButterworth(tod);
    int nn=fft_real2complex_nelem(tod->ndata);
#pragma omp for schedule(dynamic,1) 
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	act_fftw_execute_dft_r2c(tod->p_forward,tod->data[i],vec);
	for (int j=0;j<nn;j++)
	  vec[j]/=(filt[j]*tod->ndata);
	act_fftw_execute_dft_c2r(tod->p_back,vec,tod->data[i]);
      }
    }

    free(vec);
    free(filt);
  }
#else

    actComplex **data_ft=fft_all_data(tod);
#pragma omp parallel shared(tod,data_ft) default(none)
    {
      int nn=fft_real2complex_nelem(tod->ndata);
      actComplex *filt=nkMCEButterworth(tod);
#pragma omp for
      for (int i=0;i<tod->ndet;i++)
	if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	  for (int j=0;j<nn;j++)
	    data_ft[i][j]/=filt[j];  //no normalization required as the all_data transforms do it.	  
	}
      free(filt);
    }
    
    ifft_all_data(tod,data_ft);
    free(data_ft[0]);
    free(data_ft);
#endif

  
}

/*--------------------------------------------------------------------------------*/
void nkReButterworth(mbTOD *tod)
{
  assert(tod);
  assert(tod->data);

#pragma omp parallel shared(tod) default(none)
  {
    actComplex *vec=(actComplex *)malloc(tod->ndata*sizeof(actComplex));
    actComplex *filt=nkMCEButterworth(tod);
    int nn=fft_real2complex_nelem(tod->ndata);
#pragma omp for schedule(dynamic,1) 
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	act_fftw_execute_dft_r2c(tod->p_forward,tod->data[i],vec);
	for (int j=0;j<nn;j++)
	  vec[j]*=(filt[j]/tod->ndata);
	act_fftw_execute_dft_c2r(tod->p_back,vec,tod->data[i]);
      }
    }

    free(vec);
    free(filt);
  }
  
}

/*--------------------------------------------------------------------------------*/
void nkDeconvolveTimeConstants(mbTOD *tod)
{
  assert(tod);
  assert(tod->data);
  assert(tod->time_constants);
#if 0
#pragma omp parallel shared(tod) default(none)
  {
    actComplex *vec=(actComplex *)malloc(tod->ndata*sizeof(actComplex));
    actData *freqs=get_freq_vec(tod);
    int nn=fft_real2complex_nelem(tod->ndata);
#pragma omp for schedule(dynamic,1) 
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	act_fftw_execute_dft_r2c(tod->p_forward,tod->data[i],vec);
	//if ((tod->rows[i]==0)&&(tod->cols[i]==2))
	//printf("Time constant on %d %d is %14.5g\n",tod->rows[i],tod->cols[i],tod->time_constants[tod->rows[i]][tod->cols[i]]);
	actComplex myfac=2*M_PI*I*tod->time_constants[tod->rows[i]][tod->cols[i]];
	for (int j=0;j<nn;j++)
	  vec[j]*=(1+myfac*freqs[j])/((actData)tod->ndata);
	act_fftw_execute_dft_c2r(tod->p_back,vec,tod->data[i]);
      }
    }

    free(vec);
    free(freqs);
  }
#else
  actComplex **data_ft=fft_all_data(tod);
#pragma omp parallel shared(tod,data_ft) default(none)
  {
    int nn=fft_real2complex_nelem(tod->ndata);
    actData *freqs=get_freq_vec(tod);
#pragma omp for 
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
        actComplex myfac=2*M_PI*I*tod->time_constants[tod->rows[i]][tod->cols[i]];
        for (int j=0;j<nn;j++)
          data_ft[i][j]*=(1+myfac*freqs[j]);
	
      }
    }
    free(freqs);
  }
  ifft_all_data(tod,data_ft);
  free(data_ft[0]);
  free(data_ft);
  
#endif  
}
/*--------------------------------------------------------------------------------*/
void nkReconvolveTimeConstants(mbTOD *tod)
{
  assert(tod);
  assert(tod->data);
  assert(tod->time_constants);

#pragma omp parallel shared(tod) default(none)
  {
    actComplex *vec=(actComplex *)malloc(tod->ndata*sizeof(actComplex));
    actData *freqs=get_freq_vec(tod);
    int nn=fft_real2complex_nelem(tod->ndata);
#pragma omp for schedule(dynamic,1) 
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	act_fftw_execute_dft_r2c(tod->p_forward,tod->data[i],vec);
	//if ((tod->rows[i]==0)&&(tod->cols[i]==2))
	//printf("Time constant on %d %d is %14.5g\n",tod->rows[i],tod->cols[i],tod->time_constants[tod->rows[i]][tod->cols[i]]);
	actComplex myfac=2*M_PI*I*tod->time_constants[tod->rows[i]][tod->cols[i]];
	for (int j=0;j<nn;j++)
	  vec[j]*=1.0/(1+myfac*freqs[j])/((actData)tod->ndata);
	act_fftw_execute_dft_c2r(tod->p_back,vec,tod->data[i]);
      }
    }

    free(vec);
    free(freqs);
  }
  
}
/*--------------------------------------------------------------------------------*/
void destroy_noise_raw(mbNoiseVectorStruct *noise)
{
  free(noise->noises);
  free(noise);
}
/*--------------------------------------------------------------------------------*/
void destroy_tod_noise(mbTOD *tod)
{
  assert(tod);
  assert(tod->noise);
  destroy_noise_raw(tod->noise);
  tod->noise=NULL;
  
}

/*--------------------------------------------------------------------------------*/
int get_nn(int n) 
//Get the length of a real2complex fft
{
  return (n)/2+1;
}
/*--------------------------------------------------------------------------------*/
actData *get_freq_vec(mbTOD *tod)
{
  //int nn=1+tod->ndata/2;
  int nn=get_nn(tod->ndata);
  actData *vec=(actData *)malloc(sizeof(actData)*nn);
  
  actData dnu=1.0/(tod->deltat*(actData)tod->ndata);
  for (int i=0;i<nn;i++)
    vec[i]=dnu*(actData)i;
  vec[0]=0.5*dnu;
  
  return vec;
}


/*--------------------------------------------------------------------------------*/
void allocate_tod_noise_bands(mbTOD *tod,actData *bands, int nband)
{
  tod->band_noise=(mbNoiseVectorStructBands *)calloc(sizeof(mbNoiseVectorStructBands),1);
  tod->band_noise->ndet=tod->ndet;
  tod->band_noise->deltat=tod->deltat;
  
  tod->band_noise->nband=nband;
  tod->band_noise->bands=(actData *)calloc(sizeof(actData)*(nband+1),1);
  memcpy(tod->band_noise->bands,bands,sizeof(actData)*(nband+1));
  tod->band_noise->ibands=ivector(nband+1);
  
  
  actData *freqs=get_freq_vec(tod);
  int nn=get_nn(tod->ndata);

  for (int j=0;j<nband+1;j++) {
    tod->band_noise->ibands[j]=0;
    for (int i=0;i<nn;i++)
      if (freqs[i]<tod->band_noise->bands[j])
	tod->band_noise->ibands[j]=i;	
  }
  if (tod->band_noise->ibands[nband]==nn-1)
    tod->band_noise->ibands[nband]++;
  
  



  tod->band_noise->rot_mats=(actData ***)calloc(sizeof(actData **)*nband,1);
  tod->band_noise->inv_rot_mats_transpose=(actData ***)calloc(sizeof(actData **)*nband,1);
  
  tod->band_noise->noise_params=(mbNoiseParams1PixBand **)malloc(sizeof(mbNoiseParams1PixBand *)*nband);
  for (int i=0;i<nband;i++)
    tod->band_noise->noise_params[i]=(mbNoiseParams1PixBand *)calloc(sizeof(mbNoiseParams1PixBand)*tod->band_noise->ndet,1);
  tod->band_noise->do_rotations=(bool *)malloc(sizeof(bool)*nband);
  for (int i=0;i<nband;i++)
    tod->band_noise->do_rotations[i]=false;
  
}
/*--------------------------------------------------------------------------------*/

actComplex **fft_all_data_flag(mbTOD *tod, unsigned flags) 
{
  assert(tod);
  assert(tod->have_data);
  assert(tod->ndet>0);

  int nn=get_nn(tod->ndata);
  actComplex **data_fft=cmatrix(tod->ndet,nn);

  int *n=ivector(tod->ndet);
  for (int i=0;i<tod->ndet;i++)
    n[i]=tod->ndata;
  

  fftw_plan plan=fftw_plan_many_dft_r2c(1,n,tod->ndet,tod->data[0],NULL,1,tod->ndata,data_fft[0],NULL,1,nn,flags);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  free(n);

  return data_fft;
  
}
/*--------------------------------------------------------------------------------*/

actComplex **fft_all_data(mbTOD *tod) 
{
  actComplex **data_fft=fft_all_data_flag(tod,FFTW_ESTIMATE);
  return data_fft;
  
}
/*--------------------------------------------------------------------------------*/

void ifft_all_data_flag(mbTOD *tod,actComplex **data_fft,unsigned flag) 
{
  assert(tod);
  assert(tod->have_data);
  assert(tod->ndet>0);
  
  int nn=get_nn(tod->ndata);

  int *n=ivector(tod->ndet);
  for (int i=0;i<tod->ndet;i++)
    n[i]=tod->ndata;

  
#if 0
#pragma omp parallel
#pragma omp single
  fftw_plan_with_nthreads(omp_get_num_procs());
#endif
  
  fftw_plan plan=fftw_plan_many_dft_c2r(1,n,tod->ndet,data_fft[0],NULL,1,nn,tod->data[0],NULL,1,tod->ndata,flag);
  fftw_execute(plan);  
  fftw_destroy_plan(plan);
  

  actData fn=tod->ndata;
#pragma omp parallel for shared(tod,fn) default(none)  
  for (int i=0;i<tod->ndet;i++)
    for (int j=0;j<tod->ndata;j++)
      tod->data[i][j]/=fn;

  
  
  free(n);
  return;
  
}

/*--------------------------------------------------------------------------------*/

void ifft_all_data(mbTOD *tod,actComplex **data_fft) 
{
#if 1
  ifft_all_data_flag(tod,data_fft,FFTW_ESTIMATE);
#else
  assert(tod);
  assert(tod->have_data);
  assert(tod->ndet>0);
  
  int nn=get_nn(tod->ndata);

  int *n=ivector(tod->ndet);
  for (int i=0;i<tod->ndet;i++)
    n[i]=tod->ndata;

  
  fftw_plan plan=fftw_plan_many_dft_c2r(1,n,tod->ndet,data_fft[0],NULL,1,nn,tod->data[0],NULL,1,tod->ndata,FFTW_ESTIMATE);
  fftw_execute(plan);  
  fftw_destroy_plan(plan);
  

  actData fn=tod->ndata;
#pragma omp parallel for shared(tod,fn) default(none)  
  for (int i=0;i<tod->ndet;i++)
    for (int j=0;j<tod->ndata;j++)
      tod->data[i][j]/=fn;

  
  
  free(n);
#endif
  return;
  
}
/*--------------------------------------------------------------------------------*/
void act_syrk(char uplo, char trans, int n, int m, actData alpha, actData *a, int lda, actData beta, actData *b, int ldb)
{
#ifdef ACTDATA_DOUBLE
  clapack_dsyrk(uplo,trans,n,m,alpha,a,lda,beta,b,ldb);
#else
  clapack_ssyrk(uplo,trans,n,m,alpha,a,lda,beta,b,ldb);
#endif

}
/*--------------------------------------------------------------------------------*/

actData **get_banded_correlation_matrix_from_fft(mbTOD *tod, actComplex **data_fft, actData nu_min, actData nu_max)
{
  assert(tod);
  assert(sizeof(actComplex *)==sizeof(actData *));
  int nn=get_nn(tod->ndata);
  actData **mat=matrix(tod->ndet,tod->ndet);

  memset(mat[0],0,sizeof(actData)*tod->ndet*tod->ndet);

  actData *freqs=get_freq_vec(tod);

  int imin=0;
  int imax=0;

  for (int i=0;i<nn;i++) {
    if (freqs[i]<nu_min)
      imin=i;
    if (freqs[i]<nu_max)
      imax=i;
  }
  
  act_syrk('u','t',tod->ndet,2*(imax-imin+1),1.0,(double *)(&data_fft[0][imin]),2*nn,0.0,mat[0],tod->ndet);
  for (int i=0;i<tod->ndet;i++)
    for (int j=0;j<i;j++)
      mat[j][i]=mat[i][j];


  free(freqs);
  return mat;

}

/*--------------------------------------------------------------------------------*/

actData **get_banded_correlation_matrix_from_fft_int(mbTOD *tod, actComplex **data_fft, int imin, int imax)
{
  assert(tod);
  assert(sizeof(actComplex *)==sizeof(actData *));
  int nn=get_nn(tod->ndata);
  actData **mat=matrix(tod->ndet,tod->ndet);

  memset(mat[0],0,sizeof(actData)*tod->ndet*tod->ndet);
  
  act_syrk('u','t',tod->ndet,2*(imax-imin),1.0,(double *)(&data_fft[0][imin]),2*nn,0.0,mat[0],tod->ndet);
  for (int i=0;i<tod->ndet;i++)
    for (int j=0;j<i;j++)
      mat[j][i]=mat[i][j];
  return mat;

}


/*--------------------------------------------------------------------------------*/
void get_eigenvectors(actData **mat, int n) 
//transform mat into its eigenvectors.
{
  
  actData fwork;
  int liwork,lwork,iiwork;
  char jobz='v';
  char uplo='u';
  
  int info;

  actData *w=vector(n);
  
  lwork=-1;
  liwork=-1;
#ifdef ACTDATA_DOUBLE
  dsyevd_(&jobz,&uplo,&n,mat[0],&n,w,&fwork,&lwork,&iiwork,&liwork,&info);
#else
  ssyevd_(&jobz,&uplo,&n,mat[0],&n,w,&fwork,&lwork,&iiwork,&liwork,&info);
#endif
  lwork=fwork;
  liwork=iiwork;

  actData *work=vector(lwork);
  int *iwork=ivector(liwork);

  //printf("sizes in eigs are %d %d\n",lwork,liwork);
  
#ifdef ACTDATA_DOUBLE
  dsyevd_(&jobz,&uplo,&n,mat[0],&n,w,work,&lwork,iwork,&liwork,&info);
#else
  ssyevd_(&jobz,&uplo,&n,mat[0],&n,w,work,&lwork,iwork,&liwork,&info);
#endif
  
  free(w);
  free(work);
  free(iwork);
  
  return;
}

/*--------------------------------------------------------------------------------*/
bool do_I_have_rotations(mbTOD *tod)
{
  bool rotated=false;
  mbNoiseVectorStructBands *noise=tod->band_noise;
  for (int band=0;band<noise->nband;band++) 
    if (noise->do_rotations[band])
      rotated=true;
  return rotated;
}
/*--------------------------------------------------------------------------------*/
actComplex **apply_banded_rotations(mbTOD *tod, actComplex **mat_in, bool do_forward)
{
  assert(tod);
  assert(tod->band_noise);
  mbNoiseVectorStructBands *noise=tod->band_noise;

  int nn=get_nn(tod->ndata);
  actComplex **mat=cmatrix(tod->ndet,nn);
  memcpy(mat[0],mat_in[0],sizeof(actComplex)*nn*tod->ndet);

  for (int band=0;band<noise->nband;band++) {
    if (noise->do_rotations[band]) {
      //printf("doing rotation on band %d\n",band);
      actData **rotmat;
      char trans;
      if (do_forward) {
	trans='n';
	rotmat=noise->rot_mats[band];
      }
      else {
	trans='t';
	rotmat=noise->inv_rot_mats_transpose[band];
      }
      act_gemm('n',trans,2*(noise->ibands[band+1]-noise->ibands[band]),tod->ndet,tod->ndet,1.0,(actData *)(&(mat_in[0][noise->ibands[band]])),2*nn,rotmat[0],tod->ndet,0.0,(actData *)(&(mat[0][noise->ibands[band]])),2*nn);
    }
  }
  
  return mat;
}

/*--------------------------------------------------------------------------------*/
int fit_banded_noise_1det_full(mbNoiseParams1PixBand *params, actData *dat)
{
  int n=params->i_high-params->i_low;
  params->noise_data=vector(n);
  if (!params->noise_data)
    return -1;  //only should be  able to fail if malloc dies.
  for (int i=0;i<n;i++) {
    int imin=i-NOISE_FIT_WIDTH;
    int imax=i+NOISE_FIT_WIDTH;
    if (imin<0)
      imin=0;
    if (imax>n)
      imax=n;
    actData mymax=dat[imin];
    for (int ii=imin+1;ii<imax;ii++)
      if (dat[ii]>mymax)
	mymax=dat[ii];
    params->noise_data[i]=1.0/mymax/mymax;
  }
  return 0;
}
/*--------------------------------------------------------------------------------*/
void scale_banded_noise_band( mbTOD *tod,int which_band,actData fac)
//scale the weight in a band by factor fac.  
{
  assert(tod);
  assert(tod->band_noise);
  mbNoiseVectorStructBands *noise=tod->band_noise;
  for (int i=0;i<tod->ndet;i++){
    switch(noise->noise_params[which_band][i].noise_type) {
    case MBNOISE_CONSTANT: 
    case MBNOISE_INTERP: 
      noise->noise_params[which_band][i].noise_data[0]*fac;
      break;     
    case MBNOISE_FULL: {
      mbNoiseParams1PixBand *nn=&(noise->noise_params[which_band][i]);
      for (int j=0;j<nn->i_high-nn->i_low;j++)
	nn->noise_data[j]*=fac;
      break;
    }
    default:
      assert(1==0);  //break if we don't know the type.
      break;
    }
  }
}


/*--------------------------------------------------------------------------------*/
int fit_banded_noise_1det_constant( mbNoiseParams1PixBand *params, actData *dat)
{
  actData tot=0;
  int n=params->i_high-params->i_low;
  for (int i=0;i<n;i++)
    tot+=dat[i]*dat[i];
  params->noise_data=vector(1);
  params->white=sqrt(tot/((actData)n));
  //params->noise_data[0]=1.0/params->white/params->white;
  params->noise_data[0]=((actData)n)/tot;
  if (isfinite(tot))
    return 0;
  else
    return -1;
}
/*--------------------------------------------------------------------------------*/
int fit_banded_noise_1det( mbNoiseParams1PixBand *params, actComplex *dataft)
{
  int n=params->i_high-params->i_low;
  actData *myabs=vector(n);

  for (int i=0;i<n;i++)
    myabs[i]=cabs(dataft[i+params->i_low]);

  switch(params->noise_type) {
  case MBNOISE_CONSTANT:
  case MBNOISE_INTERP:
    if (fit_banded_noise_1det_constant(params,myabs)) {
      free(myabs);
      return -1;
    }
    break;
  case MBNOISE_FULL:
    if (fit_banded_noise_1det_full(params,myabs)) {
      free(myabs);
      return -1;
    }
    break;
  default:
    printf("Unrecognized type in fitting noise.\n");
    break;
    
  }
  
  free(myabs);
  return 0;

}


/*--------------------------------------------------------------------------------*/

void get_simple_banded_noise_model(mbTOD *tod, bool *do_rots, mbNoiseType *types)
{
  assert(tod);
  assert(tod->have_data);
  assert(tod->band_noise);
  mbNoiseVectorStructBands *noise=tod->band_noise;


  actComplex **data_fft=fft_all_data(tod);
  
  for (int band=0;band<noise->nband;band++) {
    noise->do_rotations[band]=do_rots[band];
    if (do_rots[band]) {
      actData **corrmat=get_banded_correlation_matrix_from_fft_int(tod,data_fft,noise->ibands[band],noise->ibands[band+1]);
      get_eigenvectors(corrmat,tod->ndet);
      noise->rot_mats[band]=corrmat;
      noise->inv_rot_mats_transpose[band]=corrmat;
    }
  }
  
  actComplex **data_rot=apply_banded_rotations(tod, data_fft, true);
  
  free(data_fft[0]);
  free(data_fft);

  //printf("fitting noise now.\n");
  for (int band=0;band<noise->nband;band++) {
#pragma omp parallel for shared(band,tod,noise,data_rot,types) default(none)
    for (int i=0;i<tod->ndet;i++) {
      noise->noise_params[band][i].noise_type=types[band];
      noise->noise_params[band][i].i_low=noise->ibands[band];
      noise->noise_params[band][i].i_high=noise->ibands[band+1];
      fit_banded_noise_1det(&(noise->noise_params[band][i]),data_rot[i]);
    }
  }
  
  free(data_rot[0]);
  free(data_rot);
}
/*--------------------------------------------------------------------------------*/
void apply_banded_noise_1det_full(mbNoiseParams1PixBand *params,actComplex *dat)
{
  int n=params->i_high-params->i_low;
  for (int i=0;i<n;i++)
    dat[i+params->i_low]*=params->noise_data[i];
  
}

/*--------------------------------------------------------------------------------*/
void apply_banded_noise_1det_interp(mbNoiseParams1PixBand *params_left, mbNoiseParams1PixBand *params, mbNoiseParams1PixBand *params_right, actComplex *dat)
{

  int n=params->i_high - params->i_low;
  int nn=n/2;


  for (int i=0;i<nn;i++) {
    actData fac=((actData)i)/((actData) nn);
    dat[i+params->i_low]*=(1.0-fac)*params_left->noise_data[0]+fac*params->noise_data[0];
  }
  for (int i=nn;i<n;i++) {
    actData fac=((actData)(i-nn))/((actData)( n-nn));
    dat[i+params->i_low]*=fac*params_right->noise_data[0]+(1-fac)*params->noise_data[0];
  }
  
}

/*--------------------------------------------------------------------------------*/
void apply_banded_noise_1det_constant(mbNoiseParams1PixBand *params,actComplex *dat)
{
  int n=params->i_high-params->i_low;
  for (int i=0;i<n;i++)
    dat[i+params->i_low]*=params->noise_data[0];

}
/*--------------------------------------------------------------------------------*/
void apply_banded_noise_complex(mbTOD *tod, actComplex **dat)
{

  mbNoiseVectorStructBands *noise=tod->band_noise;

#pragma omp parallel for shared(tod,dat,noise) default(none)
  for (int det=0;det<tod->ndet;det++) {
    int fwee;      
    for (int i=0;i<noise->nband;i++) 
      
      switch(noise->noise_params[i][det].noise_type) {
      case MBNOISE_INTERP:
	fwee=0;
	int ileft=i-1;
	if (ileft<0)
	  ileft=0;
	int iright=i+1;
	if (iright>=noise->nband)
	  iright=noise->nband-1;
	apply_banded_noise_1det_interp(&(noise->noise_params[ileft][det]),&(noise->noise_params[i][det]),&(noise->noise_params[iright][det]),dat[det]);
	break;

      case MBNOISE_FULL:
	apply_banded_noise_1det_full(&(noise->noise_params[i][det]),dat[det]);
	break;
      case MBNOISE_CONSTANT:
	apply_banded_noise_1det_constant(&(noise->noise_params[i][det]),dat[det]);
	break;
      default:
	printf("Warning - unrecognized noise type in apply_banded_noise_complex.\n");
	break;	
      }
  }
}
/*--------------------------------------------------------------------------------*/
void apply_banded_noise_model(mbTOD *tod)
{

  assert(tod->band_noise);
  
  if (do_I_have_rotations(tod)) {
    
    actComplex **data_ft=fft_all_data(tod);
    actComplex **data_rot=apply_banded_rotations(tod,data_ft,true);
    free(data_ft[0]);
    free(data_ft);
    
    apply_banded_noise_complex(tod,data_rot);
    actComplex **data_back=apply_banded_rotations(tod,data_rot,false);
    free(data_rot[0]);
    free(data_rot);
    
    ifft_all_data(tod,data_back);
    free(data_back[0]);
    free(data_back);
  }
  else
    {
      actComplex **data_ft=fft_all_data(tod);
      apply_banded_noise_complex(tod,data_ft);
      ifft_all_data(tod,data_ft);
      free(data_ft[0]);
      free(data_ft);
    }
  
}
