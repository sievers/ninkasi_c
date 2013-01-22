////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file common_mode.c
/// Defines structure and operations for computing and using common-mode signals.
/// Author: Jon Sievers, CITA.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

//#include "act.h"
//#include "mbUtils.h"
//#include "mbErrorCodes.h"
#include "mbCommon.h"
#include "mbCuts.h"
#include <math.h>
#include "sys/time.h"
//#include "cblas.h"
//#include "clapack.h"
//#if !defined(MB_SKIP_OMP)
#  include <omp.h>
//#endif
#include "mbTOD.h"
#include <assert.h>
#include <ninkasi.h>
#include <string.h>
#ifndef _SKIP_MKL
#include <mkl.h>
#else
//#include <cblas.h>
#endif
#include <nk_clapack.h>
#include "ninkasi_mathutils.h"


size_t mbFreadWrite(bool dowrite, void *ptr, size_t sz, long nobj, FILE *stream) 
{
  if (dowrite)
    return freadwrite(ptr,sz,nobj,stream,1);
  else
    return freadwrite(ptr,sz,nobj,stream,0);
}
/*---------------------------------------------------------------------------------------------------------*/
/// FORTRAN subroutine from BLAS.  See http://www.netlib.org/blas/sgemv.f
/// Computes y = alpha*A*x + beta*y  where A is a matrix; alpha,beta are scalars; x,y are vectors.

/// Wraps a BLAS call.

static void csgemv(const char trans, const int m, const int n, const float alpha, const float *a,
                   const int lda, const float *x, const int incx, const float beta, const float *y,
                   const int incy)
{
  sgemv_(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy,1);
}


static void cdgemv(const char trans, const int m, const int n, const double alpha, const double *a,
                   const int lda, const double *x, const int incx, const double beta, const double *y,
                   const int incy)
{
  dgemv_(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy,1);
}





/*---------------------------------------------------------------------------------------------------------*/
/// Free a ptr to an array of ptrs, along with what it points to.
/// \param p  The ptr to ptrs.

static inline void psFree2d(actData **p) {
  if (p)
    psFree(*p);
  psFree(p);
}

/*---------------------------------------------------------------------------------------------------------*/


void mbNoiseCommonModeFree( mbNoiseCommonMode *cm )
{
  if (cm->have_data) {
    cm->have_data=false;
    psFree2d(cm->data);
    cm->data=NULL;
  }
  if (cm->have_vecs) {
    cm->have_vecs=false;
    psFree2d(cm->vecs);
    cm->vecs=NULL;
  }
  if (cm->have_ata) {
    cm->have_ata=false;
    psFree2d(cm->ata);
    cm->ata=NULL;
  }

  psFree2d(cm->atx);
  cm->atx=NULL;
  psFree(cm->weights);
  cm->weights=NULL;
  cm->have_weights=false;

  psFree(cm->calbol_start);
  cm->calbol_start=NULL;
  psFree(cm->calbol_stop);
  cm->calbol_stop=NULL;
  psFree(cm->common_mode);
  cm->common_mode=NULL;
  psFree(cm->common_mode_nocalbol);
  cm->common_mode_nocalbol=NULL;
  psFree(cm->median_vals);
  cm->median_vals=NULL;
  psFree(cm->median_scats);
  cm->median_scats=NULL;
  psFree2d(cm->fit_params);
  cm->fit_params=NULL;
  psFree(cm->frac_errs);
  cm->frac_errs=NULL;
  psFree(cm->unsmooth_ratio);
  cm->unsmooth_ratio=NULL;
  psFree(cm->icutvec);
  cm->icutvec=NULL;
}


/*---------------------------------------------------------------------------------------------------------*/
/// Constructor for the common mode data object.
/// \param tod   The TOD being studied.
/// \param np_common The order of the ???
/// \param np_poly   The order of the polynomial fit to ??

mbNoiseCommonMode *mbAllocateCommonMode(const mbTOD *tod, int np_common, int np_poly)
{
  assert (tod != NULL);
  mbNoiseCommonMode *data=(mbNoiseCommonMode *)psAlloc(sizeof(mbNoiseCommonMode));  
  memset(data,0,sizeof(mbNoiseCommonMode));
  data->ndet=tod->ndet;
  data->ndata=tod->ndata;
  if (np_common<0)
    data->np_common=MB_NOISE_DEFAULT_NP_COMMON;
  else
    data->np_common=np_common;
  
  if (np_poly<0)
    data->np_poly=MB_NOISE_DEFAULT_NP_POLY;
  else
    data->np_poly=np_poly;
  data->ncalbol=0;
  data->icutvec=NULL;
  data->nparam=data->np_poly+data->np_common+data->ncalbol;
  
  data->common_mode=(actData *)psAlloc(data->ndata*sizeof(actData));
  data->common_mode_nocalbol=(actData *)psAlloc(data->ndata*sizeof(actData));
  data->median_vals=(actData *)psAlloc(data->ndet*sizeof(actData));
  data->median_scats=(actData *)psAlloc(data->ndet*sizeof(actData));
  data->fit_params=psAllocMatrix(data->nparam,data->ndet);
  data->calbol_pad=0.5;
  data->keep_vecs=true;
  data->data=NULL;
  data->vecs=NULL;
  data->weights=NULL;
  data->have_data=false;
  data->have_vecs=false;
  data->have_weights=false;
  data->frac_errs=(actData *)calloc(data->ndet,sizeof(actData));
  

  data->apply_common=false;
  data->common_is_applied=false;

  data->ata=NULL;
  data->have_ata=false;
  data->atx=NULL;
  data->have_atx=false;
  
  data->nsig=10.0;      // Glitch finding paramters
  data->tGlitch=1.5;   // These are set such that we find calbols
  data->tSmooth=5.0;   // t_smoooooth! Spin me some records, yo.
  data->unsmooth_ratio=(actData *)psAlloc(data->ndet*sizeof(actData));
  data->have_unsmoothed_ratio=false;

  //data->dt=(tod->dt[tod->ndata-1]-tod->dt[0])/((actData)tod->ndata);
  data->dt=tod->deltat;
 
#ifdef HAVE_FREE_FUNC
  psMemSetDeallocator( data, (psFreeFunc) mbNoiseCommonModeFree ); 
#endif
  return data;
}


/// Remove any stored data and free relevant vectors, so the object can be re-used or freed.
/// \param data  The object to clean up.

void mbCleanUpCommonMode(mbNoiseCommonMode *data)
{
  if (data->have_data) {
    psFree2d(data->data);
    data->data=NULL;
    data->have_data=false;	  
  }

  if (data->have_vecs) {
    psFree2d(data->vecs);
    data->vecs=NULL;
    data->have_vecs=false;	
  }

  if (data->have_ata) {
    psFree2d(data->ata);    
    data->ata=NULL;
    data->have_ata=false;
  }
}



/*---------------------------------------------------------------------------------------------------------*/
/// Interface to BLAS subroutine SGEMM.  See http://www.netlib.org/blas/sgemm.f
/// Performs a matrix-matrix multiply and add operation.

static void sgemm_simple(int k, int m, int n,actData *a,int alen_is_n, actData *b,int blen_is_n, actData *c)
{
#if 1
  char transa,transb;
  int alen,blen;
  if (alen_is_n) {
    transa='N';
    alen=n;
  } else {
    transa='T';
    alen=k;
  }


  if (blen_is_n) {
    transb='T';
    blen=n;
  } else {
    transb='N';
    blen=m; 
  }


  act_gemm(transa,transb,m,n,k,1.0,a,alen,b,blen,0.0,c,m);
  //  actData done=1.0;
  //  actData dzero=0.0;
  //#ifdef ACTDATA_DOUBLE
  //  dgemm_(&transa,&transb,&m,&n,&k,&done,a,&alen,b,&blen,&dzero,c,&m);
  //#else
  //  sgemm_(&transa,&transb,&m,&n,&k,&done,a,&alen,b,&blen,&dzero,c,&m);
  //#endif
  
#else
  //int transa,alen, transb,blen;
  int alen, blen;
#ifdef _SKIP_MKL
  enum
#endif
    CBLAS_TRANSPOSE transa, transb;
  if (alen_is_n) {
    transa=CblasNoTrans;
    alen=n; 
  } else {
    transa=CblasTrans;
    alen=k; 
  }

  if (blen_is_n) {
    transb=CblasTrans;
    blen=n;
  } else {
    transb=CblasNoTrans;
    blen=m; 
  }

#ifdef ACTDATA_DOUBLE
  cblas_dgemm(CblasRowMajor,transa,transb,k,m,n,1.0,a,alen,b,blen,0.0,c,m);
#else
  cblas_sgemm(CblasRowMajor,transa,transb,k,m,n,1.0,a,alen,b,blen,0.0,c,m);
#endif
#endif
}


/*---------------------------------------------------------------------------------------------------------*/
/// Compute and return the mean of a data vector.
/// \param vec  Data to be averaged.
/// \param n    Length of vec.
/// \return     The mean value.

static inline actData array_mean(const actData *vec, int n)
{
  actData sum=0;
  for (int i=0;i<n;i++)
    sum+=vec[i];
  return sum/((actData)n);
}


/*---------------------------------------------------------------------------------------------------------*/
#if 0
/// Allocate and fill a structure holding a copy of all ndet data sets from a TOD.
/// \param tod  The TOD to copy.
/// \return An mbDataCopy containing a copy of the TOD's data.

mbDataCopy *mbCopyData(const mbTOD *tod)
{
  mbDataCopy *data;
  assert(tod->ndet!=0);
  assert(tod->ndata!=0);

  data=(mbDataCopy *)psAlloc(sizeof(mbDataCopy));
  data->ndet=tod->ndet;
  data->ndata=tod->ndata;
  data->data=psAllocMatrix(tod->ndet,tod->ndata);
  memcpy(data->data[0],tod->data[0],tod->ndet*tod->ndata*sizeof(actData));

  return data;
}
#endif



/*---------------------------------------------------------------------------------------------------------*/
/// Make a scratch copy of the raw data, but rebias each detector's TOD by subtracting that
/// detector's median data value, and then rescale so that each detector has unit "scatter", meaning
/// its median absolute value is one.
/// \param tod  The TOD to study.
/// \param fit  The common mode object being updated.

void mbRescaleArray(const mbTOD *tod, mbNoiseCommonMode *fit)
{
  
  assert(tod!=NULL);
  assert(fit!=NULL);
  // If you fail ndet>1, it's because you only sent in one detector, which isn't much of a common mode.
  assert(tod->ndet>1);  
  assert(tod->ndata>1);
  assert(tod->ndata==fit->ndata);
  assert(tod->ndet==fit->ndet);

  mbTimeValue ticker;
  mbStartTime(&ticker);

  
  bool had_data=false;
  if (fit->have_data)
    had_data=true;
  else {
    fit->data=psAllocMatrix(fit->ndet,fit->ndata);
    fit->have_data=true;
  }
  
  // Make a scratch copy of the full TOD 
  memcpy(fit->data[0],tod->data[0],tod->ndet*tod->ndata*sizeof(actData));
  
#if !defined(MB_SKIP_OMP)
#pragma omp parallel for shared(tod,fit) default(none)
#endif
  // Find each detector's median data value.
  for (int i=0;i<tod->ndet;i++)
    fit->median_vals[i] = compute_median(tod->ndata,fit->data[i]);
  memcpy(fit->data[0],tod->data[0],tod->ndet*tod->ndata*sizeof(actData));
  
  psTrace("moby.pcg",3,"Took %8.5f seconds to find first medians.\n",mbElapsedTime(&ticker));

#if !defined(MB_SKIP_OMP)
#pragma omp parallel for shared(tod,fit) default(none)
#endif

  // Rebias each detector's TOD by removing that detector's median.  Then go to abs value. (?)
  for (int i=0;i<tod->ndet;i++) {
    for (int j=0;j<tod->ndata;j++) {
      fit->data[i][j]-=fit->median_vals[i];
      fit->data[i][j]=fabs(fit->data[i][j]);
    }
  }  
  psTrace("moby.pcg",3,"Made absvecs at %8.5f seconds.\n",mbElapsedTime(&ticker));  
  

#if !defined(MB_SKIP_OMP)
#pragma omp parallel for shared(tod,fit) default(none)
#endif
  // Find median absolute distance between each detector's data and its median.  Call that "scatter".
  for (int i=0;i<tod->ndet;i++)
    fit->median_scats[i] = compute_median(tod->ndata,fit->data[i]);
  

  for (int i=0;i<tod->ndet;i++)
    if (fit->median_scats[i]<=0) {
      psTrace("moby.pcg",4,"Detector index %d had zero scatter.  Should it be cut?\n",i);  
      fit->median_scats[i]=1.0;
    }
  psTrace("moby.pcg",3,"Made median scats at %8.5f seconds.\n",mbElapsedTime(&ticker));  

  // Scratch copy of data was messed up.  Refresh it from the true TOD.
  memcpy(fit->data[0],tod->data[0],tod->ndet*tod->ndata*sizeof(actData));
#if !defined(MB_SKIP_OMP)
#pragma omp parallel for shared(fit) default(none)
#endif

  // Bias the (scratch) data by the detector's median, then rescale by the detector's "scatter".
  for (int i=0;i<fit->ndet;i++)
    for (int j=0;j<fit->ndata;j++)
      fit->data[i][j]= (fit->data[i][j]-fit->median_vals[i])/fit->median_scats[i];
}



/*---------------------------------------------------------------------------------------------------------*/
/// Build list of the un-cut stretches of data associated with a single-detector TOD. If the cuts
/// are NULL, set final istop to -1.
/// \param cuts   The cuts object to check.
/// \param ndata  The number of data values in the TOD linked to the cuts.
/// \param row    The detector row number.
/// \param col    The detector column number.
/// \param istart_out  Returns array of indices for the beginnings of valid ranges.  Return val gives size.
/// \param istop_out   Returns array of indices for the ends of valid ranges.  Return val gives size.
/// \return  Number of valid segments of data.

int mbGetNoCutInds(const mbCuts *cuts, int ndata, int row, int col,
                   int **istart_out, int **istop_out)
{
  int *istart, *istop;
  bool cut_all=true;
  if (cuts)
    if (cuts->detCuts[row][col])
      if (cuts->detCuts[row][col]->ncuts>0)
	cut_all=false;

  // Handle case where cut object is not able to specify real cuts by asserting that no data are good.
  if (cut_all) {
    istart=(int *)psAlloc(sizeof(int));
    istop=(int *)psAlloc(sizeof(int));
    istart[0]=0;
    istop[0]=ndata;
    *istart_out=istart;
    *istop_out=istop;
    return 1;
  } else {

    mbCutList *mycuts=cuts->detCuts[row][col];
    int nseg=mycuts->ncuts+1;
    mbSingleCut *cur=mycuts->head;
    istart=(int *)psAlloc(nseg*sizeof(int));
    istop=(int *)psAlloc(nseg*sizeof(int));

    // Handle boundary case where first value is cut.
    if (cur->indexFirst==0) {
      nseg--;
      istart[0]=cur->indexLast;
      istop[0]=ndata;  //make sure we're ok if we return because cur->next=0;
      cur=cur->next;
    } else
      istart[0]=0;

    int i=0;
    while (cur) {
      istop[i]=cur->indexFirst;
      //fprintf(stderr,"det %2d %2d, segment %d of %d has limits %d %d.\n",row,col,i,nseg,istart[i],istop[i]);
      istart[i+1]=cur->indexLast+1;
      cur=cur->next;
      i++;
    }

    // Handle boundary case where the last value in the set is cut.
    istop[i]=ndata;
    if (istop[i]<=istart[i])
      nseg--;

    // Set up the arrays to be returned.
    *istart_out=istart;
    *istop_out=istop;
    
    return nseg;
  }
}



/*---------------------------------------------------------------------------------------------------------*/
/// Similar to mbRescaleArray, except using means instead of medians and applying cuts.
/// Make a scratch copy of the raw data, but rebias each detector's TOD by subtracting that
/// detector's mean data value, and then rescale so that each detector has unit "scatter", meaning
/// its mean squared deviation from the mean (now zero) is one.
/// \param tod  The TOD to study.
/// \param fit  The common mode object being updated.
/// \param cuts Cuts to apply when computing means.

void mbRescaleArrayMean(const mbTOD *tod, mbNoiseCommonMode *fit, const mbCuts *cuts)
{
  assert(tod!=NULL);
  assert(fit!=NULL);
  // If you fail ndet>1, it's because you only sent in one detector, which isn't much of a common mode.
  assert(tod->ndet>1);  
  assert(tod->ndata>1);
  assert(tod->ndata==fit->ndata);
  assert(tod->ndet==fit->ndet);

  mbTimeValue ticker;
  mbStartTime(&ticker);
  
  bool had_data=false;
  if (fit->have_data)
    had_data=true;
  else {
    fit->data=psAllocMatrix(fit->ndet,fit->ndata);
    fit->have_data=true;
  }
  
  // Make a scratch copy of the full TOD 
  memcpy(fit->data[0],tod->data[0],tod->ndet*tod->ndata*sizeof(actData));
  psTrace("moby.pcg",3,"Allocated space.\n");

#pragma omp parallel for shared(tod,fit,cuts) default(none)

  // Compute the mean of all not-cut data.
  for (int i=0;i<tod->ndet;i++) {
    if (mbCutsIsAlwaysCut(cuts, tod->rows[i],tod->cols[i])) {
      fit->median_vals[i]=0;
      continue;
    }
    actData sum=0;
    int nseg, *istart, *istop;
    nseg=mbGetNoCutInds(cuts,tod->ndata, tod->rows[i],tod->cols[i],&istart,&istop);
    int ndata_good=0;
    for (int j=0;j<nseg;j++) {
      for (int k=istart[j];k<istop[j];k++)
        sum+=tod->data[j][k];
      ndata_good+=(istop[j]-istart[j]);
    }
      
    psFree(istart);
    psFree(istop);
    fit->median_vals[i]=sum/(actData)ndata_good;
  }

  memcpy(fit->data[0],tod->data[0],tod->ndet*tod->ndata*sizeof(actData));  
  psTrace("moby.pcg",3,"Took %8.5f seconds to find first medians.\n",mbElapsedTime(&ticker));

#pragma omp parallel for shared(tod,fit) default(none)
  // Rebias each detector's TOD by subtracting that detector's mean value.
  for (int i=0;i<tod->ndet;i++) {
    for (int j=0;j<tod->ndata;j++)	{
      fit->data[i][j]-=fit->median_vals[i];
      // fit->data[i][j]=fabs(fit->data[i][j]); Don't need absolute value vectors since we're going
      // to be doing an RMS calculation instead of a median one.
    }
  }  
  psTrace("moby.pcg",3,"Made subtracted vecs at %8.5f seconds.\n",mbElapsedTime(&ticker));  
  

#pragma omp parallel for shared(tod,fit,cuts) default(none)
  // Find RMS deviation about the mean for each detector.
  for (int i=0;i<tod->ndet;i++) {
    psTrace("moby.pcg",6,"working on %d, which is %2d %2d\n",i,tod->rows[i],tod->cols[i]);
    if (mbCutsIsAlwaysCut(cuts, tod->rows[i],tod->cols[i])) {
      fit->median_scats[i]=1.0;
      continue;
    }

    actData sumsqr=0;
    int nseg, *istart, *istop;
    nseg=mbGetNoCutInds(cuts,tod->ndata,tod->rows[i],tod->cols[i],&istart,&istop);
    psTrace("moby.pcg",6,"nseg is %d, lims are %d %d\n",nseg,istart[0],istop[nseg-1]);
    int ndata_good=0;
    for (int j=0;j<nseg;j++) {
      for (int k=istart[j];k<istop[j];k++)
        sumsqr+=fit->data[i][k]*fit->data[i][k];
      ndata_good+=(istop[j]-istart[j]);
    }
    
    psFree(istart);
    psFree(istop);
    fit->median_scats[i]=sqrt(sumsqr/ndata_good);
  }
  
  
  for (int i=0;i<tod->ndet;i++)
    if (fit->median_scats[i]<=0) {
      psTrace("moby.pcg",1,"Detector index %d had zero scatter.  Should it be cut?\n",i);  
      fit->median_scats[i]=1.0;
    }
  psTrace("moby.pcg",3,"Made median scats at %8.5f seconds.\n",mbElapsedTime(&ticker));  
  
  

#pragma omp parallel for shared(fit) default(none)
  // Scratch copy of data was NOT messed up.  It's had the mean subtracted; now rescale to give it
  // unit RMS.
  for (int i=0;i<fit->ndet;i++)
    for (int j=0;j<fit->ndata;j++)
      fit->data[i][j] /=fit->median_scats[i];
}



/*---------------------------------------------------------------------------------------------------------*/
/// Compute the "common mode" vector as the median over all detectors of the data at each time step.
/// Note that here "the data" means data with a bias and rescaling to be zero-median, unit-"scatter".
/// \param fit  The common mode object being used and updated.

void mbCalculateCommonMode(mbNoiseCommonMode *fit)
{
  assert(fit->have_data);  /*check to make sure you pre-calculate the rescaled data.*/
  
#pragma omp parallel shared(fit) default(none)
  {
    actData *mymedian=(actData *)psAlloc(fit->ndet*sizeof(actData));
    for (int i=0;i<fit->ndet;i++)
      mymedian[i]=0;
    
#pragma omp for
    for (int i=0;i<fit->ndata;i++)  {
      for (int j=0;j<fit->ndet;j++)
        mymedian[j]=fit->data[j][i];
      
      fit->common_mode[i] = compute_median(fit->ndet,mymedian);
    }
    psFree(mymedian);
  }
}



/*---------------------------------------------------------------------------------------------------------*/
/// Compute the "common mode" vector as the mean over all not-cut detectors of the data at each time
/// step.  Note that here "the data" means data with a bias and rescaling to be zero-mean,
/// unit-RMS.
/// \param fit  The common mode object being used and updated.
/// \param tod  The TOD being averaged.
/// \param cuts Cuts to apply to the TOD.

void mbCalculateMeanCommonMode(mbNoiseCommonMode *fit, const mbTOD *tod, const mbCuts *cuts)
{
  assert(fit->have_data);  /*check to make sure you pre-calculate the rescaled data.*/
  
  if (fit->have_weights==false) { //if we don't have weights, put in ones.
    fit->weights=(actData *)psAlloc(sizeof(actData)*fit->ndet);
    for (int i=0;i<fit->ndet;i++)
      fit->weights[i]=1.0;
    fit->have_weights=true;
  }
  
  actData *weightsum=(actData *)psAlloc(sizeof(actData)*fit->ndata);
  
  for (int i =0;i<fit->ndata;i++) {
    weightsum[i]=0;
    fit->common_mode[i]=0;
  }

  for (int j=0;j<fit->ndet;j++) {
    if (cuts)
      if (mbCutsIsAlwaysCut(cuts, tod->rows[j],tod->cols[j]))
        continue;
    
    bool use_all_of_me=true;
    if (cuts) {
      mbCutList *mycuts=cuts->detCuts[tod->rows[j]][tod->cols[j]];
      if (mycuts)
        if (mycuts->ncuts>0)
          use_all_of_me=false;
    }
    if (use_all_of_me) {  // Why not take all of me?
      for (int i=0;i<fit->ndata;i++) {
        weightsum[i]+=fit->weights[j];
        fit->common_mode[i]+=fit->data[j][i];
      }

    } else {
      //fprintf(stderr,"Dealing with cuts on detector %d %d\n",tod->rows[j],tod->cols[j]);
      mbCutList *mycuts=cuts->detCuts[tod->rows[j]][tod->cols[j]];
      int istart=0;
      mbSingleCut *cur=mycuts->head;
      while (cur)	{
        for (int i=istart; i<cur->indexFirst; i++) {
          weightsum[i]+=fit->weights[j];
          fit->common_mode[i]+=fit->data[j][i];		  
        }
        istart=cur->indexLast+1;
        cur=cur->next;
      }
      for (int i=istart;i<fit->ndata;i++) {
        weightsum[i]+=fit->weights[j];
        fit->common_mode[i]+=fit->data[j][i];		  
      }
	      
    }
  }
  for (int i=0;i<fit->ndata;i++)
    fit->common_mode[i]/=weightsum[i];
  
  psFree(weightsum);
}



/*---------------------------------------------------------------------------------------------------------*/
/// Apply an anti-glitch filter to remove calbols from the common-mode data set.  Sets
/// mbNoiseCommonMode *fit->icutvec such that times with a calbol are marked with a 1 and clears out
/// evidence of calbol from fit->commom_mode_nocalbol -- uses glitch_one_detector_simple
/// \param fit  The object to update.

void mbFindCalbols(mbNoiseCommonMode *fit)
{
  bool do_cuts=true;
  bool do_smooth=false;
  bool apply_glitch=true;
  
  
  //this is something that I don't think needs to be here, but I'm getting segfaults if I don't
  memcpy(fit->common_mode_nocalbol,fit->common_mode,sizeof(actData)*fit->ndata);

  fit->icutvec=glitch_one_detector_simple(fit->common_mode_nocalbol, 
                                          fit->ndata,
                                          do_smooth,
                                          apply_glitch,
                                          do_cuts,
                                          fit->nsig,
                                          fit->tGlitch,
                                          fit->tSmooth,
                                          fit->dt,
                                          1);

#if 0
  int *icutvec=fit->icutvec;
  int ncalbol=0;
  if (do_cuts) {
      
    /*first, pad all calbol events with a pad's worth of ones, so that almost touching
      events are counted as one.*/
    int npad=fit->calbol_pad/fit->dt;
    int *icutvec_copy=(int *)psAlloc(fit->ndata*sizeof(int));      
    memcpy(icutvec_copy,icutvec,fit->ndata*sizeof(int));
    int nflip=0;
    //fprintf(stderr,"icutvec[0]=%d\n",icutvec[0]);
    int istart=-1;
    int istop=-1;
    for (int i=1;i<fit->ndata;i++)
      if (icutvec[i]!=icutvec[i-1])  {
        if (icutvec[i]) {
          istart=i;
          ncalbol++;
          int jmin=i-npad;
          if (jmin<0)
            jmin=0;
          for (int j=jmin;j<=i;j++)
            icutvec_copy[j]=1;
        } else {
          istop=i;
          int jmax=i+npad;
          if (jmax>fit->ndata)
            jmax=fit->ndata;
          for (int j=i;j<jmax;j++)
            icutvec_copy[j]=1;
        }
        fprintf(stderr,"icutvec[%d]=%d\n",i,icutvec[i]);
        nflip++;
      }
    memcpy(icutvec,icutvec_copy,fit->ndata*sizeof(int));
      

    /*OK - now check to see how many calbols we really have*/
    nflip=0;
    fprintf(stderr,"icutvec[0]=%d\n",icutvec[0]);
    istart=-1;
    istop=-1;
    ncalbol=0;
    for (int i=1;i<fit->ndata;i++)
      if (icutvec[i]!=icutvec[i-1]) {
        if (icutvec[i]) {
          istart=i;
          ncalbol++;
        } else {
          istop=i;
        }
        fprintf(stderr,"icutvec[%d]=%d\n",i,icutvec[i]);
        nflip++;
      }

      

    if ((nflip==2)&&(icutvec[0]==0))	{
      //Don't need these pads, since we already did 'em
      //istart-=npad;
      //istop+=npad;
      actData left_mean=array_mean(&(fit->common_mode_nocalbol[istart-npad]),npad);
      actData right_mean=array_mean(&(fit->common_mode_nocalbol[istop]),npad);
      actData i1=istart-npad/2;
      actData i2=istop+npad/2;
      for (int i=istart;i<istop;i++) {
        fit->common_mode_nocalbol[i]=left_mean+(right_mean-left_mean)*(i-i1)/(i2-i1);
      }
      //fprintf(stderr,"have 1 cut interval between %8.3f and %8.3f\n",tod->dt[istart],tod->dt[istop]);
    }
    psFree(icutvec);
    psFree(icutvec_copy);

      

    if (ncalbol>1){
      fprintf(stderr,"WARNING WARNING WARNING - chopping off some calbols, currently have %d.\n",ncalbol);
      ncalbol=1;

    }

  }

  //fprintf(stderr,"finished with %d calbols.\n",ncalbol);
  fit->ncalbol=ncalbol;
  fit->nparam=fit->np_common+fit->np_poly+fit->ncalbol;
  if (fit->fit_params) {
    psFree(fit->fit_params[0]);
    psFree(fit->fit_params);
    fit->fit_params=psAllocMatrix(fit->nparam,fit->ndet);
  }



  actData *ivec=(actData *)psAlloc(fit->ndata*sizeof(actData));
  for (int i=0;i<fit->ndata;i++)
    ivec[i]=((actData)i)/((actData)fit->ndata);

  actData **vecs=psAllocMatrix(fit->nparam,fit->ndata);
#endif
}



/*---------------------------------------------------------------------------------------------------------*/
/// Allocate space in the common mode object to store data about some number of calbol pulses.
/// \param fit  The object to update.
/// \param ncalbol   The number of calbol pulses found.

void mbAllocateCalbols(mbNoiseCommonMode *fit,int ncalbol)
{
  fit->ncalbol=ncalbol;

  if (ncalbol>0) {
    fit->calbol_start=(int *)psAlloc(sizeof(int)*ncalbol);
    fit->calbol_stop=(int *)psAlloc(sizeof(int)*ncalbol);
    fit->nparam=fit->np_poly+fit->np_common+fit->ncalbol;
  
    psFree(fit->fit_params[0]);
    psFree(fit->fit_params);
    fit->fit_params=psAllocMatrix(fit->nparam,fit->ndet);
  }
}



/*---------------------------------------------------------------------------------------------------------*/
/// Pre-compute vectors needed for common-mode fitting.  These include the fractional time through
/// the TOD and various powers of it, and also these same vectors times the common_mode_nocalbol.
/// \param fit  The object to update.

void mbCalculateCommonModeVecs(mbNoiseCommonMode *fit)
{
  // A vector to hold the fractional time through the TOD for each data index.
  actData *ivec=(actData *)psAlloc(fit->ndata*sizeof(actData));
  for (int i=0; i<fit->ndata; i++)
    ivec[i]=((actData)i)/((actData)fit->ndata);

  // A matrix to hold the pre-computed fractional time to various powers, up to power (np_poly-1).
  actData **vecs=psAllocMatrix(fit->nparam,fit->ndata);
  
  if (fit->np_poly>0)
#pragma omp parallel for shared(fit,ivec,vecs) default(none)
    for (int i=0; i<fit->ndata; i++) {
      vecs[0][i]=1;
      for (int j=1;j<fit->np_poly;j++)
        vecs[j][i]=vecs[j-1][i]*ivec[i];
    }

  // Matrix vecs also holds precomputed common mode function times various powers of the fractional time.
  if (fit->np_common>0)
#pragma omp parallel for shared(fit,ivec,vecs) default(none)
    for (int i=0; i<fit->ndata; i++) {
      vecs[0+fit->np_poly][i] = fit->common_mode_nocalbol[i];
      for (int j=1;j<fit->np_common;j++)
        vecs[j+fit->np_poly][i]=vecs[j+fit->np_poly-1][i]*ivec[i];
    }    
  
  // Loop over all calbol pulses.
  for (int j=0;j<fit->ncalbol;j++) {
    int ind=fit->np_poly+fit->np_common+j;
/*     fprintf(stderr,"Putting in index %d out of %d\n",ind,fit->nparam); */
    memset(vecs[ind],0,fit->ndata*sizeof(actData));

    // Compute a range [imin, imax] where the calbol is to be cut.
    int npad=fit->calbol_pad/fit->dt;
    int imin=fit->calbol_start[j]-npad;
    if (imin<0)
      imin=0;
    int imax=fit->calbol_stop[j]+npad;
    if (imax>fit->ndata)
      imax=fit->ndata;

    for (int i=imin;i<imax;i++)
      vecs[ind][i]=fit->common_mode[i]-fit->common_mode_nocalbol[i];
  }
  fit->vecs=vecs;
  fit->have_vecs=true;
  psFree(ivec);
}



/*---------------------------------------------------------------------------------------------------------*/
/// Compute the common mode fit parameters for each detector.  There will be fit->nparam per det,
/// which equals the total of the number of polynomial orders of the two fits and the number of
/// calbol pulses removed.
/// \param fit  The object to update.

void mbCalculateCommonModeFitParams(mbNoiseCommonMode *fit)
{
  assert(fit->have_data==true);  //A copy of the rescaled data should exist here.

  actData **ata=psAllocMatrix(fit->nparam,fit->nparam);


#if 1
  act_gemm('N','T',fit->nparam,fit->nparam,fit->ndata,1,fit->vecs[0],fit->ndata,fit->vecs[0],fit->ndata,0,ata[0],fit->nparam);
#else
#ifdef ACTDATA_DOUBLE
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans, fit->nparam, fit->nparam, fit->ndata,1, 
              fit->vecs[0], fit->ndata ,fit->vecs[0],fit->ndata,0,ata[0],fit->nparam);
#else
  cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasTrans, fit->nparam, fit->nparam, fit->ndata,1, 
              fit->vecs[0], fit->ndata ,fit->vecs[0],fit->ndata,0,ata[0],fit->nparam);
#endif
#endif

  //psTrace("moby.pcg",3,"Made ata at %8.5f seconds.\n",mbElapsedTime(&ticker));      
  
  if (!fit->have_ata) { //if we don't have te 
    fit->ata=psAllocMatrix(fit->nparam,fit->nparam);
    fit->have_ata=true;
    
    memcpy(fit->ata[0],ata[0],sizeof(actData)*fit->nparam*fit->nparam);
  } else
    psTrace("moby.pcg",0,"Already think I have ata in mbCalculateCommonModeFitParams.  Be careful...\n");

  // Compute Cholesky factorization of a real symmetric positive definite matrix.
  // ata[0] goes in as the real symmetric matrix to be factored and returns as the upper-triangular
  // "Cholesky square root" U, where the input matrix A = (U^T) U.
  int info=0;
#ifdef ACTDATA_DOUBLE
  clapack_dpotrf('u',fit->nparam,ata[0],fit->nparam,&info);
#else
  clapack_spotrf('u',fit->nparam,ata[0],fit->nparam,&info);
#endif
  if(info!=0) {
    psTrace("moby",0,"Failed spotrf\n");
    return;
  }

  // Compute inverse of a matrix already Cholesky factorized.
#ifdef ACTDATA_DOUBLE
  clapack_dpotri('u',fit->nparam,ata[0],fit->nparam,&info);
#else
  clapack_spotri('u',fit->nparam,ata[0],fit->nparam,&info);
#endif
  if(info!=0) {
    psTrace("moby",0,"Failed spotri\n");
    return;
  }

  // _spotri stores the result upper-diagonal.  Fill the opposite side.
  for (int i=0;i<fit->nparam;i++)
    for (int j=i+1;j<fit->nparam;j++)
      ata[i][j]=ata[j][i];

  // ???
  actData **tmp_mat=psAllocMatrix(fit->nparam,fit->ndet);

#if 1
  act_gemm('N','T',fit->nparam,fit->ndet,fit->ndata,1.0,fit->vecs[0],fit->ndata,fit->data[0],fit->ndata,0.0,tmp_mat[0],fit->ndet);
#else
#ifdef ACTDATA_DOUBLE
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans, fit->nparam, fit->ndet, fit->ndata,1, 
              fit->vecs[0], fit->ndata ,fit->data[0],fit->ndata,0,tmp_mat[0],fit->ndet);
#else
  cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasTrans, fit->nparam, fit->ndet, fit->ndata,1, 
              fit->vecs[0], fit->ndata ,fit->data[0],fit->ndata,0,tmp_mat[0],fit->ndet);
#endif
#endif
  if ((!fit->have_atx)||(fit->atx==NULL)) {
    fit->atx=psAllocMatrix(fit->nparam,fit->ndet);
    memcpy(fit->atx[0],tmp_mat[0],fit->nparam*fit->ndet*sizeof(actData));
    fit->have_atx=true;
  } else {
    fprintf(stderr,"screwed up...fit->have_atx already exists\n");
  }

  //psTrace("moby.pcg",3,"Made temp mat at %8.5f seconds.\n",mbElapsedTime(&ticker)); 

  for (int i=0;i<fit->nparam;i++)
    for (int j=0;j<fit->ndet;j++) {
      fit->fit_params[i][j]=0;
      for (int k=0;k<fit->nparam;k++)
        fit->fit_params[i][j]+=ata[i][k]*tmp_mat[k][j];
    }
  
  sgemm_simple( fit->nparam,fit->ndet,fit->nparam,ata[0],0 ,tmp_mat[0],0,fit->fit_params[0]);

  // Now apply the median scatter scaling to the fit parameters
  for (int i=0;i<fit->ndet;i++)
    for (int j=0;j<fit->nparam;j++)
      fit->fit_params[j][i]*=fit->median_scats[i];

  psFree(tmp_mat[0]);
  psFree(tmp_mat);
  psFree(ata[0]);
  psFree(ata);
}



/*---------------------------------------------------------------------------------------------------------*/
///  Remove the common mode, plus extra polynomial terms from the tod.  Optionally remove calbols.
///  Assumes that the fit parameters have already been computed.
/// \param tod   The TOD to adjust.
/// \param fit   The common mode data to be applied.
/// \param cutCalbols  Whether to cut calbol pulse periods.

void mbApplyCommonMode(mbTOD *tod, mbNoiseCommonMode *fit,bool cutCalbols)
{
  // data_fit will hold the full TOD-sized array of best fits to the common mode for each detector.
  actData **data_fit=psAllocMatrix(tod->ndet,tod->ndata);  
  if (cutCalbols)
    psTrace("moby.pcg",3,"Cutting calbols.\n");
  else
    psTrace("moby.pcg",3,"Keeping calbols.\n");

  if (cutCalbols)
    sgemm_simple( fit->ndet,fit->ndata,fit->nparam,fit->fit_params[0],0 ,fit->vecs[0],0,data_fit[0]);
  else
    sgemm_simple( fit->ndet,fit->ndata,fit->np_poly+fit->np_common,fit->fit_params[0],0 ,fit->vecs[0],0,data_fit[0]);
  //psTrace("moby.pcg",3,"Made data templates at %8.5f seconds.\n",mbElapsedTime(&ticker));        
  psTrace("moby.pcg",3,"Made data templates.\n");

#pragma omp parallel for shared(tod,data_fit,fit) default(none)
  for (int i=0;i<tod->ndet;i++)
    for (int j=0;j<tod->ndata;j++)
      tod->data[i][j] = tod->data[i][j]-(data_fit[i][j]+fit->median_vals[i]);
  //psTrace("moby.pcg",3,"Replaced with common mode at %8.5f seconds.\n",mbElapsedTime(&ticker)); 
  psTrace("moby.pcg",3,"Replaced with common mode.\n"); 
  fit->common_is_applied=true;

  //Was turned on, don't remember what it's doing...
  //cblas_saxpy(tod->ndata*tod->ndet,-1.0,data_fit[0],1,tod->data[0],1);
  //for (int i=0;i<tod->ndet;i++)
  //  fit->frac_errs->data.F32[i]=sqrt(cblas_sdot(tod->ndata,data[i],1,data[i],1)/((actData)tod->ndata));
  

  //#if !defined(MB_SKIP_OMP)
  //#pragma omp parallel for shared(tod,data_fit,fit) default(none)
  //#endif
  //for (int i=0;i<tod->ndet;i++)
   // {
  //  actData *vec=(actData *)psAloc(sizeof(actData)*fit->ndata);
  //  for (int j=0;j<fit->ndata;j++)
  //vec[j]=tod->data[i][j]-fit->median_vals[i];
  //  fit->frac_errs[i]=sqrt(cblas_sdot(tod->ndata,data[i],1,data[i],1)/((actData)tod->ndata));
  //  psFree(vec);
  // }

  psFree(data_fit[0]);
  psFree(data_fit);
}



/*---------------------------------------------------------------------------------------------------------*/
/// Compute the fractional error for each detector's common-mode fit.
/// \param tod   The TOD to check.
/// \param fit   The common mode being fit against.

void mbCalculateCommonFracErrs(const mbTOD *tod, mbNoiseCommonMode *fit)
{
#pragma omp parallel for shared(tod,fit) default(none)
  for (int i=0;i<tod->ndet;i++) {
    // If the common mode has been applied, we can go to town.  Otherwise, we need to calculate the
    // de-common-moded data before we can do the fractional errors.
    if (fit->common_is_applied) {

#if 1
      fit->frac_errs[i]=sqrt(act_dot(tod->ndata,tod->data[i],1,tod->data[i],1)/((actData)tod->ndata))/fit->median_scats[i];
#else
#ifdef ACTDATA_DOUBLE
      fit->frac_errs[i]=sqrt(cblas_ddot(tod->ndata,tod->data[i],1,tod->data[i],1)/
                             ((actData)tod->ndata))/fit->median_scats[i];
#else
      fit->frac_errs[i]=sqrt(cblas_sdot(tod->ndata,tod->data[i],1,tod->data[i],1)/
                             ((actData)tod->ndata))/fit->median_scats[i];
#endif
#endif
      continue;
    } else {

      actData *vec=(actData *)psAlloc(sizeof(actData)*fit->ndata);
      actData *params=(actData *)psAlloc(sizeof(actData)*fit->nparam);
      for (int j=0;j<fit->nparam;j++)
        params[j]=fit->fit_params[j][i];
#ifdef ACTDATA_DOUBLE
      cdgemv('n',fit->ndata,fit->nparam,1.0,fit->vecs[0],fit->ndata,params,1,0.0,vec,1);
#else
      csgemv('n',fit->ndata,fit->nparam,1.0,fit->vecs[0],fit->ndata,params,1,0.0,vec,1);
#endif
      for (int j=0;j<fit->ndata;j++)
        vec[j]=(tod->data[i][j]-vec[j]-fit->median_vals[i]);


#if 1
      fit->frac_errs[i]=sqrt(act_dot(tod->ndata,vec,1,vec,1)/((actData)tod->ndata))/fit->median_scats[i];
#else
#ifdef ACTDATA_DOUBLE
      fit->frac_errs[i]=sqrt(cblas_ddot(tod->ndata,vec,1,vec,1)/((actData)tod->ndata))/fit->median_scats[i];
#else
      fit->frac_errs[i]=sqrt(cblas_sdot(tod->ndata,vec,1,vec,1)/((actData)tod->ndata))/fit->median_scats[i];
#endif      
#endif
      
      psFree(params);
      psFree(vec);
    }
  }
}



/*---------------------------------------------------------------------------------------------------------*/
/// Read/Write the common mode to disk.  
/// \param fit_in  I/O access to pointer to the common mode object.
/// \param fname   Filename to read or write.
/// \param doWrite If true, write to a file; if false, read from a file.

static psErrorCode mbReadWriteCommon(mbNoiseCommonMode **fit_in, const char *fname, bool doWrite)
{
  FILE *iofile;
  mbNoiseCommonMode *fit;
  if (doWrite) {
    fit=*fit_in;
    iofile=fopen(fname,"w");
    if (!iofile) {
      psTrace("moby.pcg",1,"Requested file %s file not available for writing in mbReadWriteCommon.\n",fname);
      return PS_ERR_BAD_PARAMETER_VALUE;
    }
  } else {
    fit=(mbNoiseCommonMode *)psAlloc(sizeof(mbNoiseCommonMode));
    *fit_in=fit;
    iofile=fopen(fname,"r");
    if (!iofile) {
      psTrace("moby.pcg",1,"Requested file %s file not available for reading in mbReadWriteCommon.\n",fname);
      return PS_ERR_BAD_PARAMETER_VALUE;
    }      
  }  
  
  bool doRead;
  if (doWrite)
    doRead=true;
  else
    doRead=false;
  mbFreadWrite(doWrite,&fit->ndata,sizeof(int),1,iofile);
  mbFreadWrite(doWrite,&fit->ndet,sizeof(int),1,iofile);  
  mbFreadWrite(doWrite,&fit->np_poly,sizeof(int),1,iofile);
  mbFreadWrite(doWrite,&fit->np_common,sizeof(int),1,iofile);
  mbFreadWrite(doWrite,&fit->ncalbol,sizeof(int),1,iofile);
  mbFreadWrite(doWrite,&fit->nparam,sizeof(int),1,iofile);
  if (fit->ncalbol) {
    if (doRead) {
      fit->calbol_start=(int *)psAlloc(sizeof(int)*fit->ncalbol);
      fit->calbol_stop=(int *)psAlloc(sizeof(int)*fit->ncalbol);
    }
    mbFreadWrite(doWrite,fit->calbol_start,sizeof(int),fit->ncalbol,iofile);
    mbFreadWrite(doWrite,fit->calbol_stop,sizeof(int),fit->ncalbol,iofile);
  }

  if (doRead) {
    fit->common_mode=(actData *)psAlloc(sizeof(actData)*fit->ndata);
    fit->common_mode_nocalbol=(actData *)psAlloc(sizeof(actData)*fit->ndata);
    fit->median_vals=(actData *)psAlloc(sizeof(actData)*fit->ndet);
    fit->median_scats=(actData *)psAlloc(sizeof(actData)*fit->ndet);
    fit->frac_errs=(actData *)psAlloc(sizeof(actData)*fit->ndet);
    fit->fit_params=psAllocMatrix(fit->nparam,fit->ndet);
  }

  mbFreadWrite(doWrite,fit->common_mode,sizeof(actData),fit->ndata,iofile);
  mbFreadWrite(doWrite,fit->common_mode_nocalbol,sizeof(actData),fit->ndata,iofile); 
  mbFreadWrite(doWrite,fit->median_vals,sizeof(actData),fit->ndet,iofile);
  mbFreadWrite(doWrite,fit->median_scats,sizeof(actData),fit->ndet,iofile);
  mbFreadWrite(doWrite,fit->frac_errs,sizeof(actData),fit->ndet,iofile);
  mbFreadWrite(doWrite,fit->fit_params[0],sizeof(actData),fit->ndet*fit->nparam,iofile);
  mbFreadWrite(doWrite,&fit->calbol_pad,sizeof(actData),1,iofile);
  mbFreadWrite(doWrite,&fit->have_ata,sizeof(bool),1,iofile);
  if (fit->have_ata) {
    if (doRead) 
      fit->ata=psAllocMatrix(fit->nparam,fit->nparam);
    
    mbFreadWrite(doWrite,fit->ata[0],sizeof(actData),fit->nparam*fit->nparam,iofile);
  }
  mbFreadWrite(doWrite,&fit->have_atx,sizeof(bool),1,iofile);
  if (fit->have_atx) {
    if (doRead)
      fit->atx=psAllocMatrix(fit->nparam,fit->ndet);
    mbFreadWrite(doWrite,fit->atx[0],sizeof(actData),fit->ndet*fit->nparam,iofile);
  }
  mbFreadWrite(doWrite,&fit->nsig,sizeof(actData),1,iofile);
  mbFreadWrite(doWrite,&fit->tGlitch,sizeof(actData),1,iofile);
  mbFreadWrite(doWrite,&fit->tSmooth,sizeof(actData),1,iofile);
  mbFreadWrite(doWrite,&fit->dt,sizeof(actData),1,iofile);

  mbFreadWrite(doWrite,&fit->have_unsmoothed_ratio,sizeof(bool),1,iofile);
  if (doRead)  //bass-ackwards, since even though I have a flag, the alloc routine makes space for this.
    fit->unsmooth_ratio=(actData *)psAlloc(fit->ndet*sizeof(actData));
  mbFreadWrite(doWrite,fit->unsmooth_ratio,sizeof(actData),fit->ndet,iofile);
  
  
  mbFreadWrite(doWrite,&fit->apply_common,sizeof(bool),1,iofile);
  mbFreadWrite(doWrite,&fit->common_is_applied,sizeof(bool),1,iofile);
  mbFreadWrite(doWrite,&fit->keep_vecs,sizeof(bool),1,iofile);
  
  mbFreadWrite(doWrite,&fit->have_vecs,sizeof(bool),1,iofile);
  if (fit->have_vecs) {
    if (doRead) {
      fit->vecs=psAllocMatrix(fit->nparam,fit->ndata);
      // Make sure we zero since we won't read the full calbol in:
      memset(fit->vecs[0],0,sizeof(actData)*fit->nparam*fit->ndata); 
    }

    //Do a slight modification here, only save the calbols where they're non-zero
    int nfitparam=fit->np_poly+fit->np_common;
    mbFreadWrite(doWrite,fit->vecs[0],sizeof(actData),nfitparam*fit->ndata,iofile);
    for (int i=0;i<fit->ncalbol;i++)
      mbFreadWrite(doWrite,&fit->vecs[nfitparam+i][fit->calbol_start[i]],sizeof(actData),
                 fit->calbol_stop[i]-fit->calbol_start[i]+1,iofile);
  }
  
  fclose(iofile);
  return PS_ERR_NONE;
}



/// Read in a common mode file.
/// Given the filename of a common mode file (see mbWriteCommonMode), allocates and
/// Uses mbReadWriteCommon.
/// \param filename  The file to open and read.
/// \return  A pointer to the corresponding mbNoiseCommonMode structure.

mbNoiseCommonMode *mbReadCommon( const char *filename )
{
  mbNoiseCommonMode *n;
  if ( mbReadWriteCommon( &n, filename, false ) != PS_ERR_NONE ){
#ifdef HAVE_PSERR
    psError( MB_ERR_IO, 1, "mbReadCommon: Trouble reading in file %s.", filename );
#else
    fprintf(stderr,"mbReadCommon: Trouble reading in file %s.", filename );
#endif
  }
  return n;
}


/// Write a common mode file from a structure.  Given the filename of a common mode file (see
/// mbWriteCommonMode) and an mbNoiseCommonMode structure, writes the structure to a binary file.
/// Uses mbReadWriteCommon.
/// \param filename The file to open and read.

void mbWriteCommon( const char *filename, mbNoiseCommonMode *common )
/// write an mbnoiseCommonMode to a file
{
  mbNoiseCommonMode **n;
  n = &common;
  if ( mbReadWriteCommon( n, filename, true ) != PS_ERR_NONE ){
#ifdef HAVE_PSERR
    psError( MB_ERR_IO, 1, "mbWriteCommon: Trouble writing to file %s.", filename );
#else
    fprintf(stderr,"mbWriteCommon: Trouble writing to file %s.", filename );
#endif
  }
}


 
/*---------------------------------------------------------------------------------------------------------*/
/// If Common mode fits have already been computed, change the results based on applying a set of
/// cuts.
/// \param tod  TOD being studied.
/// \param fit  Common mode to be updated.
/// \param cuts Cuts to apply.

void mbApplyCutsToCommonMode(const mbTOD *tod, mbNoiseCommonMode *fit, const mbCuts *cuts)
{
  if (cuts==NULL) {
    psTrace("moby.pcg",1,"Skipping mbAplyCutsToCommonMode as cuts is NULL.\n");
    return;
  }

  assert(tod->ndet==fit->ndet);  //and more checks
#if 0
  fprintf(stderr,"nrow and ncol are %d %d\n",cuts->nrow,cuts->ncol);
  if (cuts->detCuts[15][15]) {
    fprintf(stderr,"15,15 has cuts.\n");
    mbCutList *mycuts=cuts->detCuts[15][15];
    fprintf(stderr,"ncuts is %d\n",mycuts->ncuts);
    mbSingleCut *cur=mycuts->head;
    for (int i=0;i<mycuts->ncuts;i++) {
      assert(cur!=NULL);
      fprintf(stderr,"lims are %6d %6d\n",cur->indexFirst,cur->indexLast);
      cur=cur->next;
    }
  } else
    fprintf(stderr,"missing cuts on 15,15.\n");
  for (int i=0;i<tod->ndet;i++) {
    //fprintf(stderr,"row and col are %3d %3d\n",tod->rows[i],tod->cols[i]);
  }
#endif

  actData **myata=psAllocMatrix(fit->nparam,fit->nparam);
  actData *myatx=(actData *)psAlloc(fit->nparam*sizeof(actData));
  for (int i=0;i<tod->ndet;i++) {
    mbCutList *mycuts=cuts->detCuts[tod->rows[i]][tod->cols[i]];
    if (mycuts)  {//we have some cuts to do.
      memcpy(myata[0],fit->ata[0],fit->nparam*fit->nparam*sizeof(actData));
      for (int j=0;j<fit->nparam;j++)
        myatx[j]=fit->atx[j][i];
      mbSingleCut *cur=mycuts->head;
      while (cur) {
        int jmin=cur->indexFirst;
        int jmax=cur->indexLast+1;
        //fprintf(stderr,"backing off %d %d\n",jmin,jmax);
#if 1
        for (int j=jmin; j<jmax;j++)
          for (int k=0;k<fit->nparam;k++)
            myatx[k]-=fit->vecs[k][j]*(tod->data[i][j]-fit->median_vals[i])/fit->median_scats[i];
#else
        printf("hello!\n");
#endif
	      
#if 1
        for (int m=0;m<fit->nparam;m++)
          for (int n=0;n<fit->nparam;n++)
            for (int j=jmin;j<jmax;j++)
              myata[m][n]-=fit->vecs[m][j]*fit->vecs[n][j];
#else
        printf("hello!\n");
#endif
        cur=cur->next;
      }
      mbInvertPosdefMat(myata,fit->nparam);
#ifdef ACTDATA_DOUBLE
      cdgemv('n',fit->nparam,fit->nparam,fit->median_scats[i],myata[0],fit->nparam,
             myatx,1,0.0,&fit->fit_params[0][i],fit->ndet);
#else
      csgemv('n',fit->nparam,fit->nparam,fit->median_scats[i],myata[0],fit->nparam,
             myatx,1,0.0,&fit->fit_params[0][i],fit->ndet);
#endif
      //for (int j=0;j<fit->nparam;j++)
      //  fprintf(stderr,"second: %d = %14.5e\n",j,fit->fit_params[j][i]);
      
    }
  }
  psFree(myatx);
  psFree(myata[0]);
  psFree(myata);
}



/*---------------------------------------------------------------------------------------------------------*/
/// Compute the ratio of each detector's data RMS with smoothing to without.
/// If fit object has no smoothing time, use 2.0 seconds.
/// \param tod  The data set to smooth.
/// \param fit  The common mode object.

void mbCalculateUnsmoothRatio(mbTOD *tod, mbNoiseCommonMode *fit)
{
  assert(tod->ndata==fit->ndata);
  assert(tod->ndet==fit->ndet);

  actData **dataCopy=psAllocMatrix(tod->ndet,tod->ndata); 
  memcpy(dataCopy[0],tod->data[0],tod->ndet*tod->ndata*sizeof(actData));

  if (fit->t_unsmooth<=0)
    fit->t_unsmooth=2.0;
  smooth_tod(tod,fit->t_unsmooth,dataCopy);
  
#pragma omp parallel for shared(tod,fit,dataCopy) default (none)
  for (int i=0;i<tod->ndet;i++){


#if 1
    fit->unsmooth_ratio[i]=sqrt(act_dot(tod->ndata,dataCopy[i],1,dataCopy[i],1)/act_dot(tod->ndata,tod->data[i],1,tod->data[i],1));
#else
#ifdef ACTDATA_DOUBLE
    fit->unsmooth_ratio[i]=sqrt(cblas_ddot(tod->ndata,dataCopy[i],1,dataCopy[i],1)/
                                cblas_ddot(tod->ndata,tod->data[i],1,tod->data[i],1));
#else
    fit->unsmooth_ratio[i]=sqrt(cblas_sdot(tod->ndata,dataCopy[i],1,dataCopy[i],1)/
                                cblas_sdot(tod->ndata,tod->data[i],1,tod->data[i],1));
#endif
#endif
  }
    
  psFree(dataCopy[0]);
  psFree(dataCopy);
  fit->have_unsmoothed_ratio=true;
  return;
}

/*---------------------------------------------------------------------------------------------------------*/
/// C version of the python mbCutUnsmoothDets, which is mostly an interface to 
/// mbCalculateUnsmoothRatio
/// 
/// 
void nkCutUnsmoothDets(mbTOD *tod, mbNoiseCommonMode *fit, actData maxErr, actData tUnsmooth)
{
  assert(tod);
  assert(fit);
  assert(tod->data);

  fit->t_unsmooth=tUnsmooth;
  mbCalculateUnsmoothRatio(tod,fit);
  actData *vec=vector(tod->ndata);
  int nkept=0;
  for (int i=0;i<tod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(tod->cuts, tod->rows[i],tod->cols[i])) {
      vec[nkept]=fit->unsmooth_ratio[i];
      nkept++;
    }    
  }
  actData thresh=compute_median(nkept, vec)+maxErr*compute_median_scat(nkept,vec);
  free(vec);
  int ncut=0;
  for (int i=0;i<tod->ndet;i++) {
    if (fit->unsmooth_ratio[i]>thresh) {
      if (!mbCutsIsAlwaysCut(tod->cuts, tod->rows[i],tod->cols[i])) {
	mbCutsSetAlwaysCut(tod->cuts, tod->rows[i], tod->cols[i]);
	ncut++;
	//printf("cutting detector %d %d because of unsmoothness.\n",tod->rows[i],tod->cols[i]);
      }
    }
  }
  printf("cut a total of %d additional detectors because of unsmoothness.\n",ncut);
  
}
/************************************************************************************************************/
/*!
 *  Invert a positive definite matrix.  Return value is the row at which it failed inside of lapack.  
 */
int  mbInvertPosdefMat(actData **mat, int n)
{
#if 1
  return invert_posdef_mat(mat,n);
#else
  int info;
#ifdef ACTDATA_DOUBLE
  clapack_dpotrf('u', n, mat[0], n, &info);
#else
  clapack_spotrf('u', n, mat[0], n, &info);
#endif
  if (info)
    return info;
  //assert(info==0);  //should be a bit gentler, oh well, if one needs to be...
#ifdef ACTDATA_DOUBLE
  clapack_dpotri('u', n, mat[0], n, &info);
#else
  clapack_spotri('u', n, mat[0], n, &info);
#endif
  //assert(info==0);
  if (info)
    return info;
  for (int i = 0; i < n; i++)
    for (int j = i+1; j < n; j++)
      mat[i][j] = mat[j][i];
  return info;
#endif
}

/************************************************************************************************************/
/*!
 *  Cut detectors that aren't suitably correlated with the common mode.
 */

 
void nkCutUncorrDets(mbTOD *tod, mbNoiseCommonMode *fit,actData maxErr)
{
  assert(tod);
  assert(fit);
  assert(tod->ndet==fit->ndet);
  FILE *outfile=fopen("frac_err.txt","w");
  for (int i=0;i<tod->ndet;i++)
    fprintf(outfile,"%16.8e\n",fit->frac_errs[i]);
  fclose(outfile);
  if (tod->cuts==NULL)
    return;
  actData median_err=compute_median_inplace(tod->ndet,fit->frac_errs);
  actData *vec=vector(tod->ndet);
  for (int i=0;i<tod->ndet;i++) 
    vec[i]=fabs(fit->frac_errs[i]-median_err);
  actData median_scat=compute_median(tod->ndet,vec);
  int ncut=0;
  for (int i=0;i<tod->ndet;i++) {
    if (fit->frac_errs[i]>median_err+maxErr*median_scat) {
      //printf("cutting detector %3d %2d %2d for being uncorrelated.\n",i,tod->rows[i],tod->cols[i]);
      mbCutsSetAlwaysCut(tod->cuts, tod->rows[i], tod->cols[i]);
      ncut++;
    }
  }
  printf("cut a total of %d detectors for being uncorrelated.\n",ncut);
  free(vec);
  
}

