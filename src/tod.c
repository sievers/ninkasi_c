/*!
 * \file tod.c  Supports mbSample and mbTOD objects.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifndef NO_FFTW
#include <fftw3.h>
#endif

#include "act.h"
#include "mbTOD.h"
#include "mbCuts.h"
#include "mbErrorCodes.h"

/************************************************************************************************************/
// Local functions and constants.
/************************************************************************************************************/

/*!
 * Free a ptr to an array of ptrs, along with what it points to.
 * \param p  The ptr to ptrs.
 */

static inline void psFree2d(void **p) {
  if (p)
    psFree(*p);
  psFree(p);
}



/************************************************************************************************************/
// mbSample methods.
/************************************************************************************************************/

/*!
 * Deallocator for the mbSample object.
 * \param samp  Object to be freed.
 */

static void sampleFree(mbSample *samp)
{
  psFree(samp->data);
}



/*!
 * Allocate an mbSample object with a specified data type.
 * \return The new object.
 */

mbSample *mbSampleAllocType(int nrow,   //!< Number of rows
                            int ncol,   //!< Number of columns
                            actDataTypes dataType) //<! Type code for data.
{
  mbSample *samp = psAlloc(sizeof(mbSample));

  samp->tv_sec = 0;             // time (seconds since Unix Epoch)
  samp->tv_usec = 0;            // time (microseconds; 0 < tv_usec < 1e6)
  samp->az = 0;			// azimuth of boresight; degrees
  samp->alt = 0;		// altitude of boresight; degrees
    
  samp->data = actImageAllocType(nrow, ncol, dataType);

  psMemSetDeallocator(samp, (psFreeFunc)sampleFree);

  return samp;
}



/*!
 *  Allocate an mbSample object, using the default act data type.
 * \return The new object.
 */

mbSample *mbSampleAlloc(int nrow,	//!< Number of rows
			int ncol)	//!< Number of columns
{
  return mbSampleAllocType(nrow, ncol, ACT_TYPE_DEFAULT);
}



/************************************************************************************************************/
// mbTOD methods.
/************************************************************************************************************/

/*!
 * Deallocator for the mbTOD object.
 * \param tod  Object to be freed.
 */

static void TODFree(mbTOD *tod)
{
  psFree(tod->tv_sec);
  psFree(tod->tv_usec);
  psFree(tod->dt);
  psFree(tod->lmst);
    
  psFree(tod->az);
  psFree(tod->alt);
  psFree(tod->ra);
  psFree(tod->dec);

  psFree(tod->cols);
  psFree(tod->rows);
  psFree(tod->dark);

  psFree2d((void **)(tod->dets));
  psFree2d((void **)(tod->data));

  psFree(tod->cuts);
  psFree(tod->corrs[0]);
  psFree(tod->corrs);

  psFree(tod->pointingOffset);
  psFree(tod->detPtAz);
  psFree(tod->detPtAlt);
  psFree(tod->detPtRa);
  psFree(tod->detPtDec);
}



/*!
 * Allocate an mbTOD object.
 * \return  The new object.
 */

mbTOD *mbTODAlloc(int nrow,        //!< Number of rows in TES array
                  int ncol,        //!< Number of columns
                  int ndata,       //!< Number of data points
                  int ndet,        //!< Number of pixels for which to allocate memory 
                  const int *rows, //!< List of row numbers for pixels to be stored
                  const int *cols, //!< List of column numbers for pixels to be stored
                  actFilter band)
{
  assert (ndata > 0);
  assert(ndet > 0);

  mbTOD *tod = psAlloc(sizeof(mbTOD));

  tod->band = band; 

  tod->ndata_alloc = tod->ndata = ndata;
  tod->nrow = nrow;
  tod->ncol = ncol;
  tod->ndark = 0;

  tod->times_absent = false;
  tod->point_absent = false;
  tod->tv_sec = psAlloc(ndata*sizeof(psS32));
  tod->tv_usec = psAlloc(ndata*sizeof(psS32));
  tod->dt = NULL;
  tod->sampleTime = ACT_NO_VALUE;
  tod->lmst = NULL;
  tod->lat = ACT_NO_VALUE;
  tod->lon = ACT_NO_VALUE;

  tod->az = psAlloc(ndata*sizeof(float));
  tod->alt = psAlloc(ndata*sizeof(float));

  tod->ra = NULL;
  tod->dec = NULL;

  tod->data = psAlloc(ndet*sizeof(actData *));
  tod->data[0] = psAlloc(ndet*ndata*sizeof(actData));
  for (int i = 0; i < ndet; i++){
    tod->data[i] = tod->data[0] + i*ndata;
  }

  tod->ends = psAlloc(ndata*sizeof(actData *));
  tod->ends[0] = psAlloc(3*ndata*sizeof(actData));
  for (int i = 0; i < ndata; i++){
    tod->ends[i] = tod->ends[0] + 2*i;
    tod->ends[i][0] = 0.0;
    tod->ends[i][1] = 0.0;
  }
  tod->detrended = psAlloc(ndata*sizeof(int));
  for (int i = 0; i < ndet; i++) tod->detrended[i] = 0;

  tod->dark = psAlloc(ndet*sizeof(int));
  tod->cols = psAlloc(ndet*sizeof(int));
  tod->rows = psAlloc(ndet*sizeof(int));
  tod->ndet_alloc = tod->ndet = ndet;
  tod->dets = psAlloc(nrow*sizeof(int *));
  tod->dets[0] = psAlloc(ncol*nrow*sizeof(int));
  for (int i = 0; i < nrow; i++){
    tod->dets[i] = tod->dets[0] + i*ncol;
  }
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      tod->dets[i][j] = ACT_NO_VALUE;
    }
  }
  for (int i = 0; i < ndet; i++){
    tod->cols[i] = cols[i];
    tod->rows[i] = rows[i];
    tod->dets[rows[i]][cols[i]] = i;
  }

  tod->corrs=psAllocMatrix(tod->ndet,tod->ndet);
  memset(tod->corrs[0],0,tod->ndet*tod->ndet*sizeof(actData));

  tod->cuts = (struct mbCuts*) mbCutsAlloc( tod->nrow, tod->ncol );
  tod->pointingOffset = NULL;
  tod->detPtAlt = NULL;
  tod->detPtAz = NULL;
  tod->detPtRa = NULL;
  tod->detPtDec = NULL;

  psMemSetDeallocator(tod, (psFreeFunc)TODFree);

  return tod;
}

void
mbTODDarkAlloc(mbTOD *tod, const int *dets, int ndet)
{
  tod->ndark = ndet;
  for (int i = 0; i < ndet; i++) tod->dark[i] = dets[i];
}

void
mbTODDataAlloc( mbTOD *tod )
{
    if ( tod->data != NULL || tod->ndata_alloc != 0 || tod->ndet_alloc != 0 )
        return;

    tod->data = psAlloc(tod->ndet*sizeof(actData *));
    tod->data[0] = psAlloc(tod->ndet*tod->ndata*sizeof(actData));
    for ( int i = 0; i < tod->ndet; i++ )
        tod->data[i] = tod->data[0] + i*tod->ndata;

    tod->ndet_alloc = tod->ndet;
    tod->ndata_alloc = tod->ndata;
}

void
mbTODDataFree( mbTOD *tod )
{
    psFree2d((void **)(tod->data));
    tod->data = NULL;
    tod->ndet_alloc = 0;
    tod->ndata_alloc = 0;
}

void
mbTODDataClear( mbTOD *tod )
{
    assert( tod->data != NULL && tod->ndata_alloc > 0 && tod->ndet_alloc > 0 );
    memset( tod->data[0], '\0', tod->ndata_alloc*tod->ndet_alloc*sizeof(actData) );
}

int
mbTODGetNumberOfDetectors( const mbTOD *tod )
{
    return tod->ndet;
}

void
mbTODGetDetectorRowCol( const mbTOD *tod, int idet, int *row, int *col )
{
    if ( idet < 0 || idet >= tod->ndet )
    {
        *row = ACT_NO_VALUE;
        *col = ACT_NO_VALUE;
        return;
    }

    *row = tod->rows[idet];
    *col = tod->cols[idet];
}

/*!
 * Replace any former pointing offset with a new one.
 * \param tod  The TOD to alter.
 * \param po   The PointingOffset to keep.
 */

void mbTODSetPointingOffset( mbTOD *tod, 
                             struct mbPointingOffset *po)
{
  psFree( tod->pointingOffset );
  tod->pointingOffset = psMemIncrRefCounter( po );
}



/*!
 * Return whether this TOD has a pointingOffset object associated.
 */

bool mbTODHasPointingOffset( const mbTOD *tod )
{
  if (tod->pointingOffset != NULL){
    return true;
  } else {
    return false;
  }
}
                             

/*!
 * Given a TOD and a listed subset of detectors, reallocate the TOD data
 * such that it only contains the data from that subset of detectors.
 *
*/

void mbTODShrink( mbTOD *tod,        //!< the tod to operate on
                  int ndet,          //!< number of detectors for which to allocate memory
                  const int *rows,   //!< row numbers for detectors to be stored
                  const int *cols)   //!< column numbers for detectors to be stored
{
  int ndata, nrow, ncol;
  actData **data;
  int *new_rows, *new_cols;

  nrow  = tod->nrow;
  ncol  = tod->ncol;
  ndata = tod->ndata;

  data = psAlloc(ndet*sizeof(actData *));
  data[0] = psAlloc(ndet*ndata*sizeof(actData));
  for (int i = 0; i < ndet; i++){
    data[i] = data[0] + i*ndata;
    if (tod->dets[rows[i]][cols[i]] != ACT_NO_VALUE){
      memcpy(data[i],tod->data[tod->dets[rows[i]][cols[i]]], 
             ndata*sizeof(actData));
    }
  }

  new_cols = psAlloc(ndet*sizeof(int));
  new_rows = psAlloc(ndet*sizeof(int));
  ndet = ndet;

  //Reassign detectors (don't need to de- or re-allocate)
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      tod->dets[i][j] = ACT_NO_VALUE;
    }
  }
  for (int i = 0; i < ndet; i++){
    new_cols[i] = cols[i];
    new_rows[i] = rows[i];
    tod->dets[rows[i]][cols[i]] = i;
  }

  psFree(tod->data[0]);
  psFree(tod->data);
  psFree(tod->cols);
  psFree(tod->rows);
    
  tod->data = data;
  tod->rows = new_rows;
  tod->cols = new_cols;
  tod->ndet = ndet;
}


/*!
 *  Remove a mode from a detector's TOD.
 * \param tod    The mbTOD to change.
 * \param mode   The vector to subtract from the data.
 * \param nmode  The length of the mode vector (must equal the TOD length).
 * \param row    The detector row number.
 * \param col    The detector column number.
 */

void mbTODRemoveMode( mbTOD *tod, const float *mode, int nmode, int row, int col )
{
  assert(nmode == tod->ndata);
  int det = tod->dets[row][col];
  for (int i = 0; i < tod->ndata; i++)
    tod->data[det][i] -= mode[i];
}


/*!
 *  Calibrate a detector's TOD.  This means rescale by a constant factor.
 * \param tod    The mbTOD to change.
 * \param cal    The factor to apply.
 * \param row    The detector row number.
 * \param col    The detector column number.
 */

void mbTODCalibrate( mbTOD *tod, float cal, int row, int col )
{
  int det = tod->dets[row][col];
  for (int i = 0; i < tod->ndata; i++)
    tod->data[det][i] *= cal;
}



/*!
 * Change the number of data points in an mbTOD object.
 *
 * \note{Currently this can only change the number of data points used, it cannot extend the arrays}
 */

void mbTODRealloc(mbTOD *tod,		//!< the TOD
		  const int ndata)	//!< how long the data should be said to be
{
  assert (tod != NULL);
  assert (ndata >= 0);
  if (ndata > tod->ndata_alloc) {
    psAbort("PS_FILE_LINE", "mbTODRealloc can currently only adjust ndata within allocated range\n");
  }
  tod->ndata = ndata;
}




/*--------------------------------------------------------------------------------*/

/*!
 * Remove a linear trend between the first and last data values for several detectors.
 * \param tod   The TOD to operate on.
 * \param dets  The list of detector IDs to use
 * \param ndet  The length of the list dets.
 */

void mbDetrendDataC(mbTOD *tod, const int *dets, int ndet)
{
  int win = 1000;
  if (win > tod->ndata) win = tod->ndata;

#if !defined(MB_SKIP_OMP)
#pragma omp parallel for shared(tod,dets,ndet,win)
#endif

  for (int dd=0; dd<ndet; dd++)
  {
    int det=dets[dd];
    if (tod->data[det] && tod->detrended[det] == 0)
    {
      double x0 = 0.0;
      double x1 = 0.0;
      for (int i = 0; i < win; i++) {
        x0 += (double)tod->data[det][i];
        x1 += (double)tod->data[det][tod->ndata-1 - i];
      }
      x0 /= (double)win;
      x1 /= (double)win;

      float m = (x0+x1)/2.0;
      tod->ends[det][0] = (float)x0 - m;
      tod->ends[det][1] = (float)x1 - m;

      actData dx = (actData)(x1-x0);
      actData tmax = (actData)(tod->ndata-1.0);
      actData trendSlope = dx/tmax;
    
      for (int j = 0; j < tod->ndata; j++)
        tod->data[det][j] -= (x0 - m + trendSlope*((actData) j));
      tod->detrended[det] = 1;
    }
  }
}

/*--------------------------------------------------------------------------------*/

/*!
 * Recover the original linear trend between the first and last data values for several detectors.
 * \param tod   The TOD to operate on.
 * \param dets  The list of detector IDs to use
 * \param ndet  The length of the list dets.
 */

void mbRetrendDataC(mbTOD *tod, const int *dets, int ndet)
{

#if !defined(MB_SKIP_OMP)
#pragma omp parallel for shared(tod,dets,ndet)
#endif

  for (int dd=0; dd<ndet; dd++)
  {
    int det=dets[dd];

    if (tod->data[det] && tod->detrended[det] == 1)
    {
      actData dx = tod->ends[det][1] - tod->ends[det][0];
      actData tmax = (actData)(tod->ndata-1.0);
      actData trendSlope = dx/tmax;
    
      for (int j = 0; j < tod->ndata; j++)
        tod->data[det][j] += (tod->ends[det][0] + trendSlope*((actData) j));
      tod->ends[det][0] = 0.0;
      tod->ends[det][1] = 0.0;
      tod->detrended[det] = 0;
    }
  }
}


/*--------------------------------------------------------------------------------*/
/*!
 * Filter a TOD with a give Fourier-domain transfer function.
 * \param to_filt   The data vector to be filtered.  Modified and returns the filtered result.
 * \param filt      The filter (transfer function) to apply.
 * \param fftspace  An allocated vector that returns the DFT of the filtered data.
 * \param p_forward The FFTW plan for going to the frequency domain.
 * \param p_back    The FFTW plan for going back to the time domain.
 * \param n         The length of the data set.
 */

static void filter_one_tod(actData *to_filt, const actData *filt, fftwf_complex *fftspace, 
                           fftwf_plan p_forward, fftwf_plan p_back, int n)
{
  fftwf_execute_dft_r2c(p_forward,to_filt, fftspace);
  for (int i=0;i<n/2+1;i++) {
    fftspace[i][0]*=filt[i];
    fftspace[i][1]*=filt[i];
  }
  fftwf_execute_dft_c2r(p_back,fftspace,to_filt);
  
}



/*--------------------------------------------------------------------------------*/
/*!
 * Apply a filter to a whole TOD.  If you pass in a non-NULL in to_filt_in, 
 * filter that instead and store the result into the TOD.
 * \param tod        The TOD to filter.
 * \param filt       The filter (transfer function) to apply.
 * \param to_filt_in The data set to filter in place of the TOD.  If not NULL, it must be of length
 *                   tod->ndet each entry having length tod->ndata.
 */

void filter_whole_tod(mbTOD *tod, const actData *filt, actData **to_filt_in)
{
  actData **to_filt;    // The list of vectors to filter.
  if (to_filt_in==NULL)
    to_filt=tod->data;
  else
    to_filt=to_filt_in;
  

#if !defined(MB_SKIP_OMP)
#pragma omp parallel shared(tod,to_filt,filt) default (none)
#endif
  {
    int n=tod->ndata;
    actData *myfilt=(actData *)psAlloc(n*sizeof(actData));
    memcpy(myfilt,filt,sizeof(actData)*n);
    fftwf_complex *tmpfft=(fftwf_complex *)fftwf_malloc((n/2+1)*sizeof(fftwf_complex));
    
    fftwf_plan p_forward;
    fftwf_plan p_back;
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
    {
      p_forward =fftwf_plan_dft_r2c_1d(n,to_filt[0],tmpfft,FFTW_ESTIMATE);  
      p_back    =fftwf_plan_dft_c2r_1d(n,tmpfft,to_filt[0],FFTW_ESTIMATE);  
    }
    
#if !defined(MB_SKIP_OMP)
#pragma omp for
#endif
    for (int i=0; i<tod->ndet; i++) {
      filter_one_tod(to_filt[i],filt,tmpfft,p_forward,p_back,n);
    }
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
    {
      fftwf_destroy_plan(p_forward);
      fftwf_destroy_plan(p_back);
    }
    fftwf_free(tmpfft);
    psFree(myfilt);
  }
}



/*--------------------------------------------------------------------------------*/
/*!
 * Smooth a whole TOD on time scale t_smooth.  If to_filt_in is non-NULL, smooth that and store
 * the result into the TOD.
 * \param tod        The TOD to filter.
 * \param t_smooth   The smoothing time (sigma of a Gaussian), in units of tod->dt (presumably: seconds).
 * \param to_filt_in The data set to filter in place of the TOD.  If not NULL, it must be of length
 *                   tod->ndet each entry having length tod->ndata.
 */
void smooth_tod(mbTOD *tod, actData t_smooth, actData **to_filt_in)
{
  assert(tod->dt!=NULL);
  actData dt=(tod->dt[tod->ndata-1]-tod->dt[0])/((actData)tod->ndata);
  actData *filt=calculate_glitch_filterC(0.0, t_smooth,dt,tod->ndata, 1);
  for (int i=0; i<tod->ndata; i++)
    filt[i]/=tod->ndata;

  filter_whole_tod(tod,filt,to_filt_in);
  psFree(filt);
  
}



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

float sselect(unsigned long k, unsigned long n, float *arr)
{
  unsigned long i,ir,j,l,mid;
  float a,temp;

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

/*--------------------------------------------------------------------------------*/
/*!
 * Apply a pre-calculated smoothing filter to a data vector.  Then do one of
 * (1) Just smooth the data if (do_smooth).
 * (2) Preserve the data, replacing with smoothed data only during large glitches if (apply_glitch).
 * (3) Fill a vector cutvec saying whether to cut each point due to glitches if (do_cuts).
 * You can call any combination that doesn't involve both #1 and #2.
 * If apply_glitch is true, then replace anything flagged as over the cut threshold (defined as a
 * smoothed value differing from raw by at least nsig*median absolute deviation of smoothed values)
 * with the smoothed value.  Since the glitch filter should have zero in the middle, height of
 * cosmic ray doesn't matter.
 */

static void glitch_one_detector(
  actData *mydat,        ///< Data vector to filter.  Modified and returned.
  const float *filt,     ///< The pre-calculated smoothing filter to use.
  fftwf_complex *tmpfft, ///< Pre-allocated vector of length n.  Returns DFT of smoothed data..
  float *tmpvec,         ///< Pre-allocated vector of length n.  Returns unmodified input.
  float *tmpclean,       ///< Pre-allocated vector of length n.  Returns abs(raw-smooth).
  int *cutvec,           ///< Pre-allocated vector of length n.  Returns cut list if (do_cuts).
  int n,                 ///< Length of the input data vector mydat.
  fftwf_plan p_forward,  ///< The FFTW plan for going to the frequency domain.
  fftwf_plan p_back,     ///< The FFTW plan for going back to the time domain.
  bool do_smooth,        ///< Whether to replace the entire data vector with its smoothed value.
  bool apply_glitch,     ///< Whether to replace data exceeding the cut threshold with its smoothed value.
  bool do_cuts,          ///< Whether to build the cutvec of data failing the cut threshhold test.
  actData nsig)          ///< The cut threshold is this factor times the median abs deviation of smooths.
{
  // If we have both set to true, we don't know what to put in mydat
  if(do_smooth && apply_glitch)
    psError(MB_ERR_BAD_VALUE, true, 
            "Cannot both smooth and deglitch (=not smooth) data in glitch_one_detector\n");

  // User hasn't requested any action!
  if (! (do_smooth || apply_glitch || do_cuts))
    return;
  
  // Save a copy of the input data, then apply the smoothing filter only.
  memcpy(tmpvec,mydat,n*sizeof(actData));
  filter_one_tod(mydat, filt, tmpfft,p_forward,p_back,n);

  // Joe mystified by this rescaling.
#if 1
  const actData nf=n;
  for (int i=0;i<n;i++)
    mydat[i]/=nf;
#endif

  // If we don't need to know nor apply the cuts, then we're done.
  if (do_smooth & (!do_cuts))
    return;  
  
  // tmpclean is the absolute difference between smooth and raw vectors.  Find its median (for use in cuts).
  for (int j=0;j<n;j++)
    tmpclean[j]=fabs(mydat[j]-tmpvec[j]);
  actData thresh=sselect(n/2,n,tmpclean-1)*nsig;

  if (do_cuts) {
    memset(cutvec,0,sizeof(int)*n);
    for (int j=0;j<n;j++) {
      if (fabs(mydat[j]-tmpvec[j])>thresh)
        cutvec[j]=1;
    }
  }
  if (do_smooth)
    return;

  // By this point, we know the user didn't want smooth data.  Copy data back to the input vector.
  // If apply_glitch, then mostly this will be the raw data, but we leave the smoothed data there
  // during all glitches.  Otherwise, we want the entire raw data set back (and presumably called
  // this function only to get the cuts).
  if (apply_glitch) {
    for (int j=0;j<n;j++)
      if (fabs(mydat[j]-tmpvec[j])<thresh)
        mydat[j]=tmpvec[j];
  } else
    memcpy(mydat,tmpvec,sizeof(actData)*n);
}



/*--------------------------------------------------------------------------------*/
/*!
 * Compute and return a filter appropriate for finding glitches in the TOD.  The filter is given in
 * the Fourier domain and is a Gaussian filter for smoothing, minus a second Gaussian filter for
 * glitches (provided that t_glitch > 0).
 * \note The units of t_glitch, t_smooth, and dt can be seconds or samples or Planck times or
 * anything else, provided that they are all the same units.
 * \param t_glitch   The width (Gaussian sigma) of the glitch time
 * \param t_smooth   The width (Gaussian sigma) of the smoothing time
 * \param dt         Time between data samples.
 * \param n          Length of the desired filter (number of frequency samples).
 * \param filt_type  Select filter type (currently only Gaussian filters are implemented).
 */

actData *calculate_glitch_filterC(const actData t_glitch, const actData t_smooth, const actData dt, 
                                  const int n, const int filt_type)
{
  int nodd=2*(n/2)+1;
  actData *filtvec=psAlloc(nodd*sizeof(actData));

  // All frequencies here are in units of discrete FT samples.  Thus freqSamp=1 corresponds a freq of 1/(N*dt)
  actData freqSamp;
  const actData sig_smooth = n*dt/(2*M_PI*t_smooth); // Width (in samples) of smoothing filter in freq space.

  if (t_glitch <= 0) {
    for (int i=0;i<nodd;i++) {
      if (i>nodd/2)
        freqSamp = i-nodd;
      else
        freqSamp = i;

      filtvec[i]=(exp(-0.5*freqSamp*freqSamp/sig_smooth/sig_smooth));
    }
  } else {
    const actData sig_glitch = n*dt/(2*M_PI*t_glitch); // Width (in samples) of glitch filter in freq space.
    for (int i=0;i<nodd;i++) {
      if (i>nodd/2)
        freqSamp = i-nodd;
      else
        freqSamp = i;

      // Joe--writing comments--doesn't understand this normalization at all!
      filtvec[i]=(exp(-0.5*freqSamp*freqSamp/sig_smooth/sig_smooth)*t_smooth -
                  exp(-0.5*freqSamp*freqSamp/sig_glitch/sig_glitch)*t_glitch)/(t_smooth-t_glitch);
    }
  }
  return filtvec;
}

/*--------------------------------------------------------------------------------*/
/*!
 * Apply a smoothing and de-glitching filter to a data vector.  Then do one of
 * 1. Just smooth the data if (do_smooth).
 * 2. Preserve the data, replacing with smoothed data only during large glitches if (apply_glitch).
 * 3. Fill a vector cutvec saying whether to cut each point due to glitches if (do_cuts).
 * You can call any combination that doesn't involve both #1 and #2.
 * If apply_glitch is true, then replace anything flagged as over the cut threshold (defined as a
 * smoothed value differing from raw by at least nsig*median absolute deviation of smoothed values)
 * with the smoothed value.  Since the glitch filter should have zero in the middle, height of
 * cosmic ray doesn't matter.
 * This function relies on glitch_one_detector, but it computes the requested filter rather than
 * expecing the user to pass it in.
 * \return A vector containing the list of which data samples were cut if (do_cuts).
 */

int *glitch_one_detector_simple(
  actData *mydat,        ///< Data vector to filter.  Modified and returned.
  int n,                 ///< Length of the input data vector mydat.
  bool do_smooth,        ///< Whether to replace the entire data vector with its smoothed value.
  bool apply_glitch,     ///< Whether to replace data exceeding the cut threshold with its smoothed value.
  bool do_cuts,          ///< Whether to build the cutvec of data failing the cut threshhold test.
  actData nsig,          ///< The cut threshold is this factor times the median abs deviation of smooths.
  actData t_glitch,      ///< The width (Gaussian sigma) of the glitch time
  actData t_smooth,      ///< The width (Gaussian sigma) of the smoothing time
  actData dt,            ///< Time between data samples.
  int filt_type)         ///< Select filter type (currently only Gaussian filters are implemented).
{
  // If we have both set to true, we don't know what to put in mydat
  if(do_smooth && apply_glitch)
    psError(MB_ERR_BAD_VALUE, true, 
            "Cannot both smooth and deglitch (=not smooth) data in glitch_one_detector\n");

  fftwf_complex *tmpfft=fftwf_malloc(sizeof(fftwf_complex)*n);
  actData *tmpvec=psAlloc(n*sizeof(actData));
  actData *tmpclean=psAlloc(n*sizeof(actData));
  actData *filtvec=calculate_glitch_filterC(t_glitch,t_smooth,dt,n,1);
  int *cutvec=psAlloc(n*sizeof(int));
      
  fftwf_plan p_forward;
  fftwf_plan p_back;
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
  {
    p_forward =fftwf_plan_dft_r2c_1d(n,tmpvec,tmpfft,FFTW_ESTIMATE);  
    p_back    =fftwf_plan_dft_c2r_1d(n,tmpfft,tmpvec,FFTW_ESTIMATE);  
  }
  glitch_one_detector(mydat,filtvec, tmpfft,tmpvec,tmpclean,cutvec,n, p_forward,  p_back,do_smooth,
                      apply_glitch, do_cuts,nsig);

      
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
  {
    fftwf_destroy_plan(p_forward);
    fftwf_destroy_plan(p_back);
  }
  psFree(tmpvec);
  psFree(tmpclean);
  psFree(filtvec);
  fftwf_free(tmpfft);
  
  if (do_cuts)
    return cutvec;
  else {
    psFree(cutvec);
    return NULL;
  }
}



/*--------------------------------------------------------------------------------*/
/*!
 * Remove glitches from multiple detector data vectors in a TOD using a pre-computed smoothing (and
 * glitch removal) filter. Also, if (do_cuts), then also apply all glitches as new cuts in the cuts
 * object.
 */

void mbGlitchC(
  mbTOD *tod,            ///< The TOD to modify.
  const float *filt,     ///< The pre-calculated smoothing filter to use.
  int nelem,             ///< POSSIBLE BUG: this value is ignored.
  const int *dets,       ///< List of detectors (by index) to clean.
  int ndet,              ///< Number of detectors to clean.
  bool do_smooth,        ///< Whether to replace the data vector with its smoothed value.
  bool apply_glitch,     ///< Whether to replace data exceeding the cut threshold with its smoothed value.
  mbCuts *cuts,   ///< A cuts object of the appropriate size.  Will be updated with cuts on glitches.
  actData nsig,          ///< The cut threshold is this factor times the median abs deviation of smooths.
  int maxGlitch)         ///< Maximum allowed number of glitches.  If exceeded, cut whole detector.
{
  
#if !defined(MB_SKIP_OMP)
#pragma omp parallel shared(tod,filt,nelem,dets,ndet,do_smooth,cuts,apply_glitch,nsig,maxGlitch) default(none)
#endif
  {
    fftwf_plan p_forward;
    fftwf_plan p_back;
    
    const int n=tod->ndata;
    fftwf_complex *tmpfft=fftwf_malloc(sizeof(fftwf_complex)*n);
    float *tmpvec=psAlloc(n*sizeof(actData));
    float *tmpclean=psAlloc(n*sizeof(actData));
    int *cutvec=psAlloc(n*sizeof(int));
    bool do_cuts=false;

    if (cuts)
      do_cuts=true;
    
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
    {
      p_forward =fftwf_plan_dft_r2c_1d(n,tmpvec,tmpfft,FFTW_ESTIMATE);  
      p_back    =fftwf_plan_dft_c2r_1d(n,tmpfft,tmpvec,FFTW_ESTIMATE);  
    }

#if !defined(MB_SKIP_OMP)
#pragma omp for
#endif
    for (int i=0;i<ndet;i++) {
      actData *mydat=tod->data[dets[i]];
      glitch_one_detector(mydat, filt, tmpfft,tmpvec, tmpclean, cutvec ,n,p_forward, p_back,
                          do_smooth, apply_glitch,  do_cuts, nsig);
      if (do_cuts)
        if ( mbCutsExtendByArray((mbCuts *)cuts,tod->rows[dets[i]], 
                                 tod->cols[dets[i]],cutvec,tod->ndata) > maxGlitch ){
          mbCutsSetAlwaysCut((mbCuts *)tod->cuts, tod->rows[dets[i]], tod->cols[dets[i]]);
        }
    }
    
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
    {
      fftwf_destroy_plan(p_forward);
      fftwf_destroy_plan(p_back);
      
      fftwf_free(tmpfft);
      psFree(tmpvec);
      psFree(tmpclean);
      psFree(cutvec);
      
    }
  }
}



/*--------------------------------------------------------------------------------*/
/*!
 * Return the sample number within the TOD that has either the maximum or the minimum value for that
 * detector in the entire TOD.
 * \param tod    The mbTOD to search.
 * \param row    The detector row number.
 * \param col    The detector column number.
 * \param maximum Whether we want to get the maximum
 */

int mbMaxSamp(const mbTOD *tod, int row, int col, bool maximum)
{
  int det=tod->dets[row][col];  

  if (tod->dets[row][col] == ACT_NO_VALUE) {
    return -1;
  }
  
  actData val=tod->data[det][0];
  int spot=0;
  if (maximum) {
    for (int i=0;i<tod->ndata;i++)
      if (tod->data[det][i]>val) {
        val=tod->data[det][i];
        spot=i;
      }
  } else {
    for (int i=0;i<tod->ndata;i++)
      if (tod->data[det][i]<val) {
        val=tod->data[det][i];
        spot=i;
      }      
  }
  
  return spot;

}

/*--------------------------------------------------------------------------------*/
/*!
 * Quick & dirty background subtractor.  De-glitch the tod, then smooth it. Subtract that from the
 * original TOD, and you get a background-removed TOD with cosmic rays/planets left (mostly) intact.
 * We do better things now (like common mode removal), but not a bad way to have a quick look at an
 * effectively de-common-moded single TOD, etc.
 * In the absence of glitches, this process amounts to a high-pass filter.
 */
void mbFlattenC(
  mbTOD *tod,            ///< The TOD to modify.
  const float *filt,     ///< The pre-calculated smoothing filter to use.
  int nelem,             ///< POSSIBLE BUG: this value is ignored.
  const int *dets,       ///< List of detectors (by index) to clean.
  int ndet,              ///< Number of detectors to clean.
  actData nsig)          ///< The cut threshold is this factor times the median abs deviation of smooths.
{
#if !defined(MB_SKIP_OMP)
#pragma omp parallel shared(tod,filt,nelem,dets,ndet,nsig) default(none)
#endif
  {
    fftwf_plan p_forward;
    fftwf_plan p_back;
    
    int n=tod->ndata;
    fftwf_complex *tmpfft=fftwf_malloc(sizeof(fftwf_complex)*n);
    float *tmpvec=psAlloc(n*sizeof(actData));
    float *hold_data=psAlloc(n*sizeof(actData));
    float *tmpclean=psAlloc(n*sizeof(actData));
    
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
    {
      p_forward =fftwf_plan_dft_r2c_1d(n,tmpvec,tmpfft,FFTW_ESTIMATE);  
      p_back    =fftwf_plan_dft_c2r_1d(n,tmpfft,tmpvec,FFTW_ESTIMATE);  
    }
#if !defined(MB_SKIP_OMP)
#pragma omp for
#endif
    for (int i=0;i<ndet;i++)
    {
      actData *mydat=tod->data[dets[i]];
      memcpy(hold_data,mydat,sizeof(actData)*n);
	
      // First, just apply the glitch-remover. Then smooth the result.
      glitch_one_detector(mydat, filt, tmpfft, tmpvec, tmpclean, NULL, n, p_forward, p_back,
                          false, true, false, nsig);
      glitch_one_detector(mydat, filt, tmpfft, tmpvec, tmpclean, NULL, n, p_forward, p_back,
                          true, false, false, nsig);

      // Save the raw data minus the smoothed, de-glitched data.
      for (int j=0;j<n;j++)
        mydat[j]=hold_data[j]-mydat[j];	
    }
    
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
    {
      fftwf_destroy_plan(p_forward);
      fftwf_destroy_plan(p_back);
      
      fftwf_free(tmpfft);
      psFree(tmpvec);
      psFree(hold_data);
      psFree(tmpclean);
    }
  }
}



/*!
 * Apply a time shift in the TOD.
 * \param tod        The TOD to shift.
 * \param shift      The number of samples to shift (shift to earlier times).
 */

void mbTODShiftData( mbTOD *tod, int shift )
{
  if (shift == 0) 
    return;
  if (shift < 0)
    psError(MB_ERR_BAD_VALUE, true, "mbTODShiftData cannot shift by negative number %d\n", shift);

  for ( int i = 0; i < tod->ndet; i++ ) {
    for ( int j = 0; j < tod->ndata - shift; j++ ) {
      tod->data[i][j] = tod->data[i][j+shift];
    }
  }
}




extern void mbGetDetAltAz(float *restrict detAlt, float *restrict detAz, 
                          const struct  mbPointingOffset *offset, int i, int j, 
                          float boresightAlt, float boresightAz);



/// Function for setting the detPt... information describing the pointing
/// of a particular detector in a TOD. The TOD must have a pointingOffset.
/// A psError is raised if pointingOffset is not set.

psErrorCode mbTODSetDetectorPointing( mbTOD *tod, ///<TOD to be changed
                                      int row,    ///<detector row
                                      int col )   ///<detector column
{
  // Make sure we have a pointing offset.
  if (tod->pointingOffset == NULL){
    return psError(MB_ERR_BAD_VALUE, 1, "No pointing offset information");
  }    

  // Assign row and col
  tod->detPtRow = row;
  tod->detPtCol = col;

  // Allocate space for the detector pointing arrays if not done so already
  tod->detPtAlt = (float *) psRealloc( tod->detPtAlt, tod->ndata*sizeof(float) );
  tod->detPtAz = (float *) psRealloc( tod->detPtAz, tod->ndata*sizeof(float) );
  tod->detPtRa = (float *) psRealloc( tod->detPtRa, tod->ndata*sizeof(float) );
  tod->detPtDec = (float *) psRealloc( tod->detPtDec, tod->ndata*sizeof(float) );

  // Fill in the detector's alt and az
  for (int i = 0; i < tod->ndata; i++){
    mbGetDetAltAz( &tod->detPtAlt[i], &tod->detPtAz[i],
                   tod->pointingOffset, row, col, 
                   tod->alt[i], tod->az[i] );
  }

  // Fill in the RA/Dec This can probably be done smarter..., but for now, make up a fake
  // TOD and use its operations to compute RA and dec.
  const int tmpTODcol[1] = {0};
  const int tmpTODrow[1] = {0};
    
  mbTOD *tmpTOD = mbTODAlloc(1, 1, tod->ndata, 1, tmpTODrow, tmpTODcol, tod->band);
  memcpy(tmpTOD->tv_sec, tod->tv_sec, tod->ndata*sizeof(tod->tv_sec[0]));
  memcpy(tmpTOD->tv_usec, tod->tv_usec, tod->ndata*sizeof(tod->tv_usec[0]));
  memcpy(tmpTOD->az, tod->detPtAz, tod->ndata*sizeof(tod->detPtAz[0]));
  memcpy(tmpTOD->alt, tod->detPtAlt, tod->ndata*sizeof(tod->detPtAlt[0]));

  mbTODSetRaDec(tmpTOD, ACT_SITE_LONGITUDE*M_PI/180., 
                ACT_SITE_LATITUDE*M_PI/180.);

  memcpy(tod->detPtRa, tmpTOD->ra, tod->ndata*sizeof(tod->detPtAlt[0]));
  memcpy(tod->detPtDec, tmpTOD->dec, tod->ndata*sizeof(tod->detPtAlt[0]));

  psFree(tmpTOD);
  return PS_ERR_NONE;
}




#if 0

#define MAXLINELENPO MAXLEN

/*--------------------------------------------------------------------------------
/  Code stolen from mbPointingOffset.c
/--------------------------------------------------------------------------------*/



/// Allocator for an mbPointingOffset object.
/// \param nrow  Number of rows in camera.
/// \param ncol  Number of columns in camera.
/// \param fitType  The pointing model used in the fitting.
/// \return The new object

mbPointingOffset *mbPointingOffsetAlloc(int nrow, int ncol, int fitType)  {

  assert((fitType>=0) && (fitType<MB_NPTSOLNS));
  assert((nrow > 0) && (ncol > 0));       //otherwise, there's no data

  mbPointingOffset *offset = psAlloc(sizeof(mbPointingOffset));
  offset->nrow = nrow;
  offset->ncol = ncol;
  offset->offsetAlt = psAlloc(nrow*sizeof(float *));
  offset->offsetAlt[0] = psAlloc(nrow*ncol*sizeof(float));
  offset->offsetAzCosAlt = psAlloc(nrow*sizeof(float *));
  offset->offsetAzCosAlt[0] = psAlloc(nrow*ncol*sizeof(float));
  for (int i = 0; i < nrow; i++)  {
    offset->offsetAlt[i] = offset->offsetAlt[0]+i*ncol;
    offset->offsetAzCosAlt[i] = offset->offsetAzCosAlt[0]+i*ncol;
  }
#if 0  
  offset->fitType = fitType;
  switch (fitType) {
  case MB_FIX_ALT_AZ: // fall through
  case MB_FIX_ALT_AZ_SGEOM:
    offset->nparam = AA_NPARAM;
    break;
    
  case MB_BASE_TILT:
    offset->nparam = BT_NPARAM;
    break; 
  }
  offset->fit = psAlloc((offset->nparam)*sizeof(float));
  
  //psMemSetDeallocator(offset, (psFreeFunc)mbPointingOffsetFree);
#endif
  return offset;
}



/// Read a PointingOffset from a file.
/// \param filename The file to read.
/// \return The new, loaded object.

mbPointingOffset *mbReadPointingOffset(const char *filename)  {
  char line[MAXLINELENPO], *r;
  int nrow, ncol,fitType;
  
  FILE *ifp=fopen(filename,"r");
  if(ifp == NULL)  {
    psError(MB_ERR_IO, true, "Unable to open %s to read pointing offset.",filename);
    return NULL;
  }
  
  r=fgets(line,MAXLINELENPO,ifp);
  if(!(strstr(line,"nrow =")))  {
    psError(MB_ERR_IO, true, "%s has the wrong format.  nrow = ?",filename);
    return NULL;
  }
  if(!(strstr(line,"ncol =")))  {
    psError(MB_ERR_IO, true, "%s has the wrong format.  ncol = ?",filename);
    return NULL;
  }
  
  r=strchr(line,'=');
  assert(r);
  nrow = atoi(r+1);
  r=strrchr(line,'=');
  assert(r);
  ncol = atoi(r+1);
  
  r=fgets(line,MAXLINELENPO,ifp);
  if(!(strstr(line,"fitType =")))  {
    psError(MB_ERR_IO, true, "%s has the wrong format.  fitType = ?",filename);
    return NULL;
  }
  
  r=strchr(line,'=');
  assert(r);
  fitType = atoi(r+1);
  
  mbPointingOffset *offset = mbPointingOffsetAlloc(nrow,ncol,fitType);
  
  switch(offset->fitType)  {
  case(MB_FIX_ALT_AZ): // fall through
  case(MB_FIX_ALT_AZ_SGEOM):
    for (int i=0; i<offset->nparam; i+=2)  {
      fscanf(ifp,"%e %e\n",&(offset->fit[i]),&(offset->fit[i+1]));
    }
    break;
    
  case(MB_BASE_TILT):
    for (int i=0; i<(offset->nparam-9); i+=2)  {
      fscanf(ifp,"%e %e\n",&(offset->fit[i]),&(offset->fit[i+1]));
    }
    fscanf(ifp,"%e %e %e\n",&(offset->fit[BT_11]),&(offset->fit[BT_12]),&(offset->fit[BT_13]));
    fscanf(ifp,"%e %e %e\n",&(offset->fit[BT_21]),&(offset->fit[BT_22]),&(offset->fit[BT_23]));
    fscanf(ifp,"%e %e %e\n",&(offset->fit[BT_31]),&(offset->fit[BT_32]),&(offset->fit[BT_33]));
    break;
    
  }
  
  fscanf(ifp,"\n");
  for (int i=0; i<offset->nrow; i++)  {
    for (int j=0; j<offset->ncol; j++)  {
      fscanf(ifp,"%e %e\n",&(offset->offsetAlt[i][j]),&(offset->offsetAzCosAlt[i][j]));
    }
    fscanf(ifp,"\n");
  }
  fclose(ifp);
  return offset;
}



/*--------------------------------------------------------------------------------
/  End of code stolen from mbPointingOffset.c
/--------------------------------------------------------------------------------*/
#endif
