#include <ninkasi.h>
#include <mbTOD.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "ninkasi_mathutils.h"
/*--------------------------------------------------------------------------------*/
int ind_from_rowcol(mbTOD *tod, int row, int col)
{
  for (int i=0;i<tod->ndet;i++)
    if ((tod->rows[i]==row)&&(tod->cols[i])==col)
      return i;
  return -1;
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

static void filter_one_tod(actData *to_filt, const actData *filt, actComplex *fftspace, 
                           act_fftw_plan p_forward, act_fftw_plan p_back, int n)
{
  act_fftw_execute_dft_r2c(p_forward,to_filt, fftspace);
  for (int i=0;i<n/2+1;i++) {
    fftspace[i][0]*=filt[i];
    fftspace[i][1]*=filt[i];
  }
  act_fftw_execute_dft_c2r(p_back,fftspace,to_filt);
  
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
  //assert(tod->dt!=NULL);
  //actData dt=(tod->dt[tod->ndata-1]-tod->dt[0])/((actData)tod->ndata);
  assert(tod->deltat>0);
  actData *filt=calculate_glitch_filterC(0.0, t_smooth,tod->deltat,tod->ndata, 1);
  for (int i=0; i<tod->ndata; i++)
    filt[i]/=tod->ndata;

  filter_whole_tod(tod,filt,to_filt_in);
  psFree(filt);
  
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
 * Apply a filter to a whole TOD.  If you pass in a non-NULL in to_filt_in, 
 * filter that instead and store the result into the TOD.
 * \param tod        The TOD to filter.
 * \param filt       The filter (transfer function) to apply.
 * \param to_filt_in The data set to filter in place of the TOD.  If not NULL, it must be of length
 *                   tod->ndet each entry having length tod->ndata.
 */

void filter_whole_tod(mbTOD *tod, const actData *filt, actData **to_filt_in)
{

#if 1  //use this to avoid problems below wherein it seems that some stuff is incompatible with also using fft_all_data with MKL ffts
  if (to_filt_in==NULL) {
    apply_real_filter_to_data(tod,filt);
#pragma omp parallel for shared(tod) default(none)
    for (int i=0;i<tod->ndet;i++)
      for (int j=0;j<tod->ndata;j++)
	tod->data[i][j]*=tod->ndata;
    return;
  }
#endif

  actData **to_filt;    // The list of vectors to filter.
  if (to_filt_in==NULL)
    to_filt=tod->data;
  else
    to_filt=to_filt_in;
  

#pragma omp parallel shared(tod,to_filt,filt) default (none)
  {
    int n=tod->ndata;
    actData *myfilt=(actData *)psAlloc(n*sizeof(actData));
    memcpy(myfilt,filt,sizeof(actData)*n);
    actComplex *tmpfft=(actComplex *)act_fftw_malloc((n/2+1)*sizeof(actComplex));
    
    act_fftw_plan p_forward;
    act_fftw_plan p_back;
#pragma omp critical
    {
      p_forward =act_fftw_plan_dft_r2c_1d(n,to_filt[0],tmpfft,FFTW_ESTIMATE);  
      p_back    =act_fftw_plan_dft_c2r_1d(n,tmpfft,to_filt[0],FFTW_ESTIMATE);  
    }
    
#pragma omp for
    for (int i=0; i<tod->ndet; i++) {
      filter_one_tod(to_filt[i],filt,tmpfft,p_forward,p_back,n);
    }
#pragma omp critical
    {
      act_fftw_destroy_plan(p_forward);
      act_fftw_destroy_plan(p_back);
    }
    act_fftw_free(tmpfft);
    psFree(myfilt);
#if 0
    actData nn=n;
#pragma omp for
    for (int i=0;i<tod->ndet;i++)
      for (int j=0;j<n;j++)
	tod->data[i][j]/=nn;
#endif
  }
}



/*--------------------------------------------------------------------------------*/


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
  const actData *filt,     ///< The pre-calculated smoothing filter to use.
  actComplex *tmpfft, ///< Pre-allocated vector of length n.  Returns DFT of smoothed data..
  actData *tmpvec,         ///< Pre-allocated vector of length n.  Returns unmodified input.
  actData  *tmpclean,       ///< Pre-allocated vector of length n.  Returns abs(raw-smooth).
  int *cutvec,           ///< Pre-allocated vector of length n.  Returns cut list if (do_cuts).
  int n,                 ///< Length of the input data vector mydat.
  //fftwf_plan p_forward,  ///< The FFTW plan for going to the frequency domain.
  //fftwf_plan p_back,     ///< The FFTW plan for going back to the time domain.
  act_fftw_plan p_forward,  ///< The FFTW plan for going to the frequency domain.
  act_fftw_plan p_back,     ///< The FFTW plan for going back to the time domain.
  bool do_smooth,        ///< Whether to replace the entire data vector with its smoothed value.
  bool apply_glitch,     ///< Whether to replace data exceeding the cut threshold with its smoothed value.
  bool do_cuts,          ///< Whether to build the cutvec of data failing the cut threshhold test.
  actData nsig)          ///< The cut threshold is this factor times the median abs deviation of smooths.
{
  // If we have both set to true, we don't know what to put in mydat
#ifdef HAVE_PSERR
  if(do_smooth && apply_glitch)
    psError(MB_ERR_BAD_VALUE, true, 
            "Cannot both smooth and deglitch (=not smooth) data in glitch_one_detector\n");
#endif

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
void highpass_tod(mbTOD *tod, actData nu_low, actData nu_high)
//  Filter the data, zapping frequencies below nu_low, and ramping up to no filtering at nu_high
{
  assert(tod->have_data);
  actData *filt=vector(tod->ndata);
  actData dnu=1.0/(tod->deltat*tod->ndata);
  for (int i=0;i<tod->ndata;i++) {
    actData nu=i*dnu;
    if (nu<nu_low) {
      filt[i]=0;
    }
    else {
      if (nu>nu_high)
	filt[i]=1;
      else {
	filt[i]=(nu-nu_low)/(nu_high-nu_low);
      }
      
    }
    
  }

  
  actData nn=tod->ndata;
  for (int i=0;i<tod->ndata;i++)
    filt[i]/=nn;
      
  filter_whole_tod(tod,filt,NULL);


  free(filt);
  return;
  
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
#ifdef HAVE_PSERR
  if(do_smooth && apply_glitch)
    psError(MB_ERR_BAD_VALUE, true, 
            "Cannot both smooth and deglitch (=not smooth) data in glitch_one_detector\n");
#endif

  actComplex *tmpfft=act_fftw_malloc(sizeof(actComplex)*n);
  actData *tmpvec=psAlloc(n*sizeof(actData));
  actData *tmpclean=psAlloc(n*sizeof(actData));
  actData *filtvec=calculate_glitch_filterC(t_glitch,t_smooth,dt,n,1);
  int *cutvec=psAlloc(n*sizeof(int));
      
  act_fftw_plan p_forward;
  act_fftw_plan p_back;
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
  {
    p_forward =act_fftw_plan_dft_r2c_1d(n,tmpvec,tmpfft,FFTW_ESTIMATE);  
    p_back    =act_fftw_plan_dft_c2r_1d(n,tmpfft,tmpvec,FFTW_ESTIMATE);  
  }
  glitch_one_detector(mydat,filtvec, tmpfft,tmpvec,tmpclean,cutvec,n, p_forward,  p_back,do_smooth,
                      apply_glitch, do_cuts,nsig);

      
#if !defined(MB_SKIP_OMP)
#pragma omp critical
#endif
  {
    act_fftw_destroy_plan(p_forward);
    act_fftw_destroy_plan(p_back);
  }
  psFree(tmpvec);
  psFree(tmpclean);
  psFree(filtvec);
  act_fftw_free(tmpfft);
  
  if (do_cuts)
    return cutvec;
  else {
    psFree(cutvec);
    return NULL;
  }
}



/*--------------------------------------------------------------------------------*/
actData *CopyData1Det(mbTOD *tod,int det)
{
  assert(tod->have_data);
  size_t nbyte=tod->ndata*sizeof(actData);
  actData *vec=psAlloc(nbyte);
  memcpy(vec,tod->data[det],nbyte);
  return vec;
}

/*--------------------------------------------------------------------------------*/

int *glitch_all_detectors_simple(mbTOD *tod, bool do_smooth, bool apply_glitch, bool do_cuts, actData nsig, actData t_glitch, actData t_smooth, int filt_type)
{
  assert(tod);
  assert(tod->have_data);
#pragma omp parallel for shared(tod,do_smooth,apply_glitch,do_cuts,nsig,t_glitch,t_smooth,filt_type) default(none) schedule(dynamic,1)  
  for (int i=0;i<tod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
      int *cuts=glitch_one_detector_simple(tod->data[i],tod->ndata, do_smooth, apply_glitch,do_cuts,nsig,t_glitch,t_smooth,tod->deltat,filt_type);
      if (do_cuts)
	free(cuts);
    }
  }
  return NULL;
}
