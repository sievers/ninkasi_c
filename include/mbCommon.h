////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file mbCommon.h
/// Defines structure and declares operations for computing and using common-mode signals.
/// Author: Jon Sievers, CITA.
////////////////////////////////////////////////////////////////////////////////////////////////////

#if !defined(MB_COMMON_H)
#define MB_COMMON_H

#if 0
#include "act.h"
#include "mbTOD.h"
#include "mbCuts.h"
#ifndef NO_FFTW
#include <fftw3.h>
#endif
#include "sys/time.h"
#include "mbNoise.h"
#include <nk_clapack.h>
#endif

#include <ninkasi.h>


/************************************************************************************************************/

#define MB_NOISE_DEFAULT_NP_COMMON 2  ///< 
#define MB_NOISE_DEFAULT_NP_POLY 3    ///< 

/// Structure to hold common-mode info.

typedef struct {
  actData *common_mode;   ///< The common mode vector, computed as the median/mean unit-scatter/RMS values.
  actData *common_mode_nocalbol; ///< Common mode vector with calbols smoothed over.
  actData calbol_pad;     ///< Padding (seconds) on either side of a detected calbol pulse to treat as pulse.
  actData *median_vals;   ///< Median value of incoming tods (one median per detector).
  actData *median_scats;  ///< Median scatter of a tod (abs value of its distance) about its median value
  actData **fit_params;   ///< Best-fit paramters describing the common-mode plus polynomial to the tod's
  actData *frac_errs;     ///< RMS fractional scatter after subtracting off the best-fit common mode.

  bool apply_common;      ///< Whether to apply the common-mode fit to TODs.
  bool common_is_applied; ///< Whether common-mode fit has been applied to the TOD.
  bool keep_vecs;         ///< Not currently useful.

  actData **vecs;         ///< Precomputed vectors of time and powers thereof for Common-mode fitting.
  bool have_vecs;         ///< Whether vecs are computed.

  actData **data;         ///< Spot for temporary scratch copy of the full TOD data set.
  bool have_data;         ///< Whether data is allocated and filled.

  actData **ata;          ///< Used in parameter fitting??
  bool have_ata;          ///< Whether ata is allocated and filled. 

  actData **atx;          ///< Used in parameter fitting??
  bool have_atx;          ///< Whether atx is allocated and filled. 

  actData *weights;       ///< Weight to use when combining timestreams into a common mode
  bool have_weights;

  int ndet;               ///< Number of detectors in companion TOD.
  int ndata;              ///< Number of data points in companion TOD.
  int nparam;             ///< Number of total parameters to fit.
  int np_common;          ///< Number of polynomial times common-mode coefficients to fit.
  int np_poly;            ///< Number of polynomial coefficients to fit
  int ncalbol;            ///< Number of calbol pulses found.
  int *calbol_start;      ///< Array of data indices to start of calbol pulses.
  int *calbol_stop;       ///< Array of data indices to end of calbol pulses.
  int *icutvec;           ///< Array ndata long with whether each data index was during a calbol.
  actData dt;             ///< Average time between data samples (seconds).

  actData nsig;           ///< Glitch-finding parameters for locating calbols: significance of glitch.
  actData tSmooth;        ///< Glitch-finding parameters for locating calbols: smoothing time (sec).
  actData tGlitch;        ///< Glitch-finding parameters for locating calbols: glitch time (sec).
  actData *unsmooth_ratio;///< Ratio of each detector's data RMS with smoothing to without.
  actData t_unsmooth;     ///< Time scale for smoothing (sec) when computing unsmooth ratio.
  bool have_unsmoothed_ratio;///< Whether unsmooth_ratio has been computed.

} mbNoiseCommonMode;



/************************************************************************************************************/


extern void sgemv_(const char *trans, const int *m, const int *n, const float *alpha, const float *a, 
                   const int *lda, const float *x, const int *incx, const float *beta, const float *y,
                   const int *incy, int translen);


mbNoiseCommonMode *mbAllocateCommonMode(const mbTOD *tod, int np_common, int np_poly);
void mbCleanUpCommonMode(mbNoiseCommonMode *data);
void mbCalculateCommonMode(mbNoiseCommonMode *fit);
void mbCalculateMeanCommonMode(mbNoiseCommonMode *fit, const mbTOD *tod, const mbCuts *cuts);
int mbGetNoCutInds(const mbCuts *cuts, int ndata, int row, int col, int **istart_out, int **istop_out);
void mbRescaleArray(const mbTOD *tod, mbNoiseCommonMode *fit);
void mbRescaleArrayMean(const mbTOD *tod, mbNoiseCommonMode *fit, const mbCuts *cuts);
void mbFindCalbols(mbNoiseCommonMode *fit);
void mbAllocateCalbols(mbNoiseCommonMode *fit,int ncalbol);
void mbCalculateCommonModeVecs(mbNoiseCommonMode *fit);
void mbCalculateCommonModeFitParams(mbNoiseCommonMode *fit);
void mbApplyCommonMode(mbTOD *tod, mbNoiseCommonMode *fit,bool cutCalbols);
void mbApplyCutsToCommonMode(const mbTOD *tod, mbNoiseCommonMode *fit, const mbCuts *cuts);
void mbCalculateCommonFracErrs(const mbTOD *tod, mbNoiseCommonMode *fit);
void mbCalculateUnsmoothRatio(mbTOD *tod, mbNoiseCommonMode *fit);

mbNoiseCommonMode *mbReadCommon( const char *filename );
void mbWriteCommon( const char *filename, mbNoiseCommonMode *common );

void mbNoiseCommonModeFree( mbNoiseCommonMode *cm );
int  mbInvertPosdefMat(actData **mat, int n);

void nkCutUncorrDets(mbTOD *tod, mbNoiseCommonMode *fit,actData maxErr);
void nkCutUnsmoothDets(mbTOD *tod, mbNoiseCommonMode *fit, actData maxErr, actData tUnsmooth);

#endif
