#if !defined(ACTTOD_H)
#define ACTTOD_H 1

//#include "actConsts.h"
//#include "actImage.h"

#include <ninkasi_defs.h>
//#include <ninkasi.h>
#ifdef ACTPOL
#include <actpol/actpol.h>
#endif


#include <ps_stuff.h>
#include <mbCuts.h>
#include "noise_types.h"

#include "ninkasi_mathutils.h"


typedef char actFits;


/*!
 * A single sample of data taken from one camera.
 */
#if 0
typedef struct {
  psS32 tv_sec;       //!< Time (seconds since Unix Epoch)
  psS32 tv_usec;      //!< Time (microseconds since seond; 0 <= tv_usec < 1e6)
  float az;           //!< Azimuth of boresight; degrees
  float alt;          //!< Altitude of boresight; degrees
  
  actImage *data;     //!< The actual detector data list.
} mbSample;

// From tod.c
mbSample *mbSampleAlloc(int nrow, int ncol);
mbSample *mbSampleAllocType(int nrow, int ncol, actDataTypes dataType);

#endif
/*!
 * A complete TOD; a set of mbSamples, transposed into a set of TOD streams
 */



/// A global pointing solution.
/// Contains a single list of pointing solution parameters, along with the 2D offset data of the
/// individual detectors.


typedef struct {
  int nrow;               ///< Number of detector rows.
  int ncol;               ///< Number of detector columns.
  int nparam;             ///< Number of parameters in the pointing solution model.
  int fitType;            ///< Specify which model you're using for the pointing solution.
  actData **offsetAlt;      ///< Array of angular distance offsets in Alt direction (radians).
  actData **offsetAzCosAlt; ///< Array of angular distance offsets in the Az direction (radians).
  actData *fit;             ///< Best fit parameters and errors for fitType fit.
} mbPointingOffset;

/*--------------------------------------------------------------------------------*/

typedef struct {
  actData *azvec;
  int naz;
  actData *altvec;
  int nalt;
  actData **ra_mat;
  actData **dec_mat;
  actData *ra_clock;
  actData *dec_clock;
  int nclock;

} TiledPointingFit;
/*--------------------------------------------------------------------------------*/
#ifdef ACTPOL
typedef struct {
  int nhorn;

  actData *dx;
  actData *dy;
  actData *theta;
  actData freq;
  ACTpolArray *array;
  actData alt0;
  actData az0;
  actData az_throw;
  ACTpolWeather weather;
  int *ipiv;
  int npiv;
  int dpiv;


  actData **ra_piv;
  actData **dec_piv;
  actData **sin2gamma_piv;
  actData **cos2gamma_piv;

  int n_gamma_az_coeffs;
  actData **gamma_az_cos_coeffs;
  actData *gamma_ctime_cos_coeffs;
  actData **gamma_az_sin_coeffs;
  actData *gamma_ctime_sin_coeffs;

} ACTpolPointingFit;
#endif

/*--------------------------------------------------------------------------------*/

typedef struct {
  actData **timestreams;
  int nstream;
  int n;

} BadTimestreams;

/*--------------------------------------------------------------------------------*/

typedef struct {
  PolyParams2d *ra_fit;
  PolyParams2d *dec_fit;
  actData ra_clock_rate;
  actData dec_clock_rate;
  
  actData ra_offset;  // for use in shifting TOD pointing, for sims etc.
  actData dec_offset;  

  //actData ra_scale;
  //actData dec_scale;

  //this stuff is for trying to subsample quickly
  int ncoarse;
  int *coarse_ind;


  TiledPointingFit *tiled_fit;

} PointingFit;

typedef enum { ACT_AR1=1,
	       ACT_AR2=2,
	       ACT_AR3=3} actArray;

typedef actArray actFilter;



typedef struct {
  int todtype;
  void *generic;
  actFilter band;      //!< enum: ACT_145, ACT_215, or ACT_265
  char *dirfile;

  int ndark;           //!< # of dark or dead detectors 
  int ndet;            //!< # of detectors 
  int ndata;           //!< number of data points
  int ndet_alloc;      //!< number of detectors allocated
  int ndata_alloc;     //!< number of data points allocated

  bool times_absent;   //!< Are the ctime sec/usec absent from the raw data?
  bool point_absent;   //!< Are the az/alt encoder channels absent from the raw data?
  psS32 *tv_sec;       //!< time (seconds since Unix Epoch)
  psS32 *tv_usec;      //!< time (microseconds since second; 0 <= tv_usec < 1e6)
  double *dt;           //!< time (seconds since the start of the scan) 
  float sampleTime;    //!< time (seconds since the start of the scan)
  actData *az;           //!< azimuth of boresight; radians
  actData *alt;          //!< altitude of boresight; radians
  actData *ra;           //!< RA of boresight; radians
  actData *dec;          //!< declination of boresight; radians
  int have_ra;
  int have_dec;

  double ctime;       //!< ctime of first sample
  double deltat;  //!< ctime difference between adjacent samples

  int nrow;            //!< dimensions of array; row
  int ncol;            //!< dimensions of array; column
  actData **data;      //!< the actual pixel data; use as data[det_number][time_index] 
                       //!< where det_number is the detector index into the rows, cols arrays.
  actData **calib_facs_saved;  //2D-array where calibration factors for the data are saved.
  int have_data;
  int decimate;        //!< decimate factor - for each value here, apply a factor of 2 decimation to the data.
  int n_to_window;     //!< how many samples at end to cut/window out.

  actData **corrs;     //!< the detector-detector correlation matrix.
  actData **rotmat;    //!< a rotation matrix, for use in diagonalizing the TOD.
  int *dark;           //!< detector numbers of dark or dead detectors in TOD.
  int *rows;           //!< row numbers of detectors; rows[det]
  int *cols;           //!< column numbers of detectors; cols[det]
  int **dets;          //!< det index matrix; det[row][col]
                       //!< = ACT_NO_VALUE if detector at position 
                       //!< row,col in the array is not stored in
                       //!< the TOD    
  double *lmst;        //!< local mean sidereal time
  double lat;          //!< latitude of observations
  double lon;          //!< longitude of observations
  mbCuts *cuts; //!< cuts that apply to this TOD
  mbUncut ***uncuts;  //!< uncut regions
  mbUncut ***uncuts_for_interp;  //!< uncut regions
  mbUncut ***cuts_as_uncuts;  //!< cut regions stored using uncuts
  mbCutFitParams ***cuts_fit_params;  //legendre polynomial coefficients for gapfilling across the cuts
  
  mbUncut ***cuts_as_vec;  //!< vectorized cut regions for indexing into a global array
  mbUncut ***kept_data; //have a second copy of uncuts in case we wish to project/gapfill different data
  actData **ends;        //!< values at beginning and end of TOD for detrend/retrend
  int *detrended;      //!< is 0 if detector was not detrended and 1 if detrended


  mbPointingOffset *pointingOffset; //!< detector pointing offsets from boresight
  int detPtRow;        ///< row in camera for detPtAz etc.
  int detPtCol;        ///< column in camera for detPtAz etc.
  float *detPtAz;      ///< Azimuth pointing for detector (detPtRow,detPtCol)
  float *detPtAlt;     ///< Altitude pointing for "
  float *detPtRa;      ///< Right Ascension pointing for "
  float *detPtDec;     ///< Declination pointing for "




  actData ramin;
  actData ramax;
  actData decmin;
  actData decmax;
  actData *avec;
  int have_avec;


  actData **time_constants;
  mbNoiseVectorStruct *noise;
  mbNoiseVectorStructBands *band_noise;
  mbNoiseStructBandsVecs *band_vecs_noise;


  actData *noise_amp;
  actData *noise_knee;
  actData *powlaw;
  actData *dra;
  actData *ddec;  
  
#ifdef ACTPOL
  //ACTpolArray *polarray;
  //ACTpolWeather weather;  //if there's TOD-based weather info.
  ACTpolPointingFit *actpol_pointing;
  actData *hwp;
  actData **twogamma_saved;
#endif

  PointingFit *pointing_fit;  //pointing fit, turn alt/az into ra/dec
  int **pixelization_saved;  //save a map pixelization in here.  Will break if there are multiple classes of maps with different pixelizations.
  actData **ra_saved;
  actData **dec_saved;
  actData **data_saved;
  BadTimestreams *bad_timestreams;

  bool have_plans;
  act_fftw_plan p_forward;
  act_fftw_plan p_back;
  char *froot;
  long seed;


} mbTOD;

// From tod.c

mbTOD *mbTODAlloc(int nrow, int ncol, int ndata, int ndet, const int *rows, const int *cols, actFilter band);
void mbTODRealloc(mbTOD *tod, const int ndata);

void   mbTODDarkAlloc(mbTOD *tod, const int *dets, int ndet);
void   mbTODDataAlloc( mbTOD *tod );
void   mbTODDataFree( mbTOD *tod );
void   mbTODDataClear( mbTOD *tod );
int    mbTODGetNumberOfDetectors( const mbTOD *tod );
void   mbTODGetDetectorRowCol( const mbTOD *tod, int idet, int *row, int *col );

// Routines from todUtils.c

actData **psAllocMatrix(int nrow,int ncol);  //psallocate a 2-d matrix;
psErrorCode mbTODSetPixel(mbTOD *tod,           // the TOD to set
                          const int row, const int col, // which pixel
                          actData val);                 // desired value
void mbTODClearAllData(mbTOD *tod);
//void *mbArrayFromTOD(int *n, const mbTOD *tod, const int camRow, const int camCol,
//                     psElemType type, bool copy, bool own);
//psVector *mbVectorFromTOD(const mbTOD *tod, const int row, const int col, psElemType type, bool copy);
//psErrorCode mbTODFromVector(mbTOD *tod, const int row, const int col, const psVector *vec);
psErrorCode mbCopyTODData(mbTOD *otod, int ocamRow, int ocamCol,                
                          const mbTOD *itod, int icamRow, int icamCol);

void mbTODSetDt(mbTOD *tod);            // the time-ordered-data
void mbTODSetRaDec(mbTOD *tod,          // the time-ordered-data
                   double longitude,    // longitude (radians)
                   double latitude);    // latitude (radians)

void mbTODSetLMST( mbTOD *tod, double longitude );


// I/O routines from todIO.c
#if 0
mbTOD *mbReadTOD(const char restrict *filename,  // desired dirfile, or NULL
                 int row0,                       // first row to read
                 int row1,                       // last row to read (if -ve, relative to end of file)
                 int npix,                       // number of camera pixels to be read (-1: full array)
                 const int *camrows,             // rows corresponding to camera pixels
                 const int *camcols,             // columns corresponding to camera pixels
                 actFilter band);

mbTOD *mbReadFitsTOD(actFits *restrict fits,         // file to read from
                     const char *restrict hdu,       // desired HDU, or NULL
                     int row0,                       // first row to read
                     int row1,                       // last row to read (if -ve, relative to end of file)
                     int npix,                       // number of camera pixels to be read (-1: full array)
                     const int *camrows,             // rows corresponding to camera pixels
                     const int *camcols,             // columns corresponding to camera pixels
                     actFilter band);

mbTOD *mbReadDataTOD(
  actFits *restrict fits,             //!< file to read from
  int row0,                           //!< first row to read
  int row1,                           //!< last row to read (if -ve, relative to end of file)
  int npix,                           //!< number of camera pixels to be read (-1: full array)
  const int *camrows,                       //!< rows corresponding to camera pixels
  const int *camcols,                       //!< columns corresponding to camera pixels
  actFilter band);                      //!< frequency band for data

int mbGetTODLength( 
        const char *restrict filename ); //!< dirfile we're inquiring about

#endif
// From tod.c

void mbTODShrink(  mbTOD *tod,           //!< the tod to operate on
                   int ndet,             //!< number of pixels for which to allocate memory
                   const int *rows,      //!< row numbers for pixels to be stored
                   const int *cols);     //!< column numbers for pixels to be stored

//void mbTODSetPointingOffset( mbTOD *tod,
//                             struct mbPointingOffset *po);
psErrorCode mbTODSetDetectorPointing( mbTOD *tod, int row, int col );



bool mbTODHasPointingOffset( const mbTOD *tod );

void mbTODRemoveMode( mbTOD *tod, const float *mode, int nmode, int row, int col );

void mbTODShiftData( mbTOD *tod, int shift );

void mbTODCalibrate( mbTOD *tod, float cal, int row, int col );

void mbDetrendDataC(mbTOD *tod, const int *dets, int ndet);
void mbRetrendDataC(mbTOD *tod, const int *dets, int ndet);

actData *calculate_glitch_filterC(actData t_glitch, actData t_smooth, actData dt, int n, int filt_type);
int *glitch_one_detector_simple(actData *mydat, int n,bool do_smooth, bool apply_glitch, bool do_cuts,
                                actData nsig, actData t_glitch,actData t_smooth, actData dt, int filt_type);
void mbGlitchC(mbTOD *tod, const float *filt, int nelem, const int *dets, int ndet, bool do_smooth,
               bool apply_glitch, mbCuts *cuts, actData nsig, int maxGlitch);
void mbFlattenC(mbTOD *tod, const float *filt, int nelem, const int *dets, int ndet, actData nsig);
int mbMaxSamp(const mbTOD *tod, int row, int col, bool maxmin);


void filter_whole_tod(mbTOD *tod, const actData *filt, actData **to_filt_in);
void smooth_tod(mbTOD *tod, actData t_smooth, actData **to_filt_in);
void highpass_tod(mbTOD *tod, actData nu_low, actData nu_high);
void mbTODFourierFilter(mbTOD *tod,  const float *real, const float *imag, int nelem, const float *tau, 
                        const int *dets, int ndet); 
// void mbArrayFourierFilter(float **data,  const float *real, const float *imag, int nelem, int nvec); 
// void mbArrayFourierFilter2(float *datos,  const float *real, const float *imag, int nelem); 
// void mbProcessVPhiJumps(mbTOD *tod, const int *dets, int ndet, int w1, int w2, int w3, int skip);


int *glitch_all_detectors_simple(mbTOD *tod, bool do_smooth, bool apply_glitch, bool do_cuts, actData nsig, actData t_glitch, actData t_smooth, int filt_type);
actData *CopyData1Det(mbTOD *tod,int det);
int ind_from_rowcol(mbTOD *tod, int row, int col);

#endif
