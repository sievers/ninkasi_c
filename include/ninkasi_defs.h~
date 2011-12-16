#ifndef NINKASI_DEFS_H
#define  NINKASI_DEFS_H

#include "config.h"

#include <fftw3.h>
#include <complex.h>

#define MAXLEN 256
#define DOWRITE 1
#define DOREAD 0
#define EPS 1e-4
//#define MAPS_PREALLOC
#define MAXDET 2000   //these two guys are only for setting seeds consistently if
#define MAXTOD 10000   //simulating noise internally.   Take that back - also now statically store space for tod filenames.


#define ACT_NO_VALUE -98747423

#ifndef INT_MAX
#define INT_MAX 2147483647
#endif


//sometimes this doesn't show up
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



#ifndef  RAD2DEG
#define RAD2DEG 57.2957795130823
#endif


typedef struct timeval pca_time;

typedef pca_time mbTimeValue;


#ifdef ACTDATA_DOUBLE
  typedef double actData;
  typedef fftw_plan act_fftw_plan;
  typedef fftw_complex act_fftw_complex;

#else
#define SINGLE
  typedef float actData;
  typedef fftwf_plan act_fftw_plan;
  typedef fftwf_complex act_fftw_complex;
//typedef float complex actComplex;
//typedef act_fftw_complex actComplex;
#endif

typedef act_fftw_complex actComplex;


#if 0

#ifdef ACTDATA_DOUBLE
  typedef double actData;
  typedef fftw_plan act_fftw_plan;
  typedef fftw_complex act_fftw_complex;
  typedef double complex actComplex;
#else
#define SINGLE
  typedef float actData;
  typedef fftwf_plan act_fftw_plan;
  typedef fftwf_complex act_fftw_complex;
  typedef float complex actComplex;
#endif
#endif

#endif
