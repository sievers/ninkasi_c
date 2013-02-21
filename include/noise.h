#if !defined(MB_NOISE_H)
#define MB_NOISE_H 1

#define MB_BAD_LIKE -1e30
#define MB_DEFAULT_MAX_ITER_NOISEFIT 100

#define MB_READ_NOISE 0
#define MB_WRITE_NOISE 1
#include "noise_types.h"

void act_syrk(char uplo, char trans, int n, int m, actData alpha, actData *a, int lda, actData beta, actData *b, int ldb);


int nkFitNoise1Det(NoiseParams1Pix *noise,mbTOD *tod,int mydet);
int nkFitNoise_LinearPowlaw(NoiseParams1Pix *noise,mbTOD *tod,int mydet);
void SubtractMedian(actData *vec, int n);
psErrorCode nkSetNoiseFreqRange(mbTOD *tod, NoiseParams1Pix *noise_params, int *nn_start_out, int *nn_max_out, actData *delta_out);
actData mbDot(long n,actData *x,long incx,actData *y,long incy);
psErrorCode nkCalculateOneoverFVecs(mbTOD *tod, NoiseParams1Pix *noise_params,actData ***vecs_out,int *nparam_out);
psErrorCode nkFitSpecUncorrDataQuadratic(actData *params,actData *data,long n,actData **vecs,int nparam);
psErrorCode nkFitStartingSpecUncorrData(actData *params,actData *data,int  n,actData **vecs,int nparam);
mbNoiseVectorStruct *nkFitTODNoise(mbTOD *tod, mbNoiseType  noise_type, actData minFreq, actData maxFreq, actData powlaw);
void mbVecMult(long n,actData *x,actData *y,actData *z);
actData mbVecSum(long n,actData *x);
psErrorCode mbGetCurveDerivUncorrData(actData *params, actData *data, actData *noise, actData **vecs, actData *like, actData *deriv, actData **curve, long n, int nparam);
psErrorCode mbGetShiftsUncorrData(int nparam, actData *shifts, actData **curve,actData *deriv,actData like,actData oldlike,actData *lambda,actData *max_shift);
void mbDaxpy(long n,actData a,actData *x,long incx,actData *y,long incy);
psErrorCode nkFitNoiseUncorrData(actData *params, actData *data, actData *noise, actData **vecs, long n,  int nparam,actData *like );
void cut_badnoise_dets(mbTOD *tod);
void filter_data_wnoise(mbTOD *tod);
void nkAssignNoiseKnee( mbTOD *tod, actData knee);
void add_noise_to_tod_new(mbTOD *tod);
void set_tod_noise(mbTOD *tod, actData white, actData knee, actData powlaw);

void apply_real_filter_to_data(mbTOD *tod, const actData *filt);
void apply_complex_filter_to_data(mbTOD *tod, actComplex *filt);
void nkDeButterworth(mbTOD *tod);
void nkReButterworth(mbTOD *tod);
void nkDeconvolveTimeConstants(mbTOD *tod);
void nkReconvolveTimeConstants(mbTOD *tod);
void destroy_tod_noise(mbTOD *tod);


int get_nn(int n);
actComplex **fft_all_data(mbTOD *tod);
void ifft_all_data(mbTOD *tod,actComplex **data_fft) ;
void ifft_all_data_flag(mbTOD *tod,actComplex **data_fft,unsigned flag);
actComplex **fft_all_data_flag(mbTOD *tod, unsigned flags);
actData **get_banded_correlation_matrix_from_fft(mbTOD *tod, actComplex **data_fft, actData nu_min, actData nu_max);
void act_syrk(char uplo, char trans, int n, int m, actData alpha, actData *a, int lda, actData beta, actData *b, int ldb);
void get_eigenvectors(actData **mat, int n);
void allocate_tod_noise_bands(mbTOD *tod,actData *bands, int nband);
actComplex **apply_banded_rotations(mbTOD *tod, actComplex **mat_in, bool do_forward);
void get_simple_banded_noise_model(mbTOD *tod, bool *do_rots, mbNoiseType *types);
void get_simple_banded_noise_model_onerotmat(mbTOD *tod, bool *do_rots, mbNoiseType *types);
int fit_banded_noise_1det_constant(mbNoiseParams1PixBand *params, actData *dat);
int fit_banded_noise_1det( mbNoiseParams1PixBand *params, actComplex *dataft);
void apply_banded_noise_model(mbTOD *tod);
void apply_banded_noise_complex(mbTOD *tod, actComplex **dat);
void scale_banded_noise_band( mbTOD *tod,int which_band,actData fac);
actData *get_freq_vec(mbTOD *tod);
bool do_I_have_rotations(mbTOD *tod);
void add_noise_to_tod_gaussian(mbTOD *tod);


void apply_noise_1det(mbTOD *tod, int det, actComplex *ts);
void apply_noise_1det_powlaw(mbTOD *tod, int det, act_fftw_complex *ts );
void apply_noise(mbTOD *tod);
void set_noise_powlaw(mbTOD *tod, actData *amps, actData *knees, actData *pows);

void simple_test_diag_proj_noise_inv(actData **data_in, actData **data_out, actData *noise, actData **vecs, int ndata, int ndet, int nvecs);
void apply_diag_proj_noise_inv_bands(actData **data_in, actData **data_out, actData *noise, actData **vecs, int ndata, int ndet, int nvecs, int imin, int imax);

void fill_sin_cos_mat(actData *theta, int ndata, int nterm, actData **mat);
void fill_tod_sin_cos_vec(mbTOD *tod, int nterm, actData *vec);
void fit_hwp_poly_to_data(mbTOD *tod, int nsin, int npoly, actData **fitp, actData **vecs);
void remove_hwp_poly_from_data(mbTOD *tod, int nsin, int npoly);
int get_demodulated_hwp_data(mbTOD *tod, actData hwp_freq, actComplex **tdata,actComplex **poldata);
int remodulate_hwp_data(mbTOD *tod, actData hwp_freq, actComplex **tdata,actComplex **poldata);
actData get_hwp_freq(mbTOD *tod);
void demodulate_data(mbTOD *tod, DemodData *demod);
void destroy_demod_data(DemodData *demod);
void free_demod_data(DemodData *demod);
int get_demod_nchannel(DemodData *demod);
void set_demod_freqs(DemodData *demod, actData *freqs, int nfreq);
DemodData *init_demod_data(mbTOD *tod, actData hwp_freq, actData lowpass_freq, actData lowpass_taper, actData highpass_freq,actData highpass_taper);


#endif
