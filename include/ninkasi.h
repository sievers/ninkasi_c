
#ifndef NINKASI_H
#define NINKASI_H


#include <stdio.h>
#ifndef NO_FFTW
#include <fftw3.h>
#endif
#include <sys/time.h>
#include <omp.h>
#include <ninkasi_defs.h>
#include <ps_stuff.h>

#include <stdarg.h>
#include <stdbool.h>

#include "parse_strings.h"
#include "astro.h"

#include "mbTOD.h"

#include "ninkasi_pointing.h"
#include "ninkasi_projection_types.h"

//#include <mbCuts.h>
#if 1
#include "ninkasi_types.h"
#else
/*--------------------------------------------------------------------------------*/
struct params_struct_s {
  char inname[MAXLEN],outname[MAXLEN],tempname[MAXLEN],rawname[MAXLEN];
  char datanames[MAXTOD][MAXLEN];
  int ntod;
  char pointing_file[MAXLEN];
  char altaz_file[MAXLEN];
  actData tol;
  bool do_sim;
  bool do_blank;
  int maxiter;
  int maxtod;
  bool remove_common;
  bool remove_mean;  //use this to remove the common mode from the data, but not
                    //during mapmaking.
  bool no_noise;
  long seed;
  bool add_noise;
  actData pixsize;  //map pixel size, arcmin
  bool quit;
  bool rawonly;
  bool precondition;
  bool use_input_limits;
  bool deglitch;


  bool write_pointing;

  
  int n_use_rows;
  int n_use_cols;
  int *use_rows;  //actually use these rows/columns in mapping.  separate from 
  int *use_cols;  //cuts in that non-cut detectors not in here will go into things like common mode, noise stats etc.

} params;
typedef struct params_struct_s PARAMS;
/*--------------------------------------------------------------------------------*/

struct todvec_struct_s {
  int ntod; //# of privately owned tod's
  int total_tod;  //total number of tod's
  mbTOD *tods;
  char **my_fnames;
  //char *froot;
  char froot[MAXLEN];
  
  actData ramin,ramax,decmin,decmax;  //global ra/dec limits

} todvec_struct;
typedef struct todvec_struct_s TODvec;
/*--------------------------------------------------------------------------------*/

struct map_struct_s {
  actData pixsize;
  actData ramin,ramax,decmin,decmax;
  actData *map;
  int nx,ny;
  long npix;
#ifdef ACTPOL
  int pol_state;  //flag to tell me how many polarizations I have.  0/1=I, 2=Q+U, 3=I+Q+U
#endif

  //Stuff for reductions
  int have_locks;
  omp_lock_t  *locks;  //omp_init_lock   
  int lock_len;
  int nlock;

  nkProjection *projection;
} map_struct;
typedef struct map_struct_s MAP;

/*--------------------------------------------------------------------------------*/

struct mapvec_struct_s {
  int nmap;
  MAP **maps;
} mapvec_struct;
typedef struct mapvec_struct_s MAPvec;
/*--------------------------------------------------------------------------------*/

#endif



#include "ninkasi_projection.h"









actData mbElapsedTime(pca_time *tt);
void mbStartTime(pca_time *tt);
size_t freadwrite(void *ptr, size_t sz, long nobj, FILE *stream, int dowrite);

act_fftw_complex *act_fftw_malloc(size_t n);
void act_fftw_destroy_plan(act_fftw_plan p);
void act_fftw_free(act_fftw_complex *vec);
void act_fftw_execute_dft_r2c(act_fftw_plan p, actData *r, act_fftw_complex *c);
void act_fftw_execute_dft_c2r(act_fftw_plan p,  act_fftw_complex *c, actData *r);
act_fftw_plan act_fftw_plan_dft_r2c_1d(int n, actData *vec, act_fftw_complex *vec2,unsigned flags);
act_fftw_plan act_fftw_plan_dft_c2r_1d(int n, act_fftw_complex *vec2,actData *vec, unsigned flags);


void createFFTWplans1TOD(mbTOD *mytod);
void copy_mapset2mapset(MAPvec *map2, MAPvec *map);
actData mapset_times_mapset(MAPvec *x, MAPvec *y);
void remove_common_mode(mbTOD *tod);
void readwrite_simple_map(MAP *map, char *filename, int dowrite);
void detrend_data(mbTOD *tod);
void demean_data(mbTOD *tod);
void cut_tod_ends(mbTOD *tod,actData tcut);
void set_tod_window(mbTOD *tod,actData tcut);

void window_data(mbTOD *tod);
void unwindow_data(mbTOD *tod);
void apply_preconditioner( MAPvec *maps,MAPvec *weights,PARAMS *params);
void dump_data(mbTOD *tod, char *fname);
void destroy_map(MAP *map);
FILE *fopen_safe(char *filename, char *mode);
void get_pointing_vec(const mbTOD *tod, const MAP *map, int det, int *ind);
void get_pointing_vec_new(const mbTOD *tod, const MAP *map, int det, int *ind, PointingFitScratch *scratch);
void write_tod_pointing_to_disk(const mbTOD *tod, char *fname);
void mapset2det(const MAPvec *maps, mbTOD *tod, const PARAMS *params, actData *vec, int *ind, int det, PointingFitScratch *scratch);
void map2tod(const MAP *map, mbTOD *tod,const PARAMS *params);
void polmap2tod(MAP *map, mbTOD *tod);

void clear_map(MAP *map);
void clear_tod(mbTOD *tod);
actData **matrix(long n,long m);
int **imatrix(long n,long m);
actComplex **cmatrix(long n,long m);
float **smatrix(long n,long m);
void free_matrix(actData **mat);
void free_smatrix(float **mat);
double **dmatrix(long n,long m);
actData *vector(long n);
int *ivector(long n);
float *svector(long n);
double *dvector(long n);
actComplex *cvector(long n);
int how_many_tods(char *froot, PARAMS *params);
int find_my_tods(TODvec *tods, PARAMS *params);
int read_all_tod_headers(TODvec *tods,PARAMS *params);
void set_global_radec_lims(TODvec *tods);
actData tocksilent(pca_time *tt);
void tick(pca_time *tt);
void pca_pause(actData pauselen);
void mprintf(FILE *stream, char *format, ...);


void keep_1_det(mbTOD *tod,int row, int col);

void write_timestream(actData *vec, int n, char *fname) ;
bool is_det_used(mbTOD *tod, int det);


int get_parameters(int argc, char *argv[], PARAMS *params);
void print_options(PARAMS *params);
int setup_maps(MAPvec *maps, PARAMS *params);
void clear_mapset(MAPvec *maps);
void createFFTWplans(TODvec *tod);
void run_PCG(MAPvec *maps, TODvec *tods, PARAMS *params);

void allocate_tod_storage(mbTOD *tod);
void tod2mapset(MAPvec *maps, mbTOD *tod, PARAMS *params);
void assign_tod_value(mbTOD *tod, actData val);
void free_tod_storage(mbTOD *tod);
void map2det_scaled(const MAP *map, const mbTOD *tod, actData *vec, actData scale_fac, int *ind, int det, PointingFitScratch *scratch);

int get_map_poltag(const MAP *map);
int get_npol_in_map(const MAP *map);
void set_map_polstate(MAP *map, int *pol_state);
int is_map_polarized(MAP *map);


void tod2map(MAP *map, mbTOD *tod, PARAMS *params);
void tod2polmap(MAP *map,mbTOD *tod);
void tod2polmap_copy(MAP *map,mbTOD *tod);
int *tod2map_actpol(MAP *map, mbTOD *tod, int *ipiv_proc);


actData tod_times_map(const MAP *map, const mbTOD *tod, PARAMS *params);

int read_tod_data(mbTOD *tod);
MAP *make_map_copy(MAP *map);
MAP *deres_map(MAP *map);
MAP *upres_map(MAP *map);

actData map_times_map(MAP *x, MAP *y);
void map_axpy(MAP *y, MAP *x, actData a);
bool is_det_listed(const mbTOD *tod, const PARAMS *params, int det);
void purge_cut_detectors(mbTOD *tod);
mbUncut ***get_uncut_regions(mbTOD *tod);
int get_numel_cut(mbTOD *tod);
void set_global_indexed_cuts(mbTOD *tod);
mbCutFitParams ***setup_cut_fit_params(mbTOD *tod, int *nparams_from_length);
void setup_cutsfits_precon(mbTOD *tod);
void apply_cutfits_precon(mbTOD *tod, actData *params_in, actData *params_out);


mbUncut ***get_cut_regions(mbTOD *tod);
void get_tod_cut_regions(mbTOD *tod);
mbUncut ***get_cut_regions_global_index(mbTOD *tod);
int tod2cutvec(mbTOD *tod, actData *cutvec);
int cutvec2tod(mbTOD *tod, actData *cutvec);



void get_tod_uncut_regions(mbTOD *tod);
void reverse_tod_uncut_regions(mbTOD *tod);
void decimate_uncut_regions(mbTOD *tod);
void fill_gaps_stupid(mbTOD *tod);
void clear_cut_data(mbTOD *tod);
void reverse_tod_data(mbTOD *tod);

void get_data_corrs(mbTOD *tod);
void rotate_data(mbTOD *tod, char trans, actData *rotmat);
void multiply_all_data(mbTOD *tod,actData val);

long is_tod_inbounds(const MAP *map, mbTOD *tod,const PARAMS *params);
void *malloc_retry(size_t n);

actData how_far_am_i_from_radec_radians(actData ra, actData dec, mbTOD *tod);
actData how_far_am_i_from_radec(actData rah,actData ram,actData ras,actData dd,actData dm,actData ds,mbTOD *tod);
actData  *how_far_am_i_from_radec_radians_vec( actData *ra, actData *dec, int npt, mbTOD *tod);
int tod_hits_source(actData ra, actData dec, actData dist, mbTOD *tod);
void add_src2tod(mbTOD *tod, actData ra, actData dec, actData src_amp, const actData *beam, actData dtheta, int nbeam, int oversamp);
void add_srcvec2tod(mbTOD *tod, actData *ra, actData *dec, actData *src_amp, int nsrc,const actData *beam, actData dtheta, int nbeam, int oversamp);
void tod2srcvec(actData *src_amp_out,mbTOD *tod, actData *ra_in, actData *dec_in, int nsrc_in,const actData *beam, actData dtheta, int nbeam, int oversamp);
void find_map_index_limits(MAP *map, mbTOD *tod, int *imin_out, int *imax_out);
void invert_pol_precon(MAP *map);
void apply_pol_precon(MAP *map, MAP *precon);

//from ninkasi_projection
//void get_map_projection(const mbTOD *tod, const MAP *map, int det, int *ind, PointingFitScratch *scratch);

#endif
