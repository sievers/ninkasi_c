#ifndef NINKASI_TYPES_H
#define NINKASI_TYPES_H

#include <ninkasi_defs.h>

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


#define MAX_NPOL 6
#define POL_ERROR 0
#define POL_I 1
#define POL_IQU 2
#define POL_QU 3
#define POL_IQU_PRECON 4
#define POL_QU_PRECON 5



struct map_struct_s {
  actData pixsize;
  actData ramin,ramax,decmin,decmax;
  actData *map;
  int nx,ny;
  long npix;
#ifdef ACTPOL
  //int pol_state;  //flag to tell me how many polarizations I have.  0/1=I, 2=Q+U, 3=I+Q+U
  int pol_state[MAX_NPOL];  //flag to tell me how many polarizations I have.  0/1=I, 2=Q+U, 3=I+Q+U
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


#endif
