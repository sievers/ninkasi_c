#ifndef NINKASI_POINTING_H
#define NINKASI_POINTING_H
#include "config.h"
#include "ninkasi.h"
#include "ninkasi_mathutils.h"

typedef struct {
  actData *ra;
  actData *dec;
  actData *alt;
  actData *az;
  actData *tod_alt;
  actData *tod_az;

  actData *az_coarse;
  actData *alt_coarse;
  actData *ra_coarse;
  actData *dec_coarse;
  actData *time_coarse;
  PointingFit *pointing_fit; 
} PointingFitScratch;


void destroy_pointing_fit_raw(PointingFit *fit);
void destroy_pointing_fit(mbTOD *tod);
void assign_tod_ra_dec(mbTOD *tod);
void set_tod_pointing_tiled(mbTOD *tod, actData *azvec, int naz, actData *altvec, int nalt, actData **ra_mat, actData **dec_mat, actData *ra_clock, actData *dec_clock, int nclock);

void cut_mispointed_detectors(mbTOD *tod);
mbPointingOffset *nkReadPointingOffset(const char *filename);
void find_tod_radec_lims(mbTOD *tod);

void set_tod_starting_altaz_ctime(mbTOD *tod, actData alt, actData az, double ctime);


void get_radec_from_altaz_fit_tiled(const TiledPointingFit *tile, actData *alt, actData *az, actData *ctime, actData *ra, actData *dec, int npoint);
void get_radec_from_altaz_exact_1det(const mbTOD *tod,int det,    PointingFitScratch *scratch);
void get_radec_from_altaz_fit_1det(const mbTOD *tod,int det, PointingFitScratch *scratch);
void get_radec_from_altaz_fit_1det_coarse(const mbTOD *tod, int det, PointingFitScratch *scratch);
void get_radec_from_altaz_fit_1det_coarse_exact(const mbTOD *tod, int det, PointingFitScratch *scratch);
actData find_max_pointing_err(const mbTOD *tod);

void destroy_pointing_fit_scratch(PointingFitScratch *scratch);
PointingFitScratch *allocate_pointing_fit_scratch(const mbTOD *tod);
mbPointingOffset *nkPointingOffsetAlloc(int nrow,int ncol,int fitType);

void destroy_pointing_offset(mbPointingOffset *pt);
int *find_az_turnarounds(mbTOD *tod, int *nturn);
void find_pointing_pivots(mbTOD *tod, actData tol);
void shift_tod_pointing(mbTOD *tod, actData ra_shift, actData dec_shift);
int act_observed_altaz_to_mean_radec( const Site *site, double freq_GHz,
        int n, const double ctime[], const actData alt[], const actData az[],
				      actData ra[], actData dec[] );


#endif
