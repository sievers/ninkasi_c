
#ifndef NINKASI_MATHUTILS_H
#define NINKASI_MATHUTILS_H



#include "config.h"
#include "ninkasi_defs.h"




typedef struct {
  actData xcent;  //will subtract off x/y median values for stability
  actData ycent;
  actData xwidth; //will scale by xwidth and ywidth, 
  actData ywidth;
  int nx;  //number of parameters in the x direction
  int ny;  //number of parameters in the y direction
  actData **params;
  
} PolyParams2d;

PolyParams2d *copy_2d_polyfit(PolyParams2d *fit);
int  invert_posdef_mat(actData **mat, int n);
int  invert_posdef_mat_double(double **mat, int n);
actData  compute_median(int ndata, actData *data);
actData compute_median_inplace(int ndata, actData *data);
actData compute_median_scat(int ndata, actData *data);
actData sselect(unsigned long k, unsigned long n, actData *arr);
PolyParams2d *fit_2d_poly(actData *x, actData *y, actData *d, int n, actData *errs, int nx_param, int ny_param);
void eval_2d_poly_inplace(actData *x, actData *y, int n, PolyParams2d *fit, actData *d);
void eval_2d_poly_pair_inplace(actData *x, actData *y, int n,const  PolyParams2d *fit, actData *d, const PolyParams2d *fit2, actData *d2);
actData *eval_2d_poly(actData *x, actData *y, int n, PolyParams2d *fit);

void add_outer_product(actData *col, int ncol, actData *row, int nrow, actData **mat);

actData vecmin(actData *vec, int n);
actData vecmax(actData *vec, int n);

float svecmin(float *vec, int n);
float svecmax(float *vec, int n);

int ivecmin(int *vec, int n);
int ivecmax(int *vec, int n);
void destroy_2d_poly(PolyParams2d *fit);

actData mygasdev(unsigned *seed);
actData myrand(unsigned *seed);

void act_gemm(char transa, char transb, int m, int n, int k, actData alpha, actData *a, int lda, actData *b, int ldb, actData beta, actData *c, int ldc);


/*--------------------------------------------------------------------------------*/
inline actData cos5(actData x) {
  //good to max err of ~5e-8 on (pi/2,pi/2)
  actData x2=x*x;
  return 0.999999953247608 +x2*(-0.499999050628101 +x2*( 0.0416635789306837 +x2*( -0.00138536669329216 +x2* 2.31531741531081e-05)));
}

/*--------------------------------------------------------------------------------*/
inline actData sin4(actData x)
{
  actData x2=x*x;
  actData x4=x2*x2;

  //Think this one is better, sticking w/ old for not-changingness sake.
  //return x*(0.999996601050172  -0.16664823561673*x2  + 0.00830628614181227*x4     -0.000183627485767557*x2*x4);
  return x*(0.999999241345692  -0.166656796188478*x2 + 0.00831322507990857*x4   -0.000185234483301391*x2*x4 );

}

/*--------------------------------------------------------------------------------*/
inline actData sin5( actData x)
{
  
  actData x2=x*x;
  actData x4=x2*x2;

  //Think this one is better, but sticking with old for not-changingness sake
  //return x*(0.999999976513754 -0.166666475934889*x2+ 0.00833289922282484*x4 -0.000198008653066048*x2*x4+ 2.59043003161813e-06*x4*x4); 
  return x*(0.999999995715839  -0.166666579699042*x2 + 0.00833305061731921*x4  -0.000198090463568892*x2*x4 + 2.6051662751399e-06*x4*x4);
  
}

/*--------------------------------------------------------------------------------*/



#endif
