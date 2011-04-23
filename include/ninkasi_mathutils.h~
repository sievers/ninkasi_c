
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

#endif
