
#ifndef NINKASI_PROJECTION_H
#define NINKASI_PROJECTION_H




#if 1
#include "ninkasi_projection_types.h"
#else
typedef enum {
  NK_RECT,
  NK_CEA,
  NK_TAN
} nkProjectionType;


typedef struct {
  nkProjectionType proj_type;
  
  actData radelt;  //maps to cdelt1
  actData decdelt; //maps to cdelt2

  int rapix; //maps to crpix1
  int decpix; //maps to crpix2
  actData pv;  //maps to pv2_1

  
} nkProjection;

#endif
#include "ninkasi_types.h"

void get_map_projection(const mbTOD *tod, const MAP *map, int det, int *ind, PointingFitScratch *scratch);
void get_map_projection_wchecks(const mbTOD *tod, const MAP *map, int det, int *ind, PointingFitScratch *scratch, bool *inbounds);
int set_map_projection_cea_simple( MAP *map);
int set_map_projection_tan_simple( MAP *map);
int set_map_projection_cea_simple_keeppix( MAP *map);
void radec2pix_cea(MAP *map, actData ra, actData dec, int *rapix, int *decpix);
void pix2radec_cea(MAP *map, int rapix, int decpix, actData *ra, actData *dec);
int set_map_projection_cea_predef( MAP *map,actData radelt, actData decdelt, int rapix, actData decpix, actData pv, int nra, int ndec);

void pix2radec_tan(MAP *map, int rapix, int decpix, actData *ra, actData *dec);
void set_map_projection_tan_predef(MAP *map, actData ra_cent, actData dec_cent, actData rapix, actData decpix, actData pixsize, int nra, int ndec);
void radec2xy_tan_raw(actData *x, actData *y, actData ra, actData dec, actData ra0, actData dec0);
void radec2xy_tan(actData *x, actData *y,actData ra, actData dec, nkProjection *proj);



#endif
