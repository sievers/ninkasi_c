
#ifndef NINKASI_PROJECTION_TYPES_H
#define NINKASI_PROJECTION_TYPES_H

#include <ninkasi_defs.h>

typedef enum {
  NK_RECT,
  NK_CEA,
  NK_TAN

} nkProjectionType;


typedef struct {
  nkProjectionType proj_type;
  
  actData radelt;  //maps to cdelt1
  actData decdelt; //maps to cdelt2

  actData ra_cent;  //maps to crval1
  actData dec_cent;  //maps to crval2


  actData rapix; //maps to crpix1
  actData decpix; //maps to crpix2
  actData pv;  //maps to pv2_1

  
} nkProjection;

#endif
