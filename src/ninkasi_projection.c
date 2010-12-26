//Module to turn RA/Dec into map coordinates.
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "ninkasi.h"
#include "ninkasi_projection.h"

/*--------------------------------------------------------------------------------*/
static inline actData sin5(actData x)
{
  
  actData x2=x*x;
  actData x4=x2*x2;
  
  return x*(0.999999995715839  -0.166666579699042*x2 + 0.00833305061731921*x4  -0.000198090463568892*x2*x4 + 2.6051662751399e-06*x4*x4);
  
}

/*--------------------------------------------------------------------------------*/
static inline actData sin4(actData x)
{

  actData x2=x*x;
  actData x4=x2*x2;

  return x*(0.999999241345692  -0.166656796188478*x2 + 0.00831322507990857*x4   -0.000185234483301391*x2*x4 );

}


/*--------------------------------------------------------------------------------*/

void get_map_projection(const mbTOD *tod, const MAP *map, int det, int *ind, PointingFitScratch *scratch)
{
#if 1
  get_map_projection_wchecks(tod,map,det,ind,scratch,NULL);
#else
  get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
  switch(map->projection->proj_type) {
  case(NK_RECT): 
    //printf("Doing rectangular projection.\n");
    for (int i=0;i<tod->ndata;i++)
      ind[i]=(int)((scratch->dec[i]-map->decmin)/map->pixsize)+map->ny*(int)((scratch->ra[i]-map->ramin)/map->pixsize);    
    break;
  case (NK_CEA):
    {
      //printf("Doing CEA projection.\n");
      double rafac=RAD2DEG/map->projection->radelt;  //180/pi since ra is in radians, but CEA FITS likes degrees.
      double decfac=RAD2DEG/map->projection->pv/map->projection->decdelt;
      for (int i=0;i<tod->ndata;i++) {
	int rapix=(scratch->ra[i]*rafac+map->projection->rapix-1.0+0.5);  //-1 to go from zero-offset to unity-offset maps
	int decpix=(sin(scratch->dec[i])*decfac+map->projection->decpix-1.0+0.5);  //+0.5 so pixel is centered on coordinates
	ind[i]=map->nx*decpix+rapix;
      }
    }
    break;
  default: 
    printf("unknown type in get_map_projection.\n");
    assert(1==0);  //fail if we're not set up
    break;

  }
#endif
}
/*--------------------------------------------------------------------------------*/

void get_map_projection_wchecks(const mbTOD *tod, const MAP *map, int det, int *ind, PointingFitScratch *scratch, bool *inbounds)
//if inbounds is non-null, check to see if any pixels are out of bounds.  If so, flag 'em.  The check is done once per detector, so
//this should be as fast as a non-checking version for mainline production.
{

  get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
  switch(map->projection->proj_type) {
  case(NK_RECT): 
    //printf("Doing rectangular projection.\n");
    for (int i=0;i<tod->ndata;i++)
      ind[i]=(int)((scratch->dec[i]-map->decmin)/map->pixsize)+map->ny*(int)((scratch->ra[i]-map->ramin)/map->pixsize);    
    break;
  case(NK_TAN): {
    
    actData x,y;
    for (int i=0;i<tod->ndata;i++) {
      radec2xy_tan(&x,&y,scratch->ra[i],scratch->dec[i],map->projection);
      ind[i]=(int)(x+0.5)+map->nx*((int)(y+0.5));
    }
    break;
  }
    
  case (NK_CEA):
    {
      //printf("Doing CEA projection.\n");
      double rafac=RAD2DEG/map->projection->radelt;  //180/pi since ra is in radians, but CEA FITS likes degrees.
      double decfac=RAD2DEG/map->projection->pv/map->projection->decdelt;
      for (int i=0;i<tod->ndata;i++) {
	int rapix=scratch->ra[i]*rafac+map->projection->rapix-1+0.5;
	int decpix=sin5(scratch->dec[i])*decfac+map->projection->decpix-1+0.5;  //change!  13 Aug 2010, should be faster, good to 1e-3 arcsec
	//int decpix=sin(scratch->dec[i])*decfac+map->projection->decpix-1+0.5;  //change!  13 Aug 2010, should be faster, good to 1e-3 arcsec
	ind[i]=map->nx*decpix+rapix;
      }
      if (inbounds)
	{
	  printf("doing inbounds\n");
	  double rafac=RAD2DEG/map->projection->radelt;  //180/pi since ra is in radians, but CEA FITS likes degrees.
	  double decfac=RAD2DEG/map->projection->pv/map->projection->decdelt;
	  for (int i=0;i<tod->ndata;i++) {
	    int rapix=scratch->ra[i]*rafac+map->projection->rapix-1;
	    int decpix=sin(scratch->dec[i])*decfac+map->projection->decpix-1;
	    inbounds[i]=true;
	    if (rapix<0)
	      inbounds[i]=false;
	    if (decpix<0)
	      inbounds[i]=false;
	    if (rapix>=map->nx)
	      inbounds[i]=false;
	    if (decpix>=map->ny)
	      inbounds[i]=false;
	    
	  }
	}
    }
    break;
  default: 
    printf("unknown type in get_map_projection.\n");
    assert(1==0);  //fail if we're not set up
    break;

  }
}

/*--------------------------------------------------------------------------------*/
void radec2pix_cea(MAP *map, actData ra, actData dec, int *rapix, int *decpix)
{
  nkProjection *projection=map->projection;
  double rafac=RAD2DEG/map->projection->radelt;
  double decfac=RAD2DEG/map->projection->pv/map->projection->decdelt;
  *rapix=ra*rafac+projection->rapix-1;
  *decpix=sin(dec)*decfac+projection->decpix-1;
  
}

/*--------------------------------------------------------------------------------*/
void pix2radec_cea(MAP *map, int rapix, int decpix, actData *ra, actData *dec)
{
  nkProjection *projection=map->projection;
  double rafac=RAD2DEG/map->projection->radelt;
  double decfac=RAD2DEG/map->projection->pv/map->projection->decdelt;
  *ra=(rapix-projection->rapix+1)/rafac;
  *dec=asin((decpix-projection->decpix+1)/decfac);
  
}
/*--------------------------------------------------------------------------------*/
void pix2radec_tan(MAP *map, int rapix, int decpix, actData *ra, actData *dec)
{

  printf("ra crap is %14.4e %14.4e %14.4e\n",rapix,map->projection->rapix,map->projection->radelt);

  double xx=(rapix-map->projection->rapix)*map->projection->radelt;
  double yy=(decpix-map->projection->decpix)*map->projection->decdelt;
  
  printf("xx and yy are %14.4e %14.4e\n",xx,yy);

  double rho=sqrt(xx*xx+yy*yy);
  double cc=atan(rho);

  printf("rho and cc-1 are %14.4e %14.4e\n",rho,cc-1);

  *dec=asin(cos(cc)*sin(map->projection->dec_cent)+yy*sin(cc)*cos(map->projection->dec_cent)/rho);
  *ra=map->projection->ra_cent+atan(xx*sin(cc)/(rho*cos(map->projection->dec_cent)*cos(cc)-yy*sin(map->projection->dec_cent)*sin(cc)));
    
}
/*--------------------------------------------------------------------------------*/

int set_map_projection_cea_simple( MAP *map)
{
  map->projection->proj_type=NK_CEA;

  double cos0=cos(0.5*(map->decmax+map->decmin));
  

  map->projection->decdelt=map->pixsize/cos0*RAD2DEG;
  map->projection->radelt= -map->projection->decdelt;
  map->projection->pv=cos0*cos0;
  actData ra0=map->ramax*RAD2DEG/map->projection->radelt;
  actData ra1=map->ramin*RAD2DEG/map->projection->radelt+1;  //+1 is in in case of equality
  actData dec0=sin(map->decmin)*RAD2DEG/map->projection->pv/map->projection->decdelt;
  actData dec1=sin(map->decmax)*RAD2DEG/map->projection->pv/map->projection->decdelt+1;  //+1 is in case of equality
  map->projection->rapix=-ra0;
  map->projection->decpix=-dec0;
  map->nx=ra1-ra0;
  map->ny=dec1-dec0;
  free(map->map);
  map->npix=map->nx*map->ny;
  map->map=(actData *)malloc(sizeof(actData)*map->npix);

  mprintf(stdout,"Offsets are %d %d\n",map->projection->rapix,map->projection->decpix);
  mprintf(stdout,"Pixsizes are %14.6f %14.6f\n",map->projection->radelt,map->projection->decdelt);
  mprintf(stdout,"ra/dec0/1 are %14.8f %14.8f %14.8f %14.8f\n",ra0,ra1,dec0,dec1);
  mprintf(stdout,"pv is %14.6f\n",map->projection->pv);
  mprintf(stdout,"map sizes are %d %d\n",map->nx, map->ny);

  return 0;
}


/*--------------------------------------------------------------------------------*/
void radec2xy_tan_raw(actData *x, actData *y, actData ra, actData dec, actData ra0, actData dec0)
{

  
  actData cosc=sin(dec0)*sin(dec)+cos(dec0)*cos(dec)*cos(ra-ra0);
  *x=cos(dec)*sin(ra-ra0)/cosc;
  *y=(cos(dec0)*sin(dec)-sin(dec0)*cos(dec)*cos(ra-ra0))/cosc;

}
/*--------------------------------------------------------------------------------*/
void radec2xy_tan(actData *x, actData *y,actData ra, actData dec, nkProjection *proj)
{
  //printf("hello!\n");
  radec2xy_tan_raw(x,y,ra,dec,proj->ra_cent,proj->dec_cent);
  //printf("in here, x and y are %14.4e %14.4e\n",*x,*y);
  *x=proj->rapix+*x/proj->radelt;
  *y=proj->decpix-*y/proj->radelt;
}
/*--------------------------------------------------------------------------------*/
void set_map_projection_tan_explicit(MAP *map, actData rapix, actData decpix, actData radelt, actData decdelt, actData pv, actData ra_cent, actData dec_cent, int nra, int ndec)
{
  if (map==NULL) {
    printf("appear to be missing map.\n");
    return;
  }
  if (map->projection==NULL) {
    printf("projection appears to be missing.\n");
    return;
  }
  map->projection->proj_type=NK_TAN;
  map->projection->pv=pv;
  map->projection->radelt=radelt;
  map->projection->decdelt=decdelt;
  map->projection->ra_cent=ra_cent;
  map->projection->dec_cent=dec_cent;
  map->projection->rapix=rapix;
  map->projection->decpix=decpix;
  map->nx=nra;
  map->ny=ndec;
  map->npix=nra*ndec;
  printf("doing minima.\n");
  pix2radec_cea(map,0,0,&(map->ramin),&(map->decmin));
  pix2radec_cea(map,map->nx-1,map->ny-1,&(map->ramax),&(map->decmax));
  
  if (map->ramin>map->ramax) {
    actData tmp=map->ramin;
    map->ramin=map->ramax;
    map->ramax=tmp;
  }
  
  if (map->decmin>map->decmax) {
    actData tmp=map->decmin;
    map->decmin=map->decmax;
    map->decmax=tmp;
  }

  free(map->map);
  map->map=(actData *)malloc(sizeof(actData)*map->npix);


}
/*--------------------------------------------------------------------------------*/
void set_map_projection_tan_predef(MAP *map, actData ra_cent, actData dec_cent, actData rapix, actData decpix, actData pixsize, int nra, int ndec)
{
  map->projection->proj_type=NK_TAN;
  map->projection->pv=0;
  map->projection->radelt=-pixsize;
  map->projection->decdelt=pixsize;
  map->projection->ra_cent=ra_cent;
  map->projection->dec_cent=dec_cent;
  map->projection->rapix=rapix;
  map->projection->decpix=decpix;



  printf("setting rapix/decpix to %14.4g %14.4g\n",map->projection->rapix,map->projection->decpix);
  printf("setting rapix/decpix from %14.4g %14.4g\n",rapix,decpix);

  map->nx=nra;
  map->ny=ndec;
  map->npix=nra*ndec;

  pix2radec_cea(map,0,0,&(map->ramin),&(map->decmin));
  pix2radec_cea(map,map->nx-1,map->ny-1,&(map->ramax),&(map->decmax));

  
  if (map->ramin>map->ramax) {
    actData tmp=map->ramin;
    map->ramin=map->ramax;
    map->ramax=tmp;
  }
  
  if (map->decmin>map->decmax) {
    actData tmp=map->decmin;
    map->decmin=map->decmax;
    map->decmax=tmp;
  }



  free(map->map);
  map->map=(actData *)malloc(sizeof(actData)*map->npix);

  
}
/*--------------------------------------------------------------------------------*/
int set_map_projection_tan_simple( MAP *map)
{
  map->projection->proj_type=NK_TAN;

  double dec_cent=0.5*(map->decmax+map->decmin);
  double ra_cent=0.5*(map->ramax+map->ramin);

  actData x0,x1,y0,y1;
  radec2xy_tan_raw(&x0,&y0,map->ramin,map->decmin,ra_cent,dec_cent);
  radec2xy_tan_raw(&x1,&y1,map->ramax,map->decmax,ra_cent,dec_cent);
  
  printf("ra lims are %12.4e %12.4e, dec lims are %12.4e %12.4e\n",x0,x1,y0,y1);

  return 0;

  map->projection->decdelt=map->pixsize*RAD2DEG;
  map->projection->radelt= -map->projection->decdelt;
  map->projection->pv=1.0;

  actData ra0=map->ramax*RAD2DEG/map->projection->radelt;
  actData ra1=map->ramin*RAD2DEG/map->projection->radelt+1;  //+1 is in in case of equality
  actData dec0=sin(map->decmin)*RAD2DEG/map->projection->pv/map->projection->decdelt;
  actData dec1=sin(map->decmax)*RAD2DEG/map->projection->pv/map->projection->decdelt+1;  //+1 is in case of equality
  map->projection->rapix=-ra0;
  map->projection->decpix=-dec0;
  map->nx=ra1-ra0;
  map->ny=dec1-dec0;
  free(map->map);
  map->npix=map->nx*map->ny;
  map->map=(actData *)malloc(sizeof(actData)*map->npix);

  mprintf(stdout,"Offsets are %d %d\n",map->projection->rapix,map->projection->decpix);
  mprintf(stdout,"Pixsizes are %14.6f %14.6f\n",map->projection->radelt,map->projection->decdelt);
  mprintf(stdout,"ra/dec0/1 are %14.8f %14.8f %14.8f %14.8f\n",ra0,ra1,dec0,dec1);
  mprintf(stdout,"pv is %14.6f\n",map->projection->pv);
  mprintf(stdout,"map sizes are %d %d\n",map->nx, map->ny);

  return 0;
}

/*--------------------------------------------------------------------------------*/
int set_map_projection_cea_simple_keeppix( MAP *map)
//fix the above so that the pixel size is preserved
{
  map->projection->proj_type=NK_CEA;

  double cos0=cos(0.5*(map->decmax+map->decmin));
  
  //only change in this line
  map->projection->decdelt=map->pixsize*RAD2DEG;
  map->projection->radelt= -map->projection->decdelt;
  map->projection->pv=cos0*cos0;
  actData ra0=map->ramax*RAD2DEG/map->projection->radelt;
  actData ra1=map->ramin*RAD2DEG/map->projection->radelt+1;  //+1 is in in case of equality
  actData dec0=sin(map->decmin)*RAD2DEG/map->projection->pv/map->projection->decdelt;
  actData dec1=sin(map->decmax)*RAD2DEG/map->projection->pv/map->projection->decdelt+1;  //+1 is in case of equality
  map->projection->rapix=-ra0;
  map->projection->decpix=-dec0;
  map->nx=ra1-ra0;
  map->ny=dec1-dec0;
  free(map->map);
  map->npix=map->nx*map->ny;
  map->map=(actData *)malloc(sizeof(actData)*map->npix);

  mprintf(stdout,"Offsets are %d %d\n",map->projection->rapix,map->projection->decpix);
  mprintf(stdout,"Pixsizes are %14.6f %14.6f\n",map->projection->radelt,map->projection->decdelt);
  mprintf(stdout,"ra/dec0/1 are %14.8f %14.8f %14.8f %14.8f\n",ra0,ra1,dec0,dec1);
  mprintf(stdout,"pv is %14.6f\n",map->projection->pv);
  mprintf(stdout,"map sizes are %d %d\n",map->nx, map->ny);

  return 0;
}


/*--------------------------------------------------------------------------------*/
int set_map_projection_cea_predef( MAP *map,actData radelt, actData decdelt, int rapix, actData decpix, actData pv, int nra, int ndec)
{
  map->projection->proj_type=NK_CEA;
  
  map->projection->decdelt=decdelt;
  map->projection->radelt=radelt;
  map->projection->pv=pv;
  map->projection->rapix=rapix;
  map->projection->decpix=decpix;
  
  map->nx=nra;
  map->ny=ndec;
  map->npix=nra*ndec;
  
  pix2radec_cea(map,0,0,&(map->ramin),&(map->decmin));
  pix2radec_cea(map,map->nx-1,map->ny-1,&(map->ramax),&(map->decmax));
  
  if (map->ramin>map->ramax) {
    actData tmp=map->ramin;
    map->ramin=map->ramax;
    map->ramax=tmp;
  }
  
  if (map->decmin>map->decmax) {
    actData tmp=map->decmin;
    map->decmin=map->decmax;
    map->decmax=tmp;
  }
  

  free(map->map);
  map->map=(actData *)malloc(sizeof(actData)*map->npix);

  //printf("map limits are %14.5f %14.5f %14.5f %14.5f\n",map->ramin,map->ramax,map->decmin,map->decmax);

  
  return 0;
}
