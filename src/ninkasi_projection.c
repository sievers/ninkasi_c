//Module to turn RA/Dec into map coordinates.
#ifndef MAKEFILE_HAND
#include "config.h"
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifdef USE_HEALPIX
#include "chealpix.h"
#define PI_OVER_TWO 1.5707963267948966
#endif

#include "ninkasi.h"
#include "ninkasi_projection.h"

/*--------------------------------------------------------------------------------*/
#if 0
static inline actData sin5(actData x)
{
  
  actData x2=x*x;
  actData x4=x2*x2;

  //Think this one is better, but sticking with old for not-changingness sake
  //return x*(0.999999976513754 -0.166666475934889*x2+ 0.00833289922282484*x4 -0.000198008653066048*x2*x4+ 2.59043003161813e-06*x4*x4); 
  return x*(0.999999995715839  -0.166666579699042*x2 + 0.00833305061731921*x4  -0.000198090463568892*x2*x4 + 2.6051662751399e-06*x4*x4);
  
}
#endif
/*--------------------------------------------------------------------------------*/
#if 0
static inline actData sin4(actData x)
{
  actData x2=x*x;
  actData x4=x2*x2;

  //Think this one is better, sticking w/ old for not-changingness sake.
  //return x*(0.999996601050172  -0.16664823561673*x2  + 0.00830628614181227*x4     -0.000183627485767557*x2*x4);
  return x*(0.999999241345692  -0.166656796188478*x2 + 0.00831322507990857*x4   -0.000185234483301391*x2*x4 );

}
#endif


/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/
void radecvec2cea_pix(const actData *ra, const actData *dec, int *rapix, int *decpix, int *ind, int ndata, const MAP *map)
{
  double rafac=RAD2DEG/map->projection->radelt;  //180/pi since ra is in radians, but CEA FITS likes degrees.                                                                                          
  double decfac=RAD2DEG/map->projection->pv/map->projection->decdelt;
  actData dra=map->projection->rapix-1+0.5;
  actData ddec=map->projection->decpix-1+0.5;
  if ((rapix!=NULL)&&(ind!=NULL)) {
    assert(decpix!=NULL);
    for (int i=0;i<ndata;i++) {
      rapix[i]=ra[i]*rafac+dra;
      decpix[i]=sin5(dec[i])*decfac+ddec;
      ind[i]=map->nx*decpix[i]+rapix[i];
    }
    return;
  }
  if ((rapix!=NULL)&&(ind==NULL))  {
    assert(decpix!=NULL);
    for (int i=0;i<ndata;i++) {
      rapix[i]=ra[i]*rafac+dra;
      decpix[i]=sin5(dec[i])*decfac+ddec;
    }
    return;
  }
  if ((rapix==NULL)&&(ind!=NULL)) {
    for (int i=0;i<ndata;i++) {
      int tmp_ra=ra[i]*rafac+dra;
      int tmp_dec=sin5(dec[i])*decfac+ddec;
      ind[i]=map->nx*tmp_dec+tmp_ra;
    }
    return;
  }
  assert(1==0);  //should never get here, means we didn't find an appropriate set of jobs to carry out.
  return;

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
  //fprintf(stderr,"Working on detector %d\n",det);
  if (tod->pixelization_saved) {
    memcpy(ind,tod->pixelization_saved[det],sizeof(int)*tod->ndata);
    return;
  }
  get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
#if 1
  convert_radec_to_map_pixel(scratch->ra,scratch->dec,ind,tod->ndata,map);
  if (inbounds) {  //checking map pixellization inboundness
    printf("checking pixellization.\n");
    for (int i=0;i<tod->ndata;i++) {
      if ((ind[i]<0)||(ind[i]>=map->npix)) {
	fprintf(stderr,"We are going to have a problem at ra/dec %14.6f %14.6f mapped to pixel %d\n",scratch->ra[i],scratch->dec[i],ind[i]);
	inbounds[i]=false;
      }
      else
	inbounds[i]=true;
    }
  }
  
  return;
#endif
  //printf("using old code.\n");

  //printf("back from get_radec_from_altaz_fit_1det_coarse.\n");
  switch(map->projection->proj_type) {
  case(NK_RECT): 
    //printf("Doing rectangular projection.\n");
    for (int i=0;i<tod->ndata;i++)
      ind[i]=(int)((scratch->dec[i]-map->decmin)/map->pixsize)+map->ny*(int)((scratch->ra[i]-map->ramin)/map->pixsize);    
    break;
    
  case(NK_TAN):
    {
      actData x,y;
      for (int i=0;i<tod->ndata;i++) {
	radec2xy_tan(&x,&y,scratch->ra[i],scratch->dec[i],map->projection);
	ind[i]=(int)(x+0.5)+map->nx*((int)(y+0.5));
      }
    }
    break;
#ifdef USE_HEALPIX
  case(NK_HEALPIX_RING): 
    {
      //int iter=0;
      for (int i=0;i<tod->ndata;i++) {
	/*	
	if (iter==100) {
	  iter=0;
	  fprintf(stderr,"i is %d\n",i);
	}
	if (i>331550)
	  printf("getting pixel for %14.6f %14.6f\n",scratch->dec[i],scratch->ra[i]);
	*/

	//ang2pix_ring(map->projection->nside,PI_OVER_TWO-scratch->dec[i],scratch->ra[i],&ind[i]);
	long tmp;
	ang2pix_ring(map->projection->nside,PI_OVER_TWO-scratch->dec[i],scratch->ra[i],&tmp);
	ind[i]=tmp;
	/*
	if (i>331550)
	  printf("pixel is %d\n",ind[i]);
	iter++;
	*/
      }
    }
    break;
    
  case(NK_HEALPIX_NEST): 
    for (int i=0;i<tod->ndata;i++) {
      //printf("getting pixel for %14.6f %14.6f\n",scratch->dec[i],scratch->ra[i]);

      //ang2pix_nest(map->projection->nside,PI_OVER_TWO-scratch->dec[i],scratch->ra[i],&ind[i]);
      long tmp;
      ang2pix_nest(map->projection->nside,PI_OVER_TWO-scratch->dec[i],scratch->ra[i],&tmp);
      ind[i]=tmp;
      //printf("pixel is %d\n",ind[i]);
      
    }
    break;
    
#endif
  case (NK_CEA):
    {
      //printf("Doing CEA projection.\n");
#if 1  //use the vectorized call.  This lets you calculate pixels (including x and y) in other places as well
      radecvec2cea_pix(scratch->ra,scratch->dec, NULL,NULL,ind,tod->ndata,map);
#else
      
      
      double rafac=RAD2DEG/map->projection->radelt;  //180/pi since ra is in radians, but CEA FITS likes degrees.
      double decfac=RAD2DEG/map->projection->pv/map->projection->decdelt;
      for (int i=0;i<tod->ndata;i++) {
	int rapix=scratch->ra[i]*rafac+map->projection->rapix-1+0.5;
	int decpix=sin5(scratch->dec[i])*decfac+map->projection->decpix-1+0.5;  //change!  13 Aug 2010, should be faster, good to 1e-3 arcsec
	//int decpix=sin(scratch->dec[i])*decfac+map->projection->decpix-1+0.5;  //change!  13 Aug 2010, should be faster, good to 1e-3 arcsec
	ind[i]=map->nx*decpix+rapix;
      }
#endif
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
  //printf("done.\n");
}

/*--------------------------------------------------------------------------------*/
void radec2pix_cea(MAP *map, actData ra, actData dec, int *rapix, int *decpix)
{
#if 1
  radecvec2cea_pix(&ra, &dec, rapix, decpix, NULL,1,map);
#else
  nkProjection *projection=map->projection;
  double rafac=RAD2DEG/map->projection->radelt;
  double decfac=RAD2DEG/map->projection->pv/map->projection->decdelt;
  *rapix=ra*rafac+projection->rapix-1;
  *decpix=sin(dec)*decfac+projection->decpix-1;
#endif
  
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

  printf("ra crap is %d %14.4e %14.4e\n",rapix,map->projection->rapix,map->projection->radelt);

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
  
  mprintf(stdout,"Map limits in set_map_projection_cea_simple are %14.8f %14.8f %14.8f %14.8f\n",map->ramin,map->ramax,map->decmin,map->decmax);
  map->projection->decdelt=map->pixsize/cos0*RAD2DEG;
  map->projection->radelt= -map->projection->decdelt;
  map->projection->pv=cos0*cos0;
  actData ra0=map->ramax*RAD2DEG/map->projection->radelt;
  actData ra1=map->ramin*RAD2DEG/map->projection->radelt+1;  //+1 is in in case of equality
  actData dec0=sin(map->decmin)*RAD2DEG/map->projection->pv/map->projection->decdelt;
  actData dec1=sin(map->decmax)*RAD2DEG/map->projection->pv/map->projection->decdelt+1;  //+1 is in case of equality
  mprintf(stdout,"Before rounding, ra0/dec0 are %14.8f %14.8f\n",ra0,dec0);
  map->projection->rapix=-(int)ra0;
  map->projection->decpix=-(int)dec0;
  map->nx=ra1-ra0;
  map->ny=dec1-dec0;
  free(map->map);
  map->npix=map->nx*map->ny;
  map->map=(actData *)malloc(sizeof(actData)*map->npix);

  mprintf(stdout,"Offsets are %12.4f %12.4f\n",map->projection->rapix,map->projection->decpix);
  mprintf(stdout,"Pixsizes are %14.8f %14.8f\n",map->projection->radelt,map->projection->decdelt);
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
void radec2xy_tan(actData *x, actData *y,actData ra, actData dec, const nkProjection *proj)
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

#ifdef USE_HEALPIX
int set_map_projection_healpix_ring(MAP *map, int nside) {
  int i=1;
  int npix=12*nside*nside;
  map->projection->proj_type=NK_HEALPIX_RING;
  map->projection->nside=nside;

  map->nx=npix;
  map->ny=1;
  map->npix=npix;
  map->map=(actData *)malloc(sizeof(actData)*map->npix);
  

  return 0;
}
/*--------------------------------------------------------------------------------*/
 
int set_map_projection_healpix_nest(MAP *map, int nside) {
  int i=1;
  int npix=12*nside*nside;
  map->projection->proj_type=NK_HEALPIX_NEST;
  map->projection->nside=nside;
  map->nx=npix;
  map->ny=1;
  map->npix=npix;
  map->map=(actData *)malloc(sizeof(actData)*map->npix);
  //make an unthreaded call to ang2pix_nest as there are statics to be set up
  long idum;
  double theta=0.5;
  double phi=0.5;

  ang2pix_nest(nside,theta,phi,&idum);
  pix2ang_nest(nside,idum,&theta,&phi);

 

  return 0;
}
#endif  //use_healpix
/*--------------------------------------------------------------------------------*/
int set_map_projection_cea_predef( MAP *map,actData radelt, actData decdelt, actData rapix, actData decpix, actData pv, int nra, int ndec)
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
/*--------------------------------------------------------------------------------*/
nkProjection *deres_projection(nkProjection *proj)
{
  nkProjection *proj2=(nkProjection *)malloc(sizeof(nkProjection));
  switch(proj->proj_type) {
  case(NK_CEA):
    proj2->proj_type=proj->proj_type;
    proj2->radelt=proj->radelt*2;
    proj2->decdelt=proj->decdelt*2;
    proj2->ra_cent=proj->ra_cent;
    proj2->dec_cent=proj->dec_cent;
    proj2->rapix=proj->rapix/2+0.25;
    proj2->decpix=proj->decpix/2+0.25;
;
    proj2->pv=proj->pv;
    return proj2;
    break;
  default:
    printf("Unsupported type in deres_projection.\n");
    assert(1==0);
    break;
  }
  return NULL;
}

/*--------------------------------------------------------------------------------*/
nkProjection *upres_projection(nkProjection *proj)
{
  nkProjection *proj2=(nkProjection *)malloc(sizeof(nkProjection));
  switch(proj->proj_type) {
  case(NK_CEA):
    proj2->proj_type=proj->proj_type;
    proj2->radelt=proj->radelt/2;
    proj2->decdelt=proj->decdelt/2;
    proj2->ra_cent=proj->ra_cent;
    proj2->dec_cent=proj->dec_cent;
    proj2->rapix=proj->rapix*2-0.5;
    proj2->decpix=proj->decpix*2-0.5;
    proj2->pv=proj->pv;
    return proj2;
    break;
  default:
    printf("Unsupported type in deres_projection.\n");
    assert(1==0);
    break;
  }
  return NULL;
}
/*--------------------------------------------------------------------------------*/
MAP *extract_subregion_map_cea(MAP *map, actData ramin, actData ramax, actData decmin, actData decmax, int do_copy)
{
  return NULL; //still in progress
#if 0
  int xmin,xmax,ymin,ymax;
  radec2pix_cea(map,ramin,decmin,&xmin,&ymin);
  radec2pix_cea(map,ramax,decmax,&xmax,&ymax);
  MAP *small_map=(MAP *)calloc(sizeof(MAP),1);
  small_map->pixsize=map->pixsize;
  small_map->nx=(xmax-xmin)+1;
  small_map->ny=(ymax-ymin)+1;
  small_map->npix=small_map->nx*small_map->ny;
  small_map->have_locks=0;
  small_map->projection=(nkProjection *)calloc(sizeof(nkProjection),1);
  small_map->projection->radelt=map->projection->radelt;
  small_map->projection->decdelt=map->projection->decdelt;
  small_map->projection->rapix=map->projection->rapix+ymin;
  small_map->projection->decpix=map->projection->decpix+ymin;
  small_map->projection->pv=map->projection->pv;
#ifdef ACTPOL
  small_map->pol_state=map->pol_state;
  small_map->npix*=get_npol_in_map(map);
#endif
  small_map->map=(actData *)malloc(sizeof(actData)*small_map->npix);
  if (do_copy) {
    int fac=get_npol_in_map(small_map);
    for (int i=0;i<small_map->ny;i++) {
    }
    
  }
#endif
}

/*--------------------------------------------------------------------------------*/
void convert_radec_to_map_pixel(const actData *ra, const actData *dec, int *ind, long ndata, const MAP *map)
{
  const nkProjection *proj=map->projection;
  switch(proj->proj_type) {
  case(NK_RECT): 
    //printf("Doing rectangular projection.\n");
    for (long i=0;i<ndata;i++)
      ind[i]=(int)((dec[i]-map->decmin)/map->pixsize)+map->ny*(int)((ra[i]-map->ramin)/map->pixsize);    
    break;
  case(NK_TAN):
    {
      actData x,y;
      for (long i=0;i<ndata;i++) {
        radec2xy_tan(&x,&y,ra[i],dec[i],proj);
        ind[i]=(int)(x+0.5)+map->nx*((int)(y+0.5));
      }
    }
    break;
#ifdef USE_HEALPIX
  case(NK_HEALPIX_RING):
    for (int i=0;i<ndata;i++) {
      //ang2pix_ring(map->projection->nside,PI_OVER_TWO-dec[i],ra[i],&ind[i]);
      long tmp;
      ang2pix_ring(map->projection->nside,PI_OVER_TWO-dec[i],ra[i],&tmp);
      ind[i]=tmp;
    }
    break;
  case(NK_HEALPIX_NEST):
    for (int i=0;i<ndata;i++) {
      //ang2pix_nest(map->projection->nside,PI_OVER_TWO-dec[i],ra[i],&ind[i]); 
      long tmp;
      ang2pix_nest(map->projection->nside,PI_OVER_TWO-dec[i],ra[i],&tmp); 
      ind[i]=tmp;
    }
    break;
#endif
  case (NK_CEA):
    radecvec2cea_pix(ra,dec, NULL,NULL,ind,ndata,map);
    break;
  default:
    printf("Unknown type in convert_radec_to_map_pixel.\n");
    assert(1==0);
    break;
  }
  
}

/*--------------------------------------------------------------------------------*/
void convert_saved_pointing_to_pixellization(mbTOD *tod, MAP *map)
//convert the ra/dec saved in a TOD to a map pixellization.  Free the RA/Dec.
{
  if ((!tod->ra_saved) ||(!tod->dec_saved)) {
    fprintf(stderr,"Missing ra/dec in TOD in convert_saved_pointing_to_pixellization.\n");
    return;
  }
  if (!tod->pixelization_saved)
    tod->pixelization_saved=imatrix(tod->ndet,tod->ndata);
#pragma omp parallel for shared(tod,map) default(none)
  for (int i=0;i<tod->ndet;i++)
    convert_radec_to_map_pixel(tod->ra_saved[i],tod->dec_saved[i],tod->pixelization_saved[i],tod->ndata,map);
  
  free(tod->ra_saved[0]);
  free(tod->ra_saved);
  tod->ra_saved=NULL;

  free(tod->dec_saved[0]);
  free(tod->dec_saved);
  tod->dec_saved=NULL;
  
}
