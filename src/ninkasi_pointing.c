#include "config.h"

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "ninkasi.h"
#include "ninkasi_mathutils.h"
#include "ninkasi_pointing.h"
#include "mbTOD.h"

//#define OLD_POINTING_OFFSET

#define BAD_POINTING 100
#define ALTAZ_SPACING 0.003  //sample every 0.003 radians (~10 arcmin) when doing pointing fit




/*--------------------------------------------------------------------------------*/

void
ACTSite( Site *p )
{
  p->latitude = deg2rad(-22.9585);
  p->east_longitude = deg2rad(-67.7876);
  p->elevation_m = 5188.;
  p->temperature_K = 273.;
  p->pressure_mb = 550.;
  p->relative_humidity = 0.2;
}




/*--------------------------------------------------------------------------------*/
int
act_observed_altaz_to_mean_radec( const Site *site, double freq_GHz,
        int n, const double ctime[], const actData alt[], const actData az[],
	  actData ra[], actData dec[] )
{
#ifdef ACTDATA_DOUBLE
  return dobserved_altaz_to_mean_radec(site,freq_GHz,n,ctime,alt,az,ra,dec);
#else
  return observed_altaz_to_mean_radec(site,freq_GHz,n,ctime,alt,az,ra,dec);
#endif

}



/*--------------------------------------------------------------------------------*/
actData get_az_offset(const mbTOD *tod, int det)
{
  if (tod)
    if (tod->pointingOffset)
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[det],tod->cols[det])) 
	return tod->pointingOffset->offsetAzCosAlt[tod->rows[det]][tod->cols[det]];
  return ACT_NO_VALUE;
}
/*--------------------------------------------------------------------------------*/
actData get_alt_offset(const mbTOD *tod, int det)
{
  if (tod)
    if (tod->pointingOffset)
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[det],tod->cols[det])) 
	return tod->pointingOffset->offsetAlt[tod->rows[det]][tod->cols[det]];
  return ACT_NO_VALUE;
}

/*--------------------------------------------------------------------------------*/
void destroy_pointing_offset(mbPointingOffset *pt)
{
  free_matrix(pt->offsetAlt);
  free_matrix(pt->offsetAzCosAlt);
  if (pt->nparam >0)
    free(pt->fit);
  free(pt);
}
/*--------------------------------------------------------------------------------*/

  mbPointingOffset *nkPointingOffsetAlloc(int nrow,int ncol,int fitType) 
{
  mbPointingOffset *crud= ( mbPointingOffset *)calloc(1,sizeof(  mbPointingOffset ));
  assert(crud);
  crud->offsetAlt=matrix(nrow,ncol);
  crud->offsetAzCosAlt=matrix(nrow,ncol);
  int nparam=0;
  switch(fitType) {
  case(0):
    nparam=2;
    break;
  default:
    fprintf(stderr,"Unrecognized fitType in mbPointingOffsetAlloc - type requested is %d\n",fitType);
    exit(EXIT_FAILURE);
  }
  crud->fit=vector(nparam);
  crud->nrow=nrow;
  crud->ncol=ncol;
  crud->nparam=nparam;
  crud->fitType=fitType;
  return crud;
}

/*--------------------------------------------------------------------------------*/
void cut_mispointed_detectors(mbTOD *tod)
{
  assert(tod);
  if (!tod->pointingOffset)
    return;
  for (int i=0;i<tod->ndet;i++) {
    if ((fabs(get_az_offset(tod,i))>BAD_POINTING)||(fabs(get_alt_offset(tod,i))>BAD_POINTING))
      mbCutsSetAlwaysCut(tod->cuts, tod->rows[i], tod->cols[i]);
  }
  
}

/*--------------------------------------------------------------------------------*/
mbPointingOffset *nkReadPointingOffset(const char *filename)  {
  char *bigline=read_all_stdin(filename);
  assert(bigline);  //if this fails, we were sent a bad filename.
  //printf("read bigline.\n");
  //printf("length is %ld\n",strlen(bigline));
  int argc;
  char **argv=create_argv_new(bigline,&argc," =\n");
  //printf("argc is %d\n",argc);



  int *found_list=(int *)calloc(argc,sizeof(int));
  int nrow=atoi(find_argument(argc,argv,"nrow",found_list));
  int ncol=atoi(find_argument(argc,argv,"ncol",found_list));
  //printf("nrow and ncol are %d %d\n",nrow,ncol);

  assert(nrow>0);
  assert(ncol>0);
  char *crud=find_argument(argc,argv,"fitType",found_list);
  assert(crud);
  int fitType=atoi(crud);
  mbPointingOffset *offset = nkPointingOffsetAlloc(nrow,ncol,fitType);



  
  int fitType_ind=find_tok(argc,argv,"fitType");
  int ind_start=fitType_ind+offset->nparam+2;  //plus two is since we have to skip over fitType key as well
  //printf("ind_start is %d\n",ind_start);
  assert(ind_start+2*(offset->nrow+offset->ncol)-1 <=argc);  //if this fails, the file was too short.

  for (int i=0;i<offset->nparam;i++) {
    int ii=i+fitType_ind+2;
    offset->fit[i]=atoi(argv[ii]);
    found_list[ii]=1;
    
    
  }
  
  for (int row=0;row<offset->nrow;row++)
    for (int col=0;col<offset->ncol;col++) {
      int ii=2*(col+row*(offset->ncol))+ind_start;
      offset->offsetAlt[row][col]=atof(argv[ii]);
      found_list[ii]=1;
      ii++;
      offset->offsetAzCosAlt[row][col]=atof(argv[ii]);
      found_list[ii]=1;
    }
  
  for (int i=0;i<argc;i++) {
    if (!found_list[i]) {
      fprintf(stderr,"Warning - unrecognized command %s in position %d in nkReadPointingOffset.\n",argv[i],i);
    }
  }

  free(found_list);
  free(bigline);
  free_argv(argc,argv);
  return offset;
  
}


/*--------------------------------------------------------------------------------*/
void assign_tod_ra_dec_old(mbTOD *tod)
{
  assert(tod);
  assert(tod->az);
  assert(tod->alt);
  assert(tod->ndata>0);
  assert(tod->ctime>0);
  assert(tod->deltat>0);
  double *ctime=(double *)malloc(sizeof(double)*tod->ndata);
  for (int i=0;i<tod->ndata;i++) 
    ctime[i]=tod->ctime+tod->deltat*((double)i);
  
  assert(tod->ra==NULL);  //you shouldn't have anything here coming in, most likely.  Still, if
  assert(tod->dec==NULL); //you kill these, code will overwrite if ptrs not null, else alloc.
  if (tod->ra==NULL)
    tod->ra=vector(tod->ndata);
  if (tod->dec==NULL)
    tod->dec=vector(tod->ndata);
  Site site;
  ACTSite(&site);
  act_observed_altaz_to_mean_radec(&site,150.0,tod->ndata,ctime,tod->alt,tod->az,tod->ra,tod->dec);
#if 0
  write_timestream(tod->ra,tod->ndata,"ra.out");
  write_timestream(tod->dec,tod->ndata,"dec.out");
  write_timestream(tod->az,tod->ndata,"az.out");
  write_timestream(tod->alt,tod->ndata,"alt.out");
  write_double_timestream(ctime,tod->ndata,"ctime.out");
#endif



  free(ctime);
}
/*--------------------------------------------------------------------------------*/
actData inbounds_ra_element(actData ra, actData ra0)
{

  while (ra-ra0>M_PI)
    ra-=2*M_PI;
  while (ra-ra0 < -M_PI)
    ra+= 2*M_PI;
  return ra;
}

/*--------------------------------------------------------------------------------*/
void unwrap_ra(actData *ra, int n)
//make sure that RA doesn't wrap around 0
{
  for (int i=1;i<n;i++) 
    ra[i]=inbounds_ra_element(ra[i],ra[0]);
}

/*--------------------------------------------------------------------------------*/
void set_tod_pointing_tiled(mbTOD *tod, actData *azvec, int naz, actData *altvec, int nalt, actData **ra_mat, actData **dec_mat, actData *ra_clock, actData *dec_clock, int nclock)
{
  assert(tod);
  assert(tod->pointing_fit);
  TiledPointingFit *tile=(TiledPointingFit *)malloc(sizeof(TiledPointingFit));
  tile->naz=naz;
  tile->nalt=nalt;
  tile->azvec=vector(naz);
  tile->altvec=vector(nalt);
  tile->ra_mat=matrix(naz,nalt);
  tile->dec_mat=matrix(naz,nalt);
  tile->nclock=nclock;
  tile->ra_clock=vector(nclock);
  tile->dec_clock=vector(nclock);
  //printf("naz and nalt are %d %d\n",naz,nalt);
  memcpy(tile->ra_clock,ra_clock,sizeof(actData)*nclock);
  memcpy(tile->dec_clock,dec_clock,sizeof(actData)*nclock);

  memcpy(tile->azvec,azvec,sizeof(actData)*naz);
  memcpy(tile->altvec,altvec,sizeof(actData)*nalt);
  memcpy(tile->ra_mat[0],ra_mat[0],sizeof(actData)*naz*nalt);
  memcpy(tile->dec_mat[0],dec_mat[0],sizeof(actData)*naz*nalt);

  tod->pointing_fit->tiled_fit=tile;

}
/*--------------------------------------------------------------------------------*/
void assign_tod_ra_dec(mbTOD *tod)
{

  assert(tod->pointing_fit==NULL);
  //tod->pointing_fit=(PointingFit *)malloc(sizeof(PointingFit));
  tod->pointing_fit=(PointingFit *)calloc(1,sizeof(PointingFit));


  actData altmin=vecmin(tod->alt,tod->ndata);
  actData azmin=vecmin(tod->az,tod->ndata);
  actData altmax=vecmax(tod->alt,tod->ndata);
  actData azmax=vecmax(tod->az,tod->ndata);

  actData *daz=vector(tod->ndet);
  actData *dalt=vector(tod->ndet);
  int ndet_use=0;
  int nrow=ivecmax(tod->rows,tod->ndet)+1;
  int ncol=ivecmax(tod->cols,tod->ndet)+1;

  //printf("az/alt[0] is %14.5e %14.5e\n",tod->alt[0],tod->az[0]);


  for (int i=0;i<tod->ndet;i++) {
    if ( !mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
      daz[ndet_use]=get_az_offset(tod,i);
      dalt[ndet_use]=get_alt_offset(tod,i);
      ndet_use++;
    }
  }
  actData dalt_max=vecmax(dalt,ndet_use);
  actData dalt_min=vecmin(dalt,ndet_use);
  actData daz_max=vecmax(daz,ndet_use);
  actData daz_min=vecmin(daz,ndet_use);

  //find out the grid we're going to sample the center of the array on.
  // the plus 2 is because we always include the min and max points.
  int naz=2+((int)((azmax-azmin)/ALTAZ_SPACING));
  int nalt=2+((int)((altmax-altmin)/ALTAZ_SPACING));
  int ncenters=naz*nalt;


  actData *azcent_1d=vector(naz);
  actData *altcent_1d=vector(nalt);
  for (int i=0;i<naz-1;i++)
    azcent_1d[i]=azmin+ALTAZ_SPACING*(actData)i;
  azcent_1d[naz-1]=azmax;
  
  for (int i=0;i<nalt-1;i++)
    altcent_1d[i]=altmin+ALTAZ_SPACING*(actData)i;
  altcent_1d[nalt-1]=altmax;
  
  int nsamp=ncenters*ndet_use;
  actData *azvec=vector(nsamp);
  actData *altvec=vector(nsamp);
  double *ctime=(double *)malloc(sizeof(double)*nsamp);
  int ind=0;
  for (int iaz=0;iaz<naz;iaz++) 
    for (int ialt=0;ialt<nalt;ialt++)
      for (int idet=0;idet<ndet_use;idet++) {
	altvec[ind]=altcent_1d[ialt]+dalt[idet];
#ifdef OLD_POINTING_OFFSET
	azvec[ind]=azcent_1d[iaz]+daz[idet]/cos(altvec[ind]);
#else
	azvec[ind]=azcent_1d[iaz]+daz[idet]/cos(tod->alt[tod->ndata/2]);
#endif
	ctime[ind]=tod->ctime;
	ind++;
      }
  assert(ind==nsamp);  //sanity check.

  actData *ra=vector(nsamp);
  actData *dec=vector(nsamp);
#if 1
  Site site;
  ACTSite(&site);

  act_observed_altaz_to_mean_radec(&site,150.0,nsamp,ctime,altvec,azvec,ra,dec);  
  

  if (ra[0]>4)
    ra[0]-=2*M_PI;
  unwrap_ra(ra,nsamp);

#if 0
  //hack to see if pointing matters for streaks
  tod->pointing_fit->ra_fit=fit_2d_poly(altvec,azvec,ra,nsamp,NULL,4,5);
  tod->pointing_fit->dec_fit=fit_2d_poly(altvec,azvec,dec,nsamp,NULL,4,5);
#else
  //this is the default.
  tod->pointing_fit->ra_fit=fit_2d_poly(altvec,azvec,ra,nsamp,NULL,2,3);
  tod->pointing_fit->dec_fit=fit_2d_poly(altvec,azvec,dec,nsamp,NULL,2,3);
#endif
  
  actData alt_median=compute_median_inplace(nsamp,altvec);
  actData az_median=compute_median_inplace(nsamp,azvec);
  actData ra_start, ra_stop,dec_start,dec_stop;
  act_observed_altaz_to_mean_radec(&site,150.0,1,&(tod->ctime),&alt_median,&az_median,&ra_start,&dec_start);
  ra_start=inbounds_ra_element(ra_start,ra[0]);
  //if (ra_start>5.0)
  //ra_start -=2*M_PI;
  double ctime_end=tod->ctime+tod->deltat*((double)(tod->ndata-1));
  act_observed_altaz_to_mean_radec(&site,150.0,1,&ctime_end,&alt_median,&az_median,&ra_stop,&dec_stop);
  ra_stop=inbounds_ra_element(ra_stop,ra[0]);
  //  if (ra_stop>5.0)
  //ra_start -=2*M_PI;


  tod->pointing_fit->ra_clock_rate=(ra_stop-ra_start)/((actData)(tod->ndata-1));
  tod->pointing_fit->dec_clock_rate=(dec_stop-dec_start)/((actData)(tod->ndata-1));
  
#endif

  free(azcent_1d);
  free(altcent_1d);

  free(daz);
  free(dalt);
  free(ra);
  free(dec);
  free(azvec);
  free(altvec);
  free(ctime);



  actData *ra_true=vector(tod->ndata);
  actData *dec_true=vector(tod->ndata);
  actData *ra_fit=vector(tod->ndata);    
  actData *dec_fit=vector(tod->ndata);
  actData *az=vector(tod->ndata);
  actData *alt=vector(tod->ndata);


#if 0
  
  ind=0;
  while (mbCutsIsAlwaysCut(tod->cuts,tod->rows[ind],tod->cols[ind]))
    ind++;

  pca_time tt;
  tick(&tt);
  get_radec_from_altaz_fit_1det(tod,ind,alt,az,ra_fit,dec_fit);
  printf("time to do fit:  %8.4f\n",tocksilent(&tt));
  tick(&tt);
  get_radec_from_altaz_exact_1det(tod,ind,alt,az,ra_true,dec_true);
  printf("time to do exact:  %8.4f\n",tocksilent(&tt));
  write_timestream(ra_true,tod->ndata,"myra_true.out");
  write_timestream(ra_fit,tod->ndata,"myra_fit.out");
  
  write_timestream(dec_true,tod->ndata,"mydec_true.out");
  write_timestream(dec_fit,tod->ndata,"mydec_fit.out");
  
  exit(EXIT_SUCCESS);
#endif

}
/*--------------------------------------------------------------------------------*/
int *find_az_turnarounds(mbTOD *tod, int *nturn)
{
  assert(tod);
  assert(tod->az);
  assert(tod->ndata>2);
  int *vec=(int *)calloc(tod->ndata,sizeof(int));
  assert(vec);
  vec[0]=1;
  vec[tod->ndata-1]=1;
  for (int i=1;i<tod->ndata-1;i++) {
    //Do things twice with the equality in different spots.
    //This way, a long run of identicals will have both the beginning
    //and end marked, but nothing in between.
    if ((tod->az[i]>=tod->az[i-1])&&(tod->az[i]>tod->az[i+1]))
      vec[i]=1;
    if ((tod->az[i]>tod->az[i-1])&&(tod->az[i]>=tod->az[i+1]))
      vec[i]=1;
    if ((tod->az[i]<=tod->az[i-1])&&(tod->az[i]<tod->az[i+1]))
      vec[i]=1;
    if ((tod->az[i]<tod->az[i-1])&&(tod->az[i]<=tod->az[i+1]))
      vec[i]=1;
  }
  int nkept=0;
  for (int i=0;i<tod->ndata;i++)
    if (vec[i])
      nkept++;
  int *turns=(int *)malloc(sizeof(int)*nkept);
  int cur=0;
  for (int i=0;i<tod->ndata;i++) {
    if (vec[i]) {
      turns[cur]=i;
      cur++;
    }
  }
  assert(cur==nkept);
  *nturn=nkept;
  free(vec);
  return turns;
    
}
/*--------------------------------------------------------------------------------*/
static int inbounds_index(int ind, int n)
{
  if (ind<0)
    return 0;
  if (ind>n)
    return n;
  return ind;
}

/*--------------------------------------------------------------------------------*/
#if 0
int *find_az_sample_points(mbTOD *tod, int *nsample)
{
  assert(tod);
  assert(tod->ndata>2);
  assert(tod->az);
  int nturn;
  int *turns=find_az_turnarounds(tod, &nturn);
  int coarse_sample_rate=20;
  int turn_pad=50;
  int *is_measured=(int *)calloc(tod->ndata,sizeof(int));
  is_measured[0]=1;
  is_measured[tod->ndata-1]=1;
  for (int i=0;i<tod->ndata;i+=coarse_sample_rate)
    is_measured[i]=1;
  
  for (int ii=0;ii<nturn;ii++) {
    for (int i=inbounds_ind(turns[ii]-turn_pad,tod->ndata);i<inbounds_ind(turns[ii]+turn_pad+1,tod->ndata);i++) 
      is_measured[i]=1;    
  }
  int n_to_measure=0;
  for (int i=0;i<tod->ndata;i++)
    if (is_measured[i])
      n_to_measure++;
  int *inds=(int *)malloc(sizeof(int)*n_to_measure);
  
    
  free(turns);
  free(is_measured);

}
#endif
/*--------------------------------------------------------------------------------*/
  
void get_uncut_daz_dalt(mbTOD *tod,actData **daz_out,actData **dalt_out,int *ndet)
{
  assert(tod);
  assert(tod->cuts);
  assert(tod->pointingOffset);
  assert(tod->rows);
  assert(tod->cols);
  
  actData *daz=vector(tod->ndet);
  actData *dalt=vector(tod->ndet);
  
  int ndet_good=0;
  for (int i=0;i<tod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
      daz[ndet_good]=get_az_offset(tod,i);
      dalt[ndet_good]=get_alt_offset(tod,i);
      ndet_good++;
    }
  }
  *daz_out=daz;
  *dalt_out=dalt;
  *ndet=ndet_good;
}

/*--------------------------------------------------------------------------------*/
void get_radec_from_altaz_fit_with_ind(mbTOD *tod, actData *alt, actData *az, int *samp, int nsamp, actData *ra, actData *dec)
{
  eval_2d_poly_inplace(alt,az,nsamp,tod->pointing_fit->ra_fit,ra);
  eval_2d_poly_inplace(alt,az,nsamp,tod->pointing_fit->dec_fit,dec);
  
  for (int i=0;i<nsamp;i++) {
    ra[i]+=((actData)samp[i])*tod->pointing_fit->ra_clock_rate;
    dec[i]+=((actData)samp[i])*tod->pointing_fit->dec_clock_rate;
  }

}

/*--------------------------------------------------------------------------------*/
PointingFit *copy_pointing_fit(PointingFit *fit)
{
  PointingFit *fit2=(PointingFit *)malloc(sizeof(PointingFit));
  fit2->ra_clock_rate=fit->ra_clock_rate;
  fit2->dec_clock_rate=fit->dec_clock_rate;
  fit2->ra_fit=copy_2d_polyfit(fit->ra_fit);
  fit2->dec_fit=copy_2d_polyfit(fit->dec_fit);

  fit2->ncoarse=fit->ncoarse;
  if (fit->ncoarse>0) {
    fit2->coarse_ind=(int *)malloc(sizeof(int)*fit->ncoarse);
    memcpy(fit2->coarse_ind,fit->coarse_ind,sizeof(int)*fit->ncoarse);
  }
  else
    fit2->coarse_ind=NULL;
  

  return fit2;
}
/*--------------------------------------------------------------------------------*/

PointingFitScratch *allocate_pointing_fit_scratch(const mbTOD *tod)
{
  assert(tod);
  assert(tod->ndata>0);

  PointingFitScratch *scratch=(PointingFitScratch *)malloc(sizeof(PointingFitScratch));
  scratch->ra=vector(tod->ndata);
  scratch->dec=vector(tod->ndata);
  scratch->alt=vector(tod->ndata);
  scratch->az=vector(tod->ndata);
  scratch->tod_alt=vector(tod->ndata);
  scratch->tod_az=vector(tod->ndata);

  if (tod->az) {  //breaks if the pointing model is saved.
    assert(tod->az);  
    assert(tod->alt);
    memcpy(scratch->tod_alt,tod->alt,sizeof(actData)*tod->ndata);
    memcpy(scratch->tod_az,tod->az,sizeof(actData)*tod->ndata);
  }

  if (tod->pointing_fit)
    scratch->pointing_fit=copy_pointing_fit(tod->pointing_fit);
  else
    scratch->pointing_fit=NULL;
  //scratch->pointing_fit=tod->pointing_fit;
  if (scratch->pointing_fit)
    if (scratch->pointing_fit->ncoarse>0) {
      scratch->az_coarse=vector(scratch->pointing_fit->ncoarse);
      scratch->alt_coarse=vector(scratch->pointing_fit->ncoarse);
      scratch->ra_coarse=vector(scratch->pointing_fit->ncoarse);
      scratch->dec_coarse=vector(scratch->pointing_fit->ncoarse);
      scratch->time_coarse=vector(scratch->pointing_fit->ncoarse);
    }

  
  return scratch;
  
}

/*--------------------------------------------------------------------------------*/
void destroy_pointing_fit_scratch(PointingFitScratch *scratch)
{
  free(scratch->ra);
  free(scratch->dec);
  free(scratch->alt);
  free(scratch->az);
  free(scratch->tod_alt);
  free(scratch->tod_az);
  if (scratch->pointing_fit)
    if (scratch->pointing_fit->ncoarse) {
      free(scratch->ra_coarse);
      free(scratch->dec_coarse);
      free(scratch->alt_coarse);
      free(scratch->az_coarse);
      free(scratch->time_coarse);
    }
  
  if (scratch->pointing_fit)
    destroy_pointing_fit_raw(scratch->pointing_fit);
  free(scratch);
  
}
/*--------------------------------------------------------------------------------*/
void apply_time_ramp_to_fit(const mbTOD *tod, PointingFitScratch *scratch)
{  
  for (int i=0;i<tod->ndata;i++) {
    scratch->ra[i]+=((actData)i)*tod->pointing_fit->ra_clock_rate;
    scratch->dec[i]+=((actData)i)*tod->pointing_fit->dec_clock_rate;
  }
}  

/*--------------------------------------------------------------------------------*/
void get_radec_from_altaz_fit_tiled(const TiledPointingFit *tile, actData *alt, actData *az, actData *elapsed_time, actData *ra, actData *dec, int npoint)
{
  actData daz=tile->azvec[1]-tile->azvec[0];
  actData dalt=tile->altvec[1]-tile->altvec[0];
  actData az0=tile->azvec[0];
  actData alt0=tile->altvec[0];
  for (int i=0;i<npoint;i++) {
    actData az_interp=(az[i]-az0)/daz;
    actData alt_interp=(alt[i]-alt0)/dalt;
    int az_pix=az_interp;
    int alt_pix=alt_interp;
    actData az_frac=az_interp-az_pix;
    actData alt_frac=alt_interp-alt_pix;
    actData v1,v2;
#if 0   //debug
    printf("pixes are %d %d %d %d\n",az_pix,alt_pix,tile->naz, tile->nalt);
    if (az_pix<0)
      return;
    if (az_pix>tile->naz-2)
      return;
    if (alt_pix<0)
      return;
    if (alt_pix>tile->nalt-2)
      return;

#endif

    v1=tile->ra_mat[az_pix][alt_pix]*(1-alt_frac)+alt_frac*tile->ra_mat[az_pix][alt_pix+1];
    v2=tile->ra_mat[az_pix+1][alt_pix]*(1-alt_frac)+alt_frac*tile->ra_mat[az_pix+1][alt_pix+1];
    ra[i]=v1*(1-az_frac)+v2*az_frac;
    
    v1=tile->dec_mat[az_pix][alt_pix]*(1-alt_frac)+alt_frac*tile->dec_mat[az_pix][alt_pix+1];
    v2=tile->dec_mat[az_pix+1][alt_pix]*(1-alt_frac)+alt_frac*tile->dec_mat[az_pix+1][alt_pix+1];
    dec[i]=v1*(1-az_frac)+v2*az_frac;
    
    actData tt=1;
    for (int j=0;j<tile->nclock;j++) {
      ra[i]+=tt*tile->ra_clock[j];
      dec[i]+=tt*tile->dec_clock[j];
      tt*=elapsed_time[i];      
    }
  }
}
/*--------------------------------------------------------------------------------*/
void get_radec_from_altaz_fit_1det_coarse(const mbTOD *tod, int det, PointingFitScratch *scratch)
{

#if 0
  //drop this in to do exact pointing
  get_radec_from_altaz_exact_1det(tod,det,scratch);
  return;
#endif

#if 0
  //drop this in to do approximately exact pointing

  get_radec_from_altaz_fit_1det_coarse_exact(tod,det,scratch);
  return;
#endif

  //If we have saved full pointing information, copy it into *scratch here.
  assert(tod);
  if (tod->ra_saved) {
    assert(tod->dec_saved);
    memcpy(scratch->ra,tod->ra_saved[det],tod->ndata*sizeof(actData));
    memcpy(scratch->dec,tod->dec_saved[det],tod->ndata*sizeof(actData));
    return;
  }

  assert(tod->pointing_fit);
  assert(tod->pointingOffset);
  actData mydalt=get_alt_offset(tod,det);
  actData mydaz=get_az_offset(tod,det);
  actData myalt=tod->alt[0];
  actData myaz=tod->az[0];

#ifdef OLD_POINTING_OFFSET
  actData fac=mydaz/cos(scratch->tod_alt[0]+mydalt);
#else
  actData fac=mydaz/cos(scratch->tod_alt[0]);
#endif
  int ncoarse=scratch->pointing_fit->ncoarse;
  assert(ncoarse>0);
  int *ind=scratch->pointing_fit->coarse_ind;

  for (int i=0;i<ncoarse;i++) {

#ifdef OLD_POINTING_OFFSET
    scratch->alt_coarse[i]=scratch->tod_alt[ind[i]]+mydalt;
    scratch->az_coarse[i]=scratch->tod_az[ind[i]]+mydaz/cos(scratch->alt_coarse[i]);
#else
    actData alt0=scratch->tod_alt[ind[i]];
    scratch->alt_coarse[i]=alt0+mydalt;
    scratch->az_coarse[i]=scratch->tod_az[ind[i]]+mydaz/cos(alt0);
#endif
    scratch->time_coarse[i]=tod->deltat*(ind[i]);
  }
  




  if (tod->pointing_fit->tiled_fit) {
    get_radec_from_altaz_fit_tiled(tod->pointing_fit->tiled_fit,scratch->alt_coarse, scratch->az_coarse, scratch->time_coarse,scratch->ra_coarse, scratch->dec_coarse, ncoarse);
    for (int i=0;i<ncoarse-1;i++) {
      scratch->ra[ind[i]]=scratch->ra_coarse[i];
      scratch->dec[ind[i]]=scratch->dec_coarse[i];     
      for (int j=ind[i]+1;j<ind[i+1];j++) {
	scratch->ra[j]=scratch->ra_coarse[i]+((actData)(j-ind[i]))/((actData)(ind[i+1]-ind[i]))*(scratch->ra_coarse[i+1]-scratch->ra_coarse[i]);
	scratch->dec[j]=scratch->dec_coarse[i]+((actData)(j-ind[i]))/((actData)(ind[i+1]-ind[i]))*(scratch->dec_coarse[i+1]-scratch->dec_coarse[i]);
      } 
      scratch->ra[tod->ndata-1]=scratch->ra_coarse[ncoarse-1];
      scratch->dec[tod->ndata-1]=scratch->dec_coarse[ncoarse-1];
      
    }
  }
  else {
    eval_2d_poly_pair_inplace(scratch->alt_coarse,scratch->az_coarse,ncoarse,scratch->pointing_fit->ra_fit,scratch->ra_coarse,scratch->pointing_fit->dec_fit,scratch->dec_coarse);
    
    for (int i=0;i<ncoarse-1;i++) {
      scratch->ra[ind[i]]=scratch->ra_coarse[i] +((actData)ind[i])*tod->pointing_fit->ra_clock_rate;
      scratch->dec[ind[i]]=scratch->dec_coarse[i]+((actData)ind[i])*tod->pointing_fit->dec_clock_rate;
      
      for (int j=ind[i]+1;j<ind[i+1];j++) {
	scratch->ra[j]=scratch->ra_coarse[i]+((actData)(j-ind[i]))/((actData)(ind[i+1]-ind[i]))*(scratch->ra_coarse[i+1]-scratch->ra_coarse[i])+((actData)j)*tod->pointing_fit->ra_clock_rate;
	scratch->dec[j]=scratch->dec_coarse[i]+((actData)(j-ind[i]))/((actData)(ind[i+1]-ind[i]))*(scratch->dec_coarse[i+1]-scratch->dec_coarse[i])+((actData)j)*tod->pointing_fit->dec_clock_rate;
      }
      
    }
    
    scratch->ra[tod->ndata-1]=scratch->ra_coarse[ncoarse-1]+((actData)tod->ndata-1)*tod->pointing_fit->ra_clock_rate;
    scratch->dec[tod->ndata-1]=scratch->dec_coarse[ncoarse-1]+((actData)tod->ndata-1)*tod->pointing_fit->dec_clock_rate;
  }

  
  //fprintf(stderr,"checking pointing offsets.\n");
  
#if 0
  if (det==0) {
    FILE *outfile=fopen("pointing_blob.txt","w"); 
    for (int i=0;i<ncoarse;i++) {
      fprintf(outfile,"%4d %14.6e %14.6e\n",ind[i],scratch->ra_coarse[i]+((actData)ind[i])*tod->pointing_fit->ra_clock_rate,scratch->dec_coarse[i]+((actData)ind[i])*tod->pointing_fit->dec_clock_rate);
    }
    fclose(outfile);
  }
#endif
  
  
  if (tod->pointing_fit->ra_offset!=0) {
    for (int i=0;i<tod->ndata;i++)
      scratch->ra[i]+=tod->pointing_fit->ra_offset;
  }
  
  if (tod->pointing_fit->dec_offset!=0) {
    //printf("adjusting dec.\n");
    for (int i=0;i<tod->ndata;i++)
      scratch->dec[i]+=tod->pointing_fit->dec_offset;
  }

  

  
  //  if (tod->pointing_fit->dec_offset!=0) {
  //for (int i=0;i<ncoarse;i++)
  //  scratch->dec_coarse[i]+=tod->pointing_fit->dec_offset;
  //}



  //apply_time_ramp_to_fit(tod,scratch);
  return;  

}

/*--------------------------------------------------------------------------------*/
void get_radec_from_altaz_fit_1det_coarse_exact(const mbTOD *tod, int det, PointingFitScratch *scratch)
//do the full evaluation on the pivot points.  Should help a lot with speed, keep accuracy to ~1"
{

  //drop this in to do exact pointing
#if 0
  get_radec_from_altaz_exact_1det(tod,det,scratch);
  return;
#endif

  assert(tod);
  assert(tod->pointing_fit);
  assert(tod->pointingOffset);
  actData mydalt=get_alt_offset(tod,det);
  actData mydaz=get_az_offset(tod,det);
  actData myalt=tod->alt[0];
  actData myaz=tod->az[0];

#ifdef OLD_POINTING_OFFSET
  actData fac=mydaz/cos(scratch->tod_alt[0]+mydalt);
#else
  actData fac=mydaz/cos(scratch->tod_alt[0]);
#endif
  int ncoarse=scratch->pointing_fit->ncoarse;
  assert(ncoarse>0);
  int *ind=scratch->pointing_fit->coarse_ind;

  double *ctime=dvector(ncoarse);
  for (int i=0;i<ncoarse;i++) {
    ctime[i]=tod->ctime+((double)ind[i])*tod->deltat;
#ifdef OLD_POINTING_OFFSET
    scratch->alt_coarse[i]=scratch->tod_alt[ind[i]]+mydalt;
    scratch->az_coarse[i]=scratch->tod_az[ind[i]]+mydaz/cos(scratch->alt_coarse[i]);
#else
    actData alt0=scratch->tod_alt[ind[i]];
    scratch->alt_coarse[i]=alt0+mydalt;
    scratch->az_coarse[i]=scratch->tod_az[ind[i]]+mydaz/cos(alt0);
#endif
  }
  Site site;
  ACTSite (&site);
  act_observed_altaz_to_mean_radec(&site,150.0,ncoarse,ctime,scratch->alt_coarse,scratch->az_coarse,scratch->ra_coarse,scratch->dec_coarse);
  free(ctime);

  //eval_2d_poly_pair_inplace(scratch->alt_coarse,scratch->az_coarse,ncoarse,scratch->pointing_fit->ra_fit,scratch->ra_coarse,scratch->pointing_fit->dec_fit,scratch->dec_coarse);
  //fprintf(stderr,"checking pointing offsets.\n");



#if 0
  if (det==0) {
    FILE *outfile=fopen("pointing_blob.txt","w"); 
    for (int i=0;i<ncoarse;i++) {
      fprintf(outfile,"%4d %14.6e %14.6e\n",ind[i],scratch->ra_coarse[i]+((actData)ind[i])*tod->pointing_fit->ra_clock_rate,scratch->dec_coarse[i]+((actData)ind[i])*tod->pointing_fit->dec_clock_rate);
    }
    fclose(outfile);
  }
#endif
  for (int i=0;i<ncoarse-1;i++) {
    scratch->ra[ind[i]]=scratch->ra_coarse[i];// +((actData)ind[i])*tod->pointing_fit->ra_clock_rate;
    scratch->dec[ind[i]]=scratch->dec_coarse[i];//+((actData)ind[i])*tod->pointing_fit->dec_clock_rate;
    
    for (int j=ind[i]+1;j<ind[i+1];j++) {
      scratch->ra[j]=scratch->ra_coarse[i]+((actData)(j-ind[i]))/((actData)(ind[i+1]-ind[i]))*(scratch->ra_coarse[i+1]-scratch->ra_coarse[i]);//+((actData)j)*tod->pointing_fit->ra_clock_rate;
      scratch->dec[j]=scratch->dec_coarse[i]+((actData)(j-ind[i]))/((actData)(ind[i+1]-ind[i]))*(scratch->dec_coarse[i+1]-scratch->dec_coarse[i]);//+((actData)j)*tod->pointing_fit->dec_clock_rate;
    }
    
  }

  scratch->ra[tod->ndata-1]=scratch->ra_coarse[ncoarse-1]+((actData)tod->ndata-1)*tod->pointing_fit->ra_clock_rate;
  scratch->dec[tod->ndata-1]=scratch->dec_coarse[ncoarse-1]+((actData)tod->ndata-1)*tod->pointing_fit->dec_clock_rate;


  if (tod->pointing_fit->ra_offset!=0) {
    for (int i=0;i<tod->ndata;i++)
      scratch->ra[i]+=tod->pointing_fit->ra_offset;
  }

  if (tod->pointing_fit->dec_offset!=0) {
    //printf("adjusting dec.\n");
    for (int i=0;i<tod->ndata;i++)
      scratch->dec[i]+=tod->pointing_fit->dec_offset;
  }

  //  if (tod->pointing_fit->dec_offset!=0) {
  //for (int i=0;i<ncoarse;i++)
  //  scratch->dec_coarse[i]+=tod->pointing_fit->dec_offset;
  //}



  //apply_time_ramp_to_fit(tod,scratch);
  return;  

}

/*--------------------------------------------------------------------------------*/


void get_radec_from_altaz_fit_1det(const mbTOD *tod,int det,  PointingFitScratch *scratch)
//set alt and az, then evaluate the polynomial fit.
{
  assert(tod);
  assert(tod->pointing_fit);
  assert(tod->pointingOffset);
  actData mydalt=get_alt_offset(tod,det);
  actData mydaz=get_az_offset(tod,det);
  actData myalt=tod->alt[0];
  actData myaz=tod->az[0];
#ifdef OLD_POINTING_OFFSET
  actData fac=mydaz/cos(scratch->tod_alt[0]+mydalt);
#else
  actData fac=mydaz/cos(scratch->tod_alt[0]);
#endif
  for (int i=0;i<tod->ndata;i++) {
    scratch->alt[i]=scratch->tod_alt[i]+mydalt;
#ifdef OLD_POINTING_OFFSET
    scratch->az[i]=scratch->tod_az[i]+mydaz/cos(scratch->alt[i]);
#else
    scratch->az[i]=scratch->tod_az[i]+mydaz/cos(scratch->tod_alt[i]);
#endif
    //az[i]=tod->az[i]+fac;
    //alt[i]=myalt;
    //az[i]=myaz+fac;
  }
  eval_2d_poly_pair_inplace(scratch->alt,scratch->az,tod->ndata,scratch->pointing_fit->ra_fit,scratch->ra,scratch->pointing_fit->dec_fit,scratch->dec);
  
  //eval_2d_poly_inplace(alt,az,tod->ndata,tod->pointing_fit->ra_fit,ra);
  //eval_2d_poly_inplace(alt,az,tod->ndata,tod->pointing_fit->dec_fit,dec);

  apply_time_ramp_to_fit(tod,scratch);

}
/*--------------------------------------------------------------------------------*/
actData get_interp_segment_az_max_err(actData *az, int i1, int i2)
{
  actData max_err=0;
  actData idelt=i2-i1;
  actData fac=(az[i2]-az[i1])/idelt;
  for (int i=i1+1;i<i2;i++) {
    actData di=i-i1;
    actData err=fabs(az[i1]+di*fac-az[i]);
    if (err>max_err)
      max_err=err;
  }
  return max_err;
}
/*--------------------------------------------------------------------------------*/
void find_pointing_pivots(mbTOD *tod, actData tol)
{
  //printf("pivot tolerance is %8.3g\n",tol);
  if (tol<=0)
    tol=1.0;  //1" threshold.
  tol=tol/3600.0/180.0*M_PI/cos(tod->alt[0]);  //turn arcsecond pointing tolerance into Radians
  //printf("tol is %14.4e\n",tol);

#if 1
  int n=tod->ndata;
  int *am_i_pivot=(int *)calloc(n,sizeof(int));
  am_i_pivot[0]=1;
  am_i_pivot[n-1]=1;
  int i1=1;
  int i2=2;

  while (1) {
    while ((get_interp_segment_az_max_err(tod->az, i1, i2)<tol)&&(i2<n)) {
      //printf("i2 is %d\n",i2);
      //printf("maxerr is %14.5e on %d %d \n",get_interp_segment_az_max_err(tod->az, i1, i2),i1,i2);
      i2++;
      //if (i2==500)
      //exit(EXIT_SUCCESS);
    
    }
    if (i2==n)
      break;
    i1=i2-1;
    am_i_pivot[i1]=1;
    //printf("added pivot at %d\n",i1);
  }
  int npivot=0;
  for (int i=0;i<n;i++)
    if (am_i_pivot[i])
      npivot++;
  int *coarse_ind=(int *)malloc(sizeof(int)*npivot);
  int icur=0;
  for (int i=0;i<n;i++)
    if (am_i_pivot[i]) {
      coarse_ind[icur]=i;
      icur++;
    }
  tod->pointing_fit->ncoarse=npivot;
  tod->pointing_fit->coarse_ind=coarse_ind;
  free(am_i_pivot);
#if 0
  FILE *outfile=fopen("my_pivots.txt","w");
  for (int i=0;i<npivot;i++) {
    fprintf(outfile,"%5d\n",coarse_ind[i]);
  }
  fclose(outfile);
#endif
  return;
#else
  int samp_size=15;
  int nsamp=tod->ndata/samp_size;
  if (nsamp*samp_size<tod->ndata)
    nsamp++;
  int *coarse_ind=(int *)malloc(sizeof(int)*nsamp);
  for (int i=0;i<nsamp;i++)
    coarse_ind[i]=i*samp_size;
  coarse_ind[nsamp-1]=tod->ndata-1;
  
  tod->pointing_fit->ncoarse=nsamp;
  tod->pointing_fit->coarse_ind=coarse_ind;

#endif

}
/*--------------------------------------------------------------------------------*/
actData find_max_pointing_err(const mbTOD *tod)
{
  actData max_err=0;
  int ndet_used=0;
  actData tot_err=0;
#pragma omp parallel shared(max_err, tod) reduction(+:ndet_used, tot_err) default(none)
  {
    
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    PointingFitScratch *scratch_exact=allocate_pointing_fit_scratch(tod);
    actData my_max_err=0;
#pragma omp for schedule(dynamic, 1)
    for (int i=0;i<tod->ndet;i++) 
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	ndet_used++;
	get_radec_from_altaz_fit_1det_coarse(tod,i,scratch);
	get_radec_from_altaz_exact_1det(tod,i,scratch_exact);
	actData mycos=cos(scratch_exact->dec[0]);
	for (int j=0;j<tod->ndata;j++) {
	  actData dra=(scratch->ra[i]-scratch_exact->ra[i]);
	  if (dra>M_PI)
	    dra-= 2*M_PI;
	  if (dra<-M_PI)
	    dra+= 2*M_PI;
	  dra *= mycos;
	  actData ddec=scratch->dec[i]-scratch_exact->dec[i];
	  actData err=sqrt(dra*dra+ddec*ddec);
	  tot_err+=err;
	  if (err>my_max_err)
	    my_max_err=err;
#if 0
	  if (err>max_err)
#pragma omp critical
	    {
	      if (err>max_err) {
		max_err=err;
	      }
	      
	    }
#endif
	}
      }
#pragma omp critical
    if (my_max_err>max_err)
      max_err=my_max_err;

    
  }
  actData ndata=ndet_used*tod->ndata;
  printf("used %d detectors, average err is %12.4e.\n",ndet_used,tot_err/ndata);
  return max_err;
}


/*--------------------------------------------------------------------------------*/
void get_radec_from_altaz_exact_1det(const mbTOD *tod,int det,    PointingFitScratch *scratch)
//set alt and az, create the ctime, then call slalib to get ra/dec
{
  
  assert(tod);
  assert(tod->pointingOffset);
  actData mydalt=get_alt_offset(tod,det);
  actData mydaz=get_az_offset(tod,det);
  double *ctime=dvector(tod->ndata);
  //printf("Starting altaz\n");
  for (int i=0;i<tod->ndata;i++) {
    scratch->alt[i]=scratch->tod_alt[i]+mydalt;
#ifdef OLD_POINTING_OFFSET
    scratch->az[i]=scratch->tod_az[i]+mydaz/cos(scratch->alt[i]);
#else
    scratch->az[i]=scratch->tod_az[i]+mydaz/cos(scratch->tod_alt[i]);
#endif
    ctime[i]=tod->ctime+((double)i)*tod->deltat;
  }
  //printf("Finished altaz\n");
  Site site;
  ACTSite(&site);
  //printf("got site.\n");
  act_observed_altaz_to_mean_radec(&site,150.0,tod->ndata,ctime,scratch->alt,scratch->az,scratch->ra,scratch->dec);
  //printf("finished calculation.\n");
  free(ctime);

}


/*--------------------------------------------------------------------------------*/

void find_tod_radec_lims(mbTOD *tod) 
{
  assert(tod);
  assert(tod->pointing_fit);
  assert(tod->az);
  assert(tod->alt);
  actData azmin=vecmin(tod->az,tod->ndata);
  actData azmax=vecmax(tod->az,tod->ndata);
  actData altmin=vecmin(tod->alt,tod->ndata);
  actData altmax=vecmax(tod->alt,tod->ndata);
  
  actData *daz, *dalt;
  int ndet_good;
  get_uncut_daz_dalt(tod,&daz,&dalt,&ndet_good);
  
  actData *ra=vector(8*ndet_good);
  actData *dec=vector(8*ndet_good);
  actData *az=vector(8*ndet_good);
  actData *alt=vector(8*ndet_good);
  int *ind=(int *)malloc(sizeof(int)*8*ndet_good);
  int ii=0;
  
  actData az_corners[4]={azmin, azmin, azmax ,azmax};
  actData alt_corners[4]={altmin, altmax ,altmin, altmax};
  
  for (int i=0;i<ndet_good;i++) 
    for (int j=0;j<4;j++) {
      alt[ii]=alt_corners[j]+dalt[i];
#ifdef OLD_POINTING_OFFSET
      az[ii]=az_corners[j]+daz[i]/cos(alt[ii]);
#else
      az[ii]=az_corners[j]+daz[i]/cos(tod->alt[tod->ndata/2]);
#endif
      ind[ii]=0;
      ii++;
      alt[ii]=alt[ii-1];
      az[ii]=az[ii-1];
      ind[ii]=tod->ndata-1;
      ii++;
    }
  assert(ii==8*ndet_good); //sanity check.
  get_radec_from_altaz_fit_with_ind(tod,alt,az,ind,ii,ra,dec);
  tod->ramin=vecmin(ra,ii);
  tod->ramax=vecmax(ra,ii);
  tod->decmin=vecmin(dec,ii);
  tod->decmax=vecmax(dec,ii);
  
  free(ind);
  free(daz);
  free(dalt);
  free(az);
  free(alt);

}

/*--------------------------------------------------------------------------------*/

void set_tod_starting_altaz_ctime(mbTOD *tod, actData alt, actData az, double ctime)
//overwrite the starting locations with something else.  useful for doing things like
//simulations where you can just pull in a single TOD w/ its scan pattern and adjust.
//if the ctime is not positive, quit and go home.
{
  if (ctime<=0)
    return;
  assert(tod);
  assert(tod->alt);
  assert(tod->az);

  assert(tod->ndata>0);
  actData az0=tod->az[0];
  actData alt0=tod->alt[0];
  for (int i=0;i<tod->ndata;i++) {
    tod->alt[i]=alt+(tod->alt[i]-alt0);
    tod->az[i]=az+(tod->az[i]-az0);    
  }
  tod->ctime=ctime;

}

/*--------------------------------------------------------------------------------*/

void shift_tod_pointing(mbTOD *tod, actData ra_shift, actData dec_shift)
{
  assert(tod);
  assert(tod->pointing_fit);
  tod->pointing_fit->ra_offset+=ra_shift;
  tod->pointing_fit->dec_offset+=dec_shift;

}

/*--------------------------------------------------------------------------------*/

void destroy_pointing_fit_raw(PointingFit *fit)
{
  destroy_2d_poly(fit->ra_fit);
  destroy_2d_poly(fit->dec_fit);
  if (fit->ncoarse>0)
    free(fit->coarse_ind);
  free(fit);
}
/*--------------------------------------------------------------------------------*/

void destroy_pointing_fit(mbTOD *tod)
{
  if (tod->pointing_fit==NULL)
    return;
  destroy_pointing_fit_raw(tod->pointing_fit);
  //destroy_2d_poly(tod->pointing_fit->ra_fit);
  //destroy_2d_poly(tod->pointing_fit->dec_fit);
  //free(tod->pointing_fit);
  tod->pointing_fit=NULL;
}

