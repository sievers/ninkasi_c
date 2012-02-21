
#ifndef MAKEFILE_HAND
#include "config.h"
#endif

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <getopt.h>
#include <omp.h>
//#include <cpgplot.h>

#include "ninkasi.h"

#ifdef HAVE_MPI
#  include <mpi.h>
#  ifndef ACTDATA_DOUBLE
#    define MPI_NType MPI_FLOAT
#  else
#    define MPI_NType MPI_DOUBLE
#  endif
#endif

#include "dirfile.h"
#include "astro.h"
#include "mbCommon.h"
#include "mbCuts.h"
#include "mbTOD.h"
#include "noise.h"
#include "ninkasi_pointing.h"


#define ALTAZ_PER_LINE 3

/*================================================================================*/
 
 int main(int argc, char *argv[])
 {


#if 0
   Site site;
   site.latitude=0.670784476117809;
   site.east_longitude=-1.39346796591019;
   site.elevation_m=1000.0;
   site.temperature_K=290.0;
   site.pressure_mb=900.0;
   site.relative_humidity=0.5;
   struct timeval tv;
   gettimeofday(&tv,NULL);
   printf("ctime is %d\n",tv.tv_sec);
   double ctime_1, ctime_2;
   ctime_1=1223150033.0+75600;  //1PM EDT, Sunday October 5, 2008
   ctime_2=ctime_1+3600.0*3.5;
   //az=58, el=70 in degrees @start


   float *ra=vector(2);
   float *dec=vector(2);
   float *alt=vector(2);
   float *az=vector(2);
   double *tt2=dvector(2);
   tt2[0]=ctime_1;
   tt2[1]=tt2[0];
   alt[0]=mydeg2rad(atof(argv[1]));
   alt[1]=alt[0];
   az[0]=mydeg2rad(atof(argv[2]));
   az[1]=az[0]+mydeg2rad(78.0/3600/cos(alt[1]));
   //az[1]=deg2rad(az[0]+78.0/3600/cos(alt[1]));
   
   observed_altaz_to_mean_radec(&site,30.0,2,tt2,alt,az,ra,dec);
   printf("ra/dec are %14.6f %14.6f\n",ra[0]*180.0/3.14159265/15.0,dec[0]*180.0/3.14159265);
   printf("ra/dec are %14.6f %14.6f\n",ra[1]*180.0/3.14159265/15.0,dec[1]*180.0/3.14159265);
   printf("dra and ddec are %14.6f %14.6f\n",(ra[1]-ra[0])*180.0/3.14159265/15.0*3600.0,(dec[1]-dec[0])*180.0/3.14159265*3600.0);
   exit(EXIT_SUCCESS);
#endif   

   TODvec tods;
   MAPvec maps;
   PARAMS params;
   memset(&params,0,sizeof(PARAMS));
   
   
#ifdef HAVE_MPI
   
   int ierr=MPI_Init(&argc, &argv);
   int myrank, nproc;
  
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

  get_parameters(argc, argv, &params);
#ifdef HAVE_MPI
  if (myrank==0)
#endif
    print_options(&params);
  if (params.quit)
    exit(EXIT_SUCCESS);  
  
  
#if 0
  for (int i=0;i<1500;i++) {
    printf("i is %d\n",i);
    mbTOD *tmp=read_dirfile_tod_header(params.datanames[0]);
    printf("file %s has ndata and ndet %d %d\n",params.datanames[0],tmp->ndata,tmp->ndet);
    tmp->cuts=mbCutsAlloc(tmp->nrow,tmp->ncol);
    tmp->pointingOffset=nkReadPointingOffset(params.pointing_file);

    assign_tod_ra_dec(tmp);
    //destroy_pointing_fit(tmp);

    destroy_pointing_offset(tmp->pointingOffset);
    free(tmp->alt);
    free(tmp->az);

  }
  exit(EXIT_SUCCESS);
#endif

  tods.total_tod=how_many_tods(tods.froot,&params);
  mprintf(stdout,"have %d tods.\n",tods.total_tod);
  find_my_tods(&tods,&params);
  
  for (int i=0;i<tods.ntod;i++)
    mprintf(stdout,"I own %s\n",tods.my_fnames[i]);
  
  read_all_tod_headers(&tods,&params);
  set_global_radec_lims(&tods);
  mprintf(stdout,"global limits are %12.5f %12.5f %12.5f %12.5f\n",tods.ramin,tods.ramax,tods.decmin,tods.decmax);


#if 0
  mbTOD *tod=&(tods.tods[0]);
  
  pca_time tt;
  tick(&tt);
  
  find_pointing_pivots(tod,1.0);
  printf("have %d samples.\n",tod->pointing_fit->ncoarse);
  //exit(EXIT_SUCCESS);
#pragma omp parallel shared(tod) default(none)
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
#pragma omp for schedule(dynamic, 1)
    for (int i=0;i<tod->ndet;i++) 
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	get_radec_from_altaz_fit_1det_coarse(tod,i,scratch);
      }
    destroy_pointing_fit_scratch(scratch);
  }
  printf("Time to do pointing is %10.4f seconds.\n",tocksilent(&tt));
  
  printf("Max pointing eror is %14.3f\n",find_max_pointing_err(&(tods.tods[0]))*180*3600/M_PI);
  
  exit(EXIT_SUCCESS);
#endif

  maps.nmap=1;
  {
    MAP *mapvec=(MAP *)malloc(maps.nmap*sizeof(MAP));
    for (int i=0;i<maps.nmap;i++) {
      mapvec[i].projection=(nkProjection *)malloc(sizeof(nkProjection));
      mapvec[i].projection->proj_type=NK_RECT;
      mapvec[i].have_locks=0;
    }
    maps.maps=&mapvec;
  }

  maps.maps[0]->pixsize=params.pixsize;  //30 arcsecond pixels
	      
  maps.maps[0]->ramin=tods.ramin;
  maps.maps[0]->ramax=tods.ramax;
  maps.maps[0]->decmin=tods.decmin;
  maps.maps[0]->decmax=tods.decmax;
  setup_maps(&maps,&params);
  clear_mapset(&maps);
  createFFTWplans(&tods);


  run_PCG(&maps,&tods,&params);
  readwrite_simple_map(maps.maps[0],params.outname,DOWRITE);

  exit(EXIT_SUCCESS);
  run_PCG(&maps,&tods,&params);
  readwrite_simple_map(maps.maps[0],"test_crap.map",DOWRITE);  

    
}
