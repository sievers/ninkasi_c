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
#ifndef NO_FFTW
#include <fftw3.h>
#endif

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
#include "readtod.h"
#include "astro.h"
#include "mbCommon.h"
#include "mbCuts.h"
#include "mbTOD.h"
#include "noise.h"
#include "ninkasi_pointing.h"
#include "nk_clapack.h"
#include "ninkasi_mathutils.h"

#define ALTAZ_PER_LINE 3
/*--------------------------------------------------------------------------------*/
void pca_pause(actData pauselen)
{
  pca_time tt;
  tick(&tt);
  while (tocksilent(&tt)<pauselen)
    ;

}
/*--------------------------------------------------------------------------------*/

void write_timestream(actData *vec, int n, char *fname) 
{
  FILE *outfile=fopen(fname,"w");
  if (!outfile) {
    fprintf(stderr,"unable to open %s for writing in write_timestream.\n",fname);
    return;
  }
  fwrite(&n,1,sizeof(int),outfile);
  fwrite(vec,n,sizeof(actData),outfile);
  fclose(outfile);
}
/*--------------------------------------------------------------------------------*/

void write_double_timestream(double *vec, int n, char *fname) 
{
  FILE *outfile=fopen(fname,"w");
  if (!outfile) {
    fprintf(stderr,"unable to open %s for writing in write_timestream.\n",fname);
    return;
  }
  n*=-1;
  fwrite(&n,1,sizeof(int),outfile);
  n*=-1; 
  fwrite(vec,n,sizeof(double),outfile);
  fclose(outfile);

}

/*--------------------------------------------------------------------------------*/
void mprintf(FILE *stream, char *format, ...)
/*In an mpi world, only the master will print stuff out.*/
{
#ifdef HAVE_MPI
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#else
  int myid=0;
#endif
  if (myid==0) {
    va_list args;
      va_start(args,format);
      vfprintf(stream,format,args);
      va_end(args);
      
  }
}

/*--------------------------------------------------------------------------------*/

void tick(pca_time *tt)
{
  gettimeofday(tt,NULL);
}

/*--------------------------------------------------------------------------------*/


actData tocksilent(pca_time *tt)
{
  pca_time tnow;
  gettimeofday(&tnow,NULL);
  return (actData)((tnow.tv_usec-tt->tv_usec)/1.0e6+(tnow.tv_sec-tt->tv_sec));  
  
}
/*--------------------------------------------------------------------------------*/
actData mbElapsedTime(pca_time *tt)
{
  return (actData)tocksilent(tt);
}
/*--------------------------------------------------------------------------------*/
void mbStartTime(pca_time *tt)
{
  tick(tt);
}

/*--------------------------------------------------------------------------------*/
void *malloc_retry(size_t n)
{
  int retry=10; //try this many times
  actData pause_len=0.1;  //wait this long if we fail
  for (int i=0;i<retry;i++) {
    void *vec=malloc(n);
    if (vec)
      return vec;
    fprintf(stderr,"Malloc failure when asking for %ld bytes..  Retrying %d...\n",n,i);
    pca_pause(pause_len);
  }
  return NULL;
}
/*--------------------------------------------------------------------------------*/

actData *vector(long n)
{
  actData *data;
  data=(actData *)malloc_retry(sizeof(actData)*n);
  assert(data!=NULL);
   
   return data;  
}
/*--------------------------------------------------------------------------------*/

int *ivector(long n)
{
  int *data;
  data=(int *)malloc_retry(sizeof(int)*n);
  assert(data!=NULL);
   
   return data;  
}
/*--------------------------------------------------------------------------------*/

float *svector(long n)
{
  float *data;
  data=(float *)malloc_retry(sizeof(float)*n);
  assert(data!=NULL);
   
   return data;  
}

/*--------------------------------------------------------------------------------*/

double *dvector(long n)
{
  double *data;
  data=(double *)malloc_retry(sizeof(double)*n);
  assert(data!=NULL);
   
   return data;  
}

/*--------------------------------------------------------------------------------*/

actComplex *cvector(long n)
{
  actComplex *data;
   data=(actComplex *)malloc_retry(sizeof(actComplex)*n);
   assert(data!=NULL);
   
   return data;  
}

/*--------------------------------------------------------------------------------*/

actData **matrix(long n,long m)
{
  actData *data, **ptrvec;
  data=(actData *)malloc_retry(sizeof(actData)*n*m);
  assert(data!=NULL);

  ptrvec=(actData **)malloc_retry(sizeof(actData *)*n);
  assert(ptrvec!=NULL);
  
  for (int i=0;i<n;i++)
    ptrvec[i]=&data[i*m];
   
   return ptrvec;
}
/*--------------------------------------------------------------------------------*/

int **imatrix(long n,long m)
{
  int *data, **ptrvec;
  data=(int *)malloc_retry(sizeof(int)*n*m);
  assert(data!=NULL);

  ptrvec=(int **)malloc_retry(sizeof(int *)*n);
  assert(ptrvec!=NULL);
  
  for (int i=0;i<n;i++)
    ptrvec[i]=&data[i*m];
   
   return ptrvec;
}
/*--------------------------------------------------------------------------------*/

actComplex **cmatrix(long n,long m)
{
  actComplex *data, **ptrvec;
  data=(actComplex *)malloc_retry(sizeof(actComplex)*n*m);
  assert(data!=NULL);
  ptrvec=(actComplex **)malloc_retry(sizeof(actComplex *)*n);
   assert(ptrvec!=NULL);
   for (int i=0;i<n;i++)
     ptrvec[i]=&data[i*m];
   
   return ptrvec;
}
/*--------------------------------------------------------------------------------*/

float **smatrix(long n,long m)
{
  float *data, **ptrvec;
  data=(float *)malloc_retry(sizeof(float)*n*m);
  assert(data!=NULL);
  ptrvec=(float **)malloc_retry(sizeof(float *)*n);
   assert(ptrvec!=NULL);
   for (int i=0;i<n;i++)
     ptrvec[i]=&data[i*m];
   
   return ptrvec;
}
/*--------------------------------------------------------------------------------*/
void free_matrix(actData **mat)
{
  free(mat[0]);
  free(mat);
}
/*--------------------------------------------------------------------------------*/
void free_smatrix(float **mat)
{
  free(mat[0]);
  free(mat);
}

/*--------------------------------------------------------------------------------*/

double **dmatrix(long n,long m)
{
  double *data, **ptrvec;
  data=(double *)malloc_retry(sizeof(double)*n*m);
  assert(data!=NULL);
  ptrvec=(double **)malloc_retry(sizeof(double *)*n);
   assert(ptrvec!=NULL);
   for (int i=0;i<n;i++)
     ptrvec[i]=&data[i*m];
   
   return ptrvec;
}


/*--------------------------------------------------------------------------------*/
void act_fftw_free(act_fftw_complex *vec)
{
#ifndef ACTDATA_DOUBLE
  fftwf_free(vec);
#else
  fftw_free(vec);
#endif
}

/*--------------------------------------------------------------------------------*/
act_fftw_complex *act_fftw_malloc(size_t n)
{
#ifndef ACTDATA_DOUBLE
  return fftwf_malloc(n);
#else
  return fftw_malloc(n);
#endif
}

/*--------------------------------------------------------------------------------*/
void act_fftw_destroy_plan(act_fftw_plan p)
{
#ifndef ACTDATA_DOUBLE
  fftwf_destroy_plan(p);
#else
  fftw_destroy_plan(p);
#endif
}

/*--------------------------------------------------------------------------------*/
void act_fftw_execute_dft_r2c(act_fftw_plan p, actData *r, act_fftw_complex *c)
{
#ifndef ACTDATA_DOUBLE
  fftwf_execute_dft_r2c(p,r,c);
#else
  fftw_execute_dft_r2c(p,r,c);
#endif

}
/*--------------------------------------------------------------------------------*/
void act_fftw_execute_dft_c2r(act_fftw_plan p,  act_fftw_complex *c, actData *r)
{
#ifndef ACTDATA_DOUBLE
  fftwf_execute_dft_c2r(p,c,r);
#else
  fftw_execute_dft_c2r(p,c,r);
#endif

}
/*--------------------------------------------------------------------------------*/
act_fftw_plan act_fftw_plan_dft_r2c_1d(int n, actData *vec, act_fftw_complex *vec2,unsigned flags)
{
#ifndef ACTDATA_DOUBLE  
  act_fftw_plan p=fftwf_plan_dft_r2c_1d(n,vec,vec2,flags);
#else
  act_fftw_plan p=fftw_plan_dft_r2c_1d(n,vec,vec2,flags);
#endif
  return p;
}
/*--------------------------------------------------------------------------------*/

act_fftw_plan act_fftw_plan_dft_c2r_1d(int n, act_fftw_complex *vec2,actData *vec, unsigned flags)
{
#ifndef ACTDATA_DOUBLE  
  act_fftw_plan p=fftwf_plan_dft_c2r_1d(n,vec2,vec,flags);
#else
  act_fftw_plan p=fftw_plan_dft_c2r_1d(n,vec2,vec,flags);
#endif
  return p;
}
			
/*--------------------------------------------------------------------------------*/
FILE *fopen_safe(char *filename, char *mode)
//retry a few times if we have a problem.
{
  int maxtry=10;
  actData dt=0.5;  //try maxtry times, with a spacing of dt
  pca_time tt;

  for (int i=0;i<maxtry;i++)
    {
      tick(&tt);
      FILE *file=fopen(filename,mode);
      if (file)
	return file;
      actData elapsed=tocksilent(&tt);
      while (elapsed<dt)
	elapsed=tocksilent(&tt);
    }
  //if we got here, it's 'cause we couldn't open a file.
  return NULL;
}

/*--------------------------------------------------------------------------------*/

int how_many_tods(char *froot, PARAMS *params)
{
  if (strlen(params->altaz_file)) {
    int myargc;
    char *line_in=read_all_stdin(params->altaz_file);
    char **myargv=create_argv_new(line_in,&myargc," \n");
    int nfiles=myargc/ALTAZ_PER_LINE;
    if (myargc>params->ntod) {
      for (int i=params->ntod;i<nfiles;i++) {
	strncpy(params->datanames[i],params->datanames[i-1],MAXLEN-1);	
      }
      params->ntod=nfiles;
      free(line_in);
      free_argv(myargc,myargv);
    }
    
  }
  return params->ntod;
}

/*--------------------------------------------------------------------------------*/
int find_my_tods(TODvec *tods, PARAMS *params)
{
#ifdef HAVE_MPI
  int ierr,myid,nproc;
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  assert(ierr==0);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  assert(ierr==0);
#else
  int myid=0;  //update these when ready
  int nproc=1;
#endif

  tods->ntod=0;


  //first, do a simple loop to figure out how many I own
  for (int i=myid;i<tods->total_tod;i+=nproc) {
    tods->ntod++;
  }
  tods->my_fnames=(char **)malloc_retry(sizeof(char *)*tods->ntod);
  tods->tods=(mbTOD *)calloc(sizeof(mbTOD),tods->ntod);
  int ii=0;
  for (int i=myid;i<tods->total_tod;i+=nproc) {
    tods->tods[ii].seed=((1+fabs(params->seed))*MAXTOD+i)*MAXDET;
    tods->my_fnames[ii]=strdup(params->datanames[i]);
    
    ii++;

  }
  return 0;
  
}
/*--------------------------------------------------------------------------------*/
actData *read_1d_datafile(char *fname, int *n)
//read a data file that has n elements of actData pass in 0 in *n to assign the value
//
{
  FILE *infile=fopen_safe(fname,"r");
  if (!infile)
    return NULL;
  int n_check;
  size_t crap; //just pulling the values of fread so we don't get warnings.
  crap=fread(&n_check,1,sizeof(int),infile);
  if (*n!=0) 
    assert(*n==n_check);
  else
    *n=n_check;
  actData *vec=vector(*n);
  crap=fread(vec,*n,sizeof(actData),infile);
  fclose(infile);
  return vec;
}

/*--------------------------------------------------------------------------------*/
double *read_1d_datafile_double(char *fname, int *n)
//read a data file that has n elements of actData pass in 0 in *n to assign the value
//
{
  FILE *infile=fopen_safe(fname,"r");
  if (!infile)
    return NULL;
  int n_check;
  size_t crap;  //pull return value of fread just to avoid a warning
  crap=fread(&n_check,1,sizeof(int),infile);
  if (*n!=0) 
    assert(*n==n_check);
  else
    *n=n_check;
  double *vec=(double *)malloc_retry(sizeof(double)*(*n));
  crap=fread(vec,*n,sizeof(double),infile);
  fclose(infile);
  return vec;
}
/*--------------------------------------------------------------------------------*/

void purge_vector(mbTOD *tod, actData *vec, int ngood)
{
  //printf("shrinking from %d to %d\n",tod->ndet,ngood);
  actData *tmp=vector(ngood);
  int icur=0;
  for (int i=0;i<tod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
      tmp[icur]=vec[i];
      icur++;
    }
  }
  actData *crud=(actData *)realloc(vec,sizeof(actData)*ngood);
  assert(crud==vec);
  memcpy(vec,tmp,sizeof(actData)*ngood);
  free(tmp);
  
  //free(vec);
  //*vec_in=tmp;
}
/*--------------------------------------------------------------------------------*/
void purge_cut_detectors(mbTOD *tod)
{
  int ngood=0;
  for (int i=0;i<tod->ndet;i++)
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))
      ngood++;

#ifdef ACTPOL
  if (tod->actpol_pointing) {
    ACTpolPointingFit *fit=tod->actpol_pointing;
    if (fit->dx) 
      purge_vector(tod,(tod->actpol_pointing->dx),ngood);
    if (fit->dy)
      purge_vector(tod,(tod->actpol_pointing->dy),ngood);
    if (fit->theta)
      purge_vector(tod,(tod->actpol_pointing->theta),ngood);
  }
#endif

  int *rows=(int *)malloc_retry(sizeof(int)*ngood);
  int *cols=(int *)malloc_retry(sizeof(int)*ngood);
  int icur=0;
  for (int i=0;i<tod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
      rows[icur]=tod->rows[i];
      cols[icur]=tod->cols[i];
      if (!(tod->dets==NULL))
	if (!(tod->dets[tod->rows[i]]==NULL))
	  tod->dets[tod->rows[i]][tod->cols[i]]=icur;
      icur++;
    }
    else {
      if (!(tod->dets==NULL))
	if (!(tod->dets[tod->rows[i]]==NULL))
	  tod->dets[tod->rows[i]][tod->cols[i]]=ACT_NO_VALUE;
    }
    
  }
  assert(icur==ngood);  //simple sanity check
  free(tod->rows);
  free(tod->cols);
  tod->rows=rows;
  tod->cols=cols;
  tod->ndet=ngood;

}

/*--------------------------------------------------------------------------------*/
int get_numel_cut(mbTOD *tod) //return how many elements have been cut in a TOD
{
  if (tod->cuts_fit_params) {
    int dd=tod->ndet-1;
    while (tod->cuts_fit_params[tod->rows[dd]][tod->cols[dd]]->nregions==0) {
      dd--;
      if (dd<0)
	assert(1==0);
    }
    mbCutFitParams *params=tod->cuts_fit_params[tod->rows[dd]][tod->cols[dd]];
    return params->starting_param[params->nregions-1]+params->nparams[params->nregions-1];
  }
  if (tod->cuts_as_vec) {
    int det=tod->ndet-1;
    while (1) {
      int row=tod->rows[det];
      int col=tod->cols[det];
      mbUncut *ptr=tod->cuts_as_vec[row][col];
      if (ptr)
	if (ptr->nregions>0)
	  return ptr->indexLast[ptr->nregions-1];
      det--;
      if (det==0)
	return 0;
    
    }
  }
  return 0;
    
}

/*--------------------------------------------------------------------------------*/
void setup_cutsfits_precon(mbTOD *tod)
{
  assert(tod);
  assert(tod->cuts_fit_params);
  actData *winvec=(actData *)malloc(sizeof(actData)*tod->ndata);
  for (int i=0;i<tod->ndata;i++)
    winvec[i]=1.0;
  int nsamp=tod->n_to_window;
  
  //#define CUTFITS_DEBUG

#if 0
  if (nsamp>0) {
    for (int i=0;i<nsamp;i++) {
      winvec[i]=0.5-0.5*cos(M_PI*i/(nsamp+0.0));
      winvec[i]=winvec[i]*winvec[i];
      //window2[i]=0.5+0.5*cos(M_PI*i/(nsamp+0.0));
      winvec[tod->ndata-i-1]=winvec[i];
    }    
  }
#endif

#pragma omp parallel for shared(tod,winvec,nsamp) default(none)
  for (int det=0;det<tod->ndet;det++) {
    mbCutFitParams *params=tod->cuts_fit_params[tod->rows[det]][tod->cols[det]];
    mbUncut *cut=tod->cuts_as_uncuts[tod->rows[det]][tod->cols[det]];
    params->precon=(actData ***)malloc(sizeof(actData **)*params->nregions);
    for (int i=0;i<params->nregions;i++) {
      int nelem=cut->indexLast[i]-cut->indexFirst[i];
      actData **mat=legendre_mat(nelem,params->nparams[i]);
#ifdef CUTFITS_DEBUG
      for (int ii=0;ii<params->nparams[i];ii++)
	for (int jj=0;jj<nelem;jj++) 
	  if (!isfinite(mat[ii][jj])) {
	    printf("Have a not-finite element on matrix element %d %d, segment  %d on detector %d %d with length %d\n",ii,jj,i,tod->rows[det],tod->cols[det],nelem);
	    return;
	  }
#endif

#ifdef CUTFITS_DEBUG
      if (i==2) {
	FILE *outfile=fopen("mat_temp.txt","w");
	for (int ii=0;ii<nelem;ii++) {
	  for (int jj=0;jj<params->nparams[i];jj++)
	    fprintf(outfile,"%13.5e ",mat[jj][ii]);
	  fprintf(outfile,"  %13.5e\n",winvec[ii+cut->indexFirst[i]]);
	}
	fclose(outfile);
	return;
      }
#endif

      actData **precon=matrix(params->nparams[i],params->nparams[i]);
      for (int ii=0;ii<params->nparams[i];ii++)
	for (int jj=ii;jj<params->nparams[i];jj++) {
	  precon[ii][jj]=0;
	  for (int kk=0;kk<nelem;kk++)
	    precon[ii][jj]+=mat[ii][kk]*mat[jj][kk]*winvec[kk+cut->indexFirst[i]];
	  precon[jj][ii]=precon[ii][jj];
	}
      
#ifdef CUTFITS_DEBUG
      for (int ii=0;ii<params->nparams[i];ii++)
	for (int jj=0;jj<params->nparams[i];jj++) 
	  if (!isfinite(precon[ii][jj])) {
	    printf("Have a not-finite precon element on matrix element %d %d, segment  %d on detector %d %d with length %d\n",ii,jj,i,tod->rows[det],tod->cols[det],nelem);
	    return;
	  }

      for (int ii=0;ii<params->nparams[i];ii++) {
	for (int jj=0;jj<params->nparams[i];jj++)
	  printf("%13.5e ",precon[ii][jj]);
	printf("\n");
      }
#endif
      int info=invert_posdef_mat(precon,params->nparams[i]);
      if (info) {
	printf("error inverting segment %d on detector %d %d %d\n",i,det,tod->rows[det],tod->cols[det]);
	//return;
      }
      for (int ii=0;ii<params->nparams[i];ii++)
	for (int jj=0;jj<params->nparams[i];jj++) 
	  if (!isfinite(precon[ii][jj])) {
	    printf("Have a not-finite element on segment %d on detector %d %d with length %d\n",i,tod->rows[det],tod->cols[det],nelem);
	    //return;
	  }
      params->precon[i]=precon;
      free(mat[0]);
      free(mat);
    }
  }
  free(winvec);
}
/*--------------------------------------------------------------------------------*/
void apply_cutfits_precon(mbTOD *tod, actData *params_in, actData *params_out)
{
  assert(tod);
  assert(tod->cuts_fit_params);

#pragma omp parallel for shared(tod,params_in,params_out) default(none) 
  for (int det=0;det<tod->ndet;det++) {
    mbCutFitParams *params=tod->cuts_fit_params[tod->rows[det]][tod->cols[det]];
    for (int i=0;i<params->nregions;i++) {
      for (int ii=0;ii<params->nparams[i];ii++) 
	for (int jj=0;jj<params->nparams[i];jj++) 
	  params_out[params->starting_param[i]+jj]+=params_in[params->starting_param[i]+jj]*params->precon[i][ii][jj];
    }
  }
  
}
/*--------------------------------------------------------------------------------*/
mbUncut ***get_cut_regions_global_index(mbTOD *tod) 
//make the globally referenced cuts indices, for a 1-d array for gapfilling.

{
  assert(tod);
  assert(tod->cuts_as_uncuts);
  
  mbUncut **vec=(mbUncut **)malloc_retry(sizeof(mbUncut *)*tod->nrow*tod->ncol);
  mbUncut ***mat=(mbUncut ***)malloc_retry(sizeof(mbUncut **)*tod->nrow);
  memset(vec,0,sizeof(mbUncut**)*tod->nrow*tod->ncol);
  for (int i=0;i<tod->nrow;i++)
    mat[i]=vec+i*tod->ncol;
  int **detsums=imatrix(tod->nrow,tod->ncol);
  memset(detsums[0],0,sizeof(int)*tod->nrow*tod->ncol);



  //sum up total # of cuts for each detector
  for (int det=0;det<tod->ndet;det++) {
    int row=tod->rows[det];
    int col=tod->cols[det];
    if (!mbCutsIsAlwaysCut(tod->cuts,row,col)) {
      mbUncut *cut=tod->cuts_as_uncuts[row][col];
      int tot=0;
      for (int i=0;i<cut->nregions;i++)
	tot+=(cut->indexLast[i]-cut->indexFirst[i]);
      detsums[row][col]=tot;
    }
  }

  //figure out total used to-date, shifting over one
  int tmp=0;
  for (int det=0;det<tod->ndet;det++) {
    int row=tod->rows[det];
    int col=tod->cols[det];
    int tmp2=detsums[row][col];
    detsums[row][col]=tmp;
    tmp+=tmp2;
  }
  
  
  
  //for (int row=0;row<tod->nrow;row++) 
  // for (int col=0;col<tod->ncol;col++) 
  for (int det=0;det<tod->ndet;det++) {
    int row=tod->rows[det];
    int col=tod->cols[det];
    if (!mbCutsIsAlwaysCut(tod->cuts,row,col))
      if (tod->cuts_as_uncuts[row][col]) {   //if the cuts are null, assuming nothing is cut.
	//printf("working on %d %d\n",row,col);

	
	mbUncut *todcut=tod->cuts_as_uncuts[row][col];

	mbUncut *cut=(mbUncut *)malloc(sizeof(mbUncut));
	cut->nregions=todcut->nregions;
	cut->indexFirst=(int *)malloc(sizeof(int)*todcut->nregions);
	cut->indexLast=(int *)malloc(sizeof(int)*todcut->nregions);
	mat[row][col]=cut;

	
	int cumsum=0;
	for (int i=0;i<todcut->nregions;i++) {
	  //printf("region is %d\n",i);
	  //printf("limits are %d %d\n",todcut->indexFirst[i],todcut->indexLast[i]);
	  cut->indexFirst[i]=cumsum+detsums[row][col];
	  cut->indexLast[i]=cut->indexFirst[i]+(todcut->indexLast[i]-todcut->indexFirst[i]);
	  //cut->indexLast[i]=cumsum+(todcut->indexLast[i]-todcut->indexFirst[i]);
	  cumsum+=cut->indexLast[i]-cut->indexFirst[i];
	}	
      }
  }
  return mat;
}
/*--------------------------------------------------------------------------------*/

void set_global_indexed_cuts(mbTOD *tod)
{
  assert(tod);
  assert(tod->cuts);
  if (!tod->uncuts) {
    //fprintf(stderr,"Setting uncuts in set_global_indexed_cuts.  Maybe not what you want to do.\n");
    get_tod_uncut_regions(tod);
  }
  if (!tod->cuts_as_uncuts)
    get_tod_cut_regions(tod);
  tod->cuts_as_vec=get_cut_regions_global_index(tod);
}
/*--------------------------------------------------------------------------------*/
mbUncut ***get_cut_regions(mbTOD *tod) 
//calculate the cut regions a la uncut regions and return a matrix to 'em.
//preferred form for gap-filling
{
  assert(tod);
  assert(tod->uncuts);
  
  mbUncut **vec=(mbUncut **)malloc_retry(sizeof(mbUncut *)*tod->nrow*tod->ncol);
  mbUncut ***mat=(mbUncut ***)malloc_retry(sizeof(mbUncut **)*tod->nrow);
  memset(vec,0,sizeof(mbUncut**)*tod->nrow*tod->ncol);

  for (int i=0;i<tod->nrow;i++)
    mat[i]=vec+i*tod->ncol;

  for (int i=0;i<tod->nrow;i++)
    for (int j=0;j<tod->ncol;j++) {
      //printf("working on detector %d %d\n",i,j);
      if (tod->uncuts[i][j]==NULL)
	printf("uncuts is null.\n");
      mat[i][j]=mbCutsInvertUncut(tod->uncuts[i][j],tod->ndata);
    }
  return mat;
}
/*--------------------------------------------------------------------------------*/
mbCutFitParams ***setup_cut_fit_params(mbTOD *tod, int *nparams_from_length)
//nparams_from_length contains the mapping from cut length to # of polynomial parameters
{
  assert(tod);
  assert(tod->cuts_as_uncuts);

  mbCutFitParams **vec=(mbCutFitParams **)malloc_retry(sizeof(mbCutFitParams *)*tod->nrow*tod->ncol);
  mbCutFitParams ***mat=(mbCutFitParams ***)malloc_retry(sizeof(mbCutFitParams **)*tod->nrow);
  memset(vec,0,sizeof(mbCutFitParams **)*tod->nrow*tod->ncol);

  for (int i=0;i<tod->nrow;i++)
    mat[i]=vec+i*tod->ncol;

  int global_accum=0;
#if 0
  for (int i=0;i<tod->nrow;i++)
    for (int j=0;j<tod->ncol;j++) {
      //printf("working on detector %d %d\n",i,j);                                                                                                                                           
      if (tod->cuts_as_uncuts[i][j]==NULL)
        printf("cuts_as_uncuts is null.\n");
      mat[i][j]=(mbCutFitParams *)malloc_retry(sizeof(mbCutFitParams));
      mat[i][j]->nregions=tod->cuts_as_uncuts[i][j]->nregions;
      mat[i][j]->starting_param=(int *)malloc(mat[i][j]->nregions*sizeof(int));
      mat[i][j]->nparams=(int *)malloc(mat[i][j]->nregions*sizeof(int));
      for (int k=0;k<tod->cuts_as_uncuts[i][j]->nregions;k++) {
	mat[i][j]->nparams[k]=nparams_from_length[tod->cuts_as_uncuts[i][j]->indexLast[k]-tod->cuts_as_uncuts[i][j]->indexFirst[k]];
	mat[i][j]->starting_param[k]=global_accum;
	global_accum+=mat[i][j]->nparams[k];
      }
    }
#else
  for (int det=0;det<tod->ndet;det++) {
    int i=tod->rows[det];
    int j=tod->cols[det];
    //printf("working on %4d %2d %2d with accum %8d\n",det,i,j,global_accum);
    if (tod->cuts_as_uncuts[i][j]==NULL)
      printf("cuts_as_uncuts is null.\n");
    mat[i][j]=(mbCutFitParams *)malloc_retry(sizeof(mbCutFitParams));
    mat[i][j]->nregions=tod->cuts_as_uncuts[i][j]->nregions;
    mat[i][j]->starting_param=(int *)malloc(mat[i][j]->nregions*sizeof(int));
    mat[i][j]->nparams=(int *)malloc(mat[i][j]->nregions*sizeof(int));
    for (int k=0;k<tod->cuts_as_uncuts[i][j]->nregions;k++) {
      int cutlen=tod->cuts_as_uncuts[i][j]->indexLast[k]-tod->cuts_as_uncuts[i][j]->indexFirst[k];
      //printf("cutlen is %5d, nparams is %d\n",cutlen,nparams_from_length[cutlen]);
      mat[i][j]->nparams[k]=nparams_from_length[tod->cuts_as_uncuts[i][j]->indexLast[k]-tod->cuts_as_uncuts[i][j]->indexFirst[k]];
      mat[i][j]->starting_param[k]=global_accum;
      global_accum+=mat[i][j]->nparams[k];
    }
  }
  
#endif
  return mat;
  
}
/*--------------------------------------------------------------------------------*/
int cutpolys2tod(mbTOD *tod,actData *cutvec)
{
  if (!tod->have_data) {
    fprintf(stderr,"missing data in tod in cutvec2tod\n");
    return 1;
  }
  if (tod->cuts_as_vec==NULL) {
    fprintf(stderr,"Error, tod does not contain cuts_as_vec in cutpolys2tod.\n");
    return 1;
  }
  if (tod->cuts_as_uncuts==NULL) {
    fprintf(stderr,"Error, tod does not contain cuts_as_uncut in cutpolys2tod.\n");
    return 1;    
  }
  if (tod->cuts_fit_params==NULL) {
    fprintf(stderr,"Error - tod does not contain cut fit parameters in cutpolys2tod.\n");
    return 1;
  }
#pragma omp parallel for shared(tod,cutvec) default(none)
  for (int det=0;det<tod->ndet;det++) {
    mbCutFitParams *params=tod->cuts_fit_params[tod->rows[det]][tod->cols[det]];
    mbUncut *cut=tod->cuts_as_uncuts[tod->rows[det]][tod->cols[det]];
    for (int i=0;i<params->nregions;i++) {
      legendre_eval(tod->data[det]+cut->indexFirst[i],cut->indexLast[i]-cut->indexFirst[i],cutvec+params->starting_param[i],params->nparams[i]);
    }
  }
  
  return 0;
}
/*--------------------------------------------------------------------------------*/
int tod2cutpolys(mbTOD *tod,actData *cutvec)
{
  if (!tod->have_data) {
    fprintf(stderr,"missing data in tod in cutvec2tod\n");
    return 1;
  }
  if (tod->cuts_as_vec==NULL) {
    fprintf(stderr,"Error, tod does not contain cuts_as_vec in cutpolys2tod.\n");
    return 1;
  }
  if (tod->cuts_as_uncuts==NULL) {
    fprintf(stderr,"Error, tod does not contain cuts_as_uncut in cutpolys2tod.\n");
    return 1;    
  }
  if (tod->cuts_fit_params==NULL) {
    fprintf(stderr,"Error - tod does not contain cut fit parameters in cutpolys2tod.\n");
    return 1;
  }
#pragma omp parallel for shared(tod,cutvec) default(none)
  for (int det=0;det<tod->ndet;det++) {
    mbCutFitParams *params=tod->cuts_fit_params[tod->rows[det]][tod->cols[det]];
    mbUncut *cut=tod->cuts_as_uncuts[tod->rows[det]][tod->cols[det]];
    for (int i=0;i<params->nregions;i++) {
      legendre_project(tod->data[det]+cut->indexFirst[i],cut->indexLast[i]-cut->indexFirst[i],cutvec+params->starting_param[i],params->nparams[i]);
    }
  }
  
  return 0;
}
/*--------------------------------------------------------------------------------*/
mbUncut ***get_uncut_regions(mbTOD *tod) 
//calculate the uncut regions and return a matrix to 'em.
{
  mbUncut **vec=(mbUncut **)malloc_retry(sizeof(mbUncut *)*tod->nrow*tod->ncol);
  mbUncut ***mat=(mbUncut ***)malloc_retry(sizeof(mbUncut **)*tod->nrow);
  for (int i=0;i<tod->nrow;i++)
    mat[i]=vec+i*tod->ncol;
  //tod->uncuts=mat;

  //int fac=1;
  //for (int i=0;i<tod->decimate;i++)
  //fac*=2;
  for (int i=0;i<tod->nrow;i++)
    for (int j=0;j<tod->ncol;j++)
      mat[i][j]=mbCutsGetUncut(tod->cuts,i,j,0,tod->ndata);
  //decimate_uncut_regions(tod);
  return mat;
}
/*--------------------------------------------------------------------------------*/
void get_tod_cut_regions(mbTOD *tod)
{
  if (tod->uncuts==NULL)
    get_tod_uncut_regions(tod);
  tod->cuts_as_uncuts=get_cut_regions(tod);
}


/*--------------------------------------------------------------------------------*/
void get_tod_uncut_regions(mbTOD *tod)
{
#if 1
  tod->uncuts=get_uncut_regions(tod);
#else  
  mbUncut **vec=(mbUncut **)malloc_retry(sizeof(mbUncut *)*tod->nrow*tod->ncol);
  mbUncut ***mat=(mbUncut ***)malloc_retry(sizeof(mbUncut **)*tod->nrow);
  for (int i=0;i<tod->nrow;i++)
    mat[i]=vec+i*tod->ncol;
  tod->uncuts=mat;

  //int fac=1;
  //for (int i=0;i<tod->decimate;i++)
  //fac*=2;
  for (int i=0;i<tod->nrow;i++)
    for (int j=0;j<tod->ncol;j++)
      tod->uncuts[i][j]=mbCutsGetUncut(tod->cuts,i,j,0,tod->ndata);
  //decimate_uncut_regions(tod);
#endif
}
/*--------------------------------------------------------------------------------*/
void reverse_tod_uncut_regions(mbTOD *tod)
//time-reflect the uncut regions.  So if the first ten seconds were uncut, make it so the last
//ten seconds are uncut.
{
  for (int i=0;i<tod->nrow;i++)
    for (int j=0;j<tod->nrow;j++) 
      if (!mbCutsIsAlwaysCut(tod->cuts,i,j))
	if (tod->uncuts[i][j]->nregions>0) {
	  //printf("reversing %2d %2d\n",i,j);
	  int nreg=tod->uncuts[i][j]->nregions;
	  int *new_first=(int *)malloc_retry(sizeof(int)*nreg);
	  int *new_last=(int *)malloc_retry(sizeof(int)*nreg);
	  
	  for (int reg=0;reg<nreg;reg++) {
	    new_first[reg]=tod->ndata-tod->uncuts[i][j]->indexLast[nreg-reg-1];
	    new_last[reg]=tod->ndata-tod->uncuts[i][j]->indexFirst[nreg-reg-1];	  
	    //new_first[reg]=tod->uncuts[i][j]->indexLast[nreg-reg-1];
	    //new_last[reg]=tod->uncuts[i][j]->indexFirst[nreg-reg-1];	  
	  }
	  memcpy(tod->uncuts[i][j]->indexFirst,new_first,sizeof(int)*nreg);
	  memcpy(tod->uncuts[i][j]->indexLast,new_last,sizeof(int)*nreg);
	  

	  free(new_first);
	  free(new_last);
	  //printf("finished.\n");
	}
  
  
}
/*--------------------------------------------------------------------------------*/
void decimate_uncut_regions(mbTOD *tod)
{
  for (int ii=0;ii<tod->decimate;ii++) {
    for (int row=0;row<tod->nrow;row++) 
      for (int col=0;col<tod->ncol;col++) {
	if (tod->uncuts[row][col]->nregions>0)
	  for (int i=0;i<tod->uncuts[row][col]->nregions;i++) {
	    tod->uncuts[row][col]->indexFirst[i]=(tod->uncuts[row][col]->indexFirst[i]+1)/2;
	    tod->uncuts[row][col]->indexLast[i]=(tod->uncuts[row][col]->indexLast[i])/2;
	  }
      }
  }
}
/*--------------------------------------------------------------------------------*/
void fill_gaps_stupid(mbTOD *tod)
{
  
  if (!(tod->have_data)) {
    fprintf(stderr,"Warning - attempting to gapfill non-existent data.\n");
    return;
  }
  if (tod->uncuts==NULL) {
    fprintf(stderr,"Warning - attempting to gapfill without uncuts.\n");
    return;
  }
#pragma omp parallel for shared(tod) default(none)
  for (int det=0;det<tod->ndet;det++) {
    mbUncut *uncut;
    uncut=tod->uncuts[tod->rows[det]][tod->cols[det]];
    if (tod->uncuts_for_interp)
      uncut=tod->uncuts_for_interp[tod->rows[det]][tod->cols[det]];
    int nreg=uncut->nregions;
    if (nreg>0) {
      if (uncut->indexFirst[0]>0)
	for (int i=0;i<uncut->indexFirst[0];i++)
	  tod->data[det][i]=tod->data[det][uncut->indexFirst[0]];
      if (uncut->indexLast[nreg-1]<tod->ndata)
	for (int i=uncut->indexLast[nreg-1];i<tod->ndata;i++)
	  tod->data[det][i]=tod->data[det][uncut->indexLast[nreg-1]-1];
      for (int reg=0;reg<nreg-1;reg++) {
	int i1=uncut->indexLast[reg];
	int i2=uncut->indexFirst[reg+1];	
	for (int i=i1;i<i2;i++)
	  tod->data[det][i]=tod->data[det][i1-1]+(tod->data[det][i2]-tod->data[det][i1-1])*(i-i1+1.0)/(i2-i1+1.0);
      }
    }
  }
  
}

/*--------------------------------------------------------------------------------*/

void clear_cut_data(mbTOD *tod)
{
  
  if (!(tod->have_data)) {
    fprintf(stderr,"Warning - attempting to clear non-existent data.\n");
    return;
  }
  if (tod->uncuts==NULL) {
    fprintf(stderr,"Warning - attempting to clear without uncuts.\n");
    return;
  }

#pragma omp parallel for shared(tod) default(none)
  for (int det=0;det<tod->ndet;det++) {
    mbUncut *uncut=tod->uncuts[tod->rows[det]][tod->cols[det]];
    int nreg=uncut->nregions;
    if (nreg>0) {
      if (uncut->indexFirst[0]>0)
	for (int i=0;i<uncut->indexFirst[0];i++)
	  tod->data[det][i]=0;
      if (uncut->indexLast[nreg-1]<tod->ndata)
	for (int i=uncut->indexLast[nreg-1];i<tod->ndata;i++)
	  tod->data[det][i]=0;
      for (int reg=0;reg<nreg-1;reg++) {
	int i1=uncut->indexLast[reg];
	int i2=uncut->indexFirst[reg+1];	
	for (int i=i1;i<i2;i++)
	  tod->data[det][i]=0;
      }
    }
  }
  
}
/*--------------------------------------------------------------------------------*/
bool is_det_used(mbTOD *tod, int det)
{
  if (tod) {
    assert(tod->rows);  //if we have a tod, make sure that it has row & column info
    assert(tod->cols);
    if (!tod->cuts) {
      if (det<tod->ndet)
	return true;
      else
	return false;
    }
    else
      return !mbCutsIsAlwaysCut(tod->cuts,tod->rows[det],tod->cols[det]);
  }
  return false;  //non-existent tod
}

/*--------------------------------------------------------------------------------*/
int get_starting_altaz_from_file(char *fname, actData **az_out, actData **alt_out, double **ctime_out)
{

  actData *az=NULL, *alt=NULL;
  double *ctime=NULL;
  
  int myid, nproc;      
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#else
  myid=0;
  nproc=1;
#endif
  int myargc;
  char *line_in=read_all_stdin(fname);
  char **myargv=create_argv_new(line_in,&myargc," \n");
  int nfiles=myargc/ALTAZ_PER_LINE;      
  assert(myargc==nfiles*ALTAZ_PER_LINE);  //if this fails, there's cruft in the parameter file
  free(line_in);
  
  int my_nfiles=0;
  for (int i=myid;i<nfiles;i+=nproc)
    my_nfiles++;
  
  if (my_nfiles==0) {
    free_argv(myargc,myargv);
    return 0;
  }
  az=vector(my_nfiles);
  alt=vector(my_nfiles);
  ctime=dvector(my_nfiles);
  for (int i=0;i<my_nfiles;i++) {
    int ind=myid+i*nproc;
    alt[i]=atof(myargv[ind*ALTAZ_PER_LINE])*M_PI/180.0;
    az[i]=atof(myargv[ind*ALTAZ_PER_LINE+1])*M_PI/180.0;
    ctime[i]=atof(myargv[ind*ALTAZ_PER_LINE+2]);
  }
  free_argv(myargc,myargv);
  
  *az_out=az;
  *alt_out=alt;
  *ctime_out=ctime;
  
  return my_nfiles;
  
}


/*--------------------------------------------------------------------------------*/
int read_all_tod_headers(TODvec *tods,PARAMS *params)
{

  mprintf(stdout,"Reading TOD headers now.\n");
  actData *az=NULL,*alt=NULL;
  double *ctime=NULL;
  int my_naltaz=0;
  if (strlen(params->altaz_file)) {
    my_naltaz=get_starting_altaz_from_file(params->altaz_file,&az,&alt,&ctime);
  }
  
  for (int i=0;i<tods->ntod;i++) {
    mbTOD *mytod=&(tods->tods[i]);
    char *myfroot=tods->my_fnames[i];

    long seed=mytod->seed;
    //printf("reading tod header.\n");
    mbTOD *tmp=read_dirfile_tod_header(myfroot);
    //printf("read.\n");

    memcpy(mytod,tmp,sizeof(mbTOD));
    mytod->seed=seed;
    mytod->cuts=mbCutsAlloc(mytod->nrow,mytod->ncol);
    
    mprintf(stdout,"file %s had %d detectors and %d data elements, rows and cols are %d %d.  dt=%12.4e\n",myfroot,mytod->ndet,mytod->ndata,mytod->nrow, mytod->ncol,mytod->deltat);
    mytod->pointingOffset=nkReadPointingOffset(params->pointing_file);  
    //printf("read pointing offsets.\n");
    cut_mispointed_detectors(mytod);
    //printf("cut mispointed detectors.\n");
    if (i<my_naltaz) {
      set_tod_starting_altaz_ctime(mytod,alt[i],az[i],ctime[i]);
      mprintf(stdout,"starting altaz/ctime are %10.5f %10.5f %12.2f on file %d\n",mytod->alt[0],mytod->az[0],mytod->ctime,i);
    }
    
    assign_tod_ra_dec(mytod);
    //currently not used - lives inside of read_tod_header_c.cpp
    find_pointing_pivots(mytod,0.5);
    printf("got pivots inside ninkasi .\n");
    //printf("ra/dec are assigned.\n");
    find_tod_radec_lims(mytod);    
    int myid=0;
#ifdef HAVE_MPI 
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
    //mprintf(stdout,"limits are %12.5f %12.5f %12.5f %12.5f\n",mytod->ramin,mytod->ramax,mytod->decmin,mytod->decmax);
    mprintf(stdout,"limits are %4d %3d %12.5f %12.5f %12.5f %12.5f %14.2f\n",i,myid,mytod->ramin,mytod->ramax,mytod->decmin,mytod->decmax,mytod->ctime);
    
  }
  if (my_naltaz) {
    free(alt);
    free(az);
    free(ctime);
  }
  return 0;
}
/*--------------------------------------------------------------------------------*/
void set_global_radec_lims(TODvec *tods)
{
  assert(tods->ntod>0);
  tods->ramin=tods->tods[0].ramin;
  tods->ramax=tods->tods[0].ramax;
  tods->decmin=tods->tods[0].decmin;
  tods->decmax=tods->tods[0].decmax;

  for (int i=1;i<tods->ntod;i++) {
    if (tods->tods[i].ramin<tods->ramin)
      tods->ramin=tods->tods[i].ramin;
    if (tods->tods[i].ramax>tods->ramax)
      tods->ramax=tods->tods[i].ramax;
    if (tods->tods[i].decmin<tods->decmin)
      tods->decmin=tods->tods[i].decmin;
    if (tods->tods[i].decmax>tods->decmax)
      tods->decmax=tods->tods[i].decmax;
  }
  tods->ramin-=EPS;
  tods->ramax+=EPS;
  tods->decmin-=EPS;
  tods->decmax+=EPS;


#ifdef HAVE_MPI  
  actData junk;
  MPI_Allreduce(&tods->ramin,&junk,1,MPI_NType,MPI_MIN,MPI_COMM_WORLD);
  tods->ramin=junk;
  MPI_Allreduce(&tods->decmin,&junk,1,MPI_NType,MPI_MIN,MPI_COMM_WORLD);
  tods->decmin=junk;
  MPI_Allreduce(&tods->ramax,&junk,1,MPI_NType,MPI_MAX,MPI_COMM_WORLD);
  tods->ramax=junk;
  MPI_Allreduce(&tods->decmax,&junk,1,MPI_NType,MPI_MAX,MPI_COMM_WORLD);
  tods->decmax=junk;
#endif

  //do mpi stuff here
    
}
/*--------------------------------------------------------------------------------*/
int setup_maps(MAPvec *maps, PARAMS *params)
{
  assert(maps->nmap>0);
  if (params) {
    if (params->use_input_limits) {
      MAP simmap;
      readwrite_simple_map(&simmap,params->inname,DOREAD);
      MAP *mymap=maps->maps[0];
      mymap->ramin=simmap.ramin;
      mymap->ramax=simmap.ramax;
      mymap->decmin=simmap.decmin;
      mymap->decmax=simmap.decmax;
      mymap->pixsize=simmap.pixsize;
      destroy_map(&simmap);
    }
  }
  for (int i=0;i<maps->nmap;i++) {
    MAP *mymap=maps->maps[i];
    assert(mymap->pixsize>0);
    assert(mymap->ramax>mymap->ramin);
    assert(mymap->decmax>mymap->decmin);
    
    mymap->nx=(mymap->ramax-mymap->ramin)/mymap->pixsize+1;
    mymap->ny=(mymap->decmax-mymap->decmin)/mymap->pixsize+1;
    mymap->npix=mymap->nx*mymap->ny;
    mymap->map=vector(mymap->npix);
  }
    
  return 0;
}

/*--------------------------------------------------------------------------------*/
MAP *make_blank_map_copy(const MAP *map)
{
  MAP *map_copy;
  map_copy=(MAP *)malloc_retry(sizeof(MAP));
  
  map_copy->pixsize=map->pixsize;
  map_copy->ramin=map->ramin;
  map_copy->ramax=map->ramax;
  map_copy->decmin=map->decmin;
  map_copy->decmax=map->decmax;
  map_copy->nx=map->nx;
  map_copy->ny=map->ny;
  map_copy->npix=map->npix;
  map_copy->have_locks=0;  //don't recycle locks.  will create them as needed, if needed.
#ifdef ACTPOL
  memcpy(map_copy->pol_state,map->pol_state,MAX_NPOL*sizeof(map->pol_state[0]));
#endif
  map_copy->map=(actData *)malloc_retry(sizeof(actData)*map_copy->npix*get_npol_in_map(map));
  //memcpy(map_copy->map,map->map,sizeof(actData)*map_copy->npix);
  clear_map(map_copy);
  return map_copy;
}

/*--------------------------------------------------------------------------------*/
MAP *deres_map(MAP *map)
{
  MAP *map_copy;
  map_copy=(MAP *)malloc_retry(sizeof(MAP));
  
  assert(map->projection->proj_type!=NK_HEALPIX_RING);  //if healpix, have to handle things a bit differently
  assert(map->projection->proj_type!=NK_HEALPIX_NEST);
  map_copy->pixsize=map->pixsize*2;
  map_copy->ramin=map->ramin;
  map_copy->ramax=map->ramax;
  map_copy->decmin=map->decmin;
  map_copy->decmax=map->decmax;
  map_copy->nx=map->nx/2;
  map_copy->ny=map->ny/2;
  map_copy->npix=map_copy->nx*map_copy->ny;
  //map_copy->projection=(nkProjection *)malloc_retry(sizeof(nkProjection));
  //memcpy(map_copy->projection,map->projection,sizeof(nkProjection));
  map_copy->projection=deres_projection(map->projection);
  map_copy->have_locks=0;  //don't recycle locks.  will create them as needed, if needed.
  map_copy->map=(actData *)malloc_retry(sizeof(actData)*map_copy->npix);
  //memcpy(map_copy->map,map->map,sizeof(actData)*map_copy->npix);
  
  memset(map_copy->map,0,sizeof(actData)*map_copy->npix);
  for (int i=0;i<map_copy->ny*2;i++) 
    for (int j=0;j<map_copy->nx*2;j++)
      map_copy->map[(i/2)*map_copy->nx+(j/2)]+=map->map[i*map->nx+j];
  
  return map_copy;
  
}
/*--------------------------------------------------------------------------------*/
MAP *upres_map(MAP *map)
{
  MAP *map_copy;
  map_copy=(MAP *)malloc_retry(sizeof(MAP));
  
  assert(map->projection->proj_type!=NK_HEALPIX_RING);  //if healpix, have to handle things a bit differently
  assert(map->projection->proj_type!=NK_HEALPIX_NEST);
  map_copy->pixsize=map->pixsize/2;
  map_copy->ramin=map->ramin;
  map_copy->ramax=map->ramax;
  map_copy->decmin=map->decmin;
  map_copy->decmax=map->decmax;
  map_copy->nx=map->nx*2;
  map_copy->ny=map->ny*2;
  map_copy->npix=map_copy->nx*map_copy->ny;
  //map_copy->projection=(nkProjection *)malloc_retry(sizeof(nkProjection));
  //memcpy(map_copy->projection,map->projection,sizeof(nkProjection));
  map_copy->projection=upres_projection(map->projection);
  map_copy->have_locks=0;  //don't recycle locks.  will create them as needed, if needed.
  map_copy->map=(actData *)malloc_retry(sizeof(actData)*map_copy->npix);
  //memcpy(map_copy->map,map->map,sizeof(actData)*map_copy->npix);
  for (int i=0;i<map_copy->ny;i++) 
    for (int j=0;j<map_copy->nx;j++)
      map_copy->map[i*map_copy->nx+j]=map->map[(i/2)*map->nx+(j/2)];
  
  return map_copy;
  
}
/*--------------------------------------------------------------------------------*/

MAP *make_map_copy(MAP *map)
{
  MAP *map_copy;
  //map_copy=(MAP *)malloc_retry(sizeof(MAP));
  map_copy=(MAP *)calloc(sizeof(MAP),1);

  map_copy->pixsize=map->pixsize;
  map_copy->ramin=map->ramin;
  map_copy->ramax=map->ramax;
  map_copy->decmin=map->decmin;
  map_copy->decmax=map->decmax;
  map_copy->nx=map->nx;
  map_copy->ny=map->ny;
  map_copy->npix=map->npix;
  map_copy->projection=(nkProjection *)malloc_retry(sizeof(nkProjection));
  memcpy(map_copy->projection,map->projection,sizeof(nkProjection));
  map_copy->have_locks=0;  //don't recycle locks.  will create them as needed, if needed.
  map_copy->map=(actData *)malloc_retry(sizeof(actData)*map_copy->npix*get_npol_in_map(map));
#ifdef ACTPOL
  memcpy(map_copy->pol_state,map->pol_state,MAX_NPOL*sizeof(map->pol_state[0]));
#endif
  memcpy(map_copy->map,map->map,sizeof(actData)*map_copy->npix*get_npol_in_map(map));
  return map_copy;
}
/*--------------------------------------------------------------------------------*/
void destroy_map(MAP *map)
{
  free(map->map);
  if (map->have_locks)
    free(map->locks);
  //free(map->projection);
  free(map);
}
/*--------------------------------------------------------------------------------*/
void destroy_mapset(MAPvec *maps)
{
  for (int i=0;i<maps->nmap;i++)
    destroy_map(maps->maps[i]);
  free(maps->maps);
  free(maps);
}
/*--------------------------------------------------------------------------------*/
MAPvec *make_mapset_copy(MAPvec *maps)
{
  MAPvec *maps_copy;
  assert(maps->nmap>0);
  maps_copy=(MAPvec *)malloc_retry(sizeof(MAPvec));
  maps_copy->maps=(MAP **)malloc_retry(sizeof(MAP *));
  maps_copy->nmap=maps->nmap;
  for (int i=0;i<maps->nmap;i++)
    maps_copy->maps[i]=make_map_copy((maps->maps[i]));
  return maps_copy;
  
}
/*--------------------------------------------------------------------------------*/
void clear_map(MAP *map)
{
  memset(map->map,0,sizeof(actData)*map->npix*get_npol_in_map(map));
}
/*--------------------------------------------------------------------------------*/
void clear_mapset(MAPvec *maps)
{
  for (int i=0;i<maps->nmap;i++)
    clear_map(maps->maps[i]);
}
/*--------------------------------------------------------------------------------*/
bool is_map_blank(MAP *map)
//return true if all elements of a map are 0.
{
  for (int i=0;i<map->npix;i++)
    if (map->map[i])
      return false;
  return true;
}
/*--------------------------------------------------------------------------------*/
bool is_mapset_blank(MAPvec *maps)
//return true if all maps in a mapset are blank (0).
{
  assert(maps);  
  for (int i=0;i<maps->nmap;i++)
    if (!is_map_blank(maps->maps[i]))
      return false;
  return true;
}
/*--------------------------------------------------------------------------------*/


#ifdef HAVE_MPI
int mpi_reduce_map(MAP *map)
{

  int ierr;
  actData *vec=vector(map->npix);
  memset(vec,0,sizeof(actData)*map->npix);
  ierr=MPI_Allreduce(map->map,vec,map->npix,MPI_NType,MPI_SUM,MPI_COMM_WORLD);
  memcpy(map->map,vec,sizeof(actData)*map->npix);
  free(vec);
  return ierr;
  
}

/*--------------------------------------------------------------------------------*/
int  mpi_reduce_mapset(MAPvec *maps)
{
  int ierr;
  for (int i=0;i<maps->nmap;i++) {
    ierr=mpi_reduce_map(maps->maps[i]);
    assert(ierr==0);
  }
  return ierr;
}
#endif
/*--------------------------------------------------------------------------------*/
#ifdef HAVE_MPI
int mpi_broadcast_map(MAP *map, int master)
{
  int ierr=MPI_Bcast(map->map,map->npix,MPI_NType,master,MPI_COMM_WORLD);
  return ierr;
}

/*--------------------------------------------------------------------------------*/
int mpi_broadcast_mapset(MAPvec *maps,int master)
{
  int ierr;
  for (int i=0;i<maps->nmap;i++) {
    ierr=mpi_broadcast_map(maps->maps[i],master);
    assert(ierr==0);
  }
  return ierr;
  
}
#endif
/*--------------------------------------------------------------------------------*/

int read_tod_data(mbTOD *tod)
{

  if (tod->have_data==0)
    tod->data=matrix(tod->ndet,tod->ndata);
  tod->have_data=1;  
  //printf("clearing tod.\n");
  clear_tod(tod);
  //printf("reading tod.\n");
  read_dirfile_tod_data (tod);
  return 0;
}
/*--------------------------------------------------------------------------------*/

void reverse_tod_data(mbTOD *tod)
//time-reverse all the data in a tod
{
  assert(tod);
  assert(tod->have_data);
  int nn=tod->ndata;
#pragma omp parallel for shared(tod,nn) default(none)
  for (int i=0;i<tod->ndet;i++) {
    for (int j=0;j<nn/2;j++) {
      actData tmp=tod->data[i][j];
      tod->data[i][j]=tod->data[i][nn-1-j];
      tod->data[i][nn-1-j]=tmp;
    }
    
  }
}
/*--------------------------------------------------------------------------------*/
void createFFTWplans1TOD(mbTOD *mytod)
{  
  if (mytod->have_plans==false) {
    actData *vec=vector(mytod->ndata);
    act_fftw_complex *vec2=(act_fftw_complex *)act_fftw_malloc(sizeof(act_fftw_complex)*mytod->ndata);
    mytod->p_forward=act_fftw_plan_dft_r2c_1d(mytod->ndata,vec,vec2,FFTW_ESTIMATE);
    mytod->p_back=act_fftw_plan_dft_c2r_1d(mytod->ndata,vec2,vec,FFTW_ESTIMATE);
    free(vec);
    mytod->have_plans=true;
  }
  
}
/*--------------------------------------------------------------------------------*/

void createFFTWplans(TODvec *tod)
{
  for (int i=0;i<tod->ntod;i++)
    {
      createFFTWplans1TOD(&(tod->tods[i]));
    }
}

/*--------------------------------------------------------------------------------*/
void filter_data(mbTOD *tod)
{
  filter_data_wnoise(tod);
}



/*--------------------------------------------------------------------------------*/
void add_noise_to_tod(mbTOD *tod, PARAMS *params)
{

#if 1
  add_noise_to_tod_new(tod);  
#else
  assert(tod->have_avec);
  assert(tod->have_data);
#pragma omp parallel shared(tod) default(none)
  {
    act_fftw_complex *vec=(act_fftw_complex *)act_fftw_malloc(sizeof(act_fftw_complex)*tod->ndata);
#pragma omp for schedule(dynamic,1)    
    for (int i=0;i<tod->ndet;i++) {
      act_fftw_execute_dft_r2c(tod->p_forward,tod->data[i],vec);
      unsigned idum=tod->seed+i;  //tod->seed should already be set such that this guarantees uniqueness.
      for (int j=0;j<tod->ndata/2+1;j++) { 
	if (j>0) {//don't add an overall offset
	  vec[j][0] +=sqrt(tod->avec[j])*mygasdev(&idum);
	  vec[j][1] +=sqrt(tod->avec[j])*mygasdev(&idum);
	}
	else {
	  vec[j][0]=0;
	  vec[j][1]=0;
	}

	//vec[i][0]/=(tod->avec[i]*tod->n);
	//vec[i][1]/=(tod->avec[i]*tod->n);
	vec[j][0]/=(actData)(tod->ndata);
	vec[j][1]/=(actData)(tod->ndata);
	

      }
      act_fftw_execute_dft_c2r(tod->p_back,vec,tod->data[i]);
      
    }
    act_fftw_free(vec); 
  }
  //detrend_data(tod);
#endif
}


/*--------------------------------------------------------------------------------*/
void get_pointing_vec(const mbTOD *tod, const MAP *map, int det, int *ind)
{
 #if 1
  PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
  get_pointing_vec_new(tod,map,det,ind,scratch);
  destroy_pointing_fit_scratch(scratch);
#else

 //printf("limits are %14.5e %14.5e %14.5e %14.5e %d %d\n",map->ramin,map->ramax,map->decmin,map->decmax,map->nx,map->ny);
  for (int i=0;i<tod->ndata;i++) {
    actData x=tod->ra[i]+tod->dra[det];
    actData y=tod->dec[i]+tod->ddec[det];
    ind[i]=(int)((y-map->decmin)/map->pixsize)+map->ny*(int)((x-map->ramin)/map->pixsize);
  }
#endif
}

/*--------------------------------------------------------------------------------*/
void get_pointing_vec_new(const mbTOD *tod, const MAP *map, int det, int *ind, PointingFitScratch *scratch)
//az,alt,ra,dec are scratch spaces of length(tod->ndata)
//for performance reasons, please check things before here.
{

#if 1
  get_map_projection(tod,map,det,ind,scratch);
#else
  
  get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
  //get_radec_from_altaz_fit_1det(tod,det,scratch);

  for (int i=0;i<tod->ndata;i++) {
    ind[i]=(int)((scratch->dec[i]-map->decmin)/map->pixsize)+map->ny*(int)((scratch->ra[i]-map->ramin)/map->pixsize);    
  }
#endif 
  
}

/*--------------------------------------------------------------------------------*/

void write_tod_pointing_to_disk(const mbTOD *tod, char *fname)
/*Dump ra/dec of each sample to disk.  For use in external analysis.*/
{
  FILE *outfile=fopen(fname,"w");
  assert(outfile);
  
  int ngood=0;
  for (int det=0;det<tod->ndet;det++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[det],tod->cols[det]))
      ngood++;
  }
  PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
  fwrite(&ngood,sizeof(int),1,outfile);
  fwrite(&tod->ndata,sizeof(int),1,outfile);
  for (int det=0;det<tod->ndet;det++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[det],tod->cols[det])) {
      get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
      fwrite(scratch->ra,sizeof(actData),tod->ndata,outfile);
      fwrite(scratch->dec,sizeof(actData),tod->ndata,outfile);
    }    
  }
  fclose(outfile);
  
  
}

/*--------------------------------------------------------------------------------*/
void setup_omp_locks(MAP *map) 
{
  if (map->have_locks) {
    mprintf(stdout,"skipping lock creation as they already exist.\n");
    return;
  }
  if (map->nx>map->ny) {
    map->nlock=map->ny;
    map->lock_len=map->nx;
  }
  else {
    map->nlock=map->nx;
    map->lock_len=map->ny;    
  }
  map->locks=(omp_lock_t *)malloc_retry(sizeof(omp_lock_t)*map->nlock);
  for (int i=0;i<map->nlock;i++) {
    omp_init_lock(&(map->locks[i]));    
  }
  map->have_locks=1;
  
}
/*--------------------------------------------------------------------------------*/
void omp_reduce_map(MAP *map,MAP *mymap)
{
  assert(map->have_locks);  //do not setup locks here - you might do this in parallel accidentally.  Could be ugly.
  //  if (!map->have_locks) {
  //  setup_omp_locks(map);    
  //}
  for (int i=0;i<map->nlock;i++) {
    omp_set_lock(&map->locks[i]);
    for (int j=0;j<map->lock_len;j++) {
      map->map[i*map->lock_len+j]+= mymap->map[i*map->lock_len+j];
    }
    omp_unset_lock(&map->locks[i]);
  }

}

/*--------------------------------------------------------------------------------*/
bool is_det_listed(const mbTOD *tod, const PARAMS *params, int det) 
//checks parameter file to see if detector is in the included list.
{

  if (params==NULL)
    return true;

  bool use_row;
  bool use_col;
  if (params->n_use_rows==0) {
    use_row=true;
  }
  else {
    use_row=false;
    for (int i=0;i<params->n_use_rows;i++) {
      if (tod->rows[det]==params->use_rows[i])
	use_row=true;
    }
  }
  
  if (params->n_use_cols==0) {
    use_col=true;
  }
  else {
    use_col=false;
    for (int i=0;i<params->n_use_cols;i++) {
      if (tod->cols[det]==params->use_cols[i])
	use_col=true;
    }
  }
  return (use_row&&use_col);


}
/*--------------------------------------------------------------------------------*/
actData how_far_am_i_from_radec(actData rah,actData ram,actData ras,actData dd,actData dm,actData ds,mbTOD *tod)
{
  actData ra=(rah+ram/60*ras/3600)*15*M_PI/180.0;
  actData dec=(dd+dm/60+ds/3600)*M_PI/180.0;

  return how_far_am_i_from_radec_radians(ra,dec,tod)*180.0/M_PI;
}
/*--------------------------------------------------------------------------------*/
actData how_far_am_i_from_radec_radians(actData ra, actData dec, mbTOD *tod)
{
  actData bigmin=1000;
#pragma omp parallel shared(tod,ra,dec,bigmin)  default(none)
  {
    actData mindist=bigmin;
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    actData mycos=cos(dec);
#pragma omp for schedule(dynamic,1)
    for (int i=0;i<tod->ndet;i++) 
      {
	if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])))   { //&&(is_det_listed(tod,params,i))) {
	  get_radec_from_altaz_fit_1det_coarse(tod,i,scratch);
	  for (int j=0;j<tod->ndata;j++) {
	    actData dx=(ra-scratch->ra[j])*mycos;
	    actData dy=dec-scratch->dec[j];
	    actData dd=dx*dx+dy*dy;
	    if (dd<mindist)
	      mindist=dd;
	  }
	}
      }
    destroy_pointing_fit_scratch(scratch);
#pragma omp critical
    if (mindist<bigmin)
      bigmin=mindist;
  }
  
  return sqrt(bigmin);
}
/*--------------------------------------------------------------------------------*/
actData  *how_far_am_i_from_radec_radians_vec( actData *ra, actData *dec, int npt, mbTOD *tod)
{
  actData *dists=vector(npt);
  actData bigmin=1000;
  for (int i=0;i<npt;i++)
    dists[i]=bigmin;
#pragma omp parallel shared(tod,ra,dec,bigmin,npt,dists)  default(none)
  {
    actData *mindist=vector(npt);
    for (int i=0;i<npt;i++)
      mindist[i]=bigmin;
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    //actData mycos=cos(dec);
    actData *mycos=vector(npt);
    for (int k=0;k<npt;k++)
      mycos[k]=cos(dec[k]);
#pragma omp for schedule(dynamic,1)
    for (int i=0;i<tod->ndet;i++) 
      {
	if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])))   { //&&(is_det_listed(tod,params,i))) {
	  get_radec_from_altaz_fit_1det_coarse(tod,i,scratch);
	  for (int j=0;j<tod->ndata;j++) {
	    for (int k=0;k<npt;k++) {
	      actData dx=(ra[k]-scratch->ra[j])*mycos[k];
	      actData dy=dec[k]-scratch->dec[j];
	      actData dd=dx*dx+dy*dy;
	      if (dd<mindist[k])
		mindist[k]=dd;
	    }
	  }
	}
      }
    destroy_pointing_fit_scratch(scratch);
#pragma omp critical
    for (int k=0;k<npt;k++)
      if (mindist[k]<dists[k])
	dists[k]=mindist[k];
  }
  for (int k=0;k<npt;k++)
    dists[k]=sqrt(dists[k]);
  
  return dists;
}


/*--------------------------------------------------------------------------------*/
void generate_index_mapping(MAP *map, mbTOD *tod, int **inds)
{
  assert(map);
  assert(map->projection);
#pragma omp parallel shared(map,tod,inds)
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
#pragma omp for schedule(dynamic,4)
    for (int i=0;i<tod->ndet;i++) 
      if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,NULL,i)))
	get_pointing_vec_new(tod,map,i,inds[i],scratch);  
    destroy_pointing_fit_scratch(scratch);   
   }
}
/*--------------------------------------------------------------------------------*/
void find_map_index_limits(MAP *map, mbTOD *tod, int *imin_out, int *imax_out)
{
  int imin=map->npix+2;
  int imax=-1;
#pragma omp parallel shared(map,tod,imin,imax) default(none)
  {
    int myimin=imin;
    int myimax=imax;

    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    int *itmp=(int *)malloc(sizeof(int)*tod->ndata);
#pragma omp for schedule(dynamic,4)
    for (int i=0;i<tod->ndet;i++) 
      if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,NULL,i))) {
        get_pointing_vec_new(tod,map,i,itmp,scratch);
	if (tod->kept_data) {
	  int row=tod->rows[i];
	  int col=tod->cols[i];
	  mbUncut *uncut=tod->kept_data[row][col];
	  for (int region=0;region<uncut->nregions;region++)
	    for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++) {
	      if (itmp[j]>myimax)
		myimax=itmp[j];
	      if (itmp[j]<myimin)
		myimin=itmp[j];
	    }
	  
	}
	else {
	  if (tod->uncuts) {
	    int row=tod->rows[i];
	    int col=tod->cols[i];
	    mbUncut *uncut=tod->uncuts[row][col];
	    for (int region=0;region<uncut->nregions;region++)
	      for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++) {
		if (itmp[j]>myimax)
		  myimax=itmp[j];
		if (itmp[j]<myimin)
		  myimin=itmp[j];
	      }
	  }
	}

      }
#pragma omp critical
    {
      if (myimin<imin)
	imin=myimin;
      if (myimax>imax)
	imax=myimax;
    }
    destroy_pointing_fit_scratch(scratch);
    free(itmp);
  }
  *imin_out=imin;
  *imax_out=imax;
  
}
/*--------------------------------------------------------------------------------*/
int *tod2map_actpol(MAP *map, mbTOD *tod, int *ipiv_proc)
//project a tod into a map, hopefully with better ompness.  It wants to know the
//index limits for processes when mapping TODs into a map.  If ipiv_proc comes in as null,
//it will calculate some for you.
{

  //float **sin2gamma=fmatrix(tod->ndet,tod->ndata);  //eventually polarization goes in.
  assert(map);
  assert(map->projection);
  if (ipiv_proc==NULL) {
    printf("Nope, pivot finding not yet implemented.\n");
    return NULL;    
  }
    
  int **inds=imatrix(tod->ndet,tod->ndata);
  generate_index_mapping(map,tod, inds);
  

  
#pragma omp parallel shared(map,tod,inds,ipiv_proc) default(none)
  {
    int myid=omp_get_thread_num();
    
#pragma omp for
    for (int i=0;i<tod->ndet;i++) {
      if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,NULL,i))) {      
	if (tod->kept_data) {
	  int row=tod->rows[i];
	  int col=tod->cols[i];
	  mbUncut *uncut=tod->kept_data[row][col];
	  for (int region=0;region<uncut->nregions;region++)
	    for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++) {
	      int itmp=inds[i][j];
	      if ((itmp>=ipiv_proc[myid])&&(itmp<ipiv_proc[myid+1]))
		map->map[itmp]+=tod->data[i][j];
	    }
	}
	else {
	  if (tod->uncuts) {
	    int row=tod->rows[i];
	    int col=tod->cols[i];
	    mbUncut *uncut=tod->uncuts[row][col];
	    for (int region=0;region<uncut->nregions;region++)
	      for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++) {
		int itmp=inds[i][j];
		if ((itmp>=ipiv_proc[myid])&&(itmp<ipiv_proc[myid+1]))
		  map->map[itmp]+=tod->data[i][j];
	      }
	  }
	  else
	    for (int j=0;j<tod->ndata;j++)
	      map->map[inds[i][j]]+=tod->data[i][j];
	}
      }
    }
  }
  free(inds[0]);
  free(inds);

  return NULL;

}

/*--------------------------------------------------------------------------------*/
void tod2map_nocopy(MAP *map, mbTOD *tod,PARAMS *params)
{
  int **inds=imatrix(tod->ndet,tod->ndata);

  assert(map);
  assert(map->projection);
#pragma omp parallel shared(map,tod,inds)
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
#pragma omp for schedule(dynamic,4)
    for (int i=0;i<tod->ndet;i++) { 
      if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,params,i)))
	get_pointing_vec_new(tod,map,i,inds[i],scratch);  
    }
    destroy_pointing_fit_scratch(scratch);
    
   }

  for (int i=0;i<tod->ndet;i++) {
    if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,params,i))) {
      
      if (tod->kept_data) {
	int row=tod->rows[i];
	int col=tod->cols[i];
	mbUncut *uncut=tod->kept_data[row][col];
	for (int region=0;region<uncut->nregions;region++)
	  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
	    map->map[inds[i][j]]+=tod->data[i][j];
      }
      else {
	if (tod->uncuts) {
	  int row=tod->rows[i];
	  int col=tod->cols[i];
	  mbUncut *uncut=tod->uncuts[row][col];
	  for (int region=0;region<uncut->nregions;region++)
	    for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
	      map->map[inds[i][j]]+=tod->data[i][j];
	}
	else
	  for (int j=0;j<tod->ndata;j++)
	    map->map[inds[i][j]]+=tod->data[i][j];
      }
    }
  }
  free(inds[0]);
  free(inds);


      
}
/*--------------------------------------------------------------------------------*/
int is_map_polarized(MAP *map) 
{
#ifdef ACTPOL
  //return (map->pol_state>1);
  //If anything past the first polarization slot is populated, we're polarized
  for (int i=1;i<MAX_NPOL;i++)
    if (map->pol_state[i])
      return 1;
  return 0;
#else
  return 0;
#endif
}
/*--------------------------------------------------------------------------------*/
void tod2map(MAP *map, mbTOD *tod, PARAMS *params)
{


  assert(map);
  assert(map->projection);

#ifdef ACTPOLFWEEE
  if (is_map_polarized(map)) {       
    return;
  }
#endif


  int nproc;
  //printf("projecting TOD into map.\n");
#pragma omp parallel shared(nproc) default(none)
#pragma omp single
  nproc=omp_get_num_threads();
  
  if (nproc*map->npix*sizeof(actData)>tod->ndata*tod->ndet*sizeof(int)) {
    //printf("doing index-saving projection.\n");
    tod2map_nocopy(map,tod,params);
    return;
  }
  //else
  //printf("doing old projection.\n");



  if (!map->have_locks) {
    //printf("setting up locks.\n");
    setup_omp_locks(map);
    //printf("set them up.\n");
    
  }
#pragma omp parallel shared(tod,map,params) default(none)
  { 
    MAP *mymap=make_blank_map_copy(map);    
    int *ind=(int *)malloc_retry(sizeof(int)*tod->ndata);

    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
#pragma omp for nowait
    for (int i=0;i<tod->ndet;i++) { 
      if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,params,i))) {
	//printf("doing first stuff.\n");
	get_pointing_vec_new(tod,map,i,ind,scratch);
	//printf("got pointing vec.\n");
	if (tod->kept_data) {
	  int row=tod->rows[i];
	  int col=tod->cols[i];
	  mbUncut *uncut=tod->kept_data[row][col];
	  for (int region=0;region<uncut->nregions;region++) 
	    for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
	      mymap->map[ind[j]]+=tod->data[i][j];	  
	}
	else {
	  if (tod->uncuts) {
	    //printf("doing uncuts.\n");
	    int row=tod->rows[i];
	    int col=tod->cols[i];
	    mbUncut *uncut=tod->uncuts[row][col];
	    for (int region=0;region<uncut->nregions;region++) 
	      for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
		mymap->map[ind[j]]+=tod->data[i][j];
	  }
	  else 
	    for (int j=0;j<tod->ndata;j++) 
	      mymap->map[ind[j]]+=tod->data[i][j];
	}
      }
    }
    
    //printf("ready to accumulate.\n");
    
    omp_reduce_map(map,mymap);
    free(ind);
    destroy_pointing_fit_scratch(scratch);
    destroy_map(mymap);
  } 
}
/*--------------------------------------------------------------------------------*/
void polmap2tod_old(MAP *map, mbTOD *tod)
{
  assert(tod);
  assert(tod->data);
  assert(tod->pixelization_saved);
  int cur_pol=0;
  //loop through possible polarization states and project the ones we find.
  for (int pol_ind=0;pol_ind<MAX_NPOL;pol_ind++)
    if (map->pol_state[pol_ind]) {
      if (pol_ind==0) {  //I 
#pragma omp parallel shared(pol_ind,map,tod,cur_pol)
	for (int det=0;det<tod->ndet;det++) {
	  int row=tod->rows[det];
	  int col=tod->cols[det];
	  mbUncut *uncut=tod->uncuts[row][col];
	  for (int region=0;region<uncut->nregions;region++)
	    for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
	      tod->data[det][j]+=map->map[tod->pixelization_saved[det][j]];
	}
	cur_pol++;
      }
      if (pol_ind==1) {  //Q
#pragma omp parallel shared(pol_ind,map,tod,cur_pol)
        for (int det=0;det<tod->ndet;det++) {
          int row=tod->rows[det];
          int col=tod->cols[det];
          mbUncut *uncut=tod->uncuts[row][col];
          for (int region=0;region<uncut->nregions;region++)
            for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
              tod->data[det][j]+=map->map[tod->pixelization_saved[det][j]+cur_pol*map->npix]*sin(tod->twogamma_saved[det][j]);
	}
	cur_pol++;
      }
      if (pol_ind==2) {  //U
#pragma omp parallel shared(pol_ind,map,tod,cur_pol)
        for (int det=0;det<tod->ndet;det++) {
          int row=tod->rows[det];
          int col=tod->cols[det];
          mbUncut *uncut=tod->uncuts[row][col];
          for (int region=0;region<uncut->nregions;region++)
            for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
              tod->data[det][j]+=map->map[tod->pixelization_saved[det][j]+cur_pol*map->npix]*cos(tod->twogamma_saved[det][j]);
	}
	cur_pol++;
      }
      
    }
}
/*--------------------------------------------------------------------------------*/

void polmap2tod(MAP *map, mbTOD *tod)
{
  assert(tod);
  assert(tod->data);
  assert(tod->pixelization_saved);
  assert(tod->uncuts);


  const int npol=get_npol_in_map(map);
  const int poltag=get_map_poltag(map);
  const int npix=map->npix;
  if (poltag==POL_ERROR) {
    fprintf(stderr,"Error - unrecognized combination in polmap2tod.\n");
    return;
  }
#pragma omp parallel shared(map,tod) default(none)
  {
    
    const ACTpolPointingFit *pfit=tod->actpol_pointing;
    const actData *az=tod->az;
    actData ninv=1.0/tod->ndata;


    const actData *mymap=map->map;
    switch(poltag){
    case POL_I:
#pragma omp for
      for (int det=0;det<tod->ndet;det++) {
	int row=tod->rows[det];
	int col=tod->cols[det];
	mbUncut *uncut=tod->uncuts[row][col];
	for (int region=0;region<uncut->nregions;region++) {
	  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
	    tod->data[det][j]+=mymap[tod->pixelization_saved[det][j]];
	}
      }
      break;
    case POL_IQU:
#pragma omp for
      for (int det=0;det<tod->ndet;det++) {
	actData ctime_sin=pfit->gamma_ctime_sin_coeffs[det]*ninv;
	actData ctime_cos=pfit->gamma_ctime_cos_coeffs[det]*ninv;
	const actData *az_sin=pfit->gamma_az_sin_coeffs[det];
	const actData *az_cos=pfit->gamma_az_cos_coeffs[det];

	int row=tod->rows[det];
	int col=tod->cols[det];
	mbUncut *uncut=tod->uncuts[row][col];
	for (int region=0;region<uncut->nregions;region++) {
	  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++){

	    actData mysin,mycos;
	    actData aa=az[j];
	    mysin=az_sin[3]+aa*(az_sin[2]+aa*(az_sin[1]+aa*(az_sin[0])))+ctime_sin*j;
	    mycos=az_cos[3]+aa*(az_cos[2]+aa*(az_cos[1]+aa*(az_cos[0])))+ctime_cos*j;

	    int jj=tod->pixelization_saved[det][j]*npol;
	    tod->data[det][j]+=mymap[jj];
	    tod->data[det][j]+=mymap[jj+1]*mycos;
	    tod->data[det][j]+=mymap[jj+2]*mysin;
	  }
	}
      }
      break;
    case POL_QU:
#pragma omp for
      for (int det=0;det<tod->ndet;det++) {
	int row=tod->rows[det];
	int col=tod->cols[det];
	mbUncut *uncut=tod->uncuts[row][col];
	for (int region=0;region<uncut->nregions;region++) {
	  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++){
	    actData mycos=cos7_pi(tod->twogamma_saved[det][j]);
	    actData mysin=sin7_pi(tod->twogamma_saved[det][j]);
	    int jj=tod->pixelization_saved[det][j]*npol;
	    tod->data[det][j]+=mymap[jj]*mycos;
	    tod->data[det][j]+=mymap[jj+1]*mysin;
	  }
	}
      }
      break;
    default:
      printf("Error - unsupported poltag in polmap2tod.\n");
      break;
    }
  }
}

/*--------------------------------------------------------------------------------*/
void tod2polmap(MAP *map,mbTOD *tod)
{
  assert(tod);
  assert(tod->data);
  assert(tod->pixelization_saved);
  assert(tod->uncuts);
  int cur_pol=0;
  //loop through possible polarization states and project the ones we find.
  for (int pol_ind=0;pol_ind<MAX_NPOL;pol_ind++)
    if (map->pol_state[pol_ind]) {
      if (pol_ind==0) {  //I 
#pragma omp parallel for shared(pol_ind,map,tod,cur_pol)
	for (int det=0;det<tod->ndet;det++) {
	  int row=tod->rows[det];
	  int col=tod->cols[det];
	  mbUncut *uncut=tod->uncuts[row][col];
	  for (int region=0;region<uncut->nregions;region++) {
	    for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
#pragma omp atomic
	      map->map[tod->pixelization_saved[det][j]]+=tod->data[det][j];
	  }
	}
	cur_pol++;
      }
      if (pol_ind==1) {  //Q
#pragma omp parallel shared(pol_ind,map,tod,cur_pol)
        for (int det=0;det<tod->ndet;det++) {
          int row=tod->rows[det];
          int col=tod->cols[det];
          mbUncut *uncut=tod->uncuts[row][col];
          for (int region=0;region<uncut->nregions;region++)
            for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
#pragma omp atomic
	      map->map[tod->pixelization_saved[det][j]+cur_pol*map->npix]+=tod->data[det][j]*sin(tod->twogamma_saved[det][j]);
	}
	cur_pol++;
      }
      if (pol_ind==2) {  //U
#pragma omp parallel shared(pol_ind,map,tod,cur_pol)
        for (int det=0;det<tod->ndet;det++) {
          int row=tod->rows[det];
          int col=tod->cols[det];
          mbUncut *uncut=tod->uncuts[row][col];
          for (int region=0;region<uncut->nregions;region++)
            for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
#pragma omp atomic
	      map->map[tod->pixelization_saved[det][j]+cur_pol*map->npix]+=tod->data[det][j]*cos(tod->twogamma_saved[det][j]);
	}
	cur_pol++;
      }
      if (pol_ind==3) {  //Q^2
#pragma omp parallel shared(pol_ind,map,tod,cur_pol)
        for (int det=0;det<tod->ndet;det++) {
          int row=tod->rows[det];
          int col=tod->cols[det];
          mbUncut *uncut=tod->uncuts[row][col];
          for (int region=0;region<uncut->nregions;region++)
            for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++) {
	      actData mycos=cos(tod->twogamma_saved[det][j]);
#pragma omp atomic
	      map->map[tod->pixelization_saved[det][j]+cur_pol*map->npix]+=tod->data[det][j]*mycos*mycos;
	    }
	}
	cur_pol++;
      }

      if (pol_ind==4) {  //Q*U
#pragma omp parallel shared(pol_ind,map,tod,cur_pol)
        for (int det=0;det<tod->ndet;det++) {
          int row=tod->rows[det];
          int col=tod->cols[det];
          mbUncut *uncut=tod->uncuts[row][col];
          for (int region=0;region<uncut->nregions;region++)
            for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++) {
	      actData mycos=cos(tod->twogamma_saved[det][j]);
	      actData mysin=sin(tod->twogamma_saved[det][j]);
#pragma omp atomic
	      map->map[tod->pixelization_saved[det][j]+cur_pol*map->npix]+=tod->data[det][j]*mycos*mysin;
	    }
	}
	cur_pol++;
      }
      
      if (pol_ind==5) {  //U^2
#pragma omp parallel shared(pol_ind,map,tod,cur_pol)
        for (int det=0;det<tod->ndet;det++) {
          int row=tod->rows[det];
          int col=tod->cols[det];
          mbUncut *uncut=tod->uncuts[row][col];
          for (int region=0;region<uncut->nregions;region++)
            for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++) {
	      actData mysin=sin(tod->twogamma_saved[det][j]);
#pragma omp atomic
	      map->map[tod->pixelization_saved[det][j]+cur_pol*map->npix]+=tod->data[det][j]*mysin*mysin;
	    }
	}
	cur_pol++;
      }
      
    }
}
/*--------------------------------------------------------------------------------*/
static inline void get_twogamma_sincos(actData *sinval, actData *cosval, const actData *sin_azparams, const actData *cos_azparams, int nazparams, actData az,actData sin_tvec_param, actData cos_tvec_param,actData tfrac)
{
  *sinval=0;
  *cosval=0;
  for (int i=0;i<4;i++) {
    *sinval=sin_azparams[i]+az*(*sinval);
    *cosval=cos_azparams[i]+az*(*cosval);
  }
  *sinval+=sin_tvec_param*tfrac;
  *cosval+=cos_tvec_param*tfrac;
}
/*--------------------------------------------------------------------------------*/

void tod2polmap_copy(MAP *map,mbTOD *tod) 
//copy refers to each thread having a copy of the map.  
{
  assert(tod);
  assert(tod->data);
  assert(tod->pixelization_saved);
  assert(tod->uncuts);

  const int npol=get_npol_in_map(map);
  const int poltag=get_map_poltag(map);
  const int npix=map->npix;

  if (poltag==POL_ERROR) {
    fprintf(stderr,"Error - unrecognized combination in tod2polmap_copy.\n");
    return;
  }
#pragma omp parallel shared(map,tod) default(none)
  {
    actData mysin,mycos;
    const ACTpolPointingFit *pfit=tod->actpol_pointing;
    const actData *az=tod->az;
    actData *mymap=vector(npol*map->npix);
    actData ninv=1.0/tod->ndata;
    memset(mymap,0,npol*map->npix*sizeof(actData));
    switch(poltag){
    case POL_I: 
#pragma omp for
      for (int det=0;det<tod->ndet;det++) {
	int row=tod->rows[det];
	int col=tod->cols[det];
	mbUncut *uncut=tod->uncuts[row][col];
	for (int region=0;region<uncut->nregions;region++) {
	  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
	    mymap[tod->pixelization_saved[det][j]]+=tod->data[det][j];	  
	}
      }
      break;
    case POL_IQU:
#pragma omp for
      for (int det=0;det<tod->ndet;det++) {

	actData ctime_sin=pfit->gamma_ctime_sin_coeffs[det]*ninv;
	actData ctime_cos=pfit->gamma_ctime_cos_coeffs[det]*ninv;
	const actData *az_sin=pfit->gamma_az_sin_coeffs[det];
	const actData *az_cos=pfit->gamma_az_cos_coeffs[det];

	int row=tod->rows[det];
	int col=tod->cols[det];
	mbUncut *uncut=tod->uncuts[row][col];
	for (int region=0;region<uncut->nregions;region++) {
	  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++){
	    
	    //get_twogamma_sincos(&mysin,&mycos,pfit->gamma_az_sin_coeffs[det],pfit->gamma_az_cos_coeffs[det],pfit->n_gamma_az_coeffs,tod->az[j],pfit->gamma_ctime_sin_coeffs[det],pfit->gamma_ctime_cos_coeffs[det],j*ninv);
	    actData aa=az[j];
	    mysin=az_sin[3]+aa*(az_sin[2]+aa*(az_sin[1]+aa*(az_sin[0])))+ctime_sin*j;
	    mycos=az_cos[3]+aa*(az_cos[2]+aa*(az_cos[1]+aa*(az_cos[0])))+ctime_cos*j;

	    int jj=tod->pixelization_saved[det][j]*npol;
	    mymap[jj]+=tod->data[det][j];
	    mymap[jj+1]+=tod->data[det][j]*mycos;
	    mymap[jj+2]+=tod->data[det][j]*mysin;

	  }
	}
      }
      break;

    case POL_QU:
#pragma omp for
      for (int det=0;det<tod->ndet;det++) {
	int row=tod->rows[det];
	int col=tod->cols[det];
	mbUncut *uncut=tod->uncuts[row][col];
	for (int region=0;region<uncut->nregions;region++) {
	  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++){
#if 1
	    actData mycos=cos7_pi(tod->twogamma_saved[det][j]);
	    actData mysin=sin7_pi(tod->twogamma_saved[det][j]);
	    int jj=tod->pixelization_saved[det][j]*npol;
	    mymap[jj]+=tod->data[det][j]*mycos;
	    mymap[jj+1]+=tod->data[det][j]*mysin;
	    
#else
	    mymap[tod->pixelization_saved[det][j]]+=tod->data[det][j]*cos(tod->twogamma_saved[det][j]);
	    mymap[tod->pixelization_saved[det][j]+npix]+=tod->data[det][j]*sin(tod->twogamma_saved[det][j]);
#endif
	  }
	}
      }
      break;

    case POL_IQU_PRECON:
#pragma omp for
      for (int det=0;det<tod->ndet;det++) {
	actData ctime_sin=pfit->gamma_ctime_sin_coeffs[det]*ninv;
	actData ctime_cos=pfit->gamma_ctime_cos_coeffs[det]*ninv;
	const actData *az_sin=pfit->gamma_az_sin_coeffs[det];
	const actData *az_cos=pfit->gamma_az_cos_coeffs[det];

	int row=tod->rows[det];
	int col=tod->cols[det];
	mbUncut *uncut=tod->uncuts[row][col];
	for (int region=0;region<uncut->nregions;region++) {
	  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++){

	    //get_twogamma_sincos(&mysin,&mycos,pfit->gamma_az_sin_coeffs[det],pfit->gamma_az_cos_coeffs[det],pfit->n_gamma_az_coeffs,tod->az[j],pfit->gamma_ctime_sin_coeffs[det],pfit->gamma_ctime_cos_coeffs[det],j*ninv);
	    actData aa=az[j];
	    mysin=az_sin[3]+aa*(az_sin[2]+aa*(az_sin[1]+aa*(az_sin[0])))+ctime_sin*j;
	    mycos=az_cos[3]+aa*(az_cos[2]+aa*(az_cos[1]+aa*(az_cos[0])))+ctime_cos*j;
	    
	    actData mycos=cos7_pi(tod->twogamma_saved[det][j]);
	    actData mysin=sin7_pi(tod->twogamma_saved[det][j]);

	    int jj=tod->pixelization_saved[det][j]*npol;
	    
	    mymap[jj]+=tod->data[det][j];
	    mymap[jj+1]+=tod->data[det][j]*mycos;
	    mymap[jj+2]+=tod->data[det][j]*mysin;
	    mymap[jj+3]+=tod->data[det][j]*mycos*mycos;
	    mymap[jj+4]+=tod->data[det][j]*mycos*mysin;
	    mymap[jj+5]+=tod->data[det][j]*mysin*mysin;

	  }
	}
      }
      break;
    default:  //should never get here, as error should have already triggered above.
      printf("Error - Don't know how to deal with map polarization in tod2polmap_copy.\n");
      break;
    }
#pragma omp critical(reduce_first)
    for (int i=0;i<npix;i++) {
      map->map[i]+=mymap[i];
    }
    if (npol>=2) {
#pragma omp critical(reduce_second)
      for (int i=npix;i<2*npix;i++) 
	map->map[i]+=mymap[i];
    }    
    if (npol>=3) {
#pragma omp critical(reduce_third)
      for (int i=2*npix;i<3*npix;i++) 
	map->map[i]+=mymap[i];
    }    
    if (npol>=4) {
#pragma omp critical(reduce_fourth)
      for (int i=3*npix;i<4*npix;i++) 
	map->map[i]+=mymap[i];
    }    
    if (npol>=5) {
#pragma omp critical(reduce_fifth)
      for (int i=4*npix;i<5*npix;i++) 
	map->map[i]+=mymap[i];
    }    
    if (npol>=6) {
#pragma omp critical(reduce_sixth)
      for (int i=5*npix;i<6*npix;i++) 
	map->map[i]+=mymap[i];
    }    
    
    free(mymap);
  }
}
/*--------------------------------------------------------------------------------*/

actData tod_times_map(const MAP *map, const mbTOD *tod, PARAMS *params)
{


  assert(map);
  assert(map->projection);

#if 0
  int nproc;
#pragma omp parallel shared(nproc) default(none)
#pragma omp single
  nproc=omp_get_num_threads();

  if (nproc*map->npix*sizeof(actData)>tod->ndata*tod->ndet*sizeof(int)) {
    //printf("doing index-saving projection.\n");
    tod2map_nocopy(map,tod,params);
    return;
  }
  //else
  // printf("doing old projection.\n");



  if (!map->have_locks) {
    //printf("setting up locks.\n");
    setup_omp_locks(map);
    //printf("set them up.\n");
    
  }
#endif
  actData tot=0;
#pragma omp parallel shared(tod,map,params) reduction(+:tot) default(none)
  { 

    //MAP *mymap=make_blank_map_copy(map);
    int *ind=(int *)malloc_retry(sizeof(int)*tod->ndata);
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);

#pragma omp for nowait
    for (int i=0;i<tod->ndet;i++) { 
      if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,params,i))) {
	get_pointing_vec_new(tod,map,i,ind,scratch);
	if (tod->kept_data) {
	  int row=tod->rows[i];
	  int col=tod->cols[i];
	  mbUncut *uncut=tod->kept_data[row][col];
	  for (int region=0;region<uncut->nregions;region++) 
	    for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
	      tot+=map->map[ind[j]]*tod->data[i][j];	  
	}
	else {
	  if (tod->uncuts) {
	    int row=tod->rows[i];
	    int col=tod->cols[i];
	    mbUncut *uncut=tod->uncuts[row][col];
	    for (int region=0;region<uncut->nregions;region++) 
	      for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
		tot+=map->map[ind[j]]*tod->data[i][j];
	  }
	  else 
	    for (int j=0;j<tod->ndata;j++) 
	      tot+=map->map[ind[j]]*tod->data[i][j];
	}
      }
    }
    
    free(ind);
    destroy_pointing_fit_scratch(scratch);
  } 

  return tot;
}
/*--------------------------------------------------------------------------------*/
void tod2mapset(MAPvec *maps, mbTOD *tod, PARAMS *params)
{
  for (int i=0;i<maps->nmap;i++)
    tod2map(maps->maps[i],tod,params);
}


/*--------------------------------------------------------------------------------*/
void tod2map_prealloc(MAP *map, mbTOD *tod, PARAMS *params)
/*use this version if each process already has its own mapset around.*/
{
  
  assert(1==0);  //don't think I'm using this.  
  int *ind=(int *)malloc_retry(sizeof(int)*tod->ndata);
  
  //#pragma omp parallel for
#pragma omp for schedule(dynamic,1) nowait
  for (int i=0;i<tod->ndet;i++) { 
    get_pointing_vec(tod,map,i,ind); 
    //printf("pixel limits are %d %d, out of %ld\n",ivecmin(ind,tod->ndata),ivecmax(ind,tod->ndata),mymap->npix);
    for (int j=0;j<tod->ndata;j++)
      map->map[ind[j]]+=tod->data[i][j];
  }
  free(ind);
}
/*--------------------------------------------------------------------------------*/

void tod2mapset_prealloc(MAPvec **maps, mbTOD *tod, PARAMS *params)
/* Each thread should already have its own private mapset to use this routine.
 The reduction *MUST* be done elsewhere (since this will be looped over lots of tods).*/
{
#pragma omp parallel shared(maps,tod,params)  
  {
    int myid=omp_get_thread_num();
    for (int i=0;i<maps[myid]->nmap;i++)
      tod2map_prealloc(maps[myid]->maps[i],tod,params);    
  }
}
/*--------------------------------------------------------------------------------*/


void calculate_avec(mbTOD *tod, int mydet, actData *avec)
     //calculate the noise here.  Simplistically, using a single known 1/f filter for all data.
{
  return;
}

/*--------------------------------------------------------------------------------*/
int cutvec2tod(mbTOD *tod, actData *cutvec)
{
  if (tod->cuts_fit_params) {
    //if cuts are paramaterized by polynomials instead of full cutvecs, behave properly.
    return cutpolys2tod(tod,cutvec);
  }
  if (!tod->have_data) {
    fprintf(stderr,"missing data in tod in cutvec2tod\n");
    return 1;
  }
  if (tod->cuts_as_vec==NULL) {
    fprintf(stderr,"Error, tod does not contain cuts_as_vec in cutvec2tod.\n");
    return 1;
  }
  if (tod->cuts_as_uncuts==NULL) {
    fprintf(stderr,"Error, tod does not contain cuts_as_uncut in cutvec2tod.\n");
    return 1;    
  }
  
  int ncopy=0;
  for (int det=0;det<tod->ndet;det++) 
    if (tod->cuts_as_vec[tod->rows[det]][tod->cols[det]]!=NULL) {
      mbUncut *cut=tod->cuts_as_vec[tod->rows[det]][tod->cols[det]];
      mbUncut *loc_cut=tod->cuts_as_uncuts[tod->rows[det]][tod->cols[det]];
      for (int i=0;i<cut->nregions;i++)  {
	memcpy(tod->data[det]+loc_cut->indexFirst[i],cutvec+cut->indexFirst[i],sizeof(actData)*(cut->indexLast[i]-cut->indexFirst[i])); 
	ncopy+=cut->indexLast[i]-cut->indexFirst[i];
      }
      
    }
  return 0;
}

/*--------------------------------------------------------------------------------*/
int tod2cutvec(mbTOD *tod, actData *cutvec)
{

  if (tod->cuts_fit_params) {
    //if cuts are paramaterized by polynomials instead of full cutvecs, behave properly.
    return tod2cutpolys(tod,cutvec);
  }
  //printf("greetings from tod2cutvec.\n");
  if (!tod->have_data) {
    fprintf(stderr,"missing data in tod in tod2cutvec\n");
    return 1;
  }
  if (tod->cuts_as_vec==NULL) {
    fprintf(stderr,"Error, tod does not contain cuts_as_vec in tod2cutvec.\n");
    return 1;
  }
  if (tod->cuts_as_uncuts==NULL) {
    fprintf(stderr,"Error, tod does not contain cuts_as_uncut in tod2cutvec.\n");
    return 1;    
  }
  
  for (int det=0;det<tod->ndet;det++) 
    if (tod->cuts_as_vec[tod->rows[det]][tod->cols[det]]!=NULL) {
      //printf("det is %d\n",det);
      mbUncut *cut=tod->cuts_as_vec[tod->rows[det]][tod->cols[det]];
      mbUncut *loc_cut=tod->cuts_as_uncuts[tod->rows[det]][tod->cols[det]];
      for (int i=0;i<cut->nregions;i++) 
	memcpy(cutvec+cut->indexFirst[i],tod->data[det]+loc_cut->indexFirst[i],sizeof(actData)*(cut->indexLast[i]-cut->indexFirst[i]));
      
    }
  return 0;
}
/*--------------------------------------------------------------------------------*/
int _get_npol_in_state(const int *pol_state)
{
#ifdef ACTPOL
  int tot=0;
  for (int i=0;i<MAX_NPOL;i++)
    if (pol_state[i])
      tot++;
  //if (tot==0)
  // tot=1;  //We probably meant to have TT
 
  return tot;
#else
  return 1;
#endif
  
  
}
/*--------------------------------------------------------------------------------*/
int get_npol_in_map(const MAP *map)
{
#ifdef ACTPOL
  const int *ptr=map->pol_state;//&(map->pol_state[0]);
  return _get_npol_in_state(ptr);
#else
  return 1;
#endif
}
/*--------------------------------------------------------------------------------*/
int get_map_poltag(const MAP *map)
{
#ifdef ACTPOL
  if (map->pol_state==NULL)
    return POL_I;



  const int *ptr=map->pol_state;
  
  int pol_i[MAX_NPOL]={1,0,0,0,0,0};
  if (memcmp(ptr,pol_i,sizeof(int)*MAX_NPOL)==0)
    return POL_I;
  
  int pol_iqu[MAX_NPOL]={1,1,1,0,0,0};
  if (memcmp(ptr,pol_iqu,sizeof(int)*MAX_NPOL)==0)
    return POL_IQU;

  int pol_qu[MAX_NPOL]={0,1,1,0,0,0};
  if (memcmp(ptr,pol_qu,sizeof(int)*MAX_NPOL)==0)
    return POL_QU;

  int pol_iqu_precon[MAX_NPOL]={1,1,1,1,1,1};
  if (memcmp(ptr,pol_iqu_precon,sizeof(int)*MAX_NPOL)==0)
    return POL_IQU_PRECON;

  int pol_qu_precon[MAX_NPOL]={0,0,0,1,1,1};
  if (memcmp(ptr,pol_qu_precon,sizeof(int)*MAX_NPOL)==0)
    return POL_QU_PRECON;
  
#else
  return POL_I;
#endif
  printf("Unknown polarization configuration in get_map_poltag.\n");
  return POL_ERROR;
}
/*--------------------------------------------------------------------------------*/
#ifdef ACTPOL
void set_map_polstate(MAP *map, int *pol_state)
{
  int old_npol=get_npol_in_map(map);
  int new_npol=_get_npol_in_state(pol_state);
  
  actData **new_map=matrix(new_npol,map->npix);
  memset(new_map[0],0,sizeof(actData)*new_npol*map->npix);
  if (map->map) {
    int old_ind=0;
    int new_ind=0;
    for (int i=0;i<MAX_NPOL;i++) {
      if ((map->pol_state[i])&&(pol_state[i])) 
	memcpy(new_map[new_ind],&(map->map[old_ind*map->npix]),map->npix*sizeof(actData));
      if (map->pol_state[i])
	old_ind++;
      if (pol_state[i])
	new_ind++;
    }
    free(map->map);
    map->map=new_map[0];
    free(new_map);  //only freeing the pointer to the pointers, keeping the actual address saved.
    memcpy(map->pol_state,pol_state,sizeof(pol_state[0])*MAX_NPOL);
  }
}
#endif
/*--------------------------------------------------------------------------------*/
void map2det(const MAP *map, const mbTOD *tod, actData *vec, int *ind, int det, PointingFitScratch *scratch)
//add a map into a vector.
{
  //get_pointing_vec(tod,map,det,ind);
  get_pointing_vec_new(tod,map,det,ind,scratch);
  //#ifdef ACTPOL
#if 0
  switch(map->pol_state) {
  case 0:
  case 1:
    for (int j=0;j<tod->ndata;j++) {
      vec[j]+=map->map[ind[j]];      
    }
    break;
  case 2:  //maps are only (Q,U)
    for (int j=0;j<tod->ndata;j++) {
      vec[j]+=map->map[2*ind[j]]*scratch->cos2gamma[j];
      vec[j]+=map->map[2*ind[j]+1]*scratch->sin2gamma[j];
    }
    break;
  case 3:  //maps are (I,Q,U)
    for (int j=0;j<tod->ndata;j++) {
      vec[j]+=map->map[3*ind[j]];
      vec[j]+=map->map[3*ind[j]+1]*scratch->cos2gamma[j];
      vec[j]+=map->map[3*ind[j]+2]*scratch->sin2gamma[j];      
      
    }
    break;
  default:
    assert(1==0);  //this means we had an unrecognized map type.
    break;
  }
#else
  for (int j=0;j<tod->ndata;j++) {
    vec[j]+=map->map[ind[j]];
  }
#endif

  
} 
/*--------------------------------------------------------------------------------*/
void map2det_scaled(const MAP *map, const mbTOD *tod, actData *vec, actData scale_fac, int *ind, int det, PointingFitScratch *scratch)
//add a map into a vector.
{
  //get_pointing_vec(tod,map,det,ind);
  get_pointing_vec_new(tod,map,det,ind,scratch);
  
  for (int j=0;j<tod->ndata;j++) {
    vec[j]+=map->map[ind[j]]*scale_fac;
  }
  
} 

/*--------------------------------------------------------------------------------*/
long is_tod_inbounds(const MAP *map, mbTOD *tod,const PARAMS *params)
{
  long nbad=0;
#pragma omp parallel shared(tod,map) reduction(+:nbad) default(none) 
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    int *ind=(int *)malloc_retry(sizeof(int)*tod->ndata);
    bool *inbounds=(bool *)malloc_retry(sizeof(bool)*tod->ndata);
    
#pragma omp for schedule(dynamic,1)
   for (int i=0;i<tod->ndet;i++) {
     if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
       get_map_projection_wchecks(tod,map,i,ind,scratch,inbounds);
       for (int j=0;j<tod->ndata;j++)
	 if (!inbounds[j])
	   nbad++;
     }
   }
   destroy_pointing_fit_scratch(scratch);   
   free(ind);
   free(inbounds);
  }
  return nbad;
}

/*--------------------------------------------------------------------------------*/
void save_tod_projection(const MAP *map, mbTOD *tod,const PARAMS *params)
{
  int **proj=imatrix(tod->ndet,tod->ndata);
#pragma omp parallel shared(tod,map,proj) default(none) 
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    //int *ind=(int *)malloc_retry(sizeof(int)*tod->ndata);
#pragma omp for schedule(dynamic,1)
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	get_pointing_vec_new(tod,map,i,proj[i],scratch);	
      }
      
    }
  }
  tod->pixelization_saved=proj;
}

/*--------------------------------------------------------------------------------*/
void map2tod(const MAP *map, mbTOD *tod,const PARAMS *params)
{
  
  assert(tod->have_data);

  actData scale_fac=1.0;
  if (params)
    scale_fac=*((actData *)params);

#pragma omp parallel shared(tod,map,stderr,scale_fac) default(none) 
 {
   PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
   int *ind=(int *)malloc_retry(sizeof(int)*tod->ndata);

#pragma omp for schedule(dynamic,1)
   for (int i=0;i<tod->ndet;i++) {
     //fprintf(stderr,"i is %d on %d\n",i,omp_get_thread_num());
     if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
       if (scale_fac==1)
	 map2det(map,tod,tod->data[i],ind,i, scratch);
       else
	 map2det_scaled(map,tod,tod->data[i],scale_fac,ind,i, scratch);
     }
   }
   destroy_pointing_fit_scratch(scratch);   
   free(ind);


 }
 
}
/*--------------------------------------------------------------------------------*/
void mapset2det(const MAPvec *maps, mbTOD *tod, const PARAMS *params, actData *vec, int *ind, int det, PointingFitScratch *scratch)
{
  //no asserts here for speed reasons - please make sure you've done them earlier.
  for (int map=0;map<maps->nmap;map++) 
    map2det(maps->maps[map],tod,vec,ind,det,scratch);

}
/*--------------------------------------------------------------------------------*/
void add_mapset2tod(const MAPvec *maps, mbTOD *tod, const PARAMS *params, actData fac)
//project a mapset into a TOD, scaled by factor fac (so to subtract, send in fac=-1)
{
#pragma omp parallel shared(maps,tod,params,fac) default(none) 
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    int *ind=(int *)malloc_retry(sizeof(int)*tod->ndata);
    assert(ind);
    actData *vec=vector(tod->ndata);
#pragma omp for schedule(dynamic,1) 
    for (int i=0;i<tod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
	
	memset(vec,0,sizeof(actData)*tod->ndata);
	mapset2det(maps,tod,params,vec,ind,i,scratch);
	for (int j=0;j<tod->ndata;j++)
	  tod->data[i][j]+=fac*vec[j];
      }
    }
    free(ind);
    free(vec);
    destroy_pointing_fit_scratch((PointingFitScratch *)scratch);
  }
}

/*--------------------------------------------------------------------------------*/
actData inline calc_srcamp(actData ra, actData dec, actData srcra, actData srcdec, const actData *beam, actData dtheta, int nbeam, actData cosdec)
{
  actData maxdist=dtheta*(nbeam-1);
  actData ddec=srcdec-dec;
  if (fabs(ddec)<maxdist) {
    actData dra=srcra-ra;
    if (fabs(dra)*cosdec<maxdist) {
#if 1
      actData ss=sin(0.5*ddec);
      ss=ss*ss;
      actData ss2=sin(0.5*(dra));
      actData cc=cos(dec)*cosdec*ss2*ss2;
      actData mydist=2*asin(sqrt(ss+cc));

      //now do a linear interpolation                                                                                                             
      if (mydist<maxdist) {
	int ibin=mydist/dtheta;
	actData frac=mydist-ibin*dtheta;
	return (beam[ibin]*(1-frac)+beam[ibin+1]*frac);
      }

#else
      actData ss=sin5(0.5*ddec);
      ss=ss*ss;
      actData ss2=sin5(0.5*(dra));
      actData cc=cos5(dec)*cosdec*ss2*ss2;
      actData x=sqrt(ss+cc);
      //actData mydist=2*asin(sqrt(ss+cc));
      actData mydist=2*(x+x*x*x/6);
      //now do a linear interpolation                                                                                                             
      if (mydist<maxdist) {
	int ibin=mydist/dtheta;
	actData frac=mydist-ibin*dtheta;
	return (beam[ibin]*(1-frac)+beam[ibin+1]*frac);
      }
#endif
    }

  }
  return 0;
}
/*--------------------------------------------------------------------------------*/
actData inline calc_srcamp_oversamp(int i, actData *ra, actData *dec,  actData srcra, actData srcdec, const actData *beam, actData dtheta, int nbeam, actData cosdec, int ndata, int oversamp)
{
  if ((oversamp<=1)||(i==0)||(i==ndata-1)) {
    return calc_srcamp(ra[i],dec[i],srcra,srcdec,beam,dtheta,nbeam,cosdec);
  }
  actData fac=0;
  actData ra1=0.5*(ra[i]+ra[i-1]);
  actData dec1=0.5*(dec[i]+dec[i-1]);
  actData ra2=0.5*(ra[i]+ra[i+1]);
  actData dec2=0.5*(dec[i]+dec[i+1]);
  actData dra=(ra2-ra1)/(oversamp-1);
  actData ddec=(dec2-dec1)/(oversamp-1);

  for (int j=0;j<oversamp;j++) {
    fac+=calc_srcamp(ra1+dra*j,dec1+ddec*j,srcra,srcdec,beam,dtheta,nbeam,cosdec);
  }
  return fac/oversamp;

}
/*--------------------------------------------------------------------------------*/
void add_src2tod(mbTOD *tod, actData ra, actData dec, actData src_amp, const actData *beam, actData dtheta, int nbeam, int oversamp)
//Add a source with amplitude src_amp at (ra,dec) into the timestreams in TOD
//eventually, oversamp will evaluate the beam at many places per sample
{
  if (!tod->have_data) {
    fprintf(stderr,"TOD mising data in add_src2tod.\n");
    return;
  }
  
  actData tpos=0,tdist=0;

#pragma omp parallel shared(tod,ra,dec,src_amp,beam,dtheta,nbeam,oversamp) reduction(+:tpos,tdist) default(none) 
  {
    actData maxdist=dtheta*(nbeam-1);
    actData maxra=maxdist/cos(dec);
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    actData tt;
    actData cosdec=cos(dec);
#pragma omp for schedule(dynamic,1)
    for (int det=0;det<tod->ndet;det++) {
      if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[det],tod->cols[det])) {
	tt=omp_get_wtime();
	get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
	tpos+=omp_get_wtime()-tt;
	
	tt=omp_get_wtime();
	for (int i=0;i<tod->ndata;i++) {
	  
	  actData fac=0;
	  if ((i==0)||(i==tod->ndata-1)||(oversamp<=1))
	    tod->data[det][i]+=src_amp*calc_srcamp(scratch->ra[i],scratch->dec[i],ra,dec,beam,dtheta,nbeam,cosdec);
	  else {
	    actData fac=0;	    
	    actData ra1=0.5*(scratch->ra[i]+scratch->ra[i-1]);
	    actData dec1=0.5*(scratch->dec[i]+scratch->dec[i-1]);
	    actData ra2=0.5*(scratch->ra[i]+scratch->ra[i+1]);
	    actData dec2=0.5*(scratch->dec[i]+scratch->dec[i+1]);
	    actData dra=(ra2-ra1)/(oversamp-1);
	    actData ddec=(dec2-dec1)/(oversamp-1);
	    
	    for (int j=0;j<oversamp;j++) {
	      fac+=calc_srcamp(ra1+dra*j,dec1+ddec*j,ra,dec,beam,dtheta,nbeam,cosdec);	      
	    }
	    tod->data[det][i]+=src_amp*fac/oversamp;
	  }
	}
	tdist+=omp_get_wtime()-tt;
	
      }
    }
    destroy_pointing_fit_scratch(scratch);    
  }
  printf("distance and position times were %8.4f %8.4f\n",tpos,tdist);
  
}
/*--------------------------------------------------------------------------------*/
int tod_hits_source(actData ra, actData dec, actData dist, mbTOD *tod)
{
  actData cosdec=cos(0.5*(tod->decmin+tod->decmax));
  if (tod->ramin>ra+dist/cosdec)
    return 0;
  if (tod->ramax<ra-dist/cosdec)
    return 0;
  if (tod->decmin>dec+dist)
    return 0;
  if (tod->decmax<dec-dist)
    return 0;
  return 1;

}
/*--------------------------------------------------------------------------------*/
void add_srcvec2tod(mbTOD *tod, actData *ra_in, actData *dec_in, actData *src_amp_in, int nsrc_in,const actData *beam, actData dtheta, int nbeam, int oversamp)
//Add a source with amplitude src_amp at (ra,dec) into the timestreams in TOD
//eventually, oversamp will evaluate the beam at many places per sample
{
  if (!tod->have_data) {
    fprintf(stderr,"TOD mising data in add_src2tod.\n");
    return;
  }
  
  actData *ra=vector(nsrc_in);
  actData *dec=vector(nsrc_in);
  actData *src_amp=vector(nsrc_in);

  int nsrc=0;
  
  for (int i=0;i<nsrc_in;i++) 
    if (src_amp_in[i]!=0)
      if (tod_hits_source(ra_in[i],dec_in[i],nbeam*dtheta,tod)) {
	ra[nsrc]=ra_in[i];
	dec[nsrc]=dec_in[i];
	src_amp[nsrc]=src_amp_in[i];
	nsrc++;
      }
  
  actData tpos=0,tdist=0;
  
  if (nsrc>0)     
#pragma omp parallel shared(tod,ra,dec,src_amp,beam,dtheta,nbeam,oversamp,nsrc) reduction(+:tpos,tdist) default(none) 
    {
      actData tt;

      actData maxdist=dtheta*(nbeam-1);
      actData maxra=maxdist/cos(0.5*(tod->decmin+tod->decmax));
      PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
      actData *cosdec=vector(nsrc);
      for (int i=0;i<nsrc;i++)
	cosdec[i]=cos(dec[i]);
#pragma omp for schedule(dynamic,1)
      for (int det=0;det<tod->ndet;det++) {
	if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[det],tod->cols[det])) {
	  tt=omp_get_wtime();
	  get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
	  tpos+=omp_get_wtime()-tt;


	  tt=omp_get_wtime();

	  for (int src=0;src<nsrc;src++) {
	    for (int i=0;i<tod->ndata;i++)
	      tod->data[det][i]+=src_amp[src]*calc_srcamp_oversamp(i,scratch->ra,scratch->dec,ra[src],dec[src],beam,dtheta,nbeam,cosdec[src],tod->ndata,oversamp); 	  
	  }
	}
      }
      free(cosdec);
      destroy_pointing_fit_scratch(scratch);    
    }
  //printf("distance and position times were %8.4f %8.4f\n",tpos,tdist);

  free(ra);
  free(dec);
  free(src_amp);
  
}

/*--------------------------------------------------------------------------------*/

void tod2srcvec(actData *src_amp_out,mbTOD *tod, actData *ra_in, actData *dec_in, int nsrc_in,const actData *beam, actData dtheta, int nbeam, int oversamp)
//Add a source with amplitude src_amp at (ra,dec) into the timestreams in TOD
//eventually, oversamp will evaluate the beam at many places per sample
{
  if (!tod->have_data) {
    fprintf(stderr,"TOD mising data in add_src2tod.\n");
    return;
  }
  
  actData *ra=vector(nsrc_in);
  actData *dec=vector(nsrc_in);
  int *indvec=(int *)malloc(sizeof(int)*nsrc_in);


  int nsrc=0;
  
  for (int i=0;i<nsrc_in;i++) 
    if (tod_hits_source(ra_in[i],dec_in[i],nbeam*dtheta,tod)) {
      ra[nsrc]=ra_in[i];
      dec[nsrc]=dec_in[i];
      indvec[nsrc]=i;
      nsrc++;
    }
  
  if (nsrc>0)     
#pragma omp parallel shared(tod,ra,dec,indvec,src_amp_out,beam,dtheta,nbeam,oversamp,nsrc) default(none) 
    {
      actData tt;
      actData *src_amp=vector(nsrc);
      memset(src_amp,0,sizeof(actData)*nsrc);


      actData maxdist=dtheta*(nbeam-1);
      actData maxra=maxdist/cos(0.5*(tod->decmin+tod->decmax));
      PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
      actData *cosdec=vector(nsrc);
      for (int i=0;i<nsrc;i++)
	cosdec[i]=cos(dec[i]);
#pragma omp for schedule(dynamic,1)
      for (int det=0;det<tod->ndet;det++) {
	if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[det],tod->cols[det])) {
	  get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
	  if (tod->kept_data) {
	    int row=tod->rows[det];
	    int col=tod->cols[det];
	    mbUncut *uncut=tod->kept_data[row][col];

	    for (int region=0;region<uncut->nregions;region++)
	      for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
		if (tod->data[det][j]!=0)
		  for (int src=0;src<nsrc;src++)
		    src_amp[src]+=tod->data[det][j]*calc_srcamp_oversamp(j,scratch->ra,scratch->dec,ra[src],dec[src],beam,dtheta,nbeam,cosdec[src],tod->ndata,oversamp);
	  }
	  else {
	    if (tod->uncuts) {
	      int row=tod->rows[det];
	      int col=tod->cols[det];	      
	      mbUncut *uncut=tod->uncuts[row][col];
	      for (int src=0;src<nsrc;src++)
		for (int region=0;region<uncut->nregions;region++)
		  for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++)
		    src_amp[src]+=tod->data[det][j]*calc_srcamp_oversamp(j,scratch->ra,scratch->dec,ra[src],dec[src],beam,dtheta,nbeam,cosdec[src],tod->ndata,oversamp);
	    }
	    else
	      for (int src=0;src<nsrc;src++)
		for (int j=0;j<tod->ndata;j++)
		  src_amp[src]+=tod->data[det][j]*calc_srcamp_oversamp(j,scratch->ra,scratch->dec,ra[src],dec[src],beam,dtheta,nbeam,cosdec[src],tod->ndata,oversamp);
	  }
	}
      }
      free(cosdec);
      destroy_pointing_fit_scratch(scratch);    
#pragma omp critical
      for (int src=0;src<nsrc;src++) {
	src_amp_out[indvec[src]]+=src_amp[src];
      }
    }
  
  free(ra);
  free(dec);
  free(indvec);
  
}

/*--------------------------------------------------------------------------------*/


void clear_tod(mbTOD *tod)
{
  pca_time tt;
  tick(&tt);
  assert(tod->have_data==1);
#pragma omp parallel for shared(tod)
  for (int i=0;i<tod->ndet;i++)
    memset(tod->data[i],0,sizeof(actData)*tod->ndata);
}
/*--------------------------------------------------------------------------------*/
void mapset2tod(MAPvec *maps, mbTOD *tod,PARAMS *params)
{
  clear_tod(tod);
  for (int i=0;i<maps->nmap;i++)
    map2tod(maps->maps[i],tod,params);
  if (params->remove_common)
    remove_common_mode(tod);
    
}
/*--------------------------------------------------------------------------------*/
void allocate_tod_storage(mbTOD *tod)
     /*put in storage space for the tod data*/
{
  if (tod->have_data) {
    fprintf(stderr,"Warning - attempting to reallocate data when space already exists.  Skipping.\n");
    return;
  }
  assert(tod->have_data==0);
  tod->data=matrix(tod->ndet,tod->ndata);
  assert(tod->data);
  tod->have_data=1;
}
/*--------------------------------------------------------------------------------*/
void free_tod_storage(mbTOD *tod)
     /*free temporary copy of the data*/
{
  if (!tod->have_data) {
    fprintf(stderr,"Warning - attempting to free non-existent data in TOD.  Skipping.\n");
    return;
  }
  assert(tod->have_data==1);
  free_matrix(tod->data);
  tod->have_data=0;
  tod->data=NULL;
}
/*--------------------------------------------------------------------------------*/
void rotate_matrix(actData **mat1, actData **mat2, bool do_transpose)
{
  
}
/*--------------------------------------------------------------------------------*/
void mapset2mapset(MAPvec *maps, TODvec *tods, PARAMS *params)
{
  MAPvec *maps_copy=make_mapset_copy(maps);
  clear_mapset(maps_copy);
#ifdef  MAPS_PREALLOC 
  int nproc;
  MAPvec **bigmaps;
#pragma omp parallel shared(nproc,bigmaps,maps) default(none)
  {
#pragma omp single
    {
      nproc=omp_get_num_threads();
      bigmaps=(MAPvec **)malloc_retry(nproc*sizeof(MAPvec *));
    }
    int myid=omp_get_thread_num();
    bigmaps[myid]=make_mapset_copy(maps);
    clear_mapset(bigmaps[myid]);
  }
#endif
  for (int i=0;i<tods->ntod;i++) {
    mbTOD *mytod=&(tods->tods[i]);
    allocate_tod_storage(mytod);
    mapset2tod(maps,mytod,params);
    if (!params->no_noise)
      filter_data(mytod);
#ifdef MAPS_PREALLOC
    tod2mapset_prealloc(bigmaps,mytod,params);
#else    
    tod2mapset(maps_copy,mytod,params);
#endif
    free_tod_storage(mytod);    
  }
#ifdef MAPS_PREALLOC
  for (int i=0;i<maps->nmap;i++)
    setup_omp_locks(maps_copy->maps[i]);
#pragma omp parallel shared(maps_copy,bigmaps,nproc) default(none)
  {
    int myid=omp_get_thread_num();
    for (int i=0;i<maps_copy->nmap;i++)
      omp_reduce_map(maps_copy->maps[i],bigmaps[myid]->maps[i]);
    destroy_mapset(bigmaps[myid]);
  }
  free(bigmaps);
#endif
  copy_mapset2mapset(maps,maps_copy);
  destroy_mapset(maps_copy);
#ifdef HAVE_MPI
  mpi_reduce_mapset(maps);
#endif
}
/*--------------------------------------------------------------------------------*/
void detrend_data(mbTOD *tod)
{
  assert(tod->have_data);
  if (0)
    return;  //kill-off the detrend.
  else {
#pragma omp parallel for shared(tod) default(none)
    for (int i=0;i<tod->ndet;i++) {
      actData a=tod->data[i][0];
      actData b=tod->data[i][tod->ndata-1];
      actData fac=(b-a)/(actData)(tod->ndata-1);
      for (int j=0;j<tod->ndata;j++) {
	//tod->data[i][j]=tod->data[i][j]-a-(b-a)*(actData)(j)/(actData)(tod->ndata-1);
	tod->data[i][j]=tod->data[i][j]-a-fac*(actData)(j);
      }
    }
  }
}
/*--------------------------------------------------------------------------------*/
void demean_data(mbTOD *tod)
{
  assert(tod->have_data);
#pragma omp parallel for shared(tod) default(none)
  for (int i=0;i<tod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i])) {
      actData mysum=0;
      actData myn=0;
      if (tod->uncuts) {
	  int row=tod->rows[i];
	  int col=tod->cols[i];
	  mbUncut *uncut=tod->uncuts[row][col];
	  for (int region=0;region<uncut->nregions;region++) 
	    for (int j=uncut->indexFirst[region];j<uncut->indexLast[region];j++) {
	      mysum+=tod->data[i][j];
	      myn++;
	    }
      }
      else {
	myn=tod->ndata;
	for (int j=0;j<tod->ndata;j++) {
	  mysum+=tod->data[i][j];
	}
      }
      actData mymean=mysum/myn;
      for (int j=0;j<tod->ndata;j++)
	tod->data[i][j]-=mymean;
    }
  }
}
/*--------------------------------------------------------------------------------*/
void cut_tod_ends(mbTOD *tod,actData tcut)
{
  int nsamp=tcut/tod->deltat;
  mbCutsExtendGlobal(tod->cuts,0,nsamp);
  mbCutsExtendGlobal(tod->cuts,tod->ndata-nsamp,tod->ndata);
  //tod->n_to_window=nsamp;
  set_tod_window(tod,tcut);
}
/*--------------------------------------------------------------------------------*/
void set_tod_window(mbTOD *tod,actData tcut)
{
  int nsamp=tcut/tod->deltat;
  tod->n_to_window=nsamp;
}
/*--------------------------------------------------------------------------------*/
void window_data(mbTOD *tod)
{
  assert(tod->have_data);
  int nsamp=tod->n_to_window;
  if (nsamp<=0)
    return;
  actData *window=vector(nsamp);
  actData *window2=vector(nsamp);
  for (int i=0;i<nsamp;i++) {
    window[i]=0.5-0.5*cos(M_PI*i/(nsamp+0.0));
    window2[i]=0.5+0.5*cos(M_PI*i/(nsamp+0.0));
  }
  int nn=tod->ndata-nsamp;
#pragma omp parallel for shared(tod,window,nsamp)
  for (int i=0;i<tod->ndet;i++) {
    for (int j=0;j<nsamp;j++)
      tod->data[i][j]*=window[j];
    for (int j=0;j<nsamp;j++)
      tod->data[i][nn+j]*=window2[j];    
  }  

  free(window);
  free(window2);
}
/*--------------------------------------------------------------------------------*/
void unwindow_data(mbTOD *tod)
{
  assert(tod->have_data);
  int nsamp=tod->n_to_window;
  if (nsamp<=0)
    return;
  actData *window=vector(nsamp);
  actData *window2=vector(nsamp);
  for (int i=0;i<nsamp;i++) {
    window[i]=0.5-0.5*cos(M_PI*i/(nsamp+0.0));
    window2[i]=0.5+0.5*cos(M_PI*i/(nsamp+0.0));
  }
  int nn=tod->ndata-nsamp;
#pragma omp parallel for shared(tod,window,nsamp)
  for (int i=0;i<tod->ndet;i++) {
    for (int j=1;j<nsamp;j++)
      tod->data[i][j]/=window[j];
    for (int j=0;j<nsamp-1;j++)
      tod->data[i][nn+j]/=window2[j];    
  }  

  free(window);
  free(window2);
}
/*--------------------------------------------------------------------------------*/
void remove_common_mode(mbTOD *tod)
{
  assert(tod->have_data);
  detrend_data(tod);

  mbNoiseCommonMode *common=mbAllocateCommonMode(tod,2,2); //eventually, these parameters need to actually get set.
  common->nsig=50;  //don't want to find a calbol for now...
  mbRescaleArray(tod,common);
  mbCalculateCommonMode(common);
  mbFindCalbols(common);
  mbCalculateCommonModeVecs(common);
  mbCalculateCommonModeFitParams(common);  
  mbApplyCommonMode(tod,common,false);
  
  mbCalculateCommonFracErrs(tod, common);
  nkCutUncorrDets(tod,common,2.0);
  nkCutUnsmoothDets(tod, common,2.5,2.0);
  mbNoiseCommonModeFree( common);
  psFree(common);
}
/*--------------------------------------------------------------------------------*/
void assign_tod_value(mbTOD *tod, actData val)
{
  assert(tod->have_data);
#pragma omp parallel for shared(tod,val) default(none)  
  for (int i=0;i<tod->ndet;i++)
    for (int j=0;j<tod->ndata;j++)
      tod->data[i][j]=val;
}
/*--------------------------------------------------------------------------------*/
int get_weights(MAPvec *maps, TODvec *tods, PARAMS *params)
{
  clear_mapset(maps);
  for (int i=0;i<tods->ntod;i++) {
    mbTOD *mytod=&(tods->tods[i]);
    allocate_tod_storage(mytod);
    assign_tod_value(mytod,1.0);
    tod2mapset(maps,mytod,params); 
    free_tod_storage(mytod);
  }
#ifdef HAVE_MPI
  mpi_reduce_mapset(maps);
#endif

  return 0;
}
/*--------------------------------------------------------------------------------*/
int make_initial_mapset(MAPvec *maps, TODvec *tods,PARAMS *params)
{
  mprintf(stdout,"Making initial mapset.\n");
  mprintf(stdout,"Initial limits are %10.5f %10.5f %12.6f\n",maps->maps[0]->ramin,maps->maps[0]->decmin,maps->maps[0]->pixsize);


  bool is_blank=is_mapset_blank(maps);
  MAPvec *maps_in;
  if (!is_blank)
    maps_in=make_mapset_copy(maps);


  clear_mapset(maps);
  MAP simmap; //use this guy if we are simulating data internally
  if (params->do_sim) {
    readwrite_simple_map(&simmap,params->inname,DOREAD);
    mprintf(stdout,"simmap lims are %10.5f %10.5f %10.6f\n",simmap.ramin,simmap.decmin,simmap.pixsize);
  }

#ifdef HAVE_MPI
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif


  if (is_blank)
    mprintf(stdout,"maps are blank inside initial mapset.\n");
  else
    mprintf(stdout,"maps are not blank inside initial mapset.\n");
  for (int i=0;i<tods->ntod;i++) {
    mbTOD *mytod=&(tods->tods[i]);
    mprintf(stdout,"working on %d %s\n",i,mytod->dirfile);
    //keep_1_det(mytod,10,10); //make sure to get rid of this!!!    
    if (params->do_sim) {

      allocate_tod_storage(mytod);
      if (params->do_blank)
	assign_tod_value(mytod,0.0);
      else {

	mprintf(stdout,"simulating now...\n");
	map2tod(&simmap,mytod,params);
	mprintf(stdout,"Finished.\n");
	int ii=ind_from_rowcol(mytod,15,15);
	int myid=0;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
	if (myid==0) {
	  //readwrite_simple_map(&simmap,"crappy_test.map",DOWRITE);
	  write_timestream(mytod->data[ii],mytod->ndata,"data_sim.dat");
	}
	//printf("On %d, min/max are %14.5e %14.5e\n",myid,vecmin(simmap.map,simmap.npix),vecmax(simmap.map,simmap.npix));	
	//exit(EXIT_SUCCESS);

      }
    }
    else {
      if (params->do_blank) {
	allocate_tod_storage(mytod);
	assign_tod_value(mytod,0.0);
      }
      else
	read_tod_data(mytod);

    }
    
    
    if (params->add_noise) {
      mprintf(stdout,"adding noise.\n");
      
#if 1
      set_tod_noise(mytod,1.2e-3,1,-1.5);
      add_noise_to_tod_new(mytod);
#else
      add_noise_to_tod(mytod,params);
      detrend_data(mytod);
#endif
      
    }
    
    if ((params->remove_common)||(params->remove_mean)) {
      mprintf(stdout,"removing common mode.\n");      
      remove_common_mode(mytod);
    }


    if (1)    {
      int ii=ind_from_rowcol(mytod,15,15);
      if (ii>=0) {
	write_timestream(mytod->data[ii],mytod->ndata,"data_presub.dat");
	if (!is_blank)
	  add_mapset2tod(maps_in,mytod,params,-1.0);
	write_timestream(mytod->data[ii],mytod->ndata,"data_postsub.dat");
      }
    }
    
    
    
    if (0)
      nkDeButterworth(mytod);

    if (params->deglitch) {
      pca_time tt;
      tick(&tt);
      glitch_all_detectors_simple(mytod,false,true,false,5.0,0.1,2.0,1);
      mprintf(stdout,"took %8.4f seconds to glitch everybody.\n",tocksilent(&tt));
    }
    pca_time tt;
    tick(&tt);

    bool always_filter=false; //cheezy hack for testing.

    if ((!params->no_noise)||(always_filter)) {
      mprintf(stdout,"Fitting noise.\n");


      mytod->noise=nkFitTODNoise(mytod, MBNOISE_LINEAR_POWLAW,0.5,80,-1.0);          
      cut_badnoise_dets(mytod);
      //destroy_tod_noise(mytod);
      


      //set_tod_noise(mytod,1.2e-3,1,-1.5);  //REMOVE!!!!!!
      //nkAssignNoiseKnee( mytod,4.0);  //if knee < 4, set it to 4.

      //mprintf(stdout,"Took %8.3f seconds to fit noise.\n",tocksilent(&tt));

    }
  

#if 1
    int ngood=0;
    for (int j=0;j<mytod->ndet;j++) {
      if (is_det_used(mytod,j)) {
	ngood++;
	//printf("keeping detector %d %d\n",mytod->rows[j],mytod->cols[j]);
	// if (mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[j],mytod->cols[j]))
	//printf("detector %3d %3d %3d is cut\n",j,mytod->rows[j],mytod->cols[j]);
      }
    }
    mprintf(stdout,"keeping a total of %d detectors.\n",ngood);
#endif
    if ((!params->no_noise) || always_filter) {
      mprintf(stdout,"filtering data.\n");
      filter_data(mytod);
    }

    mprintf(stdout,"calling tod2mapset.\n");
  
    tod2mapset(maps,mytod,params); 
    mprintf(stdout,"finished.\n");
    free_tod_storage(mytod);
  }

    
  if (params->do_sim)
    free(simmap.map);
#ifdef HAVE_MPI
  mpi_reduce_mapset(maps);
#endif

  if (strlen(params->rawname)) {
    readwrite_simple_map(maps->maps[0],params->rawname,DOWRITE);
  }
  
  if (!is_blank)
    destroy_mapset(maps_in);

  return 0;
}
/*--------------------------------------------------------------------------------*/
void map_axpy(MAP *y, MAP *x, actData a)
{
  assert(x->npix==y->npix);
#pragma omp parallel for shared(x,y,a) default(none)  
  for (int i=0;i<x->npix;i++) {
    y->map[i]=y->map[i]+x->map[i]*a;    
  }
}
/*--------------------------------------------------------------------------------*/
void mapset_axpy(MAPvec *y, MAPvec *x, actData a)
{
  assert(x->nmap==y->nmap);
  for (int i=0;i<x->nmap;i++)
    map_axpy(y->maps[i],x->maps[i],a);
}
/*--------------------------------------------------------------------------------*/
actData map_times_map(MAP *x, MAP *y)
{
  assert(x->npix==y->npix);
  double tot=0;
#pragma omp parallel for shared(x,y) reduction(+:tot) default(none)
  for (int i=0;i<x->npix;i++)
    tot += x->map[i]*y->map[i];

  return (actData)tot;
}
/*--------------------------------------------------------------------------------*/
actData mapset_times_mapset(MAPvec *x, MAPvec *y)
{
  assert(x->nmap==y->nmap);
  actData tot=0;
  for (int i=0;i<x->nmap;i++)
    tot += map_times_map(x->maps[i],y->maps[i]);
  return tot;
}
/*--------------------------------------------------------------------------------*/
void copy_map2map(MAP *map2, MAP *map)
{
  assert(map->npix==map2->npix);
  assert(map->npix>0);
  memcpy(map2->map,map->map,sizeof(actData)*map->npix*get_npol_in_map(map));
}
/*--------------------------------------------------------------------------------*/
void copy_mapset2mapset(MAPvec *map2, MAPvec *map)
{
  assert(map2->nmap==map->nmap);
  for (int i=0;i<map->nmap;i++) {
    copy_map2map(map2->maps[i],map->maps[i]);
  }
}
/*--------------------------------------------------------------------------------*/
actData PCGstep(MAPvec *r, MAPvec *p, MAPvec *x, TODvec *tods, MAPvec *wts, PARAMS *params)
{
  pca_time tt;
  tick(&tt);
  MAPvec *ap=make_mapset_copy(p);
  MAPvec *mr=make_mapset_copy(r);
  mapset2mapset(ap,tods,params);
  apply_preconditioner(mr,wts,params);

  actData rsqr=mapset_times_mapset(r,mr);
  destroy_mapset(mr);
  actData pap=mapset_times_mapset(p,ap);
  actData alpha_k=rsqr/pap;
  mapset_axpy(x,p,alpha_k);
  MAPvec *rk=make_mapset_copy(r);
  mapset_axpy(rk,ap,-alpha_k);  
  destroy_mapset(ap);

  MAPvec *mrk=make_mapset_copy(rk);
  apply_preconditioner(mrk,wts,params);

  actData beta_k=mapset_times_mapset(rk,mrk)/rsqr;
  MAPvec *pk=make_mapset_copy(mrk);
  mapset_axpy(pk,p,beta_k);

  //int mpi_broadcast_mapset(MAPvec *maps,master)
  
  copy_mapset2mapset(r,rk);
  copy_mapset2mapset(p,pk);
  //MPI_Barrier(MPI_COMM_WORLD);
  //mpi_broadcast_mapset(r,0);
  //mpi_broadcast_mapset(x,0);
  //mpi_broadcast_mapset(p,0);



  destroy_mapset(rk);
  destroy_mapset(pk);
  destroy_mapset(mrk);
  
  //printf("Iteration took %8.3f seconds.\n",tocksilent(&tt));
  mprintf(stdout,"Iteration took %8.3f seconds.\n",tocksilent(&tt));
  
  return rsqr;
  
}
/*--------------------------------------------------------------------------------*/
actData PCGstep_noprecon(MAPvec *r, MAPvec *p, MAPvec *x, TODvec *tods, MAPvec *wts, PARAMS *params)
{  
  pca_time tt;
  tick(&tt);
  MAPvec *ap=make_mapset_copy(p);
  apply_preconditioner(ap,wts,params);
  mapset2mapset(ap,tods,params);
  apply_preconditioner(ap,wts,params);
  actData rsqr=mapset_times_mapset(r,r);
  actData pap=mapset_times_mapset(p,ap);
  actData alpha_k=rsqr/pap;
  MAPvec *rk=make_mapset_copy(r);
  mapset_axpy(rk,ap,-alpha_k);
  actData beta_k=mapset_times_mapset(rk,rk)/rsqr;
  MAPvec *pk=make_mapset_copy(rk);
  mapset_axpy(pk,p,beta_k);
  mapset_axpy(x,p,alpha_k);
  //int mpi_broadcast_mapset(MAPvec *maps,master)
  
  copy_mapset2mapset(r,rk);
  copy_mapset2mapset(p,pk);
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  mpi_broadcast_mapset(r,0);
  mpi_broadcast_mapset(x,0);
  mpi_broadcast_mapset(p,0);
#endif


  destroy_mapset(rk);
  destroy_mapset(pk);
  destroy_mapset(ap);

  //printf("Iteration took %8.3f seconds.\n",tocksilent(&tt));
  mprintf(stdout,"Iteration took %8.3f seconds.\n",tocksilent(&tt));

  return rsqr;

}
/*--------------------------------------------------------------------------------*/
size_t freadwrite(void *ptr, size_t sz, long nobj, FILE *stream, int dowrite)
{
  if (dowrite==DOWRITE)
    return fwrite(ptr,sz,nobj,stream);
  else
    return fread(ptr,sz,nobj,stream);
}
/*--------------------------------------------------------------------------------*/
void apply_preconditioner( MAPvec *maps,MAPvec *weights,PARAMS *params)
{
  if (params->precondition) {
    MAP *map=maps->maps[0];
    MAP *wt=weights->maps[0];
#pragma omp parallel for shared(map,wt,params) default(none)
    for (int i=0;i<map->npix;i++) {
      if (wt->map[i]>0)
	map->map[i]/=wt->map[i];
    }    
  }
}
/*--------------------------------------------------------------------------------*/
void dump_data(mbTOD *tod, char *fname)
{
  FILE *outfile=fopen_safe(fname,"w");
  assert(tod->have_data);
  assert(outfile);
  fwrite(&tod->ndet,1,sizeof(int),outfile);
  fwrite(&tod->ndata,1,sizeof(int),outfile);
  fwrite(tod->data[0],tod->ndet*tod->ndata,sizeof(actData),outfile);
  fclose(outfile);
}
/*--------------------------------------------------------------------------------*/
void readwrite_simple_map(MAP *map, char *filename, int dowrite)
{
#ifdef HAVE_MPI
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  //if ((dowrite==DOREAD)||((dowrite==DOWRITE)&&(myid==0))) {  //if we're writing, only the master writes.
  if ((dowrite==DOREAD)||(myid==0)) {  //if we're writing, only the master writes.
#else
    if (1) {
#endif

      if (dowrite==DOWRITE)
	if (map->projection)
	  if (map->projection->proj_type==NK_CEA) {
	    //write_cea_map(map,filename);
	    return;
	  }

      FILE *iofile;
      if (dowrite==DOWRITE) {
	mprintf(stdout,"trying to write to %s\n",filename);
	//iofile=fopen(filename,"w");
	iofile=fopen_safe(filename,"w");
	
      }
      else {
	mprintf(stdout,"trying to read from %s\n",filename);
	//iofile=fopen(filename,"r");
	iofile=fopen_safe(filename,"r");    
      }
      
      assert(iofile);
      freadwrite(&map->nx,sizeof(int),1,iofile,dowrite);
      freadwrite(&map->ny,sizeof(int),1,iofile,dowrite);
      map->npix=map->nx*map->ny;
      freadwrite(&map->pixsize,sizeof(actData),1,iofile,dowrite);
      freadwrite(&map->ramin,sizeof(actData),1,iofile,dowrite);
      freadwrite(&map->ramax,sizeof(actData),1,iofile,dowrite);
      freadwrite(&map->decmin,sizeof(actData),1,iofile,dowrite);
      freadwrite(&map->decmax,sizeof(actData),1,iofile,dowrite);
      //mprintf(stdout,"limits on %s are %10.5f %10.5f %10.5f %10.5f %4d %4d %10.5f\n",filename,map->decmin,map->decmax,map->ramin,map->ramax,map->nx,map->ny,map->pixsize);
      if (dowrite==DOREAD)
	map->map=vector(map->npix);
      freadwrite(map->map,sizeof(actData),map->npix,iofile,dowrite);
      fclose(iofile);
      //printf("npix is %ld\n",map->npix);
    }
  
}
/*--------------------------------------------------------------------------------*/
void run_PCG(MAPvec *maps, TODvec *tods, PARAMS *params)
{

  //createFFTWplans(tod);

  bool had_maps;
  MAPvec *maps_in;
  had_maps=!(is_mapset_blank(maps));
  if (had_maps) 
    maps_in=make_mapset_copy(maps);  //save 'em, since the incoming mapset gets wiped over in make_initial_mapset

  mprintf(stdout,"Making initial mapset.\n");

  make_initial_mapset(maps,tods,params);
  
  if (params->write_pointing) {
    mprintf(stdout,"Writing pointing solution to disk.\n");
    for (int i=0;i<tods->ntod;i++) {
      char fname[512];
#ifdef HAVE_MPI
      int myid;
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      sprintf(fname,"tod_pointing_%d_%d.dat",i,myid);
#else
      sprintf(fname,"tod_pointing_%d.dat",i);
#endif
      write_tod_pointing_to_disk(&tods->tods[i],fname);
    }
    
  }

  MAPvec *weights=make_mapset_copy(maps);
  get_weights(weights,tods,params);
  char wtname[MAXLEN];
  sprintf(wtname,"%s.weights",params->outname);
  readwrite_simple_map(weights->maps[0],wtname,DOWRITE);
  if (params->rawonly) 
    exit(EXIT_SUCCESS);


  //readwrite_simple_map(weights->maps[0],"weights.dat",DOWRITE);
  //apply_preconditioner(maps,weights,params);


  mprintf(stdout,"making r,p, and x\n");
  MAPvec *r=make_mapset_copy(maps);
  MAPvec *p=make_mapset_copy(maps); 

  apply_preconditioner(p,weights,params);
  MAPvec *x=make_mapset_copy(maps);
  //clear_mapset(r);
  //clear_mapset(p);
  clear_mapset(x);
  //pca_time tt;
  
  actData residual=1e20;
  int iter=0;
  int converged=0;
  actData first_residual=0;
  while ((iter<params->maxiter)&&(converged==0))
    {
      iter++;
      //tick(&tt);
      residual=PCGstep(r,p,x,tods,weights,params);
      if (iter==1) 
	first_residual=residual;

      if (residual<params->tol*first_residual)
	converged=1;
      mprintf(stderr,"residual is %14.5e at iteration %d.\n",residual,iter);
      //fprintf(stderr,"residual is %14.5e at iteration %d.\n",residual,iter);
      //fprintf(stderr,"residual is %14.5e at iteration %d.  Step took %8.3f seconds.\n",residual,iter,tocksilent(&tt));            
#ifdef HAVE_MPI
      int ierr,myid;
      ierr=MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      if (myid==0) {
	char outname[512];
	sprintf(outname,"%s_%d.out",params->tempname,iter);
	//sprintf(outname,"temporary_map_commonsub_%d.out",iter);
	readwrite_simple_map(x->maps[0],outname,DOWRITE);
      }
      
#else
      //displayMap(x->maps[0]);
#endif
    }
  copy_mapset2mapset(maps,x);
  destroy_mapset(x);
  destroy_mapset(r);
  destroy_mapset(p);

  if (had_maps) {
    //add mapset back in.
    mapset_axpy(maps,maps_in,1.0);
    destroy_mapset(maps_in);
    
  }
}
/*--------------------------------------------------------------------------------*/
void print_options(PARAMS *params)
{
  if (params->do_sim) {
    printf("Going to simulate data internally from map %s\n",params->inname);
    if (params->do_blank)
      printf("   Scratch that, going to use zero for the input data.\n");
  }
  else
    printf("Going to read data from disk.\n");
  if (params->remove_common)
    printf("Going to remove common mode in mapping.\n");
  else
    printf("Not going to remove common mode.\n");
  if (params->no_noise)
    printf("going to skip noise filtering %d.\n",params->no_noise);
  else
    printf("noise filtering is enabled %d.\n",params->no_noise);
  if (params->add_noise)
    printf("going to add noise to simulated data.\n");
  else
    printf("not going to add noise to simulated data.\n");
  printf("Temporary map name base is %s\n",params->tempname);
  printf("Output name base is %s\n",params->outname);
  for (int i=0;i<params->ntod;i++) {
    printf("dirfile %4d is %s\n",i,params->datanames[i]);
  }
  printf("Map pixel size is %6.3g arcmin.\n",params->pixsize);
  if (params->precondition)
    printf("Going to use the weights as a preconditioner.\n");
  else
    printf("Not going to use a preconditioner.\n");
  if (params->use_input_limits)
    printf("Going to use same shape as input map for output map.\n");
  if (params->remove_mean)
    printf("going to remove the common mode from the input data.\n");
  if (params->deglitch)
    printf("going to deglitch data.\n");
  else
    printf("not going to deglitch data.\n");

  if (params->rawonly)
    printf("going to quit after making the raw(dirty) map and weights.\n");
  else
    printf("not going to quit after making the raw(dirty) map and weights.\n");

  
  if (params->write_pointing)
    printf("Going to write ra/dec solutions to disk for all tods.\n");
  else
    printf("Not going to write ra/dec solutions to disk.\n");
  
  printf("Going to run to a fractional tolerance of %14.5g or %d iterations, whichever comes first.\n",params->tol,params->maxiter);
  
}

/*--------------------------------------------------------------------------------*/
void parse_params(int argc, char **argv, PARAMS *params)
//parse the parameter file.
{
  char *tok;
  int *found_list=(int *)calloc(argc,sizeof(int));


  if (tok=find_argument(argc,argv,"@tol",found_list)) {
    params->tol=atof(tok);
    printf("setting tolerance to %12.4e\n",params->tol);
  }
  if (exists_in_command_line(argc,argv,"@precon",found_list)) {
    params->precondition=true;
    printf("Enabling preconditioning.\n");    
  }
  if (exists_in_command_line(argc,argv,"@no_noise",found_list)) {
    printf("Not going to use noise filtering.\n");
    params->no_noise=true;
  }
  if (exists_in_command_line(argc,argv,"@quit",found_list)) {
    printf("Going to read script and quit.\n");
    params->quit=true;
  }
  if (exists_in_command_line(argc,argv,"@rawonly",found_list)) {
    printf("Going to make initial map and quit.\n");
    params->rawonly=true;
  }
  if (tok=find_argument(argc,argv,"@input",found_list)) {
    strncpy(params->inname,tok,MAXLEN);
    printf("Input map is %s\n",params->inname);
  }
  if (tok=find_argument(argc,argv,"@output",found_list)) {
    strncpy(params->outname,tok,MAXLEN);
    printf("Output map is %s\n",params->outname);
  }

  if (tok=find_argument(argc,argv,"@rawmap",found_list)) {
    strncpy(params->rawname,tok,MAXLEN);
    printf("Output raw map is %s\n",params->rawname);
  }

  char **file_argv=get_list_from_argv(argc,argv,"@data",&(params->ntod),found_list);
  assert(params->ntod>0);  //gotta have data to do anything!
  for (int i=0;i<params->ntod;i++) {
    assert(strlen(file_argv[i])<MAXLEN-1);
    strncpy(params->datanames[i],file_argv[i],MAXLEN-1);
  }
  free_argv(params->ntod,file_argv);

  if (tok=find_argument(argc,argv,"@temp",found_list)) {
    strncpy(params->tempname,tok,MAXLEN);
    printf("Temporary output maps are %s\n",params->tempname);
  }
  if (tok=find_argument(argc,argv,"@seed",found_list)) {
    params->seed=atoi(tok);
    printf("Random seed is %ld\n",params->seed);
  }
  if (tok=find_argument(argc,argv,"@maxiter",found_list)) {
    params->maxiter=atoi(tok);
    printf("Maximum number of map-making iterations is %d\n",params->maxiter);
  }
  
  if (tok=find_argument(argc,argv,"@maxtod",found_list)) {
    params->maxtod=atoi(tok);
    printf("Maximum number of allowed TOD's is %d\n",params->maxtod);
  }

  if (exists_in_command_line(argc,argv,"@common",found_list)) {
    printf("Going to remove common mode.\n");
    params->remove_common=true;
  }

  if (exists_in_command_line(argc,argv,"@mean",found_list)) {
    printf("Going to remove common mode from input data.\n");
    params->remove_mean=true;
  }

  if (tok=find_argument(argc,argv,"@pixsize",found_list)) {
    params->pixsize=atof(tok);
    printf("Output map pixel size is %12.4e\n",params->pixsize);
    params->pixsize *= M_PI/60.0/180.0;  //turn it into radians.
  }

  if (exists_in_command_line(argc,argv,"@sim",found_list)) {
    params->do_sim=true;
    printf("Going to run from simulated data.\n");
  }

  if (exists_in_command_line(argc,argv,"@add_noise",found_list)) {
    params->add_noise=true;
    printf("Going to add noise to the data.\n");
  }
  
  if (exists_in_command_line(argc,argv,"@blank",found_list)) {
    params->do_blank=true;
    printf("Input data will be blank.\n");
  }
  if (exists_in_command_line(argc,argv,"@lims",found_list)) {
    params->use_input_limits=true;
    printf("output map will have same pixelization as input map.\n");
  }
  
  if (exists_in_command_line(argc,argv,"@deglitch",found_list)) {
    params->deglitch=true;
    printf("going to deglitch data.\n");
  }

  if (tok=find_argument(argc,argv,"@pointing_offsets",found_list)) {
    strncpy(params->pointing_file,tok,MAXLEN-1);
    printf("pointing offset will be read from %s\n",params->pointing_file);
  }

  if (tok=find_argument(argc,argv,"@altaz_file",found_list)) {
    strncpy(params->altaz_file,tok,MAXLEN-1);
    printf("going to read alt/az/ctime positions from %s\n",params->altaz_file);
  }


  if (params->use_rows=get_int_list_from_argv(argc,argv,"@use_rows",&(params->n_use_rows),found_list)) {
    printf("only mapping rows: ");
    for (int i=0;i<params->n_use_rows;i++)
      printf(" %d ",params->use_rows[i]);
    printf("\n");
  }

  if (params->use_cols=get_int_list_from_argv(argc,argv,"@use_cols",&(params->n_use_cols),found_list)) {
    printf("only mapping cols: ");
    for (int i=0;i<params->n_use_cols;i++)
      printf(" %d ",params->use_cols[i]);
    printf("\n");
  }


  if (exists_in_command_line(argc,argv,"@write_pointing",found_list)) {
    params->write_pointing=true;
    printf("going to write pointing solutions to disk.\n");
  }



  
  for (int i=0;i<argc;i++)
    if (!found_list[i]) {
      printf("Ignored input argument %d, which was %s\n",i,argv[i]);
    }
  free(found_list);
  

}
/*--------------------------------------------------------------------------------*/
int get_parameters(int argc, char *argv[], PARAMS *params)
{

#ifdef HAVE_MPI
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank==0) 
#else
    if (1) 
#endif
      {
	params->do_sim=false;
	params->do_blank=false;
	sprintf(params->inname,"simple_map_input.dat");
	sprintf(params->outname,"map_fit.dat");
	sprintf(params->tempname,"temporary_map_");
	params->remove_common=false;
	params->maxiter=100;
	params->tol=1e-4;
	params->no_noise=false;
	params->add_noise=false;
	params->pixsize=0.5;
	params->quit=false;
	params->precondition=false;
	params->use_input_limits=false;
	params->seed=1;
	params->remove_mean=false;
	params->maxtod=0;  //0 for unlimited.
	params->deglitch=false;
	params->rawonly=false;

	int myargc;
	char **myargv;
	char *line_in;
	printf("argc is %d\n",argc);
	switch(argc) {
	case 1:
	  printf("Reading from stdin.\n");
	  line_in=read_all_stdin(NULL);
	  printf("line_in is %s\n",line_in);
	  myargv=create_argv_new(line_in,&myargc," \n");
	  free(line_in);
	  break;
	case 2:
	  line_in=read_all_stdin(argv[1]);
	  myargv=create_argv_new(line_in,&myargc," \n");
	  free(line_in);
	  break;
	default:
	  myargv=argv;
	  myargc=argc;
	  break;
	}
	
	parse_params(myargc,myargv,params);
	print_options(params);
	free_argv(myargc,myargv);
	//exit(EXIT_SUCCESS);	
      }
#ifdef HAVE_MPI
  MPI_Bcast(params,sizeof(PARAMS),MPI_CHAR,0,MPI_COMM_WORLD);
#endif
  return 0;
  
}
/*--------------------------------------------------------------------------------*/
void keep_1_det(mbTOD *tod,int row, int col)
{
  assert(!mbCutsIsAlwaysCut(tod->cuts,row,col));
  for (int i=0;i<tod->ndet;i++) {
    if ((tod->rows[i]!=row)||(tod->cols[i]!=col))
      mbCutsSetAlwaysCut(tod->cuts, tod->rows[i], tod->cols[i]);
  }
  
  
}

/*--------------------------------------------------------------------------------*/
void get_data_corrs(mbTOD *tod)
{
  assert(tod->have_data);
  if (tod->corrs==NULL) {
    tod->corrs=matrix(tod->ndet,tod->ndet);
    memset(tod->corrs[0],0,tod->ndet*tod->ndet*sizeof(actData));
  }
#ifdef ACTDATA_DOUBLE
  //printf("Calling dsyrk.\n");
  clapack_dsyrk('u','t',tod->ndet,tod->ndata,1.0,tod->data[0],tod->ndata,0.0,tod->corrs[0],tod->ndet);
  //printf("called.\n");
#else
  clapack_ssyrk('u','n',tod->ndet,tod->ndata,1.0,tod->data[0],tod->ndata,0.0,tod->corrs[0],tod->ndet);
#endif

}

/*--------------------------------------------------------------------------------*/
 void rotate_data(mbTOD *tod, char trans, actData *rotmat)
   /*
    *  Rotate data by rotmat.  If rotmat is null, use the matrix in tod->rotmat.  if that is empty, try tod->corrs
    *  trans tells you if rotmat should be transposed or not.
    */
 {
   assert(tod);
   assert(tod->have_data);
   if (rotmat==NULL) {
     if (tod->rotmat) 
       rotmat=tod->rotmat[0];
     else {  
       assert(tod->corrs);  //need to have something.
       rotmat=tod->corrs[0];
     }
   }
   
   actData **tmp=matrix(tod->ndet,tod->ndata);
   act_gemm('n',trans,tod->ndata,tod->ndet,tod->ndet,1.0,tod->data[0],tod->ndata,rotmat,tod->ndet,0.0,tmp[0],tod->ndata);
   
   memcpy(tod->data[0],tmp[0],tod->ndata*tod->ndet*sizeof(actData));
   
   free(tmp[0]);
   free(tmp);
 }


/*--------------------------------------------------------------------------------*/
 void multiply_all_data(mbTOD *tod,actData val)
 {
   assert(tod->have_data);
#pragma omp parallel for shared(tod,val) default(none)
   for (int i=0;i<tod->ndet;i++)
     for (int j=0;j<tod->ndata;j++)
       tod->data[i][j]*=val;


 }
/*--------------------------------------------------------------------------------*/
void invert_pol_precon(MAP *map)
{
  int poltag=get_map_poltag(map);
#pragma omp parallel shared(map,poltag) default(none)
  {
    actData *mm=map->map;
    switch(poltag){
    case POL_IQU_PRECON:  {
      actData **mymat=matrix(3,3);
#pragma omp for 
      for (int i=0;i<map->npix;i++)  {
	int ii=i*6; //6 polarization in this map
	if (mm[ii]>0) {
	  mymat[0][0]=mm[ii];
	  mymat[0][1]=mymat[1][0]=mm[ii+1]; //Q 
	  mymat[0][2]=mymat[2][0]=mm[ii+2]; //U
	  mymat[1][1]=mm[ii+3];             //Q^2
	  mymat[1][2]=mymat[2][1]=mm[ii+4]; //Q*U
	  mymat[2][2]=mm[ii+5];             //U^2
	  invert_posdef_mat(mymat,3);
	  mm[ii]=mymat[0][0];
	  mm[ii+1]=mymat[0][1];
	  mm[ii+2]=mymat[0][2];
	  mm[ii+3]=mymat[1][1];
	  mm[ii+4]=mymat[1][2];
	  mm[ii+5]=mymat[2][2];
	}
      }
      free(mymat[0]);
      free(mymat);
      break;
    }
    default:
      printf("Not set up yet in invert_pol_precon.  you may have some code to write...\n");
      break;
    }
          
  }
}

/*--------------------------------------------------------------------------------*/
void apply_pol_precon(MAP *map, MAP *precon)
{
  int poltag=get_map_poltag(map);
  int precontag=get_map_poltag(precon);
  if (poltag==POL_IQU) 
    if (precontag!=POL_IQU_PRECON) {
      fprintf(stderr,"ERror in apply_pol_precon - preconditioner is not appropriate for the IQU map.\n");
      return;
    }
  if (poltag==POL_QU) 
    if (precontag!=POL_QU_PRECON) {
      fprintf(stderr,"Error in apply_pol_precon - preconditioner is not appropriate for the QU map.\n");
      return;
    }
  
#pragma omp parallel shared(map,precon,poltag) default(none)
  {
    actData *mm=map->map;
    actData *pp=precon->map;
    switch(poltag){
    case POL_IQU:  {
#pragma omp for schedule(static,512)
      for (int i=0;i<map->npix;i++)  {
	int im=i*3; //3 pols in map;
	int ip=i*6; //6 pols in precon;
	//just multiplying a 3x3 matrix in precon by a 3 element vector in map, but precon only stores half the matrix, 
	//so indexing can look a little hairy
	actData tmp1=mm[im]*pp[ip]+mm[im+1]*pp[ip+1]+mm[im+2]*pp[ip+2];
	actData tmp2=mm[im]*pp[ip+1]+mm[im+1]*pp[ip+3]+mm[im+2]*pp[ip+4];
	actData tmp3=mm[im]*pp[ip+2]+mm[im+1]*pp[ip+4]+mm[im+2]*pp[ip+5];
	mm[im]=tmp1;
	mm[im+1]=tmp2;
	mm[im+2]=tmp3;
	
      }
      
    }
      break;
    default:
      printf("Don't have polarization written up in apply_pol_precon.  You probably have some coding to do.\n");
      break;
    }
  }
}
