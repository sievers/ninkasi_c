
#include "config.h"

#include <assert.h>
#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>


#include "dirfile.h"
#include "getdata.h"

#define ACT_ARRAY_MAX_ROWS 33
#define ACT_ARRAY_MAX_COLS 32
#define ACT_ARRAY_MAX_DETECTORS (ACT_ARRAY_MAX_ROWS*ACT_ARRAY_MAX_COLS)

#define MAX(A,B) ((A > B) ? (A) : (B))

// ----------------------------------------------------------------------------

/*
static void
dirfile_print_errstatus( int status )
{
    if ( status != GD_E_OK )
        fprintf(stderr, "*** dirfile error code: %d\n", status);
}
*/

#define dirfile_print_errstatus(STATUS) {\
    if ( STATUS != GD_E_OK ) \
        fprintf(stderr, "line %d: *** dirfile error code: %d\n", __LINE__, status); \
}

size_t
bytes_per_sample( char typechar )
{
    switch ( typechar )
    {
        case 'S':
        case 'U':
        case 'f':
            return 4;

        case 'd':
            return 8;
    }

    assert( 1 == 0 );
    return 0;
}

bool
dirfile_has_channel( const struct FormatType *F, const char *channel )
{
    int status;
    int samples_per_frame = GetSamplesPerFrame( F, channel, &status );
    if ( status != GD_E_OK || samples_per_frame <= 0 )
        return false;
    return true;
}


void *
dirfile_read_channel( char typechar, const struct FormatType *F,
        const char *channelname, int *nsamples_out )
{
    int status = 0;

    /*
    int nframes = GetNFrames( F, &status, channelname );
    if ( status != GD_E_OK )
    {
        dirfile_print_errstatus(status);
        return NULL;
    }
    assert( nframes > 0 );
    */

    int nsamples = GetNSamples( F, &status, channelname );
    if ( status != GD_E_OK )
    {
        dirfile_print_errstatus(status);
        return NULL;
    }
    assert( nsamples > 0 );

    int samples_per_frame = GetSamplesPerFrame( F, channelname, &status );
    if ( status != GD_E_OK )
    {
        dirfile_print_errstatus(status);
        return NULL;
    }
    assert( samples_per_frame > 0 );

    //int nsamples = nframes * samples_per_frame;
    size_t nbytes = nsamples * bytes_per_sample(typechar);

    void *data = malloc( nbytes );

    *nsamples_out = GetData( F, channelname, 0, 0,
            nsamples / samples_per_frame,
            nsamples % samples_per_frame,
            typechar, data, &status );

    if ( status != GD_E_OK || *nsamples_out <= 0 )
    {
        dirfile_print_errstatus( status );
        free( data );
        return NULL;
    }

    return data;
}
/*--------------------------------------------------------------------------------*/

void *dirfile_read_channel_direct(char typechar, const char *filename, const char *channelname, int *nsamples_out)
{
  int status;
  *nsamples_out=0;
  struct FormatType *format = GetFormat( filename, NULL, &status );
  if (format==NULL) {
    fprintf(stderr,"Warning - problem reading format from file .%s.\n",filename);
    return NULL;
  }
  void *data=dirfile_read_channel(typechar,format,channelname,nsamples_out);
  if (data==NULL) {
    fprintf(stderr,"Error reading channel %s from %s\n",channelname,filename);

  }

  GetDataClose(format);
  free(format);
  return data;

}

/*--------------------------------------------------------------------------------*/


int32_t *
dirfile_read_int32_channel( const struct FormatType *F,
        const char *channelname, int *nsamples )
{
    return (int32_t *) dirfile_read_channel( 'S', F, channelname, nsamples );
}

uint32_t *
dirfile_read_uint32_channel( const struct FormatType *F,
        const char *channelname, int *nsamples )
{
    return (uint32_t *) dirfile_read_channel( 'U', F, channelname, nsamples );
}

float *
dirfile_read_float_channel( const struct FormatType *F,
        const char *channelname, int *nsamples )
{
    return (float *) dirfile_read_channel( 'f', F, channelname, nsamples );
}

double *
dirfile_read_double_channel( const struct FormatType *F,
        const char *channelname, int *nsamples )
{
    return (double *) dirfile_read_channel( 'd', F, channelname, nsamples );
}

// ----------------------------------------------------------------------------

uint32_t
dirfile_read_uint32_sample( const struct FormatType *F,
        const char *channelname, int index )
{
    uint32_t sample;
    int status;
    int n = GetData( F, channelname, 0, index, 0, 1, 'U', &sample, &status );

    if ( status != GD_E_OK || n != 1 )
    {
        dirfile_print_errstatus( status );
        return 0;
    }

    return sample;
}

// ----------------------------------------------------------------------------

void
write_tes_channelname( char *s, int row, int col )
{
    snprintf( s, 14, "tesdatar%2.2dc%2.2d", row, col );
}


// ----------------------------------------------------------------------------


mbTOD *
read_dirfile_tod_header( const char *filename )
{
  //printf("reading .%s.\n",filename);
    assert( filename != NULL );

    int status, n;
    //printf("reading format from %s\n",filename);
    struct FormatType *format = GetFormat( filename, NULL, &status );
    assert( format != NULL );

    mbTOD *tod = (mbTOD *) calloc( 1,sizeof(mbTOD) );
    tod->dirfile=strdup(filename);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // alt/az


#ifdef ACTDATA_DOUBLE
    float *az = dirfile_read_float_channel( format, "Enc_Az_Deg", &n );
    tod->ndata = n;
    tod->az=(actData *)malloc(sizeof(actData)*n);
    for (int i=0;i<tod->ndata;i++)
      tod->az[i]=az[i];
    free(az);
    float *alt = dirfile_read_float_channel( format, "Enc_El_Deg", &n );
    tod->alt=(actData *)malloc(sizeof(actData)*n);
    for (int i=0;i<n;i++)
      tod->alt[i]=alt[i];
    free(alt);


#else
    tod->az = dirfile_read_float_channel( format, "Enc_Az_Deg", &n );
    tod->ndata = n;
    tod->alt = dirfile_read_float_channel( format, "Enc_El_Deg", &n );
#endif
    assert( n == tod->ndata );

    //printf("in dirfile, az/alt[0] are %12.4g %12.4g with %d elements.\n",tod->az[0],tod->alt[0],n);

    // convert deg->rad; also kuka az -> astro az
    const float degToRad = M_PI/180;
    for ( int i = 0; i < tod->ndata; i++ )
    {
        tod->az[i] = (tod->az[i] + 180.) * degToRad;
        tod->alt[i] *= degToRad;
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ctime

    uint32_t s0, us0, s1, us1;

    s0 = dirfile_read_uint32_sample( format, "cpu_s", 0 );
    s1 = dirfile_read_uint32_sample( format, "cpu_s", tod->ndata-1 );
    us0 = dirfile_read_uint32_sample( format, "cpu_us", 0 );
    us1 = dirfile_read_uint32_sample( format, "cpu_us", tod->ndata-1 );

    //printf("s0,s1, us0, us1 are %8d %8d %8d %8d\n",s0,s1,us0,us1);
    

    double ctime0 = s0 + 1e-6*us0;
    double ctime1 = s1 + 1e-6*us1;

    tod->ctime = ctime0;
    tod->deltat = (ctime1 - ctime0)/(tod->ndata-1);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // row/col info

    tod->ndet = 0;
    tod->nrow = 0;
    tod->ncol = 0;
    tod->rows = malloc( ACT_ARRAY_MAX_DETECTORS*sizeof(int) );
    tod->cols = malloc( ACT_ARRAY_MAX_DETECTORS*sizeof(int) );

    char tesfield[] = "tesdatar00c00";

    for ( int r = 0; r < ACT_ARRAY_MAX_ROWS; r++ )
    for ( int c = 0; c < ACT_ARRAY_MAX_COLS; c++ )
    {
        write_tes_channelname( tesfield, r, c );
        if ( dirfile_has_channel(format, tesfield) )
        {
            int idet = tod->ndet++;
            tod->rows[idet] = r;
            tod->cols[idet] = c;
            tod->nrow = MAX( tod->nrow, r+1 );
            tod->ncol = MAX( tod->ncol, c+1 );
        }
    }

    tod->rows = (int *) realloc( tod->rows, tod->ndet*sizeof(int) );
    tod->cols = (int *) realloc( tod->cols, tod->ndet*sizeof(int) );

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // data
    tod->data=NULL;  //should already be set thusly.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GetDataClose(format);
    free(format);

    return tod;
}

// ----------------------------------------------------------------------------

actData *decimate_vector(actData *vec, int *nn) //shrink a vector by a factor of 2, free the old vector, return the new.
{
  int n=*nn;
  int n2=(n+1)/2;  //round up if we are odd
  
  actData *vec2=(actData *)malloc(n2*sizeof(actData));
  for (int i=0;i<n2-1;i++) {
    int ii=2*i;
    vec2[i]=0.25*(vec[ii]+2*vec[ii+1]+vec[ii+2]);
  }
  //Now do ends
  if (2*n2==n)  //if we have an even number of elements
    vec2[n2-1]=0.25*(vec[n-3]+2*vec[n-2]+vec[n-1]);
  else
    vec2[n2-1]=0.5*(vec[n-2]+vec[n-1]);
  
  free(vec);
  *nn=n2;
  return vec2;
  
  
  
}

// ----------------------------------------------------------------------------

actData *decimate_vector_in_place(actData *vec, int *nn) //shrink a vector by a factor of 2 in-place, replace nn by its new length.
{

  int n=*nn;
  int n2=(n+1)/2;  //round up if we are odd

  actData tmp0=0.25*(vec[0]+2*vec[1]*vec[2]);
  actData tmp1=0.25*(vec[2]+2*vec[3]*vec[4]);
  for (int i=2;i<n2-1;i++) {
    int ii=2*i;
    vec[i]=0.25*(vec[ii]+2*vec[ii+1]+vec[ii+2]);
  }
  if (2*n2==n)
    vec[n2-1]=0.25*(vec[n-3]+2*vec[n-2]+vec[n-1]);
  else  
    vec[n2-1]=0.5*(vec[n-2]+vec[n-1]);
  vec[0]=tmp0;
  vec[1]=tmp1;
  *nn=n2;

}



// ----------------------------------------------------------------------------


mbTOD *
read_dirfile_tod_header_decimate( const char *filename , int decimate )
{
  //fprintf(stderr,"inside read_dirfile_tod_header_decimate.\n");
  mbTOD *tod=read_dirfile_tod_header(filename);

  if (decimate==0)
    return tod;
  
  int n=tod->ndata;
  for (int i=0;i<decimate;i++) {
    int n1=n;
    tod->az=decimate_vector(tod->az,&n1);
    int n2=n;
    tod->alt=decimate_vector(tod->alt,&n2);
    assert(n2==n1);
    tod->deltat*=2.0;
    n=n1;
  }
  tod->ndata=n;
  tod->decimate=decimate;

  return tod;
    
}


// ----------------------------------------------------------------------------



actData **read_dirfile_tod_data_from_rowcol_list (mbTOD *tod, int *row, int *col, int ndet, actData **data)
//if data is non-null, assume the space is allocated.
{
  
  int status;
  assert(tod!=NULL);
  //printf("reading format from %s.\n",tod->dirfile);
  struct FormatType *format = GetFormat( tod->dirfile, NULL, &status );
  //printf("got it.\n");
  assert( format != NULL );
  

  
  if (data==NULL) {
    actData *vec=(actData *)malloc(ndet*tod->ndata*sizeof(actData));
    assert(vec!=NULL);
    for (int i=0;i<ndet;i++)
      data[i]=vec+i*tod->ndata;
  }
  
  for (int idet=0;idet<ndet;idet++) {
    assert(row[idet]<33);
    assert(col[idet]<32);
    char tesfield[]="tesdatar00c00";
    sprintf(tesfield,"tesdatar%02dc%02d",row[idet],col[idet]);
    int n;
#ifndef ACTDATA_DOUBLE
    //printf("reading float channel\n");
    actData *chan = dirfile_read_float_channel( format, tesfield, &n );
#else
    //printf("reading double channel\n");
    actData *chan = dirfile_read_double_channel( format, tesfield, &n );
#endif
    //printf("finished channel.\n");
    for (int ii=0;ii<tod->decimate;ii++)
      chan=decimate_vector(chan,&n);
    assert(n==tod->ndata);
    memcpy(data[idet],chan,sizeof(actData)*tod->ndata);
    free(chan);
  }
  return data;
}



// ----------------------------------------------------------------------------

void read_dirfile_tod_data (mbTOD *tod)
{
#if 1
  assert(tod!=NULL);
  if (!tod->have_data)
    tod->data=NULL;
  tod->data=read_dirfile_tod_data_from_rowcol_list(tod,tod->rows,tod->cols,tod->ndet,tod->data);
  tod->have_data=1;
  

  
#else
  int status;
  assert(tod!=NULL);
  struct FormatType *format = GetFormat( tod->dirfile, NULL, &status );
  assert( format != NULL );

  if (tod->data==NULL) {
    assert(tod->ndata>0);  //make sure we know how big we ought to be.
    tod->data = (actData **) malloc( tod->ndet*sizeof(actData *) );
    actData *vec=(actData *)malloc(tod->ndet*tod->ndata*sizeof(actData));
    for (int i=0;i<tod->ndet;i++)
      tod->data[i]=&(vec[i*tod->ndata]);
  }
  //try skipping parallel read in case it bogs disks down.
  //#pragma omp parallel for shared(tod,format) default(none)    
  for ( int idet = 0; idet < tod->ndet; idet++ )  {
    char tesfield[] = "tesdatar00c00";
    int n;
    int row = tod->rows[idet];
    int col = tod->cols[idet];
    write_tes_channelname( tesfield, row, col );
#ifndef ACTDATA_DOUBLE
    actData *chan = dirfile_read_float_channel( format, tesfield, &n );
    for (int i=0;i<tod->decimate;i++)
      chan=decimate_vector(chan,&n);
    assert( n == tod->ndata );
    memcpy(tod->data[idet],chan,sizeof(actData)*tod->ndata);
    
#else
    double *chan =  dirfile_read_double_channel( format, tesfield, &n );
    for (int i=0;i<tod->decimate;i++)
      chan=decimate_vector(chan,&n);
    
    if (n!=tod->ndata) {
      printf("Warning - size mismatch in read_tod_data, %d vs. %d on detector %d %d, and file %s\n",n,tod->ndata,tod->rows[idet],tod->cols[idet],tod->dirfile);
    }
    assert( n == tod->ndata );
    for (int i=0;i<tod->ndata;i++)
      tod->data[idet][i]=chan[i];
    //memcpy(tod->data[idet],chan,sizeof(actData)*tod->ndata);

#endif
    free(chan);
  }
#endif
}

// ----------------------------------------------------------------------------


mbTOD *
read_dirfile_tod( const char *filename )
{
  mbTOD *tod=read_dirfile_tod_header(filename );
  read_dirfile_tod_data(tod);
  return tod;
}


// ----------------------------------------------------------------------------

static char *
enter_dirfile( const char *dirfilename )
{
    struct stat st;
    int status = stat( dirfilename, &st );
    if ( status != 0 || S_ISDIR(st.st_mode) == 0 )
        return NULL;

    char *cwd = getcwd( NULL, 0 );
    if ( cwd == NULL )
        return NULL;

    status = chdir( dirfilename );
    if ( status != 0 )
    {
        free(cwd);
        return NULL;
    }

    return cwd;
}

static void
exit_dirfile( char *oldcwd )
{
    chdir( oldcwd );
    free( oldcwd );
}

int
dirfile_create( const char *dirfilename )
{
    int status = mkdir( dirfilename,
            S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH );
    if ( status != 0 )
        return -1;

    char *oldwd = enter_dirfile( dirfilename );
    if ( oldwd == NULL )
        return -1;

    FILE *fp = fopen( "format", "w" );
    fprintf( fp, "# dirfile format file\n" );
    fclose(fp);

    exit_dirfile( oldwd );
    return 0;
}

int
dirfile_write_raw_channel( const char *dirfilename,
        const char *channelname, size_t samples_per_frame,
        size_t nframes, char typechar, const void *data )
{
    int status = -1;

    assert( data != NULL );
    size_t nbytes = bytes_per_sample(typechar) * samples_per_frame * nframes;
    assert( nbytes > 0 );

    char *oldwd = enter_dirfile( dirfilename );
    if ( oldwd == NULL )
        goto error;

    int fd = open( channelname, O_WRONLY|O_CREAT|O_EXCL, S_IRUSR|S_IRGRP|S_IROTH );
    ssize_t n = write( fd, data, nbytes );
    close( fd );
    if ( n < 0 )
        goto error;

    FILE *fp = fopen( "format", "a" );
    fprintf( fp, "%s RAW %c %d\n", channelname, typechar, samples_per_frame );
    fclose( fp );

    status = 0;
error:
    exit_dirfile( oldwd );
    return status;
}
