
#pragma once


#define ACTPOL_DIRFILE

#ifdef ACTPOL_DIRFILE
#include "actpol/dirfile.h"
#include "actpol/getdata.h"

void *dirfile_read_channel(char typechar, const ACTpolDirfile *dirfile,const char *channelname, int *nsamples_out );
float *dirfile_read_float_channel(const ACTpolDirfile *dirfile, const char *channelname, int *nsamples);
double *dirfile_read_double_channel(const ACTpolDirfile *dirfile, const char *channelname, int *nsamples);
uint32_t *dirfile_read_uint32_channel(const ACTpolDirfile *dirfile, const char *channelname, int *nsamples);
uint32_t dirfile_read_uint32_sample(const ACTpolDirfile *dirfile, const char *channelname, int index);
bool dirfile_has_channel(const ACTpolDirfile *dirfile, const char *channel);
#else
#include "dirfile.h"
#include "getdata.h"
#endif



#include "mbTOD.h"
void *dirfile_read_channel_direct(char typechar, const char *filename, const char *channelname, int *nsamples_out);


void read_dirfile_tod_data (mbTOD *tod);
actData **read_dirfile_tod_data_from_rowcol_list (mbTOD *tod, int *row, int *col, int ndet, actData **data, int *nout);

mbTOD *
read_dirfile_tod_header( const char *filename );


mbTOD *
read_dirfile_tod( const char *filename );

mbTOD *
read_dirfile_tod_header_decimate( const char *filename , int decimate );

mbTOD *
read_dirfile_tod_header_abs( const char *filename );
