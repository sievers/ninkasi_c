
#pragma once

#include "mbTOD.h"


void *dirfile_read_channel_direct(char typechar, const char *filename, const char *channelname, int *nsamples_out);


void read_dirfile_tod_data (mbTOD *tod);
actData **read_dirfile_tod_data_from_rowcol_list (mbTOD *tod, int *row, int *col, int ndet, actData **data);

mbTOD *
read_dirfile_tod_header( const char *filename );


mbTOD *
read_dirfile_tod( const char *filename );

mbTOD *
read_dirfile_tod_header_decimate( const char *filename , int decimate );

mbTOD *
read_dirfile_tod_header_abs( const char *filename );
