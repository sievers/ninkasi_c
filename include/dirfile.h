
#pragma once

#include "mbTOD.h"


void *dirfile_read_channel_direct(char typechar, const char *filename, const char *channelname, int *nsamples_out);


void read_dirfile_tod_data (mbTOD *tod);
actData **read_dirfile_tod_data_from_rowcol_list (mbTOD *tod, int *row, int *col, int ndet, actData **data);

mbTOD *
read_dirfile_tod_header( const char *filename );


int
dirfile_create( const char *dirfilename );

int
dirfile_write_raw_channel( const char *dirfilename,
        const char *channelname, size_t samples_per_frame,
        size_t nframes, char typechar, const void *data );

mbTOD *
read_dirfile_tod( const char *filename );

mbTOD *
read_dirfile_tod_header_decimate( const char *filename , int decimate );
