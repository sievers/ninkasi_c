
#pragma once

#include "mbTOD.h"

void read_dirfile_tod_data (mbTOD *tod);

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

