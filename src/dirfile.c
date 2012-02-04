
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
#include <stdbool.h>
#include <unistd.h>

#include "dirfile.h"

#include "getdata.h"

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

static size_t
bytes_per_sample( char typechar )
{
    switch ( typechar )
    {
        case 'c':
            return 1;

        case 's':
        case 'u':
            return 2;

        case 'S':
        case 'U':
        case 'i':
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

    int nframes = GetNFrames( F, &status, channelname );
    if ( status != GD_E_OK )
    {
        dirfile_print_errstatus(status);
        return NULL;
    }
    assert( nframes > 0 );

    int samples_per_frame = GetSamplesPerFrame( F, channelname, &status );
    if ( status != GD_E_OK )
    {
        dirfile_print_errstatus(status);
        return NULL;
    }
    assert( samples_per_frame > 0 );

    int nsamples = nframes * samples_per_frame;
    size_t nbytes = nsamples * bytes_per_sample(typechar);

    void *data = malloc( nbytes );

    //printf( "%d %d\n", nframes, samples_per_frame );

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

int16_t *
dirfile_read_int16_channel( const struct FormatType *F,
        const char *channelname, int *nsamples )
{
    return (int16_t *) dirfile_read_channel( 's', F, channelname, nsamples );
}

uint16_t *
dirfile_read_uint16_channel( const struct FormatType *F,
        const char *channelname, int *nsamples )
{
    return (uint16_t *) dirfile_read_channel( 'u', F, channelname, nsamples );
}

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

