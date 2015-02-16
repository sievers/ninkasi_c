
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "actpol/dirfile.h"
#include "actpol/getdata.h"

#define dirfile_print_errstatus(STATUS) {\
    if ( STATUS != GD_E_OK ) \
        fprintf(stderr, "line %d: *** dirfile error code: %d\n", __LINE__, status); \
}

ACTpolDirfile *
ACTpolDirfile_open(const char *ifilename)
{
    // remove .zip extension (if present)
    size_t len = strlen(ifilename);
    char *filename = strdup(ifilename);
    if (strcmp(filename + (len-4), ".zip") == 0) {
        filename[len-4] = '\0';
        printf("filename = %s\n", filename);
    }

    int status;
    struct FormatType *format = GetFormat(filename, NULL, &status);
    free(filename);
    if (status != 0) {
        return NULL;
    }

    assert(format);
    ACTpolDirfile *dirfile = malloc(sizeof(ACTpolDirfile));
    dirfile->format = format;
    return dirfile;
}

void
ACTpolDirfile_close(ACTpolDirfile *dirfile)
{
    GetDataClose(dirfile->format);
    free(dirfile->format);
    free(dirfile);
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
ACTpolDirfile_has_channel(const ACTpolDirfile *dirfile, const char *channel)
{
    int status;
    int samples_per_frame = GetSamplesPerFrame(dirfile->format, channel, &status);
    if ( status != GD_E_OK || samples_per_frame <= 0 )
        return false;
    return true;
}

//static void *
void *
ACTpolDirfile_read_channel( char typechar, const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples_out )
{
    int status = 0;
    const struct FormatType *F = dirfile->format;

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

    // last frame may be partial
    int nsamples = (nframes + 1) * samples_per_frame - 1;
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
ACTpolDirfile_read_int16_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples )
{
    return (int16_t *) ACTpolDirfile_read_channel( 's', dirfile, channelname, nsamples );
}

uint16_t *
ACTpolDirfile_read_uint16_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples )
{
    return (uint16_t *) ACTpolDirfile_read_channel( 'u', dirfile, channelname, nsamples );
}

int32_t *
ACTpolDirfile_read_int32_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples )
{
    return (int32_t *) ACTpolDirfile_read_channel( 'S', dirfile, channelname, nsamples );
}

uint32_t *
ACTpolDirfile_read_uint32_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples )
{
    return (uint32_t *) ACTpolDirfile_read_channel( 'U', dirfile, channelname, nsamples );
}

float *
ACTpolDirfile_read_float_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples )
{
    return (float *) ACTpolDirfile_read_channel( 'f', dirfile, channelname, nsamples );
}

double *
ACTpolDirfile_read_double_channel(const ACTpolDirfile *dirfile,
        const char *channelname, int *nsamples )
{
    return (double *) ACTpolDirfile_read_channel( 'd', dirfile, channelname, nsamples );
}

uint32_t
ACTpolDirfile_read_uint32_sample(const ACTpolDirfile *dirfile,
        const char *channelname, int index )
{
    uint32_t sample;
    int status;
    int n = GetData(dirfile->format, channelname, 0, index, 0, 1, 'U', &sample, &status);

    if ( status != GD_E_OK || n != 1 )
    {
        dirfile_print_errstatus( status );
        return 0;
    }

    return sample;
}

