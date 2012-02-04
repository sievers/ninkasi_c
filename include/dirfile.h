
#pragma once

#include <stdbool.h>
#include <stdint.h>

struct FormatType;

bool
dirfile_has_channel( const struct FormatType *F, const char *channel );

void *
dirfile_read_channel( char typechar, const struct FormatType *F,
        const char *channelname, int *nsamples_out );

int16_t *
dirfile_read_int16_channel( const struct FormatType *F,
        const char *channelname, int *nsamples );

uint16_t *
dirfile_read_uint16_channel( const struct FormatType *F,
        const char *channelname, int *nsamples );

int32_t *
dirfile_read_int32_channel( const struct FormatType *F,
        const char *channelname, int *nsamples );

uint32_t *
dirfile_read_uint32_channel( const struct FormatType *F,
        const char *channelname, int *nsamples );

float *
dirfile_read_float_channel( const struct FormatType *F,
        const char *channelname, int *nsamples );

double *
dirfile_read_double_channel( const struct FormatType *F,
        const char *channelname, int *nsamples );

uint32_t
dirfile_read_uint32_sample( const struct FormatType *F,
        const char *channelname, int index );

