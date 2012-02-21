
#pragma once

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline double mysecs2days( double s ) { return s/86400.; }
static inline double mydeg2rad( double deg ) { return deg*M_PI/180.; }
static inline double myrad2deg( double rad ) { return rad*180./M_PI; }
static inline double myarcsec2rad( double sec ) { return mydeg2rad(sec/3600.); }


typedef struct
{
    double latitude;            /* radians */
    double east_longitude;      /* radians */
    double elevation_m;
    double temperature_K;
    double pressure_mb;
    double relative_humidity;   /* percentage */
}
Site;

int
observed_altaz_to_mean_radec( const Site *site, double freq_GHz,
        int n, const double ctime[], const float alt[], const float az[],
        float ra[], float dec[] );


int
dobserved_altaz_to_mean_radec( const Site *site, double freq_ghz,
        int n, const double ctime[],
        const double alt[], const double az[],
        double ra[], double dec[] );


void
ACTSite( Site *p );
