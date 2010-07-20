/// \file slimlib.h
/// Function prototypes for the slimlib C library.
/// Include this in order to use the C wrappers for slim
/// compress/expand operations and the SLIMFILE opaque type.

#ifndef _SLIMLIB_H
#define _SLIMLIB_H

#include <stdio.h> // for size_t

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

typedef struct slimlib_file_t  SLIMFILE;

extern SLIMFILE *slimopen(const char *filename,
			  const char *modes);

extern int slimclose(SLIMFILE *sf);

extern size_t slimread(void *ptr, size_t size, size_t nmemb, SLIMFILE *sf);

extern long slimtell(SLIMFILE *sf);

extern void slimrewind(SLIMFILE *sf);

extern int slimseek(SLIMFILE *sf, long offset, int whence);

extern long slimrawsize(const char *filename);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif

#endif // #ifndef _SLIMLIB_H
