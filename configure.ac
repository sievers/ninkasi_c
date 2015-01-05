#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([ninkasi], [0.0.1], [sievers@cita.utoronto.ca,nolta@cita.utoronto.ca])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/ninkasi.c])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_ENABLE([actdata-double],
              AC_HELP_STRING([--enable-actdata-double],
                             [defines actData as double instead of float]),
              [],
              [enable_actdata_double=no])
AS_IF([test "x$enable_actdata_double" = xyes],
      AC_DEFINE([ACTDATA_DOUBLE],[1],[Define if you want actData to be double.]) )

# Checks for programs.
AC_USE_SYSTEM_EXTENSIONS([_GNU_SOURCE])
AC_PROG_CC
ACX_MPI([CC="$MPICC" LIBS="$MPILIBS $LIBS"
         AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])])
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_OFF_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_HEADER_TIME

AX_OPENMP([CFLAGS="$CFLAGS $OPENMP_CFLAGS"])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor gettimeofday memset pow sqrt strdup])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
