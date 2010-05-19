#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.63)
AC_INIT([cmc],[[latest]],[[s-umbreit@northwestern.edu]])

# We put all auxiliary scripts into common
AC_CONFIG_AUX_DIR([common])

AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([cmc_binbin.c])
AC_CONFIG_HEADER([config.h])

AS_IF([test "x$FFLAGS" == x], [FFLAGS='-g'], [])
AS_IF([test "x$CFLAGS" == x], [CFLAGS='-g'], [])

# try to find out the name of the host 
AC_ARG_VAR([HOSTNAME], [Full name of the host. Sets default values for certain variables.])
AX_WITH_PROG([HOSTNAME_PROG], [hostname], [])
AS_IF([test "x$HOSTNAME" == x ], [AS_IF([test "x$HOSTNAME_PROG"== x], [],[HOSTNAME=`$HOSTNAME_PROG`])])

# We only use libtool to "pass" the -static flag to the C/FORTRAN compiler in a portable way.
# We do not intend to build shared libraries.
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

# get the right fortran flags
AC_PROG_F77([gfortran pathf90 ifort g77])
AC_F77_LIBRARY_LDFLAGS

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# First, check if we are on a specific host
AS_IF([test "x$HOSTNAME" == "xfugu.phys.northwestern.edu"], SET_FUGU_DIRS)
AS_IF([test "x$HOSTNAME" == "xfugu.phys.northwestern.edu"], SET_FUGU_FLAGS)

# Checks for libraries.
AC_CHECK_LIB([m], [sinl], [],[AC_MSG_ERROR([cannot link with libm], [])], [])
AC_CHECK_LIB([dl], [dlopen])
AC_CHECK_LIB([z], [gzopen], [], [AC_MSG_ERROR([cannot link with zlib])])
AC_SEARCH_LIBS([crypt], [crypt])
AC_SEARCH_LIBS([getopt_long], [gnugetopt], [],[AC_MSG_ERROR([cannot link with libgnugetopt], [])], [])
AC_CHECK_LIB([cfitsio], [ffopen], [],[AC_MSG_ERROR([cannot link with cfitsio], [])], [])
AC_CHECK_LIB([gslcblas], [cblas_drot])
AC_CHECK_LIB([gsl], [gsl_root_fsolver_alloc])
ACX_PTHREAD([])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stddef.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow sqrt strchr strcspn strtol])

AC_CONFIG_FILES([bse_wrap/bse/Makefile bse_wrap/Makefile fewbody-0.24/Makefile libs/Makefile Makefile])
AC_OUTPUT

                 
