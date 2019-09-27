#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([cmc],[[latest]],[[s-umbreit@northwestern.edu]])

# We put all auxiliary scripts into common
AC_CONFIG_AUX_DIR([common])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])

AC_CONFIG_SRCDIR([cmc_binbin.c])
AC_CONFIG_HEADER([config.h])

AC_PREFIX_DEFAULT([$HOME])

AC_ARG_ENABLE(static-linking, [AC_HELP_STRING([--enable-static-linking], 
              [Link all external libraries statically])],
	[])

AC_ARG_ENABLE(experimental, [AC_HELP_STRING([--enable-experimental], 
              [Enable a more accurate (but slower) orbit determination.])],
              [], [enable_experimental=no])

AC_ARG_ENABLE([cuda], [AC_HELP_STRING([--enable-cuda],
              [Enable GPU acceleration for nVidia cards using CUDA C])],
              [], [enable_cuda=no])
AM_CONDITIONAL([USE_CUDA], [test "x$enable_cuda" == "xyes"])

AS_IF([test "x$enable_cuda" == "xyes" && test "x$enable_static_linking" == "xyes"], [
      AC_MSG_WARN([CUDA library cannot be linked statically, revert to dynamic linking])
      LIBTOOL_LDFLAGS=""
])

AC_ARG_ENABLE([utils],
              [AC_HELP_STRING([--enable-utils],
                              [Build all the utilities that generate the initial conditions])],
              [], [enable_utils=yes])

AM_CONDITIONAL([BUILD_UTILS], test "x$enable_utils" == "xyes")

AC_ARG_ENABLE([mpi], 
              [AC_HELP_STRING([--enable-mpi], [Enable MPI parallelization.])],
              [], [enable_mpi=yes])

AM_CONDITIONAL([USE_MPI], [test "x$enable_mpi" == "xyes"])

AS_IF([test "x$enable_mpi" == "xyes" && test "x$enable_static_linking" == "xyes"], [
      AC_MSG_WARN([MPI libraries cannot be linked statically, revert to dynamic linking])
      LIBTOOL_LDFLAGS=""
])

AC_ARG_ENABLE([tausworth],
              [AC_HELP_STRING([--enable-tausworth],
                              [Use the Combined Tausworth random number generator used
                               in the CMC part also for the kick routine in BSE.
                               This option is automatically set when OpenMP is enabled.])],
              [], [AS_IF([test $enable_openmp == yes],
                         [enable_tausworth=yes],[enable_tausworth=no])]
              )

AC_ARG_ENABLE([ran3],
              [AC_HELP_STRING([--enable-ran3],
                              [Use Knuth''s subtractive random number generator, ran3,
                               from Numerical Recipes in the BSE kick routine. Cannot be
                               used with OpenMP.])],
              [], [AS_IF([test "$enable_tausworth" == "no"],
                         [enable_ran3=yes],[])]
             )

AS_IF([test "$enable_tausworth" == "yes" && test "$enable_ran3" == "yes"],
      [AC_MSG_ERROR([--enable-tausworth and --enable-ran3 cannot be both set.])], [])

AC_PROG_CC([gcc pathcc icc cc cl.exe /usr/ucb/cc])
# Set the maximum optimization
AX_CC_MAXOPT

AS_IF([test "x$enable_experimental" == "xyes"], [CFLAGS="$CFLAGS -DEXPERIMENTAL"], [])

AS_IF([ test "x$enable_mpi" == "xyes" ], [
    AX_MPI([])
    NU_MPI_FLAGS([])
    CFLAGS="$CFLAGS -DUSE_MPI $MPI_CFLAGS"
    LIBS="$LIBS $MPI_CLINK"], 
[])

AS_IF([test "x$FFLAGS" == x], [FFLAGS='-g -O2'], [])
# Always add debug symbols
CFLAGS="$CFLAGS -g"

# Force the C compiler to issue all warnings
AX_CFLAGS_WARN_ALL([CFLAGS_WARN])
CFLAGS="$CFLAGS $CFLAGS_WARN"

AM_CONDITIONAL([USE_TAUS], test "$enable_tausworth" == "yes")

# Try to find out the name of the host.
AC_ARG_VAR([HOSTNAME], [Full name of the host. Sets default values for certain variables.])
AX_WITH_PROG([HOSTNAME_PROG], [hostname], [])
AS_IF([test "x$HOSTNAME" == x ], [AS_IF([test "x$HOSTNAME_PROG"== x], [],
      [HOSTNAME=`$HOSTNAME_PROG`])])

# Set make variables CMCVERSION and CMCDATE.
GET_REPO_VERSION([VER])
AS_IF([test "x$VER" == "xexported"], [CMCVERSION=`cat VERSION`], 
      [CMCVERSION="$VER"; echo "$VER">VERSION])
# When everything fails, set CMCVERSION to "unknown".
AS_IF([test "x$CMCVERSION" == "x"], [CMCVERSION="unknown"])

AC_CHECK_PROG([HAS_DATE], [date], [yes], [no])
AS_IF([test "x$HAS_DATE" == "xyes"], [CMCDATE=[`date | sed -e 's|[[:space:]]|_|g'`]], [CMCDATE="unknown"])

AC_SUBST([CMCVERSION])
AC_SUBST([CMCDATE])

AS_IF([test "x$HOSTNAME" == "xfugu.phys.northwestern.edu"],
      [AC_PROG_F77([ifort pathf90 gfortran g77])],
      [ AS_IF([test "$enable_tausworth" == "yes"],
              [AC_PROG_F77([ifort g77 gfortran pathf90])] ,[])
        AS_IF([test "$enable_ran3" == "yes"],
              [AC_PROG_F77([ifort g77 gfortran pathf90])], [])
        AS_IF([test "$enable_openmp" == "yes"], [AC_PROG_FC([]) F77=${FC}])
      ]
     )
# see if the f77 compiler can handle integer KINDs
# this is only needed for --enable-tausworth
AS_IF([test "$enable_tausworth" == "yes"],
      [ NU_FORTRAN_INT64_KIND
        AS_IF([test "x$KINDINT64" == "x-1"],
              [AS_IF([test "$enable_openmp" == "no"],
                     [AC_PROG_FC([]) F77=${FC}
                      NU_FORTRAN_INT64_KIND], [])
               AS_IF([test "x$KINDINT64" == "x-1"],
                     [AC_MSG_ERROR([Unable to find a Fortran compiler with integer
                                    KIND value corresponding to 64-bit integers])
                     ], [])
              ], [])
      ])
AS_IF([ test "x$enable_tausworth" == "xyes"],
      [AC_DEFINE([USE_TAUS], 1, [Define whether we use the Combined Tausworth random number generator])
      ],[])

# get the right fortran flags
AC_F77_LIBRARY_LDFLAGS
FILTER_VAR([FLIBS], [-lgcc* -lcrt*.o -LIST:* -LNO:* -lSystem])

AS_IF([ test x$enable_mpi == "xyes"], [
    AC_LANG_PUSH([Fortran 77])
    AX_MPI
    NU_MPI_FLAGS([])
    AC_LANG_POP([Fortran 77])

    FFLAGS="$FFLAGS $MPI_FFLAGS"
    LIBS="$LIBS $MPI_FLINK"
    REMOVE_DUPS([LIBS])
    AC_MSG_NOTICE(["The libs are now $LIBS"]) ], 
[])

#Hack to get around problem with path scale fortran and cfitsio both defining ffseek.
USE_PATHSCALE=${F77##path*}
set ${USE_PATHSCALE:=yes}
AS_IF([test "x$USE_PATHSCALE" == "xyes"], [LIBTOOL_LDFLAGS+=-Wl,-zmuldefs])

# We only use libtool to "pass" the -static flag to the C/FORTRAN compiler in a portable way.
# We do not intend to build shared libraries.
#AC_DISABLE_SHARED
AC_PROG_LIBTOOL

# Checks for programs.
AC_PROG_CXX
AS_IF([test "x$enable_cuda" == "xyes"], [], [CXX=$CC])

AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AS_IF([test "x$enable_cuda" == "xyes"], [
       AX_WITH_PROG([NVCC], [nvcc], [], [], [NVCC_PATH])
       CUDA_CONFIG_VAR([LIBRARIES])
       CUDA_CONFIG_VAR([INCLUDES])
       AC_SUBST([CUDA_CONFIG_LIBRARIES])
       AC_SUBST([CUDA_CONFIG_INCLUDES])
])

# Check if we are on a specific host
AS_IF([test "x$HOSTNAME" == "xfugu.phys.northwestern.edu"], SET_FUGU_DIRS)
AS_IF([test "x$HOSTNAME" == "xfugu.phys.northwestern.edu"], SET_FUGU_FLAGS)

AS_IF([test "x$enable_static_linking" == "xyes"], [LIBTOOL_LDFLAGS+=-all-static], [])
AC_SUBST([LIBTOOL_LDFLAGS])

# Checks for libraries.
AC_CHECK_LIB([m], [sinl], [],[AC_MSG_ERROR([cannot link with libm])], [])
AC_CHECK_LIB([dl], [dlopen])
AC_CHECK_LIB([z], [gzopen], [], [AC_MSG_ERROR([cannot link with zlib])])
AC_SEARCH_LIBS([crypt], [crypt])
AC_SEARCH_LIBS([getopt_long], [gnugetopt], [],
               [AC_MSG_ERROR([cannot link with libgnugetopt], [])], [])
AC_CHECK_LIB([cfitsio], [ffopen], [],[cfitsio_not_found=yes])
AS_IF([test "x$cfitsio_not_found" == "xyes"], [
	ADD_FLAGS_PKG_CONF([cfitsio])
	AS_UNSET([ac_cv_lib_cfitsio_ffopen])
	AC_CHECK_LIB([cfitsio], [ffopen], [], 
                     [AC_MSG_ERROR(["Unable to link against cfitsio."])])
])
	
AC_CHECK_LIB([gslcblas], [cblas_drot], [], [gsl_cblas_not_found=yes])
AS_IF([test "x$gsl_cblas_not_found" == "xyes"], [
	ADD_FLAGS_PKG_CONF([gsl])
	AS_UNSET([ac_cv_lib_gslcblas_cblas_drot])
	AC_CHECK_LIB([gslcblas], [cblas_drot], [],
                     [AC_MSG_ERROR(["Unable to link against gslcblas."])])
])

AC_CHECK_LIB([gsl], [gsl_root_fsolver_alloc], [], [gsl_not_found=yes])
AS_IF([test "x$gsl_not_found" == "xyes"], [
	ADD_FLAGS_PKG_CONF([gsl])
	AS_UNSET([ac_cv_lib_gsl_gsl_root_fsolver_alloc])
	AC_CHECK_LIB([gsl], [gsl_root_fsolver_alloc], [],
                     [AC_MSG_ERROR(["Unable to link against gsl."])])
])

ACX_PTHREAD([])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stddef.h stdlib.h string.h unistd.h])

AC_CHECK_HEADERS([fitsio.h], [], [
	ADD_FLAGS_PKG_CONF([cfitsio])
	AS_UNSET([ac_cv_header_fitsio_h])
	AC_CHECK_HEADERS([fitsio.h], [], 
                         [AC_MSG_ERROR([fitsio.h not available/usable.])])
])

AC_CHECK_HEADERS([gsl/gsl_odeiv.h], [], [
	ADD_FLAGS_PKG_CONF([gsl])
	AS_UNSET([ac_cv_header_gsl_gsl_odeiv_h])
	AC_CHECK_HEADERS([gsl/gsl_odeiv.h], [], 
                         [AC_MSG_ERROR([GSL headers not available/usable.])])
])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow sqrt strchr strcspn strtol])

AC_CONFIG_FILES([bse_wrap/bse/Makefile bse_wrap/bse/COSMIC/cosmic/src/int64.h bse_wrap/Makefile
                 fewbody-0.24/Makefile libs/Makefile Makefile])
AC_OUTPUT

