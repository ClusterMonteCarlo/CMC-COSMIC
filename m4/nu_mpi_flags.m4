# ===============================================================
#       MPI Link and Compiler Flags
# ===============================================================
#
# SYNOPSIS
#
#   NU_MPI_FLAGS
#
# DESCRIPTION
#
#   Tries to figure out the compile and link flags needed to
#   include MPI support. Requires that AX_MPI was run first.
#
#   On success, sets MPI_CFLAGS, MPI_LIBS with the corresponding
#   compile and link flags.
#
# LICENSE
#
#   Copyright (c) 2012 Stefan Umbreit <s-umbreit@northwestern.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.
#
AC_DEFUN([_NU_DROP_FIRST_ARG], [
    gawk '{for (i=2; i<=NF; i++) printf \"%s \", \[$]i}'
])

AC_DEFUN([_NU_MPI_EXTRACT_FLAGS], [
    NU_MPI_FLAGS_IMPL=$1
    MPICOMP=$$2
    _mpiwrap_unknown=no
    AS_CASE([$NU_MPI_FLAGS_IMPL], 
        [OpenMPI], [_mpiwrap_compile='-showme:compile'; 
                    _mpiwrap_link='-showme:link'],
        [Lam-MPI], [_mpiwrap_compile='-showme:compile';
                    _mpiwrap_link='-showme:link'],
        [MPICH2],  [_mpiwrap_compile="-compile-info| _NU_DROP_FIRST_ARG";
                    _mpiwrap_link="-link-info | _NU_DROP_FIRST_ARG"],
        [IMPI],  [_mpiwrap_compile="-compile-info| _NU_DROP_FIRST_ARG";
                    _mpiwrap_link="-link-info | _NU_DROP_FIRST_ARG"],
        [_mpiwrap_unknown=yes 
         AC_MSG_WARN([Unknown MPI implementation.])
    ])

    _nu_mpi_flags_comp_flags=$($SHELL -c "[$]MPICOMP $_mpiwrap_compile")
    #_nu_mpi_flags_link_flags=$($SHELL -c "[$]MPICOMP $_mpiwrap_link")

])

AC_DEFUN([_NU_MPI_EXTRACT_IMPL], [
    AC_REQUIRE([AX_MPI])
    # Here we only check for the C implementation and assume it to be the 
    # same for C++ and Fortran too.
    AC_LANG_PUSH([C])
    nu_mpi_flags_save="$CC"
    CC="$MPICC"
    AC_MSG_CHECKING([which MPI implementation we are using])
    AC_RUN_IFELSE([AC_LANG_SOURCE([
       #include <stdio.h>
       #include <mpi.h>
       #ifdef OPEN_MPI
       #define RESULT "OpenMPI"
       #elif defined MPICH2
       #define RESULT "MPICH2"
       #elif defined I_MPI_VERSION
       #define RESULT "IMPI"
       #elif defined LAM_MPI
       #define RESULT "Lam-MPI"
       #else
       #define RESULT "Unknown"
       #endif
       int main(void) {
         FILE*conftest;

         conftest=fopen("conftestval", "w");
         fprintf(conftest, "%s\n", RESULT);
         fclose(conftest);
         return(0);
       }
    ])], [NU_MPI_FLAGS_IMPL=`cat conftestval`], [])
    AC_MSG_RESULT([$NU_MPI_FLAGS_IMPL])
    AC_LANG_POP([])
    CC="$nu_mpi_flags_save"
    AC_SUBST([NU_MPI_FLAGS_IMPL])
])

AC_DEFUN([NU_MPI_FLAGS], [
    # The strategy is to include <mpi.h> and check for implementation 
    # specific #defines. This information is then used to call MPICC 
    # with the appropriate flag to extract the required link/compile flags.

    _NU_MPI_EXTRACT_IMPL([])

    AC_LANG_CASE([C], [
       AC_REQUIRE([AX_MPI])
       _NU_MPI_EXTRACT_FLAGS([$NU_MPI_FLAGS_IMPL], [MPICC])
       MPI_CFLAGS=$_nu_mpi_flags_comp_flags
       MPI_CLINK=$_nu_mpi_flags_link_flags
       AC_SUBST([MPI_CFLAGS])
       AC_SUBST([MPI_CLINK])
    ],
    [C++], [
       AC_REQUIRE([AX_MPI])
       _NU_MPI_EXTRACT_FLAGS([$NU_MPI_FLAGS_IMPL], [MPICXX])
       MPI_CXXFLAGS=$_nu_mpi_flags_comp_flags
       MPI_CXXLINK=$_nu_mpi_flags_link_flags
       AC_SUBST([MPI_CXXFLAGS])
       AC_SUBST([MPI_CXXLINK])
    ],
    [Fortran 77], [
       AC_REQUIRE([AX_MPI])
       _NU_MPI_EXTRACT_FLAGS([$NU_MPI_FLAGS_IMPL], [MPIF77])
       MPI_FFLAGS=$_nu_mpi_flags_comp_flags
       MPI_FLINK=$_nu_mpi_flags_link_flags
       AC_SUBST([MPI_FFLAGS])
       AC_SUBST([MPI_FLINK])
    ],
    [Fortran], [
       AC_REQUIRE([AX_MPI])
       _NU_MPI_EXTRACT_FLAGS([$NU_MPI_FLAGS_IMPL], [MPIF90])
       MPI_FCFLAGS=$_nu_mpi_flags_comp_flags
       MPI_FCLINK=$_nu_mpi_flags_link_flags
       AC_SUBST([MPI_FCFLAGS])
       AC_SUBST([MPI_FCLINK])
    ])

    #Check if the compiler can really compile and link a simple MPI program
    _test_source_c="int *argc; char ***argv; MPI_Init(argc, argv)"
    _test_source_c_prolog="#include <mpi.h>"
    _test_source_f="        call MPI_Init()"

    nu_mpi_flags_cflags=$CFLAGS
    nu_mpi_flags_cxxflags=$CXXFLAGS
    nu_mpi_flags_fflags=$FFLAGS
    nu_mpi_flags_fcflags=$FCFLAGS
    nu_mpi_flags_libs=$LIBS

    AC_LANG_CASE([C], [
      CFLAGS="$CFLAGS $MPI_CFLAGS"
      LIBS="$LIBS $MPI_CLINK"
      AC_MSG_NOTICE([The flags are: $MPI_CFLAGS])
      AC_MSG_NOTICE([The flags are: $MPI_CLINK])
      AC_MSG_CHECKING([if a simple C MPI program can be compiled and linked])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([$_test_source_c_prolog], [$_test_source_c])], 
        [AC_MSG_RESULT([yes])], 
        [AC_MSG_ERROR([Unable to compile and link simple C MPI program.])])
      CFLAGS=$nu_mpi_flags_cflags
      AC_SUBST([MPI_CFLAGS])
      AC_SUBST([MPI_CLINK])
    ],
    [C++], [
      CXXFLAGS="$CXXFLAGS $MPI_CXXFLAGS"
      LIBS="$LIBS $MPI_CXXLINK"
      AC_MSG_CHECKING([if a simple C++ MPI program can be compiled and linked])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([$_test_source_c_prolog], [$_test_source_c])],
        [AC_MSG_RESULT([yes])], 
        [AC_MSG_ERROR([Unable to compile and link simple C++ MPI program.])])
      CXXFLAGS=$nu_mpi_flags_cxxflags
      AC_SUBST([MPI_CXXFLAGS])
      AC_SUBST([MPI_CXXLINK])
    ],
    [Fortran 77], [
      FFLAGS="$FFLAGS $MPI_FFLAGS"
      LIBS="$LIBS $MPI_FLINK"
      AC_MSG_CHECKING([if a simple F77 MPI program can be compiled and linked])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [$_test_source_f])],
        [AC_MSG_RESULT([yes])], 
        [AC_MSG_ERROR([Unable to compile and link simple F77 MPI program.])])
      FFLAGS=$nu_mpi_flags_fflags
      AC_SUBST([MPI_FFLAGS])
      AC_SUBST([MPI_FLINK])
    ],
    [Fortran], [
      FCFLAGS="$FCFLAGS $MPI_FCFLAGS"
      LIBS="$LIBS $MPI_FCLINK"
      AC_MSG_CHECKING([if a simple F90 MPI program can be compiled and linked])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [$test_source_f])],
        [AC_MSG_RESULT([yes])], 
        [AC_MSG_ERROR([Unable to compile and link simple F90 MPI program.])])
      AC_SUBST([MPI_FCFLAGS])
      AC_SUBST([MPI_FCLINK])
      FCFLAGS=$nu_mpi_flags_fcflags
    ])
    LIBS=$nu_mpi_flags_libs
])
