#Sets default values to flags and search paths for fugu.phys.northwestern.edu

AC_DEFUN([SET_FUGU_DIRS], [
          LDFLAGS+=" -L/share/apps/cfitsio/lib"
          LDFLAGS+=" -L/share/apps/gsl/lib"
          CPPFLAGS+=" -I/share/apps/gsl/include"
          CPPFLAGS+=" -I/share/apps/cfitsio/include"
          ])

AC_DEFUN([SET_FUGU_FLAGS], [
        AS_IF([ test "x$enable_cuda" == "xno" && test "x$enable_mpi" == "xno"], [
        LIBTOOL_LDFLAGS+=" -all-static"
        ])
        ])

