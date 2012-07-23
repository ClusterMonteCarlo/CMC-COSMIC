# ===============================================================
#       Fortran 64-bit integer kind
# ===============================================================
#
# SYNOPSIS
#
#   NU_FORTRAN_INT64_KIND
#
# DESCRIPTION
#
#   Determines the integer KIND value corresponding to a 64-bit integer.
#   This function is only necessary to configure code that is supposed
#   to be compiled by both a Fortran77 compiler, that has support
#   KIND values but not the intrinsic selected_int_kind (e.g. g77),
#   and any Fortran90(95,03) compiler, that naturally supports both.
#
#   On success, the variable KINDINT64 contains the numeric KIND value.
#   If no 64-bit integer could be declared, an error is generated.
#
# LICENSE
#
#   Copyright (c) 2011 Stefan Umbreit <s-umbreit@northwestern.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.
#
AC_DEFUN([_NU_FORTRAN_INT64_KIND_LOOP], [
    KINDINT64=''
    m4_for([kval], [1], [20],[],[
      AC_RUN_IFELSE([
        program int64kind
        integer(kind=kval) i

        open(20, file='conftestval', status='unknown')
        write(20, '(I3.3)') bit_size(i)
        close(20)
        end],
        [ AS_IF([test "x$KINDINT64" == "x" &&
          test "x`cat conftestval`" == "x064"],
          [KINDINT64=kval],[])
        ], [])
      ])
])

AC_DEFUN([NU_FORTRAN_INT64_KIND], [
    AC_LANG_PUSH([Fortran 77])
    KINDINT64=''
    AC_MSG_CHECKING([which KIND value defines a 64-bit integer])
    AC_RUN_IFELSE([
        program int64kind
        integer i
        open(20, file='conftestval', status='unknown')
        write(20, '(I2.2)') selected_int_kind(18)
        close(20)
        end
    ], [KINDINT64=`cat conftestval`;
        KINDINT64=`expr 0 + $KINDINT64`], [_NU_FORTRAN_INT64_KIND_LOOP])
    AS_IF([test "x$KINDINT64" == "x" || test "x$KINDINT64" == "x-1"],
          [AC_MSG_RESULT([could not find one])
           KINDINT64=-1
           ],
          [AC_MSG_RESULT([$KINDINT64])])
    AC_SUBST([KINDINT64])
    AC_LANG_POP([])
])
