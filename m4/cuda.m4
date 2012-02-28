AC_DEFUN([CUDA_FLAGS_APPEND], [
        AC_MSG_CHECKING([for CUDA include directories])
        CUDA_CONFIG_VAR([INCLUDES])
        AC_MSG_RESULT([${CUDA_CONFIG_INCLUDES}])
        AC_MSG_CHECKING([for CUDA libraries])
        CUDA_CONFIG_VAR([LIBRARIES])
        AC_MSG_RESULT([${CUDA_CONFIG_LIBRARIES}])

        CPPFLAGS="${CPPFLAGS} ${CUDA_CONFIG_INCLUDES}"
        LDFLAGS="${LDFLAGS} ${CUDA_CONFIG_LIBRARIES}"
])

AC_DEFUN([CUDA_CONFIG_VAR], [
dnl        AC_REQUIRE([AC_PROG_SED])
              
        touch conftest_cuda.cu
        [$]{NVCC} -E conftest_cuda.cu -dryrun 2>&1 |
        [$]{SED} "s/#$\s*gcc.*//" |
        [$]{SED} "s/#$\s*//"|
        [$]{SED} 's/=\(.*\)$/=\"\1\"/' > conftest.config

        echo 'echo [$]{AS_TR_SH([$1])}'|
        cat conftest.config - > conftest2.config 

        [CUDA_CONFIG_]AS_TR_SH([$1])=`[$]{SHELL} conftest2.config`
])

