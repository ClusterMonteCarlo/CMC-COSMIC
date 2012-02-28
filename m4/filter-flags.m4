AC_DEFUN([FILTER_VAR], [
        flags=$$1
        rflags="m4_shift($@)"
        new_flags=""
        for flag in $flags
        do
          for rflag in $rflags
          do
            flag="${flag##$rflag}"
          done
          new_flags="${new_flags}${flag:+ $flag}"
        done
        $1="$new_flags"
])
