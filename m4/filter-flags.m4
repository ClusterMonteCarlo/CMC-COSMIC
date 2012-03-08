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

AC_DEFUN([REMOVE_DUPS], [
   flags=$$1
   
   $1=`echo $flags|
   gawk 'BEGIN {RS=" "} 
   !([$]1 in arr) {flag=[$]1; arr[[$flag]]=1; printf "%s ", $flag}'`
])

