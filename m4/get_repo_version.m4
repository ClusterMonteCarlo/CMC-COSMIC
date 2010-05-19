# Tries to determine the revision number of the repository.
# First it tries to execute subversion or mercurial if they are installed and 
# extract a revision number. On success the output variable is set to the 
# revision number, otherwise to 'exported'. Should be broken down into smaller 
# macros, e.g., GET_SVN_VERSION, GET_HG_VERSION with either a version number or
# "exported".

AC_DEFUN([GET_REPO_VERSION], [
         AC_CHECK_PROG([HAS_SVN], [svnversion], [yes], [no])
         AC_CHECK_PROG([HAS_HG], [hg], [yes], [no])
         VERSION="exported"
         AS_IF([test "x$HAS_SVN" == "xyes"], [VERSION=`svnversion`])
         AS_IF([test "x$VERSION" == "xexported"], [
               AS_IF([test "x$HAS_HG" == "xyes"], [VERSION=`hg id -n 2>/dev/null`])])
         AS_IF([test "x$VERSION" == "x"], [VERSION="exported"])

         $1=$VERSION
         ])

