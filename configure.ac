# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

AC_INIT([tagtoxfm_bspline.cc],[0.9],[John G. Sled jgsled@sickkids.ca])
AC_CONFIG_SRCDIR([tagtoxfm_bspline.cc])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH


# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL

# Check for lapack
AC_MSG_CHECKING([for lapack])
LIBS="-llapack $libs"
AC_TRY_LINK([dsysv_();]
	    AC_MSG_ERROR(cannot find lapack library))
AC_MSG_RESULT([yes])

# Checks for libraries.  See m4/README.
mni_REQUIRE_BICPL

#AC_LANG_PUSH(C++)
#AC_MSG_CHECKING([for blitz])
#AC_TRY_LINK([#include <mniVolume.h>],
#            [mniVolume::mniVolume();]
#            AC_MSG_ERROR(cannot find oobicpl library))
#AC_MSG_RESULT([yes])
#
#AC_LANG_POP


# Finish up by writing output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
