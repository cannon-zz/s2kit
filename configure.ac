AC_INIT([s2kit],[1.0],[kipp.cannon@ligo.org])
AC_COPYRIGHT([])
# a file whose existance can be used to use to check that we are in the
# top-level directory of the source tree
AC_CONFIG_SRCDIR([gnuscripts/ltmain.sh])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AM_INIT_AUTOMAKE([1.9 tar-ustar foreign])
AUTOMAKE_OPTIONS="${AUTOMAKE_OPTIONS} --add-missing --copy"
AC_CONFIG_FILES([ \
	Makefile \
	src/Makefile \
	src/bin/Makefile \
	src/include/Makefile \
	src/include/s2kit/Makefile \
	src/lib/Makefile \
	src/lib/s2kit.pc \
	test/Makefile
])

# set pkg-config module path
AC_SUBST([pkgconfigdir],[${libdir}/pkgconfig])

# Checks for programs
AC_PROG_CC_C99
AC_PROG_INSTALL
LT_INIT
PKG_PROG_PKG_CONFIG()

# Check for math library
AC_CHECK_LIB([m], [main], , [AC_MSG_ERROR([cannot find the math library])])

# Check for FFTW
PKG_CHECK_MODULES([FFTW], [fftw3], , [echo "Not Found!" ; exit 1])
AC_SUBST([FFTW_CFLAGS])
AC_SUBST([FFTW_LIBS])

# CFLAGS
AX_CFLAGS_WARN_ALL([AM_CFLAGS])
AC_SUBST([AM_CFLAGS])

# Generate configured output files
AC_PROG_MAKE_SET
AC_OUTPUT
