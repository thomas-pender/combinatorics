#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([combinatorics], [1.0], [thomaspender1992@gmail.com])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([src/diophantine.c])
AC_CONFIG_HEADERS([src/config.h])

AM_INIT_AUTOMAKE([subdir-objects])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AM_PROG_AR

# LTLibraries
LT_PREREQ([2.4.6])
LT_INIT([static disable-shared])

#-------------------------------------------------
# doxygen
#-------------------------------------------------
DX_INIT_DOXYGEN([$PACKAGE_NAME],[Doxyfile],doc)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

if test "x${DX_FLAG_doc}" = "x1"; then
_doxy=yes
fi

AM_CONDITIONAL([DOXY_], [test "x${_doxy}" = xyes])
#-------------------------------------------------

# Checks for libraries.
have_containers=no
AC_SEARCH_LIBS([arrays_new], [containers], [have_containers=yes])

# Checks for header files.
if test "x${have_containers}" = xyes; then
   AC_CHECK_HEADERS([containers/containers.h], [], have_containers=no)
fi

if test "x${have_containers}" = xno; then
   AC_MSG_ERROR([
----------------------------------------
The containers library must be installed
and in the search path of the compiler.
Check 'config.log' for more information.
----------------------------------------])
fi

AC_SUBST([LIBS])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile Doxyfile])
AC_OUTPUT

cat << EOF

============================================================

${PACKAGE_NAME} version ${PACKAGE_VERSION}

Prefix: ${prefix}
C compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

Package features:
EOF

if test "x${_doxy}" = xyes; then
cat << EOF
    - Doxygen documentation enabled.
EOF
else
cat << EOF
    - Doxygen documentation disabled.
EOF
fi

cat << EOF

Compile and install with 'make && sudo make install'.

============================================================

EOF
