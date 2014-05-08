
dnl If adding rcX to version, be sure to separate with a '-' 
AC_INIT(visual-image-sort, [0.1])
AC_CONFIG_SRCDIR(compute_asift_keypoints.cpp)
AM_INIT_AUTOMAKE

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES(yes)])


dnl Determine CC and preset CFLAGS
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([eigen3/Eigen/Cholesky])
AC_LANG_POP([C++])

AC_SEARCH_LIBS(png_init_io, png)

AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)

AC_OUTPUT(Makefile)

cat << EOF

Your build configuration:

	CPPFLAGS = $CPPFLAGS
	CXXFLAGS = $CXXFLAGS
	CFLAGS = $CFLAGS
	LDFLAGS = $LDFLAGS
	prefix: $prefix
	version: $VERSION

EOF

