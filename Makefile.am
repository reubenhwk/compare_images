
AUTOMAKE_OPTIONS = 1.2 foreign subdir-objects

AM_CFLAGS = @OPENMP_CFLAGS@
AM_CXXFLAGS = @OPENMP_CFLAGS@

noinst_PROGRAMS = \
	demo_ASIFT

noinst_LIBRARIES = \
	libnumerics.a \
	libmatch.a \
	libio_png.a

demo_ASIFT_LDADD = \
	libnumerics.a \
	libmatch.a \
	libio_png.a

demo_ASIFT_SOURCES = \
	compute_asift_keypoints.cpp \
	compute_asift_keypoints.h \
	compute_asift_matches.cpp \
	compute_asift_matches.h \
	demo_ASIFT.cpp \
	demo_lib_sift.cpp \
	demo_lib_sift.h \
	domain.cpp \
	domain.h \
	filter.cpp \
	filter.h \
	flimage.cpp \
	flimage.h \
	fproj.cpp \
	fproj.h \
	frot.cpp \
	frot.h \
	library.cpp \
	library.h \
	numerics1.cpp \
	numerics1.h \
	orsa.cpp \
	orsa.h \
	splines.cpp \
	splines.h

libio_png_a_SOURCES = \
	lib/io_png/io_png.c \
	lib/io_png/io_png.h

libnumerics_a_SOURCES = \
	lib/numerics/numerics.cpp \
	lib/numerics/rodrigues.h \
	lib/numerics/homography.cpp \
	lib/numerics/vector.cpp \
	lib/numerics/computeH.cpp \
	lib/numerics/rodrigues.cpp \
	lib/numerics/numerics.h \
	lib/numerics/matrix.h \
	lib/numerics/matrix.cpp \
	lib/numerics/homography.h

libmatch_a_SOURCES = \
	lib/match/match.cpp \
	lib/match/match.h

