AM_CPPFLAGS = -I../include

pkgconfig_DATA = s2kit.pc

lib_LTLIBRARIES = libs2kit.la

libs2kit_la_SOURCES = \
	cospmls.c \
	FST_semi.c \
	FST_semi_fly.c \
	FST_semi_memo.c \
	makeweights.c \
	naive_synthesis.c \
	pmls.c \
	primitive.c \
	seminaive.c
libs2kit_la_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS)
libs2kit_la_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS)

