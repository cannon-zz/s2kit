AM_CPPFLAGS = -I../src/include

noinst_PROGRAMS = \
	test_conv_semi_memo \
	test_conv_semi_fly \
	test_naive \
	test_s2_semi_fly \
	test_s2_semi_memo \
	test_s2_semi_memo_for \
	test_s2_semi_memo_inv \
	test_semi

TESTS = $(noinst_PROGRAMS)

test_conv_semi_memo_SOURCES = test_conv_semi_memo.c
test_conv_semi_memo_LDADD = ../src/lib/libs2kit.la

test_conv_semi_fly_SOURCES = test_conv_semi_fly.c
test_conv_semi_fly_LDADD = ../src/lib/libs2kit.la

test_naive_SOURCES = test_naive.c csecond.h csecond.c
test_naive_LDADD = ../src/lib/libs2kit.la

test_s2_semi_fly_SOURCES = test_s2_semi_fly.c csecond.h csecond.c
test_s2_semi_fly_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS)
test_s2_semi_fly_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS)
test_s2_semi_fly_LDADD = ../src/lib/libs2kit.la

test_s2_semi_memo_SOURCES = test_s2_semi_memo.c csecond.h csecond.c
test_s2_semi_memo_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS)
test_s2_semi_memo_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS)
test_s2_semi_memo_LDADD = ../src/lib/libs2kit.la

test_s2_semi_memo_for_SOURCES = test_s2_semi_memo_for.c csecond.h csecond.c
test_s2_semi_memo_for_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS)
test_s2_semi_memo_for_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS)
test_s2_semi_memo_for_LDADD = ../src/lib/libs2kit.la

test_s2_semi_memo_inv_SOURCES = test_s2_semi_memo_inv.c csecond.h csecond.c
test_s2_semi_memo_inv_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS)
test_s2_semi_memo_inv_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS)
test_s2_semi_memo_inv_LDADD = ../src/lib/libs2kit.la

test_semi_SOURCES = test_semi.c csecond.h csecond.c
test_semi_CFLAGS = $(AM_CFLAGS) $(FFTW_CFLAGS)
test_semi_LDFLAGS = $(AM_LDFLAGS) $(FFTW_LIBS)
test_semi_LDADD = ../src/lib/libs2kit.la
