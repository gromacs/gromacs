AUTOMAKE_OPTIONS=foreign
AM_CFLAGS = -I../sz/include -I../zlib -I../zstd $(GSL_CFLAGS) $(GSL_STATIC)
AM_FCFLAGS = -I../sz/include -I../sz -I../zlib -I../zstd $(GSL_CFLAGS) $(GSL_STATIC)

if WRITESTATS
AM_CFLAGS += -DHAVE_WRITESTATS
endif

if OPENMP
AM_CFLAGS += -fopenmp
endif

if FORTRAN
bin_PROGRAMS=testint_compress testint_decompress testfloat_compress testfloat_decompress testdouble_compress testdouble_decompress sz testdouble_compress_f testdouble_decompress_f
else
bin_PROGRAMS=testint_compress testint_decompress testfloat_compress testfloat_decompress testdouble_compress testdouble_decompress sz
endif

if TIMECMPR
bin_PROGRAMS += testfloat_compress_ts testfloat_decompress_ts testdouble_compress_ts testdouble_decompress_ts
endif

#if GSL
#bin_PROGRAMS += sz_gsl
#endif
	
testint_compress_SOURCES=testint_compress.c
testint_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testint_decompress_SOURCES=testint_decompress.c
testint_decompress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_compress_SOURCES=testfloat_compress.c
testfloat_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_decompress_SOURCES=testfloat_decompress.c
testfloat_decompress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testdouble_compress_SOURCES=testdouble_compress.c
testdouble_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testdouble_decompress_SOURCES=testdouble_decompress.c
testdouble_decompress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
if FORTRAN
testdouble_compress_f_SOURCES=testdouble_compress_f.f90
testdouble_compress_f_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testdouble_decompress_f_SOURCES=testdouble_decompress_f.f90
testdouble_decompress_f_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
endif
if TIMECMPR
testfloat_compress_ts_SOURCES=testfloat_compress_ts.c
testfloat_compress_ts_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_decompress_ts_SOURCES=testfloat_decompress_ts.c
testfloat_decompress_ts_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testdouble_compress_ts_SOURCES=testdouble_compress_ts.c
testdouble_compress_ts_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testdouble_decompress_ts_SOURCES=testdouble_decompress_ts.c
testdouble_decompress_ts_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
endif

if RANDOMACCESS
sz_SOURCES=sz_randomaccess.c
sz_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
else
if PASTRI
sz_SOURCES=sz_pastri.c
sz_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
else
if GSL
sz_SOURCES=sz_gsl.c
sz_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a $(GSL_LIBS) $(GSL_HDR) $(GSL_LIB) -lm
else
if OPENMP
sz_SOURCES=sz_openmp.c
sz_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
else
sz_SOURCES=sz.c
sz_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
endif
endif
endif
endif
