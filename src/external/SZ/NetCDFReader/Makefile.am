AUTOMAKE_OPTIONS=foreign

include_HEADERS=include/NetCDFReader.h

lib_LTLIBRARIES=libnetcdfreader.la
if NETCDF
libnetcdfreader_la_CFLAGS=-g -I./include -DHAVE_NETCDF $(NETCDF_HDR)
else
libnetcdfreader_la_CFLAGS=-g -I./include
endif

libnetcdfreader_la_LDFLAGS = -version-info  0:1:0
libnetcdfreader_la_LIDADD=$(NETCDF_LIBS) $(NETCDF_LIB)
libnetcdfreader_la_SOURCES=src/NetCDFReader.c
