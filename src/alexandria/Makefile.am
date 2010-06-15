## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

# These programs are not compiled or installed by default - you will have to 
# issue "make <program>" and copy the binary to the correct location yourself! 
# Add new entries in Makefile.am!

lib_LTLIBRARIES = libalexandria@LIBSUFFIX@.la

pkgconfigdir = ${libdir}/pkgconfig
pkgconfig_DATA = libalexandria@LIBSUFFIX@.pc

LDADD = $(lib_LTLIBRARIES)

libalexandria@LIBSUFFIX@_la_LIBADD         = ../gmxlib/libgmx@LIBSUFFIX@.la ../mdlib/libmd@LIBSUFFIX@.la ../tools/libgmxana@LIBSUFFIX@.la
libalexandria@LIBSUFFIX@_la_DEPENDENCIES   = ../gmxlib/libgmx@LIBSUFFIX@.la ../mdlib/libmd@LIBSUFFIX@.la ../tools/libgmxana@LIBSUFFIX@.la
libalexandria@LIBSUFFIX@_la_LDFLAGS        = -version-info @SHARED_VERSION_INFO@

bin_PROGRAMS 	=  	pdtest		gauss_test	slater_test

# Note: you don't have to list sources for "prog" if it is the single file prog.c

libalexandria@LIBSUFFIX@_la_SOURCES = \
			xml_util.c	xml_util.h	\
			poldata.c 	poldata.h \
			poldata_xml.c	poldata_xml.h \
			gaussian_integrals.c 	gaussian_integrals.h	\
			gmx_simple_comm.h	gmx_simple_comm.c	\
			slater_integrals.cpp 	slater_integrals.h	\
			slater_low.h \
			DSlater_1S_1S.cpp DSlater_1S_2S.cpp DSlater_1S_3S.cpp \
			DSlater_1S_4S.cpp DSlater_1S_5S.cpp DSlater_1S_6S.cpp \
			DSlater_2S_2S.cpp DSlater_2S_3S.cpp DSlater_2S_4S.cpp \
			DSlater_2S_5S.cpp DSlater_2S_6S.cpp DSlater_3S_3S.cpp \
			DSlater_3S_4S.cpp DSlater_3S_5S.cpp DSlater_3S_6S.cpp \
			DSlater_4S_4S.cpp DSlater_4S_5S.cpp DSlater_4S_6S.cpp \
			DSlater_5S_5S.cpp DSlater_5S_6S.cpp DSlater_6S_6S.cpp \
			Slater_1S_1S.cpp Slater_1S_2S.cpp Slater_1S_3S.cpp \
			Slater_1S_4S.cpp Slater_1S_5S.cpp Slater_1S_6S.cpp \
			Slater_2S_2S.cpp Slater_2S_3S.cpp Slater_2S_4S.cpp \
			Slater_2S_5S.cpp Slater_2S_6S.cpp Slater_3S_3S.cpp \
			Slater_3S_4S.cpp Slater_3S_5S.cpp Slater_3S_6S.cpp \
			Slater_4S_4S.cpp Slater_4S_5S.cpp Slater_4S_6S.cpp \
			Slater_5S_5S.cpp Slater_5S_6S.cpp Slater_6S_6S.cpp 


slater_test_SOURCES	= slater_test.c

gauss_test_SOURCES	= gauss_test.c

pdtest_SOURCES	= pdtest.c

CLEANFILES   = 	*~ \\\#*










