## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\" -I/sw/include

LDADD = ../kernel/libgmxpreprocess@LIBSUFFIX@.la ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la ../tools/libgmxana@LIBSUFFIX@.la


EXTRA_DIST              = README  

# These programs are not compiled or installed by default - you will have to 
# issue "make <program>" and copy the binary to the correct location yourself! 
# Add new entries in Makefile.am!
#extra_PROGRAMS  =       prfn	mkice	total	hexamer 	optwat do_shift g_anavel anaf pmetest

EXTRA_PROGRAMS 	=  	copyrgt 	mkyaw 		\
			options 	hrefify 	\
			addquote 	tune_dip	\
			ehole 		gen_table	\
			tune_pol	gentop

# Note: you don't have to list sources for "prog" if it is the single file prog.c

#total_SOURCES		= total.f

nodist_ehole_SOURCES		= ehole.c ehdata.c ehanal.c

#pmetest_SOURCES   	= pmetest.c

nodist_tune_pol_SOURCES	= tune_pol.c tune_pol.h pol_xml.c \
			gentop_matrix.c	gentop_matrix.h	\
			poldata.c merge_xml.c

nodist_tune_dip_SOURCES 	= tune_dip.c 	\
			gentop_matrix.c		gentop_matrix.h	\
			gentop_core.c		gentop_core.h	\
			gentop_nm2type.c 	gentop_nm2type.h	\
			gentop_qalpha.h		gentop_qalpha.c	\
			gentop_qgen.c 		gentop_qgen.h \
			gentop_eemprops.c 	gentop_eemprops.h\
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

nodist_gen_table_SOURCES       = gen_table.c slater_integrals.cpp \
			  slater_integrals.h \
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

nodist_gentop_SOURCES		= gentop.c	gentop_qgen.c	\
			gentop_matrix.c	gentop_matrix.h	\
			gentop_core.c	gentop_core.h	\
			gentop_nm2type.c gentop_nm2type.h	\
			gentop_qalpha.h	gentop_qalpha.c	\
			gentop_eemprops.c gentop_eemprops.h\
			pol_xml.c tune_pol.h poldata.c \
			molprop.c molprop.h molprop_xml.c \
			slater_integrals.cpp slater_integrals.h \
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

CLEANFILES   = 	*~ \\\#*










