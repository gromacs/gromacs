## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#
SUBDIRS = scripts

include $(srcdir)/../Makefile.inc

# This is necessary for VPATH to work, 
# but it can probably be done in a nicer way...
INCLUDES = @INCLUDES@ -I$(top_srcdir)/src/include
LDFLAGS = @LDFLAGS@ -L${top_builddir}/src/gmxlib -L${top_builddir}/src/mdlib 
LDADD = ../gmxlib/libgmxXXX_SUFFIX_XXX.la ../mdlib/libmdXXX_SUFFIX_XXX.la

EXTRA_DIST              = README  programs.txt  gromacs-3.0.spec

# NB: The programs in contrib do not get double/mpi suffixes automatically,
# to make it easier for new developers to understand this file. If you want
# to have it, append "XXX_SUFFIX_XXX" to the program name.
# These programs are not compiled or installed by default - you will have to 
# issue "make <program>" and copy the binary to the correct location yourself! 
# Add new entries in Makefile.am!

EXTRA_PROGRAMS 		= sigeps mk6_n copyrgt mkyaw prfn options hrefify mkice \
                          total hexamer optwat addquote do_shift g_anavel my_dssp anaf 

copyrgt_SOURCES 	= copyrgt.c

sigeps_SOURCES 		= sigeps.c 

mk6_n_SOURCES 		= mk6_n.c 

my_dssp_SOURCES 	= my_dssp.c dsspcore.c

hrefify_SOURCES 	= hrefify.c

total_SOURCES 		= total.f

optwat_SOURCES 		= optwat.c

hexamer_SOURCES 	= hexamer.c

do_shift_SOURCES 	= do_shift.c

anaf_SOURCES 		= anaf.c

mkyaw_SOURCES 		= mkyaw.c

mkice_SOURCES 		= mkice.c

prfn_SOURCES 		= prfn.c

options_SOURCES		= options.c

g_anavel_SOURCES 	= g_anavel.c

addquote_SOURCES 	= addquote.c

# clean things explicitly, since the target names might have changed
CLEANFILES   = 	${bin_PROGRAMS}	${EXTRA_PROGRAMS}	\
		*_mpi *_d *~ \\\#*
