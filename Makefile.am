## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#
AUTOMAKE_OPTIONS = foreign

#
# The real work is done in our subdirs
#
SUBDIRS = src top html tutor

mdrun:		
	(cd ${top_builddir}/src/gmxlib && $(MAKE))
	(cd ${top_builddir}/src/mdlib && $(MAKE))
	(cd ${top_builddir}/src/kernel && $(MAKE) mdrunXXX_SUFFIX_XXX)

install-mdrun:
	(cd ${top_builddir}/src/kernel && $(MAKE) install-mdrun)

CLEANFILES = *~ \\\#*








































































