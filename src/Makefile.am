## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

SUBDIRS = gmxlib mdlib kernel tools ngmx contrib
# don't do anything for the local directory

mdrun:		
	(cd ${top_builddir}/src/gmxlib && $(MAKE))
	(cd ${top_builddir}/src/mdlib && $(MAKE))
	(cd ${top_builddir}/src/kernel && $(MAKE) mdrunXXX_SUFFIX_XXX)

install-mdrun:
	cd ${top_builddir}/src/kernel 
	make install-mdrun

CLEANFILES = *~ \\\#*