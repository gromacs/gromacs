## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#
AUTOMAKE_OPTIONS = foreign

#
# The real work is done in our subdirs - automake iterates through these
#
SUBDIRS = src top html tutor man

#
# Most files in the config subdir is included automatically when
# we issue "make dist", but some versions of automake seem to have
# problems with it, so we include them all...

EXTRA_DIST = config/depcomp 	config/ltconfig 	config/ltcf-c.sh   \
             config/ltcf-f77.sh config/config.guess	config/config.sub  \
	     config/install-sh  config/missing          config/ltmain.sh   \
             config/mkinstalldirs

#
# This is a shortcut to construct the mdrun executable by first
# building the two libraries and then only mdrun in kernel.
# The second one installs it.
#
mdrun:		
	(cd ${top_builddir}/src/gmxlib && $(MAKE))
	(cd ${top_builddir}/src/mdlib && $(MAKE))
	(cd ${top_builddir}/src/kernel && $(MAKE) mdrunXXX_SUFFIX_XXX)

install-mdrun:
	(cd ${top_builddir}/src/kernel && $(MAKE) install-mdrun)

links:
	(cd /usr/local/bin && $(LN_S) ${bindir}/* .)

#
# Apart from normal things like .o, things matching this are removed
# (The second one removes files beginning with a #)
CLEANFILES = *~ \\\#* 
DISTCLEANFILES = libtool config.cache config.log









































































