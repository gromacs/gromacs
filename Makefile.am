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
# we issue "make dist", but not depcomp. This just includes that file...
#
EXTRA_DIST = config/depcomp


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

#
# Apart from normal things like .o, things matching this are removed
# (The second one removes files beginning with a #)
CLEANFILES = *~ \\\#*








































































