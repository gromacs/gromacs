## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is automatically from Makefile.am by automake.

SUBDIRS = gmxlib mdlib kernel tools ngmx contrib
# don't do anything for the local directory

#EXTRA_DIST = config.win

mdrun:		
	(cd $(top_builddir)/src/gmxlib && $(MAKE))
	(cd $(top_builddir)/src/mdlib && $(MAKE))
	(cd $(top_builddir)/src/kernel && $(MAKE) mdrun)

install-mdrun:
	(cd $(top_builddir)/src/gmxlib && $(MAKE) install ; exit 0)
	(cd $(top_builddir)/src/mdlib && $(MAKE) install ; exit 0)
	(cd $(top_builddir)/src/kernel && $(MAKE) install-libLTLIBRARIES ; exit 0)
	(cd $(top_builddir)/src/kernel && $(MAKE) install-mdrun)

CLEANFILES = *~ \\\#*
