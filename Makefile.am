## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is automatically from Makefile.am by automake.

# N.B.: To be able to use dependency tracking on all hosts (not only w. gmake/gcc)
# and to use the suffix in library names we need the development versions of
# automake. Presently, this means automake version 1.4f. In Groningen, you
# can find that under ~lindahl/gnu/bin , or download a copy from www.gromacs.org

AUTOMAKE_OPTIONS = foreign

SUBDIRS = src admin scripts share include man

#
# Most files in the config subdir is included automatically when
# we issue "make dist", but some versions of automake seem to have
# problems with it, so just include them all...

EXTRA_DIST = config/depcomp 	config/ltconfig 	config/ltcf-c.sh   \
	     config/ltcf-f77.sh config/config.guess	config/config.sub  \
	     config/install-sh  config/missing          config/ltmain.sh   \
	     config/mkinstalldirs                               

install-hook:
	@echo ""
	@echo "GROMACS is installed under $(prefix)."
	@echo "Make sure to update your PATH and MANPATH to find the"
	@echo "programs and unix manual pages, and possibly LD_LIBRARY_PATH"
	@echo "or /etc/ld.so.conf if you are using dynamic libraries."
	@echo ""
	@echo "If you want links to the executables in /usr/local/bin,"
	@echo "you can issue \"make links\" now."

# Shortcuts to construct and install the mdrun executable

mdrun:		
	(cd $(top_builddir)/src/gmxlib && $(MAKE) ; exit 0)
	(cd $(top_builddir)/src/mdlib && $(MAKE) ; exit 0)
	(cd $(top_builddir)/src/kernel && $(MAKE) mdrun ; exit 0)

install-mdrun:
	(cd $(top_builddir)/src/kernel && $(MAKE) install-mdrun ; exit 0)

links:
	programs=$(bindir)/* 
	(cd /usr/local/bin && for i in $(programs); do \
	   (test ! -e $$i && $(LN_S) $$i . ; exit 0); \
	done)

rpm:
	mkdir -p rpmdir/BUILD \
	  rpmdir/SPECS \
	  rpmdir/SOURCES \
	  rpmdir/SRPMS \
	  rpmdir/RPMS/i386
	cp -f $(PACKAGE)-$(VERSION).tar.gz rpmdir/SOURCES
	rpm --buildroot $$PWD/rpmdir/buildroot --clean -bb $(top_srcdir)/admin/$(PACKAGE).spec --define "_topdir $$PWD/rpmdir" 
	@echo "========================================================"
	@echo "Finished - packages are located under rpmdir/RPMS/i386 !"


# Dont try to make rpm and mpi-rpm in parallel.

mpi-rpm:	
	mkdir -p rpmdir/BUILD \
	  rpmdir/SPECS \
	  rpmdir/SOURCES \
	  rpmdir/SRPMS \
	  rpmdir/RPMS/i386
	cp -f $(PACKAGE)-$(VERSION).tar.gz rpmdir/SOURCES/$(PACKAGE)-mpi-$(VERSION).tar.gz
	rpm --buildroot $$PWD/rpmdir/buildroot-mpi --clean -bb rpm_build/SPECS/$(PACKAGE)-mpi.spec --define "_topdir $$PWD/rpmdir"
	@echo "========================================================"
	@echo "Finished - packages are located under rpmdir/RPMS/i386 !"

CLEANFILES = *~ \\\#* 

DISTCLEANFILES = libtool config.cache config.log









































































