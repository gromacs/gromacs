## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is automatically from Makefile.am by automake.

# N.B.: To be able to use dependency tracking on all hosts (not only w. gmake/gcc)
# and to use the suffix in library names we need the development versions of
# automake. Presently, this means automake version 1.4f. In Groningen, you
# can find that under ~lindahl/gnu/bin , or download a copy from www.gromacs.org

AUTOMAKE_OPTIONS = foreign

SUBDIRS = include src share admin scripts man 

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

# Shortcuts to construct and install the mdrun executable

mdrun:		
	(cd $(top_builddir)/src/gmxlib && $(MAKE))
	(cd $(top_builddir)/src/mdlib && $(MAKE))
	(cd $(top_builddir)/src/kernel && $(MAKE) mdrun)

install-mdrun:
	(cd $(top_builddir)/src/kernel && $(MAKE) install-mdrun)

links:
	programs=$(bindir)/* 
	(cd /usr/local/bin && for i in $(programs); do \
	   (test ! -e $$i && $(LN_S) $$i .); \
	done)

rpm:	
	mkdir -p rpm_buildroot/BUILD \
	  rpm_buildroot/SPECS \
	  rpm_buildroot/SOURCES \
	  rpm_buildroot/SRPMS \
	  rpm_buildroot/RPMS/i386
	cp -f $(PACKAGE)-$(VERSION).tar.gz rpm_buildroot/SOURCES
	cp -f admin/$(PACKAGE).spec rpm_buildroot/SPECS
	rpm --rcfile /usr/lib/rpm/rpmrc:/etc/rpmrc:admin/rpmrc --clean --ba $(PACKAGE).spec
	cp -f rpm_buildroot/RPMS/i386/*rpm $(top_builddir)
	rm -rf rpm_buildroot

mpi-rpm:	
	mkdir -p rpm_buildroot/BUILD \
	  rpm_buildroot/SPECS \
	  rpm_buildroot/SOURCES \
	  rpm_buildroot/SRPMS \
	  rpm_buildroot/RPMS/i386
	cp -f $(PACKAGE)-$(VERSION).tar.gz rpm_buildroot/SOURCES/$(PACKAGE)-mpi-$(VERSION).tar.gz
	cp -f admin/$(PACKAGE)-mpi.spec rpm_buildroot/SPECS
	rpm --rcfile /usr/lib/rpm/rpmrc:/etc/rpmrc:admin/rpmrc --clean --ba $(PACKAGE)-mpi.spec
	cp -f rpm_buildroot/RPMS/i386/*rpm $(top_builddir)
	rm -rf rpm_buildroot

CLEANFILES = *~ \\\#* 

DISTCLEANFILES = libtool config.cache config.log









































































