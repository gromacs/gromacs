## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is automatically from Makefile.am by automake.

# N.B.: To be able to use dependency tracking on all hosts (not only w. gmake/gcc)
# and to use the suffix in library names we need the development versions of
# automake. Presently, this means automake version 1.4f. In Groningen, you
# can find that under ~lindahl/gnu/bin , or download a copy from www.gromacs.org

AUTOMAKE_OPTIONS = foreign

SUBDIRS = include src admin scripts share man

#
# Most files in the config subdir is included automatically when
# we issue "make dist", but not depcomp

EXTRA_DIST = config/depcomp 	                               

install-exec-hook:
	@echo ""
	@echo "GROMACS is installed under $(prefix)."
	@echo "Make sure to update your PATH and MANPATH to find the"
	@echo "programs and unix manual pages, and possibly LD_LIBRARY_PATH"
	@echo "or /etc/ld.so.conf if you are using dynamic libraries."
	@echo ""
	@echo "Please run \"make tests\" now to verify your installation."
	@echo ""
	@echo "If you want links to the executables in /usr/local/bin,"
	@echo "you can issue \"make links\" now."

# Shortcuts to construct and install the mdrun executable

mdrun:		
	(cd $(top_builddir)/src/gmxlib && $(MAKE) ; exit 0)
	(cd $(top_builddir)/src/mdlib && $(MAKE) ; exit 0)
	(cd $(top_builddir)/src/kernel && $(MAKE) mdrun ; exit 0)

install-mdrun:
	(cd $(top_builddir)/src/gmxlib && $(MAKE) install ; exit 0)
	(cd $(top_builddir)/src/mdlib && $(MAKE) install ; exit 0)
	(cd $(top_builddir)/src/kernel && $(MAKE) install-mdrun ; exit 0)

links:
	cd $(bindir) && programs=`ls` && cd /usr/local/bin && \
	for i in $$programs; do \
	   (test ! -f $$i && $(LN_S) $(bindir)/$$i . ; exit 0); \
	done

tests:
	(if test -d "gmxtest"; then cd "gmxtest"; ./gmxtest.pl all; cd ..; \
	else echo "No gmxtest directory found. Please download and unpack it here.";\
	fi)


RPMARCH=`rpm --showrc | awk '/^build arch/ {print $$4}'`

#
# Use --buildroot ???
rpm:
	mkdir -p rpmdir/BUILD \
	  rpmdir/SPECS \
	  rpmdir/SOURCES \
	  rpmdir/SRPMS \
	  rpmdir/RPMS/$(RPMARCH)
	cp -f $(PACKAGE)-$(VERSION).tar.gz rpmdir/SOURCES
	rpmbuild --clean -ba $(top_srcdir)/admin/$(PACKAGE).spec --define "_topdir $$PWD/rpmdir" 
	@echo "============================================================"
	@echo "Finished - the packages are in rpmdir/RPMS and rpmdir/SRPMS!"


mpi-rpm:	
	mkdir -p rpmdir/BUILD \
	  rpmdir/SPECS \
	  rpmdir/SOURCES \
	  rpmdir/SRPMS \
	  rpmdir/RPMS/$(RPMARCH)
	cp -f $(PACKAGE)-$(VERSION).tar.gz rpmdir/SOURCES/$(PACKAGE)-$(VERSION).tar.gz
	rpmbuild --clean -bb $(top_srcdir)/admin/$(PACKAGE)-mpi.spec --define "_topdir $$PWD/rpmdir"
	@echo "============================================"
	@echo "Finished - the packages are in rpmdir/RPMS!"

CLEANFILES = *~ \\\#* 

DISTCLEANFILES = libtool config.cache config.log









































































