#!gmake
#
# Gromacs Makefile $Id$
#
ifeq ($(GMXHOME),)
all:
	echo "No GROMACS environment, Please source your GMXRC file"
else
include		$(GMXHOME)/src/makef/Makefile.def
#
ifneq ($(USE_XWIN),no)
SRCDIRS	=	include	gmxlib fftw	kernel	tools	ngmx
else
SRCDIRS	=	include	gmxlib fftw	kernel	tools
endif
LOCALDIRS =	local
ALLDIRS	=	$(SRCDIRS)	$(LOCALDIRS)
SRCFILES =	Makefile makef mkhtml mkcompl
SRCDIST	=	$(SRCDIRS)	$(SRCFILES)
AUXFILES =	NOGMX
AUXDIRS	=	top tutor html
DIST	=	$(AUXFILES) $(foreach dir,$(SRCDIST),src/$(dir)) $(AUXDIRS)
TARDIR	=	$(GMXHOME)
TARCMD  =	tar cvf 
#		/usr/local/bin/tar cvf 
#TARFLG  =	--exclude=CVS
COPYRGT	=	$(GMXHOME)/src/local/copyrgt
PROGRAMS =	(cd bin/s10; ls; cd ../..)

all:	oclean dirs
	$(foreach dir,$(SRCDIRS),cd $(dir); $(MAKE); cd ..;)
	./mkcompl

etags:
	find . -name "*.[chCHm4]" -print | etags -

it:
	$(foreach dir,$(ALLDIRS),cd $(dir); $(MAKE); cd ..;)

dist:	oclean
	cd $(GMXHOME); $(TARCMD) $(TARDIR)/gmxcore.tar $(TARFLG) $(DIST)
	cd $(GMXHOME); $(TARCMD) $(TARDIR)/local.tar $(TARFLG) src/$(LOCALDIRS)

# rule to make dirs where binaries and libary go:
dirs:
	\mkdir -p $(BINDIR) $(LIBDIR); 

man:
	./mknroff
html:
	./mkhtml
complete:
	./mkcompl

clean:
	$(foreach dir,$(SRCDIRS),cd $(dir); $(MAKE) clean; cd ..;)

copyrgt:
	$(foreach dir,$(ALLDIRS),echo $(dir); cd $(dir); $(COPYRGT) *.c; $(COPYRGT) *.h; cd $(GMXHOME)/src;)

oclean:
	$(foreach dir,$(SRCDIRS),cd $(dir); $(MAKE) oclean; cd ..;)
endif
