#
# Gromacs Makefile $Id$
#
include		$(GMXHOME)/src/makef/Makefile.def
#
ifneq ($(USE_XWIN),no)
DIRS	=	include	gmxlib fftw	kernel	tools	ngmx
else
DIRS	=	include	gmxlib fftw	kernel	tools
endif
LOCAL	=	local
LDIRS	=	$(DIRS)	$(LOCAL)
TELLDIRS=	$(LDIRS)
DDDIRS	=	$(DIRS)	Makefile makef mkhtml
FILES   =	NOGMX
SDIRS	=	$(foreach DIR,$(TELLDIRS),src/$(DIR))
AUX	=	top tutor html
DIST	=	$(FILES) $(foreach DIR,$(DDDIRS),src/$(DIR)) $(AUX)
GETIT	=	$(SDIRS)  src/makef
TDIR	=	$(GMXHOME)
TARCMD  =	tar cvf #/usr/local/bin/tar cvf 
#TARFLG  =	--exclude=CVS
COPYRGT	=	$(GMXHOME)/src/local/copyrgt
PROGRAMS=	(cd bin/s10; ls; cd ../..)

all:	oclean dirs
	$(foreach DIR,$(DIRS),cd $(DIR); $(MAKE); cd ..;)

etags:
	find . -name "*.[chCHm4]" -print | etags -

it:	
	$(foreach DIR,$(LDIRS),cd $(DIR); $(MAKE); cd ..;)

dist:	oclean
	cd $(GMXHOME); $(TARCMD) $(TDIR)/gmxcore.tar $(TARFLG) $(DIST)
	cd $(GMXHOME); $(TARCMD) $(TDIR)/local.tar $(TARFLG) src/$(LOCAL)

# rule to make dirs where binaries and libary go:
dirs:
	\mkdir -p $(BINDIR) $(LIBDIR); 

man:
	./mknroff
html:
	./mkhtml
compl:
	./mkcompl

clean:
	$(foreach dir,$(DIRS),cd $(dir); $(MAKE) clean; cd ..;)

copyrgt:
	$(foreach DIR,$(TELLDIRS),echo $(DIR); cd $(DIR); $(COPYRGT) *.c; $(COPYRGT) *.h; $(COPYRGT) *.cc; cd $(GMXHOME)/src;)

oclean: 
	$(foreach dir,$(DIRS),cd $(dir); $(MAKE) oclean; cd ..;)
