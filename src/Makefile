#
# Gromacs Makefile $Id$
#
include		$(GMXHOME)/src/makef/Makefile.def
#
#this is a dirty hack to stop gmake from trying to go to directory 'fftw'
#when it is not present. THIS IS NOT FOOL-PROOF!
ifeq ($(USE_SGI_FFT),yes)
DIRS	=	include	gmxlib kernel	tools	ngmx
else
DIRS	=	include	gmxlib fftw	kernel	tools	ngmx
endif
LOCAL	=	local
LDIRS	=	$(DIRS)	$(LOCAL)
TELLDIRS=	$(LDIRS)
DDDIRS	=	$(DIRS)	Makefile makef mkhtml
FILES   =	NOGMX
SDIRS	=	$(foreach DIR,$(TELLDIRS),src/$(DIR))
DIST	=	$(foreach DIR,$(DDDIRS),src/$(DIR))
AUX	=	top tutor
DIST	=	$(FILES) $(foreach DIR,$(DDDIRS),src/$(DIR)) $(AUX)
BENCH	=	tests/bench
TESTS	=	tests/mdrun
GETIT	=	$(SDIRS)  src/makef
TDIR	=	$(GMXHOME)
TARCMD  =	/usr/local/bin/tar cvf 
TARFLG  =	--exclude=CVS
COPYRGT	=	$(GMXHOME)/src/local/copyrgt
PROGRAMS=	(cd bin/s10; ls; cd ../..)

all:	dirs
	$(foreach DIR,$(DIRS),cd $(DIR); $(MAKE); cd ..;)

it:	
	$(foreach DIR,$(LDIRS),cd $(DIR); $(MAKE); cd ..;)

cdist:	clean
	cd $(GMXHOME); $(TARCMD) $(TDIR)/gmxcore.tar $(TARFLG) $(DIST)
	cd $(GMXHOME); $(TARCMD) $(TDIR)/local.tar $(TARFLG) src/$(LOCAL)

bdist:
	cd $(GMXHOME)/$(BENCH); cleanit
	cd $(GMXHOME); $(TARCMD) $(TDIR)/gmxbench.tar $(TARFLG) $(BENCH)

tdist:
	cd $(GMXHOME)/$(TESTS); cleanit
	cd $(GMXHOME); $(TARCMD) $(TDIR)/gmxtests.tar $(TARFLG) $(TESTS)

dist:	cdist	bdist	tdist


# rule to make dirs where binaries and libary go:
dirs:
	\mkdir -p $(BINDIR) $(LIBDIR); 

man:
	./mknroff
html:
	./mkhtml

clean:
	$(foreach dir,$(DIRS),cd $(dir); $(MAKE) clean; cd ..;)

copyrgt:
	$(foreach DIR,$(TELLDIRS),echo $(DIR); cd $(DIR); $(COPYRGT) *.c; $(COPYRGT) *.h; $(COPYRGT) *.cc; cd $(GMXHOME)/src;)

oclean: 
	$(foreach dir,$(DIRS),cd $(dir); mkdir obj_$(CPU) ; mv *.o obj_$(CPU) ; cd ..;)

oget:   
	$(foreach dir,$(DIRS),cd $(dir); mv obj_$(CPU)/*.o . ; cd ..;)
