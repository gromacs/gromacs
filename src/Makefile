#
# Gromacs Makefile $Id$
#
include		$(GMXHOME)/src/makef/Makefile.def
#
DIRS	=	gmxlib 	kernel	tools	ngmx	include
LOCAL	=	local
LDIRS	=	$(DIRS)	$(LOCAL)
TELLDIRS=	$(LDIRS)
DDDIRS	=	$(DIRS)	Makefile makef
FILES   =	NOGMX
SDIRS	=	$(foreach DIR,$(TELLDIRS),src/$(DIR))
DIST	=	$(foreach DIR,$(DDDIRS),src/$(DIR))
AUX	=	top tutor
DIST	=	$(FILES) $(foreach DIR,$(DDDIRS),src/$(DIR)) $(AUX)
BENCH	=	tests/bench
TESTS	=	tests/mdrun
GETIT	=	$(SDIRS)  src/makef
TDIR	=	$(GMXHOME)
TARCMD  =	/usr/local/freeware/bin/tar cvf 
TARFLG  =	--exclude=CVS
COPYRGT	=	$(GMXHOME)/src/local/copyrgt
PROGRAMS=	(cd bin/s10; ls; cd ../..)

all:	dirs
	$(foreach DIR,$(DIRS),cd $(DIR); $(MAKE); cd ..;)

it:	
	$(foreach DIR,$(LDIRS),cd $(DIR); $(MAKE); cd ..;)

get:	clean
	cd $(GMXHOME); $(foreach dir,$(GETIT),cd $(dir); sccs get SCCS; cd $(GMXHOME);)


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

tell:	info
info:
	$(foreach DIR,$(TELLDIRS),echo $(DIR); cd $(DIR); sccs info; cd ..;)

clean:
	$(foreach dir,$(DIRS),cd $(dir); $(MAKE) clean; cd ..;)

copyrgt:
	$(foreach DIR,$(TELLDIRS),echo $(DIR); cd $(DIR); sccs edit SCCS; $(COPYRGT) *.c; $(COPYRGT) *.h; $(COPYRGT) *.cc; sccs delget -y'Edited by copyrgt' SCCS ; cd $(GMXHOME)/src;)

oclean: 
	$(foreach dir,$(DIRS),cd $(dir); mkdir $(CPU).o ; mv *.o $(CPU).o ; cd ..;)

oget:   
	$(foreach dir,$(DIRS),cd $(dir); mv $(CPU).o/*.o . ; cd ..;)
