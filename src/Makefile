#
# Gromacs Makefile
#
DIRS	=	gmxlib 	kernel	tools	ngmx
LOCAL	=	local
LDIRS	=	$(DIRS)	$(LOCAL)
IDIRS	=	include	include/types
TELLDIRS=	$(LDIRS)	$(IDIRS)
DDDIRS	=	$(DIRS)	$(IDIRS) Makefile makef
SDIRS	=	$(foreach DIR,$(TELLDIRS),src/$(DIR))
DIST	=	$(foreach DIR,$(DDDIRS),src/$(DIR))
AUX	=	top tutor man/man1 man/sample
BENCH	=	tests/bench
TESTS	=	tests/mdrun
GETIT	=	$(SDIRS)  src/makef
TDIR	=	$(GMXHOME)
COPYRGT	=	$(GMXHOME)/src/local/copyrgt
PROGRAMS=	(cd bin/s10; ls; cd ../..)

all:	
	$(foreach DIR,$(DIRS),cd $(DIR); $(MAKE); cd ..;)

it:	
	$(foreach DIR,$(LDIRS),cd $(DIR); $(MAKE); cd ..;)

get:	clean
	cd $(GMXHOME); $(foreach dir,$(GETIT),cd $(dir); sccs get SCCS; cd $(GMXHOME);)


cdist:	get
	cd $(GMXHOME); tar cvFFf $(TDIR)/gmxcore.tar $(DIST) $(AUX)
	cd $(GMXHOME); tar cvFFf $(TDIR)/local.tar src/$(LOCAL)

bdist:
	cd $(GMXHOME)/$(BENCH); cleanit
	cd $(GMXHOME); tar cvFFf $(TDIR)/gmxbench.tar $(BENCH)

tdist:
	cd $(GMXHOME)/$(TESTS); cleanit
	cd $(GMXHOME); tar cvFFf $(TDIR)/gmxtests.tar $(TESTS)

dist:	cdist	bdist	tdist

man:
	./mknroff

tell:	info
info:
	$(foreach DIR,$(TELLDIRS),echo $(DIR); cd $(DIR); sccs info; cd ..;)

clean:
	$(foreach dir,$(DIRS),cd $(dir); $(MAKE) clean; cd ..;)

copyrgt:
	$(foreach DIR,$(TELLDIRS),echo $(DIR); cd $(DIR); sccs edit SCCS; $(COPYRGT) *.c; $(COPYRGT) *.h; $(COPYRGT) *.cc; sccs delget -y'Edited by copyrgt' SCCS ; cd $(GMXHOME)/src;)

