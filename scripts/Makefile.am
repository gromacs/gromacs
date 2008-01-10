 ## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

# These scripts don't depend on single/double/mpi, so don't suffix names
transform=

# Scripts generated at compile time
bin_SCRIPTS   = GMXRC GMXRC.csh GMXRC.bash      GMXRC.zsh                   

# Distributed scripts
dist_bin_SCRIPTS =   \
	completion.csh  completion.bash completion.zsh
	xplor2gmx.pl	demux.pl			


# Scripts that will not be installed
EXTRA_SCRIPTS = 		grompplog2top       	make_gromos_nb.pl  \
	make_gromos_rtp.py 	make_gromos_bon.pl 	make_gromos_rtp.pl 


# Extra files we need
EXTRA_DIST = $(EXTRA_SCRIPTS) \
        GMXRC.in GMXRC.csh.in GMXRC.bash.in GMXRC.zsh.in                

GMXRC:  GMXRC.in Makefile
	cat $(srcdir)/GMXRC.in | sed 's,@bindir\@,$(bindir),' > GMXRC

GMXRC.csh:  GMXRC.csh.in Makefile
	cat $(srcdir)/GMXRC.csh.in | sed 's,@libdir\@,$(libdir),' | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@datadir\@,$(datadir),' | \
	sed 's,@mandir\@,$(mandir),'  > GMXRC.csh

GMXRC.bash:  GMXRC.bash.in Makefile
	cat $(srcdir)/GMXRC.bash.in | sed 's,@libdir\@,$(libdir),' | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@datadir\@,$(datadir),' | \
	sed 's,@mandir\@,$(mandir),'  > GMXRC.bash

GMXRC.zsh:  GMXRC.zsh.in Makefile
	cat $(srcdir)/GMXRC.zsh.in | sed 's,@libdir\@,$(libdir),' | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@datadir\@,$(datadir),' | \
	sed 's,@mandir\@,$(mandir),'  > GMXRC.zsh

CLEANFILES = GMXRC GMXRC.csh GMXRC.bash GMXRC.zsh \\\#* *~
