 ## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

# These scripts don't depend on single/double/mpi, so don't suffix names
transform=


EXTRA_SCRIPTS = grompplog2top       make_gromos_nb.pl   make_gromos_rtp.py \
		make_gromos_bon.pl make_gromos_rtp.pl  

bin_SCRIPTS =   GMXRC GMXRC.csh GMXRC.bash GMXRC.zsh \
		completion.csh completion.bash completion.zsh

EXTRA_DIST = GMXRC.in GMXRC.csh.in GMXRC.bash.in GMXRC.zsh.in \
	     completion.csh completion.bash completion.zsh \
	     grompplog2top       make_gromos_nb.pl   make_gromos_rtp.py \
	     make_gromos_bon.pl make_gromos_rtp.pl  

GMXRC:  GMXRC.in
	cat $(srcdir)/GMXRC.in | sed 's,@bindir\@,$(bindir),' > GMXRC

GMXRC.csh:  GMXRC.csh.in
	cat $(srcdir)/GMXRC.csh.in | sed 's,@libdir\@,$(libdir),' | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@datadir\@,$(datadir),' | \
	sed 's,@mandir\@,$(mandir),'  > GMXRC.csh

GMXRC.bash:  GMXRC.bash.in
	cat $(srcdir)/GMXRC.bash.in | sed 's,@libdir\@,$(libdir),' | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@datadir\@,$(datadir),' | \
	sed 's,@mandir\@,$(mandir),'  > GMXRC.bash

GMXRC.zsh:  GMXRC.zsh.in
	cat $(srcdir)/GMXRC.zsh.in | sed 's,@libdir\@,$(libdir),' | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@datadir\@,$(datadir),' | \
	sed 's,@mandir\@,$(mandir),'  > GMXRC.zsh

CLEANFILES = GMXRC GMXRC.csh GMXRC.bash GMXRC.zsh \\\#* *~










