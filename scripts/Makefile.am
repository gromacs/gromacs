 ## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

EXTRA_SCRIPTS = grompplog2top       make_gromos_nb.pl   make_gromos_rtp.py \
		make_gromos_bon.pl make_gromos_rtp.pl  

bin_SCRIPTS = GMXRC completion.csh completion.bash completion.zsh

EXTRA_DIST = GMXRC.in completion.csh completion.bash completion.zsh \
	     grompplog2top       make_gromos_nb.pl   make_gromos_rtp.py \
	     make_gromos_bon.pl make_gromos_rtp.pl  

GMXRC:	GMXRC.in
	cat $(srcdir)/GMXRC.in | sed 's,@libdir\@,$(libdir),' | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@mandir\@,$(mandir),'  >$@

CLEANFILES = GMXRC \\\#* *~






