 ## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

# These scripts don't depend on single/double/mpi, so don't suffix names
transform=


EXTRA_SCRIPTS = grompplog2top       make_gromos_nb.pl   make_gromos_rtp.py \
		make_gromos_bon.pl make_gromos_rtp.pl  

dist_bin_SCRIPTS =   \
	GMXRC GMXRC.csh GMXRC.bash 	GMXRC.zsh 			\
	completion.csh 	completion.bash completion.zsh			\
	g_energy        g_wham          g_velacc        g_traj          \
	g_tcaf          g_sorient       g_sgangle       g_sas           \
	g_saltbr        g_rotacf        g_rmsf          g_rmsdist       \
	g_rms           g_rdf           g_rama          g_potential     \
	g_order         g_nmens         g_nmeig         g_msd           \
	g_morph         g_mindist       g_mdmat         g_lie           \
	g_helix         g_hbond         g_h2order       g_gyrate        \
	g_filter        g_enemat        g_dyndom        g_dist          \
	g_disre         g_dipoles       g_dih           g_dielectric    \
	g_density       g_covar         g_confrms       g_clustsize     \
	g_cluster       g_chi           g_bundle        g_bond          \
	g_angle         g_analyze       g_anaeig        


EXTRA_DIST = \
	GMXRC.in GMXRC.csh.in GMXRC.bash.in GMXRC.zsh.in \
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










