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
	completion.csh  completion.bash completion.zsh			\
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


# Scripts that will not be installed
EXTRA_SCRIPTS = 		grompplog2top       	make_gromos_nb.pl  \
	make_gromos_rtp.py 	make_gromos_bon.pl 	make_gromos_rtp.pl \
	g_anaeig.bat            g_dyndom.bat            g_potential.bat	\
	g_analyze.bat           g_enemat.bat            g_rama.bat	\
	g_angle.bat             g_energy.bat            g_rdf.bat	\
	g_bond.bat              g_filter.bat            g_rms.bat	\
	g_bundle.bat            g_gyrate.bat            g_rmsdist.bat	\
	g_chi.bat               g_h2order.bat           g_rmsf.bat	\
	g_cluster.bat           g_hbond.bat             g_rotacf.bat	\
	g_clustsize.bat         g_helix.bat             g_saltbr.bat	\
	g_confrms.bat           g_lie.bat               g_sas.bat	\
	g_covar.bat             g_mdmat.bat             g_sgangle.bat	\
	g_density.bat           g_mindist.bat           g_sorient.bat	\
	g_dielectric.bat        g_morph.bat             g_tcaf.bat	\
	g_dih.bat               g_msd.bat               g_traj.bat	\
	g_dipoles.bat           g_nmeig.bat             g_velacc.bat	\
	g_disre.bat             g_nmens.bat             g_wham.bat	\
	g_dist.bat              g_order.bat


# Extra files we need
EXTRA_DIST = $(EXTRA_SCRIPTS) \
        GMXRC.in GMXRC.csh.in GMXRC.bash.in GMXRC.zsh.in                

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










