## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

# We don't need a makefile for each tutor directory just to install files.
# Collect them all here instead.

# Main TUTOR dir:
tutordir        = ${datadir}/tutor
tutor_SCRIPTS   = cleanit


# GMXDEMO:

gmxdemodir      = ${datadir}/tutor/gmxdemo
gmxdemo_DATA    = gmxdemo/cpeptide.pdb
gmxdemo_SCRIPTS = gmxdemo/demo

# NMR1
nmr1dir         = ${datadir}/tutor/nmr1
nmr1_DATA       = conf.gro grompp.mdp pep.pdb topol.top


# NMR2
nmr2dir         = ${datadir}/tutor/nmr2
nmr2_DATA       = conf.gro genconf.gcp grompp.mdp pep.pdb topol.top


# SPEPTIDE
speptidedir     = ${datadir}/tutor/speptide
speptide_DATA   = em.mdp  full.mdp  pr.mdp  speptide.pdb


# WATER
waterdir        = ${datadir}/tutor/water
water_DATA      = oxygen.ndx  spc216.gro  spc216.pdb  water.mdp  water.top




EXTRA_DIST = ${tutor_SCRIPTS} ${gmxdemo_DATA} ${gmxdemo_SCRIPTS} ${nmr1_DATA} ${nmr2_DATA} ${speptide_DATA} ${water_DATA}

CLEANFILES = *~ \\\#* 
