## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

templatedir = ${datadir}/template

# we don't build this, it's meant for the user to play with.
# Since a Makefile is used for the installation, we call the
# template makefile Template.mak, and rename it to Makefile when installing.

template_DATA = template.c README

EXTRA_DIST = Template.mak.in ${template_DATA}

install-data-hook:
	$(mkinstalldirs) $(DESTDIR)$(templatedir)
	echo "$(INSTALL_DATA) Template.mak $(DESTDIR)$(templatedir)/Makefile"; \
	$(INSTALL_DATA) Template.mak $(DESTDIR)$(templatedir)/Makefile; 


CLEANFILES = *~ \\\#* *.ene *.hat *.tr? *.xtc *.edr *.ene *.xvg core *.log mdout.mdp *.tp? mon.out confout.gro deshuf.ndx output.* *.gz *.ps *.dat *.g87 *~ \\\#*





