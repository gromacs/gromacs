#
# RPM specification file to make gromacs packages, version 3.0
# Presently, you cannot relocate from /usr/local/gromacs.
#
# Usage:
#
# 1. Start from a gromacs distribution tarball, made
#    with "make dist". Put it in the RPM
#    source directory (usually /usr/src/redhat/SOURCES).
# 2. Edit the version and release info below (bump the
#    release every time you release a new rpm, restore it
#    to 1 for each a new version.)
# 3. Edit the files tags IF YOU MOVE OR ADD ANY FILES
#    (also if you change lib versions)
# 4. This file assumes a i686-pc-linux-gnu configuration -
#    you will have to change that for a different host,
#    since it enters in the directory names gromacs creates.
# 5. cd to /usr/src/redhat/SPECS and issue 
#    rpm -ba gromacs-3.0.spec
#
#    That's it - you should have both binary and source rpms now.
#

#
# Main package - only dynamic libs, and no header files
#
Summary: A package for molecular dynamics simulation 
Name: gromacs
Version: 3.0
Release: 1
Copyright: GPL
Group: Applications/Science
Source: http://www.gromacs.org/download/gromacs_source/gromacs-3.0.tar.gz
URL: http://www.gromacs.org
Packager: Erik Lindahl <lindahl@gromacs.org>
%description
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high 
performance means it is used also in several other field 
like polymer chemistry and solid state physics. This
version has the dynamic libs and executables; to hack new
utility programs you also need the headers and static
libs in gromacs-dev. Linux kernel 2.4 is STRONGLY
recommended on Pentium III and later processors since
GROMACS can then use assembly loops with SSE instructions.
#
# The header files and static libraries go into gromacs-devel...
#
%package devel
Summary: Header files and static libraries for GROMACS
Group: Applications/Science
Requires: gromacs = %{version}-%{release}
%description devel
This package contains header files, static libraries,
and a program example for the GROMACS molecular
dynamics software. You need it if you want to write your
own analysis programs.


%prep
%setup

%build
./configure

%install
make install
make links

%post
#
# Add our library dir to /etc/ld.so.conf if it is not already there
#
if test -z `grep /usr/local/gromacs/lib/i686-pc-linux-gnu /etc/ld.so.conf`; then
     cat >> /etc/ld.so.conf < /usr/local/gromacs/lib/i686-pc-linux-gnu
fi

# run ldconfig to update the runtime linker database with the new libraries
# (make sure /sbin is in the $PATH)
PATH="/sbin:$PATH" ldconfig

%postun
#
# Remove gromacs lib dir from /etc/ld.so.conf, since nothing else resides there
# 
grep -v /usr/local/gromacs/lib/i686-pc-linux-gnu /etc/ld.so.conf > tmpconf
mv tmpconf /etc/ld.so.conf

# after uninstall, run ldconfig to remove the libs from the linker database
PATH="/sbin:$PATH" ldconfig




%files
# binaries
/usr/local/gromacs/bin/i686-pc-linux-gnu/average    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_confrms     
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_gyrate  
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_order      
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_order
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_order
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_sorient  
/usr/local/gromacs/bin/i686-pc-linux-gnu/highway    
/usr/local/gromacs/bin/i686-pc-linux-gnu/trjconv
/usr/local/gromacs/bin/i686-pc-linux-gnu/do_dssp    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_covar       
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_h2order  
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_potential  
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_tcaf     
/usr/local/gromacs/bin/i686-pc-linux-gnu/luck       
/usr/local/gromacs/bin/i686-pc-linux-gnu/trjorder
/usr/local/gromacs/bin/i686-pc-linux-gnu/editconf   
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_density     
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_hbond    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_rama       
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_traj     
/usr/local/gromacs/bin/i686-pc-linux-gnu/make_ndx   
/usr/local/gromacs/bin/i686-pc-linux-gnu/wheel
/usr/local/gromacs/bin/i686-pc-linux-gnu/eneconv    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_dielectric  
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_helix    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_rdf        
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_velacc   
/usr/local/gromacs/bin/i686-pc-linux-gnu/mdrun      
/usr/local/gromacs/bin/i686-pc-linux-gnu/x2top
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_anaeig   
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_dih         
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_lie      
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_rms        
/usr/local/gromacs/bin/i686-pc-linux-gnu/genbox     
/usr/local/gromacs/bin/i686-pc-linux-gnu/mk_angndx  
/usr/local/gromacs/bin/i686-pc-linux-gnu/xmdrun
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_analyze  
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_dipoles     
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_mdmat    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_rmsdist    
/usr/local/gromacs/bin/i686-pc-linux-gnu/genconf    
/usr/local/gromacs/bin/i686-pc-linux-gnu/ngmx       
/usr/local/gromacs/bin/i686-pc-linux-gnu/xpm2ps
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_angle    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_disre       
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_mindist  
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_rmsf       
/usr/local/gromacs/bin/i686-pc-linux-gnu/genion     
/usr/local/gromacs/bin/i686-pc-linux-gnu/nmrun      
/usr/local/gromacs/bin/i686-pc-linux-gnu/xrama
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_bond     
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_dist        
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_morph    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_rotacf     
/usr/local/gromacs/bin/i686-pc-linux-gnu/genpr      
/usr/local/gromacs/bin/i686-pc-linux-gnu/pdb2gmx
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_bundle   
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_dyndom      
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_msd      
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_saltbr     
/usr/local/gromacs/bin/i686-pc-linux-gnu/gmxcheck   
/usr/local/gromacs/bin/i686-pc-linux-gnu/protonate
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_chi      
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_enemat      
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_nmeig    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_sas        
/usr/local/gromacs/bin/i686-pc-linux-gnu/gmxdump    
/usr/local/gromacs/bin/i686-pc-linux-gnu/tpbconv
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_cluster  
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_energy      
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_nmens    
/usr/local/gromacs/bin/i686-pc-linux-gnu/g_sgangle    
/usr/local/gromacs/bin/i686-pc-linux-gnu/grompp     
/usr/local/gromacs/bin/i686-pc-linux-gnu/trjcat
#links to /usr/local/bin
/usr/local/bin/average    
/usr/local/bin/g_confrms     
/usr/local/bin/g_gyrate  
/usr/local/bin/g_order      
/usr/local/bin/g_order
/usr/local/bin/g_order
/usr/local/bin/g_sorient  
/usr/local/bin/highway    
/usr/local/bin/trjconv
/usr/local/bin/do_dssp    
/usr/local/bin/g_covar       
/usr/local/bin/g_h2order  
/usr/local/bin/g_potential  
/usr/local/bin/g_tcaf     
/usr/local/bin/luck       
/usr/local/bin/trjorder
/usr/local/bin/editconf   
/usr/local/bin/g_density     
/usr/local/bin/g_hbond    
/usr/local/bin/g_rama       
/usr/local/bin/g_traj     
/usr/local/bin/make_ndx   
/usr/local/bin/wheel
/usr/local/bin/eneconv    
/usr/local/bin/g_dielectric  
/usr/local/bin/g_helix    
/usr/local/bin/g_rdf        
/usr/local/bin/g_velacc   
/usr/local/bin/mdrun      
/usr/local/bin/x2top
/usr/local/bin/g_anaeig   
/usr/local/bin/g_dih         
/usr/local/bin/g_lie      
/usr/local/bin/g_rms        
/usr/local/bin/genbox     
/usr/local/bin/mk_angndx  
/usr/local/bin/xmdrun
/usr/local/bin/g_analyze  
/usr/local/bin/g_dipoles     
/usr/local/bin/g_mdmat    
/usr/local/bin/g_rmsdist    
/usr/local/bin/genconf    
/usr/local/bin/ngmx       
/usr/local/bin/xpm2ps
/usr/local/bin/g_angle    
/usr/local/bin/g_disre       
/usr/local/bin/g_mindist  
/usr/local/bin/g_rmsf       
/usr/local/bin/genion     
/usr/local/bin/nmrun      
/usr/local/bin/xrama
/usr/local/bin/g_bond     
/usr/local/bin/g_dist        
/usr/local/bin/g_morph    
/usr/local/bin/g_rotacf     
/usr/local/bin/genpr      
/usr/local/bin/pdb2gmx
/usr/local/bin/g_bundle   
/usr/local/bin/g_dyndom      
/usr/local/bin/g_msd      
/usr/local/bin/g_saltbr     
/usr/local/bin/gmxcheck   
/usr/local/bin/protonate
/usr/local/bin/g_chi      
/usr/local/bin/g_enemat      
/usr/local/bin/g_nmeig    
/usr/local/bin/g_sas        
/usr/local/bin/gmxdump    
/usr/local/bin/tpbconv
/usr/local/bin/g_cluster  
/usr/local/bin/g_energy      
/usr/local/bin/g_nmens    
/usr/local/bin/g_sgangle    
/usr/local/bin/grompp     
/usr/local/bin/trjcat
# the topology library
/usr/local/gromacs/top/
/usr/local/gromacs/top/FF.dat
/usr/local/gromacs/top/ffgmx.itp
/usr/local/gromacs/top/ffgmxnb.itp
/usr/local/gromacs/top/ffgmxbon.itp
/usr/local/gromacs/top/ffgmx.atp
/usr/local/gromacs/top/ffgmx.hdb
/usr/local/gromacs/top/ffgmx.n2t
/usr/local/gromacs/top/ffgmx.rtp
/usr/local/gromacs/top/ffgmx-c.tdb
/usr/local/gromacs/top/ffgmx-n.tdb
/usr/local/gromacs/top/ffgmx2.itp
/usr/local/gromacs/top/ffgmx2nb.itp
/usr/local/gromacs/top/ffgmx2bon.itp
/usr/local/gromacs/top/ffgmx2.atp
/usr/local/gromacs/top/ffgmx2.hdb
/usr/local/gromacs/top/ffgmx2.rtp
/usr/local/gromacs/top/ffgmx2-c.tdb
/usr/local/gromacs/top/ffgmx2-n.tdb
/usr/local/gromacs/top/ffG43a1.itp
/usr/local/gromacs/top/ffG43a1nb.itp
/usr/local/gromacs/top/ffG43a1bon.itp
/usr/local/gromacs/top/ffG43a1.atp
/usr/local/gromacs/top/ffG43a1.hdb
/usr/local/gromacs/top/ffG43a1.rtp
/usr/local/gromacs/top/ffG43a1-c.tdb
/usr/local/gromacs/top/ffG43a1-n.tdb
/usr/local/gromacs/top/ffG43a2.itp
/usr/local/gromacs/top/ffG43a2nb.itp
/usr/local/gromacs/top/ffG43a2bon.itp
/usr/local/gromacs/top/ffG43a2.atp
/usr/local/gromacs/top/ffG43a2.hdb
/usr/local/gromacs/top/ffG43a2.rtp
/usr/local/gromacs/top/ffG43a2-c.tdb
/usr/local/gromacs/top/ffG43a2-n.tdb
/usr/local/gromacs/top/ffG43b1.itp
/usr/local/gromacs/top/ffG43b1nb.itp
/usr/local/gromacs/top/ffG43b1bon.itp
/usr/local/gromacs/top/ffG43b1.atp
/usr/local/gromacs/top/ffG43b1.hdb
/usr/local/gromacs/top/ffG43b1.rtp
/usr/local/gromacs/top/ffG43b1-c.tdb
/usr/local/gromacs/top/ffG43b1-n.tdb
/usr/local/gromacs/top/1mlg.itp
/usr/local/gromacs/top/2mlg.itp
/usr/local/gromacs/top/benzamide.itp
/usr/local/gromacs/top/bondadd.itp
/usr/local/gromacs/top/buck.itp
/usr/local/gromacs/top/decane.itp
/usr/local/gromacs/top/dlg.itp
/usr/local/gromacs/top/dmso.itp
/usr/local/gromacs/top/fa.itp
/usr/local/gromacs/top/ff_dum.itp
/usr/local/gromacs/top/flexspc.itp
/usr/local/gromacs/top/flexspce.itp
/usr/local/gromacs/top/flexwat-ferguson.itp
/usr/local/gromacs/top/h2p4o13.itp
/usr/local/gromacs/top/h2p8o25.itp
/usr/local/gromacs/top/h2po4.itp
/usr/local/gromacs/top/ions.itp
/usr/local/gromacs/top/methanol.itp
/usr/local/gromacs/top/spc.itp
/usr/local/gromacs/top/spce.itp
/usr/local/gromacs/top/tfe.itp
/usr/local/gromacs/top/tip3pgmx.itp
/usr/local/gromacs/top/tip4pgmx.itp
/usr/local/gromacs/top/urea.itp
/usr/local/gromacs/top/dec50.gro
/usr/local/gromacs/top/dmso.gro
/usr/local/gromacs/top/spc216.gro
/usr/local/gromacs/top/tip4p.gro
/usr/local/gromacs/top/urea+h2o.gro
/usr/local/gromacs/top/aminoacids.dat
/usr/local/gromacs/top/atommass.dat
/usr/local/gromacs/top/bromacs.dat
/usr/local/gromacs/top/ca-shift.dat
/usr/local/gromacs/top/cb-shift.dat
/usr/local/gromacs/top/co-shift.dat
/usr/local/gromacs/top/edissoc.dat
/usr/local/gromacs/top/gurgle.dat
/usr/local/gromacs/top/ha-shift.dat
/usr/local/gromacs/top/links.dat
/usr/local/gromacs/top/phbres.dat
/usr/local/gromacs/top/random.dat
/usr/local/gromacs/top/refi_aa.dat
/usr/local/gromacs/top/specbond.dat
/usr/local/gromacs/top/surface.dat
/usr/local/gromacs/top/vdwradii.dat
/usr/local/gromacs/top/xlateat.dat
/usr/local/gromacs/top/export.dlg
/usr/local/gromacs/top/bonds.dlg
/usr/local/gromacs/top/ss.map
/usr/local/gromacs/top/ps.m2p
/usr/local/gromacs/top/table6-10.xvg
/usr/local/gromacs/top/table6-11.xvg
/usr/local/gromacs/top/table6-12.xvg
/usr/local/gromacs/top/table6-8.xvg
/usr/local/gromacs/top/table6-9.xvg
# examples
/usr/local/gromacs/share/tutor/cleanit
/usr/local/gromacs/share/tutor/gmxdemo/cpeptide.pdb
/usr/local/gromacs/share/tutor/gmxdemo/demo
/usr/local/gromacs/share/tutor/gmxdemo/demo
/usr/local/gromacs/share/tutor/nmr1/conf.gro
/usr/local/gromacs/share/tutor/nmr1/grompp.mdp
/usr/local/gromacs/share/tutor/nmr1/pep.pdb
/usr/local/gromacs/share/tutor/nmr1/topol.top
/usr/local/gromacs/share/tutor/nmr2/conf.gro
/usr/local/gromacs/share/tutor/nmr2/grompp.mdp
/usr/local/gromacs/share/tutor/nmr2/pep.pdb
/usr/local/gromacs/share/tutor/nmr2/topol.top
/usr/local/gromacs/share/tutor/nmr2/genconf.gcp
/usr/local/gromacs/share/tutor/water/water.top
/usr/local/gromacs/share/tutor/water/water.mdp
/usr/local/gromacs/share/tutor/water/spc216.gro
/usr/local/gromacs/share/tutor/water/spc216.pdb
/usr/local/gromacs/share/tutor/water/oxygen.ndx
/usr/local/gromacs/share/tutor/speptide/speptide.pdb
/usr/local/gromacs/share/tutor/speptide/pr.mdp
/usr/local/gromacs/share/tutor/speptide/em.mdp
/usr/local/gromacs/share/tutor/speptide/full.mdp
# manual pages
/usr/local/gromacs/man/
/usr/local/gromacs/man/man1/
/usr/local/gromacs/man/man1/g_dih.1
/usr/local/gromacs/man/man1/g_msd.1
/usr/local/gromacs/man/man1/g_tcaf.1
/usr/local/gromacs/man/man1/nmrun.1
/usr/local/gromacs/man/man1/do_dssp.1
/usr/local/gromacs/man/man1/g_dipoles.1
/usr/local/gromacs/man/man1/g_nmeig.1
/usr/local/gromacs/man/man1/g_traj.1
/usr/local/gromacs/man/man1/pdb2gmx.1
/usr/local/gromacs/man/man1/editconf.1
/usr/local/gromacs/man/man1/g_disre.1
/usr/local/gromacs/man/man1/g_nmens.1
/usr/local/gromacs/man/man1/g_velacc.1
/usr/local/gromacs/man/man1/protonate.1
/usr/local/gromacs/man/man1/eneconv.1
/usr/local/gromacs/man/man1/g_dist.1
/usr/local/gromacs/man/man1/g_order.1
/usr/local/gromacs/man/man1/genbox.1
/usr/local/gromacs/man/man1/tpbconv.1
/usr/local/gromacs/man/man1/g_anaeig.1
/usr/local/gromacs/man/man1/g_dyndom.1
/usr/local/gromacs/man/man1/g_potential.1
/usr/local/gromacs/man/man1/genconf.1
/usr/local/gromacs/man/man1/trjcat.1
/usr/local/gromacs/man/man1/g_analyze.1
/usr/local/gromacs/man/man1/g_enemat.1
/usr/local/gromacs/man/man1/g_rama.1
/usr/local/gromacs/man/man1/genion.1
/usr/local/gromacs/man/man1/trjconv.1
/usr/local/gromacs/man/man1/g_angle.1
/usr/local/gromacs/man/man1/g_energy.1
/usr/local/gromacs/man/man1/g_rdf.1
/usr/local/gromacs/man/man1/genpr.1
/usr/local/gromacs/man/man1/trjorder.1
/usr/local/gromacs/man/man1/g_bond.1
/usr/local/gromacs/man/man1/g_gyrate.1
/usr/local/gromacs/man/man1/g_rms.1
/usr/local/gromacs/man/man1/gmxcheck.1
/usr/local/gromacs/man/man1/wheel.1
/usr/local/gromacs/man/man1/g_bundle.1
/usr/local/gromacs/man/man1/g_h2order.1
/usr/local/gromacs/man/man1/g_rmsdist.1
/usr/local/gromacs/man/man1/gmxdump.1
/usr/local/gromacs/man/man1/x2top.1
/usr/local/gromacs/man/man1/g_chi.1
/usr/local/gromacs/man/man1/g_hbond.1
/usr/local/gromacs/man/man1/g_rmsf.1
/usr/local/gromacs/man/man1/grompp.1
/usr/local/gromacs/man/man1/xpm2ps.1
/usr/local/gromacs/man/man1/g_cluster.1
/usr/local/gromacs/man/man1/g_helix.1
/usr/local/gromacs/man/man1/g_rotacf.1
/usr/local/gromacs/man/man1/highway.1
/usr/local/gromacs/man/man1/xrama.1
/usr/local/gromacs/man/man1/g_confrms.1
/usr/local/gromacs/man/man1/g_lie.1
/usr/local/gromacs/man/man1/g_saltbr.1
/usr/local/gromacs/man/man1/make_ndx.1
/usr/local/gromacs/man/man1/g_covar.1
/usr/local/gromacs/man/man1/g_mdmat.1
/usr/local/gromacs/man/man1/g_sas.1
/usr/local/gromacs/man/man1/mdrun.1
/usr/local/gromacs/man/man1/g_density.1
/usr/local/gromacs/man/man1/g_mindist.1
/usr/local/gromacs/man/man1/g_sgangle.1
/usr/local/gromacs/man/man1/mk_angndx.1
/usr/local/gromacs/man/man1/g_morph.1
/usr/local/gromacs/man/man1/g_sorient.1
/usr/local/gromacs/man/man1/ngmx.1
/usr/local/gromacs/man/man1/g_dielectric.1
# html pages
/usr/local/gromacs/html/
/usr/local/gromacs/html/gmxfaq.html
/usr/local/gromacs/html/online.html
/usr/local/gromacs/html/gif/
/usr/local/gromacs/html/gif/annealdn.gif
/usr/local/gromacs/html/gif/features.gif
/usr/local/gromacs/html/gif/flow_leftrightup.gif
/usr/local/gromacs/html/gif/flow_vrule.gif
/usr/local/gromacs/html/gif/annealup.gif
/usr/local/gromacs/html/gif/flow_down.gif
/usr/local/gromacs/html/gif/flow_leftup.gif
/usr/local/gromacs/html/gif/links.gif
/usr/local/gromacs/html/gif/articles.gif
/usr/local/gromacs/html/gif/flow_downleft.gif
/usr/local/gromacs/html/gif/flow_right+left.gif
/usr/local/gromacs/html/gif/mail.gif
/usr/local/gromacs/html/gif/bench.gif
/usr/local/gromacs/html/gif/flow_hline.gif
/usr/local/gromacs/html/gif/flow_right.gif
/usr/local/gromacs/html/gif/manual.gif
/usr/local/gromacs/html/gif/charts_down.gif
/usr/local/gromacs/html/gif/flow_left.gif
/usr/local/gromacs/html/gif/flow_rightleftdown.gif
/usr/local/gromacs/html/gif/rainbow.gif
/usr/local/gromacs/html/gif/charts_up.gif
/usr/local/gromacs/html/gif/flow_leftright.gif
/usr/local/gromacs/html/gif/flow_uprightleft.gif
/usr/local/gromacs/html/gif/software.gif
/usr/local/gromacs/html/gif/faq.gif
/usr/local/gromacs/html/gif/flow_leftrightdown.gif
/usr/local/gromacs/html/gif/flow_vline.gif
/usr/local/gromacs/html/gif/topologies.gif
/usr/local/gromacs/html/gif/plotje.gif
/usr/local/gromacs/html/gif/xvgr.gif
/usr/local/gromacs/html/online/
/usr/local/gromacs/html/online/edo.html
/usr/local/gromacs/html/online/g96.html
/usr/local/gromacs/html/online/log.html
/usr/local/gromacs/html/online/options.html
/usr/local/gromacs/html/online/tpa.html
/usr/local/gromacs/html/online/xvg.html
/usr/local/gromacs/html/online/edr.html
/usr/local/gromacs/html/online/m2p.html
/usr/local/gromacs/html/online/getting_started.html
/usr/local/gromacs/html/online/out.html
/usr/local/gromacs/html/online/tpb.html
/usr/local/gromacs/html/online/ene.html
/usr/local/gromacs/html/online/gro.html
/usr/local/gromacs/html/online/map.html
/usr/local/gromacs/html/online/tpr.html
/usr/local/gromacs/html/online/eps.html
/usr/local/gromacs/html/online/hat.html
/usr/local/gromacs/html/online/mdp.html
/usr/local/gromacs/html/online/xtc.html
/usr/local/gromacs/html/online/top.html
/usr/local/gromacs/html/online/pdb.html
/usr/local/gromacs/html/online/trj.html
/usr/local/gromacs/html/online/dat.html
/usr/local/gromacs/html/online/files.html
/usr/local/gromacs/html/online/mdp_opt.html
/usr/local/gromacs/html/online/rtp.html
/usr/local/gromacs/html/online/include_bot.html
/usr/local/gromacs/html/online/trr.html
/usr/local/gromacs/html/online/dlg.html
/usr/local/gromacs/html/online/flow.html
/usr/local/gromacs/html/online/mtx.html
/usr/local/gromacs/html/online/tex.html
/usr/local/gromacs/html/online/include_top.html
/usr/local/gromacs/html/online/xpm.html
/usr/local/gromacs/html/online/edi.html
/usr/local/gromacs/html/online/g87.html
/usr/local/gromacs/html/online/itp.html
/usr/local/gromacs/html/online/ndx.html
/usr/local/gromacs/html/online/style.css
/usr/local/gromacs/html/style.css
# dynamic libraries
/usr/local/gromacs/lib/i686-pc-linux-gnu/libgmx.so.1.0.0
/usr/local/gromacs/lib/i686-pc-linux-gnu/libmd.so.1.0.0
/usr/local/gromacs/lib/i686-pc-linux-gnu/libgmx.so.1
/usr/local/gromacs/lib/i686-pc-linux-gnu/libmd.so.1

#
# The header files and static libraries go into gromacs-dev...
#

%files dev
# include headers
/usr/local/gromacs/include/
/usr/local/gromacs/include/3dview.h
/usr/local/gromacs/include/do_md.h
/usr/local/gromacs/include/invblock.h
/usr/local/gromacs/include/nrjac.h
/usr/local/gromacs/include/rwtop.h
/usr/local/gromacs/include/tpxio.h
/usr/local/gromacs/include/assert.h
/usr/local/gromacs/include/do_nm.h
/usr/local/gromacs/include/javaio.h
/usr/local/gromacs/include/nrnb.h
/usr/local/gromacs/include/sheader.h
/usr/local/gromacs/include/transfer.h
/usr/local/gromacs/include/atomprop.h
/usr/local/gromacs/include/dummies.h
/usr/local/gromacs/include/list.h
/usr/local/gromacs/include/ns.h
/usr/local/gromacs/include/shift.h
/usr/local/gromacs/include/trnio.h
/usr/local/gromacs/include/axp_asm.h
/usr/local/gromacs/include/ebin.h
/usr/local/gromacs/include/macros.h
/usr/local/gromacs/include/nsb.h
/usr/local/gromacs/include/shift_util.h
/usr/local/gromacs/include/txtdump.h
/usr/local/gromacs/include/binio.h
/usr/local/gromacs/include/edsam.h
/usr/local/gromacs/include/magic.h
/usr/local/gromacs/include/nsgrid.h
/usr/local/gromacs/include/sim_util.h
/usr/local/gromacs/include/typedefs.h
/usr/local/gromacs/include/block_tx.h
/usr/local/gromacs/include/enxio.h
/usr/local/gromacs/include/main.h
/usr/local/gromacs/include/pbc.h
/usr/local/gromacs/include/smalloc.h
/usr/local/gromacs/include/update.h
/usr/local/gromacs/include/bondf.h
/usr/local/gromacs/include/ewald.h
/usr/local/gromacs/include/maths.h
/usr/local/gromacs/include/pdbio.h
/usr/local/gromacs/include/sortwater.h
/usr/local/gromacs/include/utils.h
/usr/local/gromacs/include/buffer.h
/usr/local/gromacs/include/ewald_util.h
/usr/local/gromacs/include/matio.h
/usr/local/gromacs/include/pdebug.h
/usr/local/gromacs/include/split.h
/usr/local/gromacs/include/vcm.h
/usr/local/gromacs/include/calcgrid.h
/usr/local/gromacs/include/fatal.h
/usr/local/gromacs/include/mdatoms.h
/usr/local/gromacs/include/physics.h
/usr/local/gromacs/include/vec.h
/usr/local/gromacs/include/calch.h
/usr/local/gromacs/include/ffscanf.h
/usr/local/gromacs/include/mdebin.h
/usr/local/gromacs/include/pme.h
/usr/local/gromacs/include/statusio.h
/usr/local/gromacs/include/viewit.h
/usr/local/gromacs/include/calcmu.h
/usr/local/gromacs/include/fftgrid.h
/usr/local/gromacs/include/mdrun.h
/usr/local/gromacs/include/pppm.h
/usr/local/gromacs/include/statutil.h
/usr/local/gromacs/include/vveclib.h
/usr/local/gromacs/include/callf77.h
/usr/local/gromacs/include/fftw_wrapper.h
/usr/local/gromacs/include/memdump.h
/usr/local/gromacs/include/princ.h
/usr/local/gromacs/include/steep.h
/usr/local/gromacs/include/wgms.h
/usr/local/gromacs/include/filenm.h
/usr/local/gromacs/include/memtab.h
/usr/local/gromacs/include/pull.h
/usr/local/gromacs/include/strdb.h
/usr/local/gromacs/include/wman.h
/usr/local/gromacs/include/comlib.h
/usr/local/gromacs/include/force.h
/usr/local/gromacs/include/memtest.h
/usr/local/gromacs/include/string2.h
/usr/local/gromacs/include/writeps.h
/usr/local/gromacs/include/complex.h
/usr/local/gromacs/include/futil.h
/usr/local/gromacs/include/metacode.h
/usr/local/gromacs/include/random.h
/usr/local/gromacs/include/struc2.h
/usr/local/gromacs/include/x86_3dnow.h
/usr/local/gromacs/include/comtest.h
/usr/local/gromacs/include/gbutil.h
/usr/local/gromacs/include/mpiio.h
/usr/local/gromacs/include/rbin.h
/usr/local/gromacs/include/superb.h
/usr/local/gromacs/include/x86_cpu.h
/usr/local/gromacs/include/tgroup.h
/usr/local/gromacs/include/general.h
/usr/local/gromacs/include/mshift.h
/usr/local/gromacs/include/rdgroup.h
/usr/local/gromacs/include/symtab.h
/usr/local/gromacs/include/x86_sse.h
/usr/local/gromacs/include/confio.h
/usr/local/gromacs/include/gmxfio.h
/usr/local/gromacs/include/mvdata.h
/usr/local/gromacs/include/rdklib.h
/usr/local/gromacs/include/sync.h
/usr/local/gromacs/include/xdrf.h
/usr/local/gromacs/include/constr.h
/usr/local/gromacs/include/grompp.h
/usr/local/gromacs/include/names.h
/usr/local/gromacs/include/readcomp.h
/usr/local/gromacs/include/synclib.h
/usr/local/gromacs/include/xtcio.h
/usr/local/gromacs/include/copyrite.h
/usr/local/gromacs/include/gstat.h
/usr/local/gromacs/include/network.h
/usr/local/gromacs/include/readinp.h
/usr/local/gromacs/include/sysstuff.h
/usr/local/gromacs/include/xvgr.h
/usr/local/gromacs/include/delay.h
/usr/local/gromacs/include/index.h
/usr/local/gromacs/include/nhash.h
/usr/local/gromacs/include/renum.h
/usr/local/gromacs/include/systest.h
/usr/local/gromacs/include/disre.h
/usr/local/gromacs/include/init.h
/usr/local/gromacs/include/nr.h
/usr/local/gromacs/include/reorder.h
/usr/local/gromacs/include/tags.h
/usr/local/gromacs/include/do_fit.h
/usr/local/gromacs/include/nrama.h
/usr/local/gromacs/include/rmpbc.h
/usr/local/gromacs/include/types/
/usr/local/gromacs/include/types/atoms.h
/usr/local/gromacs/include/types/edsams.h
/usr/local/gromacs/include/types/forcerec.h
/usr/local/gromacs/include/types/ifunc.h
/usr/local/gromacs/include/types/mdatom.h
/usr/local/gromacs/include/types/nsborder.h
/usr/local/gromacs/include/types/simple.h
/usr/local/gromacs/include/types/block.h
/usr/local/gromacs/include/types/energy.h
/usr/local/gromacs/include/types/graph.h
/usr/local/gromacs/include/types/inputrec.h
/usr/local/gromacs/include/types/nblist.h
/usr/local/gromacs/include/types/nsgrid.h
/usr/local/gromacs/include/types/symtab.h
/usr/local/gromacs/include/types/commrec.h
/usr/local/gromacs/include/types/enums.h
/usr/local/gromacs/include/types/group.h
/usr/local/gromacs/include/types/ishift.h
/usr/local/gromacs/include/types/nbslist.h
/usr/local/gromacs/include/types/parm.h
/usr/local/gromacs/include/types/topology.h
/usr/local/gromacs/include/types/drblock.h
/usr/local/gromacs/include/types/filenm.h
/usr/local/gromacs/include/types/idef.h
/usr/local/gromacs/include/types/matrix.h
/usr/local/gromacs/include/types/nrnb.h
/usr/local/gromacs/include/types/pulls.h
/usr/local/gromacs/include/types/trx.h
/usr/local/gromacs/share/template/template.c
/usr/local/gromacs/share/template/README
/usr/local/gromacs/share/template/Makefile
/usr/local/gromacs/lib/i686-pc-linux-gnu/libgmx.a
/usr/local/gromacs/lib/i686-pc-linux-gnu/libmd.a
/usr/local/gromacs/lib/i686-pc-linux-gnu/libgmx.la
/usr/local/gromacs/lib/i686-pc-linux-gnu/libmd.la
/usr/local/gromacs/lib/i686-pc-linux-gnu/libgmx.so
/usr/local/gromacs/lib/i686-pc-linux-gnu/libmd.so

