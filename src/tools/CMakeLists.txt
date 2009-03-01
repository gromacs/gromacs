
add_library(gmxana 
            autocorr.c      expfit.c        polynomials.c   levenmar.c      
            anadih.c        pp2shift.c      dlist.c         
            eigio.c         cmat.c          
            eigensolver.c   nsc.c           
            hxprops.c       fitahx.c        
            gmx_analyze.c   gmx_anaeig.c    gmx_angle.c     gmx_bond.c      
            gmx_bundle.c    gmx_chi.c       gmx_cluster.c   gmx_confrms.c   
            gmx_covar.c     gmx_current.c   
            gmx_density.c   gmx_densmap.c   gmx_dih.c       
            gmx_dielectric.c        
            gmx_kinetics.c  gmx_spatial.c
            gmx_dipoles.c   gmx_disre.c     gmx_dist.c      gmx_dyndom.c    
            gmx_enemat.c    gmx_energy.c    gmx_lie.c       gmx_filter.c    
            gmx_gyrate.c    gmx_h2order.c   gmx_hbond.c     gmx_helix.c     
            gmx_mindist.c   gmx_msd.c       gmx_morph.c     gmx_nmeig.c     
            gmx_nmens.c     gmx_order.c     gmx_principal.c 
            gmx_polystat.c  gmx_potential.c gmx_rama.c      
            gmx_rdf.c       gmx_rms.c       gmx_rmsdist.c   gmx_rmsf.c      
            gmx_rotacf.c    gmx_saltbr.c    gmx_sas.c       gmx_sdf.c       
            gmx_sgangle.c   gmx_sorient.c   gmx_spol.c      gmx_tcaf.c      
            gmx_traj.c      gmx_velacc.c    gmx_helixorient.c 
            gmx_clustsize.c gmx_mdmat.c     gmx_wham.c      
            correl.c        gmx_sham.c      gmx_nmtraj.c    
            gmx_trjconv.c   gmx_trjcat.c    gmx_trjorder.c  gmx_xpm2ps.c    
            gmx_editconf.c  gmx_genbox.c    gmx_genion.c    gmx_genconf.c   
            gmx_genpr.c     gmx_eneconv.c   gmx_vanhove.c   gmx_wheel.c     
            addconf.c       calcpot.c       edittop.c)

target_link_libraries(gmxana md ${GMX_EXTRA_LIBRARIES})

add_executable(do_dssp do_dssp.c)
target_link_libraries(do_dssp gmxana)

add_executable(editconf editconf.c)
target_link_libraries(editconf gmxana)

add_executable(eneconv eneconv.c)
target_link_libraries(eneconv gmxana)

add_executable(genbox genbox.c)
target_link_libraries(genbox gmxana)

add_executable(genconf genconf.c)
target_link_libraries(genconf gmxana)

add_executable(genrestr genrestr.c)
target_link_libraries(genrestr gmxana)

add_executable(g_nmtraj g_nmtraj.c)
target_link_libraries(g_nmtraj gmxana)

add_executable(make_ndx make_ndx.c)
target_link_libraries(make_ndx gmxana)

add_executable(mk_angndx mk_angndx.c)
target_link_libraries(mk_angndx gmxana)

add_executable(trjcat trjcat.c)
target_link_libraries(trjcat gmxana)

add_executable(trjconv trjconv.c)
target_link_libraries(trjconv gmxana)

add_executable(trjorder trjorder.c)
target_link_libraries(trjorder gmxana)

add_executable(wheel wheel.c)
target_link_libraries(wheel gmxana)

add_executable(xpm2ps xpm2ps.c)
target_link_libraries(xpm2ps gmxana)

add_executable(genion genion.c)
target_link_libraries(genion gmxana)

add_executable(anadock anadock.c)
target_link_libraries(anadock gmxana)

add_executable(make_edi make_edi.c)
target_link_libraries(make_edi gmxana)

add_executable(g_analyze g_analyze.c)
target_link_libraries(g_analyze gmxana)

add_executable(g_anaeig g_anaeig.c)
target_link_libraries(g_anaeig gmxana)

add_executable(g_angle g_angle.c)
target_link_libraries(g_angle gmxana)

add_executable(g_bond g_bond.c)
target_link_libraries(g_bond gmxana)

add_executable(g_bundle g_bundle.c)
target_link_libraries(g_bundle gmxana)

add_executable(g_chi g_chi.c)
target_link_libraries(g_chi gmxana)

add_executable(g_cluster g_cluster.c)
target_link_libraries(g_cluster gmxana)

add_executable(g_confrms g_confrms.c)
target_link_libraries(g_confrms gmxana)

add_executable(g_covar g_covar.c)
target_link_libraries(g_covar gmxana)

add_executable(g_current g_current.c)
target_link_libraries(g_current gmxana)

add_executable(g_density g_density.c)
target_link_libraries(g_density gmxana)

add_executable(g_densmap g_densmap.c)
target_link_libraries(g_densmap gmxana)

add_executable(g_dih g_dih.c)
target_link_libraries(g_dih gmxana)

add_executable(g_dielectric g_dielectric.c)
target_link_libraries(g_dielectric gmxana)

add_executable(g_helixorient g_helixorient.c)
target_link_libraries(g_helixorient gmxana)

add_executable(g_principal g_principal.c)
target_link_libraries(g_principal gmxana)

add_executable(g_dipoles g_dipoles.c)
target_link_libraries(g_dipoles gmxana)

add_executable(g_disre g_disre.c)
target_link_libraries(g_disre gmxana)

add_executable(g_dist g_dist.c)
target_link_libraries(g_dist gmxana)

add_executable(g_dyndom g_dyndom.c)
target_link_libraries(g_dyndom gmxana)

add_executable(g_enemat g_enemat.c)
target_link_libraries(g_enemat gmxana)

add_executable(g_energy g_energy.c)
target_link_libraries(g_energy gmxana)

add_executable(g_lie g_lie.c)
target_link_libraries(g_lie gmxana)

add_executable(g_filter g_filter.c)
target_link_libraries(g_filter gmxana)

add_executable(g_gyrate g_gyrate.c)
target_link_libraries(g_gyrate gmxana)

add_executable(g_h2order g_h2order.c)
target_link_libraries(g_h2order gmxana)

add_executable(g_hbond g_hbond.c)
target_link_libraries(g_hbond gmxana)

add_executable(g_helix g_helix.c)
target_link_libraries(g_helix gmxana)

add_executable(g_mindist g_mindist.c)
target_link_libraries(g_mindist gmxana)

add_executable(g_msd g_msd.c)
target_link_libraries(g_msd gmxana)

add_executable(g_morph g_morph.c)
target_link_libraries(g_morph gmxana)

add_executable(g_nmeig g_nmeig.c)
target_link_libraries(g_nmeig gmxana)

add_executable(g_nmens g_nmens.c)
target_link_libraries(g_nmens gmxana)

add_executable(g_order g_order.c)
target_link_libraries(g_order gmxana)

add_executable(g_polystat g_polystat.c)
target_link_libraries(g_polystat gmxana)

add_executable(g_potential g_potential.c)
target_link_libraries(g_potential gmxana)

add_executable(g_rama g_rama.c)
target_link_libraries(g_rama gmxana)

add_executable(g_rdf g_rdf.c)
target_link_libraries(g_rdf gmxana)

add_executable(g_rms g_rms.c)
target_link_libraries(g_rms gmxana)

add_executable(g_rmsdist g_rmsdist.c)
target_link_libraries(g_rmsdist gmxana)

add_executable(g_rmsf g_rmsf.c)
target_link_libraries(g_rmsf gmxana)

add_executable(g_rotacf g_rotacf.c)
target_link_libraries(g_rotacf gmxana)

add_executable(g_saltbr g_saltbr.c)
target_link_libraries(g_saltbr gmxana)

add_executable(g_sas g_sas.c)
target_link_libraries(g_sas gmxana)

add_executable(g_sgangle g_sgangle.c)
target_link_libraries(g_sgangle gmxana)

add_executable(g_sham g_sham.c)
target_link_libraries(g_sham gmxana)

add_executable(g_sorient g_sorient.c)
target_link_libraries(g_sorient gmxana)

add_executable(g_spol g_spol.c)
target_link_libraries(g_spol gmxana)

add_executable(g_sdf g_sdf.c)
target_link_libraries(g_sdf gmxana)

add_executable(g_spatial g_spatial.c)
target_link_libraries(g_spatial gmxana)

add_executable(g_tcaf g_tcaf.c)
target_link_libraries(g_tcaf gmxana)

add_executable(g_traj g_traj.c)
target_link_libraries(g_traj gmxana)

add_executable(g_tune_pme g_tune_pme.c gmx_tune_pme.c)
target_link_libraries(g_tune_pme gmxana)

add_executable(g_vanhove g_vanhove.c)
target_link_libraries(g_vanhove gmxana)

add_executable(g_velacc g_velacc.c)
target_link_libraries(g_velacc gmxana)

add_executable(g_clustsize g_clustsize.c)
target_link_libraries(g_clustsize gmxana)

add_executable(g_mdmat g_mdmat.c)
target_link_libraries(g_mdmat gmxana)

add_executable(g_wham g_wham.c)
target_link_libraries(g_wham gmxana)

add_executable(g_kinetics g_kinetics.c)
target_link_libraries(g_kinetics gmxana)

add_executable(sigeps sigeps.c)
target_link_libraries(sigeps gmxana)


install(TARGETS
        do_dssp         editconf        eneconv         
        genbox          genconf         genrestr        g_nmtraj        
        make_ndx        mk_angndx       trjcat          trjconv         
        trjorder        wheel           xpm2ps          genion          
        anadock         make_edi        
        g_analyze       g_anaeig        g_angle         g_bond          
        g_bundle        g_chi           g_cluster       g_confrms       
        g_covar         g_current       
        g_density       g_densmap       g_dih           
        g_dielectric    g_helixorient   g_principal 
        g_dipoles       g_disre         g_dist          g_dyndom        
        g_enemat        g_energy        g_lie           g_filter        
        g_gyrate        g_h2order       g_hbond         g_helix         
        g_mindist       g_msd           g_morph         g_nmeig         
        g_nmens         g_order         
        g_polystat      g_potential     g_rama          
        g_rdf           g_rms           g_rmsdist       g_rmsf          
        g_rotacf        g_saltbr        g_sas           g_sgangle       
        g_sham          g_sorient       g_spol          
        g_sdf           g_spatial       
        g_tcaf          g_traj          g_tune_pme   
        g_vanhove       g_velacc        
        g_clustsize     g_mdmat         g_wham          g_kinetics      
        sigeps
      	RUNTIME	       DESTINATION ${BIN_INSTALL_DIR})
