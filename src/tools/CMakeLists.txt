
add_library(gmxana 
            autocorr.c      expfit.c        polynomials.c   levenmar.c      
            anadih.c        pp2shift.c      dlist.c         
            eigio.c         cmat.c          
            eigensolver.c   nsc.c           
            hxprops.c       fitahx.c        
            geminate.c
            gmx_analyze.c   gmx_anaeig.c    gmx_angle.c     gmx_bond.c      
            gmx_bundle.c    gmx_chi.c       gmx_cluster.c   gmx_confrms.c   
            gmx_covar.c     gmx_current.c   
            gmx_density.c   gmx_densmap.c   gmx_dih.c       
            gmx_dielectric.c        
            gmx_kinetics.c  gmx_spatial.c   gmx_tune_pme.c
            gmx_dipoles.c   gmx_disre.c     gmx_dist.c      gmx_dyndom.c    
            gmx_enemat.c    gmx_energy.c    gmx_lie.c       gmx_filter.c    
            gmx_gyrate.c    gmx_h2order.c   gmx_hbond.c     gmx_helix.c     
            gmx_mindist.c   gmx_msd.c       gmx_morph.c     gmx_nmeig.c     
            gmx_nmens.c     gmx_order.c     gmx_principal.c 
            gmx_polystat.c  gmx_potential.c gmx_rama.c      
            gmx_rdf.c       gmx_rms.c       gmx_rmsf.c      
            gmx_rotacf.c    gmx_saltbr.c    gmx_sas.c              
            gmx_select.c    gmx_rmsdist.c   gmx_rotmat.c
            gmx_sgangle.c   gmx_sorient.c   gmx_spol.c      gmx_tcaf.c      
            gmx_traj.c      gmx_velacc.c    gmx_helixorient.c 
            gmx_clustsize.c gmx_mdmat.c     gmx_wham.c      
            correl.c        gmx_sham.c      gmx_nmtraj.c    
            gmx_trjconv.c   gmx_trjcat.c    gmx_trjorder.c  gmx_xpm2ps.c    
            gmx_editconf.c  gmx_genbox.c    gmx_genion.c    gmx_genconf.c   
            gmx_genpr.c     gmx_eneconv.c   gmx_vanhove.c   gmx_wheel.c     
            addconf.c       calcpot.c       edittop.c       gmx_bar.c
            gmx_membed.c    gmx_pme_error.c gmx_options.c   gmx_dos.c
            gmx_hydorder.c  gmx_densorder.c powerspect.c    dens_filter.c
            binsearch.c
            )


target_link_libraries(gmxana md gmx ${GSL_LIBRARIES})
set_target_properties(gmxana PROPERTIES OUTPUT_NAME "gmxana${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

# List of programs with single corresponding .c source file,
# used to create build rules automatically.
#
set(GMX_TOOLS_PROGRAMS
    do_dssp editconf eneconv genbox genconf genrestr g_nmtraj 
    make_ndx mk_angndx trjcat trjconv trjorder g_wheel 
    xpm2ps genion g_anadock make_edi g_analyze g_anaeig
    g_angle g_bond g_bundle g_chi g_cluster g_confrms g_covar
    g_current g_density g_densmap g_dih g_dielectric
    g_helixorient g_principal g_dipoles g_disre g_dist
    g_dyndom g_enemat g_energy g_lie g_filter g_gyrate
    g_h2order g_hbond g_helix g_mindist g_msd g_morph g_nmeig
    g_nmens g_order g_kinetics g_polystat g_potential g_rama g_rdf g_rms
    g_rmsf g_rotacf g_saltbr g_sas g_select g_sgangle g_sham g_sorient
    g_spol g_spatial g_tcaf g_traj g_tune_pme g_vanhove
    g_velacc g_clustsize g_mdmat g_wham g_sigeps g_bar
    g_membed g_pme_error g_rmsdist g_rotmat g_options
    g_dos    g_hydorder  g_densorder
    )

set(GMX_TOOLS_PROGRAMS_NOT_FOR_INSTALLATION
  # names of any executables that should be built but not installed can go here
    )


foreach(TOOL ${GMX_TOOLS_PROGRAMS} ${GMX_TOOLS_PROGRAMS_NOT_FOR_INSTALLATION})
    add_executable(${TOOL} ${TOOL}.c)
    target_link_libraries(${TOOL} gmxana)
    set_target_properties(${TOOL} PROPERTIES OUTPUT_NAME "${TOOL}${GMX_BINARY_SUFFIX}")
endforeach(TOOL ${GMX_TOOLS_PROGRAMS}) 


install(TARGETS gmxana DESTINATION ${LIB_INSTALL_DIR} COMPONENT runtime)
install(TARGETS ${GMX_TOOLS_PROGRAMS}
        DESTINATION ${BIN_INSTALL_DIR}
        COMPONENT runtime)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmxana.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmxana.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmxana.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgmxana${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
