
set(GMXPREPROCESS_SOURCES 
    add_par.c       
    compute_io.c    
    convparm.c      
    gen_ad.c        
    gen_vsite.c     
    genhydro.c   
    gpp_atomtype.c  
    gpp_bond_atomtype.c     
    h_db.c          
    hackblock.c             
    hizzie.c        
    nm2type.c
    pdb2top.c       
    pgutil.c        
    readir.c        
    readpull.c      
    resall.c        
    sorting.c       
    specbond.c      
    ter_db.c        
    tomorse.c       
    topdirs.c       
    topexcl.c       
    topio.c         
    toppush.c       
    topshake.c      
    toputil.c       
    tpbcmp.c        
    vsite_parm.c    
    xlate.c)

set(MDRUN_SOURCES 
    gctio.c    ionize.c runner.c
    do_gct.c     repl_ex.c  xutils.c
    md.c         mdrun.c    genalg.c)

add_library(gmxpreprocess ${GMXPREPROCESS_SOURCES})
target_link_libraries(gmxpreprocess md)

if(GMX_FAHCORE)
  add_library(fahcore ${MDRUN_SOURCES})
else(GMX_FAHCORE)

list(APPEND GMX_EXTRA_LIBRARIES gmxpreprocess md)

add_executable(grompp grompp.c)
target_link_libraries(grompp ${GMX_EXTRA_LIBRARIES})

add_executable(tpbconv tpbconv.c)
target_link_libraries(tpbconv ${GMX_EXTRA_LIBRARIES})

add_executable(pdb2gmx pdb2gmx.c)
target_link_libraries(pdb2gmx ${GMX_EXTRA_LIBRARIES})

add_executable(protonate protonate.c)
target_link_libraries(protonate ${GMX_EXTRA_LIBRARIES})

add_executable(luck luck.c)
target_link_libraries(luck ${GMX_EXTRA_LIBRARIES})

add_executable(gmxdump gmxdump.c)
target_link_libraries(gmxdump ${GMX_EXTRA_LIBRARIES})

add_executable(x2top x2top.c)
target_link_libraries(x2top ${GMX_EXTRA_LIBRARIES})

add_executable(gmxcheck gmxcheck.c)
target_link_libraries(gmxcheck ${GMX_EXTRA_LIBRARIES})

add_executable(mdrun ${MDRUN_SOURCES})
target_link_libraries(mdrun ${GMX_EXTRA_LIBRARIES})

install(TARGETS 
        grompp
        tpbconv
        pdb2gmx
        protonate
        luck
        gmxdump
        x2top
        gmxcheck
        mdrun
	gmxpreprocess DESTINATION ${LIB_INSTALL_DIR}
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endif(GMX_FAHCORE)



