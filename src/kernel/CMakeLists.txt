
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
    gctio.c    ionize.c
    do_gct.c     repl_ex.c  xutils.c
    md.c         mdrun.c    genalg.c)

add_library(gmxpreprocess ${GMXPREPROCESS_SOURCES})
target_link_libraries(gmxpreprocess md)

set(GMX_EXTRA_LIBRARIES ${GMX_EXTRA_LIBRARIES} gmxpreprocess)

add_executable(grompp grompp.c)
target_link_libraries(grompp md ${GMX_EXTRA_LIBRARIES})

add_executable(tpbconv tpbconv.c)
target_link_libraries(tpbconv md ${GMX_EXTRA_LIBRARIES})

add_executable(pdb2gmx pdb2gmx.c)
target_link_libraries(pdb2gmx md ${GMX_EXTRA_LIBRARIES})

add_executable(protonate protonate.c)
target_link_libraries(protonate md ${GMX_EXTRA_LIBRARIES})

add_executable(luck luck.c)
target_link_libraries(luck md ${GMX_EXTRA_LIBRARIES})

add_executable(gmxdump gmxdump.c)
target_link_libraries(gmxdump md ${GMX_EXTRA_LIBRARIES})

add_executable(x2top x2top.c)
target_link_libraries(x2top md ${GMX_EXTRA_LIBRARIES})

add_executable(gmxcheck gmxcheck.c)
target_link_libraries(gmxcheck md gmx ${GMX_EXTRA_LIBRARIES})

add_executable(mdrun ${MDRUN_SOURCES})
target_link_libraries(mdrun md ${GMX_EXTRA_LIBRARIES})




if(GMX_FAHCORE)
  add_library(fahcore ${MDRUN_SOURCES})
#  libfahcore_la_CFLAGS = $(AM_CFLAGS)
#  libfahcore_la_LDFLAGS = -all-static $(AM_LDFLAGS)
#  libfahcore_la_LIBADD = libgmxpreprocess.la
endif(GMX_FAHCORE)

