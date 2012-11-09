
set(GMXPREPROCESS_SOURCES 
    add_par.c       
    calc_verletbuf.c
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
    readadress.c      
    readrot.c
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
    fflibutil.c
    xlate.c)

set(MDRUN_SOURCES 
    gctio.c    ionize.c runner.c
    do_gct.c     repl_ex.c  xutils.c pme_loadbal.c
    md.c         mdrun.c    genalg.c membed.c
    md_openmm.c)

add_library(gmxpreprocess ${GMXPREPROCESS_SOURCES})
target_link_libraries(gmxpreprocess md)
set_target_properties(gmxpreprocess PROPERTIES OUTPUT_NAME "gmxpreprocess${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}"
    COMPILE_FLAGS "${OpenMP_C_FLAGS}")


if(GMX_GPU OR GMX_OPENMM) 
    include_directories(${CMAKE_SOURCE_DIR}/src/gmxlib/gpu_utils)
endif()

if(GMX_OPENMM) 
    include_directories(${OpenMM_INCLUDE_DIR})
    link_directories(${OpenMM_LIBRARY_DIR}) 
    # with this define no evn.var. is needed with OPENMM_PLUGIN_DIR
    # if the same OpenMM installation is used for running and building 
    add_definitions( -DOPENMM_PLUGIN_DIR="${OpenMM_PLUGIN_DIR}" ) 
    file(TO_CMAKE_PATH ${OpenMM_PLUGIN_DIR} _path)
    add_library(openmm_api_wrapper STATIC openmm_wrapper.cpp)
    target_link_libraries(openmm_api_wrapper ${OpenMM_LIBRARIES})
    set(GMX_OPENMM_LIBRARIES openmm_api_wrapper ${OpenMM_LIBRARIES})   
endif(GMX_OPENMM)

if(GMX_FAHCORE)
  add_library(fahcore ${MDRUN_SOURCES})
else(GMX_FAHCORE)

list(APPEND GMX_EXTRA_LIBRARIES gmxpreprocess md ${OpenMP_LINKER_FLAGS})

set(GMX_KERNEL_PROGRAMS
    grompp tpbconv pdb2gmx g_protonate g_luck gmxdump g_x2top gmxcheck)

foreach(PROGRAM ${GMX_KERNEL_PROGRAMS})
    add_executable(${PROGRAM} ${PROGRAM}.c)
    if (NOT ${PROGRAM} STREQUAL "g_luck")
        gmx_add_man_page(${PROGRAM})
    endif()
    target_link_libraries(${PROGRAM} ${GMX_EXTRA_LIBRARIES})
    set_target_properties(${PROGRAM} PROPERTIES OUTPUT_NAME "${PROGRAM}${GMX_BINARY_SUFFIX}")
endforeach()

add_executable(mdrun ${MDRUN_SOURCES})
gmx_add_man_page(mdrun)
target_link_libraries(mdrun ${GMX_EXTRA_LIBRARIES} ${GMX_OPENMM_LIBRARIES})
set_target_properties(mdrun PROPERTIES OUTPUT_NAME "mdrun${GMX_BINARY_SUFFIX}" COMPILE_FLAGS "${OpenMP_C_FLAGS}")

# this is to circumvent the following MSVC error: 
# warning LNK4098: defaultlib 'LIBCMT' conflicts with use of other libs
# fatal error LNK1169: one or more multiply defined symbols found
if(GMX_OPENMM AND MSVC)
    set_target_properties(mdrun PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT")
endif()

# Construct component groups for installation; note that a component may
# belong to only one group
foreach(PROGRAM ${GMX_KERNEL_PROGRAMS})
    set(CPACK_COMPONENT_${PROGRAM}_GROUP tools)
endforeach()
set(CPACK_COMPONENT_MDRUN_GROUP mdrun)

foreach(PROGRAM ${GMX_KERNEL_PROGRAMS} mdrun)
    # Manage CPack component dependencies
    set(CPACK_COMPONENT_${PROGRAM}_DEPENDS libraries libraries-gmxpreprocess)

    # Create custom install-xxxx target
    if (BUILD_SHARED_LIBS)
    # If shared libraries are used, we need to install the libraries in
    # addition to the mdrun binary.
       add_custom_target(install-${PROGRAM}
           COMMAND ${CMAKE_COMMAND} -DCOMPONENT=libraries
                   -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
           COMMAND ${CMAKE_COMMAND} -DCOMPONENT=libraries-gmxpreprocess
                   -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
           COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${PROGRAM}
                   -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
           COMMENT "Installing ${PROGRAM}")
    else()
       add_custom_target(install-${PROGRAM}
           COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${PROGRAM}
                   -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
           COMMENT "Installing ${PROGRAM}")
    endif()
    add_dependencies(install-${PROGRAM} ${PROGRAM})

    # Finally, trigger installation
    install(
      TARGETS ${PROGRAM}
      COMPONENT ${PROGRAM}
      DESTINATION ${BIN_INSTALL_DIR}
      )
endforeach()

install(TARGETS gmxpreprocess DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries-gmxpreprocess)

if (INSTALL_CUDART_LIB) #can be set manual by user
    if (GMX_OPENMM OR GMX_GPU)
        foreach(CUDA_LIB ${CUDA_LIBRARIES})
            string(REGEX MATCH "cudart" IS_CUDART ${CUDA_LIB})
            if(IS_CUDART) #libcuda should not be installed
                #install also name-links (linker uses those)
                file(GLOB CUDA_LIBS ${CUDA_LIB}*)
                install(FILES ${CUDA_LIBS} DESTINATION
                    ${LIB_INSTALL_DIR} COMPONENT libraries)
            endif()
        endforeach()
    else()
        message(WARNING "INSTALL_CUDART_LIB only makes sense with GMX_OPENMM or GMX_GPU")
    endif()
endif ()
endif(GMX_FAHCORE)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmxpreprocess.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmxpreprocess.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmxpreprocess.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgmxpreprocess${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
