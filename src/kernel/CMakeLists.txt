
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
    fflibutil.c
    xlate.c)

set(MDRUN_SOURCES 
    gctio.c    ionize.c runner.c
    do_gct.c     repl_ex.c  xutils.c
    md.c         mdrun.c    genalg.c md_openmm.c)

add_library(gmxpreprocess ${GMXPREPROCESS_SOURCES})
target_link_libraries(gmxpreprocess md)
set_target_properties(gmxpreprocess PROPERTIES OUTPUT_NAME "gmxpreprocess${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")
set_target_properties(md PROPERTIES OUTPUT_NAME "md${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")


if(GMX_OPENMM) 
    add_subdirectory(gmx_gpu_utils)
    include_directories(./gmx_gpu_utils ${OpenMM_INCLUDE_DIR})
    link_directories(${OpenMM_LIBRARY_DIR}) 
    # only define if this is a local build not a release 
    # we assume that the auto-generated version is not used && 
    # version string does not contain "-dev" => it's a release build
    if(NOT USE_VERSION_H AND NOT PROJECT_VERSION MATCHES ".*-dev.*")  
        add_definitions( -DOPENMM_PLUGIN_DIR="${OpenMM_PLUGIN_DIR}" ) 
    else()
        add_definitions( -DOPENMM_PLUGIN_DIR="" )
    endif()
    add_library(openmm_api_wrapper STATIC openmm_wrapper.cpp)
    target_link_libraries(openmm_api_wrapper gmx_gpu_utils ${OpenMM_LIBRARIES})
    set(GMX_OPENMM_LIBRARIES openmm_api_wrapper gmx_gpu_utils ${OpenMM_LIBRARIES})   
endif(GMX_OPENMM)

if(GMX_FAHCORE)
  add_library(fahcore ${MDRUN_SOURCES})
else(GMX_FAHCORE)

list(APPEND GMX_EXTRA_LIBRARIES gmxpreprocess md)

add_executable(grompp grompp.c)
target_link_libraries(grompp ${GMX_EXTRA_LIBRARIES})
set_target_properties(grompp PROPERTIES OUTPUT_NAME "grompp${GMX_BINARY_SUFFIX}")

add_executable(tpbconv tpbconv.c)
target_link_libraries(tpbconv ${GMX_EXTRA_LIBRARIES})
set_target_properties(tpbconv PROPERTIES OUTPUT_NAME "tpbconv${GMX_BINARY_SUFFIX}")

add_executable(pdb2gmx pdb2gmx.c)
target_link_libraries(pdb2gmx ${GMX_EXTRA_LIBRARIES})
set_target_properties(pdb2gmx PROPERTIES OUTPUT_NAME "pdb2gmx${GMX_BINARY_SUFFIX}")

add_executable(g_protonate g_protonate.c)
target_link_libraries(g_protonate ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_protonate PROPERTIES OUTPUT_NAME "g_protonate${GMX_BINARY_SUFFIX}")

add_executable(g_luck g_luck.c)
target_link_libraries(g_luck ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_luck PROPERTIES OUTPUT_NAME "g_luck${GMX_BINARY_SUFFIX}")

add_executable(gmxdump gmxdump.c)
target_link_libraries(gmxdump ${GMX_EXTRA_LIBRARIES})
set_target_properties(gmxdump PROPERTIES OUTPUT_NAME "gmxdump${GMX_BINARY_SUFFIX}")

add_executable(g_x2top g_x2top.c)
target_link_libraries(g_x2top ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_x2top PROPERTIES OUTPUT_NAME "g_x2top${GMX_BINARY_SUFFIX}")

add_executable(gmxcheck gmxcheck.c)
target_link_libraries(gmxcheck ${GMX_EXTRA_LIBRARIES})
set_target_properties(gmxcheck PROPERTIES OUTPUT_NAME "gmxcheck${GMX_BINARY_SUFFIX}")

add_executable(mdrun ${MDRUN_SOURCES})
target_link_libraries(mdrun ${GMX_EXTRA_LIBRARIES} ${GMX_OPENMM_LIBRARIES})
set_target_properties(mdrun PROPERTIES OUTPUT_NAME "mdrun${GMX_BINARY_SUFFIX}")

# this is to circumvent the following MSVC error: 
# warning LNK4098: defaultlib 'LIBCMT' conflicts with use of other libs
# fatal error LNK1169: one or more multiply defined symbols found
if(GMX_OPENMM AND MSVC)
    set_target_properties(mdrun PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT")
endif()


install(TARGETS 
        grompp
        tpbconv
        pdb2gmx
        g_protonate
        g_luck
        gmxdump
        g_x2top
        gmxcheck
        mdrun
	gmxpreprocess DESTINATION ${LIB_INSTALL_DIR}
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})


# if we build shared gromacs libs, when installing throught the install-mdrun target 
# these libs need to be installed as well
if(BUILD_SHARED_LIBS)
    # in MDRUN_LIBS we store the libraries MDRUN links against (NOTE: hardcoded!!!)
    set(MDRUN_LIBS gmxpreprocess md gmx)
    
    # generate install-libXXX custom target for each shared lib that mdrun links against
    foreach(_lib ${MDRUN_LIBS})
        # double-check that the type is SHARED
        get_target_property(_type ${_lib} TYPE)
        if(NOT ${_type} STREQUAL "SHARED_LIBRARY")
            message(FATAL_ERROR " Internal error: library ${_lib} is not shared so it's not supposed to be processed for installing")
        endif()

        # figure out the path and filename under which the lib will be installed
        # (libname with pre- and suffix)
        get_target_property(_lib_path ${_lib} LOCATION)
        string(REGEX REPLACE "/" ";" _lib_fname ${_lib_path})
        list(REVERSE _lib_fname)
        list(GET _lib_fname 0 _lib_fname)

        # create custom target for copying each library to the install location 
        # TODO: need to fix this to have the .so.6 form
        add_custom_target(install-${_lib}
            COMMAND ${CMAKE_COMMAND} -E copy 
                "${_lib_path}" "${LIB_INSTALL_DIR}/${_lib_fname}.${SOVERSION}"
            COMMAND ${CMAKE_COMMAND} -E create_symlink 
                "${_lib_fname}.${SOVERSION}" "${LIB_INSTALL_DIR}/${_lib_fname}"
            COMMENT "Installing library ${_lib}")
        add_dependencies(install-${_lib} ${_lib})

        # gather the custom target names in a string
        # set(_lib_install_targets "${_lib_install_targets} install-lib${_lib}")
        list(APPEND _lib_install_targets "install-${_lib}")
    endforeach(_lib)
endif(BUILD_SHARED_LIBS)

get_target_property(_mdrun_path mdrun LOCATION)
add_custom_target(install-mdrun
    COMMAND ${CMAKE_COMMAND} -E copy "${_mdrun_path}" 
        "${BIN_INSTALL_DIR}/${_mdrun_exec_name}"
    COMMENT "Installing mdrun")
add_dependencies(install-mdrun mdrun ${_lib_install_targets})

endif(GMX_FAHCORE)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmxpreprocess.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmxpreprocess.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmxpreprocess.pc DESTINATION ${LIB_INSTALL_DIR}/pkgconfig 
  RENAME "libgmxpreprocess${GMX_LIBS_SUFFIX}.pc")
