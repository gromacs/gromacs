include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

set(MDRUN_SOURCES
    do_gct.c    gctio.c         genalg.c    ionize.c
    md.c        md_openmm.c     mdrun.c     membed.c
    pme_switch.c  repl_ex.c   runner.c        xutils.c)

if(GMX_OPENMM) 
    include_directories(./gmx_gpu_utils ${OpenMM_INCLUDE_DIR})
    link_directories(${OpenMM_LIBRARY_DIR}) 
    # with this define no evn.var. is needed with OPENMM_PLUGIN_DIR
    # if the same OpenMM installation is used for running and building 
    add_definitions( -DOPENMM_PLUGIN_DIR="${OpenMM_PLUGIN_DIR}" ) 
    file(TO_CMAKE_PATH ${OpenMM_PLUGIN_DIR} _path)
    add_library(openmm_api_wrapper STATIC openmm_wrapper.cpp)
    target_link_libraries(openmm_api_wrapper gmx_gpu_utils ${OpenMM_LIBRARIES})
    set(GMX_OPENMM_LIBRARIES openmm_api_wrapper gmx_gpu_utils ${OpenMM_LIBRARIES})   
endif(GMX_OPENMM)

if(GMX_FAHCORE)
    add_library(fahcore ${MDRUN_SOURCES})
else(GMX_FAHCORE)
    add_executable(mdrun ${MDRUN_SOURCES})
    gmx_add_man_page(mdrun)
    target_link_libraries(mdrun ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_OPENMM_LIBRARIES}
        ${OpenMP_LINKER_FLAGS})
    set_target_properties(mdrun PROPERTIES OUTPUT_NAME "mdrun${GMX_BINARY_SUFFIX}"
        COMPILE_FLAGS "${OpenMP_C_FLAGS}")
    install(TARGETS mdrun DESTINATION ${BIN_INSTALL_DIR} COMPONENT mdrun)

    if(GMX_OPENMM AND MSVC)
        set_target_properties(mdrun PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT")
    endif()

    # Create the custom install-mdrun target
    if (BUILD_SHARED_LIBS)
        # If shared libraries are used, we need to install the libraries in
        # addition to the mdrun binary.
        add_custom_target(install-mdrun
            COMMAND ${CMAKE_COMMAND} -DCOMPONENT=libraries
                    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
            COMMAND ${CMAKE_COMMAND} -DCOMPONENT=mdrun
                    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
            COMMENT "Installing mdrun")
    else (BUILD_SHARED_LIBS)
        add_custom_target(install-mdrun
            COMMAND ${CMAKE_COMMAND} -DCOMPONENT=mdrun
                    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
            COMMENT "Installing mdrun")
    endif (BUILD_SHARED_LIBS)
    add_dependencies(install-mdrun mdrun)
endif(GMX_FAHCORE)
