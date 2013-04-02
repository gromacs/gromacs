include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

set(exename mdrun)
configure_file(../main.cpp main.cpp)

set(MDRUN_SOURCES
    do_gct.c      gctio.c       genalg.c    ionize.c
    md.c          mdrun.c     membed.c
    pme_loadbal.c repl_ex.c     runner.c    xutils.c
    ${CMAKE_CURRENT_BINARY_DIR}/main.cpp)

if(GMX_OPENMM)
    # Even though the OpenMM build has "moved to contrib", many things
    # have be be done from within the scope of the CMakeLists.txt that
    # builds its mdrun, and that is here
    include(../contrib/BuildMdrunOpenMM)
endif(GMX_OPENMM)

if(GMX_FAHCORE)
    add_library(fahcore ${MDRUN_SOURCES})
else(GMX_FAHCORE)
    add_executable(mdrun ${MDRUN_SOURCES})
    gmx_add_man_page(mdrun)
    target_link_libraries(mdrun ${GMX_EXTRA_LIBRARIES} libgromacs
        ${GMX_EXE_LINKER_FLAGS})
    set_target_properties(mdrun PROPERTIES OUTPUT_NAME "mdrun${GMX_BINARY_SUFFIX}"
        COMPILE_FLAGS "${OpenMP_C_FLAGS}")
    install(TARGETS mdrun DESTINATION ${BIN_INSTALL_DIR} COMPONENT mdrun)

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

if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)
