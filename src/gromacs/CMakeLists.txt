set(LIBGROMACS_SOURCES)

add_subdirectory(legacyheaders)
add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(gmxpreprocess)
add_subdirectory(analysisdata)
add_subdirectory(commandline)
add_subdirectory(linearalgebra)
add_subdirectory(onlinehelp)
add_subdirectory(options)
add_subdirectory(selection)
add_subdirectory(trajectoryanalysis)
add_subdirectory(utility)

file(GLOB LIBGROMACS_HEADERS *.h)
install(FILES ${LIBGROMACS_HEADERS} DESTINATION ${INCL_INSTALL_DIR}/gromacs
        COMPONENT development)

list(APPEND LIBGROMACS_SOURCES ${GMXLIB_SOURCES} ${MDLIB_SOURCES})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/version.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.h
    DESTINATION ${INCL_INSTALL_DIR}/gromacs
    COMPONENT development)

# apply gcc 4.4.x bug workaround
if(GMX_USE_GCC44_BUG_WORKAROUND)
   include(gmxGCC44O3BugWorkaround)
   gmx_apply_gcc44_bug_workaround("gmxlib/bondfree.c")
   gmx_apply_gcc44_bug_workaround("mdlib/force.c")
   gmx_apply_gcc44_bug_workaround("mdlib/constr.c")
endif()

add_library(libgromacs ${LIBGROMACS_SOURCES})
if (GMX_GIT_VERSION_INFO)
    add_dependencies(libgromacs gmx_version)
endif ()

target_link_libraries(libgromacs ${GMX_GPU_LIBRARIES}
                      ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES}
                      ${THREAD_LIB} ${OpenMP_SHARED_LINKER_FLAGS})
set_target_properties(libgromacs PROPERTIES
                      OUTPUT_NAME "gromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${SOVERSION}
                      INSTALL_NAME_DIR "${LIB_INSTALL_DIR}"
                      COMPILE_FLAGS "${OpenMP_C_FLAGS}")

install(TARGETS libgromacs DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgromacs.pc.cmakein
               ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgromacs${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
