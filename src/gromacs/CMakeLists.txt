add_subdirectory(legacyheaders)
add_subdirectory(mdlib)

set(LIBGROMACS_SOURCES ${MDLIB_SOURCES})

add_library(libgromacs ${LIBGROMACS_SOURCES})
target_link_libraries(libgromacs gmx
                      ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES})
set_target_properties(libgromacs PROPERTIES
                      OUTPUT_NAME "gromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${SOVERSION}
                      INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

install(TARGETS libgromacs DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgromacs.pc.cmakein
               ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgromacs${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
