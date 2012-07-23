find_package(Doxygen)

########################################################################
# Doxygen configuration
########################################################################
if (DOXYGEN_FOUND)
    FILE(GLOB NB_KERNEL_DIRS
         ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/nonbonded/nb_kernel_*/)
    LIST(REMOVE_ITEM NB_KERNEL_DIRS
         ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/nonbonded/nb_kernel_c)
    FOREACH(NB_KERNEL_DIR ${NB_KERNEL_DIRS})
        SET(NB_KERNEL_DIRS_TO_IGNORE_IN_DOXYGEN
            "${NB_KERNEL_DIRS_TO_IGNORE_IN_DOXYGEN} \\\n                         ${NB_KERNEL_DIR}")
    ENDFOREACH(NB_KERNEL_DIR)
    CONFIGURE_FILE(Doxyfile-common.cmakein Doxyfile-common)
    CONFIGURE_FILE(Doxyfile-full.cmakein Doxyfile-full)
    CONFIGURE_FILE(Doxyfile-lib.cmakein Doxyfile-lib)
    CONFIGURE_FILE(Doxyfile-user.cmakein Doxyfile-user)
    FILE(COPY index.html DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_target(doc-full
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-full
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating full documentation with Doxygen" VERBATIM)
    add_custom_target(doc-lib
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-lib
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating library documentation with Doxygen" VERBATIM)
    add_custom_target(doc-user
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-user
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating user documentation with Doxygen" VERBATIM)
    add_custom_target(doc-all)
    add_dependencies(doc-all doc-full doc-lib doc-user)
endif (DOXYGEN_FOUND)
