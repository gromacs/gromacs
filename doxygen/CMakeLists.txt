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

find_package(PythonInterp)

if (PYTHONINTERP_FOUND AND NOT PYTHON_VERSION_STRING VERSION_LESS "2.6")
    add_custom_target(depgraphs
        ${CMAKE_COMMAND}
        -DSRCDIR=${CMAKE_SOURCE_DIR}
        -DBUILDDIR=${CMAKE_BINARY_DIR}
        -DOUTDIR=${CMAKE_BINARY_DIR}/doxygen/depgraphs
        -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
        -DDOT_EXECUTABLE=${DOXYGEN_DOT_EXECUTABLE}
        -DMODE=GRAPHS
        -P ${CMAKE_SOURCE_DIR}/admin/includedeps.cmake
        COMMENT "Generating include dependency graphs" VERBATIM)
    # TODO: This should be moved away from here and preferably made into
    # a CTest test.
    # TODO: To do the above, issues found by the script should be fixed.
    # A suppression system would also be good to have, since it would be
    # nice to include also legacyheaders in the checks, but those generate
    # a lot of errors about missing documentation.  The script is also
    # currently quite strict, so I don't know whether we really want to
    # enforce all the rules.
    # TODO: Would be nice to not generate the file lists for each of these,
    # but what would be a good shared location for the intermediate files?
    add_custom_target(depcheck
        ${CMAKE_COMMAND}
        -DSRCDIR=${CMAKE_SOURCE_DIR}
        -DBUILDDIR=${CMAKE_BINARY_DIR}
        -DOUTDIR=${CMAKE_BINARY_DIR}/doxygen/depgraphs
        -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
        -DMODE=CHECK
        -P ${CMAKE_SOURCE_DIR}/admin/includedeps.cmake
        COMMENT "Checking include dependencies and documentation" VERBATIM)
endif ()
