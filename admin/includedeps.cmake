function (generate_module_file_list SRCDIR OUTFILE)
    set(PATH_LIST)
    foreach (MODULE analysisdata commandline linearalgebra onlinehelp
                    options selection trajectoryanalysis utility)
        list(APPEND PATH_LIST "${SRCDIR}/src/gromacs/${MODULE}/*.cpp")
        list(APPEND PATH_LIST "${SRCDIR}/src/gromacs/${MODULE}/*.h")
    endforeach ()
    list(APPEND PATH_LIST "${SRCDIR}/src/testutils/*.cpp")
    list(APPEND PATH_LIST "${SRCDIR}/src/testutils/*.h")
    set(FILE_LIST)
    foreach (PATH_EXPR ${PATH_LIST})
        file(GLOB_RECURSE FOUND_FILES ${PATH_EXPR})
        list(APPEND FILE_LIST ${FOUND_FILES})
    endforeach ()
    string(REPLACE ";" "\n" FILE_LIST "${FILE_LIST}")
    file(WRITE ${OUTFILE} "${FILE_LIST}")
endfunction ()

function (generate_installed_file_list SRCDIR BUILDDIR OUTFILE)
    file(GLOB_RECURSE INSTALL_FILE_LIST "${BUILDDIR}/cmake_install.cmake")
    set(MATCH_REGEX "${SRCDIR}/.*\\.h")
    set(HEADER_LIST)
    foreach (INSTALL_FILE ${INSTALL_FILE_LIST})
        file(STRINGS ${INSTALL_FILE} HEADER_LINES REGEX "${MATCH_REGEX}")
        foreach (HEADER_LINE ${HEADER_LINES})
            string (REGEX MATCH "${MATCH_REGEX}" HEADER "${HEADER_LINE}")
            list(APPEND HEADER_LIST "${HEADER}")
        endforeach ()
    endforeach ()
    string(REPLACE ";" "\n" HEADER_LIST "${HEADER_LIST}")
    file(WRITE ${OUTFILE} "${HEADER_LIST}")
endfunction ()

if (NOT DEFINED SRCDIR OR NOT DEFINED BUILDDIR OR NOT DEFINED OUTDIR)
    message(FATAL_ERROR "Required input variable (SRCDIR, BUILDDIR, OUTDIR) not set")
endif ()

if (NOT DEFINED PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE python)
endif ()

if (NOT DEFINED MODE OR MODE STREQUAL "CHECK")
    set(GRAPHOPTIONS)
elseif (MODE STREQUAL "GRAPHS")
    set(GRAPHOPTIONS
        --module-graph module-deps.dot --module-file-graphs
        -o ${OUTDIR})
else ()
    message(FATAL_ERROR "Unknown mode ${MODE}")
endif ()

file(MAKE_DIRECTORY ${OUTDIR})
generate_module_file_list(${SRCDIR} ${OUTDIR}/module-files.txt)
generate_installed_file_list(${SRCDIR} ${BUILDDIR} ${OUTDIR}/installed-headers.txt)
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SRCDIR}/admin/includedeps.py
                        -f ${OUTDIR}/module-files.txt
                        --installed ${OUTDIR}/installed-headers.txt
                        -R ${SRCDIR}/src -R ${BUILDDIR}/src
                        #-I ${SRCDIR}/src/gromacs/legacyheaders
                        -I ${BUILDDIR}/src/gromacs/utility
                        ${GRAPHOPTIONS})

if (MODE STREQUAL "GRAPHS" AND DOT_EXECUTABLE)
    file(GLOB DOT_INPUT_FILES ${OUTDIR}/*.dot)
    execute_process(COMMAND ${DOT_EXECUTABLE} -Tpng -O ${DOT_INPUT_FILES})
endif ()
