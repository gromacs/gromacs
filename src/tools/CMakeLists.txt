# List of programs with single corresponding .c source file,
# used to create build rules automatically.
#
set(GMX_TOOLS_PROGRAMS
    g_tune_pme g_pme_error
    )

foreach(TOOL ${GMX_TOOLS_PROGRAMS})
    add_executable(${TOOL} ${TOOL}.c)
    gmx_add_man_page(${TOOL})
    target_link_libraries(${TOOL} libgromacs ${GMX_EXE_LINKER_FLAGS})
    set_target_properties(${TOOL} PROPERTIES OUTPUT_NAME "${TOOL}${GMX_BINARY_SUFFIX}")
endforeach(TOOL ${GMX_TOOLS_PROGRAMS})

install(TARGETS ${GMX_TOOLS_PROGRAMS}
        DESTINATION ${BIN_INSTALL_DIR}
        COMPONENT runtime)
