if(NOT GMX_FAHCORE)

list(APPEND GMX_EXTRA_LIBRARIES libgromacs)
    readrot.c

add_executable(g_luck g_luck.c)
target_link_libraries(g_luck ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_luck PROPERTIES OUTPUT_NAME "g_luck${GMX_BINARY_SUFFIX}")

install(TARGETS 
        g_luck
        COMPONENT runtime
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endif(NOT GMX_FAHCORE)
