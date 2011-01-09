include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

if(NOT GMX_FAHCORE)

list(APPEND GMX_EXTRA_LIBRARIES libgromacs)

add_executable(tpbconv tpbconv.c)
target_link_libraries(tpbconv ${GMX_EXTRA_LIBRARIES})
set_target_properties(tpbconv PROPERTIES OUTPUT_NAME "tpbconv${GMX_BINARY_SUFFIX}")

add_executable(g_protonate g_protonate.c)
target_link_libraries(g_protonate ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_protonate PROPERTIES OUTPUT_NAME "g_protonate${GMX_BINARY_SUFFIX}")

add_executable(g_luck g_luck.c)
target_link_libraries(g_luck ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_luck PROPERTIES OUTPUT_NAME "g_luck${GMX_BINARY_SUFFIX}")

install(TARGETS 
        tpbconv
        g_protonate
        g_luck
        COMPONENT runtime
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endif(NOT GMX_FAHCORE)
