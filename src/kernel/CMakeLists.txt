include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

set(GMXPREPROCESS_SOURCES 
    sorting.c)

add_library(gmxpreprocess ${GMXPREPROCESS_SOURCES})
target_link_libraries(gmxpreprocess libgromacs)
set_target_properties(gmxpreprocess PROPERTIES OUTPUT_NAME "gmxpreprocess${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

if(NOT GMX_FAHCORE)

list(APPEND GMX_EXTRA_LIBRARIES gmxpreprocess libgromacs)

add_executable(tpbconv tpbconv.c)
target_link_libraries(tpbconv ${GMX_EXTRA_LIBRARIES})
set_target_properties(tpbconv PROPERTIES OUTPUT_NAME "tpbconv${GMX_BINARY_SUFFIX}")

add_executable(g_protonate g_protonate.c)
target_link_libraries(g_protonate ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_protonate PROPERTIES OUTPUT_NAME "g_protonate${GMX_BINARY_SUFFIX}")

add_executable(g_luck g_luck.c)
target_link_libraries(g_luck ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_luck PROPERTIES OUTPUT_NAME "g_luck${GMX_BINARY_SUFFIX}")

install(TARGETS gmxpreprocess DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)
install(TARGETS 
        tpbconv
        g_protonate
        g_luck
        COMPONENT runtime
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endif(NOT GMX_FAHCORE)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmxpreprocess.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmxpreprocess.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmxpreprocess.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgmxpreprocess${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
