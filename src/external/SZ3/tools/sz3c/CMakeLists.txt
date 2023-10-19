add_library(
        SZ3c
        src/sz3c.cpp
)

target_link_libraries(
        SZ3c
        PUBLIC SZ3
)

target_include_directories(
        SZ3c
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SZ3c>
)

install(TARGETS SZ3c EXPORT SZ3c
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SZ3c)
