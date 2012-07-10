#TODO: add check that source doesn't contain any untracked files
get_filename_component(CMAKE_BINARY_DIR ${CPACK_OUTPUT_CONFIG_FILE} PATH)
if(CPACK_SOURCE_PACKAGE_FILE_NAME AND NOT EXISTS "${CMAKE_BINARY_DIR}/man/man1/ngmx.1")
    message(FATAL_ERROR 
        "To generate correct source package all man pages need to be generated. "
        "The man pages are automatically build together with the binaries. "
        "Make sure to build all binaries (e.g. GMX_X11=on). ${CMAKE_SOURCE_DIR}/man/man1/ngmx.1")
endif()
