#TODO: add check that source doesn't contain any untracked files
if(CPACK_SOURCE_PACKAGE_FILE_NAME) #building source package
    get_filename_component(CMAKE_BINARY_DIR ${CPACK_OUTPUT_CONFIG_FILE} PATH)
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/share/man/man1/gmx-view.1")
        message(FATAL_ERROR
            "To create a complete source package all man pages need to be generated. "
            "You need to run 'make man' or set GMX_BUILD_MANPAGES=ON to get "
            "them automatically built together with the binaries."
    endif()
endif()
