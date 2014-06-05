find_package(Pandoc)

if(PANDOC_EXECUTABLE)
    # Do replacement of CMake variables for version strings, etc.
    configure_file(install-guide.md install-guide.md)

    add_custom_target(install-guide.html
        COMMAND pandoc -t html -o install-guide.html install-guide.md -s --toc
        DEPENDS install-guide.md
        )
    add_custom_target(INSTALL
        COMMAND pandoc -t plain -o ../INSTALL install-guide.md
        DEPENDS install-guide.md
        )
    add_custom_target(install-guide)
    add_dependencies(install-guide install-guide.html INSTALL)
endif()
