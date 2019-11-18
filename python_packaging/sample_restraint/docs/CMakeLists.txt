configure_file(Doxyfile.in Doxyfile)

find_package(Doxygen)
if (DOXYGEN_FOUND)
    add_custom_target(docs
                      COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating developer documentation"
                      VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen needs to be installed to generate the C++ documentation")
endif (DOXYGEN_FOUND)
