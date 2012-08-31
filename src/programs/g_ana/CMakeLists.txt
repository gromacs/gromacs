add_executable(g_ana g_ana.cpp)
target_link_libraries(g_ana libgromacs)
set_target_properties(g_ana PROPERTIES OUTPUT_NAME "g_ana${GMX_BINARY_SUFFIX}")

# FIXME: The list of files should not be hardcoded
add_custom_command(OUTPUT help-topics.txt g_ana-select.txt
                   COMMAND g_ana help -export
                   DEPENDS g_ana)
add_custom_target(man
                  rst2man.py g_ana-select.txt ${MAN1_PATH}/g_ana-select.1
                  DEPENDS g_ana-select.txt)
add_custom_target(html
                  COMMAND rst2html.py g_ana-select.txt ${CMAKE_BINARY_DIR}/g_ana-select.html
                  COMMAND rst2html.py help-topics.txt ${CMAKE_BINARY_DIR}/help-topics.html
                  DEPENDS help-topics.txt g_ana-select.txt)
add_custom_target(latex
                  COMMAND rst2latex.py g_ana-select.txt ${CMAKE_BINARY_DIR}/g_ana-select.tex
                  COMMAND rst2latex.py help-topics.txt ${CMAKE_BINARY_DIR}/help-topics.tex
                  DEPENDS help-topics.txt g_ana-select.txt)

install(TARGETS g_ana
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})
