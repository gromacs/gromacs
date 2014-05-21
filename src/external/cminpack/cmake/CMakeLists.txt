set(PKG_DESC "CMinPack")
set(PKG_EXTERNAL_DEPS "")
set(pkg_conf_file ${CMAKE_CURRENT_BINARY_DIR}/cminpack.pc)
configure_file(cminpack.pc.in ${pkg_conf_file} @ONLY)
install(FILES ${pkg_conf_file}
    DESTINATION ${CMINPACK_LIB_INSTALL_DIR}/pkgconfig/ COMPONENT pkgconfig)

install(FILES FindCMinpack.cmake DESTINATION ${CMAKE_ROOT}/Modules)
