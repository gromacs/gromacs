
set(ALEXANDRIA_SOURCES
	xml_util.c
	poldata.cpp
	poldata_xml.cpp
	gmx_simple_comm.c
	molselect.cpp
	gmx_resp.cpp
	gentop_qgen.cpp
	gentop_core.cpp
	gentop_nm2type.cpp
	gentop_vsite.cpp
	molprop.cpp
	molprop_util.cpp
	molprop_xml.cpp
        molprop_sqlite3.cpp
        molprop_tables.cpp
        mymol.cpp
        moldip.cpp
        gauss_io.cpp
        gentop.cpp
        tune_fc.cpp
        tune_pol.cpp
        tune_dip.cpp
        bastat.cpp
        gauss2molprop.cpp
        analyze_mp.cpp
        mp2csv.cpp
        merge_mp.cpp
        poldata_test.cpp
        molprop_test.cpp
        gen_ff.cpp
        gen_table.cpp
        alex_modules.cpp
	)
	
add_library(alexandria ${ALEXANDRIA_SOURCES})	
target_link_libraries(alexandria libgromacs)
set_target_properties(alexandria PROPERTIES OUTPUT_NAME "alexandria${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

list(APPEND GMX_EXTRA_LIBRARIES alexandria libgromacs)

add_executable(alex alex.cpp)
target_link_libraries(alex ${GMX_EXTRA_LIBRARIES})
set_target_properties(alex PROPERTIES OUTPUT_NAME "alex${GMX_BINARY_SUFFIX}")

install(TARGETS alexandria DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)
install(TARGETS 
        alex
        COMPONENT runtime
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libalexandria.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libalexandria${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
