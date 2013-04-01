
set(ALEXANDRIA_SOURCES
	xml_util.c
	poldata.c
	poldata_xml.c
	gmx_simple_comm.c
	molselect.c
	gmx_resp.cpp
	gentop_qgen.cpp
	gentop_core.cpp
	gentop_nm2type.cpp
	gentop_vsite.cpp
	molprop.cpp
	molprop_util.cpp
	molprop_xml.cpp
        molprop_sqlite3.cpp
        mymol2.cpp
        gauss_io.cpp
	)
	
add_library(alexandria ${ALEXANDRIA_SOURCES})	
target_link_libraries(alexandria libgromacs)
set_target_properties(alexandria PROPERTIES OUTPUT_NAME "alexandria${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

list(APPEND GMX_EXTRA_LIBRARIES alexandria libgromacs)

add_executable(g_poldata_test g_poldata_test.c)
target_link_libraries(g_poldata_test ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_poldata_test PROPERTIES OUTPUT_NAME "g_poldata_test${GMX_BINARY_SUFFIX}")

add_executable(g_gauss_test g_gauss_test.c)
target_link_libraries(g_gauss_test ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_gauss_test PROPERTIES OUTPUT_NAME "g_gauss_test${GMX_BINARY_SUFFIX}")

add_executable(g_slater_test g_slater_test.c)
target_link_libraries(g_slater_test ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_slater_test PROPERTIES OUTPUT_NAME "g_slater_test${GMX_BINARY_SUFFIX}")

add_executable(gentop gentop.cpp)
target_link_libraries(gentop ${GMX_EXTRA_LIBRARIES})
set_target_properties(gentop PROPERTIES OUTPUT_NAME "gentop${GMX_BINARY_SUFFIX}")

add_executable(gauss2molprop gauss2molprop.cpp)
target_link_libraries(gauss2molprop ${GMX_EXTRA_LIBRARIES})
set_target_properties(gauss2molprop PROPERTIES OUTPUT_NAME "gauss2molprop${GMX_BINARY_SUFFIX}")

add_executable(tune_pol tune_pol.cpp)
target_link_libraries(tune_pol ${GMX_EXTRA_LIBRARIES})
set_target_properties(tune_pol PROPERTIES OUTPUT_NAME "tune_pol${GMX_BINARY_SUFFIX}")

add_executable(tune_fc tune_fc.cpp)
target_link_libraries(tune_fc ${GMX_EXTRA_LIBRARIES})
set_target_properties(tune_fc PROPERTIES OUTPUT_NAME "tune_fc${GMX_BINARY_SUFFIX}")

add_executable(tune_dip tune_dip.cpp)
target_link_libraries(tune_dip ${GMX_EXTRA_LIBRARIES})
set_target_properties(tune_dip PROPERTIES OUTPUT_NAME "tune_dip${GMX_BINARY_SUFFIX}")

add_executable(bastat bastat.cpp)
target_link_libraries(bastat ${GMX_EXTRA_LIBRARIES})
set_target_properties(bastat PROPERTIES OUTPUT_NAME "bastat${GMX_BINARY_SUFFIX}")

add_executable(analyze_mp analyze_mp.cpp molprop_tables.cpp)
target_link_libraries(analyze_mp ${GMX_EXTRA_LIBRARIES})
set_target_properties(analyze_mp PROPERTIES OUTPUT_NAME "analyze_mp${GMX_BINARY_SUFFIX}")

#add_executable(mkice mkice.c)
#target_link_libraries(mkice ${GMX_EXTRA_LIBRARIES})
#set_target_properties(mkice PROPERTIES OUTPUT_NAME "mkice${GMX_BINARY_SUFFIX}")

add_executable(mkyaw mkyaw.c)
target_link_libraries(mkyaw ${GMX_EXTRA_LIBRARIES})
set_target_properties(mkyaw PROPERTIES OUTPUT_NAME "mkyaw${GMX_BINARY_SUFFIX}")

add_executable(gen_table gen_table.c)
target_link_libraries(gen_table ${GMX_EXTRA_LIBRARIES})
set_target_properties(gen_table PROPERTIES OUTPUT_NAME "gen_table${GMX_BINARY_SUFFIX}")

#add_executable(g_mptest g_mptest.c)
#target_link_libraries(g_mptest ${GMX_EXTRA_LIBRARIES})
#set_target_properties(g_mptest PROPERTIES OUTPUT_NAME "g_mptest${GMX_BINARY_SUFFIX}")

add_executable(g_merge_mp g_merge_mp.cpp)
target_link_libraries(g_merge_mp ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_merge_mp PROPERTIES OUTPUT_NAME "g_merge_mp${GMX_BINARY_SUFFIX}")

add_executable(g_mp2csv g_mp2csv.cpp)
target_link_libraries(g_mp2csv ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_mp2csv PROPERTIES OUTPUT_NAME "g_mp2csv${GMX_BINARY_SUFFIX}")

add_executable(gen_ff gen_ff.c)
target_link_libraries(gen_ff ${GMX_EXTRA_LIBRARIES})
set_target_properties(gen_ff PROPERTIES OUTPUT_NAME "gen_ff${GMX_BINARY_SUFFIX}")

add_executable(molprop_test molprop_test.cpp molprop.cpp
molprop_xml.cpp molprop_util.cpp molprop_tables.cpp)
target_link_libraries(molprop_test ${GMX_EXTRA_LIBRARIES})
set_target_properties(molprop_test PROPERTIES OUTPUT_NAME "molprop_test${GMX_BINARY_SUFFIX}")


install(TARGETS alexandria DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)
install(TARGETS 
        gentop
	tune_dip
	tune_pol
        tune_fc
	bastat
	analyze_mp
	mkyaw
	gen_ff
	gen_table
        g_merge_mp
        gauss2molprop
        COMPONENT runtime
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libalexandria.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libalexandria${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
