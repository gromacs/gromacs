
set(ALEXANDRIA_SOURCES
	xml_util.c
	poldata.c
	poldata_xml.c
	gaussian_integrals.c
	gmx_simple_comm.c
	gmx_resp.c
	gentop_qgen.c
	gentop_core.c
	gentop_nm2type.c
	gentop_vsite.c
	molprop.c
	molprop_util.c
	molprop_xml.c
	molselect.c
	slater_integrals.cpp
	DSlater_1S_1S.cpp
	DSlater_1S_2S.cpp 
	DSlater_1S_3S.cpp 
#	DSlater_1S_4S.cpp 
#	DSlater_1S_5S.cpp 
	DSlater_2S_2S.cpp 
	DSlater_2S_3S.cpp 
#	DSlater_2S_4S.cpp 
#	DSlater_2S_5S.cpp 
	DSlater_3S_3S.cpp 
#	DSlater_3S_4S.cpp 
#	DSlater_3S_5S.cpp
#	DSlater_4S_4S.cpp 
#	DSlater_4S_5S.cpp 
#	DSlater_5S_5S.cpp 
	Slater_1S_1S.cpp 
	Slater_1S_2S.cpp 
	Slater_1S_3S.cpp 
#	Slater_1S_4S.cpp 
#	Slater_1S_5S.cpp 
	Slater_2S_2S.cpp 
	Slater_2S_3S.cpp 
#	Slater_2S_4S.cpp 
#	Slater_2S_5S.cpp 
#	Slater_4S_4S.cpp 
#	Slater_4S_5S.cpp 
#	Slater_5S_5S.cpp 
#	DSlater_1S_6S.cpp 
#	DSlater_2S_6S.cpp 
#	DSlater_3S_6S.cpp 
#	DSlater_4S_6S.cpp 
#	DSlater_5S_6S.cpp 
#	DSlater_6S_6S.cpp 
#	Slater_1S_6S.cpp 
#	Slater_2S_6S.cpp 
	Slater_3S_3S.cpp 
#	Slater_3S_4S.cpp 
#	Slater_3S_5S.cpp 
#	Slater_3S_6S.cpp 
#	Slater_4S_6S.cpp 
#	Slater_5S_6S.cpp 
#	Slater_6S_6S.cpp
	)
	
add_library(alexandria ${ALEXANDRIA_SOURCES})	
target_link_libraries(alexandria gmxpreprocess md)
set_target_properties(alexandria PROPERTIES OUTPUT_NAME "alexandria${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

list(APPEND GMX_EXTRA_LIBRARIES alexandria gmxpreprocess md)

add_executable(g_poldata_test g_poldata_test.c)
target_link_libraries(g_poldata_test ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_poldata_test PROPERTIES OUTPUT_NAME "g_poldata_test${GMX_BINARY_SUFFIX}")

add_executable(g_gauss_test g_gauss_test.c)
target_link_libraries(g_gauss_test ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_gauss_test PROPERTIES OUTPUT_NAME "g_gauss_test${GMX_BINARY_SUFFIX}")

add_executable(g_slater_test g_slater_test.c)
target_link_libraries(g_slater_test ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_slater_test PROPERTIES OUTPUT_NAME "g_slater_test${GMX_BINARY_SUFFIX}")

add_executable(gentop gentop.c)
target_link_libraries(gentop ${GMX_EXTRA_LIBRARIES})
set_target_properties(gentop PROPERTIES OUTPUT_NAME "gentop${GMX_BINARY_SUFFIX}")

add_executable(analyze_mp analyze_mp.c molprop_tables.c)
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

add_executable(g_mptest g_mptest.c)
target_link_libraries(g_mptest ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_mptest PROPERTIES OUTPUT_NAME "g_mptest${GMX_BINARY_SUFFIX}")

add_executable(g_merge_mp g_merge_mp.c)
target_link_libraries(g_merge_mp ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_merge_mp PROPERTIES OUTPUT_NAME "g_merge_mp${GMX_BINARY_SUFFIX}")

add_executable(g_mp2csv g_mp2csv.c)
target_link_libraries(g_mp2csv ${GMX_EXTRA_LIBRARIES})
set_target_properties(g_mp2csv PROPERTIES OUTPUT_NAME "g_mp2csv${GMX_BINARY_SUFFIX}")

add_executable(addquote addquote.c)
target_link_libraries(addquote ${GMX_EXTRA_LIBRARIES})
set_target_properties(addquote PROPERTIES OUTPUT_NAME "addquote${GMX_BINARY_SUFFIX}")

install(TARGETS alexandria DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)
install(TARGETS 
        gentop
	analyze_mp
	mkyaw
	addquote
	gen_table
        COMPONENT runtime
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libalexandria.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libalexandria${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
