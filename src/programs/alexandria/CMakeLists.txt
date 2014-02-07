
set(ALEXANDRIA_SOURCES
     coulombintegrals/DSlater_1S_1S.cpp
     coulombintegrals/DSlater_1S_2S.cpp
     coulombintegrals/DSlater_1S_3S.cpp
     coulombintegrals/DSlater_1S_4S.cpp
     #coulombintegrals/DSlater_1S_5S.cpp
     #coulombintegrals/DSlater_1S_6S.cpp
     coulombintegrals/DSlater_2S_2S.cpp
     coulombintegrals/DSlater_2S_3S.cpp
     coulombintegrals/DSlater_2S_4S.cpp
     #coulombintegrals/DSlater_2S_5S.cpp
     #coulombintegrals/DSlater_2S_6S.cpp
     coulombintegrals/DSlater_3S_3S.cpp
     coulombintegrals/DSlater_3S_4S.cpp
     #coulombintegrals/DSlater_3S_5S.cpp
     #coulombintegrals/DSlater_3S_6S.cpp
     coulombintegrals/DSlater_4S_4S.cpp
     #coulombintegrals/DSlater_4S_5S.cpp
     #coulombintegrals/DSlater_4S_6S.cpp
     #coulombintegrals/DSlater_5S_5S.cpp
     #coulombintegrals/DSlater_5S_6S.cpp
     #coulombintegrals/DSlater_6S_6S.cpp
     coulombintegrals/Slater_1S_1S.cpp
     coulombintegrals/Slater_1S_2S.cpp
     coulombintegrals/Slater_1S_3S.cpp
     coulombintegrals/Slater_1S_4S.cpp
     #coulombintegrals/Slater_1S_5S.cpp
     #coulombintegrals/Slater_1S_6S.cpp
     coulombintegrals/Slater_2S_2S.cpp
     coulombintegrals/Slater_2S_3S.cpp
     coulombintegrals/Slater_2S_4S.cpp
     #coulombintegrals/Slater_2S_5S.cpp
     #coulombintegrals/Slater_2S_6S.cpp
     coulombintegrals/Slater_3S_3S.cpp
     coulombintegrals/Slater_3S_4S.cpp
     #coulombintegrals/Slater_3S_5S.cpp
     #coulombintegrals/Slater_3S_6S.cpp
     coulombintegrals/Slater_4S_4S.cpp
     #coulombintegrals/Slater_4S_5S.cpp
     #coulombintegrals/Slater_4S_6S.cpp
     #coulombintegrals/Slater_5S_5S.cpp
     #coulombintegrals/Slater_5S_6S.cpp
     #coulombintegrals/Slater_6S_6S.cpp
     coulombintegrals/gaussian_integrals.c
     coulombintegrals/slater_integrals.cpp
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
        analyze.cpp
        mp2csv.cpp
        merge_mp.cpp
        poldata_test.cpp
        molprop_test.cpp
        gen_ff.cpp
        gen_table.cpp
        alex_modules.cpp
	)
	
add_library(libalexandria ${ALEXANDRIA_SOURCES})	
target_link_libraries(libalexandria libgromacs)

set_target_properties(libalexandria PROPERTIES
                      OUTPUT_NAME "alexandria${GMX_LIBS_SUFFIX}"
                      SOVERSION ${LIBRARY_SOVERSION}
                      VERSION ${LIBRARY_VERSION}
                      COMPILE_FLAGS "${OpenMP_C_FLAGS}")

list(APPEND GMX_EXTRA_LIBRARIES libalexandria libgromacs)

add_executable(alexandria alexandria.cpp)
target_link_libraries(alexandria ${GMX_EXTRA_LIBRARIES})
set_target_properties(alexandria PROPERTIES OUTPUT_NAME "alexandria${GMX_BINARY_SUFFIX}")

install(TARGETS 
        alexandria
        COMPONENT runtime
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libalexandria.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libalexandria${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
