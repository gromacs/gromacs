#
# This file is part of the alexandria code base
#
# Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

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
	gentop_core.cpp
	gentop_vsite.cpp
        qgen_eem.cpp
	qgen_resp.cpp
        qgen_resp_atom.cpp
	gmx_simple_comm.cpp
	molprop.cpp
	molprop_util.cpp
	molprop_xml.cpp
	molselect.cpp
	poldata.cpp
	poldata_xml.cpp
	poldata-low.cpp
	xml_util.c
	mymol.cpp
        alex_modules.cpp
        analyze.cpp
        bastat.cpp
        categories.cpp
        composition.cpp
        gauss2molprop.cpp
        gauss_io.cpp
        gen_table.cpp
        gentop.cpp
        getmdlogger.cpp
        merge_mp.cpp
        moldip.cpp
        molprop_check.cpp
        molprop_sqlite3.cpp
        molprop_tables.cpp
        molprop_test.cpp
        mp2csv.cpp
        nmsimplex.c
        phase.cpp
        plistwrapper.cpp
        poldata_test.cpp
        regression.cpp
        stringutil.cpp
        tune_dip.cpp
        tune_fc.cpp
        tune_pol.cpp
	optparam.cpp
#        gen_ff.cpp
	)
	
#add_library(libalexandria ${ALEXANDRIA_SOURCES})	
#target_link_libraries(libalexandria libgromacs ${GMX_EXTRA_LIBRARIES})

add_library(alexandria_objlib OBJECT ${ALEXANDRIA_SOURCES})

#set_target_properties(alexandria_objlib PROPERTIES
#                      OUTPUT_NAME "alexandria${GMX_LIBS_SUFFIX}"
#                      COMPILE_FLAGS "${OpenMP_C_FLAGS}")

add_executable(alexandria alexandria.cpp
         $<TARGET_OBJECTS:alexandria_objlib>
         )
target_link_libraries(alexandria libgromacs ${GMX_EXTRA_LIBRARIES} ${XML_LIBRARIES}
         ${GMX_EXE_LINKER_FLAGS})
set(BINARY_NAME "alexandria${GMX_BINARY_SUFFIX}")
set_target_properties(alexandria PROPERTIES 
                      OUTPUT_NAME "${BINARY_NAME}")

install(TARGETS alexandria
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libalexandria.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libalexandria${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

