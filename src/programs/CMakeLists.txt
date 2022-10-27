#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2010- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

file(GLOB MDRUN_SOURCES mdrun/*.cpp)
# make an "object library" that we can re-use for multiple targets
add_library(mdrun_objlib OBJECT ${MDRUN_SOURCES})
gmx_target_compile_options(mdrun_objlib)
target_compile_definitions(mdrun_objlib PRIVATE HAVE_CONFIG_H TMPI_USE_VISIBILITY)
target_include_directories(mdrun_objlib SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
# Should be possible to remove this when resolving #3290
target_include_directories(mdrun_objlib SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
target_link_libraries(mdrun_objlib PRIVATE
        common
        legacy_api
        legacy_modules
        utility
        )

if(GMX_FAHCORE)
    # The lack of a real source file here alongside the object library
    # may break some generators, according to CMake documentation. If
    # so, we can consider adding some dummy file to make it work.
    add_library(fahcore $<TARGET_OBJECTS:mdrun_objlib>)
    target_link_libraries(fahcore PRIVATE ${GMX_COMMON_LIBRARIES} legacy_api)
else()
    file(GLOB GMX_MAIN_SOURCES gmx.cpp legacymodules.cpp)
    add_library(gmx_objlib OBJECT ${GMX_MAIN_SOURCES})
    target_link_libraries(gmx_objlib PRIVATE common legacy_api)
    # TODO: Remove `legacy_modules` once specific modules are explicitly linked.
    target_link_libraries(gmx_objlib PRIVATE legacy_modules)
    target_link_libraries(gmx_objlib PRIVATE utility commandline)
    target_include_directories(gmx_objlib SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
    target_include_directories(gmx_objlib SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
    target_compile_definitions(gmx_objlib PRIVATE TMPI_USE_VISIBILITY)
    add_executable(gmx
        $<TARGET_OBJECTS:gmx_objlib>
        $<TARGET_OBJECTS:mdrun_objlib>)
    add_executable(Gromacs::gmx ALIAS gmx)
    gmx_target_compile_options(gmx)
    target_compile_definitions(gmx PRIVATE HAVE_CONFIG_H)
    target_link_libraries(gmx PRIVATE
                          common
                          libgromacs
                          ${GMX_COMMON_LIBRARIES}
                          ${GMX_EXE_LINKER_FLAGS})
    # Ensure that on Windows the application code page is set to utf8 mode so the behavior is
    # like Linux. For proper UTF8 console output, the user needs to enable utf8 (chcp 65001).
    if (MSVC)
        target_sources(gmx PRIVATE utf8.manifest)
    endif()
    set(BINARY_NAME "gmx${GMX_BINARY_SUFFIX}")
    set_target_properties(gmx PROPERTIES
        OUTPUT_NAME "${BINARY_NAME}")
    # Add the executable to the INSTALL target and to the libgromacs EXPORT configuration.
    # The GROMACS installation then provides an importable (e.g. find_package(GROMACS))
    # Gromacs::gmx target, which may be used as the executable COMMAND argument for functions like
    # add_custom_command() or add_test(). The full path to the installed binary may retrieved
    # from the LOCATION property or with a $<TARGET_FILE:Gromacs:gmx> generator expression
    # (where generator expressions are allowed).
    install(TARGETS gmx
            EXPORT libgromacs
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    ########################
    # Completion generation

    include(gmxCustomCommandUtilities)

    set(COMPLETION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/completion)
    # Using GMX_BUILD_HELP here is somewhat confusing, but the conditions when
    # this can be done are exactly the same (ability to run the compiled
    # binaries).
    if (GMX_BUILD_HELP)
        gmx_add_custom_output_target(completion OUTPUT STAMP
            COMMAND ${CMAKE_COMMAND}
                -D GMX_EXECUTABLE=$<TARGET_FILE:gmx>
                -D ERRORS_ARE_FATAL=${GMX_BUILD_HELP_FORCE}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/BuildCompletions.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            DEPENDS gmx ${CMAKE_CURRENT_SOURCE_DIR}/BuildCompletions.cmake
            COMMENT "Generating command-line completions for programs")
        set_target_properties(completion PROPERTIES EXCLUDE_FROM_ALL OFF)
        set_directory_properties(PROPERTIES
            ADDITIONAL_MAKE_CLEAN_FILES "completion")
        set(COMPLETION_DIR ${CMAKE_CURRENT_BINARY_DIR}/completion)
    endif()
    if (SOURCE_IS_SOURCE_DISTRIBUTION OR GMX_BUILD_HELP)
        install(DIRECTORY ${COMPLETION_DIR}/
                DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime OPTIONAL)
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gmx-completion-${BINARY_NAME}.bash
             "complete -o nospace -F _gmx_compl ${BINARY_NAME}")
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gmx-completion-${BINARY_NAME}.bash
                DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
    endif()
    gmx_cpack_add_generated_source_directory(completion)

    if(BUILD_TESTING)
        add_subdirectory(mdrun/tests)
    endif()
endif()
