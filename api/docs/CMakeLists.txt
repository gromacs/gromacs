#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2017- The GROMACS Authors
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

################################
# C++ external API documentation
################################
set(GMX_API_DOC_OUTPUT ${CMAKE_CURRENT_BINARY_DIR})

find_package(Doxygen QUIET)
if (NOT DOXYGEN_FOUND)
    function (add_disabled_doxygen_target TARGET)
        add_custom_target(${TARGET}
            COMMAND ${CMAKE_COMMAND} -E echo
                "Cannot build/run target ${TARGET} because Doxygen is not available"
            COMMENT "Doxygen not available" VERBATIM)
    endfunction()
    add_disabled_doxygen_target(public-api-xml)
    add_disabled_doxygen_target(gmxapi-cppdocs)
    add_disabled_doxygen_target(gmxapi-cppdocs-dev)
else ()
    set(legacy_api_SOURCE_DIR ${CMAKE_SOURCE_DIR}/api/legacy)
    set(legacy_api_BINARY_DIR ${CMAKE_BINARY_DIR}/api/legacy)
    set(legacy_api_INCLUDE_DIRECTORIES ${legacy_api_SOURCE_DIR}/include ${legacy_api_BINARY_DIR}/include)
    set(GMX_API_INCLUDE_DIRECTORIES ${legacy_api_INCLUDE_DIRECTORIES})

    if (GMXAPI)
        file(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/docs/api-user" GMXAPI_USERDOC_OUTPUT)
        file(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/docs/api-dev" GMXAPI_DEVDOC_OUTPUT)

        file(MAKE_DIRECTORY ${GMXAPI_USERDOC_OUTPUT})
        file(MAKE_DIRECTORY ${GMXAPI_DEVDOC_OUTPUT})

        set(GMXAPI_SOURCE_DIR ${CMAKE_SOURCE_DIR}/api/gmxapi)
        set(GMXAPI_BINARY_DIR ${CMAKE_BINARY_DIR}/api/gmxapi)

        # Use GMXAPI_PUBLIC_HEADERS to populate source files. Make a local copy with
        # the list converted to a space-delimited string.
        set(GMXAPI_PUBLIC_HEADERS ${GMXAPI_SOURCE_DIR}/include ${GMXAPI_BINARY_DIR}/include)
        string(REPLACE ";" " " GMXAPI_PUBLIC_HEADERS "${GMXAPI_PUBLIC_HEADERS}")
        configure_file(Doxyfile-gmxapi0.cmakein Doxyfile-gmxapi0 @ONLY)

        # Use all gmxapi sources to populate developer docs.
        configure_file(Doxyfile-gmxapi0-dev.cmakein Doxyfile-gmxapi0-dev @ONLY)

        list(APPEND GMX_API_INCLUDE_DIRECTORIES ${GMXAPI_PUBLIC_HEADERS})
    endif()

    if (GMX_INSTALL_NBLIB_API)
        get_target_property(NBLIB_SOURCE_DIR nblib SOURCE_DIR)
        get_target_property(NBLIB_BINARY_DIR nblib BINARY_DIR)
        set(NBLIB_INCLUDE_DIRECTORIES ${NBLIB_SOURCE_DIR}/include ${NBLIB_BINARY_DIR}/include)
        list(APPEND GMX_API_INCLUDE_DIRECTORIES ${NBLIB_INCLUDE_DIRECTORIES})
    endif ()

    string(REPLACE ";" " " GMX_API_INCLUDE_DIRECTORIES "${GMX_API_INCLUDE_DIRECTORIES}")
    configure_file(Doxyfile.cmakein Doxyfile-api @ONLY)

    add_custom_target( public-api-xml
                       COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile-api
                       BYPRODUCTS ${GMX_API_DOC_OUTPUT}/xml
                       COMMENT "Generating API user documentation for ${GMX_API_DOC_OUTPUT}"
                       VERBATIM )

    if (GMXAPI)
        # Define targets gmxapi-cppdocs and gmxapi-cppdocs-dev when GMXAPI is enabled.
        # gmxapi-cppdocs produces public API documentation in docs/api-user
        # extracted from the gmxapi public headers files.
        # gmxapi-cppdocs-dev produces more complete documentation in docs/api-dev
        # with internal interfaces and implementation details.
        add_custom_target( gmxapi-cppdocs
            COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile-gmxapi0
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API user documentation for ${GMXAPI_USERDOC_OUTPUT}"
            VERBATIM )
        add_custom_target( gmxapi-cppdocs-dev
            COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile-gmxapi0-dev
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API developer documentation for ${GMXAPI_DEVDOC_OUTPUT}"
            VERBATIM )
    endif()
endif ()
