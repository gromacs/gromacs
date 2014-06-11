#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

set(LIBGROMACS_SOURCES)

function (gmx_install_headers DESTINATION)
    if (NOT GMX_BUILD_MDRUN_ONLY)
        if (DESTINATION)
            set(DESTINATION ${INCL_INSTALL_DIR}/gromacs/${DESTINATION})
        else()
            set(DESTINATION ${INCL_INSTALL_DIR}/gromacs)
        endif()
        install(FILES ${ARGN} DESTINATION ${DESTINATION} COMPONENT development)
    endif()
endfunction ()

add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(gmxpreprocess)
add_subdirectory(commandline)
add_subdirectory(fft)
add_subdirectory(linearalgebra)
add_subdirectory(math)
add_subdirectory(random)
add_subdirectory(onlinehelp)
add_subdirectory(options)
add_subdirectory(timing)
add_subdirectory(utility)
add_subdirectory(fileio)
add_subdirectory(swap)
add_subdirectory(essentialdynamics)
add_subdirectory(pulling)
add_subdirectory(simd)
add_subdirectory(imd)
if (NOT GMX_BUILD_MDRUN_ONLY)
    add_subdirectory(legacyheaders)
    add_subdirectory(gmxana)
    add_subdirectory(statistics)
    add_subdirectory(analysisdata)
    add_subdirectory(selection)
    add_subdirectory(trajectoryanalysis)
    add_subdirectory(tools)
endif()

list(APPEND LIBGROMACS_SOURCES ${GMXLIB_SOURCES} ${MDLIB_SOURCES})

# This would be the standard way to include thread_mpi, but
# we want libgromacs to link the functions directly
#if(GMX_THREAD_MPI)
#    add_subdirectory(thread_mpi)
#endif()
#target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES} ${THREAD_MPI_LIB})

tmpi_get_source_list(THREAD_MPI_SOURCES ${CMAKE_SOURCE_DIR}/src/external/thread_mpi/src)
list(APPEND LIBGROMACS_SOURCES ${THREAD_MPI_SOURCES})

file(GLOB LIBGROMACS_HEADERS *.h)
configure_file(version.h.cmakein version.h)
gmx_install_headers("" ${LIBGROMACS_HEADERS})
gmx_install_headers("" ${CMAKE_CURRENT_BINARY_DIR}/version.h)

# Add target that generates baseversion-gen.c every time make is run
# if git version info is requested, or create it statically.
# This code is here instead of utility/CMakeLists.txt because CMake
# ignores set_source_file_properties from subdirectories.
set(GENERATED_VERSION_FILE
    ${CMAKE_CURRENT_BINARY_DIR}/utility/baseversion-gen.c)
set(GENERATED_VERSION_FILE_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/utility/baseversion-gen.c.cmakein)
if (GMX_GIT_VERSION_INFO)
    add_custom_target(gmx-version ALL
            COMMAND ${CMAKE_COMMAND}
                -D GIT_EXECUTABLE="${GIT_EXECUTABLE}"
                -D PROJECT_VERSION="${PROJECT_VERSION}"
                -D PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                -D VERSION_CMAKEIN=${GENERATED_VERSION_FILE_SOURCE}
                -D VERSION_OUT=${GENERATED_VERSION_FILE}
                -P ${CMAKE_SOURCE_DIR}/cmake/gmxGenerateVersionInfo.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS ${GENERATED_VERSION_FILE_SOURCE}
            COMMENT "Generating git version information")
    set_source_files_properties(${GENERATED_VERSION_FILE}
                                PROPERTIES GENERATED true)
else()
    set(GMX_PROJECT_VERSION_STR ${PROJECT_VERSION})
    configure_file(${GENERATED_VERSION_FILE_SOURCE} ${GENERATED_VERSION_FILE})
endif()
list(APPEND LIBGROMACS_SOURCES ${GENERATED_VERSION_FILE})

# apply gcc 4.4.x bug workaround
if(GMX_USE_GCC44_BUG_WORKAROUND)
   include(gmxGCC44O3BugWorkaround)
   gmx_apply_gcc44_bug_workaround("gmxlib/bondfree.c")
   gmx_apply_gcc44_bug_workaround("mdlib/force.c")
   gmx_apply_gcc44_bug_workaround("mdlib/constr.c")
endif()

if(GMX_USE_TNG)
    # Set variables to control the TNG build for GROMACS purposes
    set(TNG_BUILD_FORTRAN OFF CACHE BOOL "Build Fortran compatible TNG library and examples for testing")
    set(TNG_BUILD_EXAMPLES OFF CACHE BOOL "Build examples showing usage of the TNG API")
    set(TNG_BUILD_COMPRESSION_TESTS OFF CACHE BOOL "Build tests of the TNG compression library")
    set(TNG_BUILD_DOCUMENTATION OFF CACHE BOOL "Use Doxygen to create the HTML based TNG API documentation")
    set(TNG_BUILD_TEST OFF CACHE BOOL "Build TNG testing binary.")
    # We only want to build the targets that are actually used,
    # e.g. so that the libraries are not automatically installed
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/external/tng_io tng_io EXCLUDE_FROM_ALL)
    # TNG puts these variables into the (common) cache, which is too
    # noisy for GROMACS use
    mark_as_advanced(TNG_BUILD_WITH_ZLIB)
    mark_as_advanced(TNG_BUILD_FORTRAN)
    mark_as_advanced(TNG_BUILD_EXAMPLES)
    mark_as_advanced(TNG_BUILD_COMPRESSION_TESTS)
    mark_as_advanced(TNG_BUILD_DOCUMENTATION)
    mark_as_advanced(TNG_BUILD_TEST)
    mark_as_advanced(TNG_EXAMPLE_FILES_DIR)

    # TODO simplify this when we require CMake version >= 2.8.9.
    if(BUILD_SHARED_LIBS AND NOT MSVC)
        # Windows does shared libraries in a fundamentally different
        # way from Unix. Not sure how icc on Windows behaves.
        if(CMAKE_VERSION VERSION_LESS "2.8.9")
            # TNG only compiles C code
            GMX_TEST_CFLAG(CFLAG_FPIC "-fPIC" C_FPIC_TNG_FLAG)
            if(NOT C_FPIC_TNG_FLAG)
                # gcc, icc and clang all support -fPIC, and platforms
                # not using those compilers likely won't build shared
                # libraries (e.g. xlc on BlueGene/Q), so in practice
                # this fatal error will not be reached.
                message(FATAL_ERROR "Cannot build TNG support into a GROMACS build using shared libraries with CMake version < 2.8.9 unless the compiler supports the -fPIC flag. Use a more recent CMake, a different compiler, or -DBUILD_SHARED_LIBS=OFF")
            endif()
            set_target_properties(tng_io_objlib PROPERTIES COMPILE_FLAGS "-fPIC")
            set_target_properties(tng_compress_objlib PROPERTIES COMPILE_FLAGS "-fPIC")
        else()
            set_target_properties(tng_io_objlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
            set_target_properties(tng_compress_objlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
        endif()
    endif()

    list(APPEND LIBGROMACS_SOURCES $<TARGET_OBJECTS:tng_io_objlib> $<TARGET_OBJECTS:tng_compress_objlib>)
endif()

add_library(libgromacs ${LIBGROMACS_SOURCES})
if (GMX_GIT_VERSION_INFO)
    add_dependencies(libgromacs gmx-version)
endif()

# Recent versions of gcc and clang give warnings on scanner.cpp, which
# is a generated source file. These are awkward to suppress inline, so
# we do it in the compilation command (after testing that the compiler
# supports the suppressions). Setting the properties only works after
# the related target has been created, e.g. after when the file is
# used with add_library().
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-unused-parameter HAS_NO_UNUSED_PARAMETER)
if (HAS_NO_UNUSED_PARAMETER)
    set(_scanner_cpp_compiler_flags "${_scanner_cpp_compiler_flags} -Wno-unused-parameter")
endif()
check_cxx_compiler_flag(-Wno-deprecated-register HAS_NO_DEPRECATED_REGISTER)
if (HAS_NO_DEPRECATED_REGISTER)
    set(_scanner_cpp_compiler_flags "${_scanner_cpp_compiler_flags} -Wno-deprecated-register")
else()
    check_cxx_compiler_flag(-Wno-deprecated HAS_NO_DEPRECATED)
    if (HAS_NO_DEPRECATED)
        set(_scanner_cpp_compiler_flags "${_scanner_cpp_compiler_flags} -Wno-deprecated")
    endif()
endif()
set_source_files_properties(selection/scanner.cpp PROPERTIES COMPILE_FLAGS "${_scanner_cpp_compiler_flags}")

target_link_libraries(libgromacs ${GMX_GPU_LIBRARIES}
                      ${GMX_EXTRA_LIBRARIES}
                      ${EXTRAE_LIBRARIES}
                      ${FFT_LIBRARIES} ${LINEAR_ALGEBRA_LIBRARIES}
                      ${XML_LIBRARIES}
                      ${THREAD_LIB} ${GMX_SHARED_LINKER_FLAGS})
set_target_properties(libgromacs PROPERTIES
                      OUTPUT_NAME "gromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${LIBRARY_SOVERSION}
                      VERSION ${LIBRARY_VERSION}
                      COMPILE_FLAGS "${OpenMP_C_FLAGS}")

# Only install the library in mdrun-only mode if it is actually necessary
# for the binary
if (NOT GMX_BUILD_MDRUN_ONLY OR BUILD_SHARED_LIBS)
    install(TARGETS libgromacs DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)
endif()

if (NOT GMX_BUILD_MDRUN_ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgromacs.pc.cmakein
                   ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc @ONLY)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc
            DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
            RENAME "libgromacs${GMX_LIBS_SUFFIX}.pc"
            COMPONENT development)
endif()

if (INSTALL_CUDART_LIB) #can be set manual by user
    if (GMX_GPU)
        foreach(CUDA_LIB ${CUDA_LIBRARIES})
            string(REGEX MATCH "cudart" IS_CUDART ${CUDA_LIB})
            if(IS_CUDART) #libcuda should not be installed
                #install also name-links (linker uses those)
                file(GLOB CUDA_LIBS ${CUDA_LIB}*)
                install(FILES ${CUDA_LIBS} DESTINATION
                    ${LIB_INSTALL_DIR} COMPONENT libraries)
            endif()
        endforeach()
    else()
        message(WARNING "INSTALL_CUDART_LIB only makes sense with GMX_GPU")
    endif()
endif()
