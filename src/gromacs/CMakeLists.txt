#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

if (GMX_CLANG_CUDA)
    include(gmxClangCudaUtils)
endif()

set_property(GLOBAL PROPERTY GMX_LIBGROMACS_SOURCES)
set_property(GLOBAL PROPERTY GMX_INSTALLED_HEADERS)

function (_gmx_add_files_to_property PROPERTY)
    foreach (_file ${ARGN})
        if (IS_ABSOLUTE "${_file}")
            set_property(GLOBAL APPEND PROPERTY ${PROPERTY} ${_file})
        else()
            set_property(GLOBAL APPEND PROPERTY ${PROPERTY}
                         ${CMAKE_CURRENT_LIST_DIR}/${_file})
        endif()
    endforeach()
endfunction ()

function (gmx_add_libgromacs_sources)
    _gmx_add_files_to_property(GMX_LIBGROMACS_SOURCES ${ARGN})
endfunction ()

function (gmx_install_headers)
    if (NOT GMX_BUILD_MDRUN_ONLY)
        file(RELATIVE_PATH _dest ${PROJECT_SOURCE_DIR}/src ${CMAKE_CURRENT_LIST_DIR})
        install(FILES       ${ARGN}
                DESTINATION "${INCL_INSTALL_DIR}/${_dest}"
                COMPONENT   development)
    endif()
    _gmx_add_files_to_property(GMX_INSTALLED_HEADERS ${ARGN})
endfunction ()

function (gmx_write_installed_header_list)
    get_property(_list GLOBAL PROPERTY GMX_INSTALLED_HEADERS)
    string(REPLACE ";" "\n" _list "${_list}")
    # TODO: Make this only update the file timestamp if the contents actually change.
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/installed-headers.txt "${_list}")
endfunction()

add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(applied-forces)
add_subdirectory(listed-forces)
add_subdirectory(commandline)
add_subdirectory(domdec)
add_subdirectory(ewald)
add_subdirectory(fft)
add_subdirectory(gpu_utils)
add_subdirectory(hardware)
add_subdirectory(linearalgebra)
add_subdirectory(math)
add_subdirectory(mdrunutility)
add_subdirectory(mdtypes)
add_subdirectory(onlinehelp)
add_subdirectory(options)
add_subdirectory(pbcutil)
add_subdirectory(random)
add_subdirectory(tables)
add_subdirectory(taskassignment)
add_subdirectory(timing)
add_subdirectory(topology)
add_subdirectory(trajectory)
add_subdirectory(utility)
add_subdirectory(fileio)
add_subdirectory(swap)
add_subdirectory(essentialdynamics)
add_subdirectory(pulling)
add_subdirectory(simd)
add_subdirectory(imd)
if (NOT GMX_BUILD_MDRUN_ONLY)
    add_subdirectory(gmxana)
    add_subdirectory(gmxpreprocess)
    add_subdirectory(correlationfunctions)
    add_subdirectory(statistics)
    add_subdirectory(analysisdata)
    add_subdirectory(selection)
    add_subdirectory(trajectoryanalysis)
    add_subdirectory(energyanalysis)
    add_subdirectory(tools)
endif()
if (GMX_BUILD_FDA)
    add_subdirectory(fda)
endif()

get_property(PROPERTY_SOURCES GLOBAL PROPERTY GMX_LIBGROMACS_SOURCES)
list(APPEND LIBGROMACS_SOURCES ${GMXLIB_SOURCES} ${MDLIB_SOURCES} ${PROPERTY_SOURCES})

# This would be the standard way to include thread_mpi, but
# we want libgromacs to link the functions directly
#if(GMX_THREAD_MPI)
#    add_subdirectory(thread_mpi)
#endif()
#target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES} ${THREAD_MPI_LIB})

tmpi_get_source_list(THREAD_MPI_SOURCES ${CMAKE_SOURCE_DIR}/src/external/thread_mpi/src)
list(APPEND LIBGROMACS_SOURCES ${THREAD_MPI_SOURCES})

get_lmfit_properties(LMFIT_SOURCES LMFIT_LIBRARIES_TO_LINK LMFIT_INCLUDE_DIRECTORY LMFIT_INCLUDE_DIR_ORDER)
include_directories(${LMFIT_INCLUDE_DIR_ORDER} SYSTEM "${LMFIT_INCLUDE_DIRECTORY}")
list(APPEND LIBGROMACS_SOURCES ${LMFIT_SOURCES})

configure_file(version.h.cmakein version.h)
gmx_install_headers(
    analysisdata.h
    commandline.h
    options.h
    random.h
    selection.h
    trajectoryanalysis.h
    utility.h
    ${CMAKE_CURRENT_BINARY_DIR}/version.h
    )

# This code is here instead of utility/CMakeLists.txt, because CMake
# custom commands and source file properties can only be set in the directory
# that contains the target that uses them.
# TODO: Generate a header instead that can be included from baseversion.c.
# That probably simplifies things somewhat.
set(GENERATED_VERSION_FILE utility/baseversion-gen.c)
gmx_configure_version_file(
    utility/baseversion-gen.c.cmakein ${GENERATED_VERSION_FILE}
    REMOTE_HASH)
list(APPEND LIBGROMACS_SOURCES ${GENERATED_VERSION_FILE})

# set up CUDA compilation with clang
if (GMX_CLANG_CUDA)
    foreach (_file ${LIBGROMACS_SOURCES})
        get_filename_component(_ext ${_file} EXT)
        if (${_ext} STREQUAL ".cu")
            gmx_compile_cuda_file_with_clang(${_file})
        endif()
    endforeach()
endif()

if (GMX_USE_CUDA)
    # Work around FindCUDA that prevents using target_link_libraries()
    # with keywords otherwise...
    set(CUDA_LIBRARIES PRIVATE ${CUDA_LIBRARIES})
    if (NOT GMX_CLANG_CUDA)
        cuda_add_library(libgromacs ${LIBGROMACS_SOURCES})
    else()
        add_library(libgromacs ${LIBGROMACS_SOURCES})
    endif()
    target_link_libraries(libgromacs PRIVATE ${CUDA_CUFFT_LIBRARIES})
else()
    add_library(libgromacs ${LIBGROMACS_SOURCES})
endif()

# Recent versions of gcc and clang give warnings on scanner.cpp, which
# is a generated source file. These are awkward to suppress inline, so
# we do it in the compilation command (after testing that the compiler
# supports the suppressions).
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-unused HAS_NO_UNUSED)
if (HAS_NO_UNUSED)
    set(_scanner_cpp_compiler_flags "${_scanner_cpp_compiler_flags} -Wno-unused-parameter -Wno-unused-function")
endif()
set_source_files_properties(selection/scanner.cpp PROPERTIES COMPILE_FLAGS "${_scanner_cpp_compiler_flags}")

gmx_setup_tng_for_libgromacs()

target_link_libraries(libgromacs
                      PRIVATE
                      ${EXTRAE_LIBRARIES}
                      ${GMX_EXTRA_LIBRARIES}
                      ${GMX_COMMON_LIBRARIES}
                      ${FFT_LIBRARIES} ${LINEAR_ALGEBRA_LIBRARIES}
                      ${LMFIT_LIBRARIES_TO_LINK}
                      ${THREAD_LIB} ${GMX_SHARED_LINKER_FLAGS} ${OPENCL_LIBRARIES}
                      ${GMX_STDLIB_LIBRARIES}
                      PUBLIC
                      ${GMX_PUBLIC_LIBRARIES}
                      )
set_target_properties(libgromacs PROPERTIES
                      OUTPUT_NAME "gromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${LIBRARY_SOVERSION_MAJOR}
                      VERSION ${LIBRARY_VERSION}
                      COMPILE_FLAGS "${OpenMP_C_FLAGS}")

gmx_write_installed_header_list()

# Only install the library in mdrun-only mode if it is actually necessary
# for the binary
if (NOT GMX_BUILD_MDRUN_ONLY OR BUILD_SHARED_LIBS)
    install(TARGETS libgromacs
            EXPORT libgromacs
            LIBRARY DESTINATION ${LIB_INSTALL_DIR}
            RUNTIME DESTINATION ${BIN_INSTALL_DIR}
            ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
            COMPONENT libraries)
endif()

if (NOT GMX_BUILD_MDRUN_ONLY)
    include(InstallLibInfo.cmake)
endif()

# Technically, the user could want to do this for an OpenCL build
# using the CUDA runtime, but currently there's no reason to want to
# do that.
if (INSTALL_CUDART_LIB) #can be set manual by user
    if (GMX_USE_CUDA)
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
        message(WARNING "INSTALL_CUDART_LIB only makes sense when configuring for CUDA support")
    endif()
endif()

if(GMX_USE_OPENCL)
    set(OPENCL_KERNELS ${MDLIB_OPENCL_KERNELS})

    install(FILES ${OPENCL_KERNELS} DESTINATION
        ${OCL_INSTALL_DIR} COMPONENT libraries)
endif()

if (BUILD_TESTING)
    add_subdirectory(compat/tests)
endif()
