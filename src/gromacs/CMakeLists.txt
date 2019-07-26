#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
set_property(GLOBAL PROPERTY GMX_LIBGROMACS_GPU_IMPL_SOURCES)
set_property(GLOBAL PROPERTY GMX_INSTALLED_HEADERS)
set_property(GLOBAL PROPERTY GMX_AVX_512_SOURCE)

set(libgromacs_object_library_dependencies "")
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

# TODO Reconsider this, as the CUDA driver API is probably a simpler
# approach, at least for the build system. See Redmine #2530
function (gmx_compile_cpp_as_cuda)
    _gmx_add_files_to_property(GMX_LIBGROMACS_GPU_IMPL_SOURCES ${ARGN})
endfunction ()

# Add these contents first because linking their tests can take a lot
# of time, so we want lots of parallel work still available after
# linking starts.
add_subdirectory(utility)
# Add normal contents
add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(applied_forces)
add_subdirectory(listed_forces)
add_subdirectory(nbnxm)
add_subdirectory(commandline)
add_subdirectory(domdec)
add_subdirectory(ewald)
add_subdirectory(fft)
add_subdirectory(gpu_utils)
add_subdirectory(hardware)
add_subdirectory(linearalgebra)
add_subdirectory(math)
add_subdirectory(mdrun)
add_subdirectory(mdrunutility)
add_subdirectory(mdspan)
add_subdirectory(mdtypes)
add_subdirectory(onlinehelp)
add_subdirectory(options)
add_subdirectory(pbcutil)
add_subdirectory(random)
add_subdirectory(restraint)
add_subdirectory(tables)
add_subdirectory(taskassignment)
add_subdirectory(timing)
add_subdirectory(topology)
add_subdirectory(trajectory)
add_subdirectory(swap)
add_subdirectory(essentialdynamics)
add_subdirectory(pulling)
add_subdirectory(awh)
add_subdirectory(simd)
add_subdirectory(imd)
add_subdirectory(compat)
add_subdirectory(mimic)
add_subdirectory(modularsimulator)
if (NOT GMX_BUILD_MDRUN_ONLY)
    add_subdirectory(gmxana)
    add_subdirectory(gmxpreprocess)
    add_subdirectory(correlationfunctions)
    add_subdirectory(statistics)
    add_subdirectory(analysisdata)
    add_subdirectory(coordinateio)
    add_subdirectory(trajectoryanalysis)
    add_subdirectory(energyanalysis)
    add_subdirectory(tools)
endif()

get_property(PROPERTY_SOURCES GLOBAL PROPERTY GMX_LIBGROMACS_SOURCES)
list(APPEND LIBGROMACS_SOURCES ${GMXLIB_SOURCES} ${MDLIB_SOURCES} ${PROPERTY_SOURCES})

# This would be the standard way to include thread_mpi, but
# we want libgromacs to link the functions directly
#if(GMX_THREAD_MPI)
#    add_subdirectory(thread_mpi)
#endif()
#target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES} ${THREAD_MPI_LIB})
tmpi_get_source_list(THREAD_MPI_SOURCES ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/src)
add_library(thread_mpi OBJECT ${THREAD_MPI_SOURCES})
target_compile_definitions(thread_mpi PRIVATE HAVE_CONFIG_H)
gmx_target_compile_options(thread_mpi)
if (WIN32)
    gmx_target_warning_suppression(thread_mpi /wd4996 HAS_NO_MSVC_UNSAFE_FUNCTION)
endif()
list(APPEND libgromacs_object_library_dependencies thread_mpi)

configure_file(version.h.cmakein version.h)
if(GMX_INSTALL_LEGACY_API)
  install(FILES
          ${CMAKE_CURRENT_BINARY_DIR}/version.h
          DESTINATION include/gromacs)
endif()

# This code is here instead of utility/CMakeLists.txt, because CMake
# custom commands and source file properties can only be set in the directory
# that contains the target that uses them.
# TODO: Generate a header instead that can be included from baseversion.cpp.
# That probably simplifies things somewhat.
set(GENERATED_VERSION_FILE utility/baseversion-gen.cpp)
gmx_configure_version_file(
    utility/baseversion-gen.cpp.cmakein ${GENERATED_VERSION_FILE}
    REMOTE_HASH
    EXTRA_VARS
        GMX_SOURCE_DOI
        GMX_RELEASE_HASH
        GMX_SOURCE_HASH
        )
list(APPEND LIBGROMACS_SOURCES ${GENERATED_VERSION_FILE})

# Mark some shared GPU implementation files to compile with CUDA if needed
if (GMX_USE_CUDA)
    get_property(LIBGROMACS_GPU_IMPL_SOURCES GLOBAL PROPERTY GMX_LIBGROMACS_GPU_IMPL_SOURCES)
    set_source_files_properties(${LIBGROMACS_GPU_IMPL_SOURCES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
endif()

# set up CUDA compilation with clang
if (GMX_CLANG_CUDA)
    foreach (_file ${LIBGROMACS_SOURCES})
        get_filename_component(_ext ${_file} EXT)
        get_source_file_property(_cuda_source_format ${_file} CUDA_SOURCE_PROPERTY_FORMAT)
        if ("${_ext}" STREQUAL ".cu" OR _cuda_source_format)
            gmx_compile_cuda_file_with_clang(${_file})
        endif()
    endforeach()
endif()

if (GMX_USE_CUDA)
    # Work around FindCUDA that prevents using target_link_libraries()
    # with keywords otherwise...
    set(CUDA_LIBRARIES PRIVATE ${CUDA_LIBRARIES})
    if (NOT GMX_CLANG_CUDA)
        gmx_cuda_add_library(libgromacs ${LIBGROMACS_SOURCES})
    else()
        add_library(libgromacs ${LIBGROMACS_SOURCES})
    endif()
    target_link_libraries(libgromacs PRIVATE ${CUDA_CUFFT_LIBRARIES})
else()
    add_library(libgromacs ${LIBGROMACS_SOURCES})
endif()

# Add these contents first because linking their tests can take a lot
# of time, so we want lots of parallel work still available after
# linking starts.
add_subdirectory(fileio)
add_subdirectory(selection)

# Suppress a warning about our abuse of t_inputrec
gmx_source_file_warning_suppression(mdtypes/inputrec.cpp -Wno-class-memaccess HAS_NO_CLASS_MEMACCESS)

# Handle the object libraries that contain the source file
# dependencies that need special handling because they are generated
# or external code.
foreach(object_library ${libgromacs_object_library_dependencies})
    if (BUILD_SHARED_LIBS)
        set_target_properties(${object_library} PROPERTIES POSITION_INDEPENDENT_CODE true)
    endif()
    target_include_directories(${object_library} SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)

    # Add the sources from the object libraries to the main library.
    target_sources(libgromacs PRIVATE $<TARGET_OBJECTS:${object_library}>)
endforeach()
gmx_target_compile_options(libgromacs)
target_compile_definitions(libgromacs PRIVATE HAVE_CONFIG_H)
target_include_directories(libgromacs SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)

if (GMX_USE_OPENCL)
    option(GMX_EXTERNAL_CLFFT "True if an external clFFT is required to be used" FALSE)
    mark_as_advanced(GMX_EXTERNAL_CLFFT)

    # Default to using clFFT found on the system
    # switch to quiet at the second run.
    if (DEFINED clFFT_LIBRARY)
        set (clFFT_FIND_QUIETLY TRUE)
    endif()
    find_package(clFFT)
    if (NOT clFFT_FOUND)
        if (GMX_EXTERNAL_CLFFT)
            message(FATAL_ERROR "Did not find required external clFFT library, consider setting clFFT_ROOT_DIR")
        endif()

        if(MSVC)
            message(FATAL_ERROR
"An OpenCL build was requested with Visual Studio compiler, but GROMACS
requires clFFT, which was not found on your system. GROMACS does bundle
clFFT to help with building for OpenCL, but that clFFT has not yet been
ported to the more recent versions of that compiler that GROMACS itself
requires. Thus for now, OpenCL is not available with MSVC and the internal
build of clFFT in GROMACS 2019. Either change compiler, try installing
a clFFT package, or use the latest GROMACS 2018 point release.")
        endif()

        # Fall back on the internal version
        set (_clFFT_dir ../external/clFFT/src)
        add_subdirectory(${_clFFT_dir} clFFT-build)
        target_sources(libgromacs PRIVATE
            $<TARGET_OBJECTS:clFFT>
        )
        target_include_directories(libgromacs SYSTEM PRIVATE ${_clFFT_dir}/include)
        # Use the magic variable for how to link any library needed for
        # dlopen, etc.  which is -ldl where needed, and empty otherwise
        # (e.g. Windows, BSD, Mac).
        target_link_libraries(libgromacs PRIVATE "${CMAKE_DL_LIBS}")
    else()
        target_link_libraries(libgromacs PRIVATE clFFT)
    endif()
endif()

# Permit GROMACS code to include externally developed headers, such as
# the functionality from the nonstd project that we use for
# gmx::compat::optional. These are included as system headers so that
# no warnings are issued from them.
#
# TODO Perhaps generalize this for all headers from src/external
target_include_directories(libgromacs SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)

# Temporary fix to allow external access to restraintpotentail
target_include_directories(libgromacs PUBLIC
                           $<INSTALL_INTERFACE:include>
                           )

if(SIMD_AVX_512_CXX_SUPPORTED AND NOT ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_512_KNL"))
    # Since we might be overriding -march=core-avx2, add a flag so we don't warn for this specific file.
    # On KNL this can cause illegal instruction because the compiler might use non KNL AVX instructions
    # with the SIMD_AVX_512_CXX_FLAGS flags.
    set_source_files_properties(hardware/identifyavx512fmaunits.cpp PROPERTIES COMPILE_FLAGS "${SIMD_AVX_512_CXX_FLAGS} ${CXX_NO_UNUSED_OPTION_WARNING_FLAGS}")
endif()

gmx_setup_tng_for_libgromacs()

target_link_libraries(libgromacs
                      PRIVATE
                      ${EXTRAE_LIBRARIES}
                      ${GMX_EXTRA_LIBRARIES}
                      ${GMX_COMMON_LIBRARIES}
                      ${FFT_LIBRARIES} ${LINEAR_ALGEBRA_LIBRARIES}
                      ${THREAD_LIB} ${GMX_SHARED_LINKER_FLAGS}
                      ${OpenCL_LIBRARIES}
                      $<$<PLATFORM_ID:SunOS>:socket>
                      PUBLIC
                      ${GMX_PUBLIC_LIBRARIES}
                      )
if (GMX_OPENMP)
    target_link_libraries(libgromacs PUBLIC OpenMP::OpenMP_CXX)
endif()
set_target_properties(libgromacs PROPERTIES
                      OUTPUT_NAME "gromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${LIBRARY_SOVERSION_MAJOR}
                      VERSION ${LIBRARY_VERSION}
                      )

gmx_manage_lmfit()
target_link_libraries(libgromacs PRIVATE lmfit)

# Fix everything found by the latest version of clang that we use in
# Jenkins testing. This should be updated when we update the latest
# tested version of clang.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION MATCHES "^7\.0")
   target_compile_options(libgromacs PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Weverything ${IGNORED_CLANG_ALL_WARNINGS}>)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
   target_compile_options(libgromacs PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/analyze /analyze:stacksize 70000
     #Control flow warnings are disabled because the commond line output is insufficient. There is no tool
     #to convert the xml report to e.g. HTML and even in Visual Studio the viewer doesn't work with cmake support.
     /wd6001  #unitialized memory
     /wd6011  #derefencing NULL
     /wd6053  #prior call not zero-terminate
     /wd6054  #might not be zero-terminated
     /wd6385  #reading invalid data
     /wd6386  #buffer overrun
     /wd6387  #could be '0'
     /wd28199 #uninitialized memory
     # For compile time constant (e.g. templates) the following warnings have flase postives
     /wd6239  #(<non-zero> && <expr>)
     /wd6240  #(<expr> && <non-zero>)
     /wd6294  #Ill-defined for-loop
     /wd6326  #comparison of constant with other constant
     /wd28020 #expression involving paramter is not true
     # Misc
     /wd6330  #incorrect type to function (warns for char (instead of unsigned) for isspace/isalpha/isdigit/..))
     /wd6993  #OpenMP ignored
     #TODO
     /wd6031  #return value ignored (important - mostly warnigns about sscanf)
     /wd6244  #hides declaration (known issue - we ingore similar warnings for other compilers)
     /wd6246  #hides declaration
     >
   )
endif()

if (GMX_CLANG_TIDY)
   set_target_properties(libgromacs PROPERTIES CXX_CLANG_TIDY
       "${CLANG_TIDY_EXE};-warnings-as-errors=*")
endif()

# Only install the library in mdrun-only mode if it is actually necessary
# for the binary
if (NOT GMX_BUILD_MDRUN_ONLY OR BUILD_SHARED_LIBS)
    install(TARGETS libgromacs
            EXPORT libgromacs
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION include
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
                    ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)
            endif()
        endforeach()
    else()
        message(WARNING "INSTALL_CUDART_LIB only makes sense when configuring for CUDA support")
    endif()
endif()

if(GMX_USE_OPENCL)
    # Install the utility headers
    file(GLOB OPENCL_INSTALLED_FILES
        gpu_utils/vectype_ops.clh
        gpu_utils/device_utils.clh
        )
    install(FILES ${OPENCL_INSTALLED_FILES}
        DESTINATION ${GMX_INSTALL_OCLDIR}/gromacs/gpu_utils
        COMPONENT libraries)
    file(GLOB OPENCL_INSTALLED_FILES
        pbcutil/ishift.h
        )
    install(FILES ${OPENCL_INSTALLED_FILES}
        DESTINATION ${GMX_INSTALL_OCLDIR}/gromacs/pbcutil
        COMPONENT libraries)

    # Install the NBNXM source and headers
    file(GLOB OPENCL_INSTALLED_FILES
        nbnxm/constants.h
        )
    install(FILES ${OPENCL_INSTALLED_FILES}
        DESTINATION ${GMX_INSTALL_OCLDIR}/gromacs/nbnxm
        COMPONENT libraries)
    file(GLOB OPENCL_INSTALLED_FILES
        nbnxm/opencl/nbnxm_ocl_kernels.cl
        nbnxm/opencl/nbnxm_ocl_kernel.clh
        nbnxm/opencl/nbnxm_ocl_kernel_pruneonly.clh
        nbnxm/opencl/nbnxm_ocl_kernels.clh
        nbnxm/opencl/nbnxm_ocl_kernels_fastgen.clh
        nbnxm/opencl/nbnxm_ocl_kernels_fastgen_add_twincut.clh
        nbnxm/opencl/nbnxm_ocl_kernel_utils.clh
        nbnxm/opencl/nbnxm_ocl_consts.h
        )
    install(FILES ${OPENCL_INSTALLED_FILES}
        DESTINATION ${GMX_INSTALL_OCLDIR}/gromacs/nbnxm/opencl
        COMPONENT libraries)

    # Install the PME source and headers
    file(GLOB OPENCL_INSTALLED_FILES
        ewald/pme_spread.clh
        ewald/pme_solve.clh
        ewald/pme_gather.clh
        ewald/pme_gpu_utils.clh
        ewald/pme_program.cl
        ewald/pme_gpu_types.h
        )
    install(FILES ${OPENCL_INSTALLED_FILES}
        DESTINATION ${GMX_INSTALL_OCLDIR}/gromacs/ewald
        COMPONENT libraries)
endif()
