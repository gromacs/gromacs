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

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
    message(FATAL_ERROR "CMake called with the wrong directory; please use the top GROMACS directory (${CMAKE_CURRENT_SOURCE_DIR}/../../)")
endif()

set(LIBGROMACS_SOURCES)

if (GMX_CLANG_CUDA)
    include(gmxClangCudaUtils)
endif()

set_property(GLOBAL PROPERTY GMX_LIBGROMACS_SOURCES)
set_property(GLOBAL PROPERTY CUDA_SOURCES)
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

# Permit the configuration to disable compiling the many nbnxm kernels
# and others involved in force calculations. Currently only
# short-ranged and bonded kernels are disabled this way, but in future
# others may be appropriate. Thus the cmake option is not specific to
# nbnxm module.
option(GMX_USE_SIMD_KERNELS "Whether to compile NBNXM and other SIMD kernels" ON)
mark_as_advanced(GMX_USE_SIMD_KERNELS)

# pulling is using the muparser target as well, so detect it early
gmx_manage_muparser()

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
add_subdirectory(simd)
add_subdirectory(imd)
add_subdirectory(compat)
add_subdirectory(mimic)
add_subdirectory(modularsimulator)
add_subdirectory(gmxana)
add_subdirectory(gmxpreprocess)
add_subdirectory(correlationfunctions)
add_subdirectory(statistics)
add_subdirectory(analysisdata)
add_subdirectory(coordinateio)
add_subdirectory(trajectoryanalysis)
add_subdirectory(energyanalysis)
add_subdirectory(tools)

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
target_compile_definitions(thread_mpi PRIVATE TMPI_EXPORTS TMPI_USE_VISIBILITY)
if(CYGWIN)
    # Needs POSIX-isms for strdup, not just std-isms
    target_compile_definitions(thread_mpi PRIVATE _POSIX_C_SOURCE=200809L)
endif()
gmx_target_compile_options(thread_mpi)
if (WIN32)
    gmx_target_warning_suppression(thread_mpi /wd4996 HAS_NO_MSVC_UNSAFE_FUNCTION)
endif()
list(APPEND libgromacs_object_library_dependencies thread_mpi)

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
if (GMX_GPU_CUDA)
    get_property(CUDA_SOURCES GLOBAL PROPERTY CUDA_SOURCES)
    set_source_files_properties(${CUDA_SOURCES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
endif()

if (GMX_GPU_CUDA)
    if(GMX_USE_cuFFTMp)
        # needed as we need to include cufftmp include path before CUDA include path
        cuda_include_directories(${cuFFTMp_INCLUDE_DIR})
    endif()
    # Work around FindCUDA that prevents using target_link_libraries()
    # with keywords otherwise...
    set(CUDA_LIBRARIES PRIVATE ${CUDA_LIBRARIES})

    if (NOT GMX_CLANG_CUDA AND NOT GMX_NVSHMEM)
        gmx_cuda_add_library(libgromacs ${LIBGROMACS_SOURCES})
    else()
        add_library(libgromacs ${LIBGROMACS_SOURCES})
    endif()

    #we cannot link with both cufftmp and cufft
    if(GMX_USE_cuFFTMp)
        target_link_libraries(libgromacs PRIVATE ${cuFFTMp_LIBRARY})
    else()
        target_link_libraries(libgromacs PRIVATE ${CUDA_CUFFT_LIBRARIES})
    endif()
    if (GMX_NVSHMEM)
        set_target_properties(libgromacs PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        # -fvisibility=hidden is required to be passed during device link phase to hide the
        # libnvshmem_device.a symbols in libgromacs in order to allow application/tests to link
        # libnvshmem_device.a for their nvshmem kernel as well while linking with libgromacs.
        # not setting this causes runtime failure.
        target_link_options(libgromacs PRIVATE $<DEVICE_LINK:-fvisibility=hidden>)
        set_target_properties(libgromacs PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
        target_link_libraries(libgromacs PRIVATE nvshmem_host_lib nvshmem_device_lib)
        # NVSHMEM only supports SM 60+ so we make use of filtered GMX_NVSHMEM_LINK_ARCHS
        set_target_properties(libgromacs PROPERTIES CUDA_ARCHITECTURES "${GMX_NVSHMEM_LINK_ARCHS}")
    endif()
    if (GMX_USE_cuFFTMp OR GMX_NVSHMEM)
        # We need to PUBLIC link to the stub libraries nvml/cuda to WAR an issue
        # with `bfd` linker which expects you always have dependencies on the link line.
        # Also leaking path via the -rpath does not work with a stub library as the SONAME
        # present inside these stub library (used by runtime) is *.so.1 instead of *.so as these
        # stub library are meant only for compile time and the runtime version is a driver provided
        # version which is not present on the system where driver is not installed.
        # At some future release cmake will provide a workaround for this under the hood at which
        # point we can use PRIVATE linking for these libs until then we stick with PUBLIC linking.
        target_link_libraries(libgromacs PUBLIC ${GMX_CUDA_DRV_LIB} ${GMX_NVIDIA_ML_LIB})
    endif()
else()
    add_library(libgromacs ${LIBGROMACS_SOURCES})
endif()

if (TARGET Heffte::Heffte)
    target_link_libraries(libgromacs PRIVATE Heffte::Heffte)
endif()

if(GMX_GPU_FFT_ROCFFT)
    target_link_libraries(libgromacs PRIVATE roc::rocfft)
endif()

target_link_libraries(libgromacs PRIVATE $<BUILD_INTERFACE:common>)
# As long as the libgromacs target has source files that reference headers from
# modules that don't provide CMake targets, libgromacs needs to use `src/`
# amongst its include directories (to support `#include "gromacs/module/header.h"`).
add_library(legacy_modules INTERFACE)
target_include_directories(legacy_modules INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)
target_link_libraries(libgromacs PRIVATE $<BUILD_INTERFACE:legacy_modules>)

# Add these contents first because linking their tests can take a lot
# of time, so we want lots of parallel work still available after
# linking starts.
add_subdirectory(fileio)
add_subdirectory(selection)

# Handle the object libraries that contain the source file
# dependencies that need special handling because they are generated
# or external code.
foreach(object_library ${libgromacs_object_library_dependencies})
    if (BUILD_SHARED_LIBS)
        set_target_properties(${object_library} PROPERTIES POSITION_INDEPENDENT_CODE true)
    endif()
    target_include_directories(${object_library} SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
    target_compile_definitions(${object_library} PRIVATE TMPI_USE_VISIBILITY)
    target_link_libraries(${object_library} PRIVATE common)

    # Add the sources from the object libraries to the main library.
    target_sources(libgromacs PRIVATE $<TARGET_OBJECTS:${object_library}>)
endforeach()
gmx_target_compile_options(libgromacs)
target_compile_definitions(libgromacs PRIVATE HAVE_CONFIG_H
                                              TMPI_EXPORTS
                                      PUBLIC TMPI_USE_VISIBILITY)

target_include_directories(libgromacs SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)

if (GMX_GPU_FFT_CLFFT)
    if (NOT GMX_GPU_OPENCL)
        message(FATAL_ERROR "clFFT is only supported in OpenCL builds")
    endif()
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
build of clFFT in GROMACS 2019 and newer. Either change compiler, try
installing a clFFT package, use VkFFT by setting -DGMX_GPU_FFT_LIBRARY=VkFFT,
or use GROMACS 2018.8 release (no longer maintained).")
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

if(SIMD_AVX_512_CXX_SUPPORTED AND NOT ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_512_KNL"))
    # Since we might be overriding -march=core-avx2, add a flag so we don't warn for this specific file.
    # On KNL this can cause illegal instruction because the compiler might use non KNL AVX instructions
    # with the SIMD_AVX_512_CXX_FLAGS flags.
    set_source_files_properties(hardware/identifyavx512fmaunits.cpp PROPERTIES COMPILE_FLAGS "${SIMD_AVX_512_CXX_FLAGS} ${CXX_NO_UNUSED_OPTION_WARNING_FLAGS}")
endif()

# Do any special handling needed for .cpp files that use
# CUDA runtime headers
if (GMX_GPU_CUDA AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    foreach(_compile_flag ${NVCC_CLANG_SUPPRESSIONS_CXXFLAGS})
        set(GMX_CUDA_CLANG_FLAGS "${GMX_CUDA_CLANG_FLAGS} ${_compile_flag}")
    endforeach()
    if (GMX_CLANG_CUDA)
        foreach (_file ${LIBGROMACS_SOURCES})
            get_filename_component(_ext ${_file} EXT)
            get_source_file_property(_cuda_source_format ${_file} CUDA_SOURCE_PROPERTY_FORMAT)
            if ("${_ext}" STREQUAL ".cu" OR _cuda_source_format)
                gmx_compile_cuda_file_with_clang(${_file})
            endif()
        endforeach()
    else()
        get_property(CUDA_SOURCES GLOBAL PROPERTY CUDA_SOURCES)
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES COMPILE_FLAGS ${GMX_CUDA_CLANG_FLAGS})
    endif()
endif()

# Only add the -fsycl flag to sources that really need it
if (GMX_GPU_SYCL)
    get_property(SYCL_SOURCES GLOBAL PROPERTY SYCL_SOURCES)
    add_sycl_to_target(TARGET libgromacs SOURCES ${SYCL_SOURCES})
endif()

gmx_setup_tng_for_libgromacs()

# We apply the SYCL flag explicitly just for libgromacs, since bugs in the beta versions of
# icpx/dpcpp leads to crashes if we try to link an library without any SYCL code with the
# -fsycl flag enabled. Once that bug is fixed, we should change it to simply add
# SYCL_CXX_FLAGS to GMX_SHARED_LINKER_FLAGS.
target_link_libraries(libgromacs
                      PRIVATE
                      ${EXTRAE_LIBRARIES}
                      ${GMX_EXTRA_LIBRARIES}
                      ${GMX_COMMON_LIBRARIES}
                      ${FFT_LIBRARIES} ${LINEAR_ALGEBRA_LIBRARIES}
                      ${THREAD_LIB} ${GMX_SHARED_LINKER_FLAGS}
                      ${SYCL_CXX_FLAGS}
                      ${OpenCL_LIBRARIES}
                      $<$<PLATFORM_ID:SunOS>:socket>
                      PUBLIC
                      ${GMX_PUBLIC_LIBRARIES}
                      )
target_link_libraries(libgromacs PRIVATE legacy_api)
# Dependencies from libgromacs to the modules are set up here, but
# once the add_subdirectory() commands are re-ordered then
# responsibility for setting this up will move to the respective
# modules.
target_link_libraries(libgromacs PRIVATE
                      $<BUILD_INTERFACE:analysisdata>
                      $<BUILD_INTERFACE:applied_forces>
                      $<BUILD_INTERFACE:commandline>
                      $<BUILD_INTERFACE:compat>
                      $<BUILD_INTERFACE:coordinateio>
                      $<BUILD_INTERFACE:correlationfunctions>
                      $<BUILD_INTERFACE:domdec>
                      $<BUILD_INTERFACE:energyanalysis>
                      $<BUILD_INTERFACE:essentialdynamics>
                      $<BUILD_INTERFACE:ewald>
                      $<BUILD_INTERFACE:fft>
                      $<BUILD_INTERFACE:fileio>
                      $<BUILD_INTERFACE:gmxana>
                      $<BUILD_INTERFACE:gmxlib>
                      $<BUILD_INTERFACE:gmxpreprocess>
                      $<BUILD_INTERFACE:gpu_utils>
                      $<BUILD_INTERFACE:hardware>
                      $<BUILD_INTERFACE:imd>
                      $<BUILD_INTERFACE:linearalgebra>
                      $<BUILD_INTERFACE:listed_forces>
                      $<BUILD_INTERFACE:math>
                      $<BUILD_INTERFACE:mdlib>
                      $<BUILD_INTERFACE:mdrun>
                      $<BUILD_INTERFACE:mdrunutility>
                      $<BUILD_INTERFACE:mdspan>
                      $<BUILD_INTERFACE:mdtypes>
                      $<BUILD_INTERFACE:mimic>
                      $<BUILD_INTERFACE:modularsimulator>
                      $<BUILD_INTERFACE:nbnxm>
                      $<BUILD_INTERFACE:onlinehelp>
                      $<BUILD_INTERFACE:options>
                      $<BUILD_INTERFACE:pbcutil>
                      $<BUILD_INTERFACE:pulling>
                      $<BUILD_INTERFACE:random>
                      $<BUILD_INTERFACE:restraint>
                      $<BUILD_INTERFACE:selection>
                      $<BUILD_INTERFACE:simd>
                      $<BUILD_INTERFACE:statistics>
                      $<BUILD_INTERFACE:swap>
                      $<BUILD_INTERFACE:tables>
                      $<BUILD_INTERFACE:taskassignment>
                      $<BUILD_INTERFACE:timing>
                      $<BUILD_INTERFACE:tools>
                      $<BUILD_INTERFACE:topology>
                      $<BUILD_INTERFACE:trajectory>
                      $<BUILD_INTERFACE:trajectoryanalysis>
                      $<BUILD_INTERFACE:utility>
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
target_link_libraries(libgromacs PRIVATE muparser::muparser)
gmx_manage_colvars()
target_link_libraries(libgromacs PRIVATE colvars)

# Make sure we fix "everything" found by compilers that support that
gmx_warn_on_everything(libgromacs)
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
   target_compile_options(libgromacs PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/analyze /analyze:stacksize 70000
     #Control flow warnings are disabled because the commond line output is insufficient. There is no tool
     #to convert the xml report to e.g. HTML and even in Visual Studio the viewer doesn't work with cmake support.
     /wd6001  #uninitialized memory
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
     /wd28020 #expression involving parameter is not true
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
   # Temporarily exclude Colvars headers to silence warnings from code that is not built with GROMACS
   # TODO: remove this when clang-tidy >= 14 is used
   set_target_properties(libgromacs PROPERTIES CXX_CLANG_TIDY
       "${CLANG_TIDY_EXE};-warnings-as-errors=*;-header-filter=^(?!.*src/external/colvars).*")
endif()

# TODO: Stop installing libgromacs. Possibly allow installation during deprecation period with GMX_INSTALL_LEGACY_API.
if (BUILD_SHARED_LIBS)
    install(TARGETS libgromacs
            EXPORT libgromacs
            LIBRARY
                DESTINATION ${CMAKE_INSTALL_LIBDIR}
                COMPONENT libraries
            RUNTIME
                DESTINATION ${CMAKE_INSTALL_BINDIR}
                COMPONENT libraries
            ARCHIVE
                DESTINATION ${CMAKE_INSTALL_LIBDIR}
                COMPONENT libraries
            INCLUDES DESTINATION include)
    target_compile_definitions(libgromacs PUBLIC $<INSTALL_INTERFACE:GMX_DOUBLE=${GMX_DOUBLE_VALUE}>)
    # legacy headers use c++17 features, so consumer codes need to use that standard, too
    if(GMX_INSTALL_LEGACY_API)
        target_compile_features(libgromacs INTERFACE cxx_std_${CMAKE_CXX_STANDARD})
    endif()
    # only install pkg-config and cmake target files when library is installed
    include(InstallLibInfo.cmake)
endif()
add_library(Gromacs::libgromacs ALIAS libgromacs)


# Technically, the user could want to do this for an OpenCL build
# using the CUDA runtime, but currently there's no reason to want to
# do that.
if (INSTALL_CUDART_LIB) #can be set manual by user
    if (GMX_GPU_CUDA)
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

if(GMX_GPU_OPENCL)
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
        ewald/pme_gpu_calculate_splines.clh
        ewald/pme_program.cl
        ewald/pme_gpu_types.h
        )
    install(FILES ${OPENCL_INSTALLED_FILES}
        DESTINATION ${GMX_INSTALL_OCLDIR}/gromacs/ewald
        COMPONENT libraries)
endif()

include(GenerateExportHeader)
generate_export_header(libgromacs EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/api/legacy/include/gromacs/libgromacs_export.h")
