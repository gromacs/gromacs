set(LIBGROMACS_SOURCES)

add_subdirectory(legacyheaders)
add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(gmxpreprocess)
add_subdirectory(analysisdata)
add_subdirectory(commandline)
add_subdirectory(linearalgebra)
add_subdirectory(onlinehelp)
add_subdirectory(options)
add_subdirectory(selection)
add_subdirectory(trajectoryanalysis)
add_subdirectory(utility)

file(GLOB LIBGROMACS_HEADERS *.h)
install(FILES ${LIBGROMACS_HEADERS} DESTINATION ${INCL_INSTALL_DIR}/gromacs
        COMPONENT development)

list(APPEND LIBGROMACS_SOURCES ${GMXLIB_SOURCES} ${MDLIB_SOURCES})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/version.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.h
    DESTINATION ${INCL_INSTALL_DIR}/gromacs
    COMPONENT development)

# Add target that generates gitversion.c every time make is run
# if git version info is requested
# This code is here instead of utility/CMakeLists.txt because CMake
# ignores set_source_file_properties from subdirectories.
if (GMX_GIT_VERSION_INFO)
    set(GENERATED_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/utility/gitversion.c)
    add_custom_target(gmx_version ALL
            COMMAND ${CMAKE_COMMAND}
                -D GIT_EXECUTABLE="${GIT_EXECUTABLE}"
                -D GIT_VERSION="${GIT_VERSION}"
                -D PROJECT_VERSION="${PROJECT_VERSION}"
                -D PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                -D VERSION_C_CMAKEIN="${CMAKE_CURRENT_SOURCE_DIR}/utility/gitversion.c.cmakein"
                -D VERSION_C_OUT=${GENERATED_VERSION_FILE}
                -P ${CMAKE_SOURCE_DIR}/cmake/gmxGenerateVersionInfo.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/utility/gitversion.c.cmakein
            COMMENT "Generating git version information")
    set_source_files_properties(${GENERATED_VERSION_FILE}
                                PROPERTIES GENERATED true)
    list(APPEND LIBGROMACS_SOURCES ${GENERATED_VERSION_FILE})
endif()

# apply gcc 4.4.x bug workaround
if(GMX_USE_GCC44_BUG_WORKAROUND)
   include(gmxGCC44O3BugWorkaround)
   gmx_apply_gcc44_bug_workaround("gmxlib/bondfree.c")
   gmx_apply_gcc44_bug_workaround("mdlib/force.c")
   gmx_apply_gcc44_bug_workaround("mdlib/constr.c")
endif()

if(GMX_GPU)
    include_directories(${CUDA_TOOLKIT_INCLUDE})
endif()

add_library(libgromacs ${LIBGROMACS_SOURCES})
if (GMX_GIT_VERSION_INFO)
    add_dependencies(libgromacs gmx_version)
endif ()

if(GMX_BUILD_OWN_FFTW)
    # This dependency has to be made here rather than the CMakeLists.txt that
    # does the FFTW build, because of the order in which
    # add_subdirectory() calls are made in the top-level CMakeLists.txt; the
    # md library target does not necessarily exist yet. Also enabling and
    # disabling GMX_BUILD_OWN_FFTW changes dependencies correctly.
    add_dependencies(libgromacs gmxfftw)
endif()

target_link_libraries(libgromacs ${GMX_GPU_LIBRARIES}
                      ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES}
                      ${THREAD_LIB} ${GMX_SHARED_LINKER_FLAGS})
set_target_properties(libgromacs PROPERTIES
                      OUTPUT_NAME "gromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${SOVERSION}
                      COMPILE_FLAGS "${OpenMP_C_FLAGS}")

install(TARGETS libgromacs DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgromacs.pc.cmakein
               ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgromacs${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)

if (INSTALL_CUDART_LIB) #can be set manual by user
    if (GMX_OPENMM OR GMX_GPU)
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
        message(WARNING "INSTALL_CUDART_LIB only makes sense with GMX_OPENMM or GMX_GPU")
    endif()
endif ()
