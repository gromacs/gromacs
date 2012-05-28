include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# add target that generates version.c every time a make is run
# only do this if we generate the version
if(USE_VERSION_H)
    add_custom_target(gmx_version ALL
            COMMAND ${CMAKE_COMMAND} 
                -D GIT_EXECUTABLE="${GIT_EXECUTABLE}"
                -D GIT_VERSION="${GIT_VERSION}"
                -D PROJECT_VERSION="${PROJECT_VERSION}"
                -D PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                -D VERSION_C_CMAKEIN="${CMAKE_SOURCE_DIR}/src/gmxlib/version.c.cmakein"
                -D VERSION_C_OUT="${CMAKE_CURRENT_BINARY_DIR}/version.c"
                -P ${CMAKE_SOURCE_DIR}/cmake/gmxGenerateVersionInfo.cmake 
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/gmxlib 
            DEPENDS ${CMAKE_SOURCE_DIR}/src/gmxlib/version.c.cmakein
            COMMENT "Generating version information")
endif()

# The nonbonded directory contains subdirectories that are only
# conditionally built, so we cannot use a GLOB_RECURSE here.
file(GLOB GMXLIB_SOURCES *.c selection/*.c trajana/*.c statistics/*.c)

# Sources for nonbonded interactions subdirectory. This will eventually be moved
# to that subdirectory, but using a full library that is archived/unarchived
# multiple times will be slow (and non-parallel) for hundreds of files. 
# CMake 2.8.8 supports object libraries, but that is still a bit too new, so for 
# now we have to live with it here!
file(GLOB NONBONDED_SOURCES nonbonded/*.c nonbonded/nb_kernel_c/*.c nonbonded/nb_kernel_adress_c/*.c)
if(GMX_FORTRAN)
  if (GMX_DOUBLE)
    file(GLOB NB_KERNEL_FORTRAN_SOURCES nb_kernel_f77_double/*.[cf])
  else(GMX_DOUBLE)
    file(GLOB NB_KERNEL_FORTRAN_SOURCES nb_kernel_f77_single/*.[cf])
  endif(GMX_DOUBLE)
endif(GMX_FORTRAN)

if(GMX_POWER6)
  file(GLOB NB_KERNEL_FORTRAN_SOURCES nonbonded/nb_kernel_power6/*.[cF])
endif(GMX_POWER6)

if(GMX_BLUEGENE)
  file(GLOB NB_KERNEL_BLUEGENE_SOURCES nonbonded/nb_kernel_bluegene/*.c)
endif(GMX_BLUEGENE)

# End of sources from nonbonded interaction subdirectory



if(NOT GMX_USE_PLUGINS)
  list(REMOVE_ITEM GMXLIB_SOURCES vmdio.c vmddlopen.c)
endif()

# This source file is generated
file(GLOB VERSION_SOURCE version.c)

if(VERSION_SOURCE)
	list(REMOVE_ITEM GMXLIB_SOURCES ${VERSION_SOURCE})
endif(VERSION_SOURCE)

# add version.c to the list of sources and tell cmake that it is generated
if(USE_VERSION_H)
LIST(APPEND GMXLIB_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/version.c) # auto-generated
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/version.c 
                PROPERTIES GENERATED true)
endif()

if(NOT GMX_EXTERNAL_BLAS)
  file(GLOB BLAS_SOURCES gmx_blas/*.c)
endif(NOT GMX_EXTERNAL_BLAS)

if(NOT GMX_EXTERNAL_LAPACK)
  file(GLOB LAPACK_SOURCES gmx_lapack/*.c)
endif(NOT GMX_EXTERNAL_LAPACK)

# This would be the standard way to include thread_mpi, but we want libgmx
# to link the functions directly
#if(GMX_THREAD_MPI)
#    add_subdirectory(thread_mpi)
#endif(GMX_THREAD_MPI)
#target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES} ${THREAD_MPI_LIB})

# Files called xxx_test.c are test drivers with a main() function for module xxx.c,
# so they should not be included in the library
file(GLOB_RECURSE NOT_GMXLIB_SOURCES *_test.c *\#*)
list(REMOVE_ITEM GMXLIB_SOURCES ${NOT_GMXLIB_SOURCES})  
# Selection has test_ instead of _test.
file(GLOB SELECTION_TEST selection/test*)
list(REMOVE_ITEM GMXLIB_SOURCES ${SELECTION_TEST})

add_library(gmx ${GMXLIB_SOURCES} ${BLAS_SOURCES} ${LAPACK_SOURCES} ${THREAD_MPI_SRC} ${NONBONDED_SOURCES} ${NB_KERNEL_FORTRAN_SOURCES} ${NB_KERNEL_BLUEGENE_SOURCES})

target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES}  ${THREAD_LIB})
if(USE_VERSION_H)
	add_dependencies(gmx gmx_version) 
endif()
set_target_properties(gmx PROPERTIES OUTPUT_NAME "gmx${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

install(TARGETS gmx DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmx.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmx.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmx.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgmx${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)

