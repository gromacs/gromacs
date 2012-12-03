include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(nonbonded)

# The nonbonded directory contains subdirectories that are only
# conditionally built, so we cannot use a GLOB_RECURSE here.
file(GLOB GMXLIB_SOURCES *.c *.cpp statistics/*.c)

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

# gpu utils + cuda tools module
if(GMX_GPU)
    # The log file output queries Cuda if GPU support is enabled
    add_subdirectory(cuda_tools)
    add_subdirectory(gpu_utils)   
    set(GMX_GPU_LIBRARIES ${GMX_GPU_LIBRARIES} gpu_utils cuda_tools PARENT_SCOPE)
endif()

if(GMX_USE_PLUGINS)
  set(GMXLIB_SOURCES ${GMXLIB_SOURCES} ${CMAKE_SOURCE_DIR}/src/external/vmd_molfile/vmddlopen.c)
else()
  list(REMOVE_ITEM GMXLIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/vmdio.c)
endif()

# An ugly hack to get absolute paths...
file(GLOB THREAD_MPI_SOURCES ${THREAD_MPI_SRC})

set(GMXLIB_SOURCES ${GMXLIB_SOURCES} ${THREAD_MPI_SOURCES} ${NONBONDED_SOURCES}
    PARENT_SCOPE)
