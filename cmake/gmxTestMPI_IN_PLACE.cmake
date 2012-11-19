# - Define macro to check if MPI_IN_PLACE exists
#
#  GMX_TEST_MPI_IN_PLACE(VARIABLE)
#
#  VARIABLE will be set to true if MPI_IN_PLACE exists
#

include(CheckCSourceCompiles)
MACRO(GMX_TEST_MPI_IN_PLACE VARIABLE)
  if(NOT DEFINED MPI_IN_PLACE_COMPILE_OK)
    MESSAGE(STATUS "Checking for MPI_IN_PLACE")

    set(CMAKE_REQUIRED_DEFINITIONS ${MPI_COMPILE_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
    check_c_source_compiles(
      "#include <mpi.h>
int main(void) {
  void* buf;
  MPI_Allreduce(MPI_IN_PLACE, buf, 10, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}" MPI_IN_PLACE_COMPILE_OK)

    if(MPI_IN_PLACE_COMPILE_OK)
        MESSAGE(STATUS "Checking for MPI_IN_PLACE - yes")
    else(MPI_IN_PLACE_COMPILE_OK)
        MESSAGE(STATUS "Checking for MPI_IN_PLACE - no")
    endif(MPI_IN_PLACE_COMPILE_OK)
    set(MPI_IN_PLACE_COMPILE_OK "${MPI_IN_PLACE_COMPILE_OK}" CACHE INTERNAL "Result of mpi_in_place check")
    set(CMAKE_REQUIRED_DEFINITIONS)
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()
  if (MPI_IN_PLACE_COMPILE_OK)
    set(${VARIABLE} ${MPI_IN_PLACE_COMPILE_OK} 
      "Result of test for MPI_IN_PLACE")
  endif()
ENDMACRO(GMX_TEST_MPI_IN_PLACE VARIABLE)



