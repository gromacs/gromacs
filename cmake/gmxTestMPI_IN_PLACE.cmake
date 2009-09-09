# - Define macro to check if MPI_IN_PLACE exists
#
#  GMX_TEST_MPI_IN_PLACE(VARIABLE)
#
#  VARIABLE will be set to true if MPI_IN_PLACE exists
#

MACRO(GMX_TEST_MPI_IN_PLACE VARIABLE)
    MESSAGE(STATUS "Checking for MPI_IN_PLACE")
    # First check without any special flags
    TRY_COMPILE(MPI_IN_PLACE_COMPILE_OK ${CMAKE_BINARY_DIR}
                    "${CMAKE_SOURCE_DIR}/cmake/TestMPI_IN_PLACE.c"
                    COMPILE_DEFINITIONS )

    if(MPI_IN_PLACE_COMPILE_OK)
    MESSAGE(STATUS "Checking for MPI_IN_PLACE - yes")
        set(${VARIABLE} ${MPI_IN_PLACE_COMPILE_OK} 
                "Result of test for MPI_IN_PLACE")
    else(MPI_IN_PLACE_COMPILE_OK)
        MESSAGE(STATUS "Checking for MPI_IN_PLACE - no")
    endif(MPI_IN_PLACE_COMPILE_OK)
ENDMACRO(GMX_TEST_MPI_IN_PLACE VARIABLE)



