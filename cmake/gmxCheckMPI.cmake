# Test for and warn about unsuitable MPI versions
exec_program(ompi_info
  ARGS -v ompi full
  OUTPUT_VARIABLE OPENMPI_TYPE
  RETURN_VALUE OPENMPI_EXEC_RETURN)
if(OPENMPI_EXEC_RETURN EQUAL 0)
  string(REGEX REPLACE ".*Open MPI: \([0-9]+\\.[0-9]*\\.?[0-9]*\).*" "\\1" OPENMPI_VERSION ${OPENMPI_TYPE})
  if(OPENMPI_VERSION VERSION_LESS "1.4.1")
    MESSAGE(WARNING "
             CMake found OpenMPI version ${OPENMPI_VERSION} on your system.
             There are known problems with GROMACS and OpenMPI version < 1.4.1.
             Please consider updating your OpenMPI if your MPI wrapper compilers
             are using the above OpenMPI version.")
  endif()
  unset(OPENMPI_VERSION)
else()
  # This is not OpenMPI, so give the old generic warning message
  MESSAGE(WARNING "
             There are known problems with some MPI implementations:
                      MVAPICH2 version <= 1.4.1
             Please consider updating your MPI if applicable.")
endif()
unset(OPENMPI_TYPE)
