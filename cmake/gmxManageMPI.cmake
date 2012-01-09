# Manage the MPI setup, assuming that CMAKE_C_COMPILER is an MPI
# (wrapper) compiler.
if(GMX_MPI)
  if(GMX_THREADS)
    set(GMX_THREADS OFF CACHE BOOL
      "Thread-based parallelization conflicts with MPI." FORCE)
  endif(GMX_THREADS)

  # Test the CMAKE_C_COMPILER for being an MPI (wrapper) compiler
  TRY_COMPILE(MPI_FOUND ${CMAKE_BINARY_DIR}
    "${CMAKE_SOURCE_DIR}/cmake/TestMPI.c"
    COMPILE_DEFINITIONS )

  if(MPI_FOUND)
    if(GMX_FAHCORE)
      add_definitions( -DMPI ) #for FAHCORE
    endif()
    include(gmxTestMPI_IN_PLACE)
    if (GMX_MPI_IN_PLACE)
      gmx_test_mpi_in_place(MPI_IN_PLACE_EXISTS)
    endif()

    # Test for unsuitable versions of MPI
    include(gmxCheckMPI)

  else(MPI_FOUND)
    message(FATAL_ERROR "MPI support requested, but no MPI compiler found.")
  endif(MPI_FOUND)

  include(gmxTestCatamount)
  gmx_test_catamount(GMX_CRAY_XT3)
  if(GMX_CRAY_XT3)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_CRAY_XT3")
    set(GMX_NO_SYSTEM 1)
    set(GMX_NO_NICE 1)
  endif(GMX_CRAY_XT3)

  set(GMX_LIB_MPI 1)
  set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_LIB_MPI")
endif(GMX_MPI)
