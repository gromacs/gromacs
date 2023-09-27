cmake_minimum_required(VERSION 3.1)

# C language is only needed due to regression in FindHDF5.
project(EXAMPLE LANGUAGES C CXX)

find_package(H5XX REQUIRED)

find_package(Boost 1.40.0 QUIET REQUIRED)
include_directories(SYSTEM "${Boost_INCLUDE_DIR}")

set(CMAKE_BUILD_TYPE "Release")

# --- non-MPI programs
foreach(module
  attribute
  dataset_multiarray
  slice_multiarray
  slice_array
  slice_std_vector
  slice_python
  )
  add_executable(${module}
    ${module}.cpp
  )
  target_link_libraries(${module}
    h5xx
    dl
    z
  )
endforeach()

if (MPI_FOUND)
  foreach(module
    slice_mpi
    )
    add_executable(${module}
      ${module}.cpp
    )
    target_link_libraries(${module}
      h5xx
      dl
      z
    )
  endforeach()
endif()
