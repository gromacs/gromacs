# This file is part of the Collective Variables module (Colvars).
# The original version of Colvars and its updates are located at:
# https://github.com/Colvars/colvars
# Please update all Colvars source files before making any changes.
# If you wish to distribute your changes, please submit them to the
# Colvars repository at GitHub.

function(gmx_manage_lepton)

  # Add Lepton library, which is developed and distributed as part of OpenMM:
  # https://github.com/openmm/openmm

  file(GLOB LEPTON_SOURCES ${PROJECT_SOURCE_DIR}/src/external/lepton/src/*.cpp)
  add_library(lepton OBJECT ${LEPTON_SOURCES})

  target_include_directories(lepton PRIVATE ${PROJECT_SOURCE_DIR}/src/external/lepton/include)
  target_compile_options(lepton PRIVATE -DLEPTON_BUILDING_STATIC_LIBRARY)

  # Set flags so that Colvars can leverage Lepton functionality
  target_include_directories(colvars PRIVATE ${PROJECT_SOURCE_DIR}/src/external/lepton/include)
  target_compile_options(colvars PRIVATE -DLEPTON -DLEPTON_USE_STATIC_LIBRARIES)
endfunction()
