# the name of the target operating system
set(CMAKE_SYSTEM_NAME BlueGeneP-static CACHE STRING "Cross-compiling for BlueGene/P")

# set the compiler
set(CMAKE_C_COMPILER mpixlc_r)
set(CMAKE_CXX_COMPILER mpixlcxx_r)

set(GMX_ACCELERATION "BlueGene" CACHE STRING "Forcing BlueGene acceleration when using BlueGene toolchain")
