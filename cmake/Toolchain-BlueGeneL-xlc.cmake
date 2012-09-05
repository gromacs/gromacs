# derived from http://cmake.org/Wiki/CmakeBlueGene

# the name of the target operating system
set(CMAKE_SYSTEM_NAME BlueGeneL CACHE STRING "Cross-compiling for BlueGene/L")

# adjust to suit your machine's versions
#    /bgl/BlueLight/V1R3M2_140_2007-070424/ppc/bglsys
set(BLRTS_PATH /bgl/BlueLight/V1R3M4_300_2008-080728/ppc/bglsys CACHE STRING "Path to the BlueGene/L system libraries and includes")

# set the compiler
set(CMAKE_C_COMPILER  /opt/ibmcmp/vac/bg/8.0/bin/blrts_xlc)
set(CMAKE_C_FLAGS "-O3 -qbgl -qarch=auto -qtune=auto -qnoautoconfig -qfloat=norngchk -qhot")
set(CMAKE_EXE_LINKER_FLAGS "-L${BLRTS_PATH}/lib")
set(CMAKE_CXX_COMPILER  /opt/ibmcmp/vacpp/bg/8.0/bin/blrts_xlC)

set(MPI_LIBRARY mpich.rts CACHE STRING "MPI library for BlueGene" FORCE)
set(MPI_EXTRA_LIBRARY msglayer.rts devices.rts rts.rts devices.rts CACHE STRING "Extra MPI libraries for BlueGene" FORCE)
set(MPI_INCLUDE_PATH ${BLRTS_PATH}/include  CACHE STRING "MPI include path for BlueGene" FORCE)

# This adds directories that find commands should specifically ignore for cross compiles.
# Most of these directories are the includeand lib directories for the frontend on BG/P systems.
# Not ignoring these can cause things like FindX11 to find a frontend PPC version mistakenly.
# We use this on BG instead of re-rooting because backend libraries are typically strewn about
# the filesystem, and we can't re-root ALL backend libraries to a single place.

set(CMAKE_SYSTEM_IGNORE_PATH
  /lib             /lib64             /include
  /usr/lib         /usr/lib64         /usr/include
  /usr/local/lib   /usr/local/lib64   /usr/local/include
  /usr/X11/lib     /usr/X11/lib64     /usr/X11/include
  /usr/lib/X11     /usr/lib64/X11     /usr/include/X11
  /usr/X11R6/lib   /usr/X11R6/lib64   /usr/X11R6/include
  /usr/X11R7/lib   /usr/X11R7/lib64   /usr/X11R7/include
)

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /bgl/BlueLight/ppcfloor/
    ${BLRTS_PATH}
    /opt/ibmcmp/xlmass/bg
)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search 
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

set(GMX_ACCELERATION "BlueGene" CACHE STRING "Forcing BlueGene acceleration when using BlueGene toolchain")
