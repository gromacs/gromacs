# Sources that should always be built
file(GLOB NONBONDED_SOURCES *.c nb_kernel_c/*.c nb_kernel_adress_c/*.c)

# Conditionally included ones
if(GMX_FORTRAN)
  if (GMX_DOUBLE)
    file(GLOB NB_KERNEL_FORTRAN_SOURCES nb_kernel_f77_double/*.[cf])
  else(GMX_DOUBLE)
    file(GLOB NB_KERNEL_FORTRAN_SOURCES nb_kernel_f77_single/*.[cf])
  endif(GMX_DOUBLE)
endif(GMX_FORTRAN)

if(GMX_POWER6)
  file(GLOB NB_KERNEL_FORTRAN_SOURCES nb_kernel_power6/*.[cF])
endif(GMX_POWER6)

if(GMX_BLUEGENE)
  file(GLOB NB_KERNEL_BLUEGENE_SOURCES nb_kernel_bluegene/*.c)
endif(GMX_BLUEGENE)

list(APPEND NONBONDED_SOURCES ${NB_KERNEL_FORTRAN_SOURCES} ${NB_KERNEL_BLUEGENE_SOURCES})

# These sources will be used in the parent directory's CMakeLists.txt
set(NONBONDED_SOURCES ${NONBONDED_SOURCES} PARENT_SCOPE)



