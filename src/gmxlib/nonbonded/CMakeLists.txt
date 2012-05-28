include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB NONBONDED_SOURCES *.c)
file(GLOB NB_KERNEL_C_SOURCES nb_kernel_c/*.c)
file(GLOB NB_KERNEL_ADRESS_C_SOURCES nb_kernel_adress_c/*.c)

if(GMX_FORTRAN)
  if (GMX_DOUBLE)
    file(GLOB NB_KERNEL_FORTRAN_SOURCES nb_kernel_f77_double/*.[cf])
  else(GMX_DOUBLE)
    file(GLOB NB_KERNEL_FORTRAN_SOURCES nb_kernel_f77_single/*.[cf])
  endif(GMX_DOUBLE)
endif(GMX_FORTRAN)

if(GMX_POWER6)
  file(GLOB NB_KERNEL_FORTRAN_SOURCES nonbonded/nb_kernel_power6/*.[cF])
endif(GMX_POWER6)

if(GMX_BLUEGENE)
  file(GLOB NB_KERNEL_BLUEGENE_SOURCES nonbonded/nb_kernel_bluegene/*.c)
endif(GMX_BLUEGENE)

add_library(nonbonded OBJECT ${NONBONDED_SOURCES} ${NB_KERNEL_C_SOURCES} ${NB_KERNEL_ADRESS_C_SOURCES} ${NB_KERNEL_FORTRAN_SOURCES} ${NB_KERNEL_BLUEGENE_SOURCES})

