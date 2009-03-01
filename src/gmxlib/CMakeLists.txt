
file(GLOB GMXLIB_SOURCES *.c gmx_blas/*.c gmx_lapack/*.c gmx_statistics/*.c nonbonded/*.c nonbonded/nb_kernel_c/*.c)

add_library(gmx ${GMXLIB_SOURCES})


