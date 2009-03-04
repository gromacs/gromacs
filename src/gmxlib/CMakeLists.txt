
# The nonbonded directory contains subdirectories that are only
# conditionally built, so we cannot use a GLOB_RECURSE here.

file(GLOB GMXLIB_SOURCES *.c gmx_blas/*.c gmx_lapack/*.c 
     gmx_statistics/*.c nonbonded/*.c nonbonded/nb_kernel_c/*.c)

# Files called xxx_test.c are test drivers with a main() function for module xxx.c,
# so they should not be included in the library
file(GLOB_RECURSE TEST_DRIVERS *_test.c)
list(REMOVE_ITEM GMXLIB_SOURCES ${TEST_DRIVERS})  

add_library(gmx ${GMXLIB_SOURCES})

install(TARGETS gmx DESTINATION ${LIB_INSTALL_DIR})

