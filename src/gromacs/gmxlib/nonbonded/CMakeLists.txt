# Sources that should always be built
file(GLOB NONBONDED_SOURCES *.c nb_kernel_c/*.c)

# These sources will be used in the parent directory's CMakeLists.txt
set(NONBONDED_SOURCES ${NONBONDED_SOURCES} PARENT_SCOPE)



