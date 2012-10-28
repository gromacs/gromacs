# Sources that should always be built
file(GLOB NONBONDED_SOURCES *.c nb_kernel_c/*.c)
file(GLOB TEMPLATE_SOURCES *.pre*.c nb_kernel_c/*.pre*.c)

if(TEMPLATE_SOURCES)
        list(REMOVE_ITEM NONBONDED_SOURCES ${TEMPLATE_SOURCES})
endif(TEMPLATE_SOURCES)

list(APPEND NONBONDED_SOURCES)

# These sources will be used in the parent directory's CMakeLists.txt
set(NONBONDED_SOURCES ${NONBONDED_SOURCES} PARENT_SCOPE)



