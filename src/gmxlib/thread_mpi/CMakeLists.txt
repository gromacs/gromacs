
# Note that not all .c files are compiled directly: some of them 
# are #included (some multiple times) from other source files.
add_library(thread_mpi STATIC threads.c tmpi_init.c errhandler.c 
            type.c group.c comm.c topology.c p2p.c collective.c 
            once.c list.c reduce_fast.c profile.c event.c barrier.c
            hwinfo.c)

# make it link to the threads library (-lpthreads, for example)
target_link_libraries(thread_mpi ${THREAD_LIB})


#configure_file(tmpi_config.h.cmakein  tmpi_config.h)
#add_definitions(-DHAVE_TMPI_CONFIG_H)
#include_directories( . )

# this configure file only works here. If you include the source
# files into your project, insert the contest of the tmpi_config.h.cmakein
# into your config.h.cmakein file
configure_file(tmpi_config.h.cmakein  ${CMAKE_CURRENT_BINARY_DIR}/tmpi_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-DHAVE_TMPI_CONFIG_H)

