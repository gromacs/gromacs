
# Note that not all .c files are compiled directly: some of them 
# are #included (some multiple times) from other source files.
set(THREAD_MPI_LIB_SOURCE 
    alltoall.c      reduce.c        scan.c
    barrier.c       list.c          reduce_fast.c
    bcast.c         lock.c          scatter.c
    collective.c    once.c          tmpi_init.c
    comm.c          p2p_protocol.c  topology.c
    errhandler.c    p2p_send_recv.c type.c
    event.c         p2p_wait.c      tmpi_malloc.c
    gather.c        profile.c
    group.c         numa_malloc.c   atomic.c)


if (THREAD_PTHREADS)
    set(THREAD_MPI_LIB_SOURCE pthreads.c ${THREAD_MPI_LIB_SOURCE})
else (THREAD_PTHREADS)
    if (THREAD_WINDOWS)
        set(THREAD_MPI_LIB_SOURCE winthreads.c ${THREAD_MPI_LIB_SOURCE})
    endif (THREAD_WINDOWS)
endif(THREAD_PTHREADS)


add_library(thread_mpi STATIC ${THREAD_MPI_LIB_SOURCE})

# make it link to the threads library (-lpthreads, for example)
target_link_libraries(thread_mpi ${THREAD_LIB})

if (TMPI_CXX_LIB)
    set(THREAD_MPI_CXX_LIB_SOURCE 
        system_error.cpp )
    add_library(thread_mpi_cxx STATIC ${THREAD_MPI_CXX_LIB_SOURCE})
    target_link_libraries(thread_mpi_cxx thread_mpi)
endif (TMPI_CXX_LIB)

#configure_file(tmpi_config.h.cmakein  tmpi_config.h)
#add_definitions(-DHAVE_TMPI_CONFIG_H)
#include_directories( . )

# this configure file only works here. If you include the source
# files into your project, insert the contest of the tmpi_config.h.cmakein
# into your config.h.cmakein file
configure_file(tmpi_config.h.cmakein  ${CMAKE_CURRENT_BINARY_DIR}/tmpi_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-DHAVE_TMPI_CONFIG_H)

install(TARGETS thread_mpi DESTINATION lib)
