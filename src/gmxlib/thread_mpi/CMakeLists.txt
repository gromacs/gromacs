

add_library(thread_mpi STATIC threads.c tmpi_init.c comm.c send_recv.c multicast.c)

# make it link to the threads library (-lpthreads, for example)
target_link_libraries(thread_mpi ${THREAD_LIB})
