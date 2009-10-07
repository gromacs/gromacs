
# Note that not all .c files are compiled directly: some of them 
# are #included (some multiple times) from other source files.
add_library(thread_mpi STATIC threads.c tmpi_init.c errhandler.c type.c group.c comm.c topology.c send_recv.c collective.c list.c)

# make it link to the threads library (-lpthreads, for example)
target_link_libraries(thread_mpi ${THREAD_LIB})
