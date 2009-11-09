

include(CheckIncludeFiles)

#option(THREAD_PTHREADS "Use posix threads" ON)

include(FindThreads)
if (CMAKE_USE_PTHREADS_INIT)
    check_include_files(pthread.h    HAVE_PTHREAD_H)
    #set(THREAD_PTHREADS 1)
    add_definitions(-DTHREAD_PTHREADS)
    set(THREAD_MPI_SRC thread_mpi/threads.c thread_mpi/tmpi_init.c 
                       thread_mpi/errhandler.c thread_mpi/type.c
                       thread_mpi/group.c thread_mpi/comm.c 
                       thread_mpi/topology.c thread_mpi/send_recv.c 
                       thread_mpi/collective.c)
    set(THREAD_LIB ${CMAKE_THREAD_LIBS_INIT})
else (CMAKE_USE_PTHREADS_INIT)
    if (CMAKE_USE_WIN32_THREADS_INIT)
        set(THREAD_WINDOWS 1)
        add_definitions(-DTHREAD_WINDOWS)
        set(THREAD_MPI_SRC thread_mpi/threads.c thread_mpi/tmpi_init.c 
                           thread_mpi/errhandler.c thread_mpi/type.c
                           thread_mpi/group.c thread_mpi/comm.c 
                           thread_mpi/topology.c thread_mpi/send_recv.c 
                           thread_mpi/collective.c)
        set(THREAD_LIBRARY )
    endif (CMAKE_USE_WIN32_THREADS_INIT)
endif (CMAKE_USE_PTHREADS_INIT)

# the busy waiting option
option(THREAD_MPI_BUSY_WAIT "Use busy waits for thread_mpi synchronization. Provides lower latency, but higher unneccesary CPU usage." ON)
mark_as_advanced(THREAD_MPI_BUSY_WAIT)
if (THREAD_MPI_BUSY_WAIT)
    add_definitions()
else (THREAD_MPI_BUSY_WAIT)
    add_definitions(-DTMPI_NO_BUSY_WAIT)
endif (THREAD_MPI_BUSY_WAIT)

# the copy buffer option
option(THREAD_MPI_COPY_BUFFER "Use an intermediate copy buffer for small message sizes, to allow blocking sends to return quickly." ON)
mark_as_advanced(THREAD_MPI_COPY_BUFFER)
if (THREAD_MPI_COPY_BUFFER)
    add_definitions()
else (THREAD_MPI_COPY_BUFFER)
    add_definitions(-DTMPI_NO_COPY_BUFFER)
endif (THREAD_MPI_COPY_BUFFER)

