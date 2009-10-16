

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

