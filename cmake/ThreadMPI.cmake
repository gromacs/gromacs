
include(CheckIncludeFiles)
include(CheckFunctionExists)
#include(CheckCSourceCompiles)

#option(THREAD_PTHREADS "Use posix threads" ON)

MACRO(TEST_TMPI_ATOMICS VARIABLE)
    if (NOT DEFINED TMPI_ATOMICS)
        try_compile(TEST_ATOMICS "${CMAKE_BINARY_DIR}"
                "${CMAKE_SOURCE_DIR}/cmake/TestAtomics.c"
                COMPILE_DEFINITIONS "-I${CMAKE_SOURCE_DIR}/include" )

        if (TEST_ATOMICS)
            message(STATUS "Atomics found")
            set(${VARIABLE} CACHE INTERNAL 1)
        else (TEST_ATOMICS)
            message(WARNING "Atomics not found for this compiler+cpu combination. Thread support will be unbearably slow: disable threads. Atomics should work on all but the most obscure CPU+compiler combinations; if your system is not obscure -- like, for example, x86 with gcc --  please contact the developers.")
            set(${VARIABLE} CACHE INTERNAL 0)
        endif(TEST_ATOMICS)
    endif(NOT DEFINED TMPI_ATOMICS)
ENDMACRO(TEST_TMPI_ATOMICS VARIABLE)


include(FindThreads)
if (CMAKE_USE_PTHREADS_INIT)
    check_include_files(pthread.h    HAVE_PTHREAD_H)
    set(THREAD_PTHREADS 1)
    #add_definitions(-DTHREAD_PTHREADS)
    set(THREAD_MPI_SRC 
        thread_mpi/alltoall.c      thread_mpi/p2p_protocol.c
        thread_mpi/barrier.c       thread_mpi/p2p_send_recv.c
        thread_mpi/bcast.c         thread_mpi/p2p_wait.c
        thread_mpi/collective.c    thread_mpi/profile.c
        thread_mpi/comm.c          thread_mpi/pthreads.c
        thread_mpi/errhandler.c    thread_mpi/reduce.c
        thread_mpi/event.c         thread_mpi/reduce_fast.c
        thread_mpi/gather.c        thread_mpi/scatter.c
        thread_mpi/group.c         thread_mpi/tmpi_init.c
        thread_mpi/topology.c      thread_mpi/list.c          
        thread_mpi/type.c          thread_mpi/lock.c
        thread_mpi/numa_malloc.c   thread_mpi/once.c)
    set(THREAD_LIB ${CMAKE_THREAD_LIBS_INIT})
else (CMAKE_USE_PTHREADS_INIT)
    if (CMAKE_USE_WIN32_THREADS_INIT)
        set(THREAD_WINDOWS 1)
        #add_definitions(-DTHREAD_WINDOWS)
        set(THREAD_MPI_SRC 
            thread_mpi/alltoall.c      thread_mpi/p2p_protocol.c
            thread_mpi/barrier.c       thread_mpi/p2p_send_recv.c
            thread_mpi/bcast.c         thread_mpi/p2p_wait.c
            thread_mpi/collective.c    thread_mpi/profile.c
            thread_mpi/comm.c          
            thread_mpi/errhandler.c    thread_mpi/reduce.c
            thread_mpi/event.c         thread_mpi/reduce_fast.c
            thread_mpi/gather.c        thread_mpi/scatter.c
            thread_mpi/group.c         thread_mpi/tmpi_init.c
            thread_mpi/topology.c      thread_mpi/list.c
            thread_mpi/type.c          thread_mpi/lock.c
            thread_mpi/winthreads.c    thread_mpi/once.c
            thread_mpi/numa_malloc.c)
        set(THREAD_LIBRARY )
    endif (CMAKE_USE_WIN32_THREADS_INIT)
endif (CMAKE_USE_PTHREADS_INIT)

# the spin-waiting option
option(THREAD_MPI_WAIT_FOR_NO_ONE "Use busy waits without yielding to the OS scheduler. Turning this on might improve performance (very) slightly at the cost of very poor performance if the threads are competing for CPU time." OFF)
mark_as_advanced(THREAD_MPI_WAIT_FOR_NO_ONE)
if (THREAD_MPI_WAIT_FOR_NO_ONE)
    add_definitions(-DTMPI_WAIT_FOR_NO_ONE)
else (THREAD_MPI_WAIT_FOR_NO_ONE)
    add_definitions()
endif (THREAD_MPI_WAIT_FOR_NO_ONE)


# the copy buffer option
option(THREAD_MPI_COPY_BUFFER "Use an intermediate copy buffer for small message sizes, to allow blocking sends to return quickly." ON)
mark_as_advanced(THREAD_MPI_COPY_BUFFER)
if (THREAD_MPI_COPY_BUFFER)
    add_definitions()
else (THREAD_MPI_COPY_BUFFER)
    add_definitions(-DTMPI_NO_COPY_BUFFER)
endif (THREAD_MPI_COPY_BUFFER)


# the profiling option
option(THREAD_MPI_PROFILING "Turn on simple MPI profiling." OFF)
mark_as_advanced(THREAD_MPI_PROFILING)
if (THREAD_MPI_PROFILING)
    add_definitions(-DTMPI_PROFILE)
else (THREAD_MPI_PROFILING)
    add_definitions()
endif (THREAD_MPI_PROFILING)

include(CheckCSourceCompiles)

# Windows NUMA allocator
if (THREAD_WINDOWS)
	check_c_source_compiles(
	"#include <windows.h>
	int main(void) { PROCESSOR_NUMBER a; return 0; }"
	HAVE_PROCESSOR_NUMBER)
	if(HAVE_PROCESSOR_NUMBER)
            #add_definitions(-DTMPI_WINDOWS_NUMA_API)
            set(TMPI_WINDOWS_NUMA_API 1)
	endif(HAVE_PROCESSOR_NUMBER)
endif(THREAD_WINDOWS)



include(CheckFunctionExists)
if (THREAD_PTHREADS)
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
    # check for sched_setaffinity
    check_c_source_compiles(
        "#define _GNU_SOURCE
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
int main(void) { cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(0, &set);
    pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
    return 0;
}"
        PTHREAD_SETAFFINITY
    )
    if (PTHREAD_SETAFFINITY)
        set(HAVE_PTHREAD_SETAFFINITY 1)
    endif (PTHREAD_SETAFFINITY)
endif (THREAD_PTHREADS)


# this runs on POSIX systems
check_include_files(unistd.h        HAVE_UNISTD_H)
check_include_files(sched.h         HAVE_SCHED_H)
check_include_files(sys/time.h      HAVE_SYS_TIME_H)
check_function_exists(sysconf       HAVE_SYSCONF)
# this runs on windows
#check_include_files(windows.h		HAVE_WINDOWS_H)
#check_function_exists(GetSystemInfo HAVE_SYSTEM_INFO)

test_tmpi_atomics(TMPI_ATOMICS)
