
include(CheckIncludeFiles)
include(CheckFunctionExists)
#include(CheckCSourceCompiles)

#option(THREAD_PTHREADS "Use posix threads" ON)

MACRO(TEST_TMPI_ATOMICS VARIABLE)
    if (NOT DEFINED TMPI_ATOMICS)
        try_compile(TEST_ATOMICS "${CMAKE_BINARY_DIR}"
                "${CMAKE_SOURCE_DIR}/cmake/TestAtomics.c"
                COMPILE_DEFINITIONS "-I${CMAKE_SOURCE_DIR}/src/gromacs/legacyheaders" )

        if (TEST_ATOMICS)
            message(STATUS "Atomics found")
            set(${VARIABLE} CACHE INTERNAL 1)
        else (TEST_ATOMICS)
            message(WARNING "Atomics not found for this compiler+cpu combination. Thread support will be unbearably slow: disable threads. Atomics should work on all but the most obscure CPU+compiler combinations; if your system is not obscure -- like, for example, x86 with gcc --  please contact the developers.")
            set(${VARIABLE} CACHE INTERNAL 0)
        endif(TEST_ATOMICS)
    endif(NOT DEFINED TMPI_ATOMICS)
ENDMACRO(TEST_TMPI_ATOMICS VARIABLE)

MACRO(TMPI_MAKE_CXX_LIB)
    set(TMPI_CXX_LIB 1)
ENDMACRO(TMPI_MAKE_CXX_LIB)

MACRO(TMPI_GET_SOURCE_LIST SRC_VARIABLE)
    foreach (_option IN ITEMS ${ARGN})
        if (_option STREQUAL "CXX")
            set(TMPI_CXX_LIB 1)
        elseif (_option STREQUAL "NOMPI")
            set(TMPI_NO_MPI_LIB 1)
        else ()
            message(FATAL_ERROR "Unknown thread_mpi option '${_option}'")
        endif ()
    endforeach ()
    set(${SRC_VARIABLE}
        thread_mpi/errhandler.c
        thread_mpi/tmpi_malloc.c)
    if (THREAD_PTHREADS)
        list(APPEND ${SRC_VARIABLE} thread_mpi/pthreads.c)
    elseif (THREAD_WINDOWS)
        list(APPEND ${SRC_VARIABLE} thread_mpi/winthreads.c)
    endif (THREAD_PTHREADS)
    if (TMPI_CXX_LIB)
        list(APPEND ${SRC_VARIABLE} thread_mpi/system_error.cpp)
    endif (TMPI_CXX_LIB)
    if (NOT TMPI_NO_MPI_LIB)
        list(APPEND ${SRC_VARIABLE}
             thread_mpi/alltoall.c      thread_mpi/p2p_protocol.c
             thread_mpi/barrier.c       thread_mpi/p2p_send_recv.c
             thread_mpi/bcast.c         thread_mpi/p2p_wait.c
             thread_mpi/collective.c    thread_mpi/profile.c
             thread_mpi/comm.c          thread_mpi/reduce.c
             thread_mpi/event.c         thread_mpi/reduce_fast.c
             thread_mpi/gather.c        thread_mpi/scatter.c
             thread_mpi/group.c         thread_mpi/tmpi_init.c
             thread_mpi/topology.c      thread_mpi/list.c
             thread_mpi/type.c          thread_mpi/lock.c
             thread_mpi/numa_malloc.c   thread_mpi/once.c
             thread_mpi/scan.c)
    endif()
ENDMACRO(TMPI_GET_SOURCE_LIST)

include(FindThreads)
if (CMAKE_USE_PTHREADS_INIT)
    check_include_files(pthread.h    HAVE_PTHREAD_H)
    set(THREAD_PTHREADS 1)
    #add_definitions(-DTHREAD_PTHREADS)
    set(THREAD_LIB ${CMAKE_THREAD_LIBS_INIT})
else (CMAKE_USE_PTHREADS_INIT)
    if (CMAKE_USE_WIN32_THREADS_INIT)
        set(THREAD_WINDOWS 1)
        #add_definitions(-DTHREAD_WINDOWS)
        set(THREAD_LIB)
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

# option to set affinity 
option(THREAD_MPI_SET_AFFINITY "Set thread affinity to a core if number of threads equal to number of hardware threads." ON)
mark_as_advanced(THREAD_MPI_SET_AFFINITY)
if (THREAD_MPI_SET_AFFINITY)
    add_definitions(-DTMPI_SET_AFFINITY)
else (THREAD_MPI_SET_AFFINITY)
    add_definitions()
endif (THREAD_MPI_SET_AFFINITY)

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


test_tmpi_atomics(TMPI_ATOMICS)
