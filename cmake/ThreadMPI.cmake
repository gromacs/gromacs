# This source code file is part of thread_mpi.
# Written by Sander Pronk, Erik Lindahl, and possibly others.
#
# Copyright (c) 2009, Sander Pronk, Erik Lindahl.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 1) Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 2) Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 3) Neither the name of the copyright holders nor the
# names of its contributors may be used to endorse or promote products
# derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY US ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL WE BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# If you want to redistribute modifications, please consider that
# scientific software is very special. Version control is crucial -
# bugs must be traceable. We will be happy to consider code for
# inclusion in the official distribution, but derived work should not
# be called official thread_mpi. Details are found in the README & COPYING
# files.

include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)
include(CheckCXXSourceCompiles)

# sets TMPI_ATOMICS to 1 if atomic operations are found, unset otherwise
# Options:
# include path for thread_mpi/atomic.h
function(TMPI_TEST_ATOMICS INCDIR)

    if (NOT DEFINED TMPI_ATOMICS)
        set(CMAKE_REQUIRED_INCLUDES ${INCDIR})
        check_cxx_source_compiles("
// Set a define that forces a compilation error if this platform
// is not yet supported.
#define TMPI_CHECK_ATOMICS

// Include the portable atomics implementation to test
#include \"thread_mpi/atomic.h\"

int main(void)
{
    int i;
    void *ptr;
    tMPI_Atomic_t some_atomic;
    tMPI_Atomic_ptr_t *some_atomic_ptr = NULL;
    tMPI_Spinlock_t some_spinlock;

    /* Make the compiler actually emit code for these functions, so
       that things like inability to emit inline assembly get
       tested. It is not expected that the code below can run. */
    tMPI_Atomic_memory_barrier();
    tMPI_Atomic_memory_barrier_acq();
    tMPI_Atomic_memory_barrier_rel();
    tMPI_Atomic_set(&some_atomic, 0);
    i   = tMPI_Atomic_get(&some_atomic);
    ptr = tMPI_Atomic_ptr_get(some_atomic_ptr);
    tMPI_Atomic_ptr_set(some_atomic_ptr, ptr);
    tMPI_Atomic_add_return(&some_atomic, 0);
    tMPI_Atomic_fetch_add(&some_atomic, 0);
    tMPI_Atomic_cas(&some_atomic, 0, 1);
    tMPI_Atomic_ptr_cas(some_atomic_ptr, ptr, ptr);
    tMPI_Atomic_swap(&some_atomic, 0);
    tMPI_Atomic_ptr_swap(some_atomic_ptr, ptr);
    tMPI_Spinlock_init(&some_spinlock);
    tMPI_Spinlock_lock(&some_spinlock);
    tMPI_Spinlock_trylock(&some_spinlock);
    tMPI_Spinlock_unlock(&some_spinlock);
    tMPI_Spinlock_islocked(&some_spinlock);
    tMPI_Spinlock_wait(&some_spinlock);
    return 0;
}" TEST_ATOMICS)
        if (TEST_ATOMICS)
            message(STATUS "Atomic operations found")
            # If the check fails, we want to be able to check again,
            # in case the user has been able to fix this without
            # needing to delete the cache. Thus we only cache
            # positive results.
            set(TMPI_ATOMICS ${TEST_ATOMICS} CACHE INTERNAL "Whether atomic operations are found")
            set(TMPI_ATOMICS_INCDIR ${INCDIR} CACHE INTERNAL "Atomic operations check include dir")
        else ()
            message(STATUS "Atomic operations not found")
            unset(TEST_ATOMICS CACHE)
        endif()
    endif()

endfunction()

try_compile(HAVE_PROCESSOR_NUMBER ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake/TestWinProcNum.cpp")

include(FindThreads)

if(CMAKE_USE_WIN32_THREADS_INIT AND NOT HAVE_PROCESSOR_NUMBER)
    message(WARNING "Incomplete Windows Processor Group API. If you want GROMACS to be able to set thread affinity, choose a Mingw distribution with a complete API (e.g. Mingw-w64).")
endif()

if (CMAKE_USE_WIN32_THREADS_INIT AND HAVE_PROCESSOR_NUMBER)
    set(THREAD_WINDOWS 1)
    set(THREAD_LIB)
elseif (CMAKE_USE_PTHREADS_INIT)
    check_include_file_cxx(pthread.h    HAVE_PTHREAD_H)
    set(THREAD_PTHREADS 1)
    set(THREAD_LIB ${CMAKE_THREAD_LIBS_INIT})
else()
    message(FATAL_ERROR "Thread support required")
endif ()

# Turns on thread_mpi core threading functions.
MACRO(TMPI_ENABLE_CORE INCDIR)
    TMPI_TEST_ATOMICS(${INCDIR})

# affinity checks
    include(CheckFunctionExists)
    if (THREAD_PTHREADS)
        set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
        # check for sched_setaffinity
        check_cxx_source_compiles(
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
        endif ()
        set(CMAKE_REQUIRED_LIBRARIES)
    endif ()


# this runs on POSIX systems
    check_include_file_cxx(unistd.h             HAVE_UNISTD_H)
    check_include_file_cxx(sched.h              HAVE_SCHED_H)
    check_include_file_cxx(sys/time.h           HAVE_SYS_TIME_H)
    check_cxx_symbol_exists(sysconf    unistd.h HAVE_SYSCONF)
# this runs on windows
#   check_include_files(windows.h               HAVE_WINDOWS_H)
ENDMACRO(TMPI_ENABLE_CORE)

# enable C++ library build.
set(TMPI_CXX_LIB 1)

# Turns on thread_mpi MPI functions.
MACRO(TMPI_ENABLE)
    if(NOT DEFINED TMPI_ATOMICS)
        message(WARNING "Atomic operations not found for this CPU+compiler combination. Thread support will be unbearably slow: disable threads. Atomic operations should work on all but the most obscure CPU+compiler combinations; if your system is not obscure -- like, for example, x86 with gcc --  please contact the developers.")
    endif()

    set(TMPI_ENABLED 1)

# the spin-waiting option
    option(THREAD_MPI_WAIT_FOR_NO_ONE "Use busy waits without yielding to the OS scheduler. Turning this on might improve performance (very) slightly at the cost of very poor performance if the threads are competing for CPU time." OFF)
    mark_as_advanced(THREAD_MPI_WAIT_FOR_NO_ONE)
    if (THREAD_MPI_WAIT_FOR_NO_ONE)
        set(TMPI_WAIT_FOR_NO_ONE 1)
    else ()
        set(TMPI_WAIT_FOR_NO_ONE 0)
    endif ()

# the copy buffer option
    option(THREAD_MPI_COPY_BUFFER "Use an intermediate copy buffer for small message sizes, to allow blocking sends to return quickly. Only useful in programs with relatively uncoupled threads (infrequent MPI communication)" OFF)
    mark_as_advanced(THREAD_MPI_COPY_BUFFER)
    if (THREAD_MPI_COPY_BUFFER)
        set(TMPI_COPY_BUFFER 1)
    else ()
        set(TMPI_COPY_BUFFER 0)
    endif ()

# the profiling option
    option(THREAD_MPI_PROFILING "Turn on simple MPI profiling." OFF)
    mark_as_advanced(THREAD_MPI_PROFILING)
    if (THREAD_MPI_PROFILING)
        set(TMPI_PROFILE 1)
    else ()
        set(TMPI_PROFILE 0)
    endif ()

# tmpi warnings for testing
    option(THREAD_MPI_WARNINGS "Turn thread_mpi warnings for testing." OFF)
    mark_as_advanced(THREAD_MPI_WARNINGS)
    if (THREAD_MPI_WARNINGS)
        set(TMPI_WARNINGS 1)
    else ()
        set(TMPI_WARNINGS 0)
    endif ()
ENDMACRO(TMPI_ENABLE)


MACRO(TMPI_GET_SOURCE_LIST SRC_VARIABLE SRC_ROOT)
    set(${SRC_VARIABLE}
        ${SRC_ROOT}/errhandler.cpp
        ${SRC_ROOT}/tmpi_malloc.cpp
        ${SRC_ROOT}/atomic.cpp
        ${SRC_ROOT}/lock.cpp)

    if (THREAD_PTHREADS)
        list(APPEND ${SRC_VARIABLE} ${SRC_ROOT}/pthreads.cpp)
    elseif (THREAD_WINDOWS)
        list(APPEND ${SRC_VARIABLE} ${SRC_ROOT}/winthreads.cpp)
    endif ()

    if (TMPI_ENABLED)
        list(APPEND ${SRC_VARIABLE}
             ${SRC_ROOT}/alltoall.cpp      ${SRC_ROOT}/p2p_protocol.cpp
             ${SRC_ROOT}/barrier.cpp       ${SRC_ROOT}/p2p_send_recv.cpp
             ${SRC_ROOT}/bcast.cpp         ${SRC_ROOT}/p2p_wait.cpp
             ${SRC_ROOT}/collective.cpp    ${SRC_ROOT}/profile.cpp
             ${SRC_ROOT}/comm.cpp          ${SRC_ROOT}/reduce.cpp
             ${SRC_ROOT}/event.cpp         ${SRC_ROOT}/reduce_fast.cpp
             ${SRC_ROOT}/gather.cpp        ${SRC_ROOT}/scatter.cpp
             ${SRC_ROOT}/group.cpp         ${SRC_ROOT}/tmpi_init.cpp
             ${SRC_ROOT}/topology.cpp      ${SRC_ROOT}/list.cpp
             ${SRC_ROOT}/type.cpp          ${SRC_ROOT}/scan.cpp
             ${SRC_ROOT}/numa_malloc.cpp   ${SRC_ROOT}/once.cpp)
    endif()
ENDMACRO(TMPI_GET_SOURCE_LIST)

