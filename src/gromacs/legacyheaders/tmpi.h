/*
This source code file is part of thread_mpi.  
Written by Sander Pronk, Erik Lindahl, and possibly others. 

Copyright (c) 2009, Sander Pronk, Erik Lindahl.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1) Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2) Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3) Neither the name of the copyright holders nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY US ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL WE BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

If you want to redistribute modifications, please consider that
scientific software is very special. Version control is crucial -
bugs must be traceable. We will be happy to consider code for
inclusion in the official distribution, but derived work should not
be called official thread_mpi. Details are found in the README & COPYING
files.
*/


/*! 
  \if THREAD_MPI_MAIN
  \mainpage thread_mpi thread_mpi threading library
  \else
  \page thread_mpi thread_mpi threading library
  \endif

   thread_mpi is a cross-platform threading library for applications in 
   high-performance computing. It supports:

   - Cross-platform thread primitives (thread creation, mutexes, spinlocks, 
     barriers, thread-local storage, etc.).
   - Cross-platform atomic operations (compare-and-swap, add-return, etc) for
     safe lock-free synchronization.
   - An implementation of (currently, much of) MPI, either as a drop-in
     replacement, or for use in conjunction with a networked MPI      
     implementation. 
   - Shared-memory allocation and memory management (planned, as of now). 
   - Basic lock-free data structures (planned, as of now).

   Because it can be used as a drop-in replacement for MPI, existing codes
   using MPI can start using thread_mpi without major changes in the 
   source code, assuming -- and this is a big assumption -- that the code 
   is thread-safe.

   Alternatively, networked MPI calls can be used in conjunction with    
   thread_mpi calls (simply by using     
    "#include <thread_mpi.h>"
   instead of 
    "#include <tmpi.h>"
   and pre-fixing all thread_mpi MPI-like calls with tMPI instead of MPI. 

   The availability of both MPI calls and shared-memory constructs makes it 
   possible to transition (relatively) seamlessly from an MPI-style code 
   to code that's optimal on multicore CPUs.

   Although MPI-style message passing isn't neccesarily optimal for 
   performance on shared-memory systems, the MPI communicator concept and 
   its emphasis on collective operations makes sense even when computing on
   one machine with multiple cores. The communicator forms the basis for
   the shared-memory allocation and lock-free data structure implementations
   in thread_mpi.

   Although usable as a stand-alone library, thread_mpi is designed to
   be incorporated in the code tree, eliminating any external build
   requirements. The BSD-style license that this library is distributed 
   with reflects this.

   Thread primitives and the atomic operations are cpu and operating system 
   dependent - thread_mpi attempts to make them available with the same 
   interface independently of the platform it's run on. 
   Currently the thread primitives are supported on:
   - any operating system supporting POSIX threads
   - Windows (XP and later). 

   The atomic operations (such as compare-and-swap) are supported on:
   - gcc on x86, x86_64, PowerPC and Itanium.
   - Intel compilers on x86, x86_64 and Itanium. 
   - xlc on PowerPC.
   - (partial) HP/UX compilers on Itanium.

   Detailed descriptions of the parts of the API can be found in:
   - thread_mpi/threads.h for threading fundamentals.
   - thread_mpi/atomic.h for atomic operations.
   - thread_mpi/tmpi.h for the MPI functions as tMPI_-prefixed functions.
   - thread_mpi/mpi_bindings.h for the MPI bindings.
*/


/** \file 
 *
 * \brief Convenience header file for MPI compatibility.
 * 
 * This file includes the tMPI header file thread_mpi/tmpi.h and the true
 * MPI-style bindings of thread_mpi/mpi.h, as well as thread_mpi/threads.h and
 * thread_mpi/atomic.h header files. If you'd like to use the components
 * individually, or be able to use a networked MPI together with thread_mpi,
 * include the relevant header files directly. 
 */


#include "thread_mpi/atomic.h"
#include "thread_mpi/threads.h"
#include "thread_mpi/barrier.h"
#include "thread_mpi/event.h"
#include "thread_mpi/lock.h"
#include "thread_mpi/tmpi.h"
#include "thread_mpi/collective.h"
#include "thread_mpi/hwinfo.h"

#include "thread_mpi/mpi_bindings.h"


