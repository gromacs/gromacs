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

To help us fund development, we humbly ask that you cite
any papers on the package - you can find them in the top README file.
*/

/*! \mainpage thread_mpi
 *
 * \section intro_sec 
 *
 * Threading library for cross-platform high-performance computing. 
 * Contains an implementation of MPI for use with threads, threading
 * basics, as well as atomic operations, collective shared-memory allocators 
 * and lock-free data structures. 
 * 
 * Although MPI-style message passing isn't neccesarily optimal for 
 * performance on shared-memory systems, the MPI communicator concept and 
 * its emphasis on collective operations make sense even when computing on
 * one machine with multiple cores. 
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *  
 * etc...
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

#include "thread_mpi/threads.h"
#include "thread_mpi/atomic.h"
#include "thread_mpi/tmpi.h"
#include "thread_mpi/mpi.h"


