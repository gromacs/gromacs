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


#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#else
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef THREAD_WINDOWS
#include <windows.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "impl.h"
#include "thread_mpi/hwinfo.h"

#ifdef TMPI_TRACE
#include <stdarg.h>
#endif


int tMPI_Get_hw_nthreads(void)
{
    int ret=-1;
#ifdef HAVE_SYSCONF
#if defined(_SC_NPROCESSORS_ONLN)
    ret=sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
    ret=sysconf(_SC_NPROC_ONLN);
#elif defined(_SC_NPROCESSORS_CONF)
    ret=sysconf(_SC_NPROCESSORS_CONF);
#elif defined(_SC_NPROC_CONF)
    ret=sysconf(_SC_NPROC_CONF);
#endif
#endif

#ifdef THREAD_WINDOWS
    SYSTEM_INFO sysinfo;
    GetSystemInfo( &sysinfo );

    ret=sysinfo.dwNumberOfProcessors;
#endif
    return ret;
}


int tMPI_Get_recommended_nthreads(void)
{
    int N=1; /* the default is 1 */

#ifndef TMPI_NO_ATOMICS
    N=tMPI_Get_hw_nthreads();
    if (N<1)
        N=1;
#endif
    return N;
}

