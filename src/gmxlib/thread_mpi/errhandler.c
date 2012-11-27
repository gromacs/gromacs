/*
 * This file is part of the GROMACS molecular simulation package,
 * version 4.6
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
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
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if ! (defined( _WIN32 ) || defined( _WIN64 ) )
#include <sys/time.h>
#endif

#include "impl.h"






struct tmpi_errhandler_ tmpi_errors_are_fatal = { 0, tmpi_errors_are_fatal_fn };
struct tmpi_errhandler_ tmpi_errors_return = { 0, tmpi_errors_return_fn };


tMPI_Errhandler TMPI_ERRORS_ARE_FATAL=&tmpi_errors_are_fatal;
tMPI_Errhandler TMPI_ERRORS_RETURN=&tmpi_errors_return;




/* error messages. Must match error codes in thread_mpi.h */
static const char *tmpi_errmsg[] =
{
    "No error",
    "malloc failure in tMPI (out of memory)",
    "tMPI Initialization error",
    "tMPI Finalize error",
    "Invalid tMPI_Group",
    "Invalid tMPI_Comm",
    "Invalid tMPI_Status",
    "Invalid tMPI_Group rank",
    "Invalid Cartesian topology dimensions",
    "Invalid Cartesian topology coordinates",
    "Insufficient number processes for Cartesian topology",
    "Invalid counterpart for MPI transfer",
    "Receive buffer size too small for transmission",
    "Overlapping send/receive buffers: probably due to thread-unsafe code.",
    "Invalid send destination",
    "Invalid receive source",
    "Invalid buffer (null pointer in send or receive buffer)",
    "Multicast operation mismatch (multicast not collective across comm)",
    "Invalid reduce operator",
    "Out of receive envelopes: this shouldn't happen (probably a bug).",
    "Out of receive requests: this shouldn't happen (probably a bug).",
    "Transmission failure",
    "Unknown tMPI error"
};



int tMPI_Error(tMPI_Comm comm, int tmpi_errno)
{
    if (comm)
    {
        comm->erh->err=tmpi_errno;
        comm->erh->fn(&comm, &tmpi_errno);
    }
    else
    {
        /* initialization errors have no comm */
        tmpi_errors_are_fatal_fn(NULL, &tmpi_errno);
    }
    return tmpi_errno;
}


int tMPI_Error_string(int errorcode, char *strn, size_t *resultlen)
{
    if (errorcode<0 || errorcode>=N_TMPI_ERR)
        errorcode=TMPI_ERR_UNKNOWN;

#if ! (defined( _WIN32 ) || defined( _WIN64 ) )
    strncpy(strn, tmpi_errmsg[errorcode], TMPI_MAX_ERROR_STRING);
#else
    strncpy_s(strn, TMPI_MAX_ERROR_STRING, tmpi_errmsg[errorcode], TMPI_MAX_ERROR_STRING);
#endif
    *resultlen=strlen(strn);
    return TMPI_SUCCESS;
}

int tMPI_Create_errhandler(tMPI_Errhandler_fn *function, 
                           tMPI_Errhandler *errhandler) 
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Create_errhandler(%p, %p)", function, errhandler);
#endif

    /* we don't use a special malloc here because this is the error handler
       creation function. */
    *errhandler=(tMPI_Errhandler)malloc(sizeof(struct tmpi_errhandler_));
    if (!*errhandler)
    {
        fprintf(stderr, "tMPI fatal error (%s), bailing out\n", 
                tmpi_errmsg[TMPI_ERR_MALLOC]);
        abort();
    }
    (*errhandler)->err=0;
    (*errhandler)->fn=*function;
    return TMPI_SUCCESS;
}

int tMPI_Errhandler_free(tMPI_Errhandler *errhandler)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Errhandler_free(%p)", errhandler);
#endif

    free(*errhandler);
    return TMPI_SUCCESS;
}


int tMPI_Comm_set_errhandler(tMPI_Comm comm, tMPI_Errhandler errhandler)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_set_errhandler(%p, %p)", comm, errhandler);
#endif

    comm->erh = errhandler;
    return TMPI_SUCCESS;
}

int tMPI_Comm_get_errhandler(tMPI_Comm comm, tMPI_Errhandler *errhandler)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_get_errhandler(%p, %p)", comm, errhandler);
#endif

    *errhandler=comm->erh;
    return TMPI_SUCCESS;
}

void tmpi_errors_are_fatal_fn(tMPI_Comm *comm, int *err)
{
    char errstr[TMPI_MAX_ERROR_STRING];
    size_t len;

    tMPI_Error_string(*err, errstr, &len);
    if (comm)
    {
        fprintf(stderr, "tMPI error: %s (in valid comm)\n", errstr);
    }
    else
    {
        fprintf(stderr, "tMPI error: %s\n", errstr);
    }
    abort();
    /*exit(0);*/
}

void tmpi_errors_return_fn(tMPI_Comm *comm, int *err)
{
    char errstr[TMPI_MAX_ERROR_STRING];
    size_t len;

    tMPI_Error_string(*err, errstr, &len);
    if (comm)
    {
        fprintf(stderr, "tMPI error: %s (in valid comm)\n", errstr);
    }
    else
    {
        fprintf(stderr, "tMPI error: %s\n", errstr);
    }
    return;
}

