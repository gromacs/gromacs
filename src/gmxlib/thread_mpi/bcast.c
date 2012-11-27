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
#include <stdarg.h>
#include <string.h>

#include "impl.h"
#include "collective.h"


/* broadcast */
int tMPI_Bcast(void* buffer, int count, tMPI_Datatype datatype, int root,
               tMPI_Comm comm)
{
    int synct;
    struct coll_env *cev;
    int myrank;
    int ret=TMPI_SUCCESS;
    struct tmpi_thread *cur=tMPI_Get_current();

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur); 
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Bcast(%p, %d, %p, %d, %p)", buffer, count, datatype, 
                     root, comm);
#endif

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, cur);

    /* we increase our counter, and determine which coll_env we get */
    cev=tMPI_Get_cev(comm, myrank, &synct);

    if (myrank==root)
    {
        /* first set up the data */
        tMPI_Post_multi(cev, myrank, 0, TMPI_BCAST_TAG, datatype, 
                        count*datatype->size, buffer, comm->grp.N-1, synct, -1);
        /* and wait until everybody is done copying */
        tMPI_Wait_for_others(cev, myrank);
    }
    else
    {
        size_t bufsize=count*datatype->size;
        /* wait until root becomes available */
        tMPI_Wait_for_data(cur, cev, myrank);
        tMPI_Mult_recv(comm, cev, root, 0, TMPI_BCAST_TAG, datatype, bufsize, 
                       buffer, &ret);
    }
#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Bcast); 
#endif
    return ret;
}





