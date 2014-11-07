/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
#include "gmxpre.h"

#include "basenetwork.h"

#include "config.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <exception>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/programcontext.h"

gmx_bool gmx_mpi_initialized(void)
{
#ifndef GMX_MPI
    return 0;
#else
    int n;
    MPI_Initialized(&n);

    return n;
#endif
}

int gmx_node_num(void)
{
#ifndef GMX_MPI
    return 1;
#else
#ifdef GMX_THREAD_MPI
    if (!gmx_mpi_initialized())
    {
        return 1;
    }
#endif
    int i;
    (void) MPI_Comm_size(MPI_COMM_WORLD, &i);
    return i;
#endif
}

int gmx_node_rank(void)
{
#ifndef GMX_MPI
    return 0;
#else
#ifdef GMX_THREAD_MPI
    if (!gmx_mpi_initialized())
    {
        return 0;
    }
#endif
    int i;
    (void) MPI_Comm_rank(MPI_COMM_WORLD, &i);
    return i;
#endif
}

static int mpi_hostname_hash(void)
{
    int hash_int;

#ifndef GMX_LIB_MPI
    /* We have a single physical node */
    hash_int = 0;
#else
    int  resultlen;
    char mpi_hostname[MPI_MAX_PROCESSOR_NAME];

    /* This procedure can only differentiate nodes with different names.
     * Architectures where different physical nodes have identical names,
     * such as IBM Blue Gene, should use an architecture specific solution.
     */
    MPI_Get_processor_name(mpi_hostname, &resultlen);

    /* The string hash function returns an unsigned int. We cast to an int.
     * Negative numbers are converted to positive by setting the sign bit to 0.
     * This makes the hash one bit smaller.
     * A 63-bit hash (with 64-bit int) should be enough for unique node hashes,
     * even on a million node machine. 31 bits might not be enough though!
     */
    hash_int =
        (int)gmx_string_fullhash_func(mpi_hostname, gmx_string_hash_init);
    if (hash_int < 0)
    {
        hash_int -= INT_MIN;
    }
#endif

    return hash_int;
}

#if defined GMX_LIB_MPI && defined GMX_TARGET_BGQ
#ifdef __clang__
/* IBM's declaration of this function in
 * /bgsys/drivers/V1R2M2/ppc64/spi/include/kernel/process.h
 * erroneously fails to specify __INLINE__, despite
 * /bgsys/drivers/V1R2M2/ppc64/spi/include/kernel/cnk/process_impl.h
 * specifiying __INLINE__, so bgclang thinks they are different enough
 * to complain about. */
static uint64_t Kernel_GetJobID();
#endif
#include <spi/include/kernel/location.h>

static int bgq_nodenum(void)
{
    int           hostnum;
    Personality_t personality;
    Kernel_GetPersonality(&personality, sizeof(personality));
    /* Each MPI rank has a unique coordinate in a 6-dimensional space
       (A,B,C,D,E,T), with dimensions A-E corresponding to different
       physical nodes, and T within each node. Each node has sixteen
       physical cores, each of which can have up to four hardware
       threads, so 0 <= T <= 63 (but the maximum value of T depends on
       the confituration of ranks and OpenMP threads per
       node). However, T is irrelevant for computing a suitable return
       value for gmx_hostname_num().
     */
    hostnum  = personality.Network_Config.Acoord;
    hostnum *= personality.Network_Config.Bnodes;
    hostnum += personality.Network_Config.Bcoord;
    hostnum *= personality.Network_Config.Cnodes;
    hostnum += personality.Network_Config.Ccoord;
    hostnum *= personality.Network_Config.Dnodes;
    hostnum += personality.Network_Config.Dcoord;
    hostnum *= personality.Network_Config.Enodes;
    hostnum += personality.Network_Config.Ecoord;

    if (debug)
    {
        std::fprintf(debug,
                     "Torus ID A: %d / %d B: %d / %d C: %d / %d D: %d / %d E: %d / %d\n"
                     "Node ID T: %d / %d core: %d / %d hardware thread: %d / %d\n",
                     personality.Network_Config.Acoord,
                     personality.Network_Config.Anodes,
                     personality.Network_Config.Bcoord,
                     personality.Network_Config.Bnodes,
                     personality.Network_Config.Ccoord,
                     personality.Network_Config.Cnodes,
                     personality.Network_Config.Dcoord,
                     personality.Network_Config.Dnodes,
                     personality.Network_Config.Ecoord,
                     personality.Network_Config.Enodes,
                     Kernel_ProcessorCoreID(),
                     16,
                     Kernel_ProcessorID(),
                     64,
                     Kernel_ProcessorThreadID(),
                     4);
    }
    return hostnum;
}
#endif

int gmx_physicalnode_id_hash(void)
{
    int hash;

#ifndef GMX_MPI
    hash = 0;
#else
#ifdef GMX_THREAD_MPI
    /* thread-MPI currently puts the thread number in the process name,
     * we might want to change this, as this is inconsistent with what
     * most MPI implementations would do when running on a single node.
     */
    hash = 0;
#else
#ifdef GMX_TARGET_BGQ
    hash = bgq_nodenum();
#else
    hash = mpi_hostname_hash();
#endif
#endif
#endif

    if (debug)
    {
        fprintf(debug, "In gmx_physicalnode_id_hash: hash %d\n", hash);
    }

    return hash;
}

#ifdef GMX_LIB_MPI
void gmx_abort(int errorno)
{
    const char *programName = "GROMACS";
    try
    {
        programName = gmx::getProgramContext().displayName();
    }
    catch (const std::exception &)
    {
    }
    const int nnodes   = gmx_node_num();
    const int noderank = gmx_node_rank();
    if (nnodes > 1)
    {
        std::fprintf(stderr, "Halting parallel program %s on rank %d out of %d\n",
                     programName, noderank, nnodes);
    }
    else
    {
        std::fprintf(stderr, "Halting program %s\n", programName);
    }

    MPI_Abort(MPI_COMM_WORLD, errorno);
    std::exit(errorno);
}
#endif
