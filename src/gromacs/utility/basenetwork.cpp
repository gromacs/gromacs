/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "basenetwork.h"

#include "config.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"

bool gmx_mpi_initialized()
{
#if !GMX_MPI
    return false;
#else
    int n = 0;
    MPI_Initialized(&n);

    return n != 0;
#endif
}

int gmx_node_num()
{
#if !GMX_MPI
    return 1;
#else
#    if GMX_THREAD_MPI
    if (!gmx_mpi_initialized())
    {
        return 1;
    }
#    endif
    int i = 0;
    (void)MPI_Comm_size(MPI_COMM_WORLD, &i);
    return i;
#endif
}

int gmx_node_rank()
{
#if !GMX_MPI
    return 0;
#else
#    if GMX_THREAD_MPI
    if (!gmx_mpi_initialized())
    {
        return 0;
    }
#    endif
    int i = 0;
    (void)MPI_Comm_rank(MPI_COMM_WORLD, &i);
    return i;
#endif
}

static int mpi_hostname_hash()
{
    int hash_int = 0;

#if GMX_LIB_MPI
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
    hash_int = static_cast<int>(gmx_string_fullhash_func(mpi_hostname, gmx_string_hash_init));
    if (hash_int < 0)
    {
        hash_int -= INT_MIN;
    }
#else

    /* thread-MPI currently puts the thread number in the process name,
     * we might want to change this, as this is inconsistent with what
     * most MPI implementations would do when running on a single node.
     */
    hash_int = 0;
#endif

    return hash_int;
}

int gmx_physicalnode_id_hash()
{
    int hash = 0;

    if (GMX_MPI)
    {
        hash = mpi_hostname_hash();
    }

    if (debug)
    {
        fprintf(debug, "In gmx_physicalnode_id_hash: hash %d\n", hash);
    }

    return hash;
}

#if GMX_LIB_MPI
void gmx_abort(int errorno)
{
    MPI_Abort(MPI_COMM_WORLD, errorno);
    std::abort();
}
#endif
