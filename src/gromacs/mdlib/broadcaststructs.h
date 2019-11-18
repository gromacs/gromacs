/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018,2019, by the GROMACS development team, led by
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

/*! \libinternal \file
 *
 * \brief Convenience wrappers for broadcasting structs.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_BROADCASTSTRUCTS_H
#define GMX_MDLIB_BROADCASTSTRUCTS_H

#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"

struct gmx_mtop_t;
struct t_commrec;
struct t_inputrec;
struct PartialDeserializedTprFile;
class t_state;

//! Convenience wrapper for gmx_bcast of a single value.
template<typename T>
void block_bc(const t_commrec* cr, T& data)
{
    gmx_bcast(sizeof(T), static_cast<void*>(&data), cr);
}
//! Convenience wrapper for gmx_bcast of a C-style array.
template<typename T>
void nblock_bc(const t_commrec* cr, int numElements, T* data)
{
    gmx_bcast(numElements * sizeof(T), static_cast<void*>(data), cr);
}
//! Convenience wrapper for gmx_bcast of an ArrayRef<T>
template<typename T>
void nblock_bc(const t_commrec* cr, gmx::ArrayRef<T> data)
{
    gmx_bcast(data.size() * sizeof(T), static_cast<void*>(data.data()), cr);
}
//! Convenience wrapper for allocation with snew of vectors that need allocation on non-master ranks.
template<typename T>
void snew_bc(const t_commrec* cr, T*& data, int numElements)
{
    if (!MASTER(cr))
    {
        snew(data, numElements);
    }
}
//! Convenience wrapper for gmx_bcast of a C-style array which needs allocation on non-master ranks.
template<typename T>
void nblock_abc(const t_commrec* cr, int numElements, T** v)
{
    snew_bc(cr, v, numElements);
    nblock_bc(cr, numElements, *v);
}
//! Convenience wrapper for gmx_bcast of a std::vector which needs resizing on non-master ranks.
template<typename T>
void nblock_abc(const t_commrec* cr, int numElements, std::vector<T>* v)
{
    if (!MASTER(cr))
    {
        v->resize(numElements);
    }
    gmx_bcast(numElements * sizeof(T), v->data(), cr);
}

//! \brief Broadcasts the, non-dynamic, state from the master to all ranks in cr->mpi_comm_mygroup
//
// This is intended to be used with MPI parallelization without
// domain decompostion (currently with NM and TPI).
void broadcastStateWithoutDynamics(const t_commrec* cr, t_state* state);

//! \brief Broadcast inputrec and mtop and allocate node-specific settings
void init_parallel(t_commrec*                  cr,
                   t_inputrec*                 inputrec,
                   gmx_mtop_t*                 mtop,
                   PartialDeserializedTprFile* partialDeserializedTpr);

#endif
