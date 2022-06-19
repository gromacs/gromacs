/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Declares functionality for communicators across physical nodes.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_MDTYPES_PHYSICALNODECOMMUNICATOR_H
#define GMX_MDTYPES_PHYSICALNODECOMMUNICATOR_H

#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

/*! \brief Wrapper function for RAII-style cleanup.
 *
 * This is needed to discard the return value so it can be used as a
 * deleter by a smart pointer. */
void MPI_Comm_free_wrapper(MPI_Comm* comm);

//! Make a smart pointer for MPI communicators.
using MPI_Comm_ptr = gmx::unique_cptr<MPI_Comm, MPI_Comm_free_wrapper>;

/*! \libinternal \brief Holds a communicator for the physical node of this rank
 *
 * This communicator should only be used for appropriate tasks,
 * e.g. during initialization and finalization. It can contain ranks
 * from PP, PME and multiple simulations with multisim, so is not
 * suited for general-purpose communication. */
class PhysicalNodeCommunicator
{
public:
    /*! \brief Constructor.
     *
     * Communicates within \c world to make intra-communicator \c
     * comm_ between all ranks that share \c physicalNodeId. */
    PhysicalNodeCommunicator(MPI_Comm world, int physicalNodeId);
    //! Communicator for all ranks on this physical node
    MPI_Comm comm_;
    //! Number of ranks on this physical node, corresponds to MPI_Comm_size of comm.
    int size_;
    //! Rank ID within this physical node, corresponds to MPI_Comm_rank of comm.
    int rank_;
    //! RAII handler for cleaning up \c comm_ only when appropriate.
    MPI_Comm_ptr commGuard_;
    //! Creates a barrier for all ranks on this physical node.
    void barrier() const;
};

} // namespace gmx

#endif
