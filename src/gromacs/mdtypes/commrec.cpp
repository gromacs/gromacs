/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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

#include "gromacs/mdtypes/commrec.h"

#include "config.h"

#include <utility>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"

t_commrec::t_commrec(const gmx::MpiComm& mpiComm) :
    commMySim(mpiComm), commMyGroup(mpiComm), mpiDefaultCommunicator(mpiComm), duty(DUTY_PP | DUTY_PME)
{
#if GMX_LIB_MPI
    GMX_RELEASE_ASSERT(gmx_mpi_initialized(), "Must have initialized MPI before building commrec");
#endif
}

t_commrec::~t_commrec()
{
#if GMX_MPI
    // TODO We need to be able to free communicators, but the
    // structure of the commrec and domdec initialization code makes
    // it hard to avoid both leaks and double frees.
    const bool mySimIsMyGroup = (commMySim.comm() == commMyGroup.comm());
    if (commMySim.comm() != MPI_COMM_NULL && commMySim.comm() != MPI_COMM_WORLD)
    {
        // TODO see above
        // MPI_Comm_free(&cr->mpi_comm_mysim);
    }
    if (!mySimIsMyGroup && commMyGroup.comm() != MPI_COMM_NULL && commMyGroup.comm() != MPI_COMM_WORLD)
    {
        // TODO see above
        // MPI_Comm_free(&cr->mpi_comm_mygroup);
    }
#endif
}

void t_commrec::setDD(std::unique_ptr<gmx_domdec_t>&& ddUniquePtr)
{
    ddUniquePtr_ = std::move(ddUniquePtr);
    dd           = ddUniquePtr_.get();
}

void t_commrec::destroyDD()
{
    ddUniquePtr_.reset(nullptr);
    dd = nullptr;
}
