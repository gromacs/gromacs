/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 *
 * Shared helpers for computing per-dimension/pulse halo plans
 */
#pragma once

#include "gromacs/utility/gmxmpi.h"

#include "domdec_struct.h"

namespace gmx
{

struct HaloPlan
{
    //! Number of halo zones associated with this halo exchange instance
    int numZone = 0;
    //! The atom offset for receive (x) or send (f) for dimension index and pulse corresponding to this halo exchange instance
    int atomOffset = 0;
    //! send copy size for X for this haloPlan
    int xSendSize = 0;
    //! recv copy size for X for this haloPlan
    int xRecvSize = 0;
    //! flag on whether the recieve for this halo exchange is performed in-place
    bool receiveInPlace = true;
    //! The indices to communicate for this halo exchange
    const gmx_domdec_ind_t* ind = nullptr;
};

HaloPlan computeHaloPlan(const gmx_domdec_comm_t& comm,
                         int                      dimIndex,
                         int                      pulse,
                         MPI_Comm                 mpiCommMySim,
                         int                      sendRankX,
                         int                      recvRankX);

} // namespace gmx
