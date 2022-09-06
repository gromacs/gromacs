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
/*! \internal \file
 *
 * \brief Implements backend-specific code for PME-PP communication using SYCL.
 *
 * Does not actually implement anything, since the only backend-specific part is for peer-to-peer
 * communication, which is not supported with SYCL yet.
 *
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "gromacs/utility/gmxassert.h"

#include "pme_force_sender_gpu_impl.h"

namespace gmx
{

/*! \brief Send PME synchronizer directly to the peer devices. Not implemented with SYCL. */
void PmeForceSenderGpu::Impl::sendFToPpPeerToPeer(int /*ppRank*/, int /*numAtoms*/, bool /*sendForcesDirectToPpGpu*/)
{
    GMX_RELEASE_ASSERT(false,
                       "Direct peer-to-peer communications not supported with SYCL and threadMPI.");
}

} // namespace gmx
