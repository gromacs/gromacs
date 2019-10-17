/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Declaration of class which sends PME Force from GPU memory to PP task
 *
 * \author Alan Gray <alang@nvidia.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */
#ifndef GMX_PMEFORCESENDERGPU_H
#define GMX_PMEFORCESENDERGPU_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxmpi.h"

/*! \libinternal
 * \brief Contains information about the PP ranks that partner this PME rank. */
struct PpRanks
{
    //! The MPI rank ID of this partner PP rank.
    int rankId = -1;
    //! The number of atoms to communicate with this partner PP rank.
    int numAtoms = -1;
};

namespace gmx
{

template<typename T>
class ArrayRef;

/*! \libinternal
 * \brief Manages sending forces from PME-only ranks to their PP ranks. */
class PmeForceSenderGpu
{

public:
    /*! \brief Creates PME GPU Force sender object
     * \param[in] pmeStream       CUDA stream used for PME computations
     * \param[in] comm            Communicator used for simulation
     * \param[in] ppRanks         List of PP ranks
     */
    PmeForceSenderGpu(void* pmeStream, MPI_Comm comm, gmx::ArrayRef<PpRanks> ppRanks);
    ~PmeForceSenderGpu();

    /*! \brief
     * Initialization of GPU PME Force sender
     * \param[in] d_f   force buffer in GPU memory
     */
    void sendForceBufferAddressToPpRanks(rvec* d_f);

    /*! \brief
     * Send PP data to PP rank
     * \param[in] ppRank           PP rank to receive data
     */
    void sendFToPpCudaDirect(int ppRank);

private:
    class Impl;
    gmx::PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
