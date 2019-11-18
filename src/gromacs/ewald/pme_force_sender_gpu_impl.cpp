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
/*! \internal \file
 *
 * \brief May be used to implement PME-PP GPU comm interfaces for non-GPU builds.
 *
 * Currently, reports and exits if any of the interfaces are called.
 * Needed to satisfy compiler on systems, where CUDA is not available.
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/ewald/pme_force_sender_gpu.h"
#include "gromacs/utility/arrayref.h"

#if GMX_GPU != GMX_GPU_CUDA

namespace gmx
{

/*!\brief Impl class stub. */
class PmeForceSenderGpu::Impl
{
};

/*!\brief Constructor stub. */
PmeForceSenderGpu::PmeForceSenderGpu(void gmx_unused* pmeStream,
                                     MPI_Comm gmx_unused    comm,
                                     gmx::ArrayRef<PpRanks> gmx_unused ppRanks) :
    impl_(nullptr)
{
    GMX_ASSERT(false,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

PmeForceSenderGpu::~PmeForceSenderGpu() = default;

/*!\brief init PME-PP GPU communication stub */
void PmeForceSenderGpu::sendForceBufferAddressToPpRanks(rvec gmx_unused* d_f)
{
    GMX_ASSERT(false,
               "A CPU stub for PME-PP GPU communication initialization was called instead of the "
               "correct implementation.");
}

void PmeForceSenderGpu::sendFToPpCudaDirect(int gmx_unused ppRank)
{
    GMX_ASSERT(false,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

} // namespace gmx

#endif /* GMX_GPU != GMX_GPU_CUDA */
