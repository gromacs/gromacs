/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 *  \brief
 *  Stubs of functions that must be defined by nbnxm sycl implementation.
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/nbnxm/gpu_common.h"
#include "gromacs/utility/exceptions.h"

#include "nbnxm_sycl_types.h"

namespace Nbnxm
{

// SYCL-TODO: remove when functions are properly implemented
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-noreturn"

/*! \brief
 * Launch asynchronously the download of nonbonded forces from the GPU
 * (and energies/shift forces if required).
 */
void gpu_launch_cpyback(NbnxmGpu* /*nb*/,
                        struct nbnxn_atomdata_t* /*nbatom*/,
                        const gmx::StepWorkload& /*stepWork*/,
                        const AtomLocality /*atomLocality*/)
{
    GMX_THROW(gmx::NotImplementedError("Not implemented on SYCL yet"));
}

/*! \brief Launch asynchronously the xq buffer host to device copy. */
void gpu_copy_xq_to_gpu(NbnxmGpu* /*nb*/, const nbnxn_atomdata_t* /*nbatom*/, const AtomLocality /*atomLocality*/)
{
    GMX_THROW(gmx::NotImplementedError("Not implemented on SYCL yet"));
}

void gpu_launch_kernel_pruneonly(NbnxmGpu* /*nb*/, const InteractionLocality /*iloc*/, const int /*numParts*/)
{
    GMX_THROW(gmx::NotImplementedError("Not implemented on SYCL yet"));
}

void gpu_launch_kernel(NbnxmGpu* /*nb*/,
                       const gmx::StepWorkload& /*stepWork*/,
                       const Nbnxm::InteractionLocality /*iloc*/)
{
    GMX_THROW(gmx::NotImplementedError("Not implemented on SYCL yet"));
}

#pragma clang diagnostic pop // SYCL-TODO: remove when functions above are properly implemented

} // namespace Nbnxm
