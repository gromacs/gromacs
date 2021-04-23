/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
 * \brief SYCL-specific routines for the GPU implementation of SETTLE constraints algorithm.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */

#include "settle_gpu_internal.h"

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

void launchSettleGpuKernel(const int /* numSettles */,
                           const DeviceBuffer<WaterMolecule> /* d_atomIds */,
                           const SettleParameters /* settleParameters */,
                           const DeviceBuffer<Float3> /* d_x */,
                           DeviceBuffer<Float3> /* d_xp */,
                           const bool /* updateVelocities */,
                           DeviceBuffer<Float3> /* d_v */,
                           const real /* invdt */,
                           const bool /* computeVirial */,
                           DeviceBuffer<float> /* virialScaled */,
                           const PbcAiuc /* pbcAiuc */,
                           const DeviceStream& /* deviceStream */)
{
    // SYCL_TODO
    GMX_RELEASE_ASSERT(false, "SETTLE is not yet implemented in SYCL.");

    return;
}

} // namespace gmx
