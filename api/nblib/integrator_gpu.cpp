/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief
 * Implements nblib integrator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#include "nblib/integrator_gpu.h"

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/utility/arrayref.h"

#include "integrator_gpu_internal.h"

namespace nblib
{

LeapFrogGPU::LeapFrogGPU(int                  numAtoms,
                         gmx::ArrayRef<real>  inverseMasses,
                         const DeviceContext& deviceContext,
                         const DeviceStream&  deviceStream) :
    numAtoms_(numAtoms), deviceContext_(deviceContext), deviceStream_(deviceStream)
{
    reallocateDeviceBuffer(
            &inverseMasses_, numAtoms_, &numInverseMasses_, &numInverseMassesAlloc_, deviceContext_);
    copyToDeviceBuffer(
            &inverseMasses_, inverseMasses.data(), 0, numAtoms_, deviceStream_, GpuApiCallBehavior::Sync, nullptr);
}

LeapFrogGPU::~LeapFrogGPU()
{
    freeDeviceBuffer(&inverseMasses_);
}

void LeapFrogGPU::integrate(float                dt,
                            DeviceBuffer<Float4> coordinates,
                            DeviceBuffer<Float3> velocities,
                            DeviceBuffer<Float3> forces)
{
    launchLeapFrogKernel(numAtoms_, dt, coordinates, velocities, forces, inverseMasses_, deviceStream_);
}


} // namespace nblib
