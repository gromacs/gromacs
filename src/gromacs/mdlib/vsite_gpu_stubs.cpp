/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Stub implementation of the VsiteGpu class for the HIP and SYCL builds.
 *
 * GPU virtual site construction and force spreading are currently only
 * implemented for CUDA (see vsite_gpu.cpp and vsite_gpu_internal.cu). These
 * stubs provide the VsiteGpu symbols required by update_constrain_gpu_impl.cpp
 * so the HIP and SYCL builds can link.
 *
 * \todo Remove this stub once the HIP and SYCL builds have their own proper implementations.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/gmxassert.h"

#include "vsite_gpu.h"

namespace gmx
{

VsiteGpu::VsiteGpu(const gmx_mtop_t& /* mtop */,
                   const DeviceContext& deviceContext,
                   const DeviceStream&  deviceStream) :
    deviceContext_(deviceContext), deviceStream_(deviceStream)
{
}

VsiteGpu::~VsiteGpu() = default;

void VsiteGpu::set(const InteractionDefinitions& idef)
{
    int numVsites = 0;
    for (const auto& type : c_vsiteTypesOnGpu)
    {
        numVsites += idef.il[type].size();
    }
    GMX_RELEASE_ASSERT(numVsites == 0,
                       "GPU virtual site handling is only implemented for CUDA, but a topology "
                       "with virtual sites reached the VsiteGpu stub.");
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void VsiteGpu::construct(DeviceBuffer<Float3> /* d_x */, const PbcAiuc& /* pbcAiuc */) {}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void VsiteGpu::spreadForces(DeviceBuffer<Float3> /* d_f */) const {}

} // namespace gmx
