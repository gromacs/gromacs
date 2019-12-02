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
 * \brief Stub for update and constraints class CPU implementation.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/mdlib/update_constrain_cuda.h"

#if GMX_GPU != GMX_GPU_CUDA

namespace gmx
{

class UpdateConstrainCuda::Impl
{
};

UpdateConstrainCuda::UpdateConstrainCuda(const t_inputrec& /* ir   */,
                                         const gmx_mtop_t& /* mtop */,
                                         const void* /* commandStream */,
                                         GpuEventSynchronizer* /* xUpdatedOnDevice */) :
    impl_(nullptr)
{
    GMX_ASSERT(false,
               "A CPU stub for UpdateConstrain was called instead of the correct implementation.");
}

UpdateConstrainCuda::~UpdateConstrainCuda() = default;

void UpdateConstrainCuda::integrate(GpuEventSynchronizer* /* fReadyOnDevice */,
                                    const real /* dt */,
                                    const bool /* updateVelocities */,
                                    const bool /* computeVirial */,
                                    tensor /* virialScaled */,
                                    const bool /* doTemperatureScaling */,
                                    gmx::ArrayRef<const t_grp_tcstat> /* tcstat */,
                                    const bool /* doParrinelloRahman */,
                                    const float /* dtPressureCouple */,
                                    const matrix /* prVelocityScalingMatrix*/)
{
    GMX_ASSERT(false,
               "A CPU stub for UpdateConstrain was called instead of the correct implementation.");
}

void UpdateConstrainCuda::scaleCoordinates(const matrix /* scalingMatrix */)
{
    GMX_ASSERT(false,
               "A CPU stub for UpdateConstrain was called instead of the correct implementation.");
}

void UpdateConstrainCuda::set(DeviceBuffer<float> /* d_x */,
                              DeviceBuffer<float> /* d_v */,
                              const DeviceBuffer<float> /* d_f */,
                              const t_idef& /* idef */,
                              const t_mdatoms& /* md */,
                              const int /* numTempScaleValues */)
{
    GMX_ASSERT(false,
               "A CPU stub for UpdateConstrain was called instead of the correct implementation.");
}

void UpdateConstrainCuda::setPbc(const t_pbc* /* pbc */)
{
    GMX_ASSERT(false,
               "A CPU stub for UpdateConstrain was called instead of the correct implementation.");
}

GpuEventSynchronizer* UpdateConstrainCuda::getCoordinatesReadySync()
{
    GMX_ASSERT(false,
               "A CPU stub for UpdateConstrain was called instead of the correct implementation.");
    return nullptr;
}

bool UpdateConstrainCuda::isNumCoupledConstraintsSupported(const gmx_mtop_t& /* mtop */)
{
    return false;
}

} // namespace gmx

#endif /* GMX_GPU != GMX_GPU_CUDA */
