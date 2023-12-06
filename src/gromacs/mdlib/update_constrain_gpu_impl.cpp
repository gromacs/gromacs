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
 * \brief Implements update and constraints class.
 *
 * The class combines Leap-Frog integrator with LINCS and SETTLE constraints.
 *
 * \todo The computational procedures in members should be integrated to improve
 *       computational performance.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "update_constrain_gpu_impl.h"

#include <cassert>
#include <cmath>
#include <cstdio>

#include <algorithm>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/mdlib/leapfrog_gpu.h"
#include "gromacs/mdlib/update_constrain_gpu.h"
#include "gromacs/mdlib/update_constrain_gpu_internal.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"

static constexpr bool sc_haveGpuConstraintSupport = (GMX_GPU_CUDA != 0) || (GMX_GPU_SYCL != 0);

namespace gmx
{

void UpdateConstrainGpu::Impl::integrate(GpuEventSynchronizer*             fReadyOnDevice,
                                         const real                        dt,
                                         const bool                        updateVelocities,
                                         const bool                        computeVirial,
                                         tensor                            virial,
                                         const bool                        doTemperatureScaling,
                                         gmx::ArrayRef<const t_grp_tcstat> tcstat,
                                         const bool                        doParrinelloRahman,
                                         const float                       dtPressureCouple,
                                         const Matrix3x3&                  prVelocityScalingMatrix)
{
    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);

    // Clearing virial matrix
    // TODO There is no point in having separate virial matrix for constraints
    clear_mat(virial);

    // Make sure that the forces are ready on device before proceeding with the update.
    fReadyOnDevice->enqueueWaitEvent(deviceStream_);

    if (numAtoms_ != 0)
    {
        // A copy of the current coordinates is saved into d_x0_ by integrate(), and
        // d_x_ is updated by integration and constraints.
        integrator_->integrate(
                d_x_, d_x0_, d_v_, d_f_, dt, doTemperatureScaling, tcstat, doParrinelloRahman, dtPressureCouple, prVelocityScalingMatrix);
        if (sc_haveGpuConstraintSupport)
        {
            lincsGpu_->apply(d_x0_, d_x_, updateVelocities, d_v_, 1.0 / dt, computeVirial, virial, pbcAiuc_);
            settleGpu_->apply(d_x0_, d_x_, updateVelocities, d_v_, 1.0 / dt, computeVirial, virial, pbcAiuc_);
        }

        // scaledVirial -> virial (methods above returns scaled values)
        float scaleFactor = 0.5F / (dt * dt);
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                virial[i][j] = scaleFactor * virial[i][j];
            }
        }
    }

    xUpdatedOnDeviceEvent_.markEvent(deviceStream_);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void UpdateConstrainGpu::Impl::scaleCoordinates(const Matrix3x3& scalingMatrix)
{
    if (numAtoms_ == 0)
    {
        return;
    }

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);

    ScalingMatrix mu(scalingMatrix);

    launchScaleCoordinatesKernel(numAtoms_, d_x_, mu, deviceStream_);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void UpdateConstrainGpu::Impl::scaleVelocities(const Matrix3x3& scalingMatrix)
{
    if (numAtoms_ == 0)
    {
        return;
    }

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);

    ScalingMatrix mu(scalingMatrix);

    launchScaleCoordinatesKernel(numAtoms_, d_v_, mu, deviceStream_);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

UpdateConstrainGpu::Impl::Impl(const t_inputrec&    ir,
                               const gmx_mtop_t&    mtop,
                               const int            numTempScaleValues,
                               const DeviceContext& deviceContext,
                               const DeviceStream&  deviceStream,
                               gmx_wallcycle*       wcycle) :
    deviceContext_(deviceContext), deviceStream_(deviceStream), wcycle_(wcycle)
{
    integrator_ = std::make_unique<LeapFrogGpu>(deviceContext_, deviceStream_, numTempScaleValues);
    if (sc_haveGpuConstraintSupport)
    {
        lincsGpu_ = std::make_unique<LincsGpu>(ir.nLincsIter, ir.nProjOrder, deviceContext_, deviceStream_);
        settleGpu_ = std::make_unique<SettleGpu>(mtop, deviceContext_, deviceStream_);
    }
}

UpdateConstrainGpu::Impl::~Impl()
{
    freeDeviceBuffer(&d_x0_);
    freeDeviceBuffer(&d_inverseMasses_);
}

void UpdateConstrainGpu::Impl::set(DeviceBuffer<Float3>          d_x,
                                   DeviceBuffer<Float3>          d_v,
                                   const DeviceBuffer<Float3>    d_f,
                                   const InteractionDefinitions& idef,
                                   const t_mdatoms&              md)
{
    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);

    GMX_ASSERT(d_x, "Coordinates device buffer should not be null.");
    GMX_ASSERT(d_v, "Velocities device buffer should not be null.");
    GMX_ASSERT(d_f, "Forces device buffer should not be null.");

    d_x_ = d_x;
    d_v_ = d_v;
    d_f_ = d_f;

    numAtoms_ = md.homenr;

    reallocateDeviceBuffer(&d_x0_, numAtoms_, &numXp_, &numXpAlloc_, deviceContext_);

    reallocateDeviceBuffer(
            &d_inverseMasses_, numAtoms_, &numInverseMasses_, &numInverseMassesAlloc_, deviceContext_);

    // Integrator should also update something, but it does not even have a method yet
    integrator_->set(numAtoms_, md.invmass, md.cTC);
    if (sc_haveGpuConstraintSupport)
    {
        lincsGpu_->set(idef, numAtoms_, md.invmass);
        settleGpu_->set(idef);
    }
    else
    {
        GMX_ASSERT(idef.il[F_SETTLE].empty(), "SETTLE not supported");
        GMX_ASSERT(idef.il[F_CONSTR].empty(), "LINCS not supported");
    }

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void UpdateConstrainGpu::Impl::setPbc(const PbcType pbcType, const matrix box)
{
    // TODO wallcycle
    setPbcAiuc(numPbcDimensions(pbcType), box, &pbcAiuc_);
}

GpuEventSynchronizer* UpdateConstrainGpu::Impl::xUpdatedOnDeviceEvent()
{
    return &xUpdatedOnDeviceEvent_;
}

UpdateConstrainGpu::UpdateConstrainGpu(const t_inputrec&    ir,
                                       const gmx_mtop_t&    mtop,
                                       const int            numTempScaleValues,
                                       const DeviceContext& deviceContext,
                                       const DeviceStream&  deviceStream,
                                       gmx_wallcycle*       wcycle) :
    impl_(new Impl(ir, mtop, numTempScaleValues, deviceContext, deviceStream, wcycle))
{
}

UpdateConstrainGpu::~UpdateConstrainGpu() = default;

void UpdateConstrainGpu::integrate(GpuEventSynchronizer*             fReadyOnDevice,
                                   const real                        dt,
                                   const bool                        updateVelocities,
                                   const bool                        computeVirial,
                                   tensor                            virialScaled,
                                   const bool                        doTemperatureScaling,
                                   gmx::ArrayRef<const t_grp_tcstat> tcstat,
                                   const bool                        doParrinelloRahman,
                                   const float                       dtPressureCouple,
                                   const gmx::Matrix3x3&             prVelocityScalingMatrix)
{
    impl_->integrate(fReadyOnDevice,
                     dt,
                     updateVelocities,
                     computeVirial,
                     virialScaled,
                     doTemperatureScaling,
                     tcstat,
                     doParrinelloRahman,
                     dtPressureCouple,
                     prVelocityScalingMatrix);
}

void UpdateConstrainGpu::scaleCoordinates(const gmx::Matrix3x3& scalingMatrix)
{
    impl_->scaleCoordinates(scalingMatrix);
}

void UpdateConstrainGpu::scaleVelocities(const gmx::Matrix3x3& scalingMatrix)
{
    impl_->scaleVelocities(scalingMatrix);
}

void UpdateConstrainGpu::set(DeviceBuffer<Float3>          d_x,
                             DeviceBuffer<Float3>          d_v,
                             const DeviceBuffer<Float3>    d_f,
                             const InteractionDefinitions& idef,
                             const t_mdatoms&              md)
{
    impl_->set(d_x, d_v, d_f, idef, md);
}

void UpdateConstrainGpu::setPbc(const PbcType pbcType, const matrix box)
{
    impl_->setPbc(pbcType, box);
}

GpuEventSynchronizer* UpdateConstrainGpu::xUpdatedOnDeviceEvent()
{
    return impl_->xUpdatedOnDeviceEvent();
}

bool UpdateConstrainGpu::isNumCoupledConstraintsSupported(const gmx_mtop_t& mtop)
{
    return LincsGpu::isNumCoupledConstraintsSupported(mtop);
}

bool UpdateConstrainGpu::areConstraintsSupported()
{
    return sc_haveGpuConstraintSupport;
}

} // namespace gmx
