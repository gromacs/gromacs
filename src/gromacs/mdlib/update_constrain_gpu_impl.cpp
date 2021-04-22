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
 * \brief Implements update and constraints class.
 *
 * The class combines Leap-Frog or SD integrator with LINCS and SETTLE constraints.
 *
 * \todo The computational procedures in members should be integrated to improve
 *       computational performance.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "update_constrain_gpu_impl.h"

#include <cassert>
#include <cmath>
#include <cstdio>

#include <algorithm>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/mdlib/langevin_gpu.h"
#include "gromacs/mdlib/leapfrog_gpu.h"
#include "gromacs/mdlib/update_constrain_gpu.h"
#include "gromacs/mdlib/update_constrain_gpu_internal.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"

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
                                         const Matrix3x3&                  prVelocityScalingMatrix,
                                         const int                         seed,
                                         const int64_t                     step)
{
    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuUpdateConstrain);

    GMX_ASSERT(integratorLeapFrog_ == nullptr || integratorLangevin_ == nullptr,
               "Cannot use Leap frog and Langevin (SD) integrators simultaneously.");

    // Clearing virial matrix
    // TODO There is no point in having separate virial matrix for constraints
    clear_mat(virial);

    // Make sure that the forces are ready on device before proceeding with the update.
    fReadyOnDevice->enqueueWaitEvent(deviceStream_);

    if (numAtoms_ != 0)
    {
        // A copy of the current coordinates is saved into d_x0_ by integrate(), and
        // d_x_ is updated by integration and constraints.
        if (integratorLeapFrog_ != nullptr)
        {
            integratorLeapFrog_->integrate(
                    d_x_, d_x0_, d_v_, d_f_, dt, doTemperatureScaling, tcstat, doParrinelloRahman, dtPressureCouple, prVelocityScalingMatrix);
        }
        else if (integratorLangevin_ != nullptr)
        {
            integratorLangevin_->integrate(d_x_, d_x0_, d_v_, d_f_, dt, seed, step, SDUpdate::ForcesOnly);
        }
        // Constraints need both coordinates before (d_x_) and after (d_xp_) update. However, after constraints
        // are applied, the d_x_ can be discarded. So we intentionally swap the d_x_ and d_xp_ here to avoid the
        // d_xp_ -> d_x_ copy after constraints. Note that the integrate saves them in the wrong order as well.
        if constexpr (GpuConfigurationCapabilities::Constraints)
        {
            lincsGpu_->apply(d_x0_, d_x_, updateVelocities, d_v_, 1.0 / dt, computeVirial, virial, pbcAiuc_);
            settleGpu_->apply(d_x0_, d_x_, updateVelocities, d_v_, 1.0 / dt, computeVirial, virial, pbcAiuc_);
        }

        if (integratorLangevin_ != nullptr)
        {
            integratorLangevin_->integrate(
                    d_x_, d_x0_, d_v_, d_f_, dt, seed, step, SDUpdate::FrictionAndNoiseOnly);
            /* Constrain the coordinates upd->x0 for half a time step */
            const bool computeVirialAtHalfTimeStep    = false;
            const bool updateVelocitiesAtHalfTimeStep = false;
            if constexpr (GpuConfigurationCapabilities::Constraints)
            {
                lincsGpu_->apply(d_x0_,
                                 d_x_,
                                 updateVelocitiesAtHalfTimeStep,
                                 nullptr,
                                 1.0 / (0.5 * dt),
                                 computeVirialAtHalfTimeStep,
                                 nullptr,
                                 pbcAiuc_);
                settleGpu_->apply(d_x0_,
                                  d_x_,
                                  updateVelocitiesAtHalfTimeStep,
                                  nullptr,
                                  1.0 / (0.5 * dt),
                                  computeVirialAtHalfTimeStep,
                                  nullptr,
                                  pbcAiuc_);
            }
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
    if (ir.eI == IntegrationAlgorithm::SD1)
    {
        GMX_ASSERT(GpuConfigurationCapabilities::UpdateSD && GpuConfigurationCapabilities::Constraints,
                   "GPU constraint support is required when using the GPU for updates with the SD "
                   "integrator.");
        integratorLangevin_ = std::make_unique<LangevinGpu>(
                deviceContext_, deviceStream_, numTempScaleValues, ir.delta_t, ir.opts.ref_t, ir.opts.tau_t);
        integratorLeapFrog_ = nullptr;
    }
    else
    {
        GMX_ASSERT(GpuConfigurationCapabilities::UpdateLeapfrog,
                   "GPU update for leapfrog needs to be available to run on the device");
        integratorLeapFrog_ =
                std::make_unique<LeapFrogGpu>(deviceContext_, deviceStream_, numTempScaleValues);
        integratorLangevin_ = nullptr;
    }
    if constexpr (GpuConfigurationCapabilities::Constraints)
    {
        lincsGpu_ = std::make_unique<LincsGpu>(ir.nLincsIter, ir.nProjOrder, deviceContext_, deviceStream_);
        settleGpu_ = std::make_unique<SettleGpu>(mtop, deviceContext_, deviceStream_);
    }
}

UpdateConstrainGpu::Impl::~Impl()
{
    try
    {
        freeDeviceBuffer(&d_x0_);
        freeDeviceBuffer(&d_inverseMasses_);
    }
    catch (gmx::InternalError& e)
    {
        fprintf(stderr, "Internal error in destructor of UpdateConstrainGpu: %s\n", e.what());
    }
}

void UpdateConstrainGpu::Impl::set(DeviceBuffer<Float3>          d_x,
                                   DeviceBuffer<Float3>          d_v,
                                   const DeviceBuffer<Float3>    d_f,
                                   const InteractionDefinitions& idef,
                                   const t_mdatoms&              md)
{
    wallcycle_start(wcycle_, WallCycleCounter::GpuSetConstr);

    GMX_ASSERT(d_x, "Coordinates device buffer should not be null.");
    GMX_ASSERT(d_v, "Velocities device buffer should not be null.");
    GMX_ASSERT(d_f, "Forces device buffer should not be null.");
    GMX_ASSERT(integratorLeapFrog_ == nullptr || integratorLangevin_ == nullptr,
               "Cannot use Leap frog and Langevin (SD) integrators simultaneously.");


    d_x_ = d_x;
    d_v_ = d_v;
    d_f_ = d_f;

    numAtoms_ = md.homenr;

    reallocateDeviceBuffer(&d_x0_, numAtoms_, &numXp_, &numXpAlloc_, deviceContext_);

    reallocateDeviceBuffer(
            &d_inverseMasses_, numAtoms_, &numInverseMasses_, &numInverseMassesAlloc_, deviceContext_);

    // Integrator should also update something, but it does not even have a method yet
    if (integratorLeapFrog_ != nullptr)
    {
        integratorLeapFrog_->set(numAtoms_, md.invmass, md.cTC);
    }
    else if (integratorLangevin_ != nullptr)
    {
        integratorLangevin_->set(numAtoms_, md.invmass, md.cTC);
    }
    if constexpr (GpuConfigurationCapabilities::Constraints)
    {
        wallcycle_sub_start(wcycle_, WallCycleSubCounter::GpuSetLincs);
        lincsGpu_->set(idef, numAtoms_, md.invmass);
        wallcycle_sub_stop(wcycle_, WallCycleSubCounter::GpuSetLincs);
        wallcycle_sub_start(wcycle_, WallCycleSubCounter::GpuSetSettle);
        settleGpu_->set(idef);
        wallcycle_sub_stop(wcycle_, WallCycleSubCounter::GpuSetSettle);
    }
    else
    {
        GMX_ASSERT(idef.il[InteractionFunction::SETTLE].empty(), "SETTLE not supported");
        GMX_ASSERT(idef.il[InteractionFunction::Constraints].empty(), "LINCS not supported");
    }

    wallcycle_stop(wcycle_, WallCycleCounter::GpuSetConstr);
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
                                       const int            numTempCouplGroups,
                                       const DeviceContext& deviceContext,
                                       const DeviceStream&  deviceStream,
                                       gmx_wallcycle*       wcycle) :
    impl_(new Impl(ir, mtop, numTempCouplGroups, deviceContext, deviceStream, wcycle))
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
                                   const gmx::Matrix3x3&             prVelocityScalingMatrix,
                                   const int                         seed,
                                   const int64_t                     step)
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
                     prVelocityScalingMatrix,
                     seed,
                     step);
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
    return GpuConfigurationCapabilities::Constraints;
}

} // namespace gmx
