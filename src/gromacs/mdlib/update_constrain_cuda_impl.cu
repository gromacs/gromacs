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
 * \brief Implements update and constraints class using CUDA.
 *
 * The class combines Leap-Frog integrator with LINCS and SETTLE constraints.
 *
 * \todo This class should take over the management of coordinates, velocities
 *       forces, virial, and PBC from its members (i.e. from Leap-Frog, LINCS
 *       and SETTLE).
 * \todo The computational procedures in members should be integrated to improve
 *       computational performance.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "update_constrain_cuda_impl.h"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/mdlib/leapfrog_cuda.cuh"
#include "gromacs/mdlib/lincs_cuda.cuh"
#include "gromacs/mdlib/settle_cuda.cuh"
#include "gromacs/mdlib/update_constrain_cuda.h"

namespace gmx
{

void UpdateConstrainCuda::Impl::integrate(const real                        dt,
                                          const bool                        updateVelocities,
                                          const bool                        computeVirial,
                                          tensor                            virial,
                                          const bool                        doTempCouple,
                                          gmx::ArrayRef<const t_grp_tcstat> tcstat,
                                          const bool                        doPressureCouple,
                                          const float                       dtPressureCouple,
                                          const matrix                      velocityScalingMatrix)
{
    // Clearing virial matrix
    // TODO There is no point in having separate virial matrix for constraints
    clear_mat(virial);

    integrator_->integrate(d_x_, d_xp_, d_v_, d_f_, dt,
                           doTempCouple, tcstat,
                           doPressureCouple, dtPressureCouple, velocityScalingMatrix);
    lincsCuda_->apply(d_x_, d_xp_,
                      updateVelocities, d_v_, 1.0/dt,
                      computeVirial, virial);
    settleCuda_->apply(d_x_, d_xp_,
                       updateVelocities, d_v_, 1.0/dt,
                       computeVirial, virial);

    // scaledVirial -> virial (methods above returns scaled values)
    float scaleFactor = 0.5f/(dt*dt);
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            virial[i][j] = scaleFactor*virial[i][j];
        }
    }

    return;
}

UpdateConstrainCuda::Impl::Impl(const t_inputrec  &ir,
                                const gmx_mtop_t  &mtop)
{
    // TODO When the code will be integrated into the schedule, it will be assigned non-default stream.
    stream_ = nullptr;

    integrator_ = std::make_unique<LeapFrogCuda>();
    lincsCuda_  = std::make_unique<LincsCuda>(ir.nLincsIter, ir.nProjOrder);
    settleCuda_ = std::make_unique<SettleCuda>(mtop);

}

UpdateConstrainCuda::Impl::~Impl()
{
}

void UpdateConstrainCuda::Impl::set(const t_idef    &idef,
                                    const t_mdatoms &md,
                                    const int        numTempScaleValues)
{
    numAtoms_ = md.nr;

    reallocateDeviceBuffer(&d_x_,  numAtoms_, &numX_,  &numXAlloc_,  nullptr);
    reallocateDeviceBuffer(&d_xp_, numAtoms_, &numXp_, &numXpAlloc_, nullptr);
    reallocateDeviceBuffer(&d_v_,  numAtoms_, &numV_,  &numVAlloc_,  nullptr);
    reallocateDeviceBuffer(&d_f_,  numAtoms_, &numF_,  &numFAlloc_,  nullptr);

    reallocateDeviceBuffer(&d_inverseMasses_, numAtoms_,
                           &numInverseMasses_, &numInverseMassesAlloc_, nullptr);

    // Integrator should also update something, but it does not even have a method yet
    integrator_->set(md, numTempScaleValues, md.cTC);
    lincsCuda_->set(idef, md);
    settleCuda_->set(idef, md);
}

void UpdateConstrainCuda::Impl::setPbc(const t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc_);
    integrator_->setPbc(pbc);
    lincsCuda_->setPbc(pbc);
    settleCuda_->setPbc(pbc);
}

void UpdateConstrainCuda::Impl::copyCoordinatesToGpu(const rvec *h_x)
{
    copyToDeviceBuffer(&d_x_, (float3*)h_x, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

void UpdateConstrainCuda::Impl::copyVelocitiesToGpu(const rvec *h_v)
{
    copyToDeviceBuffer(&d_v_, (float3*)h_v, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

void UpdateConstrainCuda::Impl::copyForcesToGpu(const rvec *h_f)
{
    copyToDeviceBuffer(&d_f_, (float3*)h_f, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

void UpdateConstrainCuda::Impl::copyCoordinatesFromGpu(rvec *h_xp)
{
    copyFromDeviceBuffer((float3*)h_xp, &d_xp_, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

void UpdateConstrainCuda::Impl::copyVelocitiesFromGpu(rvec *h_v)
{
    copyFromDeviceBuffer((float3*)h_v, &d_v_, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

void UpdateConstrainCuda::Impl::copyForcesFromGpu(rvec *h_f)
{
    copyFromDeviceBuffer((float3*)h_f, &d_f_, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

void UpdateConstrainCuda::Impl::setXVFPointers(rvec *d_x, rvec *d_xp, rvec *d_v, rvec *d_f)
{
    d_x_  = (float3*)d_x;
    d_xp_ = (float3*)d_xp;
    d_v_  = (float3*)d_v;
    d_f_  = (float3*)d_f;
}


UpdateConstrainCuda::UpdateConstrainCuda(const t_inputrec  &ir,
                                         const gmx_mtop_t  &mtop)
    : impl_(new Impl(ir, mtop))
{
}

UpdateConstrainCuda::~UpdateConstrainCuda() = default;

void UpdateConstrainCuda::integrate(const real                        dt,
                                    const bool                        updateVelocities,
                                    const bool                        computeVirial,
                                    tensor                            virialScaled,
                                    const bool                        doTempCouple,
                                    gmx::ArrayRef<const t_grp_tcstat> tcstat,
                                    const bool                        doPressureCouple,
                                    const float                       dtPressureCouple,
                                    const matrix                      pRVScalingMatrix)
{
    impl_->integrate(dt, updateVelocities, computeVirial, virialScaled,
                     doTempCouple, tcstat,
                     doPressureCouple, dtPressureCouple, pRVScalingMatrix);
}

void UpdateConstrainCuda::set(const t_idef    &idef,
                              const t_mdatoms &md,
                              const int        numTempScaleValues)
{
    impl_->set(idef, md, numTempScaleValues);
}

void UpdateConstrainCuda::setPbc(const t_pbc *pbc)
{
    impl_->setPbc(pbc);
}

void UpdateConstrainCuda::copyCoordinatesToGpu(const rvec *h_x)
{
    impl_->copyCoordinatesToGpu(h_x);
}

void UpdateConstrainCuda::copyVelocitiesToGpu(const rvec *h_v)
{
    impl_->copyVelocitiesToGpu(h_v);
}

void UpdateConstrainCuda::copyForcesToGpu(const rvec *h_f)
{
    impl_->copyForcesToGpu(h_f);
}

void UpdateConstrainCuda::copyCoordinatesFromGpu(rvec *h_xp)
{
    impl_->copyCoordinatesFromGpu(h_xp);
}

void UpdateConstrainCuda::copyVelocitiesFromGpu(rvec *h_v)
{
    impl_->copyVelocitiesFromGpu(h_v);
}

void UpdateConstrainCuda::copyForcesFromGpu(rvec *h_f)
{
    impl_->copyForcesFromGpu(h_f);
}

void UpdateConstrainCuda::setXVFPointers(rvec *d_x, rvec *d_xp, rvec *d_v, rvec *d_f)
{
    impl_->setXVFPointers(d_x, d_xp, d_v, d_f);
}

} //namespace gmx
