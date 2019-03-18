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
#include "gromacs/gpu_utils/devicebuffer.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/leapfrog_cuda.h"
#include "gromacs/mdlib/lincs_cuda.h"
#include "gromacs/mdlib/settle_cuda.h"
#include "gromacs/mdlib/update_constrain_cuda.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"

namespace gmx
{

/*! \brief Integrate
 *
 * Integrates the equation of motion using Leap-Frog algorithm and applies
 * LINCS and SETTLE constraints.
 * Updates d_xp_ and d_v_ fields of this object.
 *
 * \param[in] dt                Timestep
 * \param[in] updateVelocities  If the velocities should be constrained.
 * \param[in] computeVirial     If virial should be updated.
 * \param[out] virial           Place to save virial tensor.
 */
void UpdateConstrainCuda::Impl::integrate(const real  dt,
                                          const bool  updateVelocities,
                                          const bool  computeVirial,
                                          tensor      virial)
{
    // Clearing virial matrix
    // TODO There is no point in having saparate virial matrix for constraints
    clear_mat(virial);

    integrator_->integrate(dt);
    lincsCuda_->apply(updateVelocities, 1.0/dt, computeVirial, virial);
    settleCuda_->apply(updateVelocities, 1.0/dt, computeVirial, virial);

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

/*! \brief Create Update-Constrain object
 *
 * \param[in] numAtoms  Number of atoms.
 * \param[in] ir        Input record data: LINCS takes number of iterations and order of
 *                      projection from it.
 * \param[in] mtop      Topology of the system: SETTLE gets the masses for O and H atoms
 *                      and target O-H and H-H distances from this object.
 */
UpdateConstrainCuda::Impl::Impl(int                numAtoms,
                                const t_inputrec  &ir,
                                const gmx_mtop_t  &mtop)
    : numAtoms_(numAtoms)
{
    allocateDeviceBuffer(&d_x_,              numAtoms, nullptr);
    allocateDeviceBuffer(&d_xp_,             numAtoms, nullptr);
    allocateDeviceBuffer(&d_v_,              numAtoms, nullptr);
    allocateDeviceBuffer(&d_f_,              numAtoms, nullptr);
    allocateDeviceBuffer(&d_inverseMasses_,  numAtoms, nullptr);

    // TODO When the code will be integrated into the schedule, it will be assigned non-default stream.
    stream_ = nullptr;

    GMX_RELEASE_ASSERT(numAtoms == mtop.natoms, "State and topology number of atoms should be the same.");
    integrator_ = std::make_unique<LeapFrogCuda>(numAtoms);
    lincsCuda_  = std::make_unique<LincsCuda>(mtop.natoms, ir.nLincsIter, ir.nProjOrder);
    settleCuda_ = std::make_unique<SettleCuda>(mtop.natoms, mtop);

    integrator_->setXVFPointers((rvec*)d_x_, (rvec*)d_xp_, (rvec*)d_v_, (rvec*)d_f_);
    lincsCuda_->setXVPointers((rvec*)d_x_, (rvec*)d_xp_, (rvec*)d_v_);
    settleCuda_->setXVPointers((rvec*)d_x_, (rvec*)d_xp_,  (rvec*)d_v_);
}

UpdateConstrainCuda::Impl::~Impl()
{
}

/*! \brief
 * Update data-structures (e.g. after NB search step).
 *
 * \param[in] idef    System topology
 * \param[in] md      Atoms data.
 */
void UpdateConstrainCuda::Impl::set(const t_idef      &idef,
                                    const t_mdatoms   &md)
{
    // Integrator should also update something, but it does not even have a method yet
    integrator_->set(md);
    lincsCuda_->set(idef, md);
    settleCuda_->set(idef, md);
}

/*! \brief
 * Update PBC data.
 *
 * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
 *
 * \param[in] pbc The PBC data in t_pbc format.
 */
void UpdateConstrainCuda::Impl::setPbc(const t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc_);
    integrator_->setPbc(pbc);
    lincsCuda_->setPbc(pbc);
    settleCuda_->setPbc(pbc);
}

/*! \brief
 * Copy coordinates from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] h_x  CPU pointer where coordinates should be copied from.
 */
void UpdateConstrainCuda::Impl::copyCoordinatesToGpu(const rvec *h_x)
{
    copyToDeviceBuffer(&d_x_, (float3*)h_x, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy velocities from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] h_v  CPU pointer where velocities should be copied from.
 */
void UpdateConstrainCuda::Impl::copyVelocitiesToGpu(const rvec *h_v)
{
    copyToDeviceBuffer(&d_v_, (float3*)h_v, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy forces from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] h_f  CPU pointer where forces should be copied from.
 */
void UpdateConstrainCuda::Impl::copyForcesToGpu(const rvec *h_f)
{
    copyToDeviceBuffer(&d_f_, (float3*)h_f, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy coordinates from GPU to CPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[out] h_xp CPU pointer where coordinates should be copied to.
 */
void UpdateConstrainCuda::Impl::copyCoordinatesFromGpu(rvec *h_xp)
{
    copyFromDeviceBuffer((float3*)h_xp, &d_xp_, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy velocities from GPU to CPU.
 *
 * The velocities are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] h_v  Pointer to velocities data.
 */
void UpdateConstrainCuda::Impl::copyVelocitiesFromGpu(rvec *h_v)
{
    copyFromDeviceBuffer((float3*)h_v, &d_v_, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy forces from GPU to CPU.
 *
 * The forces are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] h_f  Pointer to forces data.
 */
void UpdateConstrainCuda::Impl::copyForcesFromGpu(rvec *h_f)
{
    copyFromDeviceBuffer((float3*)h_f, &d_f_, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Set the internal GPU-memory x, xprime and v pointers.
 *
 * Data is not copied. The data are assumed to be in float3/fvec format
 * (float3 is used internally, but the data layout should be identical).
 *
 * \param[in] d_x  Pointer to the coordinates for the input (on GPU)
 * \param[in] d_xp Pointer to the coordinates for the output (on GPU)
 * \param[in] d_v  Pointer to the velocities (on GPU)
 * \param[in] d_f  Pointer to the forces (on GPU)
 */
void UpdateConstrainCuda::Impl::setXVFPointers(rvec *d_x, rvec *d_xp, rvec *d_v, rvec *d_f)
{
    d_x_  = (float3*)d_x;
    d_xp_ = (float3*)d_xp;
    d_v_  = (float3*)d_v;
    d_f_  = (float3*)d_f;
}


UpdateConstrainCuda::UpdateConstrainCuda(int                numAtoms,
                                         const t_inputrec  &ir,
                                         const gmx_mtop_t  &mtop)
    : impl_(new Impl(numAtoms, ir, mtop))
{
}

UpdateConstrainCuda::~UpdateConstrainCuda() = default;

void UpdateConstrainCuda::integrate(const real  dt,
                                    const bool  updateVelocities,
                                    const bool  computeVirial,
                                    tensor      virialScaled)
{
    impl_->integrate(dt, updateVelocities, computeVirial, virialScaled);
}

void UpdateConstrainCuda::set(const t_idef               &idef,
                              const t_mdatoms gmx_unused &md)
{
    impl_->set(idef, md);
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
