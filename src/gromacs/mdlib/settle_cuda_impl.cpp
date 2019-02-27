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
#include "gmxpre.h"

#include "config.h"

#include "gromacs/mdlib/settle_cuda.h"

#if GMX_GPU != GMX_GPU_CUDA

/*!\brief Impl class stub */
class SettleCuda::Impl
{
};

/*!\brief Constructor stub */
SettleCuda::SettleCuda(const int        gmx_unused  nAtom,
                       const gmx_mtop_t gmx_unused &mtop)
    : impl_(nullptr)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

SettleCuda::SettleCuda(const int  gmx_unused nAtom,
                       const real gmx_unused mO,  const real gmx_unused mH,
                       const real gmx_unused dOH, const real gmx_unused dHH)
    : impl_(nullptr)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

/*!\brief Destructor stub */
SettleCuda::~SettleCuda() = default;

/*!\brief Apply LINCS stub */
void SettleCuda::apply(const bool gmx_unused  updateVelocities,
                       const real gmx_unused  invdt,
                       const bool gmx_unused  computeVirial,
                       tensor     gmx_unused  virialScaled)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

/*!\brief Set data structures stub */
void SettleCuda::set(const t_idef    gmx_unused  &idef,
                     const t_mdatoms gmx_unused  &md)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

/*!\brief Set PBC stub */
void SettleCuda::setPbc(t_pbc gmx_unused *pbc)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

/*! \brief
 * Copy coordinates from provided CPU location to GPU stub.
 *
 * \param[in] *x  CPU pointer where coordinates should be copied from.
 * \param[in] *xp CPU pointer where coordinates should be copied from.
 */
void SettleCuda::copyCoordinatesToGpu(rvec gmx_unused *x, rvec gmx_unused *xp)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

/*! \brief
 * Copy velocities from provided CPU location to GPU stub.
 *
 * \param[in] *v  CPU pointer where velocities should be copied from.
 */
void SettleCuda::copyVelocitiesToGpu(rvec gmx_unused *v)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

/*! \brief
 * Copy coordinates from GPU to provided CPU location stub.
 *
 * \param[out] *xp CPU pointer where coordinates should be copied to.
 */
void SettleCuda::copyCoordinatesFromGpu(rvec gmx_unused *xp)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

/*! \brief
 * Copy velocities from GPU to provided CPU location stub.
 *
 */
void SettleCuda::copyVelocitiesFromGpu(rvec gmx_unused *v)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

/*! \brief
 * Set the internal GPU-memory x, xprime and v pointers stub.
 *
 * \param[in] *xpDevice Pointer to the coordinates after integrator update, before update (on GPU)
 * \param[in] *vDevice  Pointer to the velocities before integrator update (on GPU)
 */
void SettleCuda::setXVPointers(rvec gmx_unused *xDevice, rvec gmx_unused *xpDevice, rvec gmx_unused *vDevice)
{
    GMX_ASSERT(false, "A CPU stub for SETTLE was called insted of the correct implementation.");
}

#endif /* GMX_GPU != GMX_GPU_CUDA */
