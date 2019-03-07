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
 * \brief May be used to implement Leap-Frog CUDA interfaces for non-GPU builds.
 *
 * Currently, reports and exits if any of the interfaces are called.
 * Needed to satisfy compiler on systems, where CUDA is not available.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/mdlib/leapfrog_cuda.h"

#if GMX_GPU != GMX_GPU_CUDA

namespace gmx
{

/*!\brief Impl class stub. */
class LeapFrogCuda::Impl
{
};

/*!\brief Constructor stub. */
LeapFrogCuda::LeapFrogCuda(gmx_unused int numAtoms)
    : impl_(nullptr)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

LeapFrogCuda::~LeapFrogCuda() = default;

/*!\brief Apply LINCS stub. */
void LeapFrogCuda::integrate(gmx_unused const real  dt)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*!\brief Set PBC stub. */
void LeapFrogCuda::setPbc(gmx_unused const t_pbc *pbc)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*! \brief Copy inverse masses from CPU to GPU stub function. */
void LeapFrogCuda::set(gmx_unused const t_mdatoms &md)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*! \brief Copy coordinates from provided CPU location to GPU stub. */
void LeapFrogCuda::copyCoordinatesToGpu(gmx_unused const rvec *h_x)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*! \brief Copy velocities from provided CPU location to GPU stub. */
void LeapFrogCuda::copyVelocitiesToGpu(gmx_unused const rvec *h_v)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*! \brief Copy forces from CPU to GPU stub. */
void LeapFrogCuda::copyForcesToGpu(gmx_unused const rvec *h_f)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*! \brief Copy coordinates from GPU to provided CPU location stub. */
void LeapFrogCuda::copyCoordinatesFromGpu(gmx_unused rvec *h_xp)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*! \brief Copy velocities from GPU to provided CPU location stub. */
void LeapFrogCuda::copyVelocitiesFromGpu(gmx_unused rvec *h_v)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*! \brief Copy forces from GPU to CPU stub. */
void LeapFrogCuda::copyForcesFromGpu(gmx_unused rvec *h_f)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

/*! \brief Set the internal GPU-memory x, xprime and v pointers stub. */
void LeapFrogCuda::setXVFPointers(gmx_unused rvec *d_x,
                                  gmx_unused rvec *d_xp,
                                  gmx_unused rvec *d_v,
                                  gmx_unused rvec *d_f)
{
    GMX_ASSERT(false, "A CPU stub for LeapFrog was called insted of the correct implementation.");
}

}      // namespace gmx

#endif /* GMX_GPU != GMX_GPU_CUDA */
