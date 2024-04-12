/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

/*! \libinternal \file
 * \brief Defines the host-side PME GPU data structures.
 * \todo Some renaming/refactoring, which does not impair the performance:
 * -- bringing the function names up to guidelines
 * -- PmeGpuSettings -> PmeGpuTasks
 * -- refining GPU notation application (#2053)
 * -- renaming coefficients to charges (?)
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_STAGING_H
#define GMX_EWALD_PME_GPU_STAGING_H

#include <array>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"

//! Number of FEP states
static constexpr int sc_numFepStates = 2;

/*! \internal \brief
 * The PME GPU intermediate buffers structure, included in the main PME GPU structure by value.
 * Buffers are managed by the PME GPU module.
 */
struct PmeGpuStaging
{
    //! Host-side force buffer
    gmx::PaddedHostVector<gmx::RVec> h_forces;

    /*! \brief Virial and energy intermediate host-side buffer. Size is PME_GPU_VIRIAL_AND_ENERGY_COUNT. */
    std::array<gmx::HostVector<float>, sc_numFepStates> h_virialAndEnergy;
    /*! \brief B-spline values intermediate host-side buffer. */
    std::array<gmx::HostVector<float>, sc_numFepStates> h_splineModuli;

    /*! \brief Pointer to the host memory with B-spline values. Only used for host-side gather, or unit tests */
    gmx::HostVector<float> h_theta;
    /*! \brief Pointer to the host memory with B-spline derivative values. Only used for host-side gather, or unit tests */
    gmx::HostVector<float> h_dtheta;
    /*! \brief Pointer to the host memory with ivec atom gridline indices. Only used for host-side gather, or unit tests */
    gmx::HostVector<int> h_gridlineIndices;
};

#endif
