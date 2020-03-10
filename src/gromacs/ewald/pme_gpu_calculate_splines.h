/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_EWALD_PME_GPU_UTILS_H
#define GMX_EWALD_PME_GPU_UTILS_H

/*! \internal \file
 * \brief This file defines small PME GPU inline host/device functions.
 * Note that OpenCL device-side functions can't use C++ features, so they are
 * located in a similar file pme_gpu_utils.clh.
 * Be sure to keep the logic in sync in both files when changing it!
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include <cassert>

#include "gromacs/math/vectypes.h"

struct PmeGpu;

/*! \internal \brief
 * Gets a base of the unique index to an element in a spline parameter buffer (theta/dtheta),
 * which is laid out for GPU spread/gather kernels. The base only corresponds to the atom index within the execution block.
 * Feed the result into getSplineParamIndex() to get a full index.
 * TODO: it's likely that both parameters can be just replaced with a single atom index, as they are derived from it.
 * Do that, verifying that the generated code is not bloated, and/or revise the spline indexing scheme.
 * Removing warp dependency would also be nice (and would probably coincide with removing c_pmeSpreadGatherAtomsPerWarp).
 *
 * \tparam order               PME order
 * \tparam atomsPerWarp        Number of atoms processed by a warp
 * \param[in] warpIndex        Warp index wrt the block.
 * \param[in] atomWarpIndex    Atom index wrt the warp (from 0 to atomsPerWarp - 1).
 *
 * \returns Index into theta or dtheta array using GPU layout.
 */
template<int order, int atomsPerWarp>
int inline getSplineParamIndexBase(int warpIndex, int atomWarpIndex)
{
    assert((atomWarpIndex >= 0) && (atomWarpIndex < atomsPerWarp));
    const int dimIndex    = 0;
    const int splineIndex = 0;
    // The zeroes are here to preserve the full index formula for reference
    return (((splineIndex + order * warpIndex) * DIM + dimIndex) * atomsPerWarp + atomWarpIndex);
}

/*! \internal \brief
 * Gets a unique index to an element in a spline parameter buffer (theta/dtheta),
 * which is laid out for GPU spread/gather kernels. The index is wrt to the execution block,
 * in range(0, atomsPerBlock * order * DIM).
 * This function consumes result of getSplineParamIndexBase() and adjusts it for \p dimIndex and \p splineIndex.
 *
 * \tparam order               PME order
 * \tparam atomsPerWarp        Number of atoms processed by a warp
 * \param[in] paramIndexBase   Must be result of getSplineParamIndexBase().
 * \param[in] dimIndex         Dimension index (from 0 to 2)
 * \param[in] splineIndex      Spline contribution index (from 0 to \p order - 1)
 *
 * \returns Index into theta or dtheta array using GPU layout.
 */
template<int order, int atomsPerWarp>
int inline getSplineParamIndex(int paramIndexBase, int dimIndex, int splineIndex)
{
    assert((dimIndex >= XX) && (dimIndex < DIM));
    assert((splineIndex >= 0) && (splineIndex < order));
    return (paramIndexBase + (splineIndex * DIM + dimIndex) * atomsPerWarp);
}

#endif
