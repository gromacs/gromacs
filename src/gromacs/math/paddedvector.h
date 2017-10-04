/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::PaddedRVecVector
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
 * \ingroup module_math
 */
#ifndef GMX_MATH_PADDEDVECTOR_H
#define GMX_MATH_PADDEDVECTOR_H

#include <algorithm>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"

namespace gmx
{

/*! \brief Temporary definition of a type usable for SIMD-style loads of RVec quantities.
 *
 * \note When resizing paddedRVecVector, the size should be chosen with
         paddedRVecVectorSize() to ensure correct padding.
 * \todo Consider replacing the padding applied in resizePaddedRVecVector()
 *       by automated padding on resize() of the vector.
 * \todo Undo the move of allocator.h and alignedallocator.h from the internal
 *       to be public API applied in Change-Id: Ifb8dacf, needed to use
 *       AlginedAllocationPolicy here, when replacing std::vector here.
 */
using PaddedRVecVector = std::vector < gmx::RVec, gmx::Allocator < gmx::RVec, AlignedAllocationPolicy > >;

/*! \brief Returns the padded size for PaddedRVecVector given the number of atoms
 *
 * \param[in] numAtoms  The number of atoms for which data will be stored in a PaddedRVecVector
 */
static inline size_t paddedRVecVectorSize(size_t numAtoms)
{
    /* We need one real extra for 4-wide SIMD load/store of RVec.
     * But because the vector contains RVecs, we need to add 1 RVec.
     */
    size_t simdScatterAccessSize;
    if (numAtoms > 0)
    {
        simdScatterAccessSize = numAtoms + 1;
    }
    else
    {
        simdScatterAccessSize = 0;
    }

    /* We need padding up to a multiple of the SIMD real width.
     * But because we don't want a dependence on the SIMD module,
     * we use GMX_REAL_MAX_SIMD_WIDTH here.
     */
    size_t simdFlatAccesSize  = ((numAtoms + (GMX_REAL_MAX_SIMD_WIDTH - 1))/GMX_REAL_MAX_SIMD_WIDTH)*GMX_REAL_MAX_SIMD_WIDTH;

    return std::max(simdScatterAccessSize, simdFlatAccesSize);
}

} // namespace gmx

// TODO This is a hack to avoid littering gmx:: all over code that is
// almost all destined to move into the gmx namespace at some point.
// An alternative would be about 20 files with using statements.
using gmx::PaddedRVecVector;

#endif
