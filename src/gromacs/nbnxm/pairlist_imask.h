/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief
 * Declares and defines pairlist interaction mask generation functions
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_PAIRLIST_IMASK_H
#define GMX_NBNXM_PAIRLIST_IMASK_H

namespace gmx
{

/*! \brief Returns a diagonal interaction mask with atoms j<i masked out
 *
 * \tparam T             Integer type, should have at least iClusterSize*jClusterSize bits
 * \tparam iClusterSize  The i-cluster size
 * \tparam jClusterSize  The j-cluster size
 *
 * Condition: jClusterSize <= iClusterSize
 */
template<typename T, int iClusterSize, int jClusterSize>
constexpr std::array<T, iClusterSize / jClusterSize> diagonalMaskJSmallerI()
{
    static_assert(jClusterSize <= iClusterSize);
    static_assert(iClusterSize * jClusterSize <= sizeof(T) * CHAR_BIT);

    // Call the constructor to initialize the array to zero
    // to permit the function to be declared constexpr.
    std::array<T, iClusterSize / jClusterSize> mask{};

    for (int maskIndex = 0; maskIndex < iClusterSize / jClusterSize; maskIndex++)
    {
        for (int i = 0; i < iClusterSize; i++)
        {
            for (int j = std::max(i + 1 - maskIndex * jClusterSize, 0); j < jClusterSize; j++)
            {
                mask[maskIndex] |= (T(1) << (i * jClusterSize + j));
            }
        }
    }

    return mask;
}

/*! \brief Returns a diagonal interaction mask with atoms j>i masked out
 *
 * \tparam T             Integer type, should have at least iClusterSize*jClusterSize bits
 * \tparam iClusterSize  The i-cluster size
 * \tparam jClusterSize  The j-cluster size
 *
 * Condition: jClusterSize >= iClusterSize
 */
template<typename T, int iClusterSize, int jClusterSize>
constexpr std::array<T, jClusterSize / iClusterSize> diagonalMaskJLargerI()
{
    static_assert(jClusterSize >= iClusterSize);
    static_assert(iClusterSize * jClusterSize <= sizeof(T) * CHAR_BIT);

    // Call the constructor to initialize the array to zero
    // to permit the function to be declared constexpr.
    std::array<T, jClusterSize / iClusterSize> mask{};

    for (int maskIndex = 0; maskIndex < jClusterSize / iClusterSize; maskIndex++)
    {
        for (int i = 0; i < iClusterSize; i++)
        {
            for (int j = maskIndex * iClusterSize + i + 1; j < jClusterSize; j++)
            {
                mask[maskIndex] |= (T(1) << (i * jClusterSize + j));
            }
        }
    }

    return mask;
}

/*! \brief Returns a diagonal or off-diagonal interaction mask
 *
 * \tparam iClusterSize  The i-cluster size
 * \tparam jClusterSize  The j-cluster size
 * \param[in] maskOutSubDiagonal  Whether to mask out the sub-diagonal interactions
 * \param[in] ci                  The i-cluster index
 * \param[in] cj                  The j-cluster index
 */
template<int iClusterSize, int jClusterSize>
static gmx_unused uint32_t getImask(const bool maskOutSubDiagonal, const int ci, const int cj)
{
    if constexpr (jClusterSize >= iClusterSize)
    {
        // static, so the diagonal mask is only created once
        static constexpr auto sc_diagonalMask =
                diagonalMaskJLargerI<uint32_t, iClusterSize, jClusterSize>();

        constexpr int ratio = jClusterSize / iClusterSize;
        const int     diff  = ci - cj * ratio;

        return (maskOutSubDiagonal && diff >= 0 && diff < ratio ? sc_diagonalMask[diff]
                                                                : NBNXN_INTERACTION_MASK_ALL);
    }
    else
    {
        // static, so the diagonal mask is only created once
        static constexpr auto sc_diagonalMask =
                diagonalMaskJSmallerI<uint32_t, iClusterSize, jClusterSize>();

        constexpr int ratio = iClusterSize / jClusterSize;
        const int     diff  = cj - ci * ratio;

        return (maskOutSubDiagonal && diff >= 0 && diff < ratio ? sc_diagonalMask[diff]
                                                                : NBNXN_INTERACTION_MASK_ALL);
    }
}

} // namespace gmx

#endif // GMX_NBNXM_PAIRLIST_IMASK_H
