/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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

/*! \internal
 * \ingroup __module_nbnxm
 *
 * \brief Defines a class for masking (sub-) diagonal interactions in an cluster pair
 *
 * Declares and defines a templated class DiagonalMasker with a maskArray function
 * which sets diagonal and sub-diagonal entries of a cluster-pair mask array to false
 * when the i- and j-cluster index combination falls on the diagonal.
 *
 * \author Berk Hess <hess@kth.se>
 */

#ifndef GMX_NBNXM_SIMD_DIAGONAL_MASKER_H
#define GMX_NBNXM_SIMD_DIAGONAL_MASKER_H

#include "gromacs/simd/simd.h"

#include "atomdata.h"
#include "nbnxm_simd.h"

namespace gmx
{

//! Base Coulomb calculator class, only specializations are used
template<int, KernelLayout, KernelLayoutClusterRatio>
class DiagonalMasker;

//! Returns the diagonal filter masks
template<int nR, KernelLayout kernelLayout>
inline std::array<std::array<SimdBool, nR>,
                  kernelLayoutClusterRatio<kernelLayout>() == KernelLayoutClusterRatio::JSizeEqualsISize ? 1 : 2>
generateDiagonalMasks(const nbnxn_atomdata_t::SimdMasks& simdMasks)
{
    constexpr KernelLayoutClusterRatio clusterRatio = kernelLayoutClusterRatio<kernelLayout>();

    /* Load j-i for the first i */
    SimdReal diagonalJMinusI = load<SimdReal>(kernelLayout == KernelLayout::r4xM
                                                      ? simdMasks.diagonal_4xn_j_minus_i.data()
                                                      : simdMasks.diagonal_2xnn_j_minus_i.data());

    /* With 2xMM layout there are two i-clusters in one register */
    const SimdReal iIndexIncrement(kernelLayout == KernelLayout::r4xM ? 1 : 2);
    const SimdReal zero(0.0_real);
    /* Generate all the diagonal masks as comparison results */
    std::array<std::array<SimdBool, nR>, clusterRatio == KernelLayoutClusterRatio::JSizeEqualsISize ? 1 : 2> diagonalMaskVV;
    for (int i = 0; i < nR; i++)
    {
        diagonalMaskVV[0][i] = (zero < diagonalJMinusI);
        diagonalJMinusI      = diagonalJMinusI - iIndexIncrement;
    }
    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    if constexpr (clusterRatio != KernelLayoutClusterRatio::JSizeEqualsISize)
    {
        if (clusterRatio == KernelLayoutClusterRatio::JSizeIsHalfISize)
        {
            /* Load j-i for the second half of the j-cluster */
            diagonalJMinusI = load<SimdReal>(simdMasks.diagonal_4xn_j_minus_i.data() + nR / 2);
        }
        for (int i = 0; i < nR; i++)
        {
            diagonalMaskVV[1][i] = (zero < diagonalJMinusI);
            diagonalJMinusI      = diagonalJMinusI - iIndexIncrement;
        }
    }
    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    return diagonalMaskVV;
}

//! Specialized masker for JSizeEqualsISize
template<int nR, KernelLayout kernelLayout>
class DiagonalMasker<nR, kernelLayout, KernelLayoutClusterRatio::JSizeEqualsISize>
{
public:
    inline DiagonalMasker(const nbnxn_atomdata_t::SimdMasks& simdMasks) :
        diagonalMaskV_(generateDiagonalMasks<nR, kernelLayout>(simdMasks)[0])
    {
    }

    //! Sets (sub-)diagonal entries in \p boolV to false when the cluster pair in on the diagonal
    inline void maskArray(const int iClusterIndex, const int jClusterIndex, std::array<SimdBool, nR>& boolV) const
    {
        if (jClusterIndex == iClusterIndex)
        {
            boolV = genBoolArr<nR>([&](int i) { return boolV[i] && diagonalMaskV_[i]; });
        }
    }

private:
    //! The diagonal mask array
    const std::array<SimdBool, nR> diagonalMaskV_;
};

//! Specialized masker for JSizeIsDoubleISize
template<int nR, KernelLayout kernelLayout>
class DiagonalMasker<nR, kernelLayout, KernelLayoutClusterRatio::JSizeIsDoubleISize>
{
public:
    inline DiagonalMasker(const nbnxn_atomdata_t::SimdMasks& simdMasks) :
        diagonalMaskVV_(generateDiagonalMasks<nR, kernelLayout>(simdMasks))
    {
    }

    //! Sets (sub-)diagonal entries in \p boolV to false when the cluster pair in on the diagonal
    inline void maskArray(const int iClusterIndex, const int jClusterIndex, std::array<SimdBool, nR>& boolV) const
    {
        if (jClusterIndex * 2 == iClusterIndex)
        {
            boolV = genBoolArr<nR>([&](int i) { return boolV[i] && diagonalMaskVV_[0][i]; });
        }
        else if (jClusterIndex * 2 + 1 == iClusterIndex)
        {
            boolV = genBoolArr<nR>([&](int i) { return boolV[i] && diagonalMaskVV_[1][i]; });
        }
    }

private:
    /*! \brief The diagonal mask array for:
     * j-cluster index * 2 = i-cluster index
     * j-cluster index * 2 + 1 = i-cluster index
     */
    const std::array<std::array<SimdBool, nR>, 2> diagonalMaskVV_;
};

//! Specialized masker for JSizeIsHalfISize
template<int nR, KernelLayout kernelLayout>
class DiagonalMasker<nR, kernelLayout, KernelLayoutClusterRatio::JSizeIsHalfISize>
{
public:
    inline DiagonalMasker(const nbnxn_atomdata_t::SimdMasks& simdMasks) :
        diagonalMaskVV_(generateDiagonalMasks<nR, kernelLayout>(simdMasks))
    {
    }

    //! Sets (sub-)diagonal entries in \p boolV to false when the cluster pair in on the diagonal
    inline void maskArray(const int iClusterIndex, const int jClusterIndex, std::array<SimdBool, nR>& boolV) const
    {
        if (jClusterIndex == iClusterIndex * 2)
        {
            boolV = genBoolArr<nR>([&](int i) { return boolV[i] && diagonalMaskVV_[0][i]; });
        }
        else if (jClusterIndex == iClusterIndex * 2 + 1)
        {
            boolV = genBoolArr<nR>([&](int i) { return boolV[i] && diagonalMaskVV_[1][i]; });
        }
    }

private:
    /*! \brief The diagonal mask array for:
     * j-cluster index = i-cluster index * 2
     * j-cluster index = i-cluster index * 2 + 1
     */
    const std::array<std::array<SimdBool, nR>, 2> diagonalMaskVV_;
};

} // namespace gmx

#endif // GMX_NBNXM_SIMD_DIAGONAL_MASKER_H
