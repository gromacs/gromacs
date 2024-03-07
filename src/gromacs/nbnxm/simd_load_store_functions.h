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
/*! \internal \file
 *
 * \brief
 * Defines functions to load data to and store data from SIMD registers
 * for the SIMD 4xM and 2xMM kernels.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_SIMD_LOAD_STORE_FUNCTIONS_H
#define GMX_NBNXM_SIMD_LOAD_STORE_FUNCTIONS_H

#include "config.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \brief Returns the j-cluster index for the given i-cluster index
 *
 * \tparam clusterRatio  The ratio of cluster size, supported are 0.5,1,2, checked at compile time
 * \param  iCluster      The index of the i-cluster
 * \returns the j-cluster index corresponding to \p iCluster
 */
template<KernelLayoutClusterRatio clusterRatio>
static inline int cjFromCi(const int iCluster)
{
    if constexpr (clusterRatio == KernelLayoutClusterRatio::JSizeEqualsISize)
    {
        return iCluster;
    }
    else if constexpr (clusterRatio == KernelLayoutClusterRatio::JSizeIsDoubleISize)
    {
        return (iCluster >> 1);
    }
    else if constexpr (clusterRatio == KernelLayoutClusterRatio::JSizeIsHalfISize)
    {
        return (iCluster << 1);
    }
    // Note that this code will refuse to compile when none of the branches above is taken
}

//! Load a single real for an i-atom into \p iRegister
template<KernelLayout kernelLayout>
inline std::enable_if_t<kernelLayout == KernelLayout::r4xM, SimdReal> loadIAtomData(const real* ptr,
                                                                                    const int offset,
                                                                                    const int iRegister)
{
    return SimdReal(ptr[offset + iRegister]);
}

//! Load a pair of consecutive reals for two i-atom into the respective halves of \p iRegister
template<KernelLayout kernelLayout>
inline std::enable_if_t<kernelLayout == KernelLayout::r2xMM, SimdReal>
loadIAtomData(const real* ptr, const int offset, const int iRegister)
{
    return loadU1DualHsimd(ptr + offset + iRegister * 2);
}

//! Returns a SIMD register containing GMX_SIMD_REAL_WIDTH reals loaded from ptr + offset
template<KernelLayout kernelLayout>
inline std::enable_if_t<kernelLayout == KernelLayout::r4xM, SimdReal> loadJAtomData(const real* ptr,
                                                                                    const int offset)
{
    return load<SimdReal>(ptr + offset);
}

//! Returns a SIMD register containing a duplicate sequence of GMX_SIMD_REAL_WIDTH/2 reals loaded from ptr + offset
template<KernelLayout kernelLayout>
inline std::enable_if_t<kernelLayout == KernelLayout::r2xMM, SimdReal> loadJAtomData(const real* ptr,
                                                                                     const int offset)
{
    return loadDuplicateHsimd(ptr + offset);
}

#if GMX_SIMD_HAVE_INT32_LOGICAL
//! Define SimdBitMask as an integer SIMD register
typedef SimdInt32 SimdBitMask;
#else
//! Define SimdBitMask as a real SIMD register
typedef SimdReal SimdBitMask;
#endif

//! Loads no interaction masks, returns an empty array
template<bool loadMasks, KernelLayout kernelLayout>
inline std::enable_if_t<!loadMasks, std::array<SimdBool, 0>>
loadSimdPairInteractionMasks(const int excl, SimdBitMask* filterBitMasksV)
{
    return std::array<SimdBool, 0>{};

    GMX_UNUSED_VALUE(excl);
    GMX_UNUSED_VALUE(filterBitMasksV);
}

//! Loads interaction masks for a cluster pair for 4xM kernel layout
template<bool loadMasks, KernelLayout kernelLayout>
inline std::enable_if_t<loadMasks && kernelLayout == KernelLayout::r4xM, std::array<SimdBool, sc_iClusterSize(kernelLayout)>>
loadSimdPairInteractionMasks(const int excl, SimdBitMask* filterBitMasksV)
{
    using namespace gmx;

    constexpr int c_iClusterSize = sc_iClusterSize(kernelLayout);

    std::array<SimdBool, c_iClusterSize> interactionMasksV;

#if GMX_SIMD_HAVE_INT32_LOGICAL
    /* Load integer interaction mask */
    SimdInt32 mask_pr_S(excl);
    for (int i = 0; i < c_iClusterSize; i++)
    {
        interactionMasksV[i] = cvtIB2B(testBits(mask_pr_S & filterBitMasksV[i]));
    }

#elif GMX_SIMD_HAVE_LOGICAL
    union
    {
#    if GMX_DOUBLE
        std::int64_t i;
#    else
        std::int32_t i;
#    endif
        real         r;
    } conv;

    conv.i = excl;
    SimdReal mask_pr_S(conv.r);

    for (int i = 0; i < c_iClusterSize; i++)
    {
        interactionMasksV[i] = testBits(mask_pr_S & filterBitMasksV[i]);
    }

#else
    // Note that there was exclusion support without logical support up to release-2022,
    // which could be resurrected. This is probably only efficient for 4-wide SIMD.
    static_assert(false,
                  "Need GMX_SIMD_HAVE_INT32_LOGICAL or GMX_SIMD_HAVE_LOGICAL for NBNxM kernels");
#endif

    return interactionMasksV;
}

//! Loads interaction masks for a cluster pair for 2xMM kernel layout
template<bool loadMasks, KernelLayout kernelLayout>
inline std::enable_if_t<loadMasks && kernelLayout == KernelLayout::r2xMM,
                        std::array<SimdBool, sc_iClusterSize(kernelLayout) / 2>>
loadSimdPairInteractionMasks(const int excl, SimdBitMask* filterBitMasksV)
{
    using namespace gmx;

    constexpr int c_iClusterSize = sc_iClusterSize(kernelLayout);

    std::array<SimdBool, c_iClusterSize / 2> interactionMasksV;

#if GMX_SIMD_HAVE_INT32_LOGICAL
    SimdInt32 mask_pr_S(excl);
    for (int i = 0; i < c_iClusterSize / 2; i++)
    {
        interactionMasksV[i] = cvtIB2B(testBits(mask_pr_S & filterBitMasksV[i]));
    }
#elif GMX_SIMD_HAVE_LOGICAL
    union
    {
#    if GMX_DOUBLE
        std::int64_t i;
#    else
        std::int32_t i;
#    endif
        real         r;
    } conv;

    conv.i = excl;
    SimdReal mask_pr_S(conv.r);

    for (int i = 0; i < c_iClusterSize / 2; i++)
    {
        interactionMasksV[i] = testBits(mask_pr_S & filterBitMasksV[i]);
    }
#else
    static_assert(false,
                  "Need GMX_SIMD_HAVE_INT32_LOGICAL or GMX_SIMD_HAVE_LOGICAL for NBNxM kernels");
#endif

    return interactionMasksV;
}

//! Return the number of atoms pairs that are within the cut-off distance
template<int nR>
inline int pairCountWithinCutoff(SimdReal rSquaredV[nR], SimdReal cutoffSquared)
{
    alignas(GMX_SIMD_ALIGNMENT) real tmp[GMX_SIMD_REAL_WIDTH];

    int npair = 0;
    for (int i = 0; i < nR; i++)
    {
        store(tmp, cutoffSquared - rSquaredV[i]);
        for (int j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            if (tmp[j] >= 0)
            {
                npair++;
            }
        }
    }

    return npair;
}

} // namespace gmx

#endif // GMX_NBNXM_SIMD_LOAD_STORE_FUNCTIONS_H
