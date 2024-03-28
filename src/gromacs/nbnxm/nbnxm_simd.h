/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * Defines constants used to know which nbNxM kernel flavours (4xM or 2xMM)
 * can be supported by the SIMD layer in use.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_NBNXM_SIMD_H
#define GMX_NBNXM_NBNXM_SIMD_H

#include "config.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/real.h"

#include "pairlist.h"
#include "pairlistparams.h"

//! The types of nbNxM SIMD kernel layout
enum class KernelLayout
{
    r4xM, //!< 4 'i'-registers each containing data for interaction with M j-atoms
    r2xMM //!< 2 'i'-registers each containing duplicated data, { M, M }, for interaction with M j-atoms
};

//! List of supported ratios for j-cluster size versus i-cluster sizes
enum class KernelLayoutClusterRatio
{
    JSizeEqualsISize,   //!< j-cluster size = i-cluster size
    JSizeIsDoubleISize, //!< j-cluster size = 2 * i-cluster size
    JSizeIsHalfISize    //!< j-cluster size = i-cluster size / 2
};

//! The NBNxM i-cluster size in atoms for the given NBNxM kernel layout
static constexpr int sc_iClusterSize(const KernelLayout kernelLayout)
{
    switch (kernelLayout)
    {
        case KernelLayout::r4xM: return 4;
        case KernelLayout::r2xMM: return 4;
    }
}

#if GMX_SIMD

//! The NBNxM j-cluster size in atoms for the given NBNxM kernel layout
static constexpr int sc_jClusterSize(const KernelLayout kernelLayout)
{
    switch (kernelLayout)
    {
        case KernelLayout::r4xM: return GMX_SIMD_REAL_WIDTH;
        case KernelLayout::r2xMM: return GMX_SIMD_REAL_WIDTH / 2;
    }
}

/*! \brief Returns the ratio of j-cluster size versus i-cluster size
 *
 * \tparam kernelLayout The kernel layout, supports r4xM and r2xMM, asserted at compile time
 *
 * Note that currently only cluster ratios 0.5, 1 and 2 are supported. This is checked
 * at compile time.
 */
template<KernelLayout kernelLayout>
static constexpr KernelLayoutClusterRatio kernelLayoutClusterRatio()
{
    constexpr int c_iClusterSize = sc_iClusterSize(kernelLayout);
    constexpr int c_jClusterSize = sc_jClusterSize(kernelLayout);

    static_assert(c_jClusterSize == c_iClusterSize || c_jClusterSize == 2 * c_iClusterSize
                  || 2 * c_jClusterSize == c_iClusterSize);

    if constexpr (c_jClusterSize == c_iClusterSize)
    {
        return KernelLayoutClusterRatio::JSizeEqualsISize;
    }
    else if constexpr (c_jClusterSize == 2 * c_iClusterSize)
    {
        return KernelLayoutClusterRatio::JSizeIsDoubleISize;
    }
    else if constexpr (2 * c_jClusterSize == c_iClusterSize)
    {
        return KernelLayoutClusterRatio::JSizeIsHalfISize;
    }
    // Note that this code will refuse to compile when none of the branches above is taken
}

#endif

/*! \brief The nbnxn SIMD 4xN and 2x(N+N) kernels can be added independently.
 * Currently the 2xNN SIMD kernels only make sense with:
 *  8-way SIMD: 4x4 setup, performance wise only useful on CPUs without FMA or on AMD Zen1
 * 16-way SIMD: 4x8 setup, used in single precision with 512 bit wide SIMD
 */
#if GMX_SIMD && GMX_USE_SIMD_KERNELS
#    define GMX_HAVE_NBNXM_SIMD_2XMM \
        ((GMX_SIMD_REAL_WIDTH == 8 || GMX_SIMD_REAL_WIDTH == 16) && GMX_SIMD_HAVE_HSIMD_UTIL_REAL)
#    define GMX_HAVE_NBNXM_SIMD_4XM \
        (GMX_SIMD_REAL_WIDTH == 2 || GMX_SIMD_REAL_WIDTH == 4 || GMX_SIMD_REAL_WIDTH == 8)
#else
#    define GMX_HAVE_NBNXM_SIMD_2XMM 0
#    define GMX_HAVE_NBNXM_SIMD_4XM 0
#endif

//! Whether we have support for NBNxM 2xM kernels
static constexpr bool sc_haveNbnxmSimd2xmmKernels = GMX_HAVE_NBNXM_SIMD_2XMM;
//! Whether we have support for NBNxM 4xM kernels
static constexpr bool sc_haveNbnxmSimd4xmKernels = GMX_HAVE_NBNXM_SIMD_4XM;

#if GMX_SIMD && GMX_USE_SIMD_KERNELS
// We use the FDV0 tables for width==4 (when we can load it in one go), or if we don't have any unaligned loads
#    if GMX_SIMD_REAL_WIDTH == 4 || !GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL
static constexpr bool c_useTableFormatFDV0 = true;
#    else
static constexpr bool c_useTableFormatFDV0 = false;
#    endif

/*! \brief Helper function to generate an std::array of SimdReal using a lambda
 *
 * An example of intended usage is:
 * auto rV = genArr<4>([&](int i) { return rSquaredV[i] * rInvV[i]; });
 *
 * \tparam N  The size of the array
 * \tparam F  The function type, deduced
 * \param[in] f  The function
 */
template<int N, class F>
std::array<gmx::SimdReal, N> genArr(F f)
{
    static_assert(N <= 2 || N == 4);

    if constexpr (N == 0)
    {
        return std::array<gmx::SimdReal, 0>{};
    }
    else if constexpr (N == 1)
    {
        return std::array<gmx::SimdReal, 1>{ f(0) };
    }
    else if constexpr (N == 2)
    {
        return std::array<gmx::SimdReal, 2>{ f(0), f(1) };
    }
    else
    {
        return std::array<gmx::SimdReal, 4>{ f(0), f(1), f(2), f(3) };
    }
}

/*! \brief Helper function to generate an std::array of SimdBool using a lambda
 *
 * An example of intended usage is:
 * auto withinCutoffV = genBoolArr<nR>([&](int i) { return rSquaredV[i] < rc2_S; });
 *
 * \tparam N  The size of the array
 * \tparam F  The function type, deduced
 * \param[in] f  The function
 */
template<int N, class F>
std::array<gmx::SimdBool, N> genBoolArr(F f)
{
    static_assert(N <= 2 || N == 4);

    if constexpr (N == 0)
    {
        return std::array<gmx::SimdBool, 0>{};
    }
    else if constexpr (N == 1)
    {
        return std::array<gmx::SimdBool, 1>{ f(0) };
    }
    else if constexpr (N == 2)
    {
        return std::array<gmx::SimdBool, 2>{ f(0), f(1) };
    }
    else
    {
        return std::array<gmx::SimdBool, 4>{ f(0), f(1), f(2), f(3) };
    }
}

/*! \brief Returns the sum over an array of SimdReal elements
 *
 * \tparam N  The size of the array
 * \tparam a  The elements to sum over
 */
template<size_t N>
gmx::SimdReal sumArray(const std::array<gmx::SimdReal, N>& a)
{
    static_assert(N == 2 || N == 4);

    if constexpr (N == 2)
    {
        return a[0] + a[1];
    }
    else
    {
        return a[0] + a[1] + a[2] + a[3];
    }
}

#endif // GMX_SIMD && GMX_USE_SIMD_KERNELS

#endif // GMX_NBNXM_NBNXM_SIMD
