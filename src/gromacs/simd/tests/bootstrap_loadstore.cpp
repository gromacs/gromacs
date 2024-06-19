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
#include "gmxpre.h"

/*! \internal \file
 * \brief
 * Separate test of SIMD load/store, before we use them in the SIMD test classes.
 *
 * Simple tests without using any classes/utilities, so we can use load/store
 * functions inside our test utilities after this has passed.
 *
 * This file tests both the aligned and (if available) unaligned load and store
 * operatations for SimdReal, SimdInt32, and Simd4Real.
 *
 * Note that you probably do not have to add more tests in this (complicated)
 * file; once the bootstrapping tests have passed we can use the working basic
 * load/store operations to test higher-level load/store operations too.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 * \ingroup module_simd
 */

#include "config.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#if GMX_SIMD

namespace gmx
{

namespace test
{

namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */


/*! \brief Generic routine to test load & store of SIMD, and check for side effects.
 *
 * The tests for load, store, unaligned load and unaligned store both for
 * real and int are pretty much similar, so we use a template function with
 * additional function pointers for the actual load/store calls.
 */
template<typename T, typename TSimd, int simdWidth>
void loadStoreTester(TSimd gmx_simdcall loadFn(const T* mem),
                     void gmx_simdcall  storeFn(T* mem, TSimd),
                     const int          loadOffset,
                     const int          storeOffset)
{
    /* We need simdWidth storage in the first place, another simdWidth elements
     * so we can create (deliberately) offset un-aligned pointers, and finally
     * simdWidth elements at the beginning and end
     * to test we are not polluting memory there either. Sum=4*simdWidth.
     */
    alignas(GMX_SIMD_ALIGNMENT) T src[simdWidth * 4];
    alignas(GMX_SIMD_ALIGNMENT) T dst[simdWidth * 4];

    // Make sure we have memory to check both before and after the test pointers
    T*  pCopySrc = src + simdWidth + loadOffset;
    T*  pCopyDst = dst + simdWidth + storeOffset;
    int i;

    for (i = 0; i < simdWidth * 4; i++)
    {
        src[i] = 1 + i;
        dst[i] = -1 - i;
    }

    storeFn(pCopyDst, loadFn(pCopySrc));

    for (i = 0; i < simdWidth; i++)
    {
        EXPECT_EQ(pCopySrc[i], pCopyDst[i])
                << "SIMD load or store not moving data correctly for element " << i;
    }

    for (i = 0; i < simdWidth * 4; i++)
    {
        EXPECT_EQ(src[i], T(1 + i)) << "Side effect on source memory, i = " << i;
        if (dst + i < pCopyDst || dst + i >= pCopyDst + simdWidth)
        {
            EXPECT_EQ(dst[i], T(-1 - i)) << "Side effect on destination memory, i = " << i;
        }
    }
}

/*! \brief Wrapper to handle proxy objects returned by some load functions.
 *
 * \tparam T      Type of scalar object
 * \tparam TSimd  Corresponding SIMD type
 * \param  m      Memory address to load from
 */
template<typename T, typename TSimd>
TSimd gmx_simdcall loadWrapper(const T* m)
{
    return load<TSimd>(m);
}

/*! \brief Wrapper to handle proxy objects returned by some loadU functions.
 *
 * \tparam T      Type of scalar object
 * \tparam TSimd  Corresponding SIMD type
 * \param  m      Memory address to load from
 */
template<typename T, typename TSimd>
TSimd gmx_simdcall loadUWrapper(const T* m)
{
    return loadU<TSimd>(m);
}


#    if GMX_SIMD_HAVE_REAL
TEST(SimdBootstrapTest, loadStore)
{
    loadStoreTester<real, SimdReal, GMX_SIMD_REAL_WIDTH>(loadWrapper, store, 0, 0);
}

#        if GMX_SIMD_HAVE_LOADU
TEST(SimdBootstrapTest, loadU)
{
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        loadStoreTester<real, SimdReal, GMX_SIMD_REAL_WIDTH>(loadUWrapper, store, i, 0);
    }
}
#        endif // GMX_SIMD_HAVE_LOADU

#        if GMX_SIMD_HAVE_STOREU
TEST(SimdBootstrapTest, storeU)
{
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        loadStoreTester<real, SimdReal, GMX_SIMD_REAL_WIDTH>(loadWrapper, storeU, 0, i);
    }
}
#        endif // GMX_SIMD_HAVE_STOREU

// Tests for SimdInt32 load & store operations
TEST(SimdBootstrapTest, loadStoreI)
{
    loadStoreTester<int, SimdInt32, GMX_SIMD_REAL_WIDTH>(loadWrapper, store, 0, 0);
}

#        if GMX_SIMD_HAVE_LOADU
TEST(SimdBootstrapTest, loadUI)
{
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        loadStoreTester<int, SimdInt32, GMX_SIMD_REAL_WIDTH>(loadUWrapper, store, i, 0);
    }
}
#        endif // GMX_SIMD_HAVE_LOADU

#        if GMX_SIMD_HAVE_STOREU
TEST(SimdBootstrapTest, storeUI)
{
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        loadStoreTester<int, SimdInt32, GMX_SIMD_REAL_WIDTH>(loadWrapper, storeU, 0, i);
    }
}
#        endif // GMX_SIMD_HAVE_STOREU
#    endif     // GMX_SIMD_HAVE_REAL

#    if GMX_SIMD4_HAVE_REAL
TEST(SimdBootstrapTest, simd4LoadStore)
{
    loadStoreTester<real, Simd4Real, GMX_SIMD4_WIDTH>(load4, store4, 0, 0);
}

#        if GMX_SIMD_HAVE_LOADU
TEST(SimdBootstrapTest, simd4LoadU)
{
    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        loadStoreTester<real, Simd4Real, GMX_SIMD4_WIDTH>(load4U, store4, i, 0);
    }
}
#        endif // GMX_SIMD_HAVE_LOADU

#        if GMX_SIMD_HAVE_STOREU
TEST(SimdBootstrapTest, simd4StoreU)
{
    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        loadStoreTester<real, Simd4Real, GMX_SIMD4_WIDTH>(load4, store4U, 0, i);
    }
}
#        endif // GMX_SIMD_HAVE_STOREU
#    endif     // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

} // namespace

} // namespace test

} // namespace gmx

#endif // GMX_SIMD
