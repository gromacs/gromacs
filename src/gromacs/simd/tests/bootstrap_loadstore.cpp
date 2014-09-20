/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

/*! \internal \file
 * \brief
 * Separate test of SIMD load/store, before we use them in the SIMD test classes.
 *
 * Simple tests without using any classes/utilities, so we can use load/store
 * functions inside our test utilities after this has passed.
 *
 * This file tests:
 *
 * - gmx_simd_align_r(),gmx_simd_align_i(),gmx_simd4_align_r(),
 * - gmx_simd_load_r(),gmx_simd_store_r(),gmx_simd_loadu_r(),gmx_simd_storeu_r()
 * - gmx_simd_load_i(),gmx_simd_store_i(), gmx_simd_loadu_i(),gmx_simd_storeu_i()
 * - gmx_simd4_load_r(),gmx_simd4_store_r(), gmx_simd4_loadu_r(),gmx_simd4_storeu_r()
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 * \ingroup module_simd
 */

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/real.h"

namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

TEST(SimdBootstrapTest, gmxSimdAlign)
{
#ifdef GMX_SIMD_HAVE_REAL
    real rdata[GMX_SIMD_REAL_WIDTH*2];
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_EQ(((size_t)gmx_simd_align_r(&rdata[i]) & (GMX_SIMD_REAL_WIDTH*sizeof(real)-1)), (size_t)0);
    }
#endif
#ifdef GMX_SIMD_HAVE_INT32
    int idata[GMX_SIMD_INT32_WIDTH*2];
    for (int i = 0; i < GMX_SIMD_INT32_WIDTH; i++)
    {
        EXPECT_EQ(((size_t)gmx_simd_align_i(&idata[i]) & (GMX_SIMD_INT32_WIDTH*sizeof(int)-1)), (size_t)0);
    }
#endif
}

/*! \brief Generic routine to test load & store of SIMD, and check for side effects.
 *
 * The tests for load, store, unaligned load and unaligned store both for
 * real and int are pretty much similar, so we use a template function with
 * additional function pointers for the actual load/store calls. This would
 * be more hacking to turn into a class, since the SIMD functionality uses
 * macros rather than functions that can be overloaded.
 */
template <typename T, typename TSimd> void
simdLoadStoreTester(TSimd simdLoadFn(T* mem), void simdStoreFn(T* mem, TSimd),
                    T * simdAlignFn(T *mem),
                    const int loadOffset, const int storeOffset, const int simdWidth)
{
    /* We want simdWidth elements before the data to check we are not polluting
     * memory. Then we need 2*simdWidth storage to be able to extract an aligned
     * pointer, another simdWidth elements so we can create (deliberately)
     * offset un-aligned pointers, and finally simdWidth elements at the end
     * to test we are not polluting memory there either. Sum=5*simdWidth!
     */
    std::vector<T>   src(simdWidth*5);
    std::vector<T>   dst(simdWidth*5);
    // Make sure we have memory to check both before and after the test pointers
    T *              pCopySrc = simdAlignFn(&src[0]) + simdWidth + loadOffset;
    T *              pCopyDst = simdAlignFn(&dst[0]) + simdWidth + storeOffset;
    int              i;

    for (i = 0; i < simdWidth*5; i++)
    {
        src[i] =  1+i;
        dst[i] = -1-i;
    }

    simdStoreFn(pCopyDst, simdLoadFn(pCopySrc));

    for (i = 0; i < simdWidth; i++)
    {
        EXPECT_EQ(pCopySrc[i], pCopyDst[i]) << "SIMD load or store not moving data correctly for element " << i;
    }

    for (i = 0; i < simdWidth*5; i++)
    {
        EXPECT_EQ(src[i], (T)(1+i)) << "Side effect on source memory, i = " << i;
        if (&dst[0]+i < pCopyDst || &dst[0]+i >= pCopyDst+simdWidth)
        {
            EXPECT_EQ(dst[i], (T)(-1-i)) << "Side effect on destination memory, i = " << i;
        }
    }
}

#ifdef GMX_SIMD_HAVE_REAL
//! Wrapper for SIMD macro to load aligned floating-point data.
gmx_simd_real_t wrapperSimdLoadR(real *m)
{
    return gmx_simd_load_r(m);
}
//! Wrapper for SIMD macro to store to aligned floating-point data.
void            wrapperSimdStoreR(real *m, gmx_simd_real_t s)
{
    gmx_simd_store_r(m, s);
}

TEST(SimdBootstrapTest, gmxSimdLoadStoreR)
{
    simdLoadStoreTester(wrapperSimdLoadR, wrapperSimdStoreR, gmx_simd_align_r, 0, 0, GMX_SIMD_REAL_WIDTH);
}

#    ifdef GMX_SIMD_HAVE_LOADU
//! Wrapper for SIMD macro to load unaligned floating-point data.
gmx_simd_real_t WrapperSimdLoadUR(real *m)
{
    return gmx_simd_loadu_r(m);
}

TEST(SimdBootstrapTest, gmxSimdLoadUR)
{
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        simdLoadStoreTester(WrapperSimdLoadUR, wrapperSimdStoreR, gmx_simd_align_r, i, 0, GMX_SIMD_REAL_WIDTH);
    }
}
#    endif

#    ifdef GMX_SIMD_HAVE_STOREU
//! Wrapper for SIMD macro to store to unaligned floating-point data.
void WrapperSimdStoreUR(real *m, gmx_simd_real_t s)
{
    gmx_simd_storeu_r(m, s);
}

TEST(SimdBootstrapTest, gmxSimdStoreUR)
{
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        simdLoadStoreTester(wrapperSimdLoadR, WrapperSimdStoreUR, gmx_simd_align_r, 0, i, GMX_SIMD_REAL_WIDTH);
    }
}
#    endif
#endif

#ifdef GMX_SIMD_HAVE_INT32
// Tests for gmx_simd_int32_t load & store operations

//! Wrapper for SIMD macro to load aligned integer data.
gmx_simd_int32_t wrapperSimdLoadI(int *m)
{
    return gmx_simd_load_i(m);
}
//! Wrapper for SIMD macro to store to aligned integer data.
void             wrapperSimdStoreI(int *m, gmx_simd_int32_t s)
{
    gmx_simd_store_i(m, s);
}

TEST(SimdBootstrapTest, gmxSimdLoadStoreI)
{
    simdLoadStoreTester(wrapperSimdLoadI, wrapperSimdStoreI, gmx_simd_align_i, 0, 0, GMX_SIMD_INT32_WIDTH);
}

#    ifdef GMX_SIMD_HAVE_LOADU
//! Wrapper for SIMD macro to load unaligned integer data.
gmx_simd_int32_t wrapperSimdLoadUI(int *m)
{
    return gmx_simd_loadu_i(m);
}

TEST(SimdBootstrapTest, gmxSimdLoadUI)
{
    for (int i = 0; i < GMX_SIMD_INT32_WIDTH; i++)
    {
        simdLoadStoreTester(wrapperSimdLoadUI, wrapperSimdStoreI, gmx_simd_align_i, i, 0, GMX_SIMD_INT32_WIDTH);
    }
}
#    endif

#    ifdef GMX_SIMD_HAVE_STOREU
//! Wrapper for SIMD macro to store to unaligned integer data.
void wrapperSimdStoreUI(int *m, gmx_simd_int32_t s)
{
    gmx_simd_storeu_i(m, s);
}

TEST(SimdBootstrapTest, gmxSimdStoreUI)
{
    for (int i = 0; i < GMX_SIMD_INT32_WIDTH; i++)
    {
        simdLoadStoreTester(wrapperSimdLoadI, wrapperSimdStoreUI, gmx_simd_align_i, 0, i, GMX_SIMD_INT32_WIDTH);
    }
}
#    endif
#endif

#ifdef GMX_SIMD4_HAVE_REAL
/* Tests for gmx_simd4_real_t load & store operations. Define wrapper functions
 * for the SIMD instructions that are typically implemented as macros.
 */

/*! \brief Separate load/store tester function for SIMD4.
 *
 * Due to the way SIMD variables
 * are implemented as deep internal data, some compilers treat them as
 * float/double with special prefixes. Unfortunately, this means that some C++
 * compilers think an 8-wide normal real SIMD and a 4-wide SIMD4 real type
 * cannot be overloaded (e.g. with gcc using 256-bit AVX single precision).
 */
template <typename T, typename TSimd> void
simd4LoadStoreTester(TSimd simd4LoadFn(T* mem), void simd4StoreFn(T* mem, TSimd),
                     T * simd4AlignFn(T *mem),
                     const int loadOffset, const int storeOffset)
{
    /* We want simdWidth elements before the data to check we are not polluting
     * memory. Then we need 2*simdWidth storage to be able to extract an aligned
     * pointer, another simdWidth elements so we can create (deliberately)
     * offset un-aligned pointers, and finally simdWidth elements at the end
     * to test we are not polluting memory there either. Sum=5*simdWidth!
     */
    T         src[GMX_SIMD4_WIDTH*5];
    T         dst[GMX_SIMD4_WIDTH*5];
    // Make sure we have memory to check both before and after the test pointers
    T *       pCopySrc = simd4AlignFn(src) + GMX_SIMD4_WIDTH + loadOffset;
    T *       pCopyDst = simd4AlignFn(dst) + GMX_SIMD4_WIDTH + storeOffset;
    int       i;

    for (i = 0; i < GMX_SIMD4_WIDTH*5; i++)
    {
        src[i] =  1+i;
        dst[i] = -1-i;
    }

    simd4StoreFn(pCopyDst, simd4LoadFn(pCopySrc));

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        EXPECT_EQ(pCopySrc[i], pCopyDst[i]) << "SIMD4 load or store not moving data correctly for element " << i;
    }

    for (i = 0; i < GMX_SIMD4_WIDTH*5; i++)
    {
        EXPECT_EQ(src[i], (T)(1+i)) << "Side effect on source memory, i = " << i;
        if (dst+i < pCopyDst || dst+i >= pCopyDst+GMX_SIMD4_WIDTH)
        {
            EXPECT_EQ(dst[i], (T)(-1-i)) << "Side effect on destination memory, i = " << i;
        }
    }
}

//! Wrapper for SIMD4 macro to load aligned floating-point data.
gmx_simd4_real_t wrapperSimd4LoadR(real *m)
{
    return gmx_simd4_load_r(m);
}
//! Wrapper for SIMD4 macro to store to aligned floating-point data.
void             wrapperSimd4StoreR(real *m, gmx_simd4_real_t s)
{
    gmx_simd4_store_r(m, s);
}

TEST(SimdBootstrapTest, gmxSimd4LoadStoreR)
{
    simd4LoadStoreTester(wrapperSimd4LoadR, wrapperSimd4StoreR, gmx_simd4_align_r, 0, 0);
}

#    ifdef GMX_SIMD_HAVE_LOADU
//! Wrapper for SIMD4 macro to load unaligned floating-point data.
gmx_simd4_real_t WrapperSimd4LoadUR(real *m)
{
    return gmx_simd4_loadu_r(m);
}

TEST(SimdBootstrapTest, gmxSimd4LoadUR)
{
    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        simd4LoadStoreTester(WrapperSimd4LoadUR, wrapperSimd4StoreR, gmx_simd4_align_r, i, 0);
    }
}
#    endif

#    ifdef GMX_SIMD_HAVE_STOREU
//! Wrapper for SIMD4 macro to store to unaligned floating-point data.
void WrapperSimd4StoreUR(real *m, gmx_simd4_real_t s)
{
    gmx_simd4_storeu_r(m, s);
}

TEST(SimdBootstrapTest, gmxSimd4StoreUR)
{
    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        simd4LoadStoreTester(wrapperSimd4LoadR, WrapperSimd4StoreUR, gmx_simd4_align_r, 0, i);
    }
}
#    endif
#endif

/*! \} */
/*! \endcond */

} // namespace
