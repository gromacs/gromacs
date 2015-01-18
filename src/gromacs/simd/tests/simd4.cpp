/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#include "simd4.h"

#include "gromacs/simd/simd.h"

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#ifdef GMX_SIMD4_HAVE_REAL

const gmx_simd4_real_t rSimd4_1_2_3    = setSimd4RealFrom3R(1, 2, 3);
const gmx_simd4_real_t rSimd4_4_5_6    = setSimd4RealFrom3R(4, 5, 6);
const gmx_simd4_real_t rSimd4_7_8_9    = setSimd4RealFrom3R(7, 8, 9);
const gmx_simd4_real_t rSimd4_5_7_9    = setSimd4RealFrom3R(5, 7, 9);
const gmx_simd4_real_t rSimd4_m1_m2_m3 = setSimd4RealFrom3R(-1, -2, -3);
const gmx_simd4_real_t rSimd4_3_1_4    = setSimd4RealFrom3R(3, 1, 4);
const gmx_simd4_real_t rSimd4_m3_m1_m4 = setSimd4RealFrom3R(-3, -1, -4);
const gmx_simd4_real_t rSimd4_2p25     = setSimd4RealFrom1R(2.25);
const gmx_simd4_real_t rSimd4_3p75     = setSimd4RealFrom1R(3.75);
const gmx_simd4_real_t rSimd4_m2p25    = setSimd4RealFrom1R(-2.25);
const gmx_simd4_real_t rSimd4_m3p75    = setSimd4RealFrom1R(-3.75);
const gmx_simd4_real_t rSimd4_Exp      = setSimd4RealFrom3R( 1.4055235171027452623914516e+18,
                                                             5.3057102734253445623914516e-13,
                                                             -2.1057102745623934534514516e+16);
#    if (defined GMX_SIMD_HAVE_DOUBLE) && (defined GMX_DOUBLE)
// Make sure we also test exponents outside single precision when we use double
const gmx_simd4_real_t  rSimd_ExpDouble = setSimd4RealFrom3R( 6.287393598732017379054414e+176,
                                                              8.794495252903116023030553e-140,
                                                              -3.637060701570496477655022e+202);
#    endif

::std::vector<real>
simd4Real2Vector(const gmx_simd4_real_t simd4)
{
    real                mem[GMX_SIMD4_WIDTH*2];
    real *              p = gmx_simd4_align_r(mem);

    gmx_simd4_store_r(p, simd4);
    std::vector<real>   v(p, p+GMX_SIMD4_WIDTH);

    return v;
}

gmx_simd4_real_t
vector2Simd4Real(const std::vector<real> &v)
{
    real                mem[GMX_SIMD4_WIDTH*2];
    real *              p = gmx_simd4_align_r(mem);

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        p[i] = v[i % v.size()];  // repeat vector contents to fill simd width
    }
    return gmx_simd4_load_r(p);
}

gmx_simd4_real_t
setSimd4RealFrom3R(real r0, real r1, real r2)
{
    std::vector<real> v(3);
    v[0] = r0;
    v[1] = r1;
    v[2] = r2;
    return vector2Simd4Real(v);
}

gmx_simd4_real_t
setSimd4RealFrom1R(real value)
{
    std::vector<real> v(GMX_SIMD4_WIDTH);
    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        v[i] = value;
    }
    return vector2Simd4Real(v);
}

testing::AssertionResult
Simd4Test::compareSimd4RealUlp(const char *  refExpr,     const char *  tstExpr,
                               const gmx_simd4_real_t ref, const gmx_simd4_real_t tst)
{
    return compareVectorRealUlp(refExpr, tstExpr, simd4Real2Vector(ref), simd4Real2Vector(tst));
}

testing::AssertionResult
Simd4Test::compareSimd4RealEq(const char * refExpr, const char * tstExpr,
                              const gmx_simd4_real_t ref, const gmx_simd4_real_t tst)
{
    return compareVectorEq(refExpr, tstExpr, simd4Real2Vector(ref), simd4Real2Vector(tst));
}

#endif  // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
