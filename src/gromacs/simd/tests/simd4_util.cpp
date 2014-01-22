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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "simd4_util.h"

namespace gmx
{
namespace test
{

#ifdef GMX_SIMD4_HAVE_REAL
const gmx_simd4_real_t rSimd4_1_2_3    = setSimd4From3R(1, 2, 3);
const gmx_simd4_real_t rSimd4_4_5_6    = setSimd4From3R(4, 5, 6);
const gmx_simd4_real_t rSimd4_7_8_9    = setSimd4From3R(7, 8, 9);
const gmx_simd4_real_t rSimd4_5_7_9    = setSimd4From3R(5, 7, 9);
const gmx_simd4_real_t rSimd4_m1_m2_m3 = setSimd4From3R(-1, -2, -3);
const gmx_simd4_real_t rSimd4_3_1_4    = setSimd4From3R(3, 1, 4);
const gmx_simd4_real_t rSimd4_m3_m1_m4 = setSimd4From3R(-3, -1, -4);
const gmx_simd4_real_t rSimd4_2p25     = gmx_simd4_set1_r(2.25);
const gmx_simd4_real_t rSimd4_3p75     = gmx_simd4_set1_r(3.75);
const gmx_simd4_real_t rSimd4_m2p25    = gmx_simd4_set1_r(-2.25);
const gmx_simd4_real_t rSimd4_m3p75    = gmx_simd4_set1_r(-3.75);
const gmx_simd4_real_t rSimd4_Exp      = setSimd4From3R( 1.4055235171027452623914516e+18,
                                                         5.3057102734253445623914516e-13,
                                                         -2.1057102745623934534514516e+16);
#    ifdef GMX_DOUBLE
// Magic FP numbers corresponding to specific bit patterns
const gmx_simd4_real_t rSimd4_Bits1    = gmx_simd4_set1_r(-1.07730874267432137e+236);
const gmx_simd4_real_t rSimd4_Bits2    = gmx_simd4_set1_r(-9.25596313493178307e+061);
const gmx_simd4_real_t rSimd4_Bits3    = gmx_simd4_set1_r(-8.57750588235293981e+003);
const gmx_simd4_real_t rSimd4_Bits4    = gmx_simd4_set1_r( 1.22416778341839096e-250);
const gmx_simd4_real_t rSimd4_Bits5    = gmx_simd4_set1_r(-1.15711777004554095e+294);
const gmx_simd4_real_t rSimd4_Bits6    = gmx_simd4_set1_r( 1.53063836115600621e-018);
#    else
const gmx_simd4_real_t rSimd4_Bits1    = gmx_simd4_set1_r(-5.9654142337e+29);
const gmx_simd4_real_t rSimd4_Bits2    = gmx_simd4_set1_r(-1.0737417600e+08);
const gmx_simd4_real_t rSimd4_Bits3    = gmx_simd4_set1_r(-6.0235290527e+00);
const gmx_simd4_real_t rSimd4_Bits4    = gmx_simd4_set1_r( 1.0788832913e-31);
const gmx_simd4_real_t rSimd4_Bits5    = gmx_simd4_set1_r(-1.0508719529e+37);
const gmx_simd4_real_t rSimd4_Bits6    = gmx_simd4_set1_r( 1.1488970369e-02);
#    endif
#endif  // GMX_SIMD4_HAVE_REAL

#ifdef GMX_SIMD4_HAVE_REAL
::std::vector<real>
simd4Real2Vector(const gmx_simd4_real_t simd)
{
    real                mem[GMX_SIMD4_WIDTH*2];
    real *              p = gmx_simd4_align_r(mem);

    gmx_simd4_store_r(p, simd);
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
setSimd4From3R(real r0, real r1, real r2)
{
    std::vector<real> v(3);
    v[0] = r0;
    v[1] = r1;
    v[2] = r2;
    return vector2Simd4Real(v);
}

gmx_simd4_real_t
setSimd4From1R(real value)
{
    std::vector<real> v(GMX_SIMD4_WIDTH);
    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        v[i] = value;
    }
    return vector2Simd4Real(v);
}

std::string
printSimd4Real(const gmx_simd4_real_t simd)
{
    std::stringstream   buf;
    std::vector<real>   v = simd4Real2Vector(simd);
    buf << "[ ";
    buf << std::setprecision(20);
    copy(v.begin(), v.end(), std::ostream_iterator<real>(buf, " "));
    buf << "]";
    return buf.str();
}

testing::AssertionResult
Simd4Test::compareSimd4RealUlp(const char *  refExpr,    const char *  tstExpr,
                               const gmx_simd4_real_t ref, const gmx_simd4_real_t tst)
{
    std::vector<real>             vref = simd4Real2Vector(ref);
    std::vector<real>             vtst = simd4Real2Vector(tst);
    std::vector<gmx_int64_t>      ulpDiff(GMX_SIMD4_WIDTH);
    bool                          eq     = true;
    bool                          signOk = true;
    int                           i;
#    ifdef GMX_DOUBLE
    union {
        double r; gmx_int64_t i;
    } conv0, conv1;
#    else
    union {
        float  r; int i;
    } conv0, conv1;
#    endif

    for (i = 0; i < GMX_SIMD4_WIDTH && eq == true; i++)
    {
        eq     = eq && ( fabs(vref[i]-vtst[i]) < absTol_ );
        signOk = signOk && ( vref[i]*vtst[i] >= 0 );
    }
    if (eq == true)
    {
        // Pass test if all values were within absolute tolerance
        return ::testing::AssertionSuccess();
    }
    else if (signOk == false)
    {
        return ::testing::AssertionFailure()
               << "Failing comparison between " << refExpr << " and " << tstExpr << std::endl
               << "SIMD contents differs in sign from reference: " << std::endl
               << "Ref values:  " << printSimd4Real(ref) << std::endl
               << "SIMD values: " << printSimd4Real(tst) << std::endl;
    }

    for (i = 0, eq = true; i < GMX_SIMD4_WIDTH; i++)
    {
        conv0.r    = vref[i];
        conv1.r    = vtst[i];
        ulpDiff[i] = llabs(conv0.i-conv1.i);
        eq         = eq && (ulpDiff[i] <= ulpTol_);
    }

    if (eq == true)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing ulp comparison between " << refExpr << " and " << tstExpr << std::endl
               << "Requested ulp tolerance: " << ulpTol_ << std::endl
               << "Requested abs tolerance: " << absTol_ << std::endl
               << "Ref  values: " << printSimd4Real(ref) << std::endl
               << "SIMD values: " << printSimd4Real(tst) << std::endl
               << "Ulp diff.:   " << ::testing::PrintToString(ulpDiff) << std::endl;
    }
}

testing::AssertionResult
Simd4Test::compareSimd4RealUlp(const char * refExpr,    const char * tstExpr,
                               real ref, const gmx_simd4_real_t tst)
{
    return compareSimd4RealUlp(refExpr, tstExpr, setSimd4From1R(ref), tst);
}

testing::AssertionResult
Simd4Test::compareSimd4RealEq(const char * refExpr, const char * tstExpr,
                              const gmx_simd4_real_t ref, const gmx_simd4_real_t tst)
{
    // We need to alter tolerances, so save old values first
    gmx_int64_t oldUlpTol = ulpTol_;
    real        oldAbsTol = absTol_;
    // Set tolerances that enforce exact matching
    ulpTol_ = 0L;
    absTol_ = 0;
    return compareSimd4RealUlp(refExpr, tstExpr, ref, tst);
    // Reset previous tolerances
    ulpTol_ = oldUlpTol;
    absTol_ = oldAbsTol;
}

testing::AssertionResult
Simd4Test::compareSimd4RealEq(const char * refExpr, const char * tstExpr,
                              real ref,             const gmx_simd4_real_t tst)
{
    return compareSimd4RealEq(refExpr, tstExpr, setSimd4From1R(ref), tst);
}

#endif

} // namespace
} // namespace
