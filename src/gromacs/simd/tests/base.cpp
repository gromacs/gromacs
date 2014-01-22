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

#include "testutils/testoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/basicoptions.h"

#include "base.h"

namespace gmx
{
namespace test
{
    
int  SimdBaseTest::s_nPoints    = 10000;
    
/*! \brief Command-line option to adjust the number of points used to test SIMD math functions. */
GMX_TEST_OPTIONS(SimdBaseTestOptions, options)
{
    options->addOption(::gmx::IntegerOption("npoints")
                       .store(&SimdBaseTest::s_nPoints)
                       .description("Number of points to test for SIMD math functions"));
}
    
testing::AssertionResult
SimdBaseTest::compareVectorRealUlp(const char * refExpr,   const char * tstExpr,
                                   std::vector<real> ref, std::vector<real> tst)
{
    std::vector<gmx_int64_t>      ulpDiff(tst.size());
    bool                          eq     = true;
    bool                          signOk = true;
    size_t                        i;

    union {
#ifdef GMX_DOUBLE
        double r; gmx_int64_t i;
#else
        float r; gmx_int32_t i;
#endif
    } conv0, conv1;
    
    // Internal test of the test - make sure reference and test have the same length.
    if(ref.size() != tst.size())
    {
        return ::testing::AssertionFailure()
            << "Internal test error - unequal size vectors in compareVectorRealUlp" << std::endl;
    }
    
    for (i = 0; i < tst.size() && eq == true; i++)
    {
        eq     = eq && ( fabs(ref[i]-tst[i]) < absTol_ );
        signOk = signOk && ( ref[i]*tst[i] >= 0 );
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
               << "Ref values:  " << ::testing::PrintToString(ref) << std::endl
               << "SIMD values: " << ::testing::PrintToString(tst) << std::endl;
    }

    for (i = 0, eq = true; i < tst.size(); i++)
    {
        conv0.r    = ref[i];
        conv1.r    = tst[i];
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
               << "Ref  values: " << ::testing::PrintToString(ref) << std::endl
               << "SIMD values: " << ::testing::PrintToString(tst) << std::endl
               << "Ulp diff.:   " << ::testing::PrintToString(ulpDiff) << std::endl;
    }
}

}      // namespace
}      // namespace
