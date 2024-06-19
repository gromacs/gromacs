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

#include "base.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>

#include <string>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/basedefinitions.h"

#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{

namespace
{

/*! \cond */
/*! \brief Command-line option to adjust the number of points used to test SIMD math functions. */
GMX_TEST_OPTIONS(SimdBaseTestOptions, options)
{
    options->addOption(::gmx::IntegerOption("npoints")
                               .store(&SimdBaseTest::s_nPoints)
                               .description("Number of points to test for SIMD math functions"));
}
/*! \endcond */

} // namespace

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

int SimdBaseTest::s_nPoints = 10000;

::testing::AssertionResult SimdBaseTest::compareVectorRealUlp(const char*              refExpr,
                                                              const char*              tstExpr,
                                                              const std::vector<real>& ref,
                                                              const std::vector<real>& tst) const
{
    std::vector<real>         absDiff(tst.size());
    std::vector<std::int64_t> ulpDiff(tst.size());
    bool                      allOk;
    size_t                    i;

    union
    {
#if GMX_DOUBLE
        double       r;
        std::int64_t i;
#else
        float        r;
        std::int32_t i;
#endif
    } conv0, conv1;

    // Internal test of the test - make sure reference and test have the same length.
    if (ref.size() != tst.size())
    {
        return ::testing::AssertionFailure()
               << "Internal test error - unequal size vectors in compareVectorRealUlp" << std::endl;
    }

    for (i = 0, allOk = true; i < tst.size(); i++)
    {
        absDiff[i] = std::abs(ref[i] - tst[i]);
        conv0.r    = ref[i];
        conv1.r    = tst[i];
        ulpDiff[i] = std::llabs(conv0.i - conv1.i);

        /* Use strict smaller-than for absolute tolerance check, so we disable it with absTol_=0 */
        allOk = allOk && ((absDiff[i] < absTol_) || ((ref[i] * tst[i] >= 0) && (ulpDiff[i] <= ulpTol_)));
    }

    if (allOk)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing comparison between " << refExpr << " and " << tstExpr << std::endl
               << "Requested abs tolerance: " << absTol_ << std::endl
               << "Requested ulp tolerance: " << ulpTol_ << std::endl
               << "(And values should not differ in sign unless within abs tolerance.)" << std::endl
               << "Reference values: " << ::testing::PrintToString(ref) << std::endl
               << "SIMD values:      " << ::testing::PrintToString(tst) << std::endl
               << "Abs. difference:  " << ::testing::PrintToString(absDiff) << std::endl
               << "Ulp difference:   " << ::testing::PrintToString(ulpDiff) << std::endl;
    }
}

/*! \} */
/*! \endcond */

} // namespace test
} // namespace gmx
