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

#include <cmath>
#include <cstdint>

#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/real.h"

#include "simd4.h"

#if GMX_SIMD

namespace gmx
{
namespace test
{


#    if GMX_SIMD4_HAVE_REAL

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

class Simd4MathTest : public Simd4Test
{
};

/*! \} */
/*! \endcond */

// Actual math function tests below

namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

// Presently, the only SIMD4 math function is 1/sqrt(x), which
// has a close-to-trivial implementation without different
// approximation intervals or special threshold points. To
// avoid having to re-implement the entire SIMD math function
// test infrastructure we only test these functions for a few
// values that are either special or exercise all bits.

TEST_F(Simd4MathTest, invsqrt)
{
    const real x0 = std::numeric_limits<float>::min();
    const real x1 = std::numeric_limits<float>::max();
    const real x2 = M_PI;

    GMX_EXPECT_SIMD4_REAL_NEAR(
            setSimd4RealFrom3R(1.0 / std::sqrt(x0), 1.0 / std::sqrt(x1), 1.0 / std::sqrt(x2)),
            invsqrt(setSimd4RealFrom3R(x0, x1, x2)));
}

TEST_F(Simd4MathTest, invsqrtSingleAccuracy)
{
    const real x0 = std::numeric_limits<float>::min();
    const real x1 = std::numeric_limits<float>::max();
    const real x2 = M_PI;

    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    GMX_EXPECT_SIMD4_REAL_NEAR(
            setSimd4RealFrom3R(1.0 / std::sqrt(x0), 1.0 / std::sqrt(x1), 1.0 / std::sqrt(x2)),
            invsqrtSingleAccuracy(setSimd4RealFrom3R(x0, x1, x2)));
}

/*! \} */
/*! \endcond */

} // namespace

#    endif // GMX_SIMD4_HAVE_REAL

} // namespace test
} // namespace gmx

#endif // GMX_SIMD
