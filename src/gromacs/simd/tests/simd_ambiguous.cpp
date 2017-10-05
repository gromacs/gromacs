/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"

/* Test that math functions which can be used both with scalar and SIMD
 * are ambiguous when applied to value returned from load.
 *
 * gmx::load returns a proxy/reference object which can be casted to either
 * a scalar (e.g. float) or a SIMD value (e.g. SIMDFloat). The gmx math
 * functions (e.g. sqrt) take both a scalar and a SIMD value as an argument.
 * Thus e.g. load(sqrt(m)) should be ambiguous. This test makes sure that
 * this does not compile. This got previously broken by introducing templates
 * which influenced the overload resolution.
 *
 * The test execution code in CMakeLists.txt tests that the code doesn't
 * compile with a SIMD implementation. To test that this code does correctly
 * compile besides causing the ambiguous overload error, it expects to
 * correctly compile for a a non-simd build. For such a build the
 * code is non-ambiguous because only the scalar version exists.
 *
 * The test execution code passes either float/double as TEST_PREC and the math
 * function to test as TEST_FUNC. Both are passed as compile definitions.
 * The file is compiled once for each combination when executing ctest and
 * the test fails if the file compiles.
 *
 * Possible extensions: Test all other math functions including those taking
 * multiple arguments.
 */
int main()
{
    /* We cannot check for SIMD float or double support at cmake,
     * only for general SIMD support. Therefore with SIMD, but without
     * SIMD float or double support we make the compilation fail with
     * a static_assert instead of the ambiguous overload error
     */
    constexpr bool testFloat       = std::is_same<TEST_PREC, float>::value;
    constexpr bool testDouble      = std::is_same<TEST_PREC, double>::value;
    constexpr bool haveSimdSupport = GMX_SIMD;
    constexpr bool haveSimdFloat   = GMX_SIMD_HAVE_FLOAT;
    constexpr bool haveSimdDouble  = GMX_SIMD_HAVE_DOUBLE;
    static_assert(!haveSimdSupport || !testFloat || haveSimdFloat, "Assertion failure to make test fail without SIMD float support");
    static_assert(!haveSimdSupport || !testDouble || haveSimdDouble, "Assertion failure to make test fail without SIMD double support");

    TEST_PREC  d = 0;
    TEST_PREC *m = &d;
    gmx::TEST_FUNC(gmx::load(m));
}
