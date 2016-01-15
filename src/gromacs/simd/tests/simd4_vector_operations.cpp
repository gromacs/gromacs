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

#include <cmath>

#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"

#include "simd4.h"

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

#if GMX_SIMD4_HAVE_REAL

/*! \brief Test fixture for SIMD4 vector operations (identical to the SIMD4 \ref Simd4Test) */
typedef Simd4Test Simd4VectorOperationsTest;

TEST_F(Simd4VectorOperationsTest, norm2)
{
    Simd4Real simdX  = setSimd4RealFrom3R(1, 2, 3);
    Simd4Real simdY  = setSimd4RealFrom3R(3, 0, 5);
    Simd4Real simdZ  = setSimd4RealFrom3R(4, 1, 8);
    Simd4Real simdR2 = setSimd4RealFrom3R(26, 5, 98);

    setUlpTol(2);
    GMX_EXPECT_SIMD4_REAL_NEAR(simdR2, norm2(simdX, simdY, simdZ));
}

#endif      // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace

#endif // GMX_SIMD
