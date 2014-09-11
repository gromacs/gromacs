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

#include <math.h>

#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"

#include "simd4.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#ifdef GMX_SIMD4_HAVE_REAL

/*! \brief Test fixture for SIMD4 vector operations (identical to the SIMD4 \ref Simd4Test) */
typedef Simd4Test Simd4VectorOperationsTest;

TEST_F(Simd4VectorOperationsTest, gmxSimd4CalcRsqR)
{
    gmx_simd4_real_t simdX  = setSimd4RealFrom3R(1, 2, 3);
    gmx_simd4_real_t simdY  = setSimd4RealFrom3R(3, 0, 5);
    gmx_simd4_real_t simdZ  = setSimd4RealFrom3R(4, 1, 8);
    gmx_simd4_real_t simdR2 = setSimd4RealFrom3R(26, 5, 98);

    setUlpTol(2);
    GMX_EXPECT_SIMD4_REAL_NEAR(simdR2, gmx_simd4_calc_rsq_r(simdX, simdY, simdZ));
}

#endif      // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace
