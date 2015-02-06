/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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

#include <vector>

#include "gromacs/math/utilities.h"

#include "simd.h"
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

#if defined GMX_SIMD_HAVE_REAL && defined GMX_SIMD4_HAVE_REAL

/*! \brief Test fixture for SIMD4 tranposition operations */
class Simd4TransposeTest : public SimdTest
{
    public:
        //! Constructor
        Simd4TransposeTest()
        {
            /* Make a matrix to transpose, and the expected result of the
               transposition. */
            for (int i = 0; i < GMX_SIMD_REAL_WIDTH; ++i)
            {
                for (int j = 0; j < GMX_SIMD_REAL_WIDTH; ++j)
                {
                    real value = i * GMX_SIMD_REAL_WIDTH + j;
                    inputValues_[i].push_back(value);
                    expectedResultValues_[j].push_back(value);
                }
            }
        }

        //! Input values
        std::vector<real> inputValues_[GMX_SIMD_REAL_WIDTH];
        //! Expected transposition result values
        std::vector<real> expectedResultValues_[GMX_SIMD_REAL_WIDTH];
};

#ifdef GMX_SIMD4_HAVE_SIMD_TRANSPOSE

/* Provide some precision-agnostic names, because we are unit testing
   functions internal to the SIMD module that aren't currently
   intended for general use. */
#ifdef GMX_DOUBLE
#  define gmx_simd_transpose4_r gmx_simd_transpose4_d
#  define gmx_simd4_transpose_to_simd_r gmx_simd4_transpose_to_simd_d
#  define gmx_simd_transpose_to_simd4_r gmx_simd_transpose_to_simd4_d
#else
#  define gmx_simd_transpose4_r gmx_simd_transpose4_f
#  define gmx_simd4_transpose_to_simd_r gmx_simd4_transpose_to_simd_f
#  define gmx_simd_transpose_to_simd4_r gmx_simd_transpose_to_simd4_f
#endif

TEST_F(Simd4TransposeTest, gmxSimdTranspose4R)
{
    gmx_simd_real_t rSimd[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t rSimdExpected[GMX_SIMD_REAL_WIDTH];

    /* Get the matrix into aligned memory and load into SIMD
       registers */
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; ++i)
    {
        rSimd[i]         = vector2SimdReal(inputValues_[i]);
        rSimdExpected[i] = vector2SimdReal(expectedResultValues_[i]);
    }

    gmx_simd_transpose4_r(&rSimd[0], &rSimd[1], &rSimd[2], &rSimd[3]);

    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[0], rSimd[0], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[1], rSimd[1], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[2], rSimd[2], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[3], rSimd[3], 4);
}

TEST_F(Simd4TransposeTest, gmxSimd4TransposeToSimdR)
{
    gmx_simd4_real_t rSimd4[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  rSimd[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  rSimdExpected[GMX_SIMD_REAL_WIDTH];

    /* Get the matrix into aligned memory and load into SIMD
       registers */
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; ++i)
    {
        rSimd4[i]        = vector2Simd4Real(inputValues_[i]);
        rSimdExpected[i] = vector2SimdReal(expectedResultValues_[i]);
    }

    gmx_simd4_transpose_to_simd_r(rSimd4, &rSimd[0], &rSimd[1], &rSimd[2], &rSimd[3]);

    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[0], rSimd[0], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[1], rSimd[1], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[2], rSimd[2], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[3], rSimd[3], 4);
}

TEST_F(Simd4TransposeTest, gmxSimdTransposeToSimd4R)
{
    gmx_simd4_real_t rSimd4[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  rSimd[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  rSimdExpected[GMX_SIMD_REAL_WIDTH];

    /* Get the matrix into aligned memory and load into SIMD
       registers */
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; ++i)
    {
        rSimd[i]         = vector2SimdReal(inputValues_[i]);
        rSimdExpected[i] = vector2SimdReal(expectedResultValues_[i]);
    }

    gmx_simd_transpose_to_simd4_r(rSimd[0], rSimd[1], rSimd[2], rSimd[3], rSimd4);

    /* Get the results back into SIMD registers */
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; ++i)
    {
        rSimd[i] = vector2SimdReal(simd4Real2Vector(rSimd4[i]));
    }

    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[0], rSimd[0], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[1], rSimd[1], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[2], rSimd[2], 4);
    GMX_EXPECT_SIMD_REALN_EQ(rSimdExpected[3], rSimd[3], 4);
}

#endif      // GMX_SIMD4_HAVE_SIMD_TRANSPOSE

#endif      // GMX_SIMD_HAVE_REAL && GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace
