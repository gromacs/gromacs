/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * This implements basic nblib box tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gromacs/nblib/box.h"

#include <cmath>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

using gmx::test::defaultRealTolerance;

namespace nblib
{

TEST(NBlibTest, CubicBoxCannotHaveNaN)
{
    real number = NAN;
    EXPECT_THROW(Box box(number), gmx::InvalidInputError);
}

TEST(NBlibTest, CubicBoxCannotHaveInf)
{
    real number = INFINITY;
    EXPECT_THROW(Box box(number), gmx::InvalidInputError);
}

TEST(NBlibTest, RectangularBoxCannotHaveNaN)
{
    real number = NAN;
    EXPECT_THROW(Box box(number, real(1.), real(1.)), gmx::InvalidInputError);
}

TEST(NBlibTest, RectangularBoxCannotHaveInf)
{
    real number = INFINITY;
    EXPECT_THROW(Box box(number, real(1.), real(1.)), gmx::InvalidInputError);
}

TEST(NBlibTest, CubicBoxWorks)
{
    real length = 3;
    Box::Matrix ref = {{length, 0, 0, 0, length, 0, 0, 0, length}};
    Box::Matrix probe = Box(length).matrix();

    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            EXPECT_REAL_EQ_TOL(ref(i, j), probe(i, j), defaultRealTolerance());
}

}  // namespace nblib
