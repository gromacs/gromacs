/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Tests matrix class
 *
 * \author Kevin Boyd <kevin.boyd@uconn.edu>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include <array>

#include <gtest/gtest.h>

#include "gromacs/math/matrix.h"
#include "gromacs/math/vec.h"
namespace gmx
{


TEST(Matrix3x3, CanConstruct) {
    Matrix3x3<int> x;
    x(1, 1) = 1;
}


TEST(Matrix3x3Test, CanBeStoredInVector)
{
    std::vector < Matrix3x3 < int>> v;
    v.emplace_back(1, 2, 3, 4, 5, 6, 7, 8, 9);
    v.resize(2);
    EXPECT_EQ(1, v[0][XX][XX]);
    EXPECT_EQ(2, v[0][XX][YY]);
    EXPECT_EQ(3, v[0][XX][ZZ]);
    EXPECT_EQ(4, v[0][YY][XX]);
    EXPECT_EQ(5, v[0][YY][YY]);
    EXPECT_EQ(6, v[0][YY][ZZ]);
    EXPECT_EQ(7, v[0][ZZ][XX]);
    EXPECT_EQ(8, v[0][ZZ][YY]);
    EXPECT_EQ(9, v[0][ZZ][ZZ]);
}


} // namespace gmx
