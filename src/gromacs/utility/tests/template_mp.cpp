/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#include "gromacs/utility/template_mp.h"

#include <string>

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{
namespace
{

enum class Options
{
    Op0   = 0,
    Op1   = 1,
    Op2   = 2,
    Count = 3
};

template<Options i, Options j>
static int testEnumTwoIPlusJPlusK(int k)
{
    return 2 * int(i) + int(j) + k;
}

template<bool doDoubling, Options i, Options j>
static int testBoolEnumTwoIPlusJPlusK(int k)
{
    return (doDoubling ? 2 : 1) * int(i) + int(j) + k;
}

template<bool doDoubling>
static int testBoolDoubleOrNot(int k)
{
    return (doDoubling ? 2 : 1) * k;
}


TEST(TemplateMPTest, DispatchTemplatedFunctionEnum)
{
    int five           = 5;
    int two1plus2plus5 = dispatchTemplatedFunction(
            [=](auto p1, auto p2) { return testEnumTwoIPlusJPlusK<p1, p2>(five); }, Options::Op1, Options::Op2);
    EXPECT_EQ(two1plus2plus5, 9);
}

TEST(TemplateMPTest, DispatchTemplatedFunctionBool)
{
    int five = 5;
    int double5 = dispatchTemplatedFunction([=](auto p1) { return testBoolDoubleOrNot<p1>(five); }, true);
    EXPECT_EQ(double5, 10);
    int just5 = dispatchTemplatedFunction([=](auto p1) { return testBoolDoubleOrNot<p1>(five); }, false);
    EXPECT_EQ(just5, 5);
}

TEST(TemplateMPTest, DispatchTemplatedFunctionEnumBool)
{
    int five           = 5;
    int two1plus2plus5 = dispatchTemplatedFunction(
            [=](auto p1, auto p2, auto p3) { return testBoolEnumTwoIPlusJPlusK<p1, p2, p3>(five); },
            true,
            Options::Op1,
            Options::Op2);
    EXPECT_EQ(two1plus2plus5, 9);
}

} // anonymous namespace
} // namespace test
} // namespace gmx
