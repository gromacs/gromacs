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

#include "gromacs/compat/mp11.h"

#include <cstddef>

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"

// Defining some dummy functions to use later

namespace gmx
{
namespace compat
{
namespace test
{
enum class Options
{
    Op0   = 0,
    Op1   = 1,
    Op2   = 2,
    Count = 3
};

template<int i>
static int testIncrement()
{
    return i + 1;
}

template<bool i>
static bool testNot()
{
    return !i;
}

template<Options i, Options j>
static int testEnumTwoIPlusJPlusK(int k)
{
    return 2 * int(i) + int(j) + k;
}

TEST(TemplateMPTest, MpWithIndexInt)
{
    static constexpr int maxArgValue = 4;
    int inc_0 = mp_with_index<maxArgValue>(0, [](auto i) { return testIncrement<i>(); });
    EXPECT_EQ(inc_0, 1);
    int inc_3 = mp_with_index<maxArgValue>(3, [](auto i) { return testIncrement<i>(); });
    EXPECT_EQ(inc_3, 4);
}

TEST(TemplateMPTest, MpWithIndexIntBad)
{
    static constexpr int maxArgValue = 4;
    int                  i           = maxArgValue;
    // Function requirement: i < maxArgValue
    EXPECT_THROW(mp_with_index<maxArgValue>(i, [](auto i) { return testIncrement<i>(); }),
                 gmx::InternalError);
}

TEST(TemplateMPTest, MpWithIndexBool)
{
    bool not_true = mp_with_index<2>(size_t(true), [](auto i) { return testNot<i>(); });
    EXPECT_FALSE(not_true);
    bool not_false = mp_with_index<2>(size_t(false), [](auto i) { return testNot<i>(); });
    EXPECT_TRUE(not_false);
}

TEST(TemplateMPTest, MpWithIndexEnum)
{
    int five           = 5;
    int two1plus2plus5 = mp_with_index<static_cast<size_t>(Options::Count)>(
            static_cast<size_t>(Options::Op2), [=](auto i) {
                return testEnumTwoIPlusJPlusK<Options::Op1, static_cast<Options>(size_t(i))>(five);
            });
    EXPECT_EQ(two1plus2plus5, 9);
}

} // namespace test
} // namespace compat
} // namespace gmx
