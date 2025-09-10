/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*! \internal \file
 * \brief Testing gmx::accessor_policy.
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "gmxpre.h"

#include "gromacs/mdspan/accessor_policy.h"

#include <cstddef>

#include <array>
#include <string>

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{
namespace
{

class BasicAccessorPolicy : public ::testing::Test
{
public:
    std::array<float, 3>  testdata = { { 1., 2., 3. } };
    accessor_basic<float> acc;
};

TEST_F(BasicAccessorPolicy, Decay)
{
    EXPECT_EQ(acc.decay(testdata.data()), testdata.data());
}

TEST_F(BasicAccessorPolicy, Access)
{
    for (size_t i = 0; i < testdata.size(); ++i)
    {
        EXPECT_EQ(acc.access(testdata.data(), i), testdata[i]);
    }
}

TEST_F(BasicAccessorPolicy, Offset)
{
    for (size_t i = 0; i < testdata.size(); ++i)
    {
        EXPECT_EQ(acc.offset(testdata.data(), i), testdata.data() + i);
    }
}

TEST_F(BasicAccessorPolicy, CopyAccessor)
{
    const auto newAcc = acc;

    EXPECT_EQ(acc.decay(testdata.data()), newAcc.decay(testdata.data()));
    for (size_t i = 0; i < testdata.size(); ++i)
    {
        EXPECT_EQ(acc.access(testdata.data(), i), newAcc.access(testdata.data(), i));
    }

    for (size_t i = 0; i < testdata.size(); ++i)
    {
        EXPECT_EQ(acc.offset(testdata.data(), i), newAcc.offset(testdata.data(), i));
    }
}

} // namespace
} // namespace test
} // namespace gmx
