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
/*
 * This file is a modified version of original work of Sandia Corporation.
 * In the spirit of the original code, this particular file can be distributed
 * on the terms of Sandia Corporation.
 */
/*
 *                         Kokkos v. 2.0
 *               Copyright (2014) Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Kokkos is licensed under 3-clause BSD terms of use:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Christian R. Trott (crtrott@sandia.gov)
 */
/*! \internal \file
 * \brief Testing gmx::extents.
 *
 * \author Christian Trott <crtrott@sandia.gov>
 * \author Carter Edwards <hedwards@nvidia.com>
 * \author David Hollman <dshollm@sandia.gov>
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "gmxpre.h"

#include "gromacs/mdspan/extents.h"

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

template<ptrdiff_t... E_STATIC>
class ExtentsTest
{
public:
    using extents_type = gmx::extents<E_STATIC...>;

    extents_type my_extents_explicit;
    extents_type my_extents_array;
    extents_type my_extents_copy;

    ExtentsTest()
    {
        my_extents_explicit = extents<E_STATIC...>();
        my_extents_array    = extents<E_STATIC...>(std::array<ptrdiff_t, 0>());
        my_extents_copy     = extents<E_STATIC...>(my_extents_explicit);
    }

    template<class... E>
    ExtentsTest(E... e)
    {
        my_extents_explicit = extents<E_STATIC...>(e...);
        my_extents_array    = extents<E_STATIC...>(std::array<ptrdiff_t, 2>({ { e... } }));
        my_extents_copy     = extents<E_STATIC...>(my_extents_explicit);
    }

    void check_rank(ptrdiff_t r)
    {
        EXPECT_EQ(my_extents_explicit.rank(), r);
        EXPECT_EQ(my_extents_array.rank(), r);
        EXPECT_EQ(my_extents_copy.rank(), r);
    }
    void check_rank_dynamic(ptrdiff_t r)
    {
        EXPECT_EQ(my_extents_explicit.rank_dynamic(), r);
        EXPECT_EQ(my_extents_array.rank_dynamic(), r);
        EXPECT_EQ(my_extents_copy.rank_dynamic(), r);
    }
    template<class... E>
    void check_extents(E... e)
    {
        std::array<ptrdiff_t, extents_type::rank()> s = { { E_STATIC... } };
        std::array<ptrdiff_t, extents_type::rank()> a = { { e... } };
        for (size_t r = 0; r < extents_type::rank(); r++)
        {
            EXPECT_EQ(my_extents_explicit.static_extent(r), s[r]);
            EXPECT_EQ(my_extents_explicit.extent(r), a[r]);

            EXPECT_EQ(my_extents_array.static_extent(r), s[r]);
            EXPECT_EQ(my_extents_array.extent(r), a[r]);

            EXPECT_EQ(my_extents_copy.static_extent(r), s[r]);
            EXPECT_EQ(my_extents_copy.extent(r), a[r]);
        }
        EXPECT_EQ(my_extents_explicit.static_extent(extents_type::rank() + 1), 1);
        EXPECT_EQ(my_extents_explicit.extent(extents_type::rank() + 1), 1);

        EXPECT_EQ(my_extents_array.static_extent(extents_type::rank() + 1), 1);
        EXPECT_EQ(my_extents_array.extent(extents_type::rank() + 1), 1);

        EXPECT_EQ(my_extents_copy.static_extent(extents_type::rank() + 1), 1);
        EXPECT_EQ(my_extents_copy.extent(extents_type::rank() + 1), 1);
    }
};

TEST(ExtentsTest, Construction)
{

    // setting two dynamic extents
    ExtentsTest<5, dynamic_extent, 3, dynamic_extent, 1> test(4, 2);

    test.check_rank(5);
    test.check_rank_dynamic(2);
    test.check_extents(5, 4, 3, 2, 1);
}

TEST(ExtentsTest, PurelyStatic)
{
    ExtentsTest<5, 4, 3> test;
    test.check_rank(3);
    test.check_rank_dynamic(0);
    test.check_extents(5, 4, 3);
}

TEST(ExtentsTest, RankNought)
{
    // Can construct extents of rank nought
    ExtentsTest<> test;
    test.check_rank(0);
    test.check_rank_dynamic(0);
}
TEST(ExtentsTest, Assignment)
{
    extents<5, dynamic_extent, 3, dynamic_extent, 1> e1(4, 2);
    extents<5, 4, 3, 2, 1>                           e2;
    e2 = e1;
    for (size_t r = 0; r < 5; r++)
    {
        EXPECT_EQ(e2.extent(r), e1.extent(r));
    }
    extents<dynamic_extent, dynamic_extent, dynamic_extent, dynamic_extent, dynamic_extent> e3(
            9, 8, 7, 6, 5);
    for (int r = 0; r < 5; r++)
    {
        EXPECT_EQ(e3.extent(r), 9 - r);
    }
    e3 = e1;
    for (int r = 0; r < 5; r++)
    {
        EXPECT_EQ(e3.extent(r), e1.extent(r));
    }
}

} // namespace
} // namespace test
} // namespace gmx
