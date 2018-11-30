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
/*
 *                          Kokkos v. 2.0
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
 * \brief Testing gmx::mdspan
 *
 * \author Christian Trott <crtrott@sandia.gov>
 * \author Christian Blau <cblau@gwdg.de>
 */

#include "gmxpre.h"

#include "gromacs/mdspan/subspan.h"

#include <cstdio>

#include <gtest/gtest.h>

namespace gmx
{

class SubSpanTest : public ::testing::Test
{
    protected:
        template <typename Map>
        std::array<ptrdiff_t, Map::extents_type::rank()> copyStridesFromMap(const Map &map)
        {
            std::array<ptrdiff_t, Map::extents_type::rank()> result;
            for (size_t r = 0; r < map.extents().rank(); r++)
            {
                result [r] = map.stride(r);
            }
            return result;
        }

        using extents_type           = extents<5, dynamic_extent, 3, dynamic_extent, 1>;
        using sub_extent_deduce_type = detail::subspan_deduce_extents<extents_type, ptrdiff_t, std::pair<int, int>, all_type, all_type, ptrdiff_t>;
        using new_extents_type       = sub_extent_deduce_type::extents_type;

        extents_type e      = {4, 2};
        ptrdiff_t    offset = 0;
        std::array<ptrdiff_t, extents_type::rank()> strides;
};


TEST_F(SubSpanTest, static_extent_deduction) {
    EXPECT_EQ(new_extents_type::rank(), 3);
    EXPECT_EQ(new_extents_type::rank_dynamic(), 2);
    EXPECT_EQ(new_extents_type::static_extent(0), dynamic_extent);
    EXPECT_EQ(new_extents_type::static_extent(1), 3);
    EXPECT_EQ(new_extents_type::static_extent(2), dynamic_extent);
}

TEST_F(SubSpanTest, dynamic_extent_deduction) {


    new_extents_type sub_extents = sub_extent_deduce_type::create_sub_extents(e, strides, offset, 2, std::pair<int, int>(1, 3), all_type(), all_type(), 0);
    EXPECT_EQ(sub_extents.static_extent(0), dynamic_extent);
    EXPECT_EQ(sub_extents.static_extent(1), 3);
    EXPECT_EQ(sub_extents.static_extent(2), dynamic_extent);
    EXPECT_EQ(sub_extents.extent(0), 2);
    EXPECT_EQ(sub_extents.extent(1), 3);
    EXPECT_EQ(sub_extents.extent(2), 2);
}

TEST_F(SubSpanTest, strides_deduction_layout_right) {

    layout_right::mapping<extents_type> map(e);
    strides = copyStridesFromMap(map);

    sub_extent_deduce_type::create_sub_extents(e, strides, offset, 2, std::pair<int, int>(1, 3), all_type(), all_type(), 0);
    EXPECT_EQ(strides[0], map.stride(1));
    EXPECT_EQ(strides[1], map.stride(2));
    EXPECT_EQ(strides[2], map.stride(3));
}

TEST_F(SubSpanTest, offset_deduction_layout_right) {

    layout_right::mapping<extents_type> map(e);
    strides = copyStridesFromMap(map);

    sub_extent_deduce_type::create_sub_extents(e, strides, offset, 2, std::pair<int, int>(1, 3), all_type(), all_type(), 0);
    EXPECT_EQ(offset, 2*map.stride(0)+map.stride(1));
}

TEST_F(SubSpanTest, basic_mdspan_layout_right) {

    layout_right::mapping<extents_type> map(e);
    std::vector<int>                    ptr(map.required_span_size());
    using mdspan_type = basic_mdspan < int, extents_type, layout_right, accessor_basic < int>>;

    mdspan_type a(ptr.data(), e);

    auto        sub = subspan(a, ptrdiff_t(2), std::pair<int, int>(1, 3), all_type(), all_type(), ptrdiff_t(0));
    EXPECT_EQ(sub.rank(), 3);
    EXPECT_EQ(sub.rank_dynamic(), 2);
    EXPECT_EQ(sub.static_extent(0), dynamic_extent);
    EXPECT_EQ(sub.static_extent(1), 3);
    EXPECT_EQ(sub.static_extent(2), dynamic_extent);
    EXPECT_EQ(sub.extent(0), 2);
    EXPECT_EQ(sub.extent(1), 3);
    EXPECT_EQ(sub.extent(2), 2);
    EXPECT_EQ(sub.is_contiguous() ? 1 : 0, 1);
    EXPECT_EQ(sub.stride(0), a.stride(1));
    EXPECT_EQ(sub.stride(1), a.stride(2));
    EXPECT_EQ(sub.stride(2), a.stride(3));

    for (int i0 = 0; i0 < a.extent(0); i0++)
    {
        for (int i1 = 0; i1 < a.extent(1); i1++)
        {
            for (int i2 = 0; i2 < a.extent(2); i2++)
            {
                for (int i3 = 0; i3 < a.extent(3); i3++)
                {
                    for (int i4 = 0; i4 < a.extent(4); i4++)
                    {
                        a(i0, i1, i2, i3, i4) = i0*10000+i1*1000+i2*100+i3*10+i4;
                    }
                }
            }
        }
    }

    for (int i0 = 0; i0 < sub.extent(0); i0++)
    {
        for (int i1 = 0; i1 < sub.extent(1); i1++)
        {
            for (int i2 = 0; i2 < sub.extent(2); i2++)
            {
                EXPECT_EQ(sub(i0, i1, i2), 2*10000 + (i0+1)*1000+i1*100+i2*10);
            }
        }
    }
}

} // namespace gmx
