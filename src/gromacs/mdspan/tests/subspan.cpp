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
#include <cstdio>

#include "gromacs/mdspan/subspan.h"
#include "gtest/gtest.h"

using namespace std::experimental::fundamentals_v3;

class subspan_ : public ::testing::Test
{
    protected:
        static void SetUpTestCase()
        {
        }

        static void TearDownTestCase()
        {
        }
};

TEST_F(subspan_, static_extent_deduction) {
    typedef extents<5, dynamic_extent, 3, dynamic_extent, 1> extents_type;
    typedef detail::subspan_deduce_extents<extents_type, ptrdiff_t, std::pair<int, int>, all_type, all_type, ptrdiff_t>::extents_type new_extents_type;

    ASSERT_EQ(new_extents_type::rank(), 3);
    ASSERT_EQ(new_extents_type::rank_dynamic(), 2);
    ASSERT_EQ(new_extents_type::static_extent(0), dynamic_extent);
    ASSERT_EQ(new_extents_type::static_extent(1), 3);
    ASSERT_EQ(new_extents_type::static_extent(2), dynamic_extent);
}

TEST_F(subspan_, dynamic_extent_deduction) {
    typedef extents<5, dynamic_extent, 3, dynamic_extent, 1> extents_type;
    typedef detail::subspan_deduce_extents<extents_type, ptrdiff_t, std::pair<int, int>, all_type, all_type, ptrdiff_t> sub_extent_deduce_type;
    typedef sub_extent_deduce_type::extents_type new_extents_type;

    std::array<ptrdiff_t, extents_type::rank()> strides;
    ptrdiff_t        offset;
    extents_type     e(4, 2);

    new_extents_type sub_extents = sub_extent_deduce_type::create_sub_extents(e, strides, offset, 2, std::pair<int, int>(1, 3), all_type(), all_type(), 0);
    ASSERT_EQ(sub_extents.static_extent(0), dynamic_extent);
    ASSERT_EQ(sub_extents.static_extent(1), 3);
    ASSERT_EQ(sub_extents.static_extent(2), dynamic_extent);
    ASSERT_EQ(sub_extents.extent(0), 2);
    ASSERT_EQ(sub_extents.extent(1), 3);
    ASSERT_EQ(sub_extents.extent(2), 2);
}

TEST_F(subspan_, strides_deduction_layout_right) {
    typedef extents<5, dynamic_extent, 3, dynamic_extent, 1> extents_type;
    typedef detail::subspan_deduce_extents<extents_type, ptrdiff_t, std::pair<int, int>, all_type, all_type, ptrdiff_t> sub_extent_deduce_type;
    //typedef sub_extent_deduce_type::extents_type new_extents_type;

    std::array<ptrdiff_t, extents_type::rank()> strides;
    ptrdiff_t    offset = 0;
    extents_type e(4, 2);
    layout_right::mapping<extents_type> map(e);
    for (size_t r = 0; r < e.rank(); r++)
    {
        strides[r] = map.stride(r);
    }

    sub_extent_deduce_type::create_sub_extents(e, strides, offset, 2, std::pair<int, int>(1, 3), all_type(), all_type(), 0);
    ASSERT_EQ(strides[0], map.stride(1));
    ASSERT_EQ(strides[1], map.stride(2));
    ASSERT_EQ(strides[2], map.stride(3));
}


TEST_F(subspan_, strides_deduction_layout_left) {
    typedef extents<5, dynamic_extent, 3, dynamic_extent, 1> extents_type;
    typedef detail::subspan_deduce_extents<extents_type, ptrdiff_t, std::pair<int, int>, all_type, all_type, ptrdiff_t> sub_extent_deduce_type;
    //typedef sub_extent_deduce_type::extents_type new_extents_type;

    std::array<ptrdiff_t, extents_type::rank()> strides;
    ptrdiff_t    offset = 0;
    extents_type e(4, 2);
    layout_left::mapping<extents_type> map(e);
    for (size_t r = 0; r < e.rank(); r++)
    {
        strides[r] = map.stride(r);
    }

    sub_extent_deduce_type::create_sub_extents(e, strides, offset, 2, std::pair<int, int>(1, 3), all_type(), all_type(), 0);
    ASSERT_EQ(strides[0], map.stride(1));
    ASSERT_EQ(strides[1], map.stride(2));
    ASSERT_EQ(strides[2], map.stride(3));
}

TEST_F(subspan_, offset_deduction_layout_right) {
    typedef extents<5, dynamic_extent, 3, dynamic_extent, 1> extents_type;
    typedef detail::subspan_deduce_extents<extents_type, ptrdiff_t, std::pair<int, int>, all_type, all_type, ptrdiff_t> sub_extent_deduce_type;
    //typedef sub_extent_deduce_type::extents_type new_extents_type;

    std::array<ptrdiff_t, extents_type::rank()> strides;
    ptrdiff_t    offset = 0;
    extents_type e(4, 2);
    layout_right::mapping<extents_type> map(e);
    for (size_t r = 0; r < e.rank(); r++)
    {
        strides[r] = map.stride(r);
    }

    sub_extent_deduce_type::create_sub_extents(e, strides, offset, 2, std::pair<int, int>(1, 3), all_type(), all_type(), 0);
    ASSERT_EQ(offset, 2*map.stride(0)+map.stride(1));
}

TEST_F(subspan_, offset_deduction_layout_left) {
    typedef extents<5, dynamic_extent, 3, dynamic_extent, 1> extents_type;
    typedef detail::subspan_deduce_extents<extents_type, ptrdiff_t, std::pair<int, int>, all_type, all_type, ptrdiff_t> sub_extent_deduce_type;
    //typedef sub_extent_deduce_type::extents_type new_extents_type;

    std::array<ptrdiff_t, extents_type::rank()> strides;
    ptrdiff_t    offset = 0;
    extents_type e(4, 2);
    layout_left::mapping<extents_type> map(e);
    for (size_t r = 0; r < e.rank(); r++)
    {
        strides[r] = map.stride(r);
    }

    sub_extent_deduce_type::create_sub_extents(e, strides, offset, 2, std::pair<int, int>(1, 3), all_type(), all_type(), 0);
    ASSERT_EQ(offset, 2*map.stride(0)+map.stride(1));
}

TEST_F(subspan_, basic_mdspan_layout_right) {
    typedef extents<5, dynamic_extent, 3, dynamic_extent, 1> extents_type;
    extents_type e(4, 2);
    layout_right::mapping<extents_type> map(e);
    int        * ptr = new int[map.required_span_size()];
    typedef basic_mdspan<int, extents_type, layout_right, accessor_basic<int> > mdspan_type;

    mdspan_type a(ptr, e);

    auto        sub = subspan(a, ptrdiff_t(2), std::pair<int, int>(1, 3), all_type(), all_type(), ptrdiff_t(0));
    ASSERT_EQ(sub.rank(), 3);
    ASSERT_EQ(sub.rank_dynamic(), 2);
    ASSERT_EQ(sub.static_extent(0), dynamic_extent);
    ASSERT_EQ(sub.static_extent(1), 3);
    ASSERT_EQ(sub.static_extent(2), dynamic_extent);
    ASSERT_EQ(sub.extent(0), 2);
    ASSERT_EQ(sub.extent(1), 3);
    ASSERT_EQ(sub.extent(2), 2);
    ASSERT_EQ(sub.is_contiguous() ? 1 : 0, 1);
    ASSERT_EQ(sub.stride(0), a.stride(1));
    ASSERT_EQ(sub.stride(1), a.stride(2));
    ASSERT_EQ(sub.stride(2), a.stride(3));

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
                ASSERT_EQ(sub(i0, i1, i2), 2*10000 + (i0+1)*1000+i1*100+i2*10);
            }
        }
    }
    delete [] ptr;
}

TEST_F(subspan_, basic_mdspan_layout_left) {
    typedef extents<5, dynamic_extent, 3, dynamic_extent, 1> extents_type;
    extents_type e(4, 2);
    layout_left::mapping<extents_type> map(e);
    int        * ptr = new int[map.required_span_size()];
    typedef basic_mdspan<int, extents_type, layout_left, accessor_basic<int> > mdspan_type;

    mdspan_type a(ptr, e);

    auto        sub = subspan(a, ptrdiff_t(2), std::pair<int, int>(1, 3), all_type(), all_type(), ptrdiff_t(0));
    ASSERT_EQ(sub.rank(), 3);
    ASSERT_EQ(sub.rank_dynamic(), 2);
    ASSERT_EQ(sub.static_extent(0), dynamic_extent);
    ASSERT_EQ(sub.static_extent(1), 3);
    ASSERT_EQ(sub.static_extent(2), dynamic_extent);
    ASSERT_EQ(sub.extent(0), 2);
    ASSERT_EQ(sub.extent(1), 3);
    ASSERT_EQ(sub.extent(2), 2);
    ASSERT_EQ(sub.is_contiguous() ? 1 : 0, 0);
    ASSERT_EQ(sub.stride(0), a.stride(1));
    ASSERT_EQ(sub.stride(1), a.stride(2));
    ASSERT_EQ(sub.stride(2), a.stride(3));

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
                ASSERT_EQ(sub(i0, i1, i2), 2*10000 + (i0+1)*1000+i1*100+i2*10);
            }
        }
    }
    delete [] ptr;
}

//TEST_F(subspan_,reduce_to_rank_0) {
//}
