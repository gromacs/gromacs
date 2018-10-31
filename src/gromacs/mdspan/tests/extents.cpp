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

#include "gromacs/mdspan/extents.h"
#include "gtest/gtest.h"

using std::experimental::fundamentals_v3::extents;
using std::experimental::fundamentals_v3::dynamic_extent;

class extents_ : public ::testing::Test
{
    protected:
        static void SetUpTestCase()
        {
        }

        static void TearDownTestCase()
        {
        }
};

template<ptrdiff_t ... E_STATIC>
struct test_extents {

    typedef extents<E_STATIC...> extents_type;

    extents_type my_extents_explicit, my_extents_array, my_extents_copy;

    test_extents()
    {
        my_extents_explicit = extents<E_STATIC...>();
        my_extents_array    = extents<E_STATIC...>(std::array<ptrdiff_t, 0>());
        my_extents_copy     = extents<E_STATIC...>(my_extents_explicit);
    }

    template<class ... E>
    test_extents(E ... e)
    {
        my_extents_explicit = extents<E_STATIC...>(e ...);
        my_extents_array    = extents<E_STATIC...>(std::array<ptrdiff_t, 2>({{e ...}}));
        my_extents_copy     = extents<E_STATIC...>(my_extents_explicit);
    }

    void check_rank(ptrdiff_t r)
    {
        ASSERT_EQ(my_extents_explicit.rank(), r);
        ASSERT_EQ(my_extents_array.rank(), r);
        ASSERT_EQ(my_extents_copy.rank(), r);
    }
    void check_rank_dynamic(ptrdiff_t r)
    {
        ASSERT_EQ(my_extents_explicit.rank_dynamic(), r);
        ASSERT_EQ(my_extents_array.rank_dynamic(), r);
        ASSERT_EQ(my_extents_copy.rank_dynamic(), r);
    }
    template<class ... E>
    void check_extents(E ... e)
    {
        std::array<ptrdiff_t, extents_type::rank()> s({{E_STATIC ...}});
        std::array<ptrdiff_t, extents_type::rank()> a({{e ...}});
        for (size_t r = 0; r < extents_type::rank(); r++)
        {
            ASSERT_EQ(my_extents_explicit.static_extent(r), s[r]);
            ASSERT_EQ(my_extents_explicit.extent(r), a[r]);

            ASSERT_EQ(my_extents_array.static_extent(r), s[r]);
            ASSERT_EQ(my_extents_array.extent(r), a[r]);

            ASSERT_EQ(my_extents_copy.static_extent(r), s[r]);
            ASSERT_EQ(my_extents_copy.extent(r), a[r]);
        }
        ASSERT_EQ(my_extents_explicit.static_extent(extents_type::rank()+1), 1);
        ASSERT_EQ(my_extents_explicit.extent(extents_type::rank()+1), 1);

        ASSERT_EQ(my_extents_array.static_extent(extents_type::rank()+1), 1);
        ASSERT_EQ(my_extents_array.extent(extents_type::rank()+1), 1);

        ASSERT_EQ(my_extents_copy.static_extent(extents_type::rank()+1), 1);
        ASSERT_EQ(my_extents_copy.extent(extents_type::rank()+1), 1);
    }

};

TEST_F(extents_, construction) {
    test_extents<5, dynamic_extent, 3, dynamic_extent, 1> test(4, 2);

    test.check_rank(5);
    test.check_rank_dynamic(2);
    test.check_extents(5, 4, 3, 2, 1);

}

TEST_F(extents_, static_only) {
    test_extents<5, 4, 3> test;
    test.check_rank(3);
    test.check_rank_dynamic(0);
    test.check_extents(5, 4, 3);
}

TEST_F(extents_, rank_0) {
    test_extents<> test;
    test.check_rank(0);
    test.check_rank_dynamic(0);
}
TEST_F(extents_, assignment) {
    extents<5, dynamic_extent, 3, dynamic_extent, 1> e1(4, 2);
    extents<5, 4, 3, 2, 1> e2;
    e2 = e1;
    for (size_t r = 0; r < 5; r++)
    {
        ASSERT_EQ(e2.extent(r), e1.extent(r));
    }
    extents<dynamic_extent, dynamic_extent, dynamic_extent, dynamic_extent, dynamic_extent> e3(9, 8, 7, 6, 5);
    for (int r = 0; r < 5; r++)
    {
        ASSERT_EQ(e3.extent(r), 9-r);
    }
    e3 = e1;
    for (int r = 0; r < 5; r++)
    {
        ASSERT_EQ(e3.extent(r), e1.extent(r));
    }
}
