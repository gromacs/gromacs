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

#include "gromacs/mdspan/mdspan"
#include "gtest/gtest.h"

using namespace std::experimental::fundamentals_v3;

class mdspan_ : public ::testing::Test
{
    protected:
        static void SetUpTestCase()
        {
        }

        static void TearDownTestCase()
        {
        }
};

template<class ElementType, class Mapping>
struct fill_raw_data {
    static void fill(ElementType* p, Mapping m)
    {
        typename Mapping::extents_type e = m.extents();
        for (ptrdiff_t i0 = 0; i0 < e.extent(0); i0++)
        {
            for (ptrdiff_t i1 = 0; i1 < e.extent(1); i1++)
            {
                for (ptrdiff_t i2 = 0; i2 < e.extent(2); i2++)
                {
                    for (ptrdiff_t i3 = 0; i3 < e.extent(3); i3++)
                    {
                        for (ptrdiff_t i4 = 0; i4 < e.extent(4); i4++)
                        {
                            p[i0*m.stride(0)+i1*m.stride(1)+i2*m.stride(2)+i3*m.stride(3)+i4*m.stride(4)] =
                                ElementType(i0*10000+i1*1000+i2*100+i3*10+i4);
                        }
                    }
                }
            }
        }
    }
};

template<class MDSPAN>
struct test_mdspan {

    typedef MDSPAN mdspan_type;
    typedef typename mdspan_type::element_type   element_type;
    typedef typename mdspan_type::extents_type   extents_type;
    typedef typename mdspan_type::mapping_type   mapping_type;
    typedef typename mdspan_type::accessor_type  accessor_type;
    typedef typename mdspan_type::pointer        pointer_type;

    mdspan_type my_mdspan_extents, my_mdspan_array, my_mdspan_mapping, my_mdspan_map_acc, my_mdspan_copy;

    template<class ... ED>
    test_mdspan(ED ... e)
    {
        mapping_type  map(extents_type(e ...));
        accessor_type acc;
        element_type* raw_ptr = new element_type[map.required_span_size()];
        fill_raw_data<element_type, mapping_type>::fill(raw_ptr, map);
        pointer_type  p(raw_ptr);
        my_mdspan_array   = mdspan_type(p, std::array<ptrdiff_t, sizeof ... (ED)>({{e ...}}));
        my_mdspan_mapping = mdspan_type(p, map);
        my_mdspan_map_acc = mdspan_type(p, map, acc);
        my_mdspan_extents = mdspan_type(p, e ...);
        my_mdspan_copy    = my_mdspan_mapping;
    }

    void check_rank(ptrdiff_t r)
    {
        ASSERT_EQ(my_mdspan_mapping.rank(), r);
        ASSERT_EQ(my_mdspan_map_acc.rank(), r);
        ASSERT_EQ(my_mdspan_extents.rank(), r);
        ASSERT_EQ(my_mdspan_copy.rank(), r);
    }
    void check_rank_dynamic(ptrdiff_t r)
    {
        ASSERT_EQ(my_mdspan_mapping.rank_dynamic(), r);
        ASSERT_EQ(my_mdspan_map_acc.rank_dynamic(), r);
        ASSERT_EQ(my_mdspan_extents.rank_dynamic(), r);
        ASSERT_EQ(my_mdspan_copy.rank_dynamic(), r);
    }
    template<class ... E>
    void check_extents(E ... e)
    {
        std::array<ptrdiff_t, extents_type::rank()> a({{e ...}});
        for (size_t r = 0; r < extents_type::rank(); r++)
        {
            ASSERT_EQ(my_mdspan_mapping.extent(r), a[r]);
            ASSERT_EQ(my_mdspan_map_acc.extent(r), a[r]);
            ASSERT_EQ(my_mdspan_extents.extent(r), a[r]);
            ASSERT_EQ(my_mdspan_copy.extent(r), a[r]);
        }
    }
    template<class ... E>
    void check_strides(E ... e)
    {
        std::array<ptrdiff_t, extents_type::rank()> a({{e ...}});
        for (size_t r = 0; r < extents_type::rank(); r++)
        {
            ASSERT_EQ(my_mdspan_mapping.stride(r), a[r]);
            ASSERT_EQ(my_mdspan_map_acc.stride(r), a[r]);
            ASSERT_EQ(my_mdspan_extents.stride(r), a[r]);
            ASSERT_EQ(my_mdspan_copy.stride(r), a[r]);
        }
    }

    void check_properties_internal(mdspan_type my_mdspan, bool always_unique, bool always_contiguous, bool always_strided,
                                   bool unique, bool contiguous, bool strided)
    {
        ASSERT_EQ(my_mdspan.is_always_unique() ? 1 : 0, always_unique ? 1 : 0);
        ASSERT_EQ(my_mdspan.is_always_contiguous() ? 1 : 0, always_contiguous ? 1 : 0);
        ASSERT_EQ(my_mdspan.is_always_strided() ? 1 : 0, always_strided ? 1 : 0);
        ASSERT_EQ(my_mdspan.is_unique() ? 1 : 0, unique ? 1 : 0);
        ASSERT_EQ(my_mdspan.is_contiguous() ? 1 : 0, contiguous ? 1 : 0);
        ASSERT_EQ(my_mdspan.is_strided() ? 1 : 0, strided ? 1 : 0);
    }

    void check_properties(bool always_unique, bool always_contiguous, bool always_strided,
                          bool unique, bool contiguous, bool strided)
    {
        check_properties_internal(my_mdspan_mapping, always_unique, always_contiguous, always_strided, unique, contiguous, strided);
        check_properties_internal(my_mdspan_map_acc, always_unique, always_contiguous, always_strided, unique, contiguous, strided);
        check_properties_internal(my_mdspan_extents, always_unique, always_contiguous, always_strided, unique, contiguous, strided);
        check_properties_internal(my_mdspan_copy, always_unique, always_contiguous, always_strided, unique, contiguous, strided);
    }

    void check_operator()
    {
        extents_type e = my_mdspan_mapping.extents();
        for (ptrdiff_t i0 = 0; i0 < e.extent(0); i0++)
        {
            for (ptrdiff_t i1 = 0; i1 < e.extent(1); i1++)
            {
                for (ptrdiff_t i2 = 0; i2 < e.extent(2); i2++)
                {
                    for (ptrdiff_t i3 = 0; i3 < e.extent(3); i3++)
                    {
                        for (ptrdiff_t i4 = 0; i4 < e.extent(4); i4++)
                        {
                            element_type value = i0*10000+i1*1000+i2*100+i3*10+i4;
                            ASSERT_EQ(my_mdspan_mapping(i0, i1, i2, i3, i4), value);
                            ASSERT_EQ(my_mdspan_map_acc(i0, i1, i2, i3, i4), value);
                            ASSERT_EQ(my_mdspan_extents(i0, i1, i2, i3, i4), value);
                            ASSERT_EQ(my_mdspan_copy   (i0, i1, i2, i3, i4), value);
                        }
                    }
                }
            }
        }
    }
};

TEST_F(mdspan_, construction_right) {
    typedef basic_mdspan<int, extents<5, dynamic_extent, 3, dynamic_extent, 1>,
                         layout_right, accessor_basic<int> > mdspan_type;
    test_mdspan<mdspan_type> test(4, 2);

    test.check_rank(5);
//  test.check_rank_dynamic(2);
//  test.check_extents(5,4,3,2,1);
//  test.check_strides(24,6,2,1,1);
}

TEST_F(mdspan_, construction_left) {
    typedef basic_mdspan<int, extents<5, dynamic_extent, 3, dynamic_extent, 1>,
                         layout_left, accessor_basic<int> > mdspan_type;
    test_mdspan<mdspan_type> test(4, 2);

    test.check_rank(5);
    test.check_rank_dynamic(2);
    test.check_extents(5, 4, 3, 2, 1);
    test.check_strides(1, 5, 20, 60, 120);
}


TEST_F(mdspan_, properties_right) {
    typedef basic_mdspan<int, extents<5, dynamic_extent, 3, dynamic_extent, 1>,
                         layout_right, accessor_basic<int> > mdspan_type;
    test_mdspan<mdspan_type> test(4, 2);

    test.check_properties(true, true, true, true, true, true);
}

TEST_F(mdspan_, properties_left) {
    typedef basic_mdspan<int, extents<5, dynamic_extent, 3, dynamic_extent, 1>,
                         layout_left, accessor_basic<int> > mdspan_type;
    test_mdspan<mdspan_type> test(4, 2);

    test.check_properties(true, true, true, true, true, true);
}

TEST_F(mdspan_, operator_right) {
    typedef basic_mdspan<int, extents<5, dynamic_extent, 3, dynamic_extent, 1>,
                         layout_right, accessor_basic<int> > mdspan_type;
    test_mdspan<mdspan_type> test(4, 2);

    test.check_operator();
}

TEST_F(mdspan_, operator_left) {
    typedef basic_mdspan<int, extents<5, dynamic_extent, 3, dynamic_extent, 1>,
                         layout_left, accessor_basic<int> > mdspan_type;
    test_mdspan<mdspan_type> test(4, 2);

    test.check_operator();
}
