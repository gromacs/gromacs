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
 * \brief Testing gmx::extents.
 *
 * \author Christian Trott <crtrott@sandia.gov>
 * \author David Hollman <dshollm@sandia.gov>
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "gmxpre.h"

#include "gromacs/mdspan/layouts.h"

#include <cstddef>

#include <array>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/mdspan/extents.h"

namespace gmx
{
namespace test
{
namespace
{

template<class Layout, ptrdiff_t... E_STATIC>
struct LayoutTests
{

    typedef Layout                                          layout_type;
    typedef extents<E_STATIC...>                            extents_type;
    typedef typename Layout::template mapping<extents_type> mapping_type;

    mapping_type my_mapping_explicit, my_mapping_copy;

    template<class... E>
    LayoutTests(E... e)
    {
        my_mapping_explicit = mapping_type(extents_type(e...));
        my_mapping_copy     = mapping_type(my_mapping_explicit);
    }

    void check_rank(ptrdiff_t r)
    {
        EXPECT_EQ(my_mapping_explicit.extents().rank(), r);
        EXPECT_EQ(my_mapping_copy.extents().rank(), r);
    }
    void check_rank_dynamic(ptrdiff_t r)
    {
        EXPECT_EQ(my_mapping_explicit.extents().rank_dynamic(), r);
        EXPECT_EQ(my_mapping_copy.extents().rank_dynamic(), r);
    }
    template<class... E>
    void check_extents(E... e)
    {
        std::array<ptrdiff_t, extents_type::rank()> a = { { e... } };
        for (size_t r = 0; r < extents_type::rank(); r++)
        {
            EXPECT_EQ(my_mapping_explicit.extents().extent(r), a[r]);
            EXPECT_EQ(my_mapping_copy.extents().extent(r), a[r]);
        }
    }
    template<class... E>
    void check_strides(E... e)
    {
        std::array<ptrdiff_t, extents_type::rank()> a = { { e... } };
        for (size_t r = 0; r < extents_type::rank(); r++)
        {
            EXPECT_EQ(my_mapping_explicit.stride(r), a[r]);
            EXPECT_EQ(my_mapping_copy.stride(r), a[r]);
        }
    }
    void check_required_span_size(ptrdiff_t size)
    {
        EXPECT_EQ(my_mapping_explicit.required_span_size(), size);
        EXPECT_EQ(my_mapping_copy.required_span_size(), size);
    }

    void check_properties(bool always_unique,
                          bool always_contiguous,
                          bool always_strided,
                          bool unique,
                          bool contiguous,
                          bool strided)
    {
        EXPECT_EQ(my_mapping_explicit.is_always_unique() ? 1 : 0, always_unique ? 1 : 0);
        EXPECT_EQ(my_mapping_explicit.is_always_contiguous() ? 1 : 0, always_contiguous ? 1 : 0);
        EXPECT_EQ(my_mapping_explicit.is_always_strided() ? 1 : 0, always_strided ? 1 : 0);
        EXPECT_EQ(my_mapping_explicit.is_unique() ? 1 : 0, unique ? 1 : 0);
        EXPECT_EQ(my_mapping_explicit.is_contiguous() ? 1 : 0, contiguous ? 1 : 0);
        EXPECT_EQ(my_mapping_explicit.is_strided() ? 1 : 0, strided ? 1 : 0);
        EXPECT_EQ(my_mapping_copy.is_always_unique() ? 1 : 0, always_unique ? 1 : 0);
        EXPECT_EQ(my_mapping_copy.is_always_contiguous() ? 1 : 0, always_contiguous ? 1 : 0);
        EXPECT_EQ(my_mapping_copy.is_always_strided() ? 1 : 0, always_strided ? 1 : 0);
        EXPECT_EQ(my_mapping_copy.is_unique() ? 1 : 0, unique ? 1 : 0);
        EXPECT_EQ(my_mapping_copy.is_contiguous() ? 1 : 0, contiguous ? 1 : 0);
        EXPECT_EQ(my_mapping_copy.is_strided() ? 1 : 0, strided ? 1 : 0);
    }

    template<class... E>
    void check_operator(ptrdiff_t offset, E... e)
    {
        EXPECT_EQ(my_mapping_explicit(e...), offset);
        EXPECT_EQ(my_mapping_copy(e...), offset);
    }
};

TEST(LayoutTests, LayoutRightConstruction)
{
    LayoutTests<layout_right, 5, dynamic_extent, 3, dynamic_extent, 1> test(4, 2);

    test.check_rank(5);
    test.check_rank_dynamic(2);
    test.check_extents(5, 4, 3, 2, 1);
    test.check_strides(24, 6, 2, 1, 1);
    test.check_required_span_size(120);
}

TEST(LayoutTests, LayoutRightProperties)
{
    LayoutTests<layout_right, 5, dynamic_extent, 3, dynamic_extent, 1> test(4, 2);

    test.check_properties(true, true, true, true, true, true);
}

TEST(LayoutTests, LayoutRightOperator)
{
    LayoutTests<layout_right, 5, dynamic_extent, 3, dynamic_extent, 1> test(4, 2);

    test.check_operator(107, 4, 1, 2, 1, 0);
    test.check_operator(0, 0, 0, 0, 0, 0);
}

} // namespace
} // namespace test
} // namespace gmx
