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
 * \brief Testing gmx::basic_mdspan.
 *
 * \author Christian Trott <crtrott@sandia.gov>
 * \author Carter Edwards <hedwards@nvidia.com>
 * \author David Hollman <dshollm@sandia.gov>
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "gmxpre.h"

#include "gromacs/mdspan/mdspan.h"

#include <cstddef>
#include <cstdio>

#include <array>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/mdspan/accessor_policy.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"

namespace gmx
{
namespace test
{
namespace
{


// Test basic_mdspan with a mixture of dynamic and static extents, as well as extent of one.
// The dynamic extents will be set to 4 and 2 repsectively so that we'll tests
// a multidimensional array of extent 5 x 4 x 3 x 2 x 1, i.e. 120 elements

//! View on int data with mixed static and dynamic extents
using mdspan_int =
        basic_mdspan<int, extents<5, dynamic_extent, 3, dynamic_extent, 1>, layout_right, accessor_basic<int>>;

TEST(MdSpanTest, MdSpanWrapsBasicMdSpanCorrectly)
{
    // Check that mdspan wraps basic_mdspan as expected
    ::testing::StaticAssertTypeEq<mdspan_int, mdspan<int, 5, dynamic_extent, 3, dynamic_extent, 1>>();
}

//! View on float data with mixed static and dynamic extents
using mdspan_float =
        basic_mdspan<float, extents<5, dynamic_extent, 3, dynamic_extent, 1>, layout_right, accessor_basic<float>>;
//! Types to be tested
using MdSpanTypes = ::testing::Types<mdspan_int, mdspan_float>;

template<class ElementType, class Mapping>
struct fill_raw_data
{
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
                            p[i0 * m.stride(0) + i1 * m.stride(1) + i2 * m.stride(2)
                              + i3 * m.stride(3) + i4 * m.stride(4)] =
                                    ElementType(i0 * 10000 + i1 * 1000 + i2 * 100 + i3 * 10 + i4);
                        }
                    }
                }
            }
        }
    }
};
template<class MDSPAN>
struct MdSpanTest : public ::testing::Test
{
    using mdspan_type   = MDSPAN;
    using element_type  = typename mdspan_type::element_type;
    using extents_type  = typename mdspan_type::extents_type;
    using mapping_type  = typename mdspan_type::mapping_type;
    using accessor_type = typename mdspan_type::accessor_type;
    using pointer_type  = typename mdspan_type::pointer;

    mdspan_type my_mdspan_extents;
    mdspan_type my_mdspan_array;
    mdspan_type my_mdspan_mapping;
    mdspan_type my_mdspan_map_acc;
    mdspan_type my_mdspan_copy;

    std::vector<element_type> rawData;

    template<class... ED>
    void SetUp(ED... e)
    {
        mapping_type  map{ extents_type(e...) };
        accessor_type acc;
        rawData.resize(map.required_span_size());
        fill_raw_data<element_type, mapping_type>::fill(rawData.data(), map);
        pointer_type p(rawData.data());

        my_mdspan_array   = mdspan_type(p, std::array<ptrdiff_t, sizeof...(ED)>({ { e... } }));
        my_mdspan_mapping = mdspan_type(p, map);
        my_mdspan_map_acc = mdspan_type(p, map, acc);
        my_mdspan_extents = mdspan_type(p, e...);
        my_mdspan_copy    = my_mdspan_mapping;
    }

    void check_rank(ptrdiff_t r)
    {
        EXPECT_EQ(my_mdspan_mapping.rank(), r);
        EXPECT_EQ(my_mdspan_map_acc.rank(), r);
        EXPECT_EQ(my_mdspan_extents.rank(), r);
        EXPECT_EQ(my_mdspan_copy.rank(), r);
    }
    void check_rank_dynamic(ptrdiff_t r)
    {
        EXPECT_EQ(my_mdspan_mapping.rank_dynamic(), r);
        EXPECT_EQ(my_mdspan_map_acc.rank_dynamic(), r);
        EXPECT_EQ(my_mdspan_extents.rank_dynamic(), r);
        EXPECT_EQ(my_mdspan_copy.rank_dynamic(), r);
    }
    template<class... E>
    void check_extents(E... e)
    {
        std::array<ptrdiff_t, extents_type::rank()> a{ { e... } };
        for (size_t r = 0; r < extents_type::rank(); r++)
        {
            EXPECT_EQ(my_mdspan_mapping.extent(r), a[r]);
            EXPECT_EQ(my_mdspan_map_acc.extent(r), a[r]);
            EXPECT_EQ(my_mdspan_extents.extent(r), a[r]);
            EXPECT_EQ(my_mdspan_copy.extent(r), a[r]);
        }
    }
    template<class... E>
    void check_strides(E... e)
    {
        std::array<ptrdiff_t, extents_type::rank()> a{ { e... } };
        for (size_t r = 0; r < extents_type::rank(); r++)
        {
            EXPECT_EQ(my_mdspan_mapping.stride(r), a[r]);
            EXPECT_EQ(my_mdspan_map_acc.stride(r), a[r]);
            EXPECT_EQ(my_mdspan_extents.stride(r), a[r]);
            EXPECT_EQ(my_mdspan_copy.stride(r), a[r]);
        }
    }

    void check_properties_internal(mdspan_type my_mdspan,
                                   bool        always_unique,
                                   bool        always_contiguous,
                                   bool        always_strided,
                                   bool        unique,
                                   bool        contiguous,
                                   bool        strided)
    {
        EXPECT_EQ(my_mdspan.is_always_unique(), always_unique);
        EXPECT_EQ(my_mdspan.is_always_contiguous(), always_contiguous);
        EXPECT_EQ(my_mdspan.is_always_strided(), always_strided);
        EXPECT_EQ(my_mdspan.is_unique(), unique);
        EXPECT_EQ(my_mdspan.is_contiguous(), contiguous);
        EXPECT_EQ(my_mdspan.is_strided(), strided);
    }

    void check_properties(bool always_unique,
                          bool always_contiguous,
                          bool always_strided,
                          bool unique,
                          bool contiguous,
                          bool strided)
    {
        check_properties_internal(
                my_mdspan_mapping, always_unique, always_contiguous, always_strided, unique, contiguous, strided);
        check_properties_internal(
                my_mdspan_map_acc, always_unique, always_contiguous, always_strided, unique, contiguous, strided);
        check_properties_internal(
                my_mdspan_extents, always_unique, always_contiguous, always_strided, unique, contiguous, strided);
        check_properties_internal(
                my_mdspan_copy, always_unique, always_contiguous, always_strided, unique, contiguous, strided);
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
                            element_type value = i0 * 10000 + i1 * 1000 + i2 * 100 + i3 * 10 + i4;
                            EXPECT_EQ(my_mdspan_mapping(i0, i1, i2, i3, i4), value);
                            EXPECT_EQ(my_mdspan_map_acc(i0, i1, i2, i3, i4), value);
                            EXPECT_EQ(my_mdspan_extents(i0, i1, i2, i3, i4), value);
                            EXPECT_EQ(my_mdspan_copy(i0, i1, i2, i3, i4), value);
                        }
                    }
                }
            }
        }
    }
};

TYPED_TEST_SUITE(MdSpanTest, MdSpanTypes);

TYPED_TEST(MdSpanTest, Rank)
{
    this->SetUp(4, 2);
    this->check_rank(5);
}

TYPED_TEST(MdSpanTest, DynamicRank)
{
    this->SetUp(4, 2);
    this->check_rank_dynamic(2);
}

TYPED_TEST(MdSpanTest, Extents)
{
    this->SetUp(4, 2);
    this->check_extents(5, 4, 3, 2, 1);
}

TYPED_TEST(MdSpanTest, Strides)
{
    this->SetUp(4, 2);
    this->check_strides(24, 6, 2, 1, 1);
}

TYPED_TEST(MdSpanTest, Properties)
{
    this->SetUp(4, 2);
    const bool always_unique     = true;
    const bool always_contiguous = true;
    const bool always_strided    = true;
    const bool unique            = true;
    const bool contiguous        = true;
    const bool strided           = true;
    this->check_properties(always_unique, always_contiguous, always_strided, unique, contiguous, strided);
}

TYPED_TEST(MdSpanTest, Operator)
{
    this->SetUp(4, 2);
    this->check_operator();
}

} // namespace
} // namespace test
} // namespace gmx
