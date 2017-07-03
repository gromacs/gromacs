/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Tests for GPU page-locked memory.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/pagelockedmemory.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "gputest.h"

namespace gmx
{

// Doxygen doesn't work for extern templates
#if !defined DOXYGEN

extern template PageLockedMemory::PageLockedMemory(ConstArrayRef<real>);
extern template PageLockedMemory::PageLockedMemory(ConstArrayRef<RVec>);

#endif

namespace
{

//! Typed test fixture
template <typename T>
class PageLockedMemoryTest : public test::GpuTest
{
    public:
        //! Convenience type
        using ValueType = T;
        //! Convenience type
        using AllocatorType = Allocator<ValueType, PageAlignedAllocationPolicy>;
        //! Convenience type
        using VectorType = std::vector<ValueType, AllocatorType>;
        // TODO The implementation of ArrayRef hard-codes std::vector,
        // which makes it hard to get an ArrayRef from a std::vector with
        // an allocator, because that is a different type.
        //! Convenience type
        using ConstArrayRefType = ConstArrayRef<ValueType>;
        //! Convenience conversion function.
        ConstArrayRefType asArrayRef(const VectorType &v)
        {
            return ConstArrayRefType::fromArray(v.data(), v.size());
        }
};

//! The types used in testing.
typedef ::testing::Types<real, RVec> TestTypes;

TYPED_TEST_CASE(PageLockedMemoryTest, TestTypes);

// Note that the use of TestFixture:: and this-> is sometimes required
// to get access to things in the fixture class (or its base classes).

TYPED_TEST(PageLockedMemoryTest, EmptyMemoryAlwaysWorks)
{
    typename TestFixture::ConstArrayRefType empty = EmptyArrayRef();
    PageLockedMemory lockedMemory(empty);
    EXPECT_EQ(nullptr, lockedMemory.memory());
}

TYPED_TEST(PageLockedMemoryTest, LocksAndUnlocksPageAlignedMemory)
{
    if (!this->haveValidGpus())
    {
        return;
    }

    typename TestFixture::VectorType host(1);
    PageLockedMemory                 lockedMemory(this->asArrayRef(host));
}

TYPED_TEST(PageLockedMemoryTest, CanSwapEmptyPageLockedMemory)
{
    if (!this->haveValidGpus())
    {
        return;
    }

    // Using 3-element initializer lists works for both normal and
    // BasicVector<normal> types.
    typename TestFixture::VectorType
    host {{
              2, 3, 4
          }};
    typename TestFixture::ConstArrayRefType empty = EmptyArrayRef();
    PageLockedMemory lockedMemory(this->asArrayRef(host)), initiallyEmpty(empty);
    initiallyEmpty.swap(lockedMemory);
    for (size_t i = 0; i != host.size(); ++i)
    {
        EXPECT_EQ(nullptr, lockedMemory.memory());
        EXPECT_EQ(host.data(), initiallyEmpty.memory());
    }
}

TYPED_TEST(PageLockedMemoryTest, CanSwapNonEmptyPageLockedMemory)
{
    if (!this->haveValidGpus())
    {
        return;
    }

    typename TestFixture::VectorType
    firstHost {{
                   2, 3, 4
               }};
    typename TestFixture::VectorType
    secondHost {{
                    4, 5, 6
                }};
    PageLockedMemory firstLockedMemory(this->asArrayRef(firstHost)), secondLockedMemory(this->asArrayRef(secondHost));
    secondLockedMemory.swap(firstLockedMemory);
    for (size_t i = 0; i != firstHost.size(); ++i)
    {
        EXPECT_EQ(secondHost.data(), firstLockedMemory.memory());
        EXPECT_EQ(firstHost.data(), secondLockedMemory.memory());
    }
}

} // namespace
} // namespace
