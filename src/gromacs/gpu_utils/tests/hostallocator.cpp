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
 * Tests for GPU host allocator.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/hostallocator.h"

#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "devicetransfers.h"
#include "gputest.h"

namespace gmx
{

namespace
{

//! The types used in testing.
typedef ::testing::Types<int, real, RVec> TestTypes;

//! Typed test fixture
template <typename T>
class HostAllocatorTest : public test::GpuTest
{
    public:
        //! Convenience type
        using ValueType = T;
        //! Convenience type
        using AllocatorType = HostAllocator<T>;
        //! Convenience type
        using VectorType = std::vector<ValueType, AllocatorType>;
        //! Convenience type
        using ViewType = ArrayRef<ValueType>;
        //! Convenience type
        using ConstViewType = ArrayRef<const ValueType>;
        //! Prepare contents of a VectorType.
        void fillInput(VectorType *input) const;
        //! Compares input and output vectors.
        void compareVectors(ConstViewType input,
                            ConstViewType output) const;
        //! Do some transfers and test the results.
        void runTest(ConstViewType input, ViewType output) const;
};

// Already documented
template <typename T>
void HostAllocatorTest<T>::fillInput(VectorType *input) const
{
    input->push_back(1);
    input->push_back(2);
    input->push_back(3);
}

//! Initialization specialization for RVec
template <>
void HostAllocatorTest<RVec>::fillInput(VectorType *input) const
{
    input->push_back({1, 2, 3});
}

// Already documented
template <typename T>
void HostAllocatorTest<T>::compareVectors(ConstViewType input,
                                          ConstViewType output) const
{
    for (size_t i = 0; i != input.size(); ++i)
    {
        EXPECT_EQ(input[i], output[i]) << "for index " << i;
    }
}

//! Comparison specialization for RVec
template <>
void HostAllocatorTest<RVec>::compareVectors(ConstViewType input,
                                             ConstViewType output) const
{
    for (size_t i = 0; i != input.size(); ++i)
    {
        EXPECT_EQ(input[i][XX], output[i][XX]) << "for index " << i;
        EXPECT_EQ(input[i][YY], output[i][YY]) << "for index " << i;
        EXPECT_EQ(input[i][ZZ], output[i][ZZ]) << "for index " << i;
    }
}

/*! \brief Convenience function to transform a view into one with base
 * type of (non-const) char.
 *
 * This transformation is useful for using containers with C APIs
 * where the function signature is not declared const even where the
 * semantics of the usage actually are const.
 *
 * \param[in]    data   The data pointer.
 * \param[in]    size   The size of the data pointer (in T).
 * \tparam       T      The base type of the container
 * */
template <typename T>
ArrayRef<char> charArrayRefFromArray(T *data, size_t size)
{
    // Make a type like T, but without its possible const qualifier.
    using NonConstT = typename std::remove_const<T>::type;
    return arrayRefFromArray<char>(reinterpret_cast<char *>(const_cast<NonConstT *>(data)), size * sizeof(T));
}

template <typename T>
void HostAllocatorTest<T>::runTest(ConstViewType input, ViewType output) const
{
    // We can't do a test that does a transfer unless we have a
    // compatible device.
    if (!this->haveValidGpus())
    {
        return;
    }

    // Convert the views of input and output to flat non-const chars,
    // so that there's no templating when we call doDeviceTransfers.
    auto inputRef  = charArrayRefFromArray(input.data(), input.size());
    auto outputRef = charArrayRefFromArray(output.data(), output.size());

    doDeviceTransfers(*this->gpuInfo_, inputRef, outputRef);
    this->compareVectors(input, output);
}

TYPED_TEST_CASE(HostAllocatorTest, TestTypes);

// Note that in GoogleTest typed tests, the use of TestFixture:: and
// this-> is sometimes required to get access to things in the fixture
// class (or its base classes).

// Note also that aspects of this code can be tested even when a GPU
// device is not available.

TYPED_TEST(HostAllocatorTest, EmptyMemoryAlwaysWorks)
{
    typename TestFixture::VectorType v;
}

TYPED_TEST(HostAllocatorTest, TransfersUsingDefaultHostAllocatorWork)
{
    typename TestFixture::VectorType input = {{1, 2, 3}}, output;
    output.resize(input.size());

    this->runTest(input, output);
}

TYPED_TEST(HostAllocatorTest, TransfersUsingNormalCpuHostAllocatorWork)
{
    // Make an allocator with a 'normal CPU' allocation policy. This
    // might be slower than another policy, but still works.
    using AllocatorType       = typename TestFixture::AllocatorType;
    using AllocatorPolicyType = typename AllocatorType::allocation_policy;
    AllocatorPolicyType              policy(AllocatorPolicyType::Impl::AllocateAligned);
    AllocatorType                    allocator(policy);

    typename TestFixture::VectorType input(allocator);
    this->fillInput(&input);
    typename TestFixture::VectorType output(allocator);
    output.resize(input.size());

    this->runTest(input, output);
}

TYPED_TEST(HostAllocatorTest, TransfersUsingGpuHostAllocatorWork)
{
    // Make an allocator with a 'for GPU' allocation policy. This
    // should be more efficient, but we can't test that.
    using AllocatorType       = typename TestFixture::AllocatorType;
    using AllocatorPolicyType = typename AllocatorType::allocation_policy;
    AllocatorPolicyType              policy(AllocatorPolicyType::Impl::AllocateForGpu);
    AllocatorType                    allocator(policy);

    typename TestFixture::VectorType input(allocator);
    this->fillInput(&input);
    typename TestFixture::VectorType output(allocator);
    output.resize(input.size());

    this->runTest(input, output);
}

TYPED_TEST(HostAllocatorTest, StatefulAllocatorUsesMemory)
{
    // The HostAllocator has state, so a container using it will be
    // larger than a normal vector, whose default allocator is
    // stateless.
    EXPECT_LT(sizeof(std::vector<typename TestFixture::ValueType>),
              sizeof(typename TestFixture::VectorType));
}

} // namespace
} // namespace
