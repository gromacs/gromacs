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

#include "config.h"

#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "devicetransfers.h"
#include "gputest.h"

namespace gmx
{

namespace
{

/*! \internal \brief Typed test fixture for infrastructure for
 * host-side memory used for GPU transfers. */
template <typename T>
class HostMemoryTest : public test::GpuTest
{
    public:
        //! Convenience type
        using ValueType = T;
        //! Convenience type
        using ViewType = ArrayRef<ValueType>;
        //! Convenience type
        using ConstViewType = ArrayRef<const ValueType>;
        //! Prepare contents of a VectorType.
        template <typename VectorType>
        void fillInput(VectorType *input) const;
        //! Compares input and output vectors.
        void compareVectors(ConstViewType input,
                            ConstViewType output) const;
        //! Do some transfers and test the results.
        void runTest(ConstViewType input, ViewType output) const;
};

// Already documented
template <typename T> template <typename VectorType>
void HostMemoryTest<T>::fillInput(VectorType *input) const
{
    input->resize(3);
    (*input)[0] = 1;
    (*input)[1] = 2;
    (*input)[2] = 3;
}

//! Initialization specialization for RVec
template <> template <typename VectorType>
void HostMemoryTest<RVec>::fillInput(VectorType *input) const
{
    input->reserve(3);
    input->resize(3);
    (*input)[0] = {1, 2, 3};
    (*input)[1] = {4, 5, 6};
    (*input)[2] = {7, 8, 9};
}

// Already documented
template <typename T>
void HostMemoryTest<T>::compareVectors(ConstViewType input,
                                       ConstViewType output) const
{
    for (size_t i = 0; i != input.size(); ++i)
    {
        EXPECT_EQ(input[i], output[i]) << "for index " << i;
    }
}

//! Comparison specialization for RVec
template <>
void HostMemoryTest<RVec>::compareVectors(ConstViewType input,
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
void HostMemoryTest<T>::runTest(ConstViewType input, ViewType output) const
{
    // Convert the views of input and output to flat non-const chars,
    // so that there's no templating when we call doDeviceTransfers.
    auto inputRef  = charArrayRefFromArray(input.data(), input.size());
    auto outputRef = charArrayRefFromArray(output.data(), output.size());

    doDeviceTransfers(*this->gpuInfo_, inputRef, outputRef);
    this->compareVectors(input, output);
}

//! The types used in testing.
typedef ::testing::Types<int, real, RVec> TestTypes;

//! Typed test fixture
template <typename T>
class HostAllocatorTest : public HostMemoryTest<T>
{
    public:
        //! Convenience type
        using ValueType = T;
        //! Convenience type
        using AllocatorType = HostAllocator<T>;
        //! Convenience type
        using VectorType = std::vector<ValueType, AllocatorType>;
};

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

TYPED_TEST(HostAllocatorTest, VectorsWithDefaultHostAllocatorAlwaysWorks)
{
    typename TestFixture::VectorType input = {{1, 2, 3}}, output;
    output.resize(input.size());
}

// Several tests actually do CUDA transfers. This is not necessary
// because the state of page alignment or pinning is not currently
// relevant to the success of a CUDA transfer. CUDA checks happen only
// during cudaHostRegister and cudaHostUnregister. Such tests are of
// value only when this behaviour changes, if ever.

TYPED_TEST(HostAllocatorTest, TransfersWithoutPinningWork)
{
    typename TestFixture::VectorType input;
    this->fillInput(&input);
    typename TestFixture::VectorType output;
    output.resize(input.size());

    this->runTest(input, output);
}

TYPED_TEST(HostAllocatorTest, FillInputAlsoWorksAfterCallingReserve)
{
    typename TestFixture::VectorType input;
    input.reserve(3);
    this->fillInput(&input);
}

#if GMX_GPU == GMX_GPU_CUDA

// Policy suitable for pinning is only supported for a CUDA build

TYPED_TEST(HostAllocatorTest, TransfersWithPinningWorkWithCuda)
{
    typename TestFixture::VectorType input;
    changePinningPolicy(&input, PinningPolicy::CanBePinned);
    this->fillInput(&input);
    typename TestFixture::VectorType output;
    changePinningPolicy(&output, PinningPolicy::CanBePinned);
    output.resize(input.size());

    this->runTest(input, output);
}

//! Helper function for wrapping a call to isHostMemoryPinned.
template <typename VectorType>
bool isPinned(const VectorType &v)
{
    void *data = const_cast<void *>(static_cast<const void *>(v.data()));
    return isHostMemoryPinned(data);
}

TYPED_TEST(HostAllocatorTest, ManualPinningOperationsWorkWithCuda)
{
    typename TestFixture::VectorType input;
    changePinningPolicy(&input, PinningPolicy::CanBePinned);
    EXPECT_FALSE(isPinned(input));

    // Unpin before allocation is fine, but does nothing.
    input.get_allocator().getPolicy().unpin();
    EXPECT_FALSE(isPinned(input));

    // Pin with no contents is fine, but does nothing.
    input.get_allocator().getPolicy().pin();
    EXPECT_FALSE(isPinned(input));

    // Fill some contents, which will be pinned because of the policy.
    this->fillInput(&input);
    EXPECT_TRUE(isPinned(input));

    // Unpin after pin is fine.
    input.get_allocator().getPolicy().unpin();
    EXPECT_FALSE(isPinned(input));

    // Repeated unpin should be a no-op.
    input.get_allocator().getPolicy().unpin();

    // Pin after unpin is fine.
    input.get_allocator().getPolicy().pin();
    EXPECT_TRUE(isPinned(input));

    // Repeated pin should be a no-op, and still pinned.
    input.get_allocator().getPolicy().pin();
    EXPECT_TRUE(isPinned(input));

    // Switching policy to CannotBePinned must unpin the buffer (via
    // realloc and copy).
    auto oldInputData = input.data();
    changePinningPolicy(&input, PinningPolicy::CannotBePinned);
    EXPECT_FALSE(isPinned(input));
    // These cannot be equal as both had to be allocated at the same
    // time for the contents to be able to be copied.
    EXPECT_NE(oldInputData, input.data());

    // Switching policy to CanBePinned must pin the buffer (via
    // realloc and copy).
    oldInputData = input.data();
    changePinningPolicy(&input, PinningPolicy::CanBePinned);
    EXPECT_TRUE(isPinned(input));
    // These cannot be equal as both had to be allocated at the same
    // time for the contents to be able to be copied.
    EXPECT_NE(oldInputData, input.data());
}

#else

TYPED_TEST(HostAllocatorTest, ChangingPinningPolicyRequiresCuda)
{
    typename TestFixture::VectorType input;
    EXPECT_DEATH(changePinningPolicy(&input, PinningPolicy::CanBePinned),
                 ".*A suitable build of GROMACS.* is required.*");
}

TYPED_TEST(HostAllocatorTest, ManualPinningOperationsWorkEvenWithoutCuda)
{
    typename TestFixture::VectorType input;

    // Since the buffer can't be pinned and isn't pinned, and the
    // calling code can't be unhappy about this, these are OK.
    input.get_allocator().getPolicy().pin();
    input.get_allocator().getPolicy().unpin();
}

#endif

TYPED_TEST(HostAllocatorTest, StatefulAllocatorUsesMemory)
{
    // The HostAllocator has state, so a container using it will be
    // larger than a normal vector, whose default allocator is
    // stateless.
    EXPECT_LT(sizeof(std::vector<typename TestFixture::ValueType>),
              sizeof(typename TestFixture::VectorType));
}

//! Declare allocator types to test.
using AllocatorTypesToTest = ::testing::Types<HostAllocator<real>,
                                              HostAllocator<int>,
                                              HostAllocator<RVec>
                                              >;

TYPED_TEST_CASE(AllocatorTest, AllocatorTypesToTest);

} // namespace
} // namespace

// Includes tests common to all allocation policies.
#include "gromacs/utility/tests/alignedallocator-impl.h"
