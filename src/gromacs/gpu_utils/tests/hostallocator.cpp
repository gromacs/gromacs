/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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

#include "gromacs/math/tests/testarrayrefs.h"

#include "devicetransfers.h"
#include "gputest.h"

namespace gmx
{

namespace test
{

/*! \internal \brief Typed test fixture for infrastructure for
 * host-side memory used for GPU transfers. */
template<typename T>
class HostMemoryTest : public test::GpuTest
{
public:
    //! Convenience type
    using ValueType = T;
};

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
template<typename T>
ArrayRef<char> charArrayRefFromArray(T* data, size_t size)
{
    // Make a type like T, but without its possible const qualifier.
    using NonConstT = std::remove_const_t<T>;
    return arrayRefFromArray<char>(reinterpret_cast<char*>(const_cast<NonConstT*>(data)),
                                   size * sizeof(T));
}

//! Does a device transfer of \c input to the device in \c gpuInfo, and back to \c output.
template<typename T>
void runTest(const gmx_gpu_info_t& gpuInfo, ArrayRef<T> input, ArrayRef<T> output)
{
    // Convert the views of input and output to flat non-const chars,
    // so that there's no templating when we call doDeviceTransfers.
    auto inputRef  = charArrayRefFromArray(input.data(), input.size());
    auto outputRef = charArrayRefFromArray(output.data(), output.size());

    ASSERT_EQ(inputRef.size(), outputRef.size());
    doDeviceTransfers(gpuInfo, inputRef, outputRef);
    compareViews(input, output);
}

struct MoveOnly
{
    MoveOnly(real x = 0) : x(x) {}
    MoveOnly(const MoveOnly&) = delete;
    MoveOnly(MoveOnly&&)      = default;
    MoveOnly& operator=(const MoveOnly&) = delete;
    MoveOnly& operator=(MoveOnly&&) = default;
    bool      operator==(const MoveOnly& o) const { return x == o.x; }
    real      operator*=(int scaleFactor) { return x *= scaleFactor; }
    real      x;
};

} // namespace test

namespace detail
{

template<>
struct PaddingTraits<test::MoveOnly>
{
    using SimdBaseType                          = real;
    static constexpr int maxSimdWidthOfBaseType = GMX_REAL_MAX_SIMD_WIDTH;
};

} // namespace detail

namespace test
{

//! The types used in testing of all operations.
typedef ::testing::Types<int32_t, real, RVec, test::MoveOnly> TestTypes;

//! Typed test fixture
template<typename T>
struct HostAllocatorTest : HostMemoryTest<T>
{
    using VectorType = PaddedHostVector<T>; //!< PaddedHostVector of type tested
};
TYPED_TEST_CASE(HostAllocatorTest, TestTypes);

//! Typed test fixture (no mem/gpu initializtion - much faster)
template<typename T>
struct HostAllocatorTestNoMem : ::testing::Test
{
    using VectorType = PaddedHostVector<T>; //!< PaddedHostVector of type tested
};
TYPED_TEST_CASE(HostAllocatorTestNoMem, TestTypes);

//! Typed test fixture for tests requiring a copyable type
template<typename T>
struct HostAllocatorTestNoMemCopyable : HostAllocatorTestNoMem<T>
{
};
//! The types used in testing minus move only types
using TestTypesCopyable = ::testing::Types<int32_t, real, RVec>;

TYPED_TEST_CASE(HostAllocatorTestNoMemCopyable, TestTypesCopyable);

//! Typed test fixture for tests requiring a copyable type
template<typename T>
using HostAllocatorTestCopyable = HostAllocatorTest<T>;
TYPED_TEST_CASE(HostAllocatorTestCopyable, TestTypesCopyable);

// Note that in GoogleTest typed tests, the use of TestFixture:: and
// this-> is sometimes required to get access to things in the fixture
// class (or its base classes).

// Note also that aspects of this code can be tested even when a GPU
// device is not available.

TYPED_TEST(HostAllocatorTest, EmptyMemoryAlwaysWorks)
{
    typename TestFixture::VectorType v;
}

TYPED_TEST(HostAllocatorTestCopyable, VectorsWithDefaultHostAllocatorAlwaysWorks)
{
    typename TestFixture::VectorType input(3), output;
    output.resizeWithPadding(input.size());
}

// Several tests actually do CUDA transfers. This is not necessary
// because the state of page alignment or pinning is not currently
// relevant to the success of a CUDA transfer. CUDA checks happen only
// during cudaHostRegister and cudaHostUnregister. Such tests are of
// value only when this behaviour changes, if ever.

TYPED_TEST(HostAllocatorTestCopyable, TransfersWithoutPinningWork)
{
    typename TestFixture::VectorType input;
    fillInput(&input, 1);
    typename TestFixture::VectorType output;
    output.resizeWithPadding(input.size());

    runTest(*this->gpuInfo_, makeArrayRef(input), makeArrayRef(output));
}

TYPED_TEST(HostAllocatorTestCopyable, FillInputAlsoWorksAfterCallingReserve)
{
    typename TestFixture::VectorType input;
    input.reserveWithPadding(3);
    fillInput(&input, 1);
}

TYPED_TEST(HostAllocatorTestNoMem, CreateVector)
{
    typename TestFixture::VectorType input1;
    EXPECT_FALSE(input1.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
    typename TestFixture::VectorType input2({ PinningPolicy::PinnedIfSupported });
    EXPECT_TRUE(input2.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
}

TYPED_TEST(HostAllocatorTestNoMem, MoveAssignment)
{
    typename TestFixture::VectorType input1({ PinningPolicy::PinnedIfSupported });
    input1 = typename TestFixture::VectorType();
    EXPECT_FALSE(input1.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);

    typename TestFixture::VectorType input2;
    input2 = typename TestFixture::VectorType({ PinningPolicy::PinnedIfSupported });
    EXPECT_TRUE(input2.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
}

TYPED_TEST(HostAllocatorTestNoMem, MoveConstruction)
{
    typename TestFixture::VectorType input1;
    typename TestFixture::VectorType input2(std::move(input1));
    EXPECT_FALSE(input2.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);

    typename TestFixture::VectorType input3({ PinningPolicy::PinnedIfSupported });
    typename TestFixture::VectorType input4(std::move(input3));
    EXPECT_TRUE(input4.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
}

TYPED_TEST(HostAllocatorTestNoMemCopyable, CopyAssignment)
{
    typename TestFixture::VectorType input1;
    typename TestFixture::VectorType input2({ PinningPolicy::PinnedIfSupported });
    input1 = input2;
    EXPECT_FALSE(input1.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
    EXPECT_TRUE(input2.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
    input2 = input1;
    EXPECT_FALSE(input1.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
    EXPECT_TRUE(input2.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
}

TYPED_TEST(HostAllocatorTestNoMemCopyable, CopyConstruction)
{
    typename TestFixture::VectorType input1;
    typename TestFixture::VectorType input2(input1); //NOLINT(performance-unnecessary-copy-initialization)
    EXPECT_FALSE(input2.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);

    typename TestFixture::VectorType input3({ PinningPolicy::PinnedIfSupported });
    typename TestFixture::VectorType input4(input3); //NOLINT(performance-unnecessary-copy-initialization)
    EXPECT_FALSE(input4.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
}

TYPED_TEST(HostAllocatorTestNoMem, Swap)
{
    typename TestFixture::VectorType input1;
    typename TestFixture::VectorType input2({ PinningPolicy::PinnedIfSupported });
    std::swap(input1, input2);
    EXPECT_TRUE(input1.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
    EXPECT_FALSE(input2.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
    std::swap(input2, input1);
    EXPECT_FALSE(input1.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
    EXPECT_TRUE(input2.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
}

TYPED_TEST(HostAllocatorTestNoMem, Comparison)
{
    using AllocatorType = typename TestFixture::VectorType::allocator_type;
    EXPECT_EQ(AllocatorType{}, AllocatorType{});
    // Should be false for different pinning policy
    EXPECT_NE(AllocatorType{}, AllocatorType{ PinningPolicy::PinnedIfSupported });
}

#if GMX_GPU == GMX_GPU_CUDA

// Policy suitable for pinning is only supported for a CUDA build

TYPED_TEST(HostAllocatorTestCopyable, TransfersWithPinningWorkWithCuda)
{
    if (!this->haveCompatibleGpus())
    {
        return;
    }

    typename TestFixture::VectorType input;
    changePinningPolicy(&input, PinningPolicy::PinnedIfSupported);
    fillInput(&input, 1);
    typename TestFixture::VectorType output;
    changePinningPolicy(&output, PinningPolicy::PinnedIfSupported);
    output.resizeWithPadding(input.size());

    runTest(*this->gpuInfo_, makeArrayRef(input), makeArrayRef(output));
}

//! Helper function for wrapping a call to isHostMemoryPinned.
template<typename VectorType>
bool isPinned(const VectorType& v)
{
    void* data = const_cast<void*>(static_cast<const void*>(v.data()));
    return isHostMemoryPinned(data);
}

TYPED_TEST(HostAllocatorTestCopyable, ManualPinningOperationsWorkWithCuda)
{
    if (!this->haveCompatibleGpus())
    {
        return;
    }

    typename TestFixture::VectorType input;
    changePinningPolicy(&input, PinningPolicy::PinnedIfSupported);
    EXPECT_TRUE(input.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
    EXPECT_EQ(0, input.size());
    EXPECT_EQ(0, input.paddedSize());
    EXPECT_TRUE(input.empty());
    EXPECT_FALSE(isPinned(input));

    // Fill some contents, which will be pinned because of the policy.
    fillInput(&input, 1);
    EXPECT_TRUE(isPinned(input));

    // Switching policy to CannotBePinned must unpin the buffer (via
    // realloc and copy).
    auto oldInputData = input.data();
    changePinningPolicy(&input, PinningPolicy::CannotBePinned);
    EXPECT_FALSE(isPinned(input));
    // These cannot be equal as both had to be allocated at the same
    // time for the contents to be able to be copied.
    EXPECT_NE(oldInputData, input.data());

    // Switching policy to PinnedIfSupported must pin the buffer (via
    // realloc and copy).
    oldInputData = input.data();
    changePinningPolicy(&input, PinningPolicy::PinnedIfSupported);
    EXPECT_TRUE(isPinned(input));
    // These cannot be equal as both had to be allocated at the same
    // time for the contents to be able to be copied.
    EXPECT_NE(oldInputData, input.data());
}

#endif

TYPED_TEST(HostAllocatorTest, StatefulAllocatorUsesMemory)
{
    // The HostAllocator has state, so a container using it will be
    // larger than a normal vector, whose default allocator is
    // stateless.
    EXPECT_LT(sizeof(std::vector<typename TestFixture::VectorType::value_type>),
              sizeof(typename TestFixture::VectorType));
}

TEST(HostAllocatorUntypedTest, Comparison)
{
    // Should always be true for the same policy, indpendent of value_type
    EXPECT_EQ(HostAllocator<float>{}, HostAllocator<double>{});
}

//! Declare allocator types to test.
using AllocatorTypesToTest =
        ::testing::Types<HostAllocator<real>, HostAllocator<int32_t>, HostAllocator<RVec>, HostAllocator<MoveOnly>>;

TYPED_TEST_CASE(AllocatorTest, AllocatorTypesToTest);

} // namespace test
} // namespace gmx

// Includes tests common to all allocation policies.
#include "gromacs/utility/tests/alignedallocator_impl.h"
