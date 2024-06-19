/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * \brief
 * Tests for GPU host allocator.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/hostallocator.h"

#include "config.h"

#include <cstddef>
#include <cstdint>

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/math/tests/testarrayrefs.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"

#include "devicetransfers.h"

namespace gmx
{
namespace detail
{
template<typename T>
struct PaddingTraits;
} // namespace detail

namespace test
{

/*! \internal \brief Typed test fixture for infrastructure for
 * host-side memory used for GPU transfers. */
template<typename T>
class HostMemoryTest : public ::testing::Test
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
void runTest(const DeviceContext& deviceContext, ArrayRef<T> input, ArrayRef<T> output)
{
    // Convert the views of input and output to flat non-const chars,
    // so that there's no templating when we call doDeviceTransfers.
    auto inputRef  = charArrayRefFromArray(input.data(), input.size());
    auto outputRef = charArrayRefFromArray(output.data(), output.size());

    ASSERT_EQ(inputRef.size(), outputRef.size());
    doDeviceTransfers(deviceContext, inputRef, outputRef);
    compareViews(input, output);
}

struct MoveOnly
{
    MoveOnly(real x = 0) : x_(x) {}
    MoveOnly(const MoveOnly&) = delete;
    MoveOnly(MoveOnly&&)      = default;
    MoveOnly& operator=(const MoveOnly&) = delete;
    MoveOnly& operator=(MoveOnly&&) = default;
    bool      operator==(const MoveOnly& o) const { return x_ == o.x_; }
    real      operator*=(int scaleFactor) { return x_ *= scaleFactor; }
    real      x_;
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
TYPED_TEST_SUITE(HostAllocatorTest, TestTypes);

//! Typed test fixture (no mem/gpu initializtion - much faster)
template<typename T>
struct HostAllocatorTestNoMem : ::testing::Test
{
    using VectorType = PaddedHostVector<T>; //!< PaddedHostVector of type tested
};
TYPED_TEST_SUITE(HostAllocatorTestNoMem, TestTypes);

//! Typed test fixture for tests requiring a copyable type
template<typename T>
struct HostAllocatorTestNoMemCopyable : HostAllocatorTestNoMem<T>
{
};
//! The types used in testing minus move only types
using TestTypesCopyable = ::testing::Types<int32_t, real, RVec>;

TYPED_TEST_SUITE(HostAllocatorTestNoMemCopyable, TestTypesCopyable);

//! Typed test fixture for tests requiring a copyable type
template<typename T>
using HostAllocatorTestCopyable = HostAllocatorTest<T>;
TYPED_TEST_SUITE(HostAllocatorTestCopyable, TestTypesCopyable);

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
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        testDevice->deviceContext().activate();
        typename TestFixture::VectorType input;
        resizeAndFillInput(&input, 3, 1);
        typename TestFixture::VectorType output;
        output.resizeWithPadding(input.size());

        runTest(testDevice->deviceContext(), makeArrayRef(input), makeArrayRef(output));
    }
}

TYPED_TEST(HostAllocatorTestCopyable, FillInputAlsoWorksAfterCallingReserve)
{
    typename TestFixture::VectorType input;
    input.reserveWithPadding(3);
    resizeAndFillInput(&input, 3, 1);
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

template<typename T>
struct HoldsHostVector
{
    HoldsHostVector(HostAllocationPolicy p) : v(1, p) {}
    HostVector<T> v;
};

TYPED_TEST(HostAllocatorTestNoMemCopyable, CopyConstructionOfStructHoldingAHostVectorDoesNotCopyTheAllocator)
{
    using Holder = HoldsHostVector<typename TestFixture::VectorType::value_type>;
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        testDevice->deviceContext().activate();
        SCOPED_TRACE("By default allocator does not propagate");
        {
            Holder c{ { PinningPolicy::PinnedIfSupported } };
            EXPECT_EQ(c.v.get_allocator().pinningPolicy(), PinningPolicy::PinnedIfSupported);
            std::vector<Holder> v(2, c);
            const auto          defaultPolicy = PinningPolicy::CannotBePinned;
            EXPECT_EQ(v[0].v.get_allocator().pinningPolicy(), defaultPolicy);
            EXPECT_EQ(v[1].v.get_allocator().pinningPolicy(), defaultPolicy);
        }
        SCOPED_TRACE("Allocator can propagate");
        {
            Holder c{ { PinningPolicy::PinnedIfSupported, true } };
            EXPECT_EQ(c.v.get_allocator().pinningPolicy(), PinningPolicy::PinnedIfSupported);
            std::vector<Holder> v(2, c);
            EXPECT_EQ(v[0].v.get_allocator().pinningPolicy(), c.v.get_allocator().pinningPolicy());
            EXPECT_EQ(v[1].v.get_allocator().pinningPolicy(), c.v.get_allocator().pinningPolicy());
        }
    }
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

#if GMX_GPU_CUDA || GMX_GPU_SYCL || GMX_GPU_HIP

// Policy suitable for pinning is supported for a CUDA, HIP or SYCL build

TYPED_TEST(HostAllocatorTestCopyable, TransfersWithPinningWorkWithDevice)
{
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        testDevice->deviceContext().activate();
        typename TestFixture::VectorType input;
        changePinningPolicy(&input, PinningPolicy::PinnedIfSupported);
        resizeAndFillInput(&input, 3, 1);
        typename TestFixture::VectorType output;
        changePinningPolicy(&output, PinningPolicy::PinnedIfSupported);
        output.resizeWithPadding(input.size());

        runTest(testDevice->deviceContext(), makeArrayRef(input), makeArrayRef(output));
    }
}

#endif

#if GMX_GPU_CUDA || GMX_GPU_HIP

// While we can allocate pinned memory with SYCL, we don't support isHostMemoryPinned yet. See #4522

//! Helper function for wrapping a call to isHostMemoryPinned.
template<typename VectorType>
bool isPinned(const VectorType& v)
{
    void* data = const_cast<void*>(static_cast<const void*>(v.data()));
    return isHostMemoryPinned(data);
}

TYPED_TEST(HostAllocatorTestCopyable, ManualPinningOperationsWork)
{
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        testDevice->deviceContext().activate();
        typename TestFixture::VectorType input;
        changePinningPolicy(&input, PinningPolicy::PinnedIfSupported);
        EXPECT_TRUE(input.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported);
        EXPECT_TRUE(input.empty());
        resizeAndFillInput(&input, 3, 1);
        // realloc and copy).
        auto* oldInputData = input.data();
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

TYPED_TEST_SUITE(AllocatorTest, AllocatorTypesToTest);

} // namespace test
} // namespace gmx

// Includes tests common to all allocation policies.
#include "gromacs/utility/tests/alignedallocator_impl.h"
