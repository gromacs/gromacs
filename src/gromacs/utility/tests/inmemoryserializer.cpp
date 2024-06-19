/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests for gmx::InMemorySerializer and InMemoryDeserializer.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_utility
 */

#include "gmxpre.h"

#include "gromacs/utility/inmemoryserializer.h"

#include <cstdint>

#include <string>
#include <type_traits>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace test
{
namespace
{
union IntAndFloat32
{
    std::int32_t int32Value_;
    float        floatValue_;
};

union IntAndFloat64
{
    std::int64_t int64Value_;
    double       doubleValue_;
};

//! Constants used for testing endian swap operations
//! \{
constexpr std::int16_t c_int16Value        = static_cast<std::int16_t>(0x7A2B);
constexpr std::int16_t c_int16ValueSwapped = static_cast<std::int16_t>(0x2B7A);
constexpr std::int32_t c_int32Value        = static_cast<std::int32_t>(0x78ABCDEF);
constexpr std::int32_t c_int32ValueSwapped = static_cast<std::int32_t>(0xEFCDAB78);
constexpr std::int64_t c_int64Value        = static_cast<std::int64_t>(0x78ABCDEF12345678);
constexpr std::int64_t c_int64ValueSwapped = static_cast<std::int64_t>(0x78563412EFCDAB78);

constexpr const IntAndFloat32 c_intAndFloat32{ c_int32Value };
constexpr const IntAndFloat64 c_intAndFloat64{ c_int64Value };

constexpr const IntAndFloat32 c_intAndFloat32Swapped{ c_int32ValueSwapped };
constexpr const IntAndFloat64 c_intAndFloat64Swapped{ c_int64ValueSwapped };
//! \}

//! Return the integer used for testing, depending on the size of int.
constexpr int integerSizeDependentTestingValue()
{
    return sizeof(int) == 4 ? c_int32Value : sizeof(int) == 8 ? c_int64Value : c_int16Value;
}

//! Return the endianess-swapped integer used for testing, depending on the size of int.
constexpr int integerSizeDependentTestingValueEndianessSwapped()
{
    return sizeof(int) == 4   ? c_int32ValueSwapped
           : sizeof(int) == 8 ? c_int64ValueSwapped
                              : c_int16ValueSwapped;
}

class InMemorySerializerTest : public ::testing::Test
{
public:
    struct SerializerValues
    {
        bool           boolValue_;
        unsigned char  unsignedCharValue_;
        char           charValue_;
        unsigned short unsignedShortValue_;
        std::int32_t   int32Value_;
        float          floatValue_;
        std::int64_t   int64Value_;
        double         doubleValue_;
        int            intValue_;
        real           realValue_;
    };

    static void serialize(ISerializer* serializer, SerializerValues* values)
    {
        EXPECT_FALSE(serializer->reading());
        doValues(serializer, values);
    }

    static SerializerValues deserialize(ISerializer* serializer)
    {
        EXPECT_TRUE(serializer->reading());
        SerializerValues result;
        doValues(serializer, &result);
        return result;
    }

    static void checkSerializerValuesforEquality(const SerializerValues& lhs, const SerializerValues& rhs)
    {
        EXPECT_EQ(lhs.boolValue_, rhs.boolValue_);
        EXPECT_EQ(lhs.unsignedCharValue_, rhs.unsignedCharValue_);
        EXPECT_EQ(lhs.charValue_, rhs.charValue_);
        EXPECT_EQ(lhs.intValue_, rhs.intValue_);
        EXPECT_EQ(lhs.int32Value_, rhs.int32Value_);
        EXPECT_EQ(lhs.int64Value_, rhs.int64Value_);
        EXPECT_EQ(lhs.unsignedShortValue_, rhs.unsignedShortValue_);
        EXPECT_EQ(lhs.realValue_, rhs.realValue_);
        EXPECT_EQ(lhs.floatValue_, rhs.floatValue_);
        EXPECT_EQ(lhs.doubleValue_, rhs.doubleValue_);
    }

private:
    static void doValues(ISerializer* serializer, SerializerValues* values)
    {
        serializer->doBool(&values->boolValue_);
        serializer->doUChar(&values->unsignedCharValue_);
        serializer->doChar(&values->charValue_);
        serializer->doInt(&values->intValue_);
        serializer->doInt32(&values->int32Value_);
        serializer->doInt64(&values->int64Value_);
        serializer->doUShort(&values->unsignedShortValue_);
        serializer->doReal(&values->realValue_);
        serializer->doFloat(&values->floatValue_);
        serializer->doDouble(&values->doubleValue_);
    }

protected:
    SerializerValues defaultValues_ = { true,
                                        0x78,
                                        0x78,
                                        static_cast<unsigned short>(c_int16Value),
                                        c_int32Value,
                                        c_intAndFloat32.floatValue_,
                                        c_int64Value,
                                        c_intAndFloat64.doubleValue_,
                                        integerSizeDependentTestingValue(),
                                        std::is_same_v<real, double>
                                                ? static_cast<real>(c_intAndFloat64.doubleValue_)
                                                : static_cast<real>(c_intAndFloat32.floatValue_) };

    SerializerValues endianessSwappedValues_ = {
        true,
        0x78,
        0x78,
        static_cast<unsigned short>(c_int16ValueSwapped),
        c_int32ValueSwapped,
        c_intAndFloat32Swapped.floatValue_,
        c_int64ValueSwapped,
        c_intAndFloat64Swapped.doubleValue_,
        integerSizeDependentTestingValueEndianessSwapped(),
        std::is_same_v<real, float> ? static_cast<real>(c_intAndFloat32Swapped.floatValue_)
                                    : static_cast<real>(c_intAndFloat64Swapped.doubleValue_)
    };
};

TEST_F(InMemorySerializerTest, Roundtrip)
{
    InMemorySerializer serializer;
    SerializerValues   values = defaultValues_;
    serialize(&serializer, &values);

    auto buffer = serializer.finishAndGetBuffer();

    InMemoryDeserializer deserializer(buffer, std::is_same_v<real, double>);

    SerializerValues deserialisedValues = deserialize(&deserializer);

    checkSerializerValuesforEquality(values, deserialisedValues);
}

TEST_F(InMemorySerializerTest, RoundtripWithEndianessSwap)
{
    InMemorySerializer serializerWithSwap(EndianSwapBehavior::Swap);
    SerializerValues   values = defaultValues_;
    serialize(&serializerWithSwap, &values);

    auto buffer = serializerWithSwap.finishAndGetBuffer();

    InMemoryDeserializer deserializerWithSwap(
            buffer, std::is_same_v<real, double>, EndianSwapBehavior::Swap);

    SerializerValues deserialisedValues = deserialize(&deserializerWithSwap);

    checkSerializerValuesforEquality(values, deserialisedValues);
}

TEST_F(InMemorySerializerTest, SerializerExplicitEndianessSwap)
{
    InMemorySerializer serializerWithSwap(EndianSwapBehavior::Swap);
    SerializerValues   values = defaultValues_;
    serialize(&serializerWithSwap, &values);

    auto buffer = serializerWithSwap.finishAndGetBuffer();

    InMemoryDeserializer deserializerWithOutSwap(buffer, std::is_same_v<real, double>);

    SerializerValues deserialisedValues = deserialize(&deserializerWithOutSwap);
    checkSerializerValuesforEquality(endianessSwappedValues_, deserialisedValues);
}

TEST_F(InMemorySerializerTest, DeserializerExplicitEndianessSwap)
{
    InMemorySerializer serializer;
    SerializerValues   values = defaultValues_;
    serialize(&serializer, &values);

    auto buffer = serializer.finishAndGetBuffer();

    InMemoryDeserializer deserializerWithSwap(
            buffer, std::is_same_v<real, double>, EndianSwapBehavior::Swap);

    SerializerValues deserialisedValues = deserialize(&deserializerWithSwap);
    checkSerializerValuesforEquality(endianessSwappedValues_, deserialisedValues);
}

TEST_F(InMemorySerializerTest, SizeIsCorrect)
{
    InMemorySerializer serializer;
    // These types all have well-defined widths in bytes,
    // which we can test below.
    serializer.doBool(&defaultValues_.boolValue_);     // 1 bytes
    serializer.doChar(&defaultValues_.charValue_);     // 1 bytes
    serializer.doInt32(&defaultValues_.int32Value_);   // 4 bytes
    serializer.doInt64(&defaultValues_.int64Value_);   // 8 bytes
    serializer.doFloat(&defaultValues_.floatValue_);   // 4 bytes
    serializer.doDouble(&defaultValues_.doubleValue_); // 8 bytes
    std::vector<char> charBuffer = { 'a', 'b', 'c' };
    serializer.doCharArray(charBuffer.data(), charBuffer.size()); // 3 bytes
    serializer.doOpaque(charBuffer.data(), charBuffer.size());    // 3 bytes
    std::vector<int32_t> int32Buffer = { 0x1BCDEF78, 0x654321FE };
    serializer.doInt32Array(int32Buffer.data(), int32Buffer.size()); // 8 bytes
    std::vector<int64_t> int64Buffer = { 0x1BCDEF78654321FE, 0x3726ABFEAB34716C };
    serializer.doInt64Array(int64Buffer.data(), int64Buffer.size()); // 16 bytes
    auto buffer = serializer.finishAndGetBuffer();
    EXPECT_EQ(buffer.size(), 56);
}

} // namespace
} // namespace test
} // namespace gmx
