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
 * Tests for gmx::XdrSerializer.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/xdr_serializer.h"

#include <cstdint>
#include <cstdio>

#include <filesystem>
#include <string>
#include <type_traits>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! Constants useful for observing endian-swap bugs
//! \{
constexpr std::int16_t c_int16Value = static_cast<std::int16_t>(0x7A2B);
constexpr std::int32_t c_int32Value = static_cast<std::int32_t>(0x78ABCDEF);
constexpr std::int64_t c_int64Value = static_cast<std::int64_t>(0x78ABCDEF12345678);
//! \}

struct XdrSerializerTestParameters
{
    bool writeAsDouble;
};

/*! \brief Test XDR serialization works
 *
 * We want to be sure that XDR serialization works when reading and
 * writing reals and rvecs from both double- and mixed-precision
 * builds of GROMACS. The GROMACS XDR-file headers must write a value
 * that describes the precision of the GROMACS that wrote the file,
 * which provides clarity to a reader built in either precision as to
 * what sizes of values it should expect to find in the XDR. In
 * testing, we can direct either kind of GROMACS configuration to
 * write as if it was either kind of configuration, and to read back
 * natively.
 */
class XdrSerializerTest : public ::testing::TestWithParam<XdrSerializerTestParameters>
{
public:
    // These types all have well-defined widths in bytes AFTER XDR serialization,
    // which we can test below.
    struct SerializerValues
    {
        bool                       boolValue;
        unsigned char              unsignedCharValue;
        char                       charValue;
        unsigned short             unsignedShortValue;
        int                        intValue;
        std::int32_t               int32Value;
        std::int64_t               int64Value;
        float                      floatValue;
        double                     doubleValue;
        real                       realValue;
        ivec                       legacyIVecValue;
        IVec                       ivecValue;
        rvec                       legacyRVecValue;
        RVec                       rvecValue;
        std::vector<char>          charArray;
        std::vector<unsigned char> ucharArray;
        std::vector<std::byte>     byteArray;
        std::vector<int32_t>       int32Array;
        std::vector<int64_t>       int64Array;
        std::vector<RVec>          rvecArray;
        std::string                emptyStringValue;
        std::string                stringValue;
        void                       serialize(ISerializer* serializer)
        {
            // Serialized sizes are found in comments, according to
            // the precision of the GROMACS configuration that wrote
            // the file (applicable to real and rvec).
            serializer->doBool(&boolValue);            // 4 bytes serialized
            serializer->doUChar(&unsignedCharValue);   // 4 bytes serialized
            serializer->doChar(&charValue);            // 4 bytes serialized
            serializer->doUShort(&unsignedShortValue); // 4 bytes serialized
            serializer->doInt(&intValue);              // 4 bytes serialized
            serializer->doInt32(&int32Value);          // 4 bytes serialized
            serializer->doInt64(&int64Value);          // 8 bytes serialized
            serializer->doFloat(&floatValue);          // 4 bytes serialized
            serializer->doDouble(&doubleValue);        // 8 bytes serialized
            serializer->doReal(&realValue);            // 4 or 8 bytes serialized
            serializer->doIvec(&legacyIVecValue);      // 12 bytes serialized
            serializer->doIvec(&ivecValue);            // 12 bytes serialized
            serializer->doRvec(&legacyRVecValue);      // 12 or 24 bytes serialized
            serializer->doRvec(&rvecValue);            // 12 or 24 bytes serialized

            // Ensure for reading that arrays are resized appropriately
            int size = charArray.size();
            serializer->doInt(&size); // 4 bytes serialized
            charArray.resize(size);
            serializer->doCharArray(charArray.data(), charArray.size()); // 12 bytes serialized

            size = ucharArray.size();
            serializer->doInt(&size); // 4 bytes serialized
            ucharArray.resize(size);
            serializer->doUCharArray(ucharArray.data(), ucharArray.size()); // 12 bytes serialized

            size = byteArray.size();
            serializer->doInt(&size); // 4 bytes serialized
            byteArray.resize(size);
            serializer->doOpaque(reinterpret_cast<char*>(byteArray.data()),
                                 byteArray.size()); // 4 bytes serialized

            size = int32Array.size();
            serializer->doInt(&size); // 4 bytes serialized
            int32Array.resize(size);
            serializer->doInt32Array(int32Array.data(), int32Array.size()); // 8 bytes serialized

            size = int64Array.size();
            serializer->doInt(&size); // 4 bytes serialized
            int64Array.resize(size);
            serializer->doInt64Array(int64Array.data(), int64Array.size()); // 16 bytes serialized

            size = rvecArray.size();
            serializer->doInt(&size); // 4 bytes serialized
            if (GMX_DOUBLE && serializer->reading())
            {
                // Note that default-initialization of real[DIM] is
                // not zero initialization, which might cause overflow
                // when reading and changing precision in the
                // implementation of doRvec. This scenario only
                // happens in these unit tests, as normal GROMACS
                // always writes double-precision files when compiled
                // in double precision.
                rvecArray.resize(size, { 0, 0, 0 });
            }
            else
            {
                rvecArray.resize(size);
            }
            serializer->doRvecArray(rvecArray); // 24 or 48 bytes serialized

            // 8 bytes serialized, 4 for the size, 1 for the null char
            // and some from the implementation
            serializer->doString(&emptyStringValue);
            // 16 bytes serialized, 4 for the size, 8 for the chars
            // and some from the implementation
            serializer->doString(&stringValue);
        }
    };

    TestFileManager       fileManager_;
    std::filesystem::path filename_ = fileManager_.getTemporaryFilePath("data.bin");
};

TEST_P(XdrSerializerTest, Works)
{
    SerializerValues valuesToWrite;
    valuesToWrite.boolValue          = true;
    valuesToWrite.unsignedCharValue  = 0x78;
    valuesToWrite.charValue          = 0x78;
    valuesToWrite.unsignedShortValue = static_cast<unsigned short>(c_int16Value);
    valuesToWrite.intValue           = c_int32Value;
    valuesToWrite.int32Value         = c_int32Value;
    valuesToWrite.int64Value         = c_int64Value;
    valuesToWrite.floatValue         = 3.0e3f;
    valuesToWrite.doubleValue        = 3.0e3; // Should be in range for conversion to float
    valuesToWrite.realValue =
            (std::is_same_v<real, double> ? valuesToWrite.doubleValue : valuesToWrite.floatValue);
    valuesToWrite.legacyIVecValue[XX] = 1;
    valuesToWrite.legacyIVecValue[YY] = 2;
    valuesToWrite.legacyIVecValue[ZZ] = 3;
    valuesToWrite.ivecValue           = { 4, 5, 6 };
    valuesToWrite.legacyRVecValue[XX] = 7._real;
    valuesToWrite.legacyRVecValue[YY] = 8._real;
    valuesToWrite.legacyRVecValue[ZZ] = 9._real;
    valuesToWrite.rvecValue           = { 10._real, 11._real, 12._real };
    valuesToWrite.charArray           = { 'a', 'b', 'c' };
    valuesToWrite.ucharArray          = { 'd', 'e', 'f' };
    valuesToWrite.byteArray           = { static_cast<std::byte>('g'),
                                          static_cast<std::byte>('h'),
                                          static_cast<std::byte>('i') };
    valuesToWrite.int32Array          = { 0x1BCDEF78, 0x654321FE };
    valuesToWrite.int64Array          = { 0x1BCDEF78654321FE, 0x3726ABFEAB34716C };
    valuesToWrite.rvecArray = { { 13._real, 14._real, 15._real }, { 16._real, 17._real, 18._real } };
    // This corner case is tested because GROMACS serializes an empty
    // string with a (sic) size of 1 and a null-terminated empty
    // string, which is side effect of always storing the size and the
    // null terminator.
    valuesToWrite.emptyStringValue = {};
    // This string is serialized with size 8, per the above.
    valuesToWrite.stringValue = "gromacs";

    {
        SCOPED_TRACE("Writing XDR file");
        XdrSerializer serializer(filename_, "w");
        // In GROMACS XDR files, a value describing the precision of
        // the build is written to a header, so that a reader built
        // with either precision can do the right thing. So we do
        // similarly in the tests.
        bool writeAsDouble = GetParam().writeAsDouble;
        serializer.doBool(&writeAsDouble);
        serializer.setDoublePrecision(writeAsDouble);
        valuesToWrite.serialize(&serializer);
    }
    {
        SCOPED_TRACE("Reading XDR file size");
        FILE* file = gmx_ffopen(filename_, "r");

        // Determine file size. This ensures mutual compatibility of
        // all the XDR implementations that GROMACS is built on, and
        // the inter-operability of files written and read by both
        // kinds of GROMACS precision configurations.
        gmx_fseek(file, 0, SEEK_END);
        gmx_off_t fileSize = gmx_ftell(file);
        EXPECT_EQ(fileSize, 224 + (GetParam().writeAsDouble ? 52 : 0));
        gmx_ffclose(file);
    }
    {
        SCOPED_TRACE("Reading XDR file contents");
        XdrSerializer serializer(filename_, "r");
        bool          writtenAsDouble;
        serializer.doBool(&writtenAsDouble);
        serializer.setDoublePrecision(writtenAsDouble);
        // Now real and rvec will read correctly

        SerializerValues valuesToRead;
        if (GMX_DOUBLE && !writtenAsDouble)
        {
            // Double-precision GROMACS never normally writes a
            // non-double-precision file, but we do so in these tests
            // for better test coverage within a single build
            // configuration. Because of those, we need to prevent
            // arithmetic exceptions from overflow from uninitialized
            // values being converted from double precision to float
            // by initializing the otherwise unused values.
            valuesToRead.realValue           = 0._real;
            valuesToRead.legacyRVecValue[XX] = 0._real;
            valuesToRead.legacyRVecValue[YY] = 0._real;
            valuesToRead.legacyRVecValue[ZZ] = 0._real;
            valuesToRead.rvecValue           = { 0._real, 0._real, 0._real };
        }
        valuesToRead.serialize(&serializer);
        EXPECT_EQ(valuesToRead.boolValue, valuesToWrite.boolValue);
        EXPECT_EQ(valuesToRead.unsignedCharValue, valuesToWrite.unsignedCharValue);
        EXPECT_EQ(valuesToRead.charValue, valuesToWrite.charValue);
        EXPECT_EQ(valuesToRead.unsignedShortValue, valuesToWrite.unsignedShortValue);
        EXPECT_EQ(valuesToRead.intValue, valuesToWrite.intValue);
        EXPECT_EQ(valuesToRead.int32Value, valuesToWrite.int32Value);
        EXPECT_EQ(valuesToRead.int64Value, valuesToWrite.int64Value);
        EXPECT_EQ(valuesToRead.floatValue, valuesToWrite.floatValue);
        EXPECT_EQ(valuesToRead.doubleValue, valuesToWrite.doubleValue);
        EXPECT_EQ(valuesToRead.realValue, valuesToWrite.realValue);
        EXPECT_EQ(valuesToRead.legacyIVecValue[XX], valuesToWrite.legacyIVecValue[XX]);
        EXPECT_EQ(valuesToRead.legacyIVecValue[YY], valuesToWrite.legacyIVecValue[YY]);
        EXPECT_EQ(valuesToRead.legacyIVecValue[ZZ], valuesToWrite.legacyIVecValue[ZZ]);
        EXPECT_EQ(valuesToRead.ivecValue[XX], valuesToWrite.ivecValue[XX]);
        EXPECT_EQ(valuesToRead.ivecValue[YY], valuesToWrite.ivecValue[YY]);
        EXPECT_EQ(valuesToRead.ivecValue[ZZ], valuesToWrite.ivecValue[ZZ]);
        EXPECT_EQ(valuesToRead.legacyRVecValue[XX], valuesToWrite.legacyRVecValue[XX]);
        EXPECT_EQ(valuesToRead.legacyRVecValue[YY], valuesToWrite.legacyRVecValue[YY]);
        EXPECT_EQ(valuesToRead.legacyRVecValue[ZZ], valuesToWrite.legacyRVecValue[ZZ]);
        EXPECT_EQ(valuesToRead.rvecValue[XX], valuesToWrite.rvecValue[XX]);
        EXPECT_EQ(valuesToRead.ivecValue[YY], valuesToWrite.ivecValue[YY]);
        EXPECT_EQ(valuesToRead.rvecValue[ZZ], valuesToWrite.rvecValue[ZZ]);
        EXPECT_THAT(valuesToRead.charArray, testing::Pointwise(testing::Eq(), valuesToWrite.charArray));
        EXPECT_THAT(valuesToRead.ucharArray, testing::Pointwise(testing::Eq(), valuesToWrite.ucharArray));
        EXPECT_THAT(valuesToRead.byteArray, testing::Pointwise(testing::Eq(), valuesToWrite.byteArray));
        EXPECT_THAT(valuesToRead.int32Array, testing::Pointwise(testing::Eq(), valuesToWrite.int32Array));
        EXPECT_THAT(valuesToRead.int64Array, testing::Pointwise(testing::Eq(), valuesToWrite.int64Array));
        EXPECT_THAT(valuesToRead.rvecArray, testing::Pointwise(testing::Eq(), valuesToWrite.rvecArray));
        EXPECT_EQ(valuesToRead.emptyStringValue, valuesToWrite.emptyStringValue);
        EXPECT_EQ(valuesToRead.stringValue, valuesToWrite.stringValue);
    }
}

const XdrSerializerTestParameters sc_testParameters[] = {
    { false },
    { true },
};
INSTANTIATE_TEST_SUITE_P(InDifferentPrecisionModes, XdrSerializerTest, ::testing::ValuesIn(sc_testParameters));

} // namespace
} // namespace test
} // namespace gmx
