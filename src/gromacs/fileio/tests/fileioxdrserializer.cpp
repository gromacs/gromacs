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
 * Tests for gmx::FileIOXdrSerializer.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include <cstdint>
#include <cstdio>

#include <filesystem>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"

#include "testutils/testfilemanager.h"

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
constexpr std::int16_t c_int16Value = static_cast<std::int16_t>(0x7A2B);
constexpr std::int32_t c_int32Value = static_cast<std::int32_t>(0x78ABCDEF);
constexpr std::int64_t c_int64Value = static_cast<std::int64_t>(0x78ABCDEF12345678);

constexpr const IntAndFloat32 c_intAndFloat32{ c_int32Value };
constexpr const IntAndFloat64 c_intAndFloat64{ c_int64Value };
//! \}

//! Return the integer used for testing, depending on the size of int.
constexpr int integerSizeDependentTestingValue()
{
    return sizeof(int) == 4 ? c_int32Value : sizeof(int) == 8 ? c_int64Value : c_int16Value;
}

class FileIOXdrSerializerTest : public ::testing::Test
{
public:
    ~FileIOXdrSerializerTest() override
    {
        if (file_)
        {
            gmx_fio_close(file_);
        }
    }
    struct SerializerValues
    {
        bool           boolValue_          = true;
        unsigned char  unsignedCharValue_  = 0x78;
        char           charValue_          = 0x78;
        unsigned short unsignedShortValue_ = static_cast<unsigned short>(c_int16Value);
        std::int32_t   int32Value_         = c_int32Value;
        float          floatValue_         = c_intAndFloat32.floatValue_;
        std::int64_t   int64Value_         = c_int64Value;
        double         doubleValue_        = c_intAndFloat64.doubleValue_;
        int            intValue_           = integerSizeDependentTestingValue();
        real           realValue_          = std::is_same_v<real, double>
                                                     ? static_cast<real>(c_intAndFloat64.doubleValue_)
                                                     : static_cast<real>(c_intAndFloat32.floatValue_);
    } defaultValues_;

    TestFileManager fileManager_;
    // Make sure the file extension is one that gmx_fio_open will
    // recognize to open as binary, even though we're just abusing it
    // to write arbitrary XDR output.
    std::filesystem::path filename_ = fileManager_.getTemporaryFilePath("data.edr");
    t_fileio*             file_     = nullptr;
};

TEST_F(FileIOXdrSerializerTest, SizeIsCorrect)
{
    file_ = gmx_fio_open(filename_, "w");
    FileIOXdrSerializer serializer(file_);
    // These types all have well-defined widths in bytes AFTER XDR serialization,
    // which we can test below.
    serializer.doBool(&defaultValues_.boolValue_);     // 4 bytes
    serializer.doChar(&defaultValues_.charValue_);     // 4 bytes
    serializer.doInt32(&defaultValues_.int32Value_);   // 4 bytes
    serializer.doInt64(&defaultValues_.int64Value_);   // 8 bytes
    serializer.doFloat(&defaultValues_.floatValue_);   // 4 bytes
    serializer.doDouble(&defaultValues_.doubleValue_); // 8 bytes
    std::vector<char> charBuffer = { 'a', 'b', 'c' };
    serializer.doCharArray(charBuffer.data(), charBuffer.size()); // 12 bytes
    serializer.doOpaque(charBuffer.data(), charBuffer.size());    // 4 bytes
    std::vector<int32_t> int32Buffer = { 0x1BCDEF78, 0x654321FE };
    serializer.doInt32Array(int32Buffer.data(), int32Buffer.size()); // 8 bytes
    std::vector<int64_t> int64Buffer = { 0x1BCDEF78654321FE, 0x3726ABFEAB34716C };
    serializer.doInt64Array(int64Buffer.data(), int64Buffer.size()); // 16 bytes
    gmx_fio_close(file_);

    file_ = gmx_fio_open(filename_, "r");

    // Determine file size
    gmx_fseek(gmx_fio_getfp(file_), 0, SEEK_END);
    gmx_off_t fileSize = gmx_fio_ftell(file_);
    EXPECT_EQ(fileSize, 72);
}

} // namespace
} // namespace test
} // namespace gmx
