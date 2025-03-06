/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Tests for H5MD data type utility routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_type.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture for all primitive types
template<typename ValueType>
class H5mdTypeTest : public ::testing::Test
{
};

//! \brief List of primitives for which to create tests
using PrimitiveList = ::testing::Types<int32_t, int64_t, float, double>;

// Set up suites for testing of all relevant types
TYPED_TEST_SUITE(H5mdTypeTest, PrimitiveList);

TEST(H5mdTypeTest, ValueTypesToHdf5TypesWorkForAllTypes)
{
    // This test may be called before the HDF5 library has been initialized with a call to H5open.
    // Thus, we also call it here to ensure that it is active for all tests in this file.
    H5open();

    // We test all types for match using H5Tequal, which returns a positive value if they are equal
    EXPECT_GT(H5Tequal(hdf5DataTypeFor<int32_t>(), H5T_NATIVE_INT32), 0) << "Mismatch for Int32";
    EXPECT_GT(H5Tequal(hdf5DataTypeFor<int64_t>(), H5T_NATIVE_INT64), 0) << "Mismatch for Int64";
    EXPECT_GT(H5Tequal(hdf5DataTypeFor<float>(), H5T_NATIVE_FLOAT), 0) << "Mismatch for Float";
    EXPECT_GT(H5Tequal(hdf5DataTypeFor<double>(), H5T_NATIVE_DOUBLE), 0) << "Mismatch for Double";
}


TYPED_TEST(H5mdTypeTest, ValueTypeIsDataTypeWorksForAllTypes)
{
    {
        SCOPED_TRACE("Int32");
        const bool expectedMatch = std::is_same_v<TypeParam, int32_t>;
        const bool actualMatch   = valueTypeIsDataType<TypeParam>(H5T_NATIVE_INT32);

        EXPECT_EQ(actualMatch, expectedMatch);
    }
    {
        SCOPED_TRACE("Int64");
        const bool expectedMatch = std::is_same_v<TypeParam, int64_t>;
        const bool actualMatch   = valueTypeIsDataType<TypeParam>(H5T_NATIVE_INT64);

        EXPECT_EQ(actualMatch, expectedMatch);
    }
    {
        SCOPED_TRACE("Float");
        const bool expectedMatch = std::is_same_v<TypeParam, float>;
        const bool actualMatch   = valueTypeIsDataType<TypeParam>(H5T_NATIVE_FLOAT);

        EXPECT_EQ(actualMatch, expectedMatch);
    }
    {
        SCOPED_TRACE("Double");
        const bool expectedMatch = std::is_same_v<TypeParam, double>;
        const bool actualMatch   = valueTypeIsDataType<TypeParam>(H5T_NATIVE_DOUBLE);

        EXPECT_EQ(actualMatch, expectedMatch);
    }
}

TYPED_TEST(H5mdTypeTest, ValueTypeIsDataTypeWorksForAllTypesWithBufferTypeInference)
{
    const TypeParam                 value = -1;
    const ArrayRef<const TypeParam> array = arrayRefFromArray(&value, 1);

    {
        SCOPED_TRACE("Int32");
        const bool expectedMatch = std::is_same_v<TypeParam, int32_t>;
        const bool actualMatch   = valueTypeIsDataType(H5T_NATIVE_INT32, array);

        EXPECT_EQ(actualMatch, expectedMatch);
    }
    {
        SCOPED_TRACE("Int64");
        const bool expectedMatch = std::is_same_v<TypeParam, int64_t>;
        const bool actualMatch   = valueTypeIsDataType(H5T_NATIVE_INT64, array);

        EXPECT_EQ(actualMatch, expectedMatch);
    }
    {
        SCOPED_TRACE("Float");
        const bool expectedMatch = std::is_same_v<TypeParam, float>;
        const bool actualMatch   = valueTypeIsDataType(H5T_NATIVE_FLOAT, array);

        EXPECT_EQ(actualMatch, expectedMatch);
    }
    {
        SCOPED_TRACE("Double");
        const bool expectedMatch = std::is_same_v<TypeParam, double>;
        const bool actualMatch   = valueTypeIsDataType(H5T_NATIVE_DOUBLE, array);

        EXPECT_EQ(actualMatch, expectedMatch);
    }
}

} // namespace
} // namespace test
} // namespace gmx
