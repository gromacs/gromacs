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
 * \brief Tests for H5MD utility functions.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_util.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_dataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture alias for functions checking whether objects exist.
using H5mdObjectExistsTest = H5mdTestBase;

//! \brief Test fixture alias for functions asserting HDF5 handles.
using H5mdHandleIsValidTest = H5mdTestBase;

TEST_F(H5mdObjectExistsTest, FindsObjectsInFileRoot)
{
    EXPECT_FALSE(objectExists(fileid(), "testGroup"));
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
    EXPECT_TRUE(objectExists(fileid(), "testGroup"));

    EXPECT_FALSE(objectExists(fileid(), "testDataSet"));
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<float>(fileid(), "testDataSet"));
    EXPECT_TRUE(objectExists(fileid(), "testDataSet"));
}

TEST_F(H5mdObjectExistsTest, FindsObjectsInGroups)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));

    EXPECT_FALSE(objectExists(group, "testSubGroup"));
    const auto [subGroup, subGroupGuard] = makeH5mdGroupGuard(createGroup(group, "testSubGroup"));
    EXPECT_TRUE(objectExists(group, "testSubGroup"));

    EXPECT_FALSE(objectExists(group, "testDataSet"));
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<float>(group, "testDataSet"));
    EXPECT_TRUE(objectExists(group, "testDataSet"));
}

TEST_F(H5mdObjectExistsTest, DoesNotSearchForObjectInsideSubGroups)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));

    const auto [subGroup, subGroupGuard] = makeH5mdGroupGuard(createGroup(group, "testSubGroup"));
    EXPECT_FALSE(objectExists(fileid(), "testSubGroup")) << "Must not find sub group in root";
    EXPECT_TRUE(objectExists(group, "testSubGroup")) << "Must find sub group inside parent group";

    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<float>(group, "testDataSet"));
    EXPECT_FALSE(objectExists(fileid(), "testDataSet")) << "Must not find data set in root";
    EXPECT_TRUE(objectExists(group, "testDataSet")) << "Must find data set inside parent group";
}

TEST_F(H5mdObjectExistsTest, SearchesExplicitMultiLevelPath)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "/path/to/testGroup"));
    const auto [subGroup, subGroupGuard] = makeH5mdGroupGuard(createGroup(group, "testSubGroup"));
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<float>(subGroup, "testDataSet"));

    EXPECT_TRUE(objectExists(subGroup, "testDataSet"));
    EXPECT_TRUE(objectExists(group, "testSubGroup/testDataSet"));
    EXPECT_TRUE(objectExists(fileid(), "/path/to/testGroup/testSubGroup/testDataSet"));
}

TEST_F(H5mdHandleIsValidTest, ReturnTrueForValidHandlesToObjects)
{
    EXPECT_TRUE(handleIsValid(fileid()));

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
    EXPECT_TRUE(handleIsValid(group));

    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<float>(fileid(), "testDataSet"));
    EXPECT_TRUE(handleIsValid(dataSet));

    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
    EXPECT_TRUE(handleIsValid(dataType));
}

TEST_F(H5mdHandleIsValidTest, ReturnFalseForInvalidHandles)
{
    EXPECT_FALSE(handleIsValid(H5I_INVALID_HID));
}

TEST_F(H5mdHandleIsValidTest, ReturnFalseForHandlesToClosedObjects)
{
    hid_t handleToClose = H5I_INVALID_HID;
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
        handleToClose                  = group;
        ASSERT_TRUE(handleIsValid(handleToClose))
                << "Sanity check: must be valid before scope guard exits scope";
    }
    EXPECT_FALSE(handleIsValid(handleToClose));
}

} // namespace
} // namespace test
} // namespace gmx
