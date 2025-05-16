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
 * \brief
 * Tests for H5MD scope guard routines
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_guard.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md.h"
#include "gromacs/fileio/h5md/h5md_group.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Determines whether an identifier is valid with the H5Iis_valid method
TEST(H5mdGuardTest, H5mdTypeGuardWorks)
{
    hid_t dataTypeToTest = H5I_INVALID_HID;

    {
        const auto [dataType, dataTypeGuard] = gmx::makeH5mdTypeGuard(H5Tcopy(H5T_NATIVE_INT));
        dataTypeToTest                       = dataType;

        ASSERT_GT(H5Iis_valid(dataTypeToTest), 0);
    }

    ASSERT_LE(H5Iis_valid(dataTypeToTest), 0)
            << "Guard failed: Data type handle should be invalid after exiting scope";
}

TEST(H5mdGuardTest, H5mdGroupGuardWorks)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    hid_t groupToTest = H5I_INVALID_HID;

    {
        const auto [group, groupGuard] =
                gmx::makeH5mdGroupGuard(createGroup(file.fileid(), "group"));
        groupToTest = group;

        ASSERT_GT(H5Iis_valid(groupToTest), 0);
    }

    ASSERT_LE(H5Iis_valid(groupToTest), 0)
            << "Guard failed: Group handle should be invalid after exiting scope";
}

TEST(H5mdGuardTest, H5mdDataSpaceGuardWorks)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    hid_t dataSpaceToTest = H5I_INVALID_HID;

    {
        const auto [dataSpace, dataSpaceGuard] = gmx::makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));
        dataSpaceToTest                        = dataSpace;

        ASSERT_GT(H5Iis_valid(dataSpaceToTest), 0);
    }

    ASSERT_LE(H5Iis_valid(dataSpaceToTest), 0)
            << "Guard failed: Data space handle should be invalid after exiting scope";
}

TEST(H5mdGuardTest, H5mdDataSetGuardWorks)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    hid_t dataSetToTest = H5I_INVALID_HID;

    const auto [dataSpace, dataSpaceGuard] = gmx::makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));

    {
        const auto [dataSet, dataSetGuard] = gmx::makeH5mdDataSetGuard(H5Dcreate(
                file.fileid(), "testDataSet", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        dataSetToTest                      = dataSet;

        ASSERT_GT(H5Iis_valid(dataSetToTest), 0);
    }

    ASSERT_LE(H5Iis_valid(dataSetToTest), 0)
            << "Guard failed: Data set handle should be invalid after exiting scope";
}

TEST(H5mdGuardTest, H5mdAttributeGuardWorks)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    hid_t attributeToTest = H5I_INVALID_HID;

    const auto [dataSpace, dataSpaceGuard] = gmx::makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));

    {
        const auto [attribute, attributeGuard] = gmx::makeH5mdAttributeGuard(H5Acreate(
                file.fileid(), "testAttribute", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT));
        attributeToTest                        = attribute;

        ASSERT_GT(H5Iis_valid(attributeToTest), 0);
    }

    ASSERT_LE(H5Iis_valid(attributeToTest), 0)
            << "Attribute handle should be invalid after exiting scope";
}

TEST(H5mdGuardTest, H5mdPropertyListGuardWorks)
{
    hid_t propertyListToTest = H5I_INVALID_HID;

    {
        const auto [propertyList, propertyListGuard] =
                gmx::makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));
        propertyListToTest = propertyList;

        ASSERT_GT(H5Iis_valid(propertyListToTest), 0);
    }

    ASSERT_LE(H5Iis_valid(propertyListToTest), 0)
            << "Property list handle should be invalid after exiting scope";
}

} // namespace
} // namespace test
} // namespace gmx
