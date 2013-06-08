/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Tests for analysis data functionality.
 *
 * These tests check the functionality of gmx::AnalysisData, as well as classes
 * used in its implementation: gmx::AbstractAnalysisData and
 * gmx::AnalysisDataStorage.
 * Most checking is done using gmx::test::AnalysisDataTestFixture and mock
 * modules that implement gmx::AnalysisDataModuleInterface.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/datatest.h"
#include "testutils/mock_datamodule.h"
#include "testutils/testasserts.h"

using gmx::test::AnalysisDataTestInput;
using gmx::test::MockAnalysisDataModule;
using gmx::test::MockAnalysisDataModulePointer;

namespace
{

/********************************************************************
 * Tests for gmx::AnalysisData.
 */

/*
 * Tests that simple initialization works.
 */
TEST(AnalysisDataInitializationTest, BasicInitialization)
{
    gmx::AnalysisData data;
    EXPECT_EQ(0, data.columnCount());
    EXPECT_FALSE(data.isMultipoint());
    EXPECT_EQ(0, data.frameCount());

    data.setColumnCount(1);
    EXPECT_EQ(1, data.columnCount());
    EXPECT_FALSE(data.isMultipoint());

    data.setColumnCount(3);
    data.setMultipoint(true);
    EXPECT_EQ(3, data.columnCount());
    EXPECT_TRUE(data.isMultipoint());

    data.setColumnCount(1);
    EXPECT_EQ(1, data.columnCount());
    EXPECT_TRUE(data.isMultipoint());
}

/*
 * Tests that checking for compatibility of modules with multicolumn data
 * works.
 */
TEST(AnalysisDataInitializationTest, ChecksMultiColumnModules)
{
    gmx::AnalysisData data;
    data.setColumnCount(2);

    MockAnalysisDataModulePointer mod1(new MockAnalysisDataModule(0));
    EXPECT_THROW_GMX(data.addModule(mod1), gmx::APIError);

    MockAnalysisDataModulePointer mod2(
            new MockAnalysisDataModule(gmx::AnalysisDataModuleInterface::efAllowMulticolumn));
    EXPECT_NO_THROW_GMX(data.addModule(mod2));
}

/*
 * Tests that checking for compatibility of modules with multipoint data
 * works.
 */
TEST(AnalysisDataInitializationTest, ChecksMultiPointModules)
{
    gmx::AnalysisData data;
    data.setColumnCount(1);
    data.setMultipoint(true);

    MockAnalysisDataModulePointer mod1(new MockAnalysisDataModule(0));
    EXPECT_THROW_GMX(data.addModule(mod1), gmx::APIError);

    MockAnalysisDataModulePointer mod2(
            new MockAnalysisDataModule(gmx::AnalysisDataModuleInterface::efAllowMultipoint));
    EXPECT_NO_THROW_GMX(data.addModule(mod2));
}


//! Test fixture for gmx::AnalysisData.
typedef gmx::test::AnalysisDataTestFixture AnalysisDataTest;

// Basic input data for gmx::AnalysisData tests.
class SimpleInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
            static SimpleInputData singleton;
            return singleton.data_;
        }

        SimpleInputData() : data_(3, false)
        {
            data_.addFrameWithValues(1.0,  0.0, 1.0, 2.0);
            data_.addFrameWithValues(2.0,  1.0, 1.0, 1.0);
            data_.addFrameWithValues(3.0,  2.0, 0.0, 0.0);
        }

    private:
        AnalysisDataTestInput  data_;
};

// Input data for multipoint gmx::AnalysisData tests.
class MultipointInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
            static MultipointInputData singleton;
            return singleton.data_;
        }

        MultipointInputData() : data_(3, true)
        {
            using gmx::test::AnalysisDataTestInputFrame;
            AnalysisDataTestInputFrame &frame1 = data_.addFrame(1.0);
            frame1.addPointSetWithValues(0, 0.0, 1.0, 2.0);
            frame1.addPointSetWithValues(0, 1.1, 2.1, 1.1);
            frame1.addPointSetWithValues(0, 2.2, 1.2, 0.2);
            AnalysisDataTestInputFrame &frame2 = data_.addFrame(2.0);
            frame2.addPointSetWithValues(1, 1.0, 1.0);
            frame2.addPointSetWithValues(0, 2.1, 1.1, 0.1);
            frame2.addPointSetWithValues(2, 1.2);
            AnalysisDataTestInputFrame &frame3 = data_.addFrame(3.0);
            frame3.addPointSetWithValues(0, 2.0, 0.0, 0.0);
            frame3.addPointSetWithValues(0, 3.1, 2.1);
            frame3.addPointSetWithValues(1, 2.2, 1.2);
        }

    private:
        AnalysisDataTestInput  data_;
};

/*
 * Tests that data is forwarded correctly to modules using two independent
 * modules.
 */
TEST_F(AnalysisDataTest, CallsModuleCorrectly)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

/*
 * Tests that data is forwarded correctly to modules that are added using
 * addColumnModule().
 * Uses two independent modules.
 */
TEST_F(AnalysisDataTest, CallsColumnModuleCorrectly)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    ASSERT_NO_THROW_GMX(addStaticColumnCheckerModule(input, 0, 2, &data));
    ASSERT_NO_THROW_GMX(addStaticColumnCheckerModule(input, 2, 1, &data));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

/*
 * Tests that data is forwarded correctly (in frame order) to modules when the
 * data is added through multiple handles in non-increasing order.
 */
TEST_F(AnalysisDataTest, CallsModuleCorrectlyWithOutOfOrderFrames)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addStaticColumnCheckerModule(input, 1, 2, &data));
    gmx::AnalysisDataHandle          handle1;
    gmx::AnalysisDataHandle          handle2;
    gmx::AnalysisDataParallelOptions options(2);
    ASSERT_NO_THROW_GMX(handle1 = data.startData(options));
    ASSERT_NO_THROW_GMX(handle2 = data.startData(options));
    ASSERT_NO_THROW_GMX(presentDataFrame(input, 1, handle1));
    ASSERT_NO_THROW_GMX(presentDataFrame(input, 0, handle2));
    ASSERT_NO_THROW_GMX(presentDataFrame(input, 2, handle1));
    ASSERT_NO_THROW_GMX(handle1.finishData());
    ASSERT_NO_THROW_GMX(handle2.finishData());
}

/*
 * Tests that data can be accessed correctly from a module that requests
 * storage using AbstractAnalysisData::requestStorage() with parameter -1.
 */
TEST_F(AnalysisDataTest, FullStorageWorks)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    ASSERT_NO_THROW_GMX(addStaticStorageCheckerModule(input, -1, &data));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

/*
 * Tests that a data module can be added to an AnalysisData object after data
 * has been added if all data is still available in storage.
 */
TEST_F(AnalysisDataTest, CanAddModuleAfterStoredData)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));
    ASSERT_TRUE(data.requestStorage(-1));

    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
}

/*
 * Tests that data can be accessed correctly from a module that requests
 * storage using AbstractAnalysisData::requestStorage() only for one frame.
 */
TEST_F(AnalysisDataTest, LimitedStorageWorks)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    ASSERT_NO_THROW_GMX(addStaticStorageCheckerModule(input, 1, &data));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

/*
 * Tests that multipoint data is forwarded correctly to modules using two
 * independent modules.
 */
TEST_F(AnalysisDataTest, MultipointCallsModuleCorrectly)
{
    const AnalysisDataTestInput &input = MultipointInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

/*
 * Tests that multipoint data is forwarded correctly to modules that are added
 * using addColumnModule().
 * Uses two independent modules.
 */
TEST_F(AnalysisDataTest, MultipointCallsColumnModuleCorrectly)
{
    const AnalysisDataTestInput &input = MultipointInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    ASSERT_NO_THROW_GMX(addStaticColumnCheckerModule(input, 0, 2, &data));
    ASSERT_NO_THROW_GMX(addStaticColumnCheckerModule(input, 2, 1, &data));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

} // namespace
