/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
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
 * Tests for analysis data functionality.
 *
 * These tests check the functionality of gmx::AnalysisData, as well as classes
 * used in its implementation: gmx::AbstractAnalysisData and
 * gmx::AnalysisDataStorage.
 * Most checking is done using gmx::test::AnalysisDataTestFixture and mock
 * modules that implement gmx::IAnalysisDataModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/analysisdata.h"

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/analysisdata/tests/datatest.h"
#include "gromacs/analysisdata/tests/mock_datamodule.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

using gmx::test::AnalysisDataTestInput;
using gmx::test::MockAnalysisDataModule;
using gmx::test::MockAnalysisDataModulePointer;

namespace gmx
{
namespace test
{
namespace
{

/********************************************************************
 * Tests for gmx::AnalysisData without any actual data.
 */

/*
 * Tests that simple initialization works.
 */
TEST(AnalysisDataInitializationTest, BasicInitialization)
{
    gmx::AnalysisData data;
    EXPECT_EQ(1, data.dataSetCount());
    EXPECT_EQ(0, data.columnCount(0));
    EXPECT_EQ(0, data.columnCount());
    EXPECT_FALSE(data.isMultipoint());
    EXPECT_EQ(0, data.frameCount());

    data.setColumnCount(0, 1);
    EXPECT_EQ(1, data.columnCount(0));
    EXPECT_EQ(1, data.columnCount());
    EXPECT_FALSE(data.isMultipoint());

    data.setDataSetCount(2);
    EXPECT_EQ(2, data.dataSetCount());
    data.setColumnCount(0, 3);
    EXPECT_EQ(3, data.columnCount(0));
    EXPECT_EQ(0, data.columnCount(1));
    data.setColumnCount(1, 2);
    EXPECT_EQ(3, data.columnCount(0));
    EXPECT_EQ(2, data.columnCount(1));

    data.setDataSetCount(1);
    EXPECT_EQ(1, data.dataSetCount());
    data.setMultipoint(true);
    EXPECT_EQ(3, data.columnCount());
    EXPECT_TRUE(data.isMultipoint());

    data.setColumnCount(0, 1);
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
    data.setColumnCount(0, 2);

    MockAnalysisDataModulePointer mod1(new MockAnalysisDataModule(0));
    EXPECT_THROW_GMX(data.addModule(mod1), gmx::APIError);

    MockAnalysisDataModulePointer mod2(
            new MockAnalysisDataModule(gmx::IAnalysisDataModule::efAllowMulticolumn));
    EXPECT_NO_THROW_GMX(data.addModule(mod2));
}

/*
 * Tests that checking for compatibility of modules with multipoint data
 * works.
 */
TEST(AnalysisDataInitializationTest, ChecksMultipointModules)
{
    gmx::AnalysisData data;
    data.setColumnCount(0, 1);
    data.setMultipoint(true);

    MockAnalysisDataModulePointer mod1(new MockAnalysisDataModule(0));
    EXPECT_THROW_GMX(data.addModule(mod1), gmx::APIError);

    MockAnalysisDataModulePointer mod2(
            new MockAnalysisDataModule(gmx::IAnalysisDataModule::efAllowMultipoint));
    EXPECT_NO_THROW_GMX(data.addModule(mod2));
}

#if GTEST_HAS_TYPED_TEST

/********************************************************************
 * Input data for tests below.
 */

// Basic input data for gmx::AnalysisData tests.
class SimpleInputData
{
public:
    static const AnalysisDataTestInput& get()
    {
        static SimpleInputData singleton;
        return singleton.data_;
    }

    SimpleInputData() : data_(1, false)
    {
        data_.setColumnCount(0, 3);
        data_.addFrameWithValues(1.0, 0.0, 1.0, 2.0);
        data_.addFrameWithValues(2.0, 1.0, 1.0, 1.0);
        data_.addFrameWithValues(3.0, 2.0, 0.0, 0.0);
    }

private:
    AnalysisDataTestInput data_;
};

// Input data with multiple data sets for gmx::AnalysisData tests.
class DataSetsInputData
{
public:
    static const AnalysisDataTestInput& get()
    {
        static DataSetsInputData singleton;
        return singleton.data_;
    }

    DataSetsInputData() : data_(2, false)
    {
        using gmx::test::AnalysisDataTestInputFrame;
        data_.setColumnCount(0, 3);
        data_.setColumnCount(1, 2);
        AnalysisDataTestInputFrame& frame1 = data_.addFrame(1.0);
        frame1.addPointSetWithValues(0, 0, 0.0, 1.0, 2.0);
        frame1.addPointSetWithValues(1, 0, 2.1, 1.1);
        AnalysisDataTestInputFrame& frame2 = data_.addFrame(2.0);
        frame2.addPointSetWithValues(0, 0, 1.0, 1.0, 1.0);
        frame2.addPointSetWithValues(1, 0, 0.1, 2.1);
        AnalysisDataTestInputFrame& frame3 = data_.addFrame(3.0);
        frame3.addPointSetWithValues(0, 0, 2.0, 0.0, 0.0);
        frame3.addPointSetWithValues(1, 0, 1.1, 1.1);
    }

private:
    AnalysisDataTestInput data_;
};

// Input data for multipoint gmx::AnalysisData tests.
class MultipointInputData
{
public:
    static const AnalysisDataTestInput& get()
    {
        static MultipointInputData singleton;
        return singleton.data_;
    }

    MultipointInputData() : data_(1, true)
    {
        using gmx::test::AnalysisDataTestInputFrame;
        data_.setColumnCount(0, 3);
        AnalysisDataTestInputFrame& frame1 = data_.addFrame(1.0);
        frame1.addPointSetWithValues(0, 0, 0.0, 1.0, 2.0);
        frame1.addPointSetWithValues(0, 0, 1.1, 2.1, 1.1);
        frame1.addPointSetWithValues(0, 0, 2.2, 1.2, 0.2);
        AnalysisDataTestInputFrame& frame2 = data_.addFrame(2.0);
        frame2.addPointSetWithValues(0, 1, 1.0, 1.0);
        frame2.addPointSetWithValues(0, 0, 2.1, 1.1, 0.1);
        frame2.addPointSetWithValues(0, 2, 1.2);
        AnalysisDataTestInputFrame& frame3 = data_.addFrame(3.0);
        frame3.addPointSetWithValues(0, 0, 2.0, 0.0, 0.0);
        frame3.addPointSetWithValues(0, 0, 3.1, 2.1);
        frame3.addPointSetWithValues(0, 1, 2.2, 1.2);
    }

private:
    AnalysisDataTestInput data_;
};

// Input data with multiple multipoint data sets for gmx::AnalysisData tests.
class MultipointDataSetsInputData
{
public:
    static const AnalysisDataTestInput& get()
    {
        static MultipointDataSetsInputData singleton;
        return singleton.data_;
    }

    MultipointDataSetsInputData() : data_(2, true)
    {
        using gmx::test::AnalysisDataTestInputFrame;
        data_.setColumnCount(0, 3);
        data_.setColumnCount(1, 2);
        AnalysisDataTestInputFrame& frame1 = data_.addFrame(1.0);
        frame1.addPointSetWithValues(0, 0, 0.0, 1.0, 2.0);
        frame1.addPointSetWithValues(0, 1, 2.1, 1.1);
        frame1.addPointSetWithValues(1, 0, 2.01, 1.01);
        frame1.addPointSetWithValues(1, 1, 0.11);
        AnalysisDataTestInputFrame& frame2 = data_.addFrame(2.0);
        frame2.addPointSetWithValues(0, 0, 1.0, 1.0, 1.0);
        frame2.addPointSetWithValues(0, 0, 0.1, 2.1);
        frame2.addPointSetWithValues(1, 1, 1.01);
        AnalysisDataTestInputFrame& frame3 = data_.addFrame(3.0);
        frame3.addPointSetWithValues(0, 0, 2.0, 0.0, 0.0);
        frame3.addPointSetWithValues(0, 1, 1.1);
    }

private:
    AnalysisDataTestInput data_;
};

/********************************************************************
 * Tests for gmx::AnalysisData that require data.
 */

using gmx::test::AnalysisDataTestFixture;

class AnalysisDataTest : public AnalysisDataTestFixture
{
public:
    explicit AnalysisDataTest(const AnalysisDataTestInput& input) : input_(input) {}

    void SetUp() override { ASSERT_NO_THROW_GMX(setupDataObject(input_, &data_)); }

    void addStaticCheckerModule()
    {
        AnalysisDataTestFixture::addStaticCheckerModule(input_, &data_);
    }
    void addStaticParallelCheckerModule()
    {
        AnalysisDataTestFixture::addStaticParallelCheckerModule(input_, &data_);
    }
    void addStaticColumnCheckerModule(int firstColumn, int columnCount)
    {
        AnalysisDataTestFixture::addStaticColumnCheckerModule(input_, firstColumn, columnCount, &data_);
    }
    void addStaticStorageCheckerModule(int storageCount)
    {
        AnalysisDataTestFixture::addStaticStorageCheckerModule(input_, storageCount, &data_);
    }
    void presentAllData() { AnalysisDataTestFixture::presentAllData(input_, &data_); }

    const AnalysisDataTestInput& input_;
    gmx::AnalysisData            data_;
};

template<class InputDataType>
class AnalysisDataCommonTest : public AnalysisDataTest
{
public:
    AnalysisDataCommonTest() : AnalysisDataTest(InputDataType::get()) {}
};

//! Test fixture for tests that are only applicable to simple data.
typedef AnalysisDataCommonTest<SimpleInputData> AnalysisDataSimpleTest;
//! Test fixture for tests that are only applicable to multipoint data.
typedef AnalysisDataCommonTest<MultipointInputData> AnalysisDataMultipointTest;
//! List of input data types for tests applicable to all types of data.
typedef ::testing::Types<SimpleInputData, DataSetsInputData, MultipointInputData, MultipointDataSetsInputData> AllInputDataTypes;
TYPED_TEST_SUITE(AnalysisDataCommonTest, AllInputDataTypes);

/*
 * Tests that data is forwarded correctly to modules using two independent
 * modules.
 */
TYPED_TEST(AnalysisDataCommonTest, CallsModuleCorrectly)
{
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticCheckerModule());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticCheckerModule());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentAllData());
}

/*
 * Tests that data is forwarded correctly to modules when there are only
 * parallel modules.
 */
TYPED_TEST(AnalysisDataCommonTest, CallsParallelModuleCorrectly)
{
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticParallelCheckerModule());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticParallelCheckerModule());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentAllData());
}

/*
 * Tests that data is forwarded correctly to modules when there are both
 * parallel and serial modules.
 */
TYPED_TEST(AnalysisDataCommonTest, CallsMixedModulesCorrectly)
{
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticCheckerModule());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticParallelCheckerModule());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentAllData());
}

/*
 * Tests that data is forwarded correctly to modules that are added using
 * addColumnModule().
 * Uses two independent modules.
 */
TYPED_TEST(AnalysisDataCommonTest, CallsColumnModuleCorrectly)
{
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticColumnCheckerModule(0, 2));
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticColumnCheckerModule(2, 1));
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentAllData());
}

/*
 * Tests that data is forwarded correctly (in frame order) to modules when the
 * data is added through multiple handles in non-increasing order.
 */
TYPED_TEST(AnalysisDataCommonTest, CallsModuleCorrectlyWithOutOfOrderFrames)
{
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticCheckerModule());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticParallelCheckerModule());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticColumnCheckerModule(1, 2));
    gmx::AnalysisDataHandle          handle1;
    gmx::AnalysisDataHandle          handle2;
    gmx::AnalysisDataParallelOptions options(2);
    ASSERT_NO_THROW_GMX(handle1 = this->data_.startData(options));
    ASSERT_NO_THROW_GMX(handle2 = this->data_.startData(options));
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentDataFrame(this->input_, 1, handle1));
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentDataFrame(this->input_, 0, handle2));
    ASSERT_NO_THROW_GMX(this->data_.finishFrameSerial(0));
    ASSERT_NO_THROW_GMX(this->data_.finishFrameSerial(1));
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentDataFrame(this->input_, 2, handle1));
    ASSERT_NO_THROW_GMX(this->data_.finishFrameSerial(2));
    ASSERT_NO_THROW_GMX(handle1.finishData());
    ASSERT_NO_THROW_GMX(handle2.finishData());
}

/*
 * Tests that data can be accessed correctly from a module that requests
 * storage using AbstractAnalysisData::requestStorage() with parameter -1.
 */
TYPED_TEST(AnalysisDataCommonTest, FullStorageWorks)
{
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticStorageCheckerModule(-1));
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentAllData());
}

/*
 * Tests that a data module can be added to an AnalysisData object after data
 * has been added if all data is still available in storage.
 */
TYPED_TEST(AnalysisDataCommonTest, CanAddModuleAfterStoredData)
{
    ASSERT_TRUE(this->data_.requestStorage(-1));

    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentAllData());
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticCheckerModule());
}

/*
 * Tests that data can be accessed correctly from a module that requests
 * storage using AbstractAnalysisData::requestStorage() only for one frame.
 */
TYPED_TEST(AnalysisDataCommonTest, LimitedStorageWorks)
{
    ASSERT_NO_THROW_GMX(AnalysisDataTest::addStaticStorageCheckerModule(1));
    ASSERT_NO_THROW_GMX(AnalysisDataTest::presentAllData());
}

#else

/* A dummy test that at least signals that something is missing if one runs the
 * unit test executable itself.
 */
TEST(DISABLED_AnalysisDataCommonTest, GenericTests)
{
    ADD_FAILURE() << "Tests for generic AnalysisData functionality require support for "
                  << "Google Test typed tests, which was not available when the tests "
                  << "were compiled.";
}

#endif

} // namespace
} // namespace test
} // namespace gmx
