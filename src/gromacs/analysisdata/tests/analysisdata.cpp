/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Tests for analysis data functionality.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/fatalerror/exceptions.h"

#include "datatest.h"
#include "mock_module.h"

using gmx::test::MockAnalysisModule;

namespace
{

/********************************************************************
 * Tests for gmx::AnalysisData.
 */

TEST(AnalysisDataInitializationTest, BasicInitialization)
{
    gmx::AnalysisData data;
    EXPECT_EQ(0, data.columnCount());
    EXPECT_FALSE(data.isMultipoint());
    EXPECT_EQ(0, data.frameCount());

    data.setColumns(1);
    EXPECT_EQ(1, data.columnCount());
    EXPECT_FALSE(data.isMultipoint());

    data.setColumns(3, true);
    EXPECT_EQ(3, data.columnCount());
    EXPECT_TRUE(data.isMultipoint());

    data.setColumns(1);
    EXPECT_EQ(1, data.columnCount());
    EXPECT_FALSE(data.isMultipoint());
}


TEST(AnalysisDataInitializationTest, ChecksMultiColumnModules)
{
    gmx::AnalysisData data;
    data.setColumns(2);

    std::auto_ptr<MockAnalysisModule> mod(new MockAnalysisModule(0));
    EXPECT_THROW(data.addModule(mod.release()), gmx::APIError);

    mod.reset(new MockAnalysisModule(gmx::AnalysisDataModuleInterface::efAllowMulticolumn));
    EXPECT_NO_THROW(data.addModule(mod.release()));
}


TEST(AnalysisDataInitializationTest, ChecksMultiPointModules)
{
    gmx::AnalysisData data;
    data.setColumns(1, true);

    std::auto_ptr<MockAnalysisModule> mod(new MockAnalysisModule(0));
    EXPECT_THROW(data.addModule(mod.release()), gmx::APIError);

    mod.reset(new MockAnalysisModule(gmx::AnalysisDataModuleInterface::efAllowMultipoint));
    EXPECT_NO_THROW(data.addModule(mod.release()));
}

typedef gmx::test::AnalysisDataTestFixture AnalysisDataTest;

using gmx::test::END_OF_DATA;
using gmx::test::END_OF_FRAME;
static const real inputdata[] = {
    1.0,  0.0, 1.0, 2.0, END_OF_FRAME,
    2.0,  1.0, 1.0, 1.0, END_OF_FRAME,
    3.0,  2.0, 0.0, 0.0, END_OF_FRAME,
    END_OF_DATA
};

TEST_F(AnalysisDataTest, CallsModuleCorrectly)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount());

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(presentAllData(input, &data));
}

TEST_F(AnalysisDataTest, CallsColumnModuleCorrectly)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount());

    ASSERT_NO_THROW(addStaticColumnCheckerModule(input, 0, 2, &data));
    ASSERT_NO_THROW(addStaticColumnCheckerModule(input, 2, 1, &data));
    ASSERT_NO_THROW(presentAllData(input, &data));
}

TEST_F(AnalysisDataTest, CallsModuleCorrectlyWithIndividualPoints)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount());

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addStaticColumnCheckerModule(input, 1, 2, &data));
    gmx::AnalysisDataHandle *handle = NULL;
    ASSERT_NO_THROW(handle = data.startData(NULL));
    for (int row = 0; row < input.frameCount(); ++row)
    {
        const gmx::test::AnalysisDataTestInputFrame &frame = input.frame(row);
        ASSERT_NO_THROW(handle->startFrame(row, frame.x(), frame.dx()));
        for (int column = 0; column < input.columnCount(); ++column)
        {
            ASSERT_NO_THROW(handle->addPoint(column, frame.yptr()[column]));
        }
        ASSERT_NO_THROW(handle->finishFrame());
        EXPECT_EQ(row + 1, data.frameCount());
    }
    ASSERT_NO_THROW(handle->finishData());
}

TEST_F(AnalysisDataTest, CallsModuleCorrectlyWithOutOfOrderFrames)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount());

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addStaticColumnCheckerModule(input, 1, 2, &data));
    gmx::AnalysisDataHandle *handle1 = NULL;
    gmx::AnalysisDataHandle *handle2 = NULL;
    ASSERT_NO_THROW(handle1 = data.startData(NULL));
    ASSERT_NO_THROW(handle2 = data.startData(NULL));
    ASSERT_NO_THROW(presentDataFrame(input, 1, handle1));
    ASSERT_NO_THROW(presentDataFrame(input, 0, handle2));
    ASSERT_NO_THROW(presentDataFrame(input, 2, handle1));
    ASSERT_NO_THROW(handle1->finishData());
    ASSERT_NO_THROW(handle2->finishData());
}

TEST_F(AnalysisDataTest, FullStorageWorks)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount());

    ASSERT_NO_THROW(addStaticStorageCheckerModule(input, -1, &data));
    ASSERT_NO_THROW(presentAllData(input, &data));
}

TEST_F(AnalysisDataTest, CanAddModuleAfterStoredData)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount());
    data.requestStorage(-1);

    ASSERT_NO_THROW(presentAllData(input, &data));
    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
}

TEST_F(AnalysisDataTest, LimitedStorageWorks)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount());

    ASSERT_NO_THROW(addStaticStorageCheckerModule(input, 1, &data));
    ASSERT_NO_THROW(presentAllData(input, &data));
}

} // namespace
