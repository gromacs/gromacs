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

#include "mock_module.h"

namespace
{

/********************************************************************
 * Tests for gmx::AnalysisData.
 */

TEST(AnalysisDataTest, BasicInitialization)
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


TEST(AnalysisDataTest, ChecksMultiColumnModules)
{
    gmx::AnalysisData data;
    data.setColumns(2);

    std::auto_ptr<MockModule> mod(new MockModule(0));
    EXPECT_THROW(data.addModule(mod.release()), gmx::APIError);

    mod.reset(new MockModule(gmx::AnalysisDataModuleInterface::efAllowMulticolumn));
    EXPECT_NO_THROW(data.addModule(mod.release()));
}


TEST(AnalysisDataTest, ChecksMultiPointModules)
{
    gmx::AnalysisData data;
    data.setColumns(1, true);

    std::auto_ptr<MockModule> mod(new MockModule(0));
    EXPECT_THROW(data.addModule(mod.release()), gmx::APIError);

    mod.reset(new MockModule(gmx::AnalysisDataModuleInterface::efAllowMultipoint));
    EXPECT_NO_THROW(data.addModule(mod.release()));
}


TEST(AnalysisDataTest, CallsModuleCorrectly)
{
    gmx::AnalysisData data;
    data.setColumns(1);

    std::auto_ptr<MockModule> mod(new MockModule(0));
    {
        ::testing::InSequence dummy;
        using ::testing::_;

        EXPECT_CALL(*mod, dataStarted(&data));
        EXPECT_CALL(*mod, frameStarted(1.0, 0.0));
        EXPECT_CALL(*mod, pointsAdded(1.0, 0.0, 0, 1, _, _, _));
        EXPECT_CALL(*mod, frameFinished());
        EXPECT_CALL(*mod, frameStarted(2.0, 0.0));
        EXPECT_CALL(*mod, pointsAdded(2.0, 0.0, 0, 1, _, _, _));
        EXPECT_CALL(*mod, frameFinished());
        EXPECT_CALL(*mod, frameStarted(3.0, 0.0));
        EXPECT_CALL(*mod, pointsAdded(3.0, 0.0, 0, 1, _, _, _));
        EXPECT_CALL(*mod, frameFinished());
        EXPECT_CALL(*mod, dataFinished());
    }
    ASSERT_NO_THROW(data.addModule(mod.release()));

    gmx::AnalysisDataHandle *dh = NULL;
    ASSERT_NO_THROW(dh = data.startData(NULL));

    ASSERT_NO_THROW(dh->startFrame(0, 1.0));
    EXPECT_NO_THROW(dh->addPoint(0, 1.5));
    EXPECT_NO_THROW(dh->finishFrame());
    EXPECT_EQ(1, data.frameCount());

    ASSERT_NO_THROW(dh->startFrame(1, 2.0));
    EXPECT_NO_THROW(dh->addPoint(0, 2.5));
    EXPECT_NO_THROW(dh->finishFrame());
    EXPECT_EQ(2, data.frameCount());

    ASSERT_NO_THROW(dh->startFrame(2, 3.0));
    EXPECT_NO_THROW(dh->addPoint(0, 3.5));
    EXPECT_NO_THROW(dh->finishFrame());
    EXPECT_EQ(3, data.frameCount());

    EXPECT_NO_THROW(data.finishData(dh));
}

} // namespace
