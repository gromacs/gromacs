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
 * Implements classes in datatest.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "datatest.h"

#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/fatalerror/gmxassert.h"
#include "gromacs/utility/format.h"

#include "testutils/refdata.h"

#include "mock_module.h"

namespace gmx
{

namespace test
{

/********************************************************************
 * AnalysisDataTestInputFrame
 */

AnalysisDataTestInputFrame::AnalysisDataTestInputFrame()
{
}


/********************************************************************
 * AnalysisDataTestInput
 */

namespace
{
    void checkValidDataItem(real data)
    {
        GMX_RELEASE_ASSERT(data != END_OF_DATA && data != END_OF_FRAME,
                           "Inconsistent data array");
    }
}

AnalysisDataTestInput::AnalysisDataTestInput(const real *data)
    : columnCount_(0)
{
    // Count rows and columns.
    int rows = 0, columns = -1;
    const real *dataptr = data;
    for (int i = 0; dataptr[i] != END_OF_DATA; ++i)
    {
        if (dataptr[i] == END_OF_FRAME)
        {
            GMX_RELEASE_ASSERT(i > 0, "Empty data frame");
            if (columns == -1)
            {
                columns = i - 1;
            }
            GMX_RELEASE_ASSERT(columns == i - 1,
                               "Different frames have different number of columns");
            ++rows;
            dataptr += i + 1;
            i = -1;
        }
    }
    GMX_RELEASE_ASSERT(rows > 0, "Empty data");
    columnCount_ = columns;

    // Store the data.
    frames_.resize(rows);
    dataptr = data;
    for (int r = 0; r < rows; ++r)
    {
        AnalysisDataTestInputFrame &frame = frames_[r];
        checkValidDataItem(*dataptr);
        frame.index_ = r;
        frame.x_ = *dataptr;
        ++dataptr;
        frame.y_.resize(columns);
        for (int c = 0; c < columns; ++c)
        {
            checkValidDataItem(dataptr[c]);
            frame.y_[c] = dataptr[c];
        }
        GMX_RELEASE_ASSERT(dataptr[columns] == END_OF_FRAME,
                           "Inconsistent data array");
        dataptr += columns + 1;
    }
    GMX_RELEASE_ASSERT(*dataptr == END_OF_DATA, "Inconsistent data array");
}


AnalysisDataTestInput::~AnalysisDataTestInput()
{
}


const AnalysisDataTestInputFrame &AnalysisDataTestInput::frame(int index) const
{
    GMX_RELEASE_ASSERT(index >= 0 && index < frameCount(),
                       "Out-of-range frame index");
    return frames_[index];
}


/********************************************************************
 * AnalysisDataTest
 */

AnalysisDataTestFixture::AnalysisDataTestFixture()
{
}


void AnalysisDataTestFixture::presentAllData(const AnalysisDataTestInput &input,
                                             AnalysisData *data)
{
    gmx::AnalysisDataHandle *handle = NULL;
    handle = data->startData(NULL);
    for (int row = 0; row < input.frameCount(); ++row)
    {
        presentDataFrame(input, row, handle);
        EXPECT_EQ(row + 1, data->frameCount());
    }
    handle->finishData();
}


void AnalysisDataTestFixture::presentDataFrame(const AnalysisDataTestInput &input,
                                               int row, AnalysisDataHandle *handle)
{
    const AnalysisDataTestInputFrame &frame = input.frame(row);
    handle->startFrame(row, frame.x(), frame.dx());
    handle->addPoints(0, input.columnCount(), frame.yptr(), frame.dyptr(),
                      frame.presentptr());
    handle->finishFrame();
}


void
AnalysisDataTestFixture::addStaticCheckerModule(const AnalysisDataTestInput &data,
                                                AbstractAnalysisData *source)
{
    std::auto_ptr<MockAnalysisModule> module(
            new MockAnalysisModule(gmx::AnalysisDataModuleInterface::efAllowMulticolumn));
    module->setupStaticCheck(data, source);
    source->addModule(module.release());
}


void
AnalysisDataTestFixture::addStaticColumnCheckerModule(const AnalysisDataTestInput &data,
                                                      int firstcol, int n,
                                                      AbstractAnalysisData *source)
{
    std::auto_ptr<MockAnalysisModule> module(
            new MockAnalysisModule(gmx::AnalysisDataModuleInterface::efAllowMulticolumn));
    module->setupStaticColumnCheck(data, firstcol, n, source);
    source->addColumnModule(firstcol, n, module.release());
}


void
AnalysisDataTestFixture::addStaticStorageCheckerModule(const AnalysisDataTestInput &data,
                                                       int storageCount,
                                                       AbstractAnalysisData *source)
{
    std::auto_ptr<MockAnalysisModule> module(
            new MockAnalysisModule(gmx::AnalysisDataModuleInterface::efAllowMulticolumn));
    module->setupStaticStorageCheck(data, storageCount, source);
    source->addModule(module.release());
}


void
AnalysisDataTestFixture::addReferenceCheckerModule(const TestReferenceChecker &checker,
                                                   AbstractAnalysisData *source)
{
    std::auto_ptr<MockAnalysisModule> module(
            new MockAnalysisModule(gmx::AnalysisDataModuleInterface::efAllowMulticolumn));
    module->setupReferenceCheck(checker, source);
    source->addModule(module.release());
}


void
AnalysisDataTestFixture::addReferenceCheckerModule(const char *id,
                                                   AbstractAnalysisData *source)
{
    TestReferenceChecker checker(data_.rootChecker().checkCompound("AnalysisData", id));
    addReferenceCheckerModule(checker, source);
}

} // namespace test
} // namespace gmx
