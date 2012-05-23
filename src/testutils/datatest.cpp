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
 * \ingroup module_testutils
 */
#include "datatest.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mock_datamodule.h"
#include "testutils/refdata.h"

namespace gmx
{

namespace test
{

/********************************************************************
 * AnalysisDataTestInputPointSet
 */

AnalysisDataTestInputPointSet::AnalysisDataTestInputPointSet()
{
}


/********************************************************************
 * AnalysisDataTestInputFrame
 */

AnalysisDataTestInputFrame::AnalysisDataTestInputFrame(int index, real x)
    : index_(index), x_(x)
{
}


/********************************************************************
 * AnalysisDataTestInput
 */

AnalysisDataTestInput::AnalysisDataTestInput(const real *data)
    : columnCount_(0), bMultipoint_(false)
{
    size_t columns = 0;
    const real *dataptr = data;

    while (*dataptr != END_OF_DATA)
    {
        frames_.push_back(AnalysisDataTestInputFrame(frames_.size(), *dataptr));
        AnalysisDataTestInputFrame &frame = frames_.back();
        GMX_RELEASE_ASSERT(*dataptr != END_OF_FRAME && *dataptr != MPSTOP,
                           "Empty data frame");
        while (*dataptr != END_OF_FRAME)
        {
            ++dataptr;
            frame.points_.push_back(AnalysisDataTestInputPointSet());
            AnalysisDataTestInputPointSet &points = frame.points_.back();
            while (*dataptr != MPSTOP && *dataptr != END_OF_FRAME)
            {
                GMX_RELEASE_ASSERT(*dataptr != END_OF_DATA,
                                   "Premature end of data marker");
                points.y_.push_back(*dataptr);
                ++dataptr;
            }
            size_t frameColumns = points.y_.size();
            GMX_RELEASE_ASSERT(frameColumns > 0U, "Empty data point set");
            if (*dataptr == MPSTOP)
            {
                bMultipoint_ = true;
            }
            GMX_RELEASE_ASSERT(!(!bMultipoint_ && columns > 0U && columns != frameColumns),
                               "Different frames have different number of columns");
            if (columns < frameColumns)
            {
                columns = frameColumns;
            }
        }
        ++dataptr;
    }
    GMX_RELEASE_ASSERT(frames_.size() > 0U, "Empty data");
    columnCount_ = columns;
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
    gmx::AnalysisDataParallelOptions options;
    gmx::AnalysisDataHandle handle = data->startData(options);
    for (int row = 0; row < input.frameCount(); ++row)
    {
        presentDataFrame(input, row, handle);
        EXPECT_EQ(row + 1, data->frameCount());
    }
    handle.finishData();
}


void AnalysisDataTestFixture::presentDataFrame(const AnalysisDataTestInput &input,
                                               int row, AnalysisDataHandle handle)
{
    const AnalysisDataTestInputFrame &frame = input.frame(row);
    handle.startFrame(row, frame.x(), frame.dx());
    for (int i = 0; i < frame.pointSetCount(); ++i)
    {
        const AnalysisDataTestInputPointSet &points = frame.points(i);
        for (int j = 0; j < points.size(); ++j)
        {
            handle.setPoint(j, points.y(j), points.dy(j), points.present(j));
        }
        if (input.isMultipoint())
        {
            handle.finishPointSet();
        }
    }
    handle.finishFrame();
}


void
AnalysisDataTestFixture::addStaticCheckerModule(const AnalysisDataTestInput &data,
                                                AbstractAnalysisData *source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupStaticCheck(data, source);
    source->addModule(module);
}


void
AnalysisDataTestFixture::addStaticColumnCheckerModule(const AnalysisDataTestInput &data,
                                                      int firstcol, int n,
                                                      AbstractAnalysisData *source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupStaticColumnCheck(data, firstcol, n, source);
    source->addColumnModule(firstcol, n, module);
}


void
AnalysisDataTestFixture::addStaticStorageCheckerModule(const AnalysisDataTestInput &data,
                                                       int storageCount,
                                                       AbstractAnalysisData *source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupStaticStorageCheck(data, storageCount, source);
    source->addModule(module);
}


void
AnalysisDataTestFixture::addReferenceCheckerModule(TestReferenceChecker checker,
                                                   const char *id,
                                                   AbstractAnalysisData *source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupReferenceCheck(checker.checkCompound("AnalysisData", id), source);
    source->addModule(module);
}


void
AnalysisDataTestFixture::addReferenceCheckerModule(const char *id,
                                                   AbstractAnalysisData *source)
{
    addReferenceCheckerModule(data_.rootChecker(), id, source);
}

} // namespace test
} // namespace gmx
