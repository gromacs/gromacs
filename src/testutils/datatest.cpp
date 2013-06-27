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
 * Implements classes in datatest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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

AnalysisDataTestInputPointSet::AnalysisDataTestInputPointSet(
        int index, int firstColumn)
    : index_(index), firstColumn_(firstColumn)
{
}


/********************************************************************
 * AnalysisDataTestInputFrame
 */

AnalysisDataTestInputFrame::AnalysisDataTestInputFrame(int index, real x)
    : index_(index), x_(x)
{
}

AnalysisDataTestInputPointSet &
AnalysisDataTestInputFrame::addPointSet(int firstColumn)
{
    pointSets_.push_back(AnalysisDataTestInputPointSet(pointSets_.size(), firstColumn));
    return pointSets_.back();
}

void AnalysisDataTestInputFrame::addPointSetWithValues(int firstColumn, real y1)
{
    AnalysisDataTestInputPointSet &pointSet = addPointSet(firstColumn);
    pointSet.addValue(y1);
}

void AnalysisDataTestInputFrame::addPointSetWithValues(int firstColumn, real y1,
                                                       real y2)
{
    AnalysisDataTestInputPointSet &pointSet = addPointSet(firstColumn);
    pointSet.addValue(y1);
    pointSet.addValue(y2);
}

void AnalysisDataTestInputFrame::addPointSetWithValues(int firstColumn, real y1,
                                                       real y2, real y3)
{
    AnalysisDataTestInputPointSet &pointSet = addPointSet(firstColumn);
    pointSet.addValue(y1);
    pointSet.addValue(y2);
    pointSet.addValue(y3);
}


/********************************************************************
 * AnalysisDataTestInput
 */

AnalysisDataTestInput::AnalysisDataTestInput(int columnCount, bool bMultipoint)
    : columnCount_(columnCount), bMultipoint_(bMultipoint)
{
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


AnalysisDataTestInputFrame &AnalysisDataTestInput::addFrame(real x)
{
    frames_.push_back(AnalysisDataTestInputFrame(frames_.size(), x));
    return frames_.back();
}

void AnalysisDataTestInput::addFrameWithValues(real x, real y1)
{
    AnalysisDataTestInputFrame &frame = addFrame(x);
    frame.addPointSetWithValues(0, y1);
}

void AnalysisDataTestInput::addFrameWithValues(real x, real y1, real y2)
{
    AnalysisDataTestInputFrame &frame = addFrame(x);
    frame.addPointSetWithValues(0, y1, y2);
}

void AnalysisDataTestInput::addFrameWithValues(real x, real y1, real y2, real y3)
{
    AnalysisDataTestInputFrame &frame = addFrame(x);
    frame.addPointSetWithValues(0, y1, y2, y3);
}


/********************************************************************
 * AnalysisDataTest
 */

AnalysisDataTestFixture::AnalysisDataTestFixture()
{
}


void AnalysisDataTestFixture::setupDataObject(const AnalysisDataTestInput &input,
                                              AnalysisData                *data)
{
    data->setColumnCount(input.columnCount());
    data->setMultipoint(input.isMultipoint());
}


void AnalysisDataTestFixture::presentAllData(const AnalysisDataTestInput &input,
                                             AnalysisData                *data)
{
    gmx::AnalysisDataParallelOptions options;
    gmx::AnalysisDataHandle          handle = data->startData(options);
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
        const AnalysisDataTestInputPointSet &points = frame.pointSet(i);
        for (int j = 0; j < points.size(); ++j)
        {
            handle.setPoint(j + points.firstColumn(),
                            points.y(j), points.dy(j), points.present(j));
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
                                                AbstractAnalysisData        *source)
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
                                                       int                          storageCount,
                                                       AbstractAnalysisData        *source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupStaticStorageCheck(data, storageCount, source);
    source->addModule(module);
}


void
AnalysisDataTestFixture::addReferenceCheckerModule(TestReferenceChecker  checker,
                                                   const char           *id,
                                                   AbstractAnalysisData *source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupReferenceCheck(checker.checkCompound("AnalysisData", id), source);
    source->addModule(module);
}


void
AnalysisDataTestFixture::addReferenceCheckerModule(const char           *id,
                                                   AbstractAnalysisData *source)
{
    addReferenceCheckerModule(data_.rootChecker(), id, source);
}

} // namespace test
} // namespace gmx
