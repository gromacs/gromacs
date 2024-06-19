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
 * Implements classes in datatest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "datatest.h"

#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/analysisdata/tests/mock_datamodule.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

/********************************************************************
 * AnalysisDataTestInputPointSet
 */

AnalysisDataTestInputPointSet::AnalysisDataTestInputPointSet(int index, int dataSetIndex, int firstColumn) :
    index_(index), dataSetIndex_(dataSetIndex), firstColumn_(firstColumn)
{
}


/********************************************************************
 * AnalysisDataTestInputFrame
 */

AnalysisDataTestInputFrame::AnalysisDataTestInputFrame(int index, real x) : index_(index), x_(x) {}

AnalysisDataTestInputPointSet& AnalysisDataTestInputFrame::addPointSet(int dataSet, int firstColumn)
{
    pointSets_.push_back(AnalysisDataTestInputPointSet(pointSets_.size(), dataSet, firstColumn));
    return pointSets_.back();
}

void AnalysisDataTestInputFrame::addPointSetWithValues(int dataSet, int firstColumn, real y1)
{
    AnalysisDataTestInputPointSet& pointSet = addPointSet(dataSet, firstColumn);
    pointSet.addValue(y1);
}

void AnalysisDataTestInputFrame::addPointSetWithValues(int dataSet, int firstColumn, real y1, real y2)
{
    AnalysisDataTestInputPointSet& pointSet = addPointSet(dataSet, firstColumn);
    pointSet.addValue(y1);
    pointSet.addValue(y2);
}

void AnalysisDataTestInputFrame::addPointSetWithValues(int dataSet, int firstColumn, real y1, real y2, real y3)
{
    AnalysisDataTestInputPointSet& pointSet = addPointSet(dataSet, firstColumn);
    pointSet.addValue(y1);
    pointSet.addValue(y2);
    pointSet.addValue(y3);
}

void AnalysisDataTestInputFrame::addPointSetWithValueAndError(int dataSet, int firstColumn, real y1, real e1)
{
    AnalysisDataTestInputPointSet& pointSet = addPointSet(dataSet, firstColumn);
    pointSet.addValueWithError(y1, e1);
}


/********************************************************************
 * AnalysisDataTestInput
 */

AnalysisDataTestInput::AnalysisDataTestInput(int dataSetCount, bool bMultipoint) :
    columnCounts_(dataSetCount), bMultipoint_(bMultipoint)
{
}


AnalysisDataTestInput::~AnalysisDataTestInput() {}


const AnalysisDataTestInputFrame& AnalysisDataTestInput::frame(int index) const
{
    GMX_RELEASE_ASSERT(index >= 0 && index < frameCount(), "Out-of-range frame index");
    return frames_[index];
}


void AnalysisDataTestInput::setColumnCount(int dataSet, int columnCount)
{
    GMX_RELEASE_ASSERT(dataSet >= 0 && dataSet < dataSetCount(), "Out-of-range data set index");
    columnCounts_[dataSet] = columnCount;
}


AnalysisDataTestInputFrame& AnalysisDataTestInput::addFrame(real x)
{
    frames_.push_back(AnalysisDataTestInputFrame(frames_.size(), x));
    return frames_.back();
}

void AnalysisDataTestInput::addFrameWithValues(real x, real y1)
{
    AnalysisDataTestInputFrame& frame = addFrame(x);
    frame.addPointSetWithValues(0, 0, y1);
}

void AnalysisDataTestInput::addFrameWithValues(real x, real y1, real y2)
{
    AnalysisDataTestInputFrame& frame = addFrame(x);
    frame.addPointSetWithValues(0, 0, y1, y2);
}

void AnalysisDataTestInput::addFrameWithValues(real x, real y1, real y2, real y3)
{
    AnalysisDataTestInputFrame& frame = addFrame(x);
    frame.addPointSetWithValues(0, 0, y1, y2, y3);
}

void AnalysisDataTestInput::addFrameWithValueAndError(real x, real y1, real e1)
{
    AnalysisDataTestInputFrame& frame = addFrame(x);
    frame.addPointSetWithValueAndError(0, 0, y1, e1);
}


/********************************************************************
 * AnalysisDataTest
 */

AnalysisDataTestFixture::AnalysisDataTestFixture() {}


void AnalysisDataTestFixture::setupDataObject(const AnalysisDataTestInput& input, AnalysisData* data)
{
    data->setDataSetCount(input.dataSetCount());
    for (int i = 0; i < input.dataSetCount(); ++i)
    {
        data->setColumnCount(i, input.columnCount(i));
    }
    data->setMultipoint(input.isMultipoint());
}


void AnalysisDataTestFixture::presentAllData(const AnalysisDataTestInput& input, AnalysisData* data)
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


void AnalysisDataTestFixture::presentDataFrame(const AnalysisDataTestInput& input, int row, AnalysisDataHandle handle)
{
    const AnalysisDataTestInputFrame& frame = input.frame(row);
    handle.startFrame(row, frame.x(), AnalysisDataTestInputFrame::dx());
    for (int i = 0; i < frame.pointSetCount(); ++i)
    {
        const AnalysisDataTestInputPointSet& points = frame.pointSet(i);
        handle.selectDataSet(points.dataSetIndex());
        for (int j = 0; j < points.size(); ++j)
        {
            if (points.hasError(j))
            {
                handle.setPoint(j + points.firstColumn(),
                                points.y(j),
                                points.error(j),
                                AnalysisDataTestInputPointSet::present(j));
            }
            else
            {
                handle.setPoint(j + points.firstColumn(),
                                points.y(j),
                                AnalysisDataTestInputPointSet::present(j));
            }
        }
        if (input.isMultipoint())
        {
            handle.finishPointSet();
        }
    }
    handle.finishFrame();
}


void AnalysisDataTestFixture::addStaticCheckerModule(const AnalysisDataTestInput& data,
                                                     AbstractAnalysisData*        source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupStaticCheck(data, source, false);
    source->addModule(module);
}


void AnalysisDataTestFixture::addStaticParallelCheckerModule(const AnalysisDataTestInput& data,
                                                             AbstractAnalysisData*        source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupStaticCheck(data, source, true);
    source->addModule(module);
}


void AnalysisDataTestFixture::addStaticColumnCheckerModule(const AnalysisDataTestInput& data,
                                                           int                          firstcol,
                                                           int                          n,
                                                           AbstractAnalysisData*        source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupStaticColumnCheck(data, firstcol, n, source);
    source->addColumnModule(firstcol, n, module);
}


void AnalysisDataTestFixture::addStaticStorageCheckerModule(const AnalysisDataTestInput& data,
                                                            int                   storageCount,
                                                            AbstractAnalysisData* source)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    module->setupStaticStorageCheck(data, storageCount, source);
    source->addModule(module);
}


void AnalysisDataTestFixture::addReferenceCheckerModule(const TestReferenceChecker&   checker,
                                                        const char*                   id,
                                                        AbstractAnalysisData*         source,
                                                        const FloatingPointTolerance& tolerance)
{
    MockAnalysisDataModulePointer module(new MockAnalysisDataModule(0));
    TestReferenceChecker          tmpChecker(checker);
    TestReferenceChecker          compoundChecker(tmpChecker.checkCompound("AnalysisData", id));
    compoundChecker.setDefaultTolerance(tolerance);
    module->setupReferenceCheck(compoundChecker, source);
    source->addModule(module);
}


void AnalysisDataTestFixture::addReferenceCheckerModule(const char* id, AbstractAnalysisData* source)
{
    addReferenceCheckerModule(data_.rootChecker(), id, source, defaultRealTolerance());
}

} // namespace test
} // namespace gmx
