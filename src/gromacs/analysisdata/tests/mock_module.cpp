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
 * Implements classes in mock_module.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "mock_module.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/fatalerror/gmxassert.h"
#include "gromacs/utility/format.h"

#include "testutils/refdata.h"

#include "datatest.h"
#include "mock_module-impl.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * MockAnalysisModule::Impl
 */

MockAnalysisModule::Impl::Impl(int flags)
    : flags_(flags), frameIndex_(0)
{
}


void
MockAnalysisModule::Impl::checkReferencePoints(real x, real dx, int firstcol, int n,
                                               const real *y, const real *dy,
                                               const bool *present)
{
    TestReferenceChecker frame(
        checker_->checkCompound("DataFrame", formatString("Frame%d", frameIndex_).c_str()));
    ++frameIndex_;
    frame.checkReal(x, "X");
    frame.checkSequenceArray(n, y, "Y");
}


/********************************************************************
 * MockAnalysisModule
 */

namespace
{

void checkFrame(real x, real dx, int firstcol, int n, const real *y,
                const AnalysisDataTestInputFrame &frame)
{
    EXPECT_FLOAT_EQ(frame.x(), x);
    EXPECT_FLOAT_EQ(frame.dx(), dx);
    for (int i = 0; i < n; ++i)
    {
        EXPECT_FLOAT_EQ(frame.y(firstcol + i), y[i]);
    }
}

class StaticDataPointsChecker
{
    public:
        StaticDataPointsChecker(const AnalysisDataTestInputFrame *frame,
                                int firstcol, int n)
            : frame_(frame), firstcol_(firstcol), n_(n)
        {
        }

        void operator()(real x, real dx, int firstcol, int n,
                        const real *y, const real *dy,
                        const bool *present) const
        {
            SCOPED_TRACE(formatString("Frame %d", frame_->index()));
            EXPECT_EQ(0, firstcol);
            EXPECT_EQ(n_, n);
            checkFrame(x, dx, firstcol_ + firstcol, n, y, *frame_);
        }

    private:
        const AnalysisDataTestInputFrame *frame_;
        int                     firstcol_;
        int                     n_;
};


class DataStorageRequester
{
    public:
        explicit DataStorageRequester(int count) : count_(count) {}

        void operator()(AbstractAnalysisData *data) const
        {
            data->requestStorage(count_);
        }

    private:
        int                     count_;
};


class StaticDataPointsStorageChecker
{
    public:
        StaticDataPointsStorageChecker(AbstractAnalysisData *source,
                                       const AnalysisDataTestInput *data,
                                       int frameIndex, int storageCount)
            : source_(source), data_(data),
              frameIndex_(frameIndex), storageCount_(storageCount)
        {
        }

        void operator()(real x, real dx, int firstcol, int n,
                        const real *y, const real *dy,
                        const bool *present) const
        {
            SCOPED_TRACE(formatString("Frame %d", frameIndex_));
            EXPECT_EQ(0, firstcol);
            EXPECT_EQ(data_->columnCount(), n);
            checkFrame(x, dx, firstcol, n, y, data_->frame(frameIndex_));
            for (int past = 0;
                 (storageCount_ < 0 || past <= storageCount_) && past <= frameIndex_;
                 ++past)
            {
                int   index = frameIndex_ - past;
                real  pastx, pastdx;
                const real *pasty;
                SCOPED_TRACE(formatString("Checking storage of frame %d", index));
                ASSERT_TRUE(source_->getDataWErr(index,
                                                 &pastx, &pastdx, &pasty, NULL));
                checkFrame(pastx, pastdx, 0, data_->columnCount(), pasty,
                           data_->frame(index));
                if (past > 0)
                {
                    ASSERT_TRUE(source_->getDataWErr(-past,
                                                     &pastx, &pastdx, &pasty, NULL));
                    checkFrame(pastx, pastdx, 0, data_->columnCount(), pasty,
                               data_->frame(index));
                }
            }
        }

    private:
        AbstractAnalysisData   *source_;
        const AnalysisDataTestInput *data_;
        int                     frameIndex_;
        int                     storageCount_;
};

} // anonymous namespace


MockAnalysisModule::MockAnalysisModule(int flags)
    : impl_(new Impl(flags))
{
}


MockAnalysisModule::~MockAnalysisModule()
{
    delete impl_;
}


int MockAnalysisModule::flags() const
{
    return impl_->flags_;
}


void
MockAnalysisModule::setupStaticCheck(const AnalysisDataTestInput &data,
                                     AbstractAnalysisData *source)
{
    GMX_RELEASE_ASSERT(data.columnCount() == source->columnCount(),
                       "Mismatching data column count");

    ::testing::InSequence dummy;
    using ::testing::_;
    using ::testing::Invoke;

    EXPECT_CALL(*this, dataStarted(source));
    for (int row = 0; row < data.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame &frame = data.frame(row);
        EXPECT_CALL(*this, frameStarted(frame.x(), frame.dx()));
        EXPECT_CALL(*this, pointsAdded(frame.x(), frame.dx(), 0,
                                       data.columnCount(), _, _, _))
            .WillOnce(Invoke(StaticDataPointsChecker(&frame, 0, data.columnCount())));
        EXPECT_CALL(*this, frameFinished());
    }
    EXPECT_CALL(*this, dataFinished());
}


void
MockAnalysisModule::setupStaticColumnCheck(const AnalysisDataTestInput &data,
                                           int firstcol, int n,
                                           AbstractAnalysisData *source)
{
    GMX_RELEASE_ASSERT(data.columnCount() == source->columnCount(),
                       "Mismatching data column count");
    GMX_RELEASE_ASSERT(firstcol >= 0 && n > 0 && firstcol + n <= data.columnCount(),
                       "Out-of-range columns");

    ::testing::InSequence dummy;
    using ::testing::_;
    using ::testing::Invoke;

    EXPECT_CALL(*this, dataStarted(_));
    for (int row = 0; row < data.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame &frame = data.frame(row);
        EXPECT_CALL(*this, frameStarted(frame.x(), frame.dx()));
        EXPECT_CALL(*this, pointsAdded(frame.x(), frame.dx(), 0, n, _, _, _))
            .WillOnce(Invoke(StaticDataPointsChecker(&frame, firstcol, n)));
        EXPECT_CALL(*this, frameFinished());
    }
    EXPECT_CALL(*this, dataFinished());
}


void
MockAnalysisModule::setupStaticStorageCheck(const AnalysisDataTestInput &data,
                                            int storageCount,
                                            AbstractAnalysisData *source)
{
    GMX_RELEASE_ASSERT(data.columnCount() == source->columnCount(),
                       "Mismatching data column count");

    ::testing::InSequence dummy;
    using ::testing::_;
    using ::testing::Invoke;

    EXPECT_CALL(*this, dataStarted(source))
        .WillOnce(Invoke(DataStorageRequester(storageCount)));
    for (int row = 0; row < data.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame &frame = data.frame(row);
        EXPECT_CALL(*this, frameStarted(frame.x(), frame.dx()));
        EXPECT_CALL(*this, pointsAdded(frame.x(), frame.dx(), 0,
                                       data.columnCount(), _, _, _))
            .WillOnce(Invoke(StaticDataPointsStorageChecker(source, &data, row,
                                                            storageCount)));
        EXPECT_CALL(*this, frameFinished());
    }
    EXPECT_CALL(*this, dataFinished());
}


void
MockAnalysisModule::setupReferenceCheck(const TestReferenceChecker &checker,
                                        AbstractAnalysisData *source)
{
    impl_->checker_.reset(new TestReferenceChecker(checker));
    // Google Mock does not support checking the order fully, because
    // the number of frames is not known.
    using ::testing::_;
    using ::testing::AnyNumber;
    using ::testing::Expectation;
    using ::testing::Invoke;

    Expectation dataStart = EXPECT_CALL(*this, dataStarted(source));
    Expectation frameStart = EXPECT_CALL(*this, frameStarted(_, _))
        .Times(AnyNumber())
        .After(dataStart);
    Expectation pointsAdd = EXPECT_CALL(*this, pointsAdded(_, _, _, _, _, _, _))
        .After(dataStart)
        .WillRepeatedly(Invoke(impl_, &Impl::checkReferencePoints));
    Expectation frameFinish = EXPECT_CALL(*this, frameFinished())
        .Times(AnyNumber())
        .After(dataStart);
    EXPECT_CALL(*this, dataFinished())
        .After(frameStart, pointsAdd, frameFinish);
}

} // namespace test
} // namespace gmx
