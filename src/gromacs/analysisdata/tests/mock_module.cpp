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
#include "gromacs/analysisdata/dataframe.h"
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
MockAnalysisModule::Impl::startReferenceFrame(const AnalysisDataFrameHeader &header)
{
    EXPECT_TRUE(frameChecker_.get() == NULL);
    EXPECT_EQ(frameIndex_, header.index());
    frameChecker_.reset(new TestReferenceChecker(
        rootChecker_->checkCompound("DataFrame",
                                    formatString("Frame%d", frameIndex_).c_str())));
    frameChecker_->checkReal(header.x(), "X");
}


void
MockAnalysisModule::Impl::checkReferencePoints(const AnalysisDataPointSetRef &points)
{
    EXPECT_TRUE(frameChecker_.get() != NULL);
    if (frameChecker_.get() != NULL)
    {
        // TODO: Add interface to points to make this easier.
        std::vector<real> tmp;
        tmp.reserve(points.columnCount());
        for (int i = 0; i < points.columnCount(); ++i)
        {
            tmp.push_back(points.y(i));
        }
        frameChecker_->checkSequenceArray(tmp.size(), &tmp[0], "Y");
    }
}


void
MockAnalysisModule::Impl::finishReferenceFrame(const AnalysisDataFrameHeader &header)
{
    EXPECT_TRUE(frameChecker_.get() != NULL);
    EXPECT_EQ(frameIndex_, header.index());
    ++frameIndex_;
    frameChecker_.reset();
}


/********************************************************************
 * MockAnalysisModule
 */

namespace
{

void checkHeader(const AnalysisDataFrameHeader &header,
                 const AnalysisDataTestInputFrame &refFrame)
{
    EXPECT_EQ(refFrame.index(), header.index());
    EXPECT_FLOAT_EQ(refFrame.x(), header.x());
    EXPECT_FLOAT_EQ(refFrame.dx(), header.dx());
}

void checkPoints(const AnalysisDataPointSetRef &points,
                 const AnalysisDataTestInputPointSet &refPoints,
                 int columnOffset)
{
    for (int i = 0; i < points.columnCount(); ++i)
    {
        EXPECT_FLOAT_EQ(refPoints.y(points.firstColumn() + columnOffset + i),
                        points.y(i))
            << "  Column: " << i << " (+" << points.firstColumn() << ") / "
            << points.columnCount();
    }
}

void checkFrame(const AnalysisDataFrameRef &frame,
                const AnalysisDataTestInputFrame &refFrame)
{
    checkHeader(frame.header(), refFrame);
    checkPoints(frame.points(), refFrame.points(), 0);
}

class StaticDataFrameHeaderChecker
{
    public:
        StaticDataFrameHeaderChecker(const AnalysisDataTestInputFrame *frame,
                                     int index)
            : frame_(frame), index_(index)
        {
        }

        void operator()(const AnalysisDataFrameHeader &header) const
        {
            SCOPED_TRACE(formatString("Frame %d", index_));
            EXPECT_EQ(index_, header.index());
            checkHeader(header, *frame_);
        }

    private:
        const AnalysisDataTestInputFrame *frame_;
        int                               index_;
};

class StaticDataPointsChecker
{
    public:
        StaticDataPointsChecker(const AnalysisDataTestInputFrame *frame,
                                const AnalysisDataTestInputPointSet *points,
                                int firstcol, int n)
            : frame_(frame), points_(points), firstcol_(firstcol), n_(n)
        {
        }

        void operator()(const AnalysisDataPointSetRef &points) const
        {
            SCOPED_TRACE(formatString("Frame %d", frame_->index()));
            EXPECT_EQ(0, points.firstColumn());
            EXPECT_EQ(n_, points.columnCount());
            checkHeader(points.header(), *frame_);
            checkPoints(points, *points_, firstcol_);
        }

    private:
        const AnalysisDataTestInputFrame *frame_;
        const AnalysisDataTestInputPointSet *points_;
        int                     firstcol_;
        int                     n_;
};


class DataStorageRequester
{
    public:
        explicit DataStorageRequester(int count) : count_(count) {}

        void operator()(AbstractAnalysisData *data) const
        {
            EXPECT_TRUE(data->requestStorage(count_));
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

        void operator()(const AnalysisDataPointSetRef &points) const
        {
            SCOPED_TRACE(formatString("Frame %d", frameIndex_));
            EXPECT_EQ(0, points.firstColumn());
            EXPECT_EQ(data_->columnCount(), points.columnCount());
            checkHeader(points.header(), data_->frame(frameIndex_));
            checkPoints(points, data_->frame(frameIndex_).points(), 0);
            for (int past = 1;
                 (storageCount_ < 0 || past <= storageCount_) && past <= frameIndex_;
                 ++past)
            {
                int   index = frameIndex_ - past;
                SCOPED_TRACE(formatString("Checking storage of frame %d", index));
                ASSERT_NO_THROW({
                    AnalysisDataFrameRef frame = source_->getDataFrame(index);
                    ASSERT_TRUE(frame.isValid());
                    checkFrame(frame, data_->frame(index));
                });
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
    impl_->flags_ |= efAllowMulticolumn | efAllowMultipoint;

    ::testing::InSequence dummy;
    using ::testing::_;
    using ::testing::Invoke;

    EXPECT_CALL(*this, dataStarted(source));
    for (int row = 0; row < data.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame &frame = data.frame(row);
        EXPECT_CALL(*this, frameStarted(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame, row)));
        for (int ps = 0; ps < frame.pointSetCount(); ++ps)
        {
            const AnalysisDataTestInputPointSet &points = frame.points(ps);
            EXPECT_CALL(*this, pointsAdded(_))
                .WillOnce(Invoke(StaticDataPointsChecker(&frame, &points, 0,
                                                         data.columnCount())));
        }
        EXPECT_CALL(*this, frameFinished(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame, row)));
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
    impl_->flags_ |= efAllowMulticolumn | efAllowMultipoint;

    ::testing::InSequence dummy;
    using ::testing::_;
    using ::testing::Invoke;

    EXPECT_CALL(*this, dataStarted(_));
    for (int row = 0; row < data.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame &frame = data.frame(row);
        EXPECT_CALL(*this, frameStarted(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame, row)));
        for (int ps = 0; ps < frame.pointSetCount(); ++ps)
        {
            const AnalysisDataTestInputPointSet &points = frame.points(ps);
            EXPECT_CALL(*this, pointsAdded(_))
                .WillOnce(Invoke(StaticDataPointsChecker(&frame, &points, firstcol, n)));
        }
        EXPECT_CALL(*this, frameFinished(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame, row)));
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
    GMX_RELEASE_ASSERT(!data.isMultipoint() && !source->isMultipoint(),
                       "Storage testing not implemented for multipoint data");
    impl_->flags_ |= efAllowMulticolumn;

    ::testing::InSequence dummy;
    using ::testing::_;
    using ::testing::Invoke;

    EXPECT_CALL(*this, dataStarted(source))
        .WillOnce(Invoke(DataStorageRequester(storageCount)));
    for (int row = 0; row < data.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame &frame = data.frame(row);
        EXPECT_CALL(*this, frameStarted(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame, row)));
        EXPECT_CALL(*this, pointsAdded(_))
            .WillOnce(Invoke(StaticDataPointsStorageChecker(source, &data, row,
                                                            storageCount)));
        EXPECT_CALL(*this, frameFinished(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame, row)));
    }
    EXPECT_CALL(*this, dataFinished());
}


void
MockAnalysisModule::setupReferenceCheck(const TestReferenceChecker &checker,
                                        AbstractAnalysisData *source)
{
    impl_->flags_ |= efAllowMulticolumn | efAllowMultipoint;

    impl_->rootChecker_.reset(new TestReferenceChecker(checker));
    // Google Mock does not support checking the order fully, because
    // the number of frames is not known.
    // Order of the frameStarted(), pointsAdded() and frameFinished()
    // calls is checked using Google Test assertions in the invoked methods.
    using ::testing::_;
    using ::testing::AnyNumber;
    using ::testing::Expectation;
    using ::testing::Invoke;

    Expectation dataStart = EXPECT_CALL(*this, dataStarted(source));
    Expectation frameStart = EXPECT_CALL(*this, frameStarted(_))
        .After(dataStart)
        .WillRepeatedly(Invoke(impl_, &Impl::startReferenceFrame));
    Expectation pointsAdd = EXPECT_CALL(*this, pointsAdded(_))
        .After(dataStart)
        .WillRepeatedly(Invoke(impl_, &Impl::checkReferencePoints));
    Expectation frameFinish = EXPECT_CALL(*this, frameFinished(_))
        .After(dataStart)
        .WillRepeatedly(Invoke(impl_, &Impl::finishReferenceFrame));
    EXPECT_CALL(*this, dataFinished())
        .After(frameStart, pointsAdd, frameFinish);
}

} // namespace test
} // namespace gmx
