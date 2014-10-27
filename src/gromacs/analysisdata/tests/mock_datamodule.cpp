/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * Implements classes in mock_datamodule.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "mock_datamodule.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "gromacs/analysisdata/tests/datatest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * MockAnalysisDataModule::Impl
 */

/*! \internal \brief
 * Private implementation class for gmx::test::MockAnalysisDataModule.
 *
 * \ingroup module_analysisdata
 */
class MockAnalysisDataModule::Impl
{
    public:
        //! Initializes a mock object with the given flags.
        explicit Impl(int flags);

        /*! \brief
         * Callback used to initialize reference data checks
         *
         * Called in response to dataStarted().
         * Records the source data for later use (for access to data properties).
         */
        void startReferenceData(AbstractAnalysisData *data);
        /*! \brief
         * Callback used to check frame start against reference data.
         *
         * Called to check parameters and order of calls to frameStarted().
         * In addition to reference data checks, this method checks statically
         * that the new frame matches \a frameIndex_.
         */
        void startReferenceFrame(const AnalysisDataFrameHeader &header);
        /*! \brief
         * Callback used to check frame points against reference data.
         *
         * Called to check parameters and order of calls to pointsAdded().
         */
        void checkReferencePoints(const AnalysisDataPointSetRef &points);
        /*! \brief
         * Callback used to check frame finish against reference data.
         *
         * Called to check parameters and order of calls to frameFinished().
         */
        void finishReferenceFrame(const AnalysisDataFrameHeader &header);
        /*! \brief
         * Callback used to check serial frame finish with reference data.
         *
         * Called to check parameters and order of calls to
         * frameFinishedSerial().
         * \a frameIndex_ is incremented here.
         */
        void finishReferenceFrameSerial(int frameIndex);

        /*! \brief
         * Reference data checker to use for checking frames.
         *
         * Must be non-NULL if startReferenceFrame() is called.
         */
        boost::scoped_ptr<TestReferenceChecker>  rootChecker_;
        /*! \brief
         * Reference data checker to use to check the current frame.
         *
         * Non-NULL between startReferenceFrame() and finishReferenceFrame()
         * calls.
         */
        boost::scoped_ptr<TestReferenceChecker>  frameChecker_;
        //! Source data.
        const AbstractAnalysisData              *source_;
        //! Flags that will be returned by the mock module.
        int                                      flags_;
        //! Index of the current/next frame.
        int                                      frameIndex_;
};

namespace
{

/*! \internal \brief
 * Checks a single AnalysisDataValue.
 *
 * \ingroup module_analysisdata
 */
void checkReferenceDataPoint(TestReferenceChecker    *checker,
                             const AnalysisDataValue &value)
{
    TestReferenceChecker compound(checker->checkCompound("DataValue", NULL));
    compound.checkReal(value.value(), "Value");
    if (compound.checkPresent(value.hasError(), "Error"))
    {
        compound.checkReal(value.error(), "Error");
    }
    if (compound.checkPresent(!value.isPresent(), "Present"))
    {
        compound.checkBoolean(value.isPresent(), "Present");
    }
}

}       // namespace

MockAnalysisDataModule::Impl::Impl(int flags)
    : source_(NULL), flags_(flags), frameIndex_(0)
{
}


void
MockAnalysisDataModule::Impl::startReferenceData(AbstractAnalysisData *data)
{
    source_ = data;
}


void
MockAnalysisDataModule::Impl::startReferenceFrame(
        const AnalysisDataFrameHeader &header)
{
    GMX_RELEASE_ASSERT(rootChecker_.get() != NULL,
                       "Root checker not set, but reference data used");
    EXPECT_TRUE(frameChecker_.get() == NULL);
    EXPECT_EQ(frameIndex_, header.index());
    frameChecker_.reset(new TestReferenceChecker(
                                rootChecker_->checkCompound("DataFrame",
                                                            formatString("Frame%d", frameIndex_).c_str())));
    frameChecker_->checkReal(header.x(), "X");
}


void
MockAnalysisDataModule::Impl::checkReferencePoints(
        const AnalysisDataPointSetRef &points)
{
    EXPECT_TRUE(frameChecker_.get() != NULL);
    if (frameChecker_.get() != NULL)
    {
        TestReferenceChecker checker(
                frameChecker_->checkCompound("DataValues", NULL));
        checker.checkInteger(points.columnCount(), "Count");
        if (checker.checkPresent(source_->dataSetCount() > 1, "DataSet"))
        {
            checker.checkInteger(points.dataSetIndex(), "DataSet");
        }
        const int  sourceColumnCount = source_->columnCount(points.dataSetIndex());
        const bool bAllColumns       = (points.firstColumn() == 0
                                        && points.columnCount() == sourceColumnCount);
        if (checker.checkPresent(!bAllColumns, "FirstColumn"))
        {
            checker.checkInteger(points.firstColumn(), "FirstColumn");
            checker.checkInteger(points.lastColumn(),  "LastColumn");
        }

        AnalysisDataValuesRef::const_iterator value;
        for (value = points.values().begin(); value != points.values().end(); ++value)
        {
            checkReferenceDataPoint(&checker, *value);
        }
    }
}


void
MockAnalysisDataModule::Impl::finishReferenceFrame(
        const AnalysisDataFrameHeader &header)
{
    EXPECT_TRUE(frameChecker_.get() != NULL);
    EXPECT_EQ(frameIndex_, header.index());
    frameChecker_.reset();
}


void
MockAnalysisDataModule::Impl::finishReferenceFrameSerial(int frameIndex)
{
    EXPECT_TRUE(frameChecker_.get() == NULL);
    EXPECT_EQ(frameIndex_, frameIndex);
    ++frameIndex_;
}


/********************************************************************
 * MockAnalysisDataModule
 */

namespace
{

/*! \brief
 * Helper function for checking the data frame header against static data.
 *
 * \param[in] header    Frame header to check.
 * \param[in] refFrame  Data to check against.
 */
void checkHeader(const AnalysisDataFrameHeader    &header,
                 const AnalysisDataTestInputFrame &refFrame)
{
    EXPECT_EQ(refFrame.index(), header.index());
    EXPECT_FLOAT_EQ(refFrame.x(), header.x());
    EXPECT_FLOAT_EQ(refFrame.dx(), header.dx());
}

/*! \brief
 * Helper function for checking a point set against static data.
 *
 * \param[in] points       Point set to check.
 * \param[in] refPoints    Data to check against.
 * \param[in] columnOffset Offset of first column of \p points in \p refPoints.
 */
void checkPoints(const AnalysisDataPointSetRef       &points,
                 const AnalysisDataTestInputPointSet &refPoints,
                 int                                  columnOffset)
{
    for (int i = 0; i < points.columnCount(); ++i)
    {
        const int column = points.firstColumn() - refPoints.firstColumn() + i + columnOffset;
        EXPECT_FLOAT_EQ(refPoints.y(column),
                        points.y(i))
        << "  Column: " << i+1 << " / " << points.columnCount()
        << " (+" << points.firstColumn() << ")\n"
        << "Ref. col: " << column+1 << " / " << refPoints.size()
        << " (+" << refPoints.firstColumn() << ", offs " << columnOffset << ")";
    }
}

/*! \brief
 * Helper function for checking a full frame against static data.
 *
 * \param[in] frame     Frame to check.
 * \param[in] refFrame  Data to check against.
 */
void checkFrame(const AnalysisDataFrameRef       &frame,
                const AnalysisDataTestInputFrame &refFrame)
{
    checkHeader(frame.header(), refFrame);
    ASSERT_EQ(refFrame.pointSetCount(), frame.pointSetCount());
    for (int i = 0; i < frame.pointSetCount(); ++i)
    {
        const AnalysisDataPointSetRef       &points    = frame.pointSet(i);
        const AnalysisDataTestInputPointSet &refPoints = refFrame.pointSet(i);
        EXPECT_EQ(refPoints.firstColumn(), points.firstColumn());
        checkPoints(points, refPoints, 0);
    }
}

/*! \brief
 * Functor for checking data frame header against static test input data.
 *
 * This functor is designed to be invoked as a handled for
 * AnalysisDataModuleInterface::frameStarted().
 */
class StaticDataFrameHeaderChecker
{
    public:
        /*! \brief
         * Constructs a checker against a given input data frame.
         *
         * \param[in] frame Frame to check against.
         *
         * \p frame must exist for the lifetime of this object.
         */
        StaticDataFrameHeaderChecker(const AnalysisDataTestInputFrame *frame)
            : frame_(frame)
        {
        }

        //! Function call operator for the functor.
        void operator()(const AnalysisDataFrameHeader &header) const
        {
            SCOPED_TRACE(formatString("Frame %d", frame_->index()));
            checkHeader(header, *frame_);
        }

    private:
        const AnalysisDataTestInputFrame *frame_;
};

/*! \brief
 * Functor for checking data frame points against static test input data.
 *
 * This functor is designed to be invoked as a handled for
 * AnalysisDataModuleInterface::pointsAdded().
 */
class StaticDataPointsChecker
{
    public:
        /*! \brief
         * Constructs a checker against a given input data frame and point set.
         *
         * \param[in] frame    Frame to check against.
         * \param[in] points   Point set in \p frame to check against.
         * \param[in] firstcol Expected first column.
         * \param[in] n        Expected number of columns.
         *
         * \p firstcol and \p n are used to create a checker that only expects
         * to be called for a subset of columns.
         * \p frame and \p points must exist for the lifetime of this object.
         */
        StaticDataPointsChecker(const AnalysisDataTestInputFrame *frame,
                                const AnalysisDataTestInputPointSet *points,
                                int firstcol, int n)
            : frame_(frame), points_(points), firstcol_(firstcol), n_(n)
        {
        }

        //! Function call operator for the functor.
        void operator()(const AnalysisDataPointSetRef &points) const
        {
            SCOPED_TRACE(formatString("Frame %d, point set %d",
                                      frame_->index(), points_->index()));
            EXPECT_EQ(points_->dataSetIndex(), points.dataSetIndex());
            const int expectedFirstColumn
                = std::max(0, points_->firstColumn() - firstcol_);
            const int expectedLastColumn
                = std::min(n_ - 1, points_->lastColumn() - firstcol_);
            EXPECT_EQ(expectedFirstColumn, points.firstColumn());
            EXPECT_EQ(expectedLastColumn,  points.lastColumn());
            checkHeader(points.header(), *frame_);
            checkPoints(points, *points_, firstcol_);
        }

    private:
        const AnalysisDataTestInputFrame    *frame_;
        const AnalysisDataTestInputPointSet *points_;
        int                                  firstcol_;
        int                                  n_;
};

/*! \brief
 * Functor for requesting data storage.
 *
 * This functor is designed to be invoked as a handled for
 * AnalysisDataModuleInterface::dataStarted().
 */
class DataStorageRequester
{
    public:
        /*! \brief
         * Constructs a functor that requests the given amount of storage.
         *
         * \param[in] count  Number of frames of storage to request, or
         *      -1 for all frames.
         *
         * \see AbstractAnalysisData::requestStorage()
         */
        explicit DataStorageRequester(int count) : count_(count) {}

        //! Function call operator for the functor.
        void operator()(AbstractAnalysisData *data) const
        {
            EXPECT_TRUE(data->requestStorage(count_));
        }

    private:
        int                     count_;
};

/*! \brief
 * Functor for checking data frame points and storage against static test input
 * data.
 *
 * This functor is designed to be invoked as a handled for
 * AnalysisDataModuleInterface::pointsAdded().
 */
class StaticDataPointsStorageChecker
{
    public:
        /*! \brief
         * Constructs a checker for a given frame.
         *
         * \param[in] source     Data object that is being checked.
         * \param[in] data       Test input data to check against.
         * \param[in] frameIndex Frame index for which this functor expects
         *      to be called.
         * \param[in] pointSetIndex Point set for which this functor expects
         *      to be called.
         * \param[in] storageCount How many past frames should be checked for
         *      storage (-1 = check all frames).
         *
         * This checker works as StaticDataPointsChecker, but additionally
         * checks that previous frames can be accessed using access methods
         * in AbstractAnalysisData and that correct data is returned.
         *
         * \p source and \p data must exist for the lifetime of this object.
         */
        StaticDataPointsStorageChecker(AbstractAnalysisData        *source,
                                       const AnalysisDataTestInput *data,
                                       int frameIndex, int pointSetIndex,
                                       int storageCount)
            : source_(source), data_(data),
              frameIndex_(frameIndex), pointSetIndex_(pointSetIndex),
              storageCount_(storageCount)
        {
        }

        //! Function call operator for the functor.
        void operator()(const AnalysisDataPointSetRef &points) const
        {
            SCOPED_TRACE(formatString("Frame %d", frameIndex_));
            const AnalysisDataTestInputFrame    &refFrame  = data_->frame(frameIndex_);
            const AnalysisDataTestInputPointSet &refPoints = refFrame.pointSet(pointSetIndex_);
            EXPECT_EQ(refPoints.firstColumn(), points.firstColumn());
            EXPECT_EQ(refPoints.size(), points.columnCount());
            checkHeader(points.header(), refFrame);
            checkPoints(points, refPoints, 0);
            for (int past = 1;
                 (storageCount_ < 0 || past <= storageCount_) && past <= frameIndex_;
                 ++past)
            {
                int   index = frameIndex_ - past;
                SCOPED_TRACE(formatString("Checking storage of frame %d", index));
                ASSERT_NO_THROW_GMX({
                                        AnalysisDataFrameRef frame = source_->getDataFrame(index);
                                        ASSERT_TRUE(frame.isValid());
                                        checkFrame(frame, data_->frame(index));
                                    });
            }
        }

    private:
        AbstractAnalysisData        *source_;
        const AnalysisDataTestInput *data_;
        int                          frameIndex_;
        int                          pointSetIndex_;
        int                          storageCount_;
};

/*! \brief
 * Sets the mock object expectation to mimick AnalysisDataModuleSerial.
 *
 * Makes MockAnalysisDataModule::parallelDataStarted() behave as if the mock
 * object was an AnalysisDataModuleSerial object: forward the call to
 * MockAnalysisDataModule::dataStarted() and return false.
 */
void setSerialExpectationForParallelDataStarted(MockAnalysisDataModule *mock)
{
    using ::testing::_;
    using ::testing::AtMost;
    using ::testing::DoAll;
    using ::testing::Invoke;
    using ::testing::Return;
    using ::testing::WithArg;
    EXPECT_CALL(*mock, parallelDataStarted(_, _))
        .Times(AtMost(1))
        .WillOnce(DoAll(WithArg<0>(Invoke(mock, &MockAnalysisDataModule::dataStarted)),
                        Return(false)));
}

}       // anonymous namespace


MockAnalysisDataModule::MockAnalysisDataModule(int flags)
    : impl_(new Impl(flags))
{
}


MockAnalysisDataModule::~MockAnalysisDataModule()
{
}


int MockAnalysisDataModule::flags() const
{
    return impl_->flags_;
}


void
MockAnalysisDataModule::setupStaticCheck(const AnalysisDataTestInput &data,
                                         AbstractAnalysisData        *source,
                                         bool                         bParallel)
{
    impl_->flags_ |= efAllowMulticolumn | efAllowMultipoint | efAllowMultipleDataSets;

    using ::testing::_;
    using ::testing::Invoke;
    using ::testing::Property;
    using ::testing::Return;

    if (bParallel)
    {
        ::testing::Expectation init =
            EXPECT_CALL(*this, parallelDataStarted(source, _))
                .WillOnce(Return(true));
        ::testing::ExpectationSet framesFinished;
        ::testing::Expectation    prevFinish;
        for (int row = 0; row < data.frameCount(); ++row)
        {
            ::testing::InSequence frameSequence;
            const AnalysisDataTestInputFrame &frame = data.frame(row);
            EXPECT_CALL(*this, frameStarted(Property(&AnalysisDataFrameHeader::index, row)))
                .After(init)
                .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame)));
            for (int ps = 0; ps < frame.pointSetCount(); ++ps)
            {
                const AnalysisDataTestInputPointSet &points = frame.pointSet(ps);
                StaticDataPointsChecker              checker(&frame, &points, 0,
                                                             data.columnCount(points.dataSetIndex()));
                EXPECT_CALL(*this, pointsAdded(Property(&AnalysisDataPointSetRef::frameIndex, row)))
                    .WillOnce(Invoke(checker));
            }
            EXPECT_CALL(*this, frameFinished(Property(&AnalysisDataFrameHeader::index, row)))
                .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame)));
            ::testing::Expectation finish;
            if (row > 0)
            {
                finish = EXPECT_CALL(*this, frameFinishedSerial(row))
                        .After(prevFinish);
            }
            else
            {
                finish = EXPECT_CALL(*this, frameFinishedSerial(row));
            }
            framesFinished += finish;
            prevFinish      = finish;
        }
        EXPECT_CALL(*this, dataFinished())
            .After(framesFinished);
    }
    else
    {
        ::testing::InSequence dummy;
        setSerialExpectationForParallelDataStarted(this);
        EXPECT_CALL(*this, dataStarted(source));
        for (int row = 0; row < data.frameCount(); ++row)
        {
            const AnalysisDataTestInputFrame &frame = data.frame(row);
            EXPECT_CALL(*this, frameStarted(_))
                .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame)));
            for (int ps = 0; ps < frame.pointSetCount(); ++ps)
            {
                const AnalysisDataTestInputPointSet &points = frame.pointSet(ps);
                StaticDataPointsChecker              checker(&frame, &points, 0,
                                                             data.columnCount(points.dataSetIndex()));
                EXPECT_CALL(*this, pointsAdded(_)).WillOnce(Invoke(checker));
            }
            EXPECT_CALL(*this, frameFinished(_))
                .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame)));
            EXPECT_CALL(*this, frameFinishedSerial(row));
        }
        EXPECT_CALL(*this, dataFinished());
    }
}


void
MockAnalysisDataModule::setupStaticColumnCheck(
        const AnalysisDataTestInput &data,
        int firstcol, int n, AbstractAnalysisData * /*source*/)
{
    impl_->flags_ |= efAllowMulticolumn | efAllowMultipoint | efAllowMultipleDataSets;

    ::testing::InSequence dummy;
    using ::testing::_;
    using ::testing::Invoke;

    setSerialExpectationForParallelDataStarted(this);
    EXPECT_CALL(*this, dataStarted(_));
    for (int row = 0; row < data.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame &frame = data.frame(row);
        EXPECT_CALL(*this, frameStarted(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame)));
        for (int ps = 0; ps < frame.pointSetCount(); ++ps)
        {
            const AnalysisDataTestInputPointSet &points = frame.pointSet(ps);
            if (points.lastColumn() >= firstcol
                && points.firstColumn() <= firstcol + n - 1)
            {
                EXPECT_CALL(*this, pointsAdded(_))
                    .WillOnce(Invoke(StaticDataPointsChecker(&frame, &points, firstcol, n)));
            }
        }
        EXPECT_CALL(*this, frameFinished(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame)));
        EXPECT_CALL(*this, frameFinishedSerial(row));
    }
    EXPECT_CALL(*this, dataFinished());
}


void
MockAnalysisDataModule::setupStaticStorageCheck(
        const AnalysisDataTestInput &data,
        int storageCount, AbstractAnalysisData *source)
{
    GMX_RELEASE_ASSERT(data.isMultipoint() == source->isMultipoint(),
                       "Mismatching multipoint properties");
    impl_->flags_ |= efAllowMulticolumn | efAllowMultipoint | efAllowMultipleDataSets;

    ::testing::InSequence dummy;
    using ::testing::_;
    using ::testing::Invoke;

    setSerialExpectationForParallelDataStarted(this);
    EXPECT_CALL(*this, dataStarted(source))
        .WillOnce(Invoke(DataStorageRequester(storageCount)));
    for (int row = 0; row < data.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame &frame = data.frame(row);
        EXPECT_CALL(*this, frameStarted(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame)));
        for (int pointSet = 0; pointSet < frame.pointSetCount(); ++pointSet)
        {
            StaticDataPointsStorageChecker checker(source, &data, row, pointSet,
                                                   storageCount);
            EXPECT_CALL(*this, pointsAdded(_)).WillOnce(Invoke(checker));
        }
        EXPECT_CALL(*this, frameFinished(_))
            .WillOnce(Invoke(StaticDataFrameHeaderChecker(&frame)));
        EXPECT_CALL(*this, frameFinishedSerial(row));
    }
    EXPECT_CALL(*this, dataFinished());
}


void
MockAnalysisDataModule::setupReferenceCheck(const TestReferenceChecker &checker,
                                            AbstractAnalysisData       *source)
{
    impl_->flags_ |= efAllowMulticolumn | efAllowMultipoint | efAllowMissing
        | efAllowMultipleDataSets;

    impl_->rootChecker_.reset(new TestReferenceChecker(checker));
    // Google Mock does not support checking the order fully, because
    // the number of frames is not known.
    // Order of the frameStarted(), pointsAdded() and frameFinished()
    // calls is checked using Google Test assertions in the invoked methods.
    using ::testing::_;
    using ::testing::AnyNumber;
    using ::testing::Expectation;
    using ::testing::Invoke;

    setSerialExpectationForParallelDataStarted(this);
    Expectation dataStart = EXPECT_CALL(*this, dataStarted(source))
            .WillOnce(Invoke(impl_.get(), &Impl::startReferenceData));
    Expectation frameStart = EXPECT_CALL(*this, frameStarted(_))
            .After(dataStart)
            .WillRepeatedly(Invoke(impl_.get(), &Impl::startReferenceFrame));
    Expectation pointsAdd = EXPECT_CALL(*this, pointsAdded(_))
            .After(dataStart)
            .WillRepeatedly(Invoke(impl_.get(), &Impl::checkReferencePoints));
    Expectation frameFinish = EXPECT_CALL(*this, frameFinished(_))
            .After(dataStart)
            .WillRepeatedly(Invoke(impl_.get(), &Impl::finishReferenceFrame));
    Expectation frameFinishSerial = EXPECT_CALL(*this, frameFinishedSerial(_))
            .After(dataStart)
            .WillRepeatedly(Invoke(impl_.get(), &Impl::finishReferenceFrameSerial));
    EXPECT_CALL(*this, dataFinished())
        .After(frameStart, pointsAdd, frameFinish, frameFinishSerial);
}

} // namespace test
} // namespace gmx
