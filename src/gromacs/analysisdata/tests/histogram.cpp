/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Tests for functionality of analysis data histogram modules.
 *
 * These tests check that classes in histogram.h compute histograms correctly
 * with simple input data.  Also different ways of initializing the histograms
 * are tested.
 * Checking is done using gmx::test::AnalysisDataTestFixture and reference
 * data.  Also the input data is written to the reference data to catch
 * out-of-date reference.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/modules/histogram.h"

#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"

#include "gromacs/analysisdata/tests/datatest.h"
#include "testutils/testasserts.h"

using gmx::test::AnalysisDataTestInput;

namespace
{

/********************************************************************
 * Tests for gmx::AnalysisHistogramSettings.
 *
 * These tests check that gmx::AnalysisHistogramSettings objects can be
 * initialized from various types of values, and that the bin positions are
 * computed correctly based on the input values.
 */

TEST(AnalysisHistogramSettingsTest, InitializesFromBins)
{
    gmx::AnalysisHistogramSettings settings(
            gmx::histogramFromBins(1.0, 5, 0.5));
    EXPECT_FLOAT_EQ(1.0, settings.firstEdge());
    EXPECT_EQ(5, settings.binCount());
    EXPECT_FLOAT_EQ(0.5, settings.binWidth());
    EXPECT_FLOAT_EQ(3.5, settings.lastEdge());
}


TEST(AnalysisHistogramSettingsTest, InitializesFromBinsWithIntegerBins)
{
    gmx::AnalysisHistogramSettings settings(
            gmx::histogramFromBins(1.0, 5, 0.5).integerBins());
    EXPECT_FLOAT_EQ(0.75, settings.firstEdge());
    EXPECT_EQ(5, settings.binCount());
    EXPECT_FLOAT_EQ(0.5, settings.binWidth());
    EXPECT_FLOAT_EQ(3.25, settings.lastEdge());
}


TEST(AnalysisHistogramSettingsTest, InitializesFromRangeWithBinCount)
{
    gmx::AnalysisHistogramSettings settings(
            gmx::histogramFromRange(1.0, 4.0).binCount(6));
    EXPECT_FLOAT_EQ(1.0, settings.firstEdge());
    EXPECT_FLOAT_EQ(4.0, settings.lastEdge());
    EXPECT_EQ(6, settings.binCount());
    EXPECT_FLOAT_EQ(0.5, settings.binWidth());
}


TEST(AnalysisHistogramSettingsTest, InitializesFromRangeWithBinWidth)
{
    gmx::AnalysisHistogramSettings settings(
            gmx::histogramFromRange(1.0, 4.0).binWidth(0.5));
    EXPECT_FLOAT_EQ(1.0, settings.firstEdge());
    EXPECT_FLOAT_EQ(4.0, settings.lastEdge());
    EXPECT_FLOAT_EQ(0.5, settings.binWidth());
    EXPECT_EQ(6, settings.binCount());
}


TEST(AnalysisHistogramSettingsTest, InitializesFromRangeWithBinCountAndIntegerBins)
{
    gmx::AnalysisHistogramSettings settings(
            gmx::histogramFromRange(1.0, 4.0).binCount(7).integerBins());
    EXPECT_FLOAT_EQ(0.75, settings.firstEdge());
    EXPECT_FLOAT_EQ(4.25, settings.lastEdge());
    EXPECT_EQ(7, settings.binCount());
    EXPECT_FLOAT_EQ(0.5, settings.binWidth());
}


TEST(AnalysisHistogramSettingsTest, InitializesFromRangeWithBinWidthAndIntegerBins)
{
    gmx::AnalysisHistogramSettings settings(
            gmx::histogramFromRange(1.0, 4.0).binWidth(0.5).integerBins());
    EXPECT_FLOAT_EQ(0.75, settings.firstEdge());
    EXPECT_FLOAT_EQ(4.25, settings.lastEdge());
    EXPECT_FLOAT_EQ(0.5, settings.binWidth());
    EXPECT_EQ(7, settings.binCount());
}


TEST(AnalysisHistogramSettingsTest, InitializesFromRangeWithRoundedRange)
{
    gmx::AnalysisHistogramSettings settings(
            gmx::histogramFromRange(1.2, 3.8).binWidth(0.5).roundRange());
    EXPECT_FLOAT_EQ(1.0, settings.firstEdge());
    EXPECT_FLOAT_EQ(4.0, settings.lastEdge());
    EXPECT_FLOAT_EQ(0.5, settings.binWidth());
    EXPECT_EQ(6, settings.binCount());
}


/********************************************************************
 * Tests for gmx::AnalysisDataSimpleHistogramModule.
 */

//! Test fixture for gmx::AnalysisDataSimpleHistogramModule.
typedef gmx::test::AnalysisDataTestFixture SimpleHistogramModuleTest;

// Input data for gmx::AnalysisDataSimpleHistogramModule tests.
class SimpleInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
#ifndef STATIC_ANON_NAMESPACE_BUG
            static SimpleInputData singleton;
            return singleton.data_;
#else
            static SimpleInputData singleton_histogram;
            return singleton_histogram.data_;
#endif
        }

        SimpleInputData() : data_(1, true)
        {
            using gmx::test::AnalysisDataTestInputFrame;
            data_.setColumnCount(0, 1);
            AnalysisDataTestInputFrame &frame1 = data_.addFrame(1.0);
            frame1.addPointSetWithValues(0, 0, 0.7);
            frame1.addPointSetWithValues(0, 0, 1.1);
            frame1.addPointSetWithValues(0, 0, 2.3);
            frame1.addPointSetWithValues(0, 0, 2.9);
            AnalysisDataTestInputFrame &frame2 = data_.addFrame(2.0);
            frame2.addPointSetWithValues(0, 0, 1.3);
            frame2.addPointSetWithValues(0, 0, 2.2);
            AnalysisDataTestInputFrame &frame3 = data_.addFrame(3.0);
            frame3.addPointSetWithValues(0, 0, 3.3);
            frame3.addPointSetWithValues(0, 0, 1.2);
            frame3.addPointSetWithValues(0, 0, 1.3);
        }

    private:
        AnalysisDataTestInput  data_;
};

TEST_F(SimpleHistogramModuleTest, ComputesCorrectly)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataSimpleHistogramModulePointer module(
            new gmx::AnalysisDataSimpleHistogramModule(
                    gmx::histogramFromRange(1.0, 3.0).binCount(4)));
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Histogram", module.get()));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("HistogramAverage",
                                                  &module->averager()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
    ASSERT_NO_THROW_GMX(module->averager().done());
}


TEST_F(SimpleHistogramModuleTest, ComputesCorrectlyWithAll)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataSimpleHistogramModulePointer module(
            new gmx::AnalysisDataSimpleHistogramModule(
                    gmx::histogramFromRange(1.0, 3.0).binCount(4).includeAll()));
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Histogram", module.get()));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("HistogramAverage",
                                                  &module->averager()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
    ASSERT_NO_THROW_GMX(module->averager().done());
}


/********************************************************************
 * Tests for gmx::AnalysisDataWeightedHistogramModule.
 */

//! Test fixture for gmx::AnalysisDataWeightedHistogramModule.
typedef gmx::test::AnalysisDataTestFixture WeightedHistogramModuleTest;

// Input data for both weighted histogram and bin average module tests.
class WeightedSimpleInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
#ifndef STATIC_ANON_NAMESPACE_BUG
            static WeightedSimpleInputData singleton;
            return singleton.data_;
#else
            static WeightedSimpleInputData singleton_histogram;
            return singleton_histogram.data_;
#endif
        }

        WeightedSimpleInputData() : data_(1, true)
        {
            using gmx::test::AnalysisDataTestInputFrame;
            data_.setColumnCount(0, 2);
            AnalysisDataTestInputFrame &frame1 = data_.addFrame(1.0);
            frame1.addPointSetWithValues(0, 0, 0.7, 0.5);
            frame1.addPointSetWithValues(0, 0, 1.1, 1.0);
            frame1.addPointSetWithValues(0, 0, 2.3, 1.0);
            frame1.addPointSetWithValues(0, 0, 2.9, 2.0);
            AnalysisDataTestInputFrame &frame2 = data_.addFrame(2.0);
            frame2.addPointSetWithValues(0, 0, 1.3, 1.0);
            frame2.addPointSetWithValues(0, 0, 2.2, 3.0);
            AnalysisDataTestInputFrame &frame3 = data_.addFrame(3.0);
            frame3.addPointSetWithValues(0, 0, 3.3, 0.5);
            frame3.addPointSetWithValues(0, 0, 1.2, 2.0);
            frame3.addPointSetWithValues(0, 0, 1.3, 1.0);
        }

    private:
        AnalysisDataTestInput  data_;
};

// Input data for both weighted histogram and bin average module tests.
class WeightedDataSetInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
#ifndef STATIC_ANON_NAMESPACE_BUG
            static WeightedDataSetInputData singleton;
            return singleton.data_;
#else
            static WeightedDataSetInputData singleton_histogram;
            return singleton_histogram.data_;
#endif
        }

        WeightedDataSetInputData() : data_(2, true)
        {
            using gmx::test::AnalysisDataTestInputFrame;
            data_.setColumnCount(0, 2);
            data_.setColumnCount(1, 2);
            AnalysisDataTestInputFrame &frame1 = data_.addFrame(1.0);
            frame1.addPointSetWithValues(0, 0, 0.7, 0.5);
            frame1.addPointSetWithValues(0, 0, 1.1, 1.0);
            frame1.addPointSetWithValues(1, 0, 2.3, 1.0);
            frame1.addPointSetWithValues(1, 0, 2.9, 2.0);
            AnalysisDataTestInputFrame &frame2 = data_.addFrame(2.0);
            frame2.addPointSetWithValues(0, 0, 1.3, 1.0);
            frame2.addPointSetWithValues(1, 0, 2.2, 3.0);
            AnalysisDataTestInputFrame &frame3 = data_.addFrame(3.0);
            frame3.addPointSetWithValues(0, 0, 3.3, 0.5);
            frame3.addPointSetWithValues(0, 0, 1.2, 2.0);
            frame3.addPointSetWithValues(1, 0, 1.3, 1.0);
        }

    private:
        AnalysisDataTestInput  data_;
};

TEST_F(WeightedHistogramModuleTest, ComputesCorrectly)
{
    const AnalysisDataTestInput &input = WeightedSimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataWeightedHistogramModulePointer module(
            new gmx::AnalysisDataWeightedHistogramModule(
                    gmx::histogramFromRange(1.0, 3.0).binCount(4)));
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Histogram", module.get()));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("HistogramAverage",
                                                  &module->averager()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
    ASSERT_NO_THROW_GMX(module->averager().done());
}


TEST_F(WeightedHistogramModuleTest, ComputesCorrectlyWithAll)
{
    const AnalysisDataTestInput &input = WeightedSimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataWeightedHistogramModulePointer module(
            new gmx::AnalysisDataWeightedHistogramModule(
                    gmx::histogramFromRange(1.0, 3.0).binCount(4).includeAll()));
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Histogram", module.get()));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("HistogramAverage",
                                                  &module->averager()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
    ASSERT_NO_THROW_GMX(module->averager().done());
}


TEST_F(WeightedHistogramModuleTest, HandlesMultipleDataSets)
{
    const AnalysisDataTestInput &input = WeightedDataSetInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataWeightedHistogramModulePointer module(
            new gmx::AnalysisDataWeightedHistogramModule(
                    gmx::histogramFromRange(1.0, 3.0).binCount(4)));
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Histogram", module.get()));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("HistogramAverage",
                                                  &module->averager()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
    ASSERT_NO_THROW_GMX(module->averager().done());
}


/********************************************************************
 * Tests for gmx::AnalysisDataBinAverageModule.
 */

//! Test fixture for gmx::AnalysisDataBinAverageModule.
typedef gmx::test::AnalysisDataTestFixture BinAverageModuleTest;

TEST_F(BinAverageModuleTest, ComputesCorrectly)
{
    const AnalysisDataTestInput &input = WeightedSimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataBinAverageModulePointer module(
            new gmx::AnalysisDataBinAverageModule(
                    gmx::histogramFromRange(1.0, 3.0).binCount(4)));
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("HistogramAverage", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}


TEST_F(BinAverageModuleTest, ComputesCorrectlyWithAll)
{
    const AnalysisDataTestInput &input = WeightedSimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataBinAverageModulePointer module(
            new gmx::AnalysisDataBinAverageModule(
                    gmx::histogramFromRange(1.0, 3.0).binCount(4).includeAll()));
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("HistogramAverage", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}


TEST_F(BinAverageModuleTest, HandlesMultipleDataSets)
{
    const AnalysisDataTestInput &input = WeightedDataSetInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataBinAverageModulePointer module(
            new gmx::AnalysisDataBinAverageModule(
                    gmx::histogramFromRange(1.0, 3.0).binCount(4)));
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("HistogramAverage", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}


/********************************************************************
 * Tests for gmx::AbstractAverageHistogram.
 *
 * This class derives from gmx::AbstractAnalysisArrayData, and is tested using
 * corresponding facilities in gmx::test::AnalysisDataTestFixture.
 */

//! Test fixture for gmx::AbstractAverageHistogram.
typedef gmx::test::AnalysisDataTestFixture AbstractAverageHistogramTest;

// Input data for gmx::AbstractAverageHistogram tests.
class AverageInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
#ifndef STATIC_ANON_NAMESPACE_BUG
            static AverageInputData singleton;
            return singleton.data_;
#else
            static AverageInputData singleton_histogram;
            return singleton_histogram.data_;
#endif
        }

        AverageInputData() : data_(1, false)
        {
            data_.setColumnCount(0, 1);
            data_.addFrameWithValueAndError(1.0,  2.0, 1.0);
            data_.addFrameWithValueAndError(1.5,  1.0, 1.0);
            data_.addFrameWithValueAndError(2.0,  3.0, 2.0);
            data_.addFrameWithValueAndError(2.5,  4.0, 2.0);
            data_.addFrameWithValueAndError(3.0,  2.0, 1.0);
            data_.addFrameWithValueAndError(3.5,  0.0, 3.0);
            data_.addFrameWithValueAndError(4.0,  1.0, 3.0);
        }

    private:
        AnalysisDataTestInput  data_;
};

/*! \brief
 * Mock object for testing gmx::AbstractAverageHistogram.
 *
 * Exposes necessary methods from gmx::AbstractAverageHistogram to use with
 * gmx::test::AnalysisDataTestFixture::setupArrayData().
 *
 * \ingroup module_analysisdata
 */
class MockAverageHistogram : public gmx::AbstractAverageHistogram
{
    public:
        //! Creates a histogram module with defined bin parameters.
        explicit MockAverageHistogram(const gmx::AnalysisHistogramSettings &settings)
            : AbstractAverageHistogram(settings)
        {
        }

        using AbstractAverageHistogram::init;
        using AbstractAverageHistogram::setColumnCount;
        using AbstractAverageHistogram::setRowCount;
        using AbstractAverageHistogram::allocateValues;
        using AbstractAverageHistogram::value;
};


TEST_F(AbstractAverageHistogramTest, ClonesCorrectly)
{
    const AnalysisDataTestInput &input = AverageInputData::get();
    MockAverageHistogram         data(
            gmx::histogramFromBins(1.0, input.frameCount(), 0.5).integerBins());
    setupArrayData(input, &data);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    gmx::AverageHistogramPointer copy(data.clone());
    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, copy.get()));
    ASSERT_NO_THROW_GMX(copy->done());
    ASSERT_NO_THROW_GMX(data.done());
    gmx::AverageHistogramPointer copy2(data.clone());
    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, copy2.get()));
    ASSERT_NO_THROW_GMX(copy2->done());
}


TEST_F(AbstractAverageHistogramTest, ComputesCumulativeHistogram)
{
    const AnalysisDataTestInput &input = AverageInputData::get();
    MockAverageHistogram         data(
            gmx::histogramFromBins(1.0, input.frameCount(), 0.5).integerBins());
    setupArrayData(input, &data);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(data.done());

    gmx::AverageHistogramPointer cumulative(data.clone());
    cumulative->makeCumulative();
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("CumulativeHistogram", cumulative.get()));
    ASSERT_NO_THROW_GMX(cumulative->done());
}


TEST_F(AbstractAverageHistogramTest, ResamplesAtDoubleBinWidth)
{
    const AnalysisDataTestInput &input = AverageInputData::get();
    MockAverageHistogram         data(
            gmx::histogramFromBins(1.0, input.frameCount(), 0.5).integerBins());
    setupArrayData(input, &data);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    gmx::AverageHistogramPointer resampled(data.resampleDoubleBinWidth(false));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("ResampledHistogram", resampled.get()));
    ASSERT_NO_THROW_GMX(data.done());
    ASSERT_NO_THROW_GMX(resampled->done());
}


TEST_F(AbstractAverageHistogramTest, ResamplesAtDoubleBinWidthWithIntegerBins)
{
    const AnalysisDataTestInput &input = AverageInputData::get();
    MockAverageHistogram         data(
            gmx::histogramFromBins(1.0, input.frameCount(), 0.5).integerBins());
    setupArrayData(input, &data);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    gmx::AverageHistogramPointer resampled(data.resampleDoubleBinWidth(true));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("ResampledHistogram", resampled.get()));
    ASSERT_NO_THROW_GMX(data.done());
    ASSERT_NO_THROW_GMX(resampled->done());
}

} // namespace
