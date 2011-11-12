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
 * Tests for functionality of analysis data histogram modules.
 *
 * These tests check that classes in histogram.h compute histograms correctly
 * with simple input data.  Also different ways of initializing the histograms
 * are tested.
 * Checking is done using gmx::test::AnalysisDataTestFixture and reference
 * data.  Also the input data is written to the reference data to catch
 * out-of-date reference.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/histogram.h"

#include "datatest.h"

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

typedef gmx::test::AnalysisDataTestFixture SimpleHistogramModuleTest;

// Input data for the tests below.
using gmx::test::END_OF_DATA;
using gmx::test::END_OF_FRAME;
using gmx::test::MPSTOP;
static const real simpleinputdata[] = {
    1.0,  0.7, MPSTOP, 1.1, MPSTOP, 2.3, MPSTOP, 2.9, END_OF_FRAME,
    2.0,  1.3, MPSTOP, 2.2, END_OF_FRAME,
    3.0,  3.3, MPSTOP, 1.2, MPSTOP, 1.3, END_OF_FRAME,
    END_OF_DATA
};

TEST_F(SimpleHistogramModuleTest, ComputesCorrectly)
{
    gmx::test::AnalysisDataTestInput input(simpleinputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount(), true);
    gmx::AnalysisDataSimpleHistogramModule *module =
        new gmx::AnalysisDataSimpleHistogramModule(
                gmx::histogramFromRange(1.0, 3.0).binCount(4));
    data.addModule(module);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("Histogram", module));
    ASSERT_NO_THROW(addReferenceCheckerModule("HistogramAverage",
                                              module->averager()));
    ASSERT_NO_THROW(presentAllData(input, &data));
    ASSERT_NO_THROW(module->averager()->done());
}


TEST_F(SimpleHistogramModuleTest, ComputesCorrectlyWithAll)
{
    gmx::test::AnalysisDataTestInput input(simpleinputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount(), true);
    gmx::AnalysisDataSimpleHistogramModule *module =
        new gmx::AnalysisDataSimpleHistogramModule(
                gmx::histogramFromRange(1.0, 3.0).binCount(4).includeAll());
    data.addModule(module);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("Histogram", module));
    ASSERT_NO_THROW(addReferenceCheckerModule("HistogramAverage",
                                              module->averager()));
    ASSERT_NO_THROW(presentAllData(input, &data));
    ASSERT_NO_THROW(module->averager()->done());
}


/********************************************************************
 * Tests for gmx::AnalysisDataWeightedHistogramModule.
 */

typedef gmx::test::AnalysisDataTestFixture WeightedHistogramModuleTest;

// Input data for the tests below (both weighted and bin average modules).
static const real weightedinputdata[] = {
    1.0,  0.7, 0.5, MPSTOP, 1.1, 1.0, MPSTOP, 2.3, 1.0, MPSTOP, 2.9, 2.0, END_OF_FRAME,
    2.0,  1.3, 1.0, MPSTOP, 2.2, 3.0, END_OF_FRAME,
    3.0,  3.3, 0.5, MPSTOP, 1.2, 2.0, MPSTOP, 1.3, 1.0, END_OF_FRAME,
    END_OF_DATA
};

TEST_F(WeightedHistogramModuleTest, ComputesCorrectly)
{
    gmx::test::AnalysisDataTestInput input(weightedinputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount(), true);
    gmx::AnalysisDataWeightedHistogramModule *module =
        new gmx::AnalysisDataWeightedHistogramModule(
                gmx::histogramFromRange(1.0, 3.0).binCount(4));
    data.addModule(module);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("Histogram", module));
    ASSERT_NO_THROW(addReferenceCheckerModule("HistogramAverage",
                                              module->averager()));
    ASSERT_NO_THROW(presentAllData(input, &data));
    ASSERT_NO_THROW(module->averager()->done());
}


TEST_F(WeightedHistogramModuleTest, ComputesCorrectlyWithAll)
{
    gmx::test::AnalysisDataTestInput input(weightedinputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount(), true);
    gmx::AnalysisDataWeightedHistogramModule *module =
        new gmx::AnalysisDataWeightedHistogramModule(
                gmx::histogramFromRange(1.0, 3.0).binCount(4).includeAll());
    data.addModule(module);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("Histogram", module));
    ASSERT_NO_THROW(addReferenceCheckerModule("HistogramAverage",
                                              module->averager()));
    ASSERT_NO_THROW(presentAllData(input, &data));
    ASSERT_NO_THROW(module->averager()->done());
}


/********************************************************************
 * Tests for gmx::AnalysisDataBinAverageModule.
 */

typedef gmx::test::AnalysisDataTestFixture BinAverageModuleTest;

TEST_F(BinAverageModuleTest, ComputesCorrectly)
{
    gmx::test::AnalysisDataTestInput input(weightedinputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount(), true);
    gmx::AnalysisDataBinAverageModule *module =
        new gmx::AnalysisDataBinAverageModule(
                gmx::histogramFromRange(1.0, 3.0).binCount(4));
    data.addModule(module);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("HistogramAverage", module));
    ASSERT_NO_THROW(presentAllData(input, &data));
}


TEST_F(BinAverageModuleTest, ComputesCorrectlyWithAll)
{
    gmx::test::AnalysisDataTestInput input(weightedinputdata);
    gmx::AnalysisData data;
    data.setColumns(input.columnCount(), true);
    gmx::AnalysisDataBinAverageModule *module =
        new gmx::AnalysisDataBinAverageModule(
                gmx::histogramFromRange(1.0, 3.0).binCount(4).includeAll());
    data.addModule(module);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("HistogramAverage", module));
    ASSERT_NO_THROW(presentAllData(input, &data));
}

} // namespace
