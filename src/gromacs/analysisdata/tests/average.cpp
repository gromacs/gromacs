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
 * Tests for functionality of analysis data averaging modules.
 *
 * These tests check that gmx::AnalysisDataAverageModule and
 * gmx::AnalysisDataFrameAverageModule compute averages correctly with simple
 * input data.
 * Checking is done using gmx::test::AnalysisDataTestFixture and reference
 * data.  Also the input data is written to the reference data to catch
 * out-of-date reference.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/modules/average.h"

#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"

#include "gromacs/analysisdata/tests/datatest.h"
#include "testutils/testasserts.h"

using gmx::test::AnalysisDataTestInput;

namespace
{

// Simple input data for gmx::AnalysisDataAverageModule tests.
class SimpleInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
#ifndef STATIC_ANON_NAMESPACE_BUG
            static SimpleInputData singleton;
            return singleton.data_;
#else
            static SimpleInputData singleton_average;
            return singleton_average.data_;
#endif
        }

        SimpleInputData() : data_(1, false)
        {
            data_.setColumnCount(0, 3);
            data_.addFrameWithValues(1.0,  0.0, 1.0, 2.0);
            data_.addFrameWithValues(2.0,  1.0, 1.0, 1.0);
            data_.addFrameWithValues(3.0,  2.0, 0.0, 0.0);
        }

    private:
        AnalysisDataTestInput  data_;
};

// Multipoint input data for gmx::AnalysisDataAverageModule tests.
class MultipointInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
#ifndef STATIC_ANON_NAMESPACE_BUG
            static MultipointInputData singleton;
            return singleton.data_;
#else
            static MultipointInputData singleton_average;
            return singleton_average.data_;
#endif
        }

        MultipointInputData() : data_(1, true)
        {
            using gmx::test::AnalysisDataTestInputFrame;
            data_.setColumnCount(0, 3);
            AnalysisDataTestInputFrame &frame1 = data_.addFrame(1.0);
            frame1.addPointSetWithValues(0, 0, 0.0, 1.0, 2.0);
            frame1.addPointSetWithValues(0, 0, 1.0, 0.0);
            frame1.addPointSetWithValues(0, 0, 2.0);
            AnalysisDataTestInputFrame &frame2 = data_.addFrame(2.0);
            frame2.addPointSetWithValues(0, 0, 1.0, 1.0);
            frame2.addPointSetWithValues(0, 0, 2.0);
            AnalysisDataTestInputFrame &frame3 = data_.addFrame(3.0);
            frame3.addPointSetWithValues(0, 0, 2.0, 0.0, 0.0);
        }

    private:
        AnalysisDataTestInput  data_;
};

// Input data with multiple data sets for gmx::AnalysisDataAverageModule tests.
class MultiDataSetInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
#ifndef STATIC_ANON_NAMESPACE_BUG
            static MultiDataSetInputData singleton;
            return singleton.data_;
#else
            static MultiDataSetInputData singleton_average;
            return singleton_average.data_;
#endif
        }

        MultiDataSetInputData() : data_(2, true)
        {
            using gmx::test::AnalysisDataTestInputFrame;
            data_.setColumnCount(0, 3);
            data_.setColumnCount(1, 2);
            AnalysisDataTestInputFrame &frame1 = data_.addFrame(1.0);
            frame1.addPointSetWithValues(0, 0, 0.0, 1.0, 2.0);
            frame1.addPointSetWithValues(0, 0, 1.0, 0.0);
            frame1.addPointSetWithValues(1, 0, 2.0, 1.0);
            frame1.addPointSetWithValues(1, 1, 2.0);
            AnalysisDataTestInputFrame &frame2 = data_.addFrame(2.0);
            frame2.addPointSetWithValues(0, 0, 1.0, 1.0);
            frame2.addPointSetWithValues(0, 2, 2.0);
            frame2.addPointSetWithValues(1, 0, 1.0, 0.0);
            AnalysisDataTestInputFrame &frame3 = data_.addFrame(3.0);
            frame3.addPointSetWithValues(0, 0, 2.0, 0.0, 0.0);
            frame3.addPointSetWithValues(1, 0, 0.0, 2.0);
        }

    private:
        AnalysisDataTestInput  data_;
};


/********************************************************************
 * Tests for gmx::AnalysisDataAverageModule.
 */

//! Test fixture for gmx::AnalysisDataAverageModule.
typedef gmx::test::AnalysisDataTestFixture AverageModuleTest;

TEST_F(AverageModuleTest, BasicTest)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataAverageModulePointer module(
            new gmx::AnalysisDataAverageModule);
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Average", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

TEST_F(AverageModuleTest, HandlesMultipointData)
{
    const AnalysisDataTestInput &input = MultipointInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataAverageModulePointer module(
            new gmx::AnalysisDataAverageModule);
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Average", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

TEST_F(AverageModuleTest, HandlesMultipleDataSets)
{
    const AnalysisDataTestInput &input = MultiDataSetInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataAverageModulePointer module(
            new gmx::AnalysisDataAverageModule);
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Average", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

TEST_F(AverageModuleTest, HandlesDataSetAveraging)
{
    const AnalysisDataTestInput &input = MultiDataSetInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataAverageModulePointer module(
            new gmx::AnalysisDataAverageModule);
    module->setAverageDataSets(true);
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Average", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

TEST_F(AverageModuleTest, CanCustomizeXAxis)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataAverageModulePointer module(new gmx::AnalysisDataAverageModule());
    data.addModule(module);
    module->setXAxis(0.5, 0.5);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Average", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

TEST_F(AverageModuleTest, CanCustomizeNonUniformXAxis)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataAverageModulePointer module(new gmx::AnalysisDataAverageModule());
    data.addModule(module);
    module->setXAxisValue(0, 2.0);
    module->setXAxisValue(1, 3.0);
    module->setXAxisValue(2, 5.0);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("Average", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

/********************************************************************
 * Tests for gmx::AnalysisDataFrameAverageModule.
 */

//! Test fixture for gmx::AnalysisDataFrameAverageModule.
typedef gmx::test::AnalysisDataTestFixture FrameAverageModuleTest;

TEST_F(FrameAverageModuleTest, BasicTest)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataFrameAverageModulePointer module(
            new gmx::AnalysisDataFrameAverageModule);
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("FrameAverage", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

TEST_F(FrameAverageModuleTest, HandlesMultipleDataSets)
{
    const AnalysisDataTestInput &input = MultiDataSetInputData::get();
    gmx::AnalysisData            data;
    ASSERT_NO_THROW_GMX(setupDataObject(input, &data));

    gmx::AnalysisDataFrameAverageModulePointer module(
            new gmx::AnalysisDataFrameAverageModule);
    data.addModule(module);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW_GMX(addReferenceCheckerModule("FrameAverage", module.get()));
    ASSERT_NO_THROW_GMX(presentAllData(input, &data));
}

} // namespace
