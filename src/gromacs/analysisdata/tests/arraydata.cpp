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
 * Tests for gmx::AnalysisArrayData functionality.
 *
 * These tests check the functionality of gmx::AnalysisArrayData and its base
 * class gmx::AbstractAnalysisArrayData.
 * Checking is done using gmx::test::AnalysisDataTestFixture and mock
 * modules that implement gmx::AnalysisDataModuleInterface.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/arraydata.h"

#include <gtest/gtest.h>

#include "gromacs/analysisdata/tests/datatest.h"
#include "testutils/testasserts.h"

using gmx::test::AnalysisDataTestInput;

namespace
{

/********************************************************************
 * Tests for gmx::AnalysisArrayData.
 */

//! Test fixture for gmx::AnalysisArrayData.
typedef gmx::test::AnalysisDataTestFixture AnalysisArrayDataTest;

// Input data for gmx::AnalysisArrayData tests.
class SimpleInputData
{
    public:
        static const AnalysisDataTestInput &get()
        {
#ifndef STATIC_ANON_NAMESPACE_BUG
            static SimpleInputData singleton;
            return singleton.data_;
#else
            static SimpleInputData singleton_arraydata;
            return singleton_arraydata.data_;
#endif
        }

        SimpleInputData() : data_(1, false)
        {
            data_.setColumnCount(0, 3);
            data_.addFrameWithValues(1.0,  0.0, 1.0, 2.0);
            data_.addFrameWithValues(2.0,  1.0, 1.0, 1.0);
            data_.addFrameWithValues(3.0,  2.0, 0.0, 0.0);
            data_.addFrameWithValues(4.0,  3.0, 2.0, 1.0);
        }

    private:
        AnalysisDataTestInput  data_;
};

TEST_F(AnalysisArrayDataTest, CallsModuleCorrectly)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisArrayData       data;
    data.setXAxis(1.0, 1.0);
    setupArrayData(input, &data);

    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW_GMX(data.valuesReady());
}

TEST_F(AnalysisArrayDataTest, StorageWorks)
{
    const AnalysisDataTestInput &input = SimpleInputData::get();
    gmx::AnalysisArrayData       data;
    data.setXAxis(1.0, 1.0);
    setupArrayData(input, &data);

    ASSERT_NO_THROW_GMX(addStaticStorageCheckerModule(input, -1, &data));
    ASSERT_NO_THROW_GMX(data.valuesReady());
}

TEST_F(AnalysisArrayDataTest, CanSetXAxis)
{
    gmx::AnalysisArrayData       data;
    data.setRowCount(5);
    data.setXAxis(1.0, 1.0);
    EXPECT_FLOAT_EQ(1.0, data.xvalue(0));
    EXPECT_FLOAT_EQ(3.0, data.xvalue(2));
    EXPECT_FLOAT_EQ(5.0, data.xvalue(4));
    data.setXAxisValue(0, 3.0);
    data.setXAxisValue(2, 1.0);
    EXPECT_FLOAT_EQ(3.0, data.xvalue(0));
    EXPECT_FLOAT_EQ(2.0, data.xvalue(1));
    EXPECT_FLOAT_EQ(1.0, data.xvalue(2));
    EXPECT_FLOAT_EQ(4.0, data.xvalue(3));
}

TEST_F(AnalysisArrayDataTest, CanSetXAxisBeforeRowCount)
{
    {
        gmx::AnalysisArrayData       data;
        data.setXAxis(1.0, 1.0);
        data.setRowCount(5);
        EXPECT_FLOAT_EQ(1.0, data.xvalue(0));
        EXPECT_FLOAT_EQ(3.0, data.xvalue(2));
        EXPECT_FLOAT_EQ(5.0, data.xvalue(4));
    }
    {
        gmx::AnalysisArrayData       data;
        data.setXAxisValue(0, 2.0);
        data.setXAxisValue(1, 3.0);
        data.setXAxisValue(2, 5.0);
        data.setRowCount(3);
        EXPECT_FLOAT_EQ(2.0, data.xvalue(0));
        EXPECT_FLOAT_EQ(3.0, data.xvalue(1));
        EXPECT_FLOAT_EQ(5.0, data.xvalue(2));
    }
}

} // namespace
