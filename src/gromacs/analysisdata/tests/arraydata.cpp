/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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
#include <gtest/gtest.h>

#include "gromacs/analysisdata/arraydata.h"

#include "testutils/datatest.h"

namespace
{

/********************************************************************
 * Tests for gmx::AnalysisArrayData.
 */

//! Test fixture for gmx::AnalysisArrayData.
typedef gmx::test::AnalysisDataTestFixture AnalysisArrayDataTest;

using gmx::test::END_OF_FRAME;
//! Input data for gmx::AnalysisArrayData tests.
const real inputdata[] = {
    1.0,  0.0, 1.0, 2.0, END_OF_FRAME,
    2.0,  1.0, 1.0, 1.0, END_OF_FRAME,
    3.0,  2.0, 0.0, 0.0, END_OF_FRAME,
    4.0,  3.0, 2.0, 1.0, END_OF_FRAME
};

TEST_F(AnalysisArrayDataTest, CallsModuleCorrectly)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisArrayData           data;
    data.setXAxis(1.0, 1.0);
    setupArrayData(input, &data);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(data.valuesReady());
}

TEST_F(AnalysisArrayDataTest, StorageWorks)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisArrayData           data;
    data.setXAxis(1.0, 1.0);
    setupArrayData(input, &data);

    ASSERT_NO_THROW(addStaticStorageCheckerModule(input, -1, &data));
    ASSERT_NO_THROW(data.valuesReady());
}

} // namespace
