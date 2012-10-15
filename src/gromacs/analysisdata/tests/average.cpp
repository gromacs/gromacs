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
 * Tests for functionality of gmx::AnalysisDataAverageModule
 *
 * These tests check that gmx::AnalysisDataAverageModule computes averages
 * correctly with simple input data.
 * Checking is done using gmx::test::AnalysisDataTestFixture and reference
 * data.  Also the input data is written to the reference data to catch
 * out-of-date reference.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include <gtest/gtest.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"

#include "testutils/datatest.h"

namespace
{

/********************************************************************
 * Tests for gmx::AnalysisDataAverageModule.
 */

typedef gmx::test::AnalysisDataTestFixture AverageModuleTest;

// Input data for the tests below.
using gmx::test::END_OF_FRAME;
static const real inputdata[] = {
    1.0,  0.0, 1.0, 2.0, END_OF_FRAME,
    2.0,  1.0, 1.0, 1.0, END_OF_FRAME,
    3.0,  2.0, 0.0, 0.0, END_OF_FRAME
};

TEST_F(AverageModuleTest, BasicTest)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumnCount(input.columnCount());
    gmx::AnalysisDataAverageModulePointer module(
            new gmx::AnalysisDataAverageModule);
    data.addModule(module);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("Average", module.get()));
    ASSERT_NO_THROW(presentAllData(input, &data));
}

TEST_F(AverageModuleTest, CanCustomizeXAxis)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisData data;
    data.setColumnCount(input.columnCount());
    gmx::AnalysisDataAverageModulePointer module(new gmx::AnalysisDataAverageModule());
    data.addModule(module);
    module->setXAxis(0.5, 0.5);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("InputData", &data));
    ASSERT_NO_THROW(addReferenceCheckerModule("Average", module.get()));
    ASSERT_NO_THROW(presentAllData(input, &data));
}

} // namespace
