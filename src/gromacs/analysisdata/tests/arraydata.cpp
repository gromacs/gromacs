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
 * Tests for gmx::AnalysisArrayData functionality.
 *
 * These tests check the functionality of gmx::AnalysisArrayData and its base
 * class gmx::AbstractAnalysisArrayData.
 * Checking is done using gmx::test::AnalysisDataTestFixture and mock
 * modules that implement gmx::AnalysisDataModuleInterface.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
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
    gmx::AnalysisArrayData data;
    data.setXAxis(1.0, 1.0);
    setupArrayData(input, &data);

    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(addStaticCheckerModule(input, &data));
    ASSERT_NO_THROW(data.valuesReady());
}

TEST_F(AnalysisArrayDataTest, StorageWorks)
{
    gmx::test::AnalysisDataTestInput input(inputdata);
    gmx::AnalysisArrayData data;
    data.setXAxis(1.0, 1.0);
    setupArrayData(input, &data);

    ASSERT_NO_THROW(addStaticStorageCheckerModule(input, -1, &data));
    ASSERT_NO_THROW(data.valuesReady());
}

} // namespace
