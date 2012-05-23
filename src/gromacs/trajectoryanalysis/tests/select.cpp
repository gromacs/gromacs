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
 * Tests for functionality of the "select" trajectory analysis module.
 *
 * These tests test most of the functionality of the module, but currently
 * missing are:
 *  - Tests related to -oc output.  This would require a more complex input
 *    structure for reasonable testing (the same structure could also be used
 *    in selection unit tests for 'insolidangle' keyword).
 *  - Tests for XVG labels.  This is a limitation of the current testing
 *    framework.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#include <gtest/gtest.h>

#include "gromacs/trajectoryanalysis/modules/select.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;

/********************************************************************
 * Tests for gmx::analysismodules::Select.
 */

typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::Select>
        SelectModuleTest;

TEST_F(SelectModuleTest, BasicTest)
{
    const char *const cmdline[] = {
        "select",
        "-select", "y < 2.5", "resname RA"
    };
    setTopology("simple.gro");
    setOutputFile("-oi", "index.dat");
    setOutputFile("-on", "index.ndx");
    excludeDataset("cfrac");
    runTest(CommandLine::create(cmdline));
}

TEST_F(SelectModuleTest, HandlesDumpOption)
{
    const char *const cmdline[] = {
        "select",
        "-select", "y < 2.5",
        "-dump"
    };
    setTopology("simple.gro");
    setOutputFile("-oi", "index.dat");
    includeDataset("index");
    runTest(CommandLine::create(cmdline));
}

TEST_F(SelectModuleTest, NormalizesSizes)
{
    const char *const cmdline[] = {
        "select",
        "-select", "y < 2.5", "resname RA and y < 2.5", "resname RA",
        "-norm"
    };
    setTopology("simple.gro");
    includeDataset("size");
    runTest(CommandLine::create(cmdline));
}

TEST_F(SelectModuleTest, WritesResidueNumbers)
{
    const char *const cmdline[] = {
        "select",
        "-select", "res_com of resname RA RD"
    };
    setTopology("simple.gro");
    includeDataset("index");
    runTest(CommandLine::create(cmdline));
}

TEST_F(SelectModuleTest, WritesResidueIndices)
{
    const char *const cmdline[] = {
        "select",
        "-select", "res_com of resname RA RD",
        "-resnr", "index"
    };
    setTopology("simple.gro");
    includeDataset("index");
    runTest(CommandLine::create(cmdline));
}

} // namespace
