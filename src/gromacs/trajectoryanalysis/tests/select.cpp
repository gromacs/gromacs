/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/select.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/cmdlinetest.h"
#include "testutils/textblockmatchers.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::ExactTextMatch;

/********************************************************************
 * Tests for gmx::analysismodules::Select.
 */

//! Test fixture for the select analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::SelectInfo> SelectModuleTest;

TEST_F(SelectModuleTest, BasicTest)
{
    const char* const cmdline[] = { "select", "-select", "y < 2.5", "resname RA" };
    setTopology("simple.gro");
    setTrajectory("simple.gro");
    setOutputFile("-oi", "index.dat", ExactTextMatch());
    setOutputFile("-on", "index.ndx", ExactTextMatch());
    excludeDataset("cfrac");
    runTest(CommandLine(cmdline));
}

TEST_F(SelectModuleTest, HandlesPDBOutputWithNonPDBInput)
{
    const char* const cmdline[] = { "select", "-select", "resname RA RD and y < 2.5" };
    setTopology("simple.gro");
    setTrajectory("simple.gro");
    includeDataset("occupancy");
    setOutputFile("-ofpdb", "occupancy.pdb", ExactTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(SelectModuleTest, HandlesPDBOutputWithPDBInput)
{
    const char* const cmdline[] = { "select", "-select", "resname RA RD and y < 2.5" };
    setTopology("simple.pdb");
    setTrajectory("simple.gro");
    includeDataset("occupancy");
    setOutputFile("-ofpdb", "occupancy.pdb", ExactTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(SelectModuleTest, HandlesMaxPDBOutput)
{
    const char* const cmdline[] = { "select",        "-select",   "resname RA RD and y < 2.5",
                                    "resname RA RB", "-pdbatoms", "maxsel" };
    setTopology("simple.pdb");
    setTrajectory("simple.gro");
    includeDataset("occupancy");
    setOutputFile("-ofpdb", "occupancy.pdb", ExactTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(SelectModuleTest, HandlesSelectedPDBOutput)
{
    const char* const cmdline[] = { "select",        "-select",   "resname RA RD and y < 2.5",
                                    "resname RA RB", "-pdbatoms", "selected" };
    setTopology("simple.pdb");
    setTrajectory("simple.gro");
    includeDataset("occupancy");
    setOutputFile("-ofpdb", "occupancy.pdb", ExactTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(SelectModuleTest, NormalizesSizes)
{
    const char* const cmdline[] = { "select",     "-select", "y < 2.5", "resname RA and y < 2.5",
                                    "resname RA", "-norm" };
    setTopology("simple.gro");
    includeDataset("size");
    runTest(CommandLine(cmdline));
}

TEST_F(SelectModuleTest, WritesResidueNumbers)
{
    const char* const cmdline[] = { "select", "-select", "res_com of resname RA RD" };
    setTopology("simple.gro");
    includeDataset("index");
    runTest(CommandLine(cmdline));
}

TEST_F(SelectModuleTest, WritesResidueIndices)
{
    const char* const cmdline[] = {
        "select", "-select", "res_com of resname RA RD", "-resnr", "index"
    };
    setTopology("simple.gro");
    includeDataset("index");
    runTest(CommandLine(cmdline));
}

} // namespace
} // namespace test
} // namespace gmx
