/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * Tests for functionality of the "rdf" trajectory analysis module.
 *
 * These tests are essentially regression tests for the actual RDF calculation
 * from very small configurations.  Exclusions are not tested (since the input
 * does not contain any), nor is the effect of -cut of -rmax.
 * At the moment, they do not test the final normalization, but only the pair
 * counts calculated for each frame.  Tests for the final normalization should
 * be added once related TODOs in the implementation/framework have been
 * resolved (currently, the test framework does not see the final output data).
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/rdf.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::NoTextMatch;

/********************************************************************
 * Tests for gmx::analysismodules::Rdf.
 */

//! Test fixture for the `rdf` analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::RdfInfo> RdfModuleTest;

TEST_F(RdfModuleTest, BasicTest)
{
    const char* const cmdline[] = { "rdf",     "-bin", "0.05",    "-ref",
                                    "name OW", "-sel", "name OW", "not name OW" };
    setTopology("spc216.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    excludeDataset("pairdist");
    runTest(CommandLine(cmdline));
}

TEST_F(RdfModuleTest, SelectionsSolelyFromIndexFileWork)
{
    const char* const cmdline[] = { "rdf",
                                    "-bin",
                                    "0.05",
                                    // Use selection that names a group in the index file
                                    "-ref",
                                    "name_OW",
                                    // Use selections that name groups in the index file
                                    "-sel",
                                    "name_OW",
                                    "not_name_OW" };
    // Note not supplying a topology file to -s
    setTrajectory("spc216.gro");
    setInputFile("-n", "index.ndx");
    setOutputFile("-o", ".xvg", NoTextMatch());
    excludeDataset("pairdist");
    runTest(CommandLine(cmdline));
}

TEST_F(RdfModuleTest, SelectionsFromBothTopologyFileAndIndexFileWork)
{
    const char* const cmdline[] = { "rdf",
                                    "-bin",
                                    "0.05",
                                    // Use selection whose parsing requires topology file
                                    "-ref",
                                    "name OW",
                                    // Use selections that name groups in the index file
                                    "-sel",
                                    "name_OW",
                                    "not_name_OW" };
    // Note supplying a topology file to -s
    setTopology("spc216.gro");
    setInputFile("-n", "index.ndx");
    setOutputFile("-o", ".xvg", NoTextMatch());
    excludeDataset("pairdist");
    runTest(CommandLine(cmdline));
}

TEST_F(RdfModuleTest, CalculatesSurf)
{
    const char* const cmdline[] = { "rdf",
                                    "-bin",
                                    "0.05",
                                    "-surf",
                                    "res",
                                    "-ref",
                                    "within 0.5 of (resnr 1 and name OW)",
                                    "-sel",
                                    "name OW",
                                    "not name OW" };
    setTopology("spc216.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    excludeDataset("pairdist");
    runTest(CommandLine(cmdline));
}

TEST_F(RdfModuleTest, CalculatesXY)
{
    const char* const cmdline[] = { "rdf",     "-bin", "0.05",    "-xy",        "-ref",
                                    "name OW", "-sel", "name OW", "not name OW" };
    setTopology("spc216.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    excludeDataset("pairdist");
    // TODO: Consider if it is possible to get a more reproducible result
    // and/or a stricter tolerance (e.g., by checking that the sum of
    // neighboring values still stays constant).
    setDatasetTolerance("paircount", gmx::test::absoluteTolerance(2.5));
    runTest(CommandLine(cmdline));
}

} // namespace
} // namespace test
} // namespace gmx
