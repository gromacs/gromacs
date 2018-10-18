/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Tests for functionality of the "scattering" trajectory analysis module.
 *
 * \author Joe Jordan <e.jjordan12@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/scattering.h"

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;
using gmx::test::XvgMatch;

/********************************************************************
 * Tests for gmx::analysismodules::Scattering.
 */

//! Test fixture for the scattering module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::ScatteringInfo>
    ScatteringModuleTest;

TEST_F(ScatteringModuleTest, ComputesNeutronScattering)
{
    const char *const cmdline[] = {
        "scattering",
        "-sel", "name CA",
        "-norm", "yes",
        "-scatter-type", "neutron",
        "-numq", "11",
    };
    setTopology("lysozyme-water.tpr");
    setTrajectory("lysozyme-water.xtc");
    setOutputFile("-o", ".xvg", XvgMatch());
    setDatasetTolerance("scattering", gmx::test::absoluteTolerance(0.000001));
    runTest(CommandLine(cmdline));
}

TEST_F(ScatteringModuleTest, ComputesNeutronScatteringMC)
{
    const char *const cmdline[] = {
        "scattering",
        "-sel", "name CA",
        "-norm", "yes",
        "-scatter-type", "neutron",
        "-numq", "11",
        "-mc",
    };
    setTopology("lysozyme-water.tpr");
    setTrajectory("lysozyme-water.xtc");
    setOutputFile("-o", ".xvg", XvgMatch());
    setDatasetTolerance("scattering", gmx::test::absoluteTolerance(0.000001));
    runTest(CommandLine(cmdline));
}

TEST_F(ScatteringModuleTest, ComputesXrayScattering)
{
    const char *const cmdline[] = {
        "scattering",
        "-sel", "name CA",
        "-norm", "yes",
        "-scatter-type", "xray",
        "-numq", "11",
    };
    setTopology("lysozyme-water.tpr");
    setTrajectory("lysozyme-water.xtc");
    setOutputFile("-o", ".xvg", XvgMatch());
    setDatasetTolerance("scattering", gmx::test::absoluteTolerance(0.000001));
    runTest(CommandLine(cmdline));
}

TEST_F(ScatteringModuleTest, ComputesXrayScatteringMC)
{
    const char *const cmdline[] = {
        "scattering",
        "-sel", "name CA",
        "-norm", "yes",
        "-scatter-type", "xray",
        "-numq", "11",
        "-mc",
    };
    setTopology("lysozyme-water.tpr");
    setTrajectory("lysozyme-water.xtc");
    setOutputFile("-o", ".xvg", XvgMatch());
    setDatasetTolerance("scattering", gmx::test::absoluteTolerance(0.000001));
    runTest(CommandLine(cmdline));
}

} // namespace
