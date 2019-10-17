/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
 * Tests for the mdrun multi-simulation functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "multisimtest.h"

#include <cmath>

#include <algorithm>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"
#include "terminationhelper.h"

namespace gmx
{
namespace test
{

MultiSimTest::MultiSimTest() :
    size_(gmx_node_num()),
    rank_(gmx_node_rank()),
    mdrunCaller_(new CommandLine)

{
    const char* directoryNameFormat = "sim_%d";

    // Modify the file manager to have a temporary directory unique to
    // each simulation. No need to have a mutex on this, nobody else
    // can access the fileManager_ yet because we only just
    // constructed it.
    std::string originalTempDirectory = fileManager_.getOutputTempDirectory();
    std::string newTempDirectory =
            Path::join(originalTempDirectory, formatString(directoryNameFormat, rank_));
    Directory::create(newTempDirectory);
    fileManager_.setOutputTempDirectory(newTempDirectory);

    mdrunCaller_->append("mdrun");
    mdrunCaller_->addOption("-multidir");
    for (int i = 0; i != size_; ++i)
    {
        mdrunCaller_->append(Path::join(originalTempDirectory, formatString(directoryNameFormat, i)));
    }
}

void MultiSimTest::organizeMdpFile(SimulationRunner* runner, const char* controlVariable, int numSteps)
{
    const real  baseTemperature = 298;
    const real  basePressure    = 1;
    std::string mdpFileContents = formatString(
            "nsteps = %d\n"
            "nstlog = 1\n"
            "nstcalcenergy = 1\n"
            "tcoupl = v-rescale\n"
            "tc-grps = System\n"
            "tau-t = 1\n"
            "ref-t = %f\n"
            // pressure coupling (if active)
            "tau-p = 1\n"
            "ref-p = %f\n"
            "compressibility = 4.5e-5\n"
            // velocity generation
            "gen-vel = yes\n"
            "gen-temp = %f\n"
            // control variable specification
            "%s\n",
            numSteps, baseTemperature + 0.0001 * rank_, basePressure * std::pow(1.01, rank_),
            /* Set things up so that the initial KE decreases with
               increasing replica number, so that the (identical)
               starting PE decreases on the first step more for the
               replicas with higher number, which will tend to force
               replica exchange to occur. */
            std::max(baseTemperature - 10 * rank_, real(0)), controlVariable);
    runner->useStringAsMdpFile(mdpFileContents);
}

void MultiSimTest::runExitsNormallyTest()
{
    if (size_ <= 1)
    {
        /* Can't test multi-sim without multiple ranks. */
        return;
    }

    SimulationRunner runner(&fileManager_);
    runner.useTopGroAndNdxFromDatabase("spc2");

    const char* pcoupl = GetParam();
    organizeMdpFile(&runner, pcoupl);
    /* Call grompp on every rank - the standard callGrompp() only runs
       grompp on rank 0. */
    EXPECT_EQ(0, runner.callGromppOnThisRank());

    ASSERT_EQ(0, runner.callMdrun(*mdrunCaller_));
}

void MultiSimTest::runMaxhTest()
{
    if (size_ <= 1)
    {
        /* Can't test replica exchange without multiple ranks. */
        return;
    }

    SimulationRunner runner(&fileManager_);
    runner.useTopGroAndNdxFromDatabase("spc2");

    TerminationHelper helper(&fileManager_, mdrunCaller_.get(), &runner);
    // Make sure -maxh has a chance to propagate
    int numSteps = 100;
    organizeMdpFile(&runner, "pcoupl = no", numSteps);
    /* Call grompp on every rank - the standard callGrompp() only runs
       grompp on rank 0. */
    EXPECT_EQ(0, runner.callGromppOnThisRank());

    helper.runFirstMdrun(runner.cptFileName_);
    helper.runSecondMdrun();
}

} // namespace test
} // namespace gmx
