/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * \brief Defines functionality used to test mdrun termination
 * functionality under different conditions
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "terminationhelper.h"

#include <filesystem>

#include <gtest/gtest.h>

#include "gromacs/utility/path.h"

#include "testutils/testfilemanager.h"

#include "programs/mdrun/tests/moduletest.h"

namespace gmx
{
namespace test
{

TerminationHelper::TerminationHelper(CommandLine* mdrunCaller, SimulationRunner* runner) :
    mdrunCaller_(mdrunCaller), runner_(runner)
{
    runner_->useTopGroAndNdxFromDatabase("spc2");
}

void TerminationHelper::runFirstMdrun(const std::string& expectedCptFileName)
{
    CommandLine firstPart(*mdrunCaller_);
    // Stop after 0.036 ms, which should be short enough that
    // numSteps isn't reached first.
    firstPart.addOption("-maxh", 1e-7);
    firstPart.addOption("-nstlist", 1);
    ASSERT_EQ(0, runner_->callMdrun(firstPart));
    EXPECT_EQ(true, File::exists(expectedCptFileName, File::returnFalseOnError))
            << expectedCptFileName << " was not found";
}

void TerminationHelper::runSecondMdrun()
{
    CommandLine secondPart(*mdrunCaller_);
    secondPart.addOption("-cpi", runner_->cptOutputFileName_);
    secondPart.addOption("-nsteps", 2);
    ASSERT_EQ(0, runner_->callMdrun(secondPart));
}

void TerminationHelper::runSecondMdrunWithNoAppend()
{
    CommandLine secondPart(*mdrunCaller_);
    secondPart.addOption("-cpi", runner_->cptOutputFileName_);
    secondPart.addOption("-nsteps", 2);
    secondPart.append("-noappend");
    ASSERT_EQ(0, runner_->callMdrun(secondPart));
}

} // namespace test
} // namespace gmx
