/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * \brief Defines functionality used to test mdrun termination
 * functionality under different conditions
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "terminationhelper.h"

#include <gtest/gtest.h>

#include "gromacs/utility/path.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

TerminationHelper::TerminationHelper(TestFileManager  *fileManager,
                                     CommandLine      *mdrunCaller,
                                     SimulationRunner *runner)
    : mdrunCaller_(mdrunCaller), runner_(runner)
{
    runner_->cptFileName_ = fileManager->getTemporaryFilePath(".cpt");
    runner_->useTopGroAndNdxFromDatabase("spc2");
}

void TerminationHelper::runFirstMdrun(const std::string &expectedCptFileName)
{
    CommandLine firstPart(*mdrunCaller_);
    // Stop after 0.036 ms, which should be short enough that
    // numSteps isn't reached first.
    firstPart.addOption("-maxh", 1e-7);
    firstPart.addOption("-nstlist", 1);
    firstPart.addOption("-cpo", runner_->cptFileName_);
    ASSERT_EQ(0, runner_->callMdrun(firstPart));
    EXPECT_EQ(true, File::exists(expectedCptFileName, File::returnFalseOnError)) << expectedCptFileName << " was not found";
}

void TerminationHelper::runSecondMdrun()
{
    CommandLine secondPart(*mdrunCaller_);
    secondPart.addOption("-cpi", runner_->cptFileName_);
    secondPart.addOption("-nsteps", 2);
    ASSERT_EQ(0, runner_->callMdrun(secondPart));
}

} // namespace
} // namespace
