/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Tests utilities for "Computational Electrophysiology" setups.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/xvgtest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

class CompelTest : public MdrunTestFixture
{
};

/* This test ensures that the compel protocol can be run, that all of
 * the swapcoords parameters from the .mdp file are understood, that
 * the output .xvg file is correct, and that the swap state variables
 * can be written to and read from checkpoint.
 *
 * The setup checks whether the swap-coordinates protocol performs
 * ion/water exchanges in the z direction.  The system is a typical
 * computational electrophysiology setup with two octanol membranes in
 * the x-y plane with water and ions inbetween. At the start, there
 * are 9 anions and 9 cations in one compartment (A) and 10 anions and
 * 10 cations in the other (B). The swap parameters are set such that
 * all ions should end up in compartment A in the swap attempted
 * during the first step, and remain there.
 *
 * The former swap_z regressiontest did the first 2 steps of this test,
 * ie. without the checkpoint restart. */
TEST_F(CompelTest, SwapCanRun)
{
    runner_.useTopGroAndNdxFromDatabase("OctaneSandwich-z");
    const std::string mdpContents = R"(
        dt                       = 0.004
        nsteps                   = 2
        tcoupl                   = V-rescale
        tc-grps                  = System
        tau-t                    = 0.5
        ref-t                    = 300
        constraints              = all-bonds
        cutoff-scheme            = Verlet
        swapcoords               = Z
        swap_frequency           = 1
        split_group0             = Ch0
        split_group1             = Ch1
        massw_split0             = yes
        massw_split1             = no
        solvent_group            = SOL
        cyl0_r                   = 0.9
        cyl0_up                  = 0.75
        cyl0_down                = 0.75
        cyl1_r                   = 0.9
        cyl1_up                  = 0.75
        cyl1_down                = 0.75
        coupl_steps              = 1
        iontypes                 = 2
        iontype0-name            = NA+
        iontype0-in-A            = 19
        iontype0-in-B            = 0
        iontype1-name            = CL-
        iontype1-in-A            = 19
        iontype1-in-B            = 0
        threshold                = 1
     )";

    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    runner_.swapFileName_ = fileManager_.getTemporaryFilePath("swap.xvg").string();

    ::gmx::test::CommandLine swapCaller;
    swapCaller.addOption("-swap", runner_.swapFileName_);

    // Do an initial mdrun that writes a checkpoint file
    ::gmx::test::CommandLine firstCaller(swapCaller);
    ASSERT_EQ(0, runner_.callMdrun(firstCaller));

    // Test the swap output xvg file
    const auto           xvgTolerance = relativeToleranceAsFloatingPoint(5.0, 1e-4);
    TestReferenceData    data;
    TestReferenceChecker checker = data.rootChecker();
    if (File::exists(runner_.swapFileName_, File::returnFalseOnError))
    {
        TestReferenceChecker fileChecker(checker.checkCompound("File", "swap output after 2 steps"));
        TextInputFile    swapXvgFileStream(runner_.swapFileName_);
        XvgMatchSettings matchSettings;
        matchSettings.tolerance = xvgTolerance;
        checkXvgFile(&swapXvgFileStream, &fileChecker, matchSettings);
    }

    // Continue mdrun from that checkpoint file
    ::gmx::test::CommandLine secondCaller(swapCaller);
    secondCaller.addOption("-cpi", runner_.cptOutputFileName_);
    runner_.nsteps_ = 2;
    ASSERT_EQ(0, runner_.callMdrun(secondCaller));

    // Test the updated swap output file
    if (File::exists(runner_.swapFileName_, File::returnFalseOnError))
    {
        TestReferenceChecker fileChecker(checker.checkCompound("File", "swap output after 4 steps"));
        TextInputFile    swapXvgFileStream(runner_.swapFileName_);
        XvgMatchSettings matchSettings;
        matchSettings.tolerance = xvgTolerance;
        checkXvgFile(&swapXvgFileStream, &fileChecker, matchSettings);
    }
}


/*! \todo Add other tests for the compel module, e.g.
 *
 *  - a test that checks that actually ion/water swaps have been done, by
 *    calling gmxcheck on the swap output file and a reference file
 */


} // namespace test
} // namespace gmx
