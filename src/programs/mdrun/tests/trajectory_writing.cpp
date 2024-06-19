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
 * Tests for the .mdp nst*out functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testfilemanager.h"

#include "moduletest.h"

namespace
{

// TODO configure these tests so they test all formats for mdrun trajectory writing
#if GMX_USE_TNG

//! Test fixture for mdrun trajectory writing
class TrajectoryWritingTest :
    public gmx::test::MdrunTestFixture,
    public ::testing::WithParamInterface<std::string>
{
public:
    //! The file name of the MDP file
    std::string theMdpFile;

    //! Execute the trajectory writing test
    void runTest(int maxwarn)
    {
        runner_.useStringAsMdpFile(theMdpFile);
        runner_.useTopGroAndNdxFromDatabase("spc-and-methanol");
        runner_.setMaxWarn(maxwarn);
        EXPECT_EQ(0, runner_.callGrompp());

        runner_.fullPrecisionTrajectoryFileName_ =
                fileManager_.getTemporaryFilePath("spc-and-methanol.tng").string();
        runner_.reducedPrecisionTrajectoryFileName_ =
                fileManager_.getTemporaryFilePath("spc-and-methanol-reduced.tng").string();
        ASSERT_EQ(0, runner_.callMdrun());
        // TODO When there is a way to sense something like the
        // output of gmx check, compare the result with that from
        // writing .trr and .xtc and assert the behaviour is
        // correct. Note that TNG will always write the box, even
        // when constant - this will be a source of
        // trajectory-file differences.
    }
};

//! Helper typedef for naming test cases like sentences
typedef TrajectoryWritingTest Trajectories;

/* This test ensures mdrun can write various quantities at various
   frequencies */
TEST_P(Trajectories, ThatDifferInNstxout)
{
    const auto nstxout = GetParam();
    theMdpFile         = gmx::formatString(
            "integrator = md\n"
            "nsteps = 6\n"
            "nstxout = %s\n"
            "nstvout = 2\n"
            "nstfout = 4\n"
            "nstxout-compressed = 5\n"
            "tcoupl = v-rescale\n"
            "tc-grps = System\n"
            "tau-t = 1\n"
            "ref-t = 298\n"
            "compressed-x-grps = Sol\n",
            nstxout.c_str());
    runTest(0);
}

//! Helper typedef for naming test cases like sentences
typedef TrajectoryWritingTest NptTrajectories;

/* This test ensures mdrun can write trajectories in TNG format from NPT ensembles. */
TEST_P(NptTrajectories, WithDifferentPcoupl)
{
    const auto& pcouple = GetParam();
    int         maxwarn = (pcouple == "Berendsen") ? 1 : 0;
    theMdpFile          = gmx::formatString(
            "integrator = md\n"
            "nsteps = 2\n"
            "nstxout = 2\n"
            "nstvout = 1\n"
            "pcoupl = %s\n"
            "tau-p = 1\n"
            "ref-p = 1\n"
            "compressibility = 4.5e-5\n"
            "tcoupl = v-rescale\n"
            "tc-grps = System\n"
            "tau-t = 1\n"
            "ref-t = 298\n",
            pcouple.c_str());
    runTest(maxwarn);
}

// TODO Consider spamming more of the parameter space when we don't
// have to write .mdp and .tpr files to do it.
INSTANTIATE_TEST_SUITE_P(MdrunCanWrite, Trajectories, ::testing::Values("1", "2", "3"));

INSTANTIATE_TEST_SUITE_P(MdrunCanWrite,
                         NptTrajectories,
                         ::testing::Values("no", "Berendsen", "Parrinello-Rahman"));

#endif

} // namespace
