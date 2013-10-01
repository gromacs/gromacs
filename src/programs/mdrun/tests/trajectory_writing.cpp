/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
 * Tests for the .mdp nst*out functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#include <gtest/gtest.h>
#include "gromacs/legacyheaders/string2.h"
#include "moduletest.h"
#include "gromacs/options/filenameoption.h"

namespace
{

//! Max string length
const int g_maxString = 4096;

//! Test fixture for mdrun trajectory writing
class TrajectoryWritingTest :
    public gmx::test::MdrunTestFixture,
    public ::testing::WithParamInterface<const char *>
{
    public:
        //! The file name of the MDP file
        char theMdpFile[g_maxString];

        //! Execute the trajectory writing test
        void runTest()
        {
            useStringAsMdpFile(theMdpFile);
            useTopAndGroFromDatabase("spc-and-methanol");
            EXPECT_EQ(0, callGrompp());

            fullPrecisionTrajectoryFileName    = fileManager_.getTemporaryFilePath("spc-and-methanol.tng");
            reducedPrecisionTrajectoryFileName = fileManager_.getTemporaryFilePath("spc-and-methanol-reduced.tng");
            ASSERT_EQ(0, callMdrun());
            // TODO compare the result with that from writing .trr and
            // .xtc and assert the behaviour is correct. TNG will
            // always write the box, even when constant - this will be
            // a source of trajectory-file differences.
        }
};

//! TODO
typedef TrajectoryWritingTest Trajectories;

/* This test ensures mdrun can write various quantities at various
   frequencies */
TEST_P(Trajectories, ThatDifferInNstxout)
{
    snprintf(theMdpFile, g_maxString,
             "integrator = md\n"
             "nsteps = 6\n"
             "nstxout = %s\n"
             "nstvout = 2\n"
             "nstfout = 4\n"
             "nstxtcout = 5\n"
             "xtc-grps = Sol\n",
             GetParam());
    runTest();
}

// TODO consider spamming more of the parameter space when we don't
// have to write .mdp and .tpr files to do it.
INSTANTIATE_TEST_CASE_P(CanWrite,
                        Trajectories,
                            ::testing::Values("1", "2", "3"));

//! TODO
typedef TrajectoryWritingTest NptTrajectories;

/* This test ensures mdrun can write trajectories in TNG format from NPT ensembles. */
TEST_P(NptTrajectories, WithDifferentPcoupl)
{
    snprintf(theMdpFile, g_maxString,
             "integrator = md\n"
             "nsteps = 2\n"
             "nstxout = 2\n"
             "nstvout = 1\n"
             "pcoupl = %s\n"
             "tau-p = 1\n"
             "ref-p = 1\n"
             "compressibility = 4.5e-5\n",
             GetParam());
    runTest();
}

INSTANTIATE_TEST_CASE_P(CanWrite,
                        NptTrajectories,
                            ::testing::Values("no", "Berendsen", "Parrinello-Rahman"));

/*! \todo
 *
 */

} // namespace
