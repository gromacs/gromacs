/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Tests for the mdrun -rerun functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#include <gtest/gtest.h>
#include "moduletest.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/utility/file.h"

int gmx_cmain(int argc, char *argv[]);

namespace
{

//! Test fixture for mdrun -rerun
typedef gmx::test::MdrunTestFixture RerunTest;

/* Among other things, this test ensures mdrun can read a trajectory. */
TEST_F(RerunTest, RerunExitsNormally)
{
    gmx::File::writeFileFromString(mdpFileName,"");
    topFileName = fileManager.getInputFilePath("spc2.top");
    groFileName = fileManager.getInputFilePath("spc2.gro");
    ASSERT_EQ(numThreads, 1) << "The toy system is so small only one thread can be used";
    ASSERT_EQ(numOpenMPThreads, 1) << "The toy system is so small only one OpenMP thread can be used";
    EXPECT_EQ(0, callGrompp());

    std::string rerunTrrFileName = fileManager.getInputFilePath("spc2.trr");
    setProgramOption("mdrun", "rerun", rerunTrrFileName);
    ASSERT_EQ(0, callMdrun());

    // TODO find a way to fake stdin
    //    createProgramCaller("gmx traj");
    /*
    createProgramCaller("gmx traj", gmx_cmain);
    setProgramOption<gmx::FileNameOption,std::string>("gmx traj", "f", rerunTrrFileName);
    setProgramOption<gmx::FileNameOption,std::string>("gmx traj", "s", tprFileName);
    runProgram("gmx traj");
    */
}

/* TODO add other tests for mdrun -rerun, e.g.
 *
 * - RerunReproducesRunWhenRunOnlyWroteEnergiesOnNeighborSearchSteps
 *   (e.g. do such a run, do a rerun, call gmxcheck)
 */

//! Test fixture for tpi integrator (which uses the re-run machinery)
typedef gmx::test::MdrunTestFixture TpiTest;

TEST_F(TpiTest, TpiExitsNormally)
{
    gmx::File::writeFileFromString(mdpFileName,
                                   "integrator = tpi\n"
                                   "tcoupl = vrescale\n"
                                   "tc-grps = System\n"
                                   "tau-t = 1\n"
                                   "rlist = 0.9\n"
                                   "rcoulomb = 0.9\n"
                                   "rvdw = 0.9\n"
                                   "ref-t = 298\n");

    //    topFileName = fileManager.getInputFilePath("spc2tpi.top");
    //    groFileName = fileManager.getInputFilePath("spc2tpi.gro");
    topFileName = fileManager.getInputFilePath("spc216tpi.top");
    groFileName = fileManager.getInputFilePath("spc216tpi.gro");
    ASSERT_EQ(numThreads, 1) << "The toy system is so small only one thread can be used";
    ASSERT_EQ(numOpenMPThreads, 1) << "The toy system is so small only one OpenMP thread can be used";
    EXPECT_EQ(0, callGrompp());

    std::string rerunTrrFileName = fileManager.getInputFilePath("spc216.gro");
    setProgramOption("mdrun", "rerun", rerunTrrFileName);
    std::string tpiFileName = fileManager.getTemporaryFilePath("tpi.xvg");
    setProgramOption("mdrun", "tpi", tpiFileName);
    setProgramOption("mdrun", "tpid", fileManager.getTemporaryFilePath("tpidist.xvg"));
    ASSERT_EQ(0, callMdrun());
}

} // namespace
