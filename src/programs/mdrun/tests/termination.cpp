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
 * \brief
 * Tests for the mdrun termination functionality
 *
 * \todo This approach is not very elegant, but "stuff doesn't
 * segfault or give a fatal error" is a useful result. We can improve
 * it when we can mock out more do_md() functionality. Before that,
 * we'd probably prefer not to run this test case in per-patchset
 * verification, but this is the best we can do for now.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <cstdio>

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "moduletest.h"
#include "terminationhelper.h"

namespace gmx
{
namespace test
{

//! Build a simple .mdp file
static void organizeMdpFile(SimulationRunner* runner, int nsteps = 2)
{
    // Make sure -maxh has a chance to propagate
    runner->useStringAsMdpFile(
            formatString("nsteps = %d\n"
                         "tcoupl = v-rescale\n"
                         "tc-grps = System\n"
                         "tau-t = 1\n"
                         "ref-t = 298\n",
                         nsteps));
}

//! Convenience typedef
typedef MdrunTestFixture MdrunTerminationTest;

TEST_F(MdrunTerminationTest, CheckpointRestartAppendsByDefault)
{
    runner_.cptOutputFileName_ = fileManager_.getTemporaryFilePath(".cpt").string();

    runner_.useTopGroAndNdxFromDatabase("spc2");
    organizeMdpFile(&runner_);
    EXPECT_EQ(0, runner_.callGrompp());

    SCOPED_TRACE("Running the first simulation part");
    {
        CommandLine firstPart;
        firstPart.append("mdrun");
        ASSERT_EQ(0, runner_.callMdrun(firstPart));
        ASSERT_TRUE(File::exists(runner_.cptOutputFileName_, File::returnFalseOnError))
                << runner_.cptOutputFileName_ << " was not found and should be";
    }
    SCOPED_TRACE("Running the second simulation part with default appending behavior");
    {
        runner_.changeTprNsteps(4);

        CommandLine secondPart;
        secondPart.append("mdrun");
        secondPart.addOption("-cpi", runner_.cptOutputFileName_);
        ASSERT_EQ(0, runner_.callMdrun(secondPart));

        auto logFileContents = TextReader::readFileToString(runner_.logFileName_);
        EXPECT_NE(std::string::npos, logFileContents.find("Restarting from checkpoint, appending to previous log file"))
                << "appending was not detected";
    }
}

TEST_F(MdrunTerminationTest, WritesCheckpointAfterMaxhTerminationAndThenRestarts)
{
    runner_.useTopGroAndNdxFromDatabase("spc2");
    organizeMdpFile(&runner_, 100);
    EXPECT_EQ(0, runner_.callGrompp());

    SCOPED_TRACE("Running the first simulation part with -maxh");
    {
        CommandLine firstPart;
        firstPart.append("mdrun");
        // Ensure maxh will trigger the halt, and that the signal will
        // have time to be propagated.
        //
        // TODO It would be nicer to set nstlist in the .mdp file, but
        // then it is not a command.
        firstPart.addOption("-maxh", 1e-7);
        firstPart.addOption("-nstlist", 1);
        ASSERT_EQ(0, runner_.callMdrun(firstPart));
        EXPECT_EQ(true, File::exists(runner_.cptOutputFileName_, File::returnFalseOnError))
                << runner_.cptOutputFileName_ << " was not found";
    }

    SCOPED_TRACE("Running the second simulation part");
    {
        runner_.changeTprNsteps(102);

        CommandLine secondPart;
        secondPart.append("mdrun");
        secondPart.addOption("-cpi", runner_.cptOutputFileName_);
        ASSERT_EQ(0, runner_.callMdrun(secondPart));

        auto logFileContents = TextReader::readFileToString(runner_.logFileName_);
        EXPECT_NE(std::string::npos, logFileContents.find("Writing checkpoint, step 102"))
                << "completion of restarted simulation was not detected";
    }
}

TEST_F(MdrunTerminationTest, CheckpointRestartWithNoAppendWorksAndCannotLaterAppend)
{
    runner_.useTopGroAndNdxFromDatabase("spc2");
    organizeMdpFile(&runner_);
    EXPECT_EQ(0, runner_.callGrompp());

    SCOPED_TRACE("Running the first simulation part");
    {
        CommandLine firstPart;
        firstPart.append("mdrun");
        ASSERT_EQ(0, runner_.callMdrun(firstPart));
        EXPECT_EQ(true, File::exists(runner_.cptOutputFileName_, File::returnFalseOnError))
                << runner_.cptOutputFileName_ << " was not found";
    }

    SCOPED_TRACE("Running the second simulation part with -noappend");
    {
        runner_.changeTprNsteps(4);

        CommandLine secondPart;
        secondPart.append("mdrun");
        secondPart.addOption("-cpi", runner_.cptOutputFileName_);
        secondPart.append("-noappend");
        ASSERT_EQ(0, runner_.callMdrun(secondPart));

        auto expectedLogFileName = fileManager_.getTemporaryFilePath(".part0002.log");
        ASSERT_EQ(true, File::exists(expectedLogFileName, File::returnFalseOnError))
                << expectedLogFileName << " was not found";
        auto expectedEdrFileName = fileManager_.getTemporaryFilePath(".part0002.edr");
        ASSERT_EQ(true, File::exists(expectedEdrFileName, File::returnFalseOnError))
                << expectedEdrFileName << " was not found";
    }

    SCOPED_TRACE("Running the third simulation part with -append, which will fail");
    runner_.logFileName_ = fileManager_.getTemporaryFilePath(".part0002.log").string();
    runner_.changeTprNsteps(6);

    {
        CommandLine thirdPart;
        thirdPart.append("mdrun");
        thirdPart.addOption("-cpi", runner_.cptOutputFileName_);
        thirdPart.append("-append");
        EXPECT_THROW_GMX(runner_.callMdrun(thirdPart), InconsistentInputError);
    }
    SCOPED_TRACE("Running the third simulation part with -noappend");
    {
        CommandLine thirdPart;
        thirdPart.append("mdrun");
        thirdPart.addOption("-cpi", runner_.cptOutputFileName_);
        thirdPart.append("-noappend");
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(".part0003.edr").string();
        ASSERT_EQ(0, runner_.callMdrun(thirdPart));

        auto expectedLogFileName = fileManager_.getTemporaryFilePath(".part0003.log");
        EXPECT_EQ(true, File::exists(expectedLogFileName, File::returnFalseOnError))
                << expectedLogFileName << " was not found";
        auto expectedEdrFileName = fileManager_.getTemporaryFilePath(".part0003.edr");
        ASSERT_EQ(true, File::exists(expectedEdrFileName, File::returnFalseOnError))
                << expectedEdrFileName << " was not found";
    }
    SCOPED_TRACE("Running the fourth simulation part with default appending");
    runner_.changeTprNsteps(8);
    {
        CommandLine fourthPart;
        fourthPart.append("mdrun");
        fourthPart.addOption("-cpi", runner_.cptOutputFileName_);
        // TODO this is necessary, but ought not be. Is this the issue in Issue #2804?
        fourthPart.append("-noappend");
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(".part0004.edr").string();
        runner_.logFileName_ = fileManager_.getTemporaryFilePath(".part0004.log").string();
        ASSERT_EQ(0, runner_.callMdrun(fourthPart));

        auto expectedLogFileName = fileManager_.getTemporaryFilePath(".part0004.log");
        ASSERT_EQ(true, File::exists(expectedLogFileName, File::returnFalseOnError))
                << expectedLogFileName << " was not found";
        auto expectedEdrFileName = fileManager_.getTemporaryFilePath(".part0004.edr");
        ASSERT_EQ(true, File::exists(expectedEdrFileName, File::returnFalseOnError))
                << expectedEdrFileName << " was not found";
    }
    SCOPED_TRACE("Running the fifth simulation part with no extra steps");
    {
        CommandLine fifthPart;
        fifthPart.append("mdrun");
        fifthPart.addOption("-cpi", runner_.cptOutputFileName_);
        // TODO this is necessary, but ought not be. Is this the issue in Issue #2804?
        fifthPart.append("-noappend");
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(".part0005.edr").string();
        runner_.logFileName_ = fileManager_.getTemporaryFilePath(".part0005.log").string();
        ASSERT_EQ(0, runner_.callMdrun(fifthPart));

        auto expectedLogFileName = fileManager_.getTemporaryFilePath(".part0005.log");
        ASSERT_EQ(true, File::exists(expectedLogFileName, File::returnFalseOnError))
                << expectedLogFileName << " was not found";
        auto expectedEdrFileName = fileManager_.getTemporaryFilePath(".part0005.edr");
        ASSERT_EQ(true, File::exists(expectedEdrFileName, File::returnFalseOnError))
                << expectedEdrFileName << " was not found";
    }
}

TEST_F(MdrunTerminationTest, CheckpointRestartWorksEvenWithMissingCheckpointFile)
{
    runner_.useTopGroAndNdxFromDatabase("spc2");
    organizeMdpFile(&runner_);
    EXPECT_EQ(0, runner_.callGrompp());

    SCOPED_TRACE("Running the first simulation part");
    {
        CommandLine firstPart;
        firstPart.append("mdrun");
        ASSERT_EQ(0, runner_.callMdrun(firstPart));
        EXPECT_EQ(true, File::exists(runner_.cptOutputFileName_, File::returnFalseOnError))
                << runner_.cptOutputFileName_ << " was not found";
    }

    SCOPED_TRACE("Running the second simulation part after deleting the checkpoint file");
    {
        runner_.changeTprNsteps(4);

        CommandLine secondPart;
        secondPart.append("mdrun");
        secondPart.addOption("-cpi", runner_.cptOutputFileName_);

        // Remove the checkpoint, so technically this can no longer be
        // a restart. But it starts again from the beginning anyway.
        //
        // TODO what do we want the behaviour to be?
        std::remove(runner_.cptOutputFileName_.c_str());

        ASSERT_EQ(0, runner_.callMdrun(secondPart));
        auto logFileContents = TextReader::readFileToString(runner_.logFileName_);
        EXPECT_EQ(std::string::npos, logFileContents.find("Restarting from checkpoint, appending to previous log file"))
                << "appending was not detected";
    }
}

TEST_F(MdrunTerminationTest, CheckpointRestartWorksEvenWithAppendAndMissingCheckpointFile)
{
    runner_.useTopGroAndNdxFromDatabase("spc2");
    organizeMdpFile(&runner_);
    EXPECT_EQ(0, runner_.callGrompp());

    SCOPED_TRACE("Running the first simulation part");
    {
        CommandLine firstPart;
        firstPart.append("mdrun");
        ASSERT_EQ(0, runner_.callMdrun(firstPart));
        EXPECT_EQ(true, File::exists(runner_.cptOutputFileName_, File::returnFalseOnError))
                << runner_.cptOutputFileName_ << " was not found";
    }

    SCOPED_TRACE(
            "Running the second simulation part with -append after deleting the checkpoint file");
    {
        runner_.changeTprNsteps(4);

        CommandLine secondPart;
        secondPart.append("mdrun");
        secondPart.addOption("-cpi", runner_.cptOutputFileName_);
        secondPart.append("-append");

        // Remove the checkpoint, so this can no longer be a
        // restart.
        std::remove(runner_.cptOutputFileName_.c_str());

        EXPECT_THROW_GMX(runner_.callMdrun(secondPart), InconsistentInputError);
    }
}

TEST_F(MdrunTerminationTest, RunWithNoAppendCreatesPartFiles)
{
    runner_.useTopGroAndNdxFromDatabase("spc2");
    organizeMdpFile(&runner_);
    EXPECT_EQ(0, runner_.callGrompp());

    SCOPED_TRACE("Running the first simulation part with -noappend");
    {
        CommandLine firstPart;
        firstPart.append("mdrun");
        firstPart.append("-noappend");
        ASSERT_EQ(0, runner_.callMdrun(firstPart));
        auto expectedLogFileName = fileManager_.getTemporaryFilePath(".part0001.log");
        ASSERT_EQ(true, File::exists(expectedLogFileName, File::returnFalseOnError))
                << expectedLogFileName << " was not found";
        auto expectedEdrFileName = fileManager_.getTemporaryFilePath(".part0001.edr");
        ASSERT_EQ(true, File::exists(expectedEdrFileName, File::returnFalseOnError))
                << expectedEdrFileName << " was not found";
        EXPECT_EQ(true, File::exists(runner_.cptOutputFileName_, File::returnFalseOnError))
                << runner_.cptOutputFileName_ << " was not found";
    }

    SCOPED_TRACE("Running the second simulation part with -noappend");
    {
        runner_.changeTprNsteps(4);

        CommandLine secondPart;
        secondPart.append("mdrun");
        secondPart.addOption("-cpi", runner_.cptOutputFileName_);
        secondPart.append("-noappend");
        ASSERT_EQ(0, runner_.callMdrun(secondPart));

        auto expectedLogFileName = fileManager_.getTemporaryFilePath(".part0002.log");
        ASSERT_EQ(true, File::exists(expectedLogFileName, File::returnFalseOnError))
                << expectedLogFileName << " was not found";
        auto expectedEdrFileName = fileManager_.getTemporaryFilePath(".part0002.edr");
        ASSERT_EQ(true, File::exists(expectedEdrFileName, File::returnFalseOnError))
                << expectedEdrFileName << " was not found";
    }
}

} // namespace test
} // namespace gmx
