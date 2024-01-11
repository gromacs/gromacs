/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Tests for functionality of the "msd" trajectory analysis module.
 *
 * \author Kevin Boyd <kevin44boyd@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/msd.h"

#include <gromacs/commandline/cmdlineoptionsmodule.h>
#include <gromacs/trajectoryanalysis/cmdlinerunner.h>
#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

#include "moduletest.h"

namespace gmx::test
{

namespace
{

using gmx::test::CommandLine;

/*! \brief Returns whether or not we care about a header line in an xvg file, for matching purposes.
 *
 * \todo This is mostly taken from xvgtest.cpp. We could probably create some modular checker
 * functionality where each line is compared against a set of subcheckers automatically. Then we
 * could build matchers like this out of the modular components.
 */
bool isRelevantXvgHeader(const std::string& line)
{
    return startsWith(line, "@")
           && (contains(line, " title ") || contains(line, " subtitle ") || contains(line, " label ")
               || contains(line, "@TYPE ") || contains(line, " legend \""));
}

//! Helper function to check a single xvg value in a sequence.
void checkXvgDataPoint(TestReferenceChecker* checker, const std::string& value)
{
    checker->checkRealFromString(value, nullptr);
}

/*! \brief MsdMatcher is effectively an extension of XvgMatcher for gmx msd results.
 *
 * In addition to the usual fields XvgMatcher checks, MsdMatcher checks for properly reported
 * diffusion coefficients.
 */
class MsdMatcher : public ITextBlockMatcher
{
public:
    MsdMatcher() = default;

    void checkStream(TextInputStream* stream, TestReferenceChecker* checker) override
    {
        checker->setDefaultTolerance(
                gmx::test::FloatingPointTolerance(gmx::test::absoluteTolerance(1.0e-5)));
        TestReferenceChecker dCoefficientChecker(
                checker->checkCompound("XvgLegend", "DiffusionCoefficient"));
        TestReferenceChecker legendChecker(checker->checkCompound("XvgLegend", "Legend"));
        TestReferenceChecker dataChecker(checker->checkCompound("XvgData", "Data"));
        std::string          line;
        int                  rowCount = 0;
        while (stream->readLine(&line))
        {
            // Legend and titles.
            if (isRelevantXvgHeader(line))
            {
                legendChecker.checkString(stripString(line.substr(1)), nullptr);
                continue;
            }
            // All other comment lines we don't care about.
            if (startsWith(line, "#") || startsWith(line, "@"))
            {
                continue;
            }
            // Actual data.
            const std::vector<std::string> columns = splitString(line);
            const std::string              id      = formatString("Row%d", rowCount);
            dataChecker.checkSequence(columns.begin(), columns.end(), id.c_str(), &checkXvgDataPoint);
            rowCount++;
        }
    }
};

class MsdMatch : public ITextBlockMatcherSettings
{
public:
    [[nodiscard]] TextBlockMatcherPointer createMatcher() const override
    {
        return std::make_unique<MsdMatcher>();
    }
};

class MsdModuleTest : public gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::MsdInfo>
{
public:
    MsdModuleTest() { setOutputFile("-o", "msd.xvg", MsdMatch()); }
    // Creates a TPR for the given starting structure and topology. Builds an mdp in place prior
    // to calling grompp. sets the -s input to the generated tpr
    void createTpr(const std::string& structure, const std::string& topology, const std::string& index)
    {
        std::string tpr             = fileManager().getTemporaryFilePath(".tpr").u8string();
        std::string mdp             = fileManager().getTemporaryFilePath(".mdp").u8string();
        std::string mdpFileContents = gmx::formatString(
                "cutoff-scheme = verlet\n"
                "rcoulomb      = 0.85\n"
                "rvdw          = 0.85\n"
                "rlist         = 0.85\n");
        gmx::TextWriter::writeFileFromString(mdp, mdpFileContents);

        // Prepare a .tpr file
        CommandLine caller;
        const auto  simDB = gmx::test::TestFileManager::getTestSimulationDatabaseDirectory();
        caller.append("grompp");
        caller.addOption("-maxwarn", 0);
        caller.addOption("-f", mdp.c_str());
        auto gro = std::filesystem::path(simDB).append(structure);
        caller.addOption("-c", gro.u8string().c_str());
        auto top = std::filesystem::path(simDB).append(topology);
        caller.addOption("-p", top.u8string().c_str());
        auto ndx = std::filesystem::path(simDB).append(index);
        caller.addOption("-n", ndx.u8string().c_str());
        caller.addOption("-o", tpr.c_str());
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));

        // setInputFile() doesn't like the temporary tpr path.
        CommandLine& cmdline = commandLine();
        cmdline.addOption("-s", tpr.c_str());
    }

    // Convenience function to set input trajectory, tpr, and index, if all of the input files
    // share a common prefix.
    void setAllInputs(const std::string& prefix)
    {
        setInputFile("-f", prefix + ".xtc");
        setInputFile("-n", prefix + ".ndx");
        createTpr(prefix + ".gro", prefix + ".top", prefix + ".ndx");
    }
};

// ----------------------------------------------
// These tests check the basic MSD and diffusion coefficient reporting capabilities on toy systems.
// Trestart is set to larger than the size of the trajectory so that all frames are only compared
// against the first frame and diffusion coefficients can be predicted exactly.
// ----------------------------------------------

// for 3D, (8 + 4 + 0) / 3 should yield 4 cm^2 / s
TEST_F(MsdModuleTest, threeDimensionalDiffusion)
{
    setInputFile("-f", "msd_traj.xtc");
    setInputFile("-s", "msd_coords.gro");
    setInputFile("-n", "msd.ndx");
    const char* const cmdline[] = {
        "-trestart",
        "200",
        "-sel",
        "0",
    };
    runTest(CommandLine(cmdline));
}

// for lateral z, (8 + 4) / 2 should yield 6 cm^2 /s
TEST_F(MsdModuleTest, twoDimensionalDiffusion)
{
    setInputFile("-f", "msd_traj.xtc");
    setInputFile("-s", "msd_coords.gro");
    setInputFile("-n", "msd.ndx");
    const char* const cmdline[] = { "-trestart", "200", "-lateral", "z", "-sel", "0" };
    runTest(CommandLine(cmdline));
}

// for type x, should yield 8 cm^2 / s
TEST_F(MsdModuleTest, oneDimensionalDiffusion)
{
    setInputFile("-f", "msd_traj.xtc");
    setInputFile("-s", "msd_coords.gro");
    setInputFile("-n", "msd.ndx");
    const char* const cmdline[] = { "-trestart", "200", "-type", "x", "-sel", "0" };
    runTest(CommandLine(cmdline));
}

TEST_F(MsdModuleTest, oneDimensionalDiffusionWithMaxTau)
{
    setInputFile("-f", "msd_traj.xtc");
    setInputFile("-s", "msd_coords.gro");
    setInputFile("-n", "msd.ndx");
    const char* const cmdline[] = { "-trestart", "200", "-type", "x", "-sel", "0", "-maxtau", "5" };
    runTest(CommandLine(cmdline));
}

TEST_F(MsdModuleTest, roundingFail)
{
    // Check that a proper exception is throws when a trajectory is saved too often, #4694
    setInputFile("-f", "msd_traj_rounding_fail.xtc");
    setInputFile("-s", "msd_coords.gro");
    setInputFile("-n", "msd.ndx");
    CommandLine cmdline = commandLine();
    cmdline.addOption("-sel", "0");

    ICommandLineOptionsModulePointer runner(
            TrajectoryAnalysisCommandLineRunner::createModule(createModule()));

    EXPECT_THROW(CommandLineTestHelper::runModuleDirect(std::move(runner), &cmdline), gmx::ToleranceError);
}


// -------------------------------------------------------------------------
// These tests operate on a more realistic trajectory, with a solvated protein,
// with 20 frames at a 2 ps dt. Note that this box is also non-square, so we're validating
// non-trivial pbc removal.

// Group 1 = protein, 2 = all waters, 3 = a subset of 6 water molecules
// ------------------------------------------------------------------------
TEST_F(MsdModuleTest, multipleGroupsWork)
{
    setAllInputs("alanine_vsite_solvated");
    // Restart every frame, select protein and water separately. Note that the reported diffusion
    // coefficient for protein is not well-sampled and doesn't correspond to anything physical.
    const char* const cmdline[] = { "-trestart", "2", "-sel", "1;2" };
    runTest(CommandLine(cmdline));
}

TEST_F(MsdModuleTest, trestartLessThanDt)
{
    setAllInputs("alanine_vsite_solvated");
    const char* const cmdline[] = { "-trestart", "1", "-sel", "2" };
    runTest(CommandLine(cmdline));
}

TEST_F(MsdModuleTest, trestartGreaterThanDt)
{
    setAllInputs("alanine_vsite_solvated");
    const char* const cmdline[] = { "-trestart", "10", "-sel", "2" };
    runTest(CommandLine(cmdline));
}

TEST_F(MsdModuleTest, molTest)
{
    setAllInputs("alanine_vsite_solvated");
    setOutputFile("-mol", "diff_mol.xvg", MsdMatch());
    const char* const cmdline[] = { "-trestart", "10", "-sel", "3" };
    runTest(CommandLine(cmdline));
}

TEST_F(MsdModuleTest, beginFit)
{
    setAllInputs("alanine_vsite_solvated");
    const char* const cmdline[] = { "-trestart", "2", "-sel", "3", "-beginfit", "20" };
    runTest(CommandLine(cmdline));
}

TEST_F(MsdModuleTest, endFit)
{
    setAllInputs("alanine_vsite_solvated");
    const char* const cmdline[] = { "-trestart", "2", "-sel", "3", "-endfit", "25" };
    runTest(CommandLine(cmdline));
}

TEST_F(MsdModuleTest, notEnoughPointsForFitErrorEstimate)
{
    setAllInputs("alanine_vsite_solvated");
    const char* const cmdline[] = { "-trestart", "2",        "-beginfit", "5",    "-endfit",
                                    "9",         "-lateral", "x",         "-sel", "all" };
    runTest(CommandLine(cmdline));
}

} // namespace

} // namespace gmx::test
