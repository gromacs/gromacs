/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * Test routines that parse mdp fields from grompp input and writes
 * mdp back out.
 *
 * In particular these will provide test coverage as we refactor to
 * use a new Options-based key-value-style mdp implementation to
 * support a more modular mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/readir.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/warninp.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
namespace gmx
{
namespace test
{

class GetIrTest : public ::testing::Test
{
public:
    GetIrTest()
    {
        snew(opts_.include, STRLEN);
        snew(opts_.define, STRLEN);
    }
    ~GetIrTest() override
    {
        done_inputrec_strings();
        sfree(opts_.include);
        sfree(opts_.define);
    }

    //! Tells whether warnings and/or errors are expected from inputrec parsing and checking, and whether we should compare the output
    enum class TestBehavior
    {
        NoErrorAndCompareOutput,      //!< Expect no warnings/error and compare output
        NoErrorAndDoNotCompareOutput, //!< Expect no warnings/error and do not compare output
        ErrorAndCompareOutput,        //!< Expect at least one warning/error and compare output
        ErrorAndDoNotCompareOutput //!< Expect at least one warning/error and do not compare output
    };

    /*! \brief Test mdp reading and writing
     *
     * \todo Modernize read_inp and write_inp to use streams,
     * which will make these tests run faster, because they don't
     * use disk files. */
    void runTest(const std::string& inputMdpFileContents,
                 const TestBehavior testBehavior = TestBehavior::NoErrorAndCompareOutput,
                 const bool         setGenVelSeedToKnownValue = true)
    {
        const bool expectError = testBehavior == TestBehavior::ErrorAndCompareOutput
                                 || testBehavior == TestBehavior::ErrorAndDoNotCompareOutput;
        const bool compareOutput = testBehavior == TestBehavior::ErrorAndCompareOutput
                                   || testBehavior == TestBehavior::NoErrorAndCompareOutput;

        WarningHandler wi{ false, 0 };
        std::string    inputMdpFilename = fileManager_.getTemporaryFilePath("input.mdp").u8string();
        std::string    outputMdpFilename;
        if (compareOutput)
        {
            outputMdpFilename = fileManager_.getTemporaryFilePath("output.mdp").u8string();
        }
        if (setGenVelSeedToKnownValue)
        {
            TextWriter::writeFileFromString(inputMdpFilename,
                                            inputMdpFileContents + "\n gen-seed = 256\n");
        }
        else
        {
            TextWriter::writeFileFromString(inputMdpFilename, inputMdpFileContents);
        }

        get_ir(inputMdpFilename.c_str(),
               outputMdpFilename.empty() ? nullptr : outputMdpFilename.c_str(),
               &mdModules_,
               &ir_,
               &opts_,
               WriteMdpHeader::no,
               &wi);

        check_ir(inputMdpFilename.c_str(), mdModules_.notifiers(), &ir_, &opts_, &wi);
        // Now check
        bool failure = warning_errors_exist(wi);
        EXPECT_EQ(failure, expectError);

        if (compareOutput)
        {
            TestReferenceData    data;
            TestReferenceChecker checker(data.rootChecker());
            checker.checkBoolean(failure, "Error parsing mdp file");

            auto outputMdpContents = TextReader::readFileToString(outputMdpFilename);
            checker.checkString(outputMdpContents, "OutputMdpFile");
        }
    }

    TestFileManager fileManager_;
    t_inputrec      ir_;
    MDModules       mdModules_;
    t_gromppopts    opts_;
};

TEST_F(GetIrTest, HandlesDifferentKindsOfMdpLines)
{
    const char*    inputMdpFile[] = { "; File to run my simulation",
                                   "title = simulation",
                                   "define = -DBOOLVAR -DVAR=VALUE",
                                   ";",
                                   "xtc_grps = System ; was Protein",
                                   "include = -I/home/me/stuff",
                                   "",
                                   "tau-t = 0.1 0.3",
                                   "ref-t = ;290 290",
                                   "tinit = 0.3",
                                   "init_step = 0",
                                   "nstcomm = 100",
                                   "integrator = steep" };
    WarningHandler wi{ false, 0 };
    runTest(joinStrings(inputMdpFile, "\n"));
}

TEST_F(GetIrTest, RejectsNonCommentLineWithNoEquals)
{
    const char* inputMdpFile = "title simulation";

    GMX_EXPECT_DEATH_IF_SUPPORTED(runTest(inputMdpFile), "No '=' to separate");
}

TEST_F(GetIrTest, AcceptsKeyWithoutValue)
{
    // Users are probably using lines like this
    const char* inputMdpFile = "xtc_grps = ";
    runTest(inputMdpFile);
}

TEST_F(GetIrTest, RejectsValueWithoutKey)
{
    const char* inputMdpFile = "= -I/home/me/stuff";
    GMX_EXPECT_DEATH_IF_SUPPORTED(runTest(inputMdpFile), "No .mdp parameter name was found");
}

TEST_F(GetIrTest, RejectsEmptyKeyAndEmptyValue)
{
    const char* inputMdpFile = " = ";
    GMX_EXPECT_DEATH_IF_SUPPORTED(runTest(inputMdpFile),
                                  "No .mdp parameter name or value was found");
}

TEST_F(GetIrTest, AcceptsDefineParametersWithValuesIncludingAssignment)
{
    const char* inputMdpFile[] = {
        "define = -DBOOL -DVAR=VALUE",
    };
    runTest(joinStrings(inputMdpFile, "\n"));
}

TEST_F(GetIrTest, AcceptsEmptyLines)
{
    const char* inputMdpFile = "";
    runTest(inputMdpFile);
}

TEST_F(GetIrTest, MtsCheckNstcalcenergy)
{
    const char* inputMdpFile[] = {
        "mts = yes", "mts-levels = 2", "mts-level2-factor = 2", "nstcalcenergy = 5"
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, MtsCheckNstenergy)
{
    const char* inputMdpFile[] = {
        "mts = yes", "mts-levels = 2", "mts-level2-factor = 2", "nstenergy = 5"
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, MtsCheckNstpcouple)
{
    const char* inputMdpFile[] = { "mts = yes",
                                   "mts-levels = 2",
                                   "mts-level2-factor = 2",
                                   "pcoupl = Berendsen",
                                   "nstpcouple = 5" };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, MtsCheckNstdhdl)
{
    const char* inputMdpFile[] = {
        "mts = yes", "mts-level2-factor = 2", "free-energy = yes", "nstdhdl = 5"
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, MtsCheckSDNotSupported)
{
    const char* inputMdpFile[] = { "mts = yes", "integrator = sd" };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

// These tests observe how the electric-field keys behave, since they
// are currently the only ones using the new Options-style handling.
TEST_F(GetIrTest, AcceptsElectricField)
{
    const char* inputMdpFile = "electric-field-x = 1.2 0 0 0";
    runTest(inputMdpFile);
}

TEST_F(GetIrTest, AcceptsElectricFieldPulsed)
{
    const char* inputMdpFile = "electric-field-y = 3.7 2.0 6.5 1.0";
    runTest(inputMdpFile);
}

TEST_F(GetIrTest, AcceptsElectricFieldOscillating)
{
    const char* inputMdpFile = "electric-field-z = 3.7 7.5 0 0";
    runTest(inputMdpFile);
}

TEST_F(GetIrTest, RejectsDuplicateOldAndNewKeys)
{
    const char* inputMdpFile[] = { "verlet-buffer-drift = 1.3", "verlet-buffer-tolerance = 2.7" };
    GMX_EXPECT_DEATH_IF_SUPPORTED(runTest(joinStrings(inputMdpFile, "\n")),
                                  "A parameter is present with both");
}

TEST_F(GetIrTest, AcceptsImplicitSolventNo)
{
    const char* inputMdpFile = "implicit-solvent = no";
    runTest(inputMdpFile);
}

TEST_F(GetIrTest, RejectsImplicitSolventYes)
{
    const char* inputMdpFile = "implicit-solvent = yes";
    GMX_EXPECT_DEATH_IF_SUPPORTED(runTest(inputMdpFile), "Invalid enum");
}

TEST_F(GetIrTest, AcceptsMimic)
{
    const char* inputMdpFile[] = { "integrator = mimic", "QMMM-grps = QMatoms" };
    runTest(joinStrings(inputMdpFile, "\n"));
}

#if HAVE_MUPARSER

TEST_F(GetIrTest, AcceptsTransformationCoord)
{
    const char* inputMdpFile[] = {
        "pull = yes",
        "pull-ngroups = 2",
        "pull-ncoords = 2",
        "pull-coord1-geometry = distance",
        "pull-coord1-groups = 1 2",
        "pull-coord1-k = 1",
        "pull-coord2-geometry = transformation",
        "pull-coord2-expression = 1/(1+x1)",
        "pull-coord2-k = 10",
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::NoErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, InvalidTransformationCoordWithConstraint)
{
    const char* inputMdpFile[] = {
        "pull = yes",
        "pull-ncoords = 1",
        "pull-coord1-geometry = transformation",
        "pull-coord1-type = constraint", // INVALID
        "pull-coord1-expression = 10",
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, InvalidPullCoordWithConstraintInTransformationExpression)
{
    const char* inputMdpFile[] = {
        "pull = yes",
        "pull-ngroups = 2",
        "pull-ncoords = 2",
        "pull-coord1-geometry = distance",
        "pull-coord1-type = constraint", // INVALID
        "pull-coord1-groups = 1 2",
        "pull-coord2-geometry = transformation",
        "pull-coord2-expression = x1",
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, InvalidTransformationCoordDxValue)
{
    const char* inputMdpFile[] = {
        "pull = yes",
        "pull-ncoords = 1",
        "pull-coord1-geometry = transformation",
        "pull-coord1-expression = 10",
        "pull-coord1-dx = 0", // INVALID
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, MissingTransformationCoordExpression)
{
    const char* inputMdpFile[] = {
        "pull = yes",
        "pull-ncoords = 1",
        "pull-coord1-geometry = transformation",
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_SC_And_ExactlyAsManyStep)
{
    // 1e5 steps and delta lambda 1e-5, should not warn (exactly right)
    const char* inputMdpFile[] = {
        "nsteps        = 100000",
        "nstdhdl       = 1",
        "nstcalcenergy       = 1",
        "free-energy   = yes",
        "init-lambda   = 0",
        "delta-lambda  = 1e-05",
        "sc-alpha      = 0.3",
        "sc-sigma      = 0.25",
        "; decoupled VdW to avoid warning about",
        "; - not the point of this test",
        "couple-lambda0 = q",
        "sc-power      = 1",
        "sc-coul       = yes",

    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::NoErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_SC_And_ExactlyAsManyStep_negativeDelta)
{
    // 1e5 steps and delta lambda -1e-5, should not warn (exactly right)
    const char* inputMdpFile[] = {
        "nsteps        = 100000",
        "nstdhdl       = 1",
        "nstcalcenergy       = 1",
        "free-energy   = yes",
        "init-lambda   = 1.0",
        "delta-lambda  = -1e-05",
        "sc-alpha      = 0.3",
        "sc-sigma      = 0.25",
        "; decoupled VdW to avoid warning about",
        "; - not the point of this test",
        "couple-lambda0 = q",
        "sc-power      = 1",
        "sc-coul       = yes",

    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::NoErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_NoSC_And_ExactlyAsManyStep)
{
    // 1e5 steps and delta lambda 1e-5, should not warn (exactly right)
    // Should not warn without softcore too
    const char* inputMdpFile[] = { "nsteps        = 100000",  "nstdhdl       = 1",
                                   "nstcalcenergy       = 1", "free-energy   = yes",
                                   "init-lambda   = 0",       "delta-lambda  = 1e-05" };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::NoErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_NoSC_And_ExactlyAsManyStep_negativeDelta)
{
    // 1e5 steps and delta lambda -1e-5, should not warn (exactly right)
    // Should not warn without softcore too
    const char* inputMdpFile[] = { "nsteps        = 100000",  "nstdhdl       = 1",
                                   "nstcalcenergy       = 1", "free-energy   = yes",
                                   "init-lambda   = 1",       "delta-lambda  = -1e-05" };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::NoErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_SC_And_OneStepTooMuch)
{
    // 1e5+1 steps and delta lambda 1e-5, should error (will be over 1 in the very last step
    // and this is unsupported by softcore)
    const char* inputMdpFile[] = {
        "nsteps        = 100001",
        "nstdhdl       = 1",
        "nstcalcenergy       = 1",
        "free-energy   = yes",
        "init-lambda   = 0",
        "delta-lambda  = 1e-05",
        "sc-alpha      = 0.3",
        "sc-sigma      = 0.25",
        "; decoupled VdW to avoid warning about",
        "; - not the point of this test",
        "couple-lambda0 = q",
        "sc-power      = 1",
        "sc-coul       = yes",
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_SC_And_OneStepTooMuch_negativeDelta)
{
    // 1e5+1 steps and delta lambda -1e-5, should error (will be under 0 in the very last step
    // and this is unsupported by softcore)
    const char* inputMdpFile[] = {
        "nsteps        = 100001",
        "nstdhdl       = 1",
        "nstcalcenergy       = 1",
        "free-energy   = yes",
        "init-lambda   = 1",
        "delta-lambda  = -1e-05",
        "sc-alpha      = 0.3",
        "sc-sigma      = 0.25",
        "; decoupled VdW to avoid warning about",
        "; - not the point of this test",
        "couple-lambda0 = q",
        "sc-power      = 1",
        "sc-coul       = yes",
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_NoSC_And_OneStepTooMuch_negativeDelta)
{
    // 1e5+1 steps and delta lambda -1e-5, should warn (will be under 0 in the very last step)
    // Without softcore, this is still a warning
    const char* inputMdpFile[] = {
        "nsteps        = 100001", "nstdhdl       = 1", "nstcalcenergy       = 1",
        "free-energy   = yes",    "init-lambda   = 1", "delta-lambda  = -1e-05",
    };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_LambdaVector_And_OneStepTooMuch)
{
    // 1e5+1 steps and delta lambda 1e-5, with lambda vector, should warn (will be capped in the very last step)
    const char* inputMdpFile[] = { "nsteps        = 100001",  "nstcalcenergy       = 1",
                                   "nstdhdl       = 1",       "free-energy   = yes",
                                   "delta-lambda  = 1e-05",   "init-lambda-state = 0",
                                   "fep_lambdas =  0 0.5 1.0" };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_LambdaVector_And_OneStepTooMuch_negativeDelta)
{
    // 1e5+1 steps and delta lambda -1e-5, with lambda vector, should warn (will be capped in the very last step)
    const char* inputMdpFile[] = { "nsteps        = 100001",  "nstdhdl       = 1",
                                   "nstcalcenergy       = 1", "free-energy   = yes",
                                   "delta-lambda  = -1e-05",  "init-lambda-state = 2",
                                   "fep_lambdas =  0 0.5 1.0" };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::ErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_LambdaVector_And_ExactlyAsManyStep)
{
    // 1e5 steps and delta lambda 1e-5, with lambda vector, should not warn of lambda capping
    const char* inputMdpFile[] = { "nsteps        = 100000",  "nstdhdl       = 1",
                                   "nstcalcenergy       = 1", "free-energy   = yes",
                                   "delta-lambda  = 1e-05",   "init-lambda-state = 0",
                                   "fep_lambdas =  0 0.5 1.0" };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::NoErrorAndDoNotCompareOutput);
}

TEST_F(GetIrTest, lambdaOverOneCheck_LambdaVector_And_ExactlyAsManyStep_negativeDelta)
{
    // 1e5 steps and delta lambda -1e-5, with lambda vector, should not warn of lambda capping
    const char* inputMdpFile[] = { "nsteps        = 100000",  "nstdhdl       = 1",
                                   "nstcalcenergy       = 1", "free-energy   = yes",
                                   "delta-lambda  = -1e-05",  "init-lambda-state = 2",
                                   "fep_lambdas =  0 0.5 1.0" };
    runTest(joinStrings(inputMdpFile, "\n"), TestBehavior::NoErrorAndDoNotCompareOutput);
}

#endif // HAVE_MUPARSER

} // namespace test
} // namespace gmx
