/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * Tests AWH runs
 *
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
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/xvgtest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

class AwhTest : public MdrunTestFixture
{
};

namespace
{

const FloatingPointTolerance pullXTolerance = relativeToleranceAsFloatingPoint(100.0, 1e-4);
const FloatingPointTolerance pullFTolerance = relativeToleranceAsFloatingPoint(100.0, 1e-3);

//! Check xvg file contents against reference data
void checkXvgOutputFile(TestReferenceChecker*         checker,
                        const std::filesystem::path&  filename,
                        const char*                   description,
                        const FloatingPointTolerance& tolerance)
{
    if (File::exists(filename, File::returnFalseOnError))
    {
        TestReferenceChecker fileChecker(checker->checkCompound("File", description));
        TextInputFile        xvgFileStream(filename);
        XvgMatchSettings     matchSettings;
        matchSettings.tolerance = tolerance;
        checkXvgFile(&xvgFileStream, &fileChecker, matchSettings);
    }
}

} // namespace

// Tests biasing of the two-dimensional dihedral angle coordinate
// (phi, psi) for alanine dipeptide in vacuum. This also tests the
// periodic AWH grid.  Note that this test can be unstable. It
// compares different bias variables on the grid which is in this case
// 2D. As MD is chaotic, it can happen rather soon that runs with
// different setup, e.g. number of MPI ranks, visit at least one
// different grid point. Therefore this test is over only 20 steps.

TEST_F(AwhTest, SingleBiasMultiDimCanRun)
{
    runner_.useTopGroAndNdxFromDatabase("alanine_vacuo");
    MdpFieldValues mdpFieldValues = prepareMdpFieldValues(
            "alanine_vacuo", "md", "v-rescale", "No", MdpParameterDatabase::Default);
    mdpFieldValues["nsteps"]    = "20";
    mdpFieldValues["nstenergy"] = "10";

    mdpFieldValues["pull"]                           = "yes";
    mdpFieldValues["pull-cylinder-r"]                = "1.5";
    mdpFieldValues["pull-constr-tol"]                = "1e-06";
    mdpFieldValues["pull-print-com"]                 = "no";
    mdpFieldValues["pull-print-ref-value"]           = "no";
    mdpFieldValues["pull-print-components"]          = "no";
    mdpFieldValues["pull-nstxout"]                   = "20";
    mdpFieldValues["pull-nstfout"]                   = "20";
    mdpFieldValues["pull-ngroups"]                   = "5";
    mdpFieldValues["pull-ncoords"]                   = "2";
    mdpFieldValues["pull-group1-name"]               = "C_&_r_1";
    mdpFieldValues["pull-group1-pbcatom"]            = "0";
    mdpFieldValues["pull-group2-name"]               = "N_&_r_2";
    mdpFieldValues["pull-group2-pbcatom"]            = "0";
    mdpFieldValues["pull-group3-name"]               = "CA";
    mdpFieldValues["pull-group3-pbcatom"]            = "0";
    mdpFieldValues["pull-group4-name"]               = "C_&_r_2";
    mdpFieldValues["pull-group4-pbcatom"]            = "0";
    mdpFieldValues["pull-group5-name"]               = "N_&_r_3";
    mdpFieldValues["pull-group5-pbcatom"]            = "0";
    mdpFieldValues["pull-coord1-type"]               = "external-potential";
    mdpFieldValues["pull-coord1-potential-provider"] = "awh";
    mdpFieldValues["pull-coord1-geometry"]           = "dihedral";
    mdpFieldValues["pull-coord1-groups"]             = "1 2 2 3 3 4";
    mdpFieldValues["pull-coord1-dim"]                = "Y Y Y";
    mdpFieldValues["pull-coord1-origin"]             = "0.0 0.0 0.0";
    mdpFieldValues["pull-coord1-vec"]                = "0.0 0.0 0.0";
    mdpFieldValues["pull-coord1-start"]              = "no";
    mdpFieldValues["pull-coord1-init"]               = "0";
    mdpFieldValues["pull-coord1-rate"]               = "0";
    mdpFieldValues["pull-coord1-k"]                  = "4000";
    mdpFieldValues["pull-coord1-kB"]                 = "1000";
    mdpFieldValues["pull-coord2-type"]               = "external-potential";
    mdpFieldValues["pull-coord2-potential-provider"] = "awh";
    mdpFieldValues["pull-coord2-geometry"]           = "dihedral";
    mdpFieldValues["pull-coord2-groups"]             = "2 3 3 4 4 5";
    mdpFieldValues["pull-coord2-dim"]                = "Y Y Y";
    mdpFieldValues["pull-coord2-origin"]             = "0.0 0.0 0.0";
    mdpFieldValues["pull-coord2-vec"]                = "0.0 0.0 0.0";
    mdpFieldValues["pull-coord2-start"]              = "no";
    mdpFieldValues["pull-coord2-init"]               = "0";
    mdpFieldValues["pull-coord2-rate"]               = "0";
    mdpFieldValues["pull-coord2-k"]                  = "4000";
    mdpFieldValues["pull-coord2-kB"]                 = "1000";

    mdpFieldValues["awh"]                        = "yes";
    mdpFieldValues["awh-potential"]              = "convolved";
    mdpFieldValues["awh-seed"]                   = "955431239";
    mdpFieldValues["awh-nstout"]                 = "10";
    mdpFieldValues["awh-nstsample"]              = "10";
    mdpFieldValues["awh-nsamples-update"]        = "1";
    mdpFieldValues["awh-share-multisim"]         = "no";
    mdpFieldValues["awh-nbias"]                  = "1";
    mdpFieldValues["awh1-error-init"]            = "10";
    mdpFieldValues["awh1-growth"]                = "exp-linear";
    mdpFieldValues["awh1-growth-factor"]         = "3";
    mdpFieldValues["awh1-equilibrate-histogram"] = "no";
    mdpFieldValues["awh1-target"]                = "constant";
    mdpFieldValues["awh1-target-beta-scaling"]   = "0";
    mdpFieldValues["awh1-target-cutoff"]         = "0";
    mdpFieldValues["awh1-user-data"]             = "no";
    mdpFieldValues["awh1-share-group"]           = "0";
    mdpFieldValues["awh1-ndim"]                  = "2";
    mdpFieldValues["awh1-dim1-coord-provider"]   = "pull";
    mdpFieldValues["awh1-dim1-coord-index"]      = "2";
    mdpFieldValues["awh1-dim1-start"]            = "150";
    mdpFieldValues["awh1-dim1-end"]              = "180";
    mdpFieldValues["awh1-dim1-force-constant"]   = "4000";
    mdpFieldValues["awh1-dim1-diffusion"]        = "0.1";
    mdpFieldValues["awh1-dim1-cover-diameter"]   = "0";
    mdpFieldValues["awh1-dim2-coord-provider"]   = "pull";
    mdpFieldValues["awh1-dim2-coord-index"]      = "1";
    mdpFieldValues["awh1-dim2-start"]            = "178";
    mdpFieldValues["awh1-dim2-end"]              = "-178";
    mdpFieldValues["awh1-dim2-force-constant"]   = "4000";
    mdpFieldValues["awh1-dim2-diffusion"]        = "0.1";
    mdpFieldValues["awh1-dim2-cover-diameter"]   = "0";

    const std::string mdpContents = prepareMdpFileContents(mdpFieldValues);
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    ::gmx::test::CommandLine awhCaller;

    // Do an initial mdrun that writes a checkpoint file
    ::gmx::test::CommandLine firstCaller(awhCaller);
    auto                     pullXFileName = fileManager_.getTemporaryFilePath("pullx.xvg");
    auto                     pullFFileName = fileManager_.getTemporaryFilePath("pullf.xvg");
    firstCaller.addOption("-px", pullXFileName);
    firstCaller.addOption("-pf", pullFFileName);
    ASSERT_EQ(0, runner_.callMdrun(firstCaller));

    TestReferenceData    data;
    TestReferenceChecker checker = data.rootChecker();

    // Test the AWH output files
    checkXvgOutputFile(&checker, pullXFileName, "AWH pull x output after 20 steps", pullXTolerance);
    checkXvgOutputFile(&checker, pullFFileName, "AWH pull f output after 20 steps", pullFTolerance);
}

// Tests biasing of two one-dimensional dihedral angle coordinates
// (phi and psi) for alanine dipeptide in vacuum. This also tests the
// periodic AWH grid. Furthermore, it checks checkpoint writing and
// continuation, which are the major components of the AWH module that
// are not covered by module tests.

TEST_F(AwhTest, MultiBiasCanRunAndRestartFromCheckpoint)
{
    runner_.useTopGroAndNdxFromDatabase("alanine_vacuo");
    MdpFieldValues mdpFieldValues = prepareMdpFieldValues(
            "alanine_vacuo", "md", "v-rescale", "No", MdpParameterDatabase::Default);
    mdpFieldValues["nsteps"]    = "40";
    mdpFieldValues["nstenergy"] = "10";

    mdpFieldValues["pull"]                           = "yes";
    mdpFieldValues["pull-cylinder-r"]                = "1.5";
    mdpFieldValues["pull-constr-tol"]                = "1e-06";
    mdpFieldValues["pull-print-com"]                 = "no";
    mdpFieldValues["pull-print-ref-value"]           = "no";
    mdpFieldValues["pull-print-components"]          = "no";
    mdpFieldValues["pull-nstxout"]                   = "20";
    mdpFieldValues["pull-nstfout"]                   = "20";
    mdpFieldValues["pull-ngroups"]                   = "5";
    mdpFieldValues["pull-ncoords"]                   = "2";
    mdpFieldValues["pull-group1-name"]               = "C_&_r_1";
    mdpFieldValues["pull-group1-pbcatom"]            = "0";
    mdpFieldValues["pull-group2-name"]               = "N_&_r_2";
    mdpFieldValues["pull-group2-pbcatom"]            = "0";
    mdpFieldValues["pull-group3-name"]               = "CA";
    mdpFieldValues["pull-group3-pbcatom"]            = "0";
    mdpFieldValues["pull-group4-name"]               = "C_&_r_2";
    mdpFieldValues["pull-group4-pbcatom"]            = "0";
    mdpFieldValues["pull-group5-name"]               = "N_&_r_3";
    mdpFieldValues["pull-group5-pbcatom"]            = "0";
    mdpFieldValues["pull-coord1-type"]               = "external-potential";
    mdpFieldValues["pull-coord1-potential-provider"] = "awh";
    mdpFieldValues["pull-coord1-geometry"]           = "dihedral";
    mdpFieldValues["pull-coord1-groups"]             = "1 2 2 3 3 4";
    mdpFieldValues["pull-coord1-dim"]                = "Y Y Y";
    mdpFieldValues["pull-coord1-origin"]             = "0.0 0.0 0.0";
    mdpFieldValues["pull-coord1-vec"]                = "0.0 0.0 0.0";
    mdpFieldValues["pull-coord1-start"]              = "no";
    mdpFieldValues["pull-coord1-init"]               = "0";
    mdpFieldValues["pull-coord1-rate"]               = "0";
    mdpFieldValues["pull-coord1-k"]                  = "4000";
    mdpFieldValues["pull-coord1-kB"]                 = "1000";
    mdpFieldValues["pull-coord2-type"]               = "external-potential";
    mdpFieldValues["pull-coord2-potential-provider"] = "awh";
    mdpFieldValues["pull-coord2-geometry"]           = "dihedral";
    mdpFieldValues["pull-coord2-groups"]             = "2 3 3 4 4 5";
    mdpFieldValues["pull-coord2-dim"]                = "Y Y Y";
    mdpFieldValues["pull-coord2-origin"]             = "0.0 0.0 0.0";
    mdpFieldValues["pull-coord2-vec"]                = "0.0 0.0 0.0";
    mdpFieldValues["pull-coord2-start"]              = "no";
    mdpFieldValues["pull-coord2-init"]               = "0";
    mdpFieldValues["pull-coord2-rate"]               = "0";
    mdpFieldValues["pull-coord2-k"]                  = "4000";
    mdpFieldValues["pull-coord2-kB"]                 = "1000";

    mdpFieldValues["awh"]                        = "yes";
    mdpFieldValues["awh-potential"]              = "convolved";
    mdpFieldValues["awh-seed"]                   = "955431239";
    mdpFieldValues["awh-nstout"]                 = "10";
    mdpFieldValues["awh-nstsample"]              = "10";
    mdpFieldValues["awh-nsamples-update"]        = "1";
    mdpFieldValues["awh-share-multisim"]         = "no";
    mdpFieldValues["awh-nbias"]                  = "2";
    mdpFieldValues["awh1-error-init"]            = "10";
    mdpFieldValues["awh1-growth"]                = "exp-linear";
    mdpFieldValues["awh1-growth-factor"]         = "3";
    mdpFieldValues["awh1-equilibrate-histogram"] = "no";
    mdpFieldValues["awh1-target"]                = "constant";
    mdpFieldValues["awh1-target-beta-scaling"]   = "0";
    mdpFieldValues["awh1-target-cutoff"]         = "0";
    mdpFieldValues["awh1-user-data"]             = "no";
    mdpFieldValues["awh1-share-group"]           = "0";
    mdpFieldValues["awh1-ndim"]                  = "1";
    mdpFieldValues["awh1-dim1-coord-provider"]   = "pull";
    mdpFieldValues["awh1-dim1-coord-index"]      = "2";
    mdpFieldValues["awh1-dim1-start"]            = "150";
    mdpFieldValues["awh1-dim1-end"]              = "180";
    mdpFieldValues["awh1-dim1-force-constant"]   = "4000";
    mdpFieldValues["awh1-dim1-diffusion"]        = "0.1";
    mdpFieldValues["awh1-dim1-cover-diameter"]   = "0";
    mdpFieldValues["awh2-error-init"]            = "10";
    mdpFieldValues["awh2-growth"]                = "exp-linear";
    mdpFieldValues["awh2-growth-factor"]         = "3";
    mdpFieldValues["awh2-equilibrate-histogram"] = "no";
    mdpFieldValues["awh2-target"]                = "constant";
    mdpFieldValues["awh2-target-beta-scaling"]   = "0";
    mdpFieldValues["awh2-target-cutoff"]         = "0";
    mdpFieldValues["awh2-user-data"]             = "no";
    mdpFieldValues["awh2-share-group"]           = "0";
    mdpFieldValues["awh2-ndim"]                  = "1";
    mdpFieldValues["awh2-dim1-coord-provider"]   = "pull";
    mdpFieldValues["awh2-dim1-coord-index"]      = "1";
    mdpFieldValues["awh2-dim1-start"]            = "178";
    mdpFieldValues["awh2-dim1-end"]              = "-178";
    mdpFieldValues["awh2-dim1-force-constant"]   = "4000";
    mdpFieldValues["awh2-dim1-diffusion"]        = "0.1";
    mdpFieldValues["awh2-dim1-cover-diameter"]   = "0";

    const std::string mdpContents = prepareMdpFileContents(mdpFieldValues);
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    ::gmx::test::CommandLine awhCaller;

    // Do an initial mdrun that writes a checkpoint file
    ::gmx::test::CommandLine firstCaller(awhCaller);
    auto                     pullXFileName = fileManager_.getTemporaryFilePath("pullx.xvg");
    auto                     pullFFileName = fileManager_.getTemporaryFilePath("pullf.xvg");
    firstCaller.addOption("-px", pullXFileName);
    firstCaller.addOption("-pf", pullFFileName);
    runner_.nsteps_ = 20;
    ASSERT_EQ(0, runner_.callMdrun(firstCaller));

    TestReferenceData    data;
    TestReferenceChecker checker = data.rootChecker();

    // Test the AWH output files
    checkXvgOutputFile(&checker, pullXFileName, "AWH pull x output after 20 steps", pullXTolerance);
    checkXvgOutputFile(&checker, pullFFileName, "AWH pull f output after 20 steps", pullFTolerance);

#if GMX_LIB_MPI
    // Don't overwrite pullXFileName / pullFFileName before all ranks had the chance to check them
    MPI_Barrier(MdrunTestFixtureBase::s_communicator);
#endif

    // Continue mdrun from that checkpoint file
    ::gmx::test::CommandLine secondCaller(awhCaller);
    secondCaller.addOption("-cpi", runner_.cptOutputFileName_);
    secondCaller.addOption("-px", pullXFileName);
    secondCaller.addOption("-pf", pullFFileName);
    ASSERT_EQ(0, runner_.callMdrun(secondCaller));

    // Test the updated AWH output files
    checkXvgOutputFile(&checker, pullXFileName, "AWH pull x output after 40 steps", pullXTolerance);
    checkXvgOutputFile(&checker, pullFFileName, "AWH pull f output after 40 steps", pullFTolerance);
}

} // namespace test
} // namespace gmx
