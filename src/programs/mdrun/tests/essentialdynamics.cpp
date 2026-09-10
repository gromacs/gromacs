/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for mdrun essential dynamics (ED) and flooding functionality.
 *
 * These tests check the mdrun results for cases where essential dynamics
 * sampling is enabled. Essential dynamics is a method to enhance sampling
 * along specific degrees of freedom (eigenvectors from PCA) or to apply
 * flooding potentials.
 *
 * Most tests use a minimal 4-atom system (butane-like) in a 10 nm box to verify:
 * - Linear acceleration (linacc): Accelerate along an ED vector
 * - Linear fixed (linfix): Keep system fixed along an ED vector
 * - Radial acceleration (radacc): Radial acceleration in ED space
 * - Radial constraint (radcon): Constrain radial distance in ED space
 * - Radial fixed (radfix): Keep system at fixed radius in ED space
 *
 * Flooding tests use the original 4076-atom protein system (guanylin in water):
 * - Flooding1: Flooding potential initialization (0 steps)
 * - Flooding2: Flooding potential evolution (50 steps)
 *
 * All tests run on a single rank (4-atom system too small for domain decomposition,
 * 4076-atom flooding system matches original regressiontests single-rank behavior).
 *
 * EDI files were generated using gmx make_edi with eigenvectors from
 * essentialdynamics-eigenvec.trr (computed via PCA on the 4-atom system).
 * The EDI files match those validated in the original regressiontests repository.
 *
 * Each test verifies that:
 * 1. The simulation runs without errors
 * 2. Energy terms are within expected tolerances
 * 3. EDSAM output (eigenvector projections) matches reference data
 *
 * Note: Trajectory coordinates are NOT checked because ED simulations are
 * chaotic - tiny numerical differences lead to completely divergent trajectories.
 *
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <filesystem>
#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/xvgtest.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"
#include "trajectorycomparison.h"

namespace gmx
{
namespace test
{
namespace
{

//! Test parameters: ED type and simulation steps
typedef std::tuple<std::string, int> EDParameters;

/*! \brief Test fixture for essential dynamics simulations
 *
 * This test runs essential dynamics simulations with various ED modes
 * and verifies energies and EDSAM output against reference data.
 */
class EssentialDynamicsTest : public MdrunTestFixture, public ::testing::WithParamInterface<EDParameters>
{
};

//! MDP settings common to all ED tests (without nsteps, set per test)
const std::string g_commonMdpSettings = R"(
integrator               = md
dt                       = 0.001
comm-mode                = linear
nstcomm                  = 10
ld-seed                  = 1993
cutoff-scheme            = Verlet
coulombtype              = Cut-off
coulomb-modifier         = None
rcoulomb                 = 2.5
vdwtype                  = Cut-off
vdw-modifier             = None
rvdw                     = 2.5
rlist                    = 2.5
verlet-buffer-tolerance  = -1
nstlist                  = 10
nstxout                  = 2
nstvout                  = 0
nstfout                  = 0
nstlog                   = 10
nstenergy                = 2
nstcalcenergy            = 2
pbc                      = xyz
)";

//! Formatter for int parameters (converts to string)
static std::string intToString(int value)
{
    return std::to_string(value);
}

//! Tuple of formatters to name the parameterized test cases (e.g., "linacc_50")
const NameOfTestFromTuple<EDParameters> sc_testNamer{ std::make_tuple(useString, intToString) };

//! Helper to make reference data filenames match test names (includes both ED type and nsteps)
const RefDataFilenameMaker<EDParameters> sc_refDataFilenameMaker{ std::make_tuple(useString, intToString) };

TEST_P(EssentialDynamicsTest, WithinTolerances)
{
    GMX_MPI_TEST(RequireRankCount<1>);

    auto params = GetParam();
    auto edType = std::get<0>(params);
    auto nsteps = std::get<1>(params);

    SCOPED_TRACE(formatString("Testing essential dynamics mode '%s' with %d steps", edType.c_str(), nsteps));

    // Flooding tests use a different system (4076-atom protein) and MDP settings
    bool isFloodingTest = (edType == "flooding1" || edType == "flooding2");
    std::string systemBaseName = isFloodingTest ? "essentialdynamics-flooding" : "essentialdynamics";
    std::string mdpBaseName = isFloodingTest ? "essentialdynamics-flooding" : "essentialdynamics";

    // Load MDP file from database and replace nsteps value
    auto mdpPath = gmx::test::TestFileManager::getTestSimulationDatabaseDirectory()
                   / (mdpBaseName + ".mdp");
    std::string mdpContents = TextReader::readFileToString(mdpPath.string());
    // Replace nsteps line with the value needed for this test
    size_t nstepsPos = mdpContents.find("nsteps");
    if (nstepsPos != std::string::npos)
    {
        size_t lineEnd = mdpContents.find('\n', nstepsPos);
        mdpContents.replace(nstepsPos, lineEnd - nstepsPos, formatString("nsteps = %d", nsteps));
    }
    runner_.useStringAsMdpFile(mdpContents);

    // Set up the coordinate and topology files
    // Flooding tests use a 4076-atom protein system, other ED tests use 4-atom system
    runner_.useTopGroAndNdxFromDatabase(systemBaseName);

    // Construct path to EDI file in the simulation database
    std::filesystem::path ediPath = gmx::test::TestFileManager::getTestSimulationDatabaseDirectory()
                                    / formatString("essentialdynamics-%s.edi", edType.c_str());

    // Call grompp with -maxwarn to allow the Berendsen thermostat/barostat warnings.
    // The original MDP files use Berendsen, which is deprecated but was the
    // standard when these tests were created. Flooding tests have both thermostat
    // and barostat warnings (maxwarn 2), while other ED tests only have the
    // thermostat warning (maxwarn 1).
    CommandLine gromppCommandLine;
    int         maxwarn = isFloodingTest ? 2 : 1;
    gromppCommandLine.addOption("-maxwarn", std::to_string(maxwarn));
    ASSERT_EQ(0, runner_.callGrompp(gromppCommandLine));

    // Set up mdrun to use the EDI file and write EDSAM output to a temporary file
    std::filesystem::path edsamFileName = fileManager_.getTemporaryFilePath("edsam.xvg");
    CommandLine           mdrunCommandLine;
    mdrunCommandLine.addOption("-ei", ediPath.string());
    mdrunCommandLine.addOption("-eo", edsamFileName.string());
    ASSERT_EQ(0, runner_.callMdrun(mdrunCommandLine));

    // Compare results to reference data using descriptive filename based on ED type
    TestReferenceData refData(sc_refDataFilenameMaker(GetParam()));
    auto              checker = refData.rootChecker();

    // Check energies with very loose tolerances appropriate for chaotic ED trajectories.
    // ED simulations are highly sensitive to numerical precision - small differences in
    // compiler optimizations, FFT libraries, or floating-point operations lead to
    // divergent trajectories and energy values. We use 1% relative tolerance to
    // accommodate these differences across build configurations while still catching
    // major regressions.
    {
        auto energyTolerance = relativeToleranceAsFloatingPoint(100.0, 0.01);

        EnergyTermsToCompare energyTermsToCompare{
            { { interaction_function[InteractionFunction::PotentialEnergy].longname, energyTolerance },
              { interaction_function[InteractionFunction::KineticEnergy].longname, energyTolerance },
              { interaction_function[InteractionFunction::TotalEnergy].longname, energyTolerance } }
        };
        checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
    }

    // Check EDSAM output file (eigenvector projections, flooding potentials, etc.)
    // This contains the ED-specific data that the original regressiontests validated
    {
        TextInputFile    edsamFile(edsamFileName.string());
        XvgMatchSettings settings;
        // Use same loose tolerance as for energies due to chaotic nature of ED
        settings.tolerance = relativeToleranceAsFloatingPoint(100.0, 0.01);
        checkXvgFile(&edsamFile, &checker, settings);
    }
}

// All tests run on single rank (small systems)
INSTANTIATE_TEST_SUITE_P(EssentialDynamicsTests,
                         EssentialDynamicsTest,
                         ::testing::Values(EDParameters{ "linacc", 50 },
                                           EDParameters{ "linfix", 50 },
                                           EDParameters{ "radacc", 50 },
                                           EDParameters{ "radcon", 50 },
                                           EDParameters{ "radfix", 50 },
                                           EDParameters{ "flooding1", 0 },
                                           EDParameters{ "flooding2", 50 }),
                         sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
