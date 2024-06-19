/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests utilities for "Densityfitting" setups.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <filesystem>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{

/*! \internal
 * \brief End to end test for the density fitting functionality in mdrun.
 *
 * This test checks whether the density fitting related .mdp input options are
 * understood by grompp and that the density fitting module in mdrun computes
 * the correct density fitting energies and stores them in the .edr output file.
 */
class DensityFittingTest : public MdrunTestFixture
{
public:
    DensityFittingTest()
    {
        runner_.useTopGroAndNdxFromDatabase("argon12");
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(".edr").string();
    }

    //! Check the output of mdrun
    void checkMdrun(real energyTermMagnitude)
    {
        const FloatingPointTolerance energyTermTolerance =
                relativeToleranceAsFloatingPoint(energyTermMagnitude, 1e-4);

        EnergyTermsToCompare energyTermsToCompare{
            { { interaction_function[F_DENSITYFITTING].longname, energyTermTolerance },
              { interaction_function[F_EPOT].longname, energyTermTolerance } }
        };

        TestReferenceData refData;
        auto              checker = refData.rootChecker();
        checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
    }

    //! Mdp values for steepest-decent energy minimization with default density fitting parameters.
    const std::string mdpEminDensfitYesUnsetValues = formatString(
            "integrator                       = steep\n"
            "nsteps                           = 2\n"
            "cutoff-scheme                    = verlet\n"
            "density-guided-simulation-active = yes\n"
            "density-guided-simulation-group  = FirstThreeOfTwelve\n"
            "density-guided-simulation-reference-density-filename = %s\n",
            TestFileManager::getInputFilePath("ellipsoid-density.mrc").string().c_str());

    //! Mdp values for md integrator with default density fitting parameters.
    const std::string mdpMdDensfitYesUnsetValues = formatString(
            "integrator                       = md\n"
            "nsteps                           = 2\n"
            "cutoff-scheme                    = verlet\n"
            "density-guided-simulation-active = yes\n"
            "density-guided-simulation-group  = FirstThreeOfTwelve\n"
            "density-guided-simulation-reference-density-filename = %s\n",
            TestFileManager::getInputFilePath("ellipsoid-density.mrc").string().c_str());

    //! Mdp values for steepest-decent energy minimization with density fitting values set to non-defaults.
    const std::string mdpDensiftAllDefaultsChanged_ = formatString(
            "density-guided-simulation-similarity-measure = relative-entropy\n"
            "density-guided-simulation-atom-spreading-weight = mass\n"
            "density-guided-simulation-force-constant = -1\n"
            "density-guided-simulation-gaussian-transform-spreading-width = 0.8\n"
            "density-guided-simulation-gaussian-transform-spreading-range-in-multiples-of-width = "
            "6\n"
            "density-guided-simulation-normalize-densities = false\n");
    //! Set mdp values so that energy calculation interval and density guided simulation interval mismatch.
    const std::string mdpEnergyAndDensityfittingIntervalMismatch_ = formatString(
            "nstcalcenergy = 7\n"
            "density-guided-simulation-nst = 3\n");
    //! Set mdp values so that we skip every other step
    const std::string mdpSkipDensityfittingEveryOtherStep_ = formatString(
            "nstenergy = 2\n"
            "density-guided-simulation-nst = 2\n");
    //! A properly set shift vector
    const std::string mdpTranslationSet_ =
            formatString("density-guided-simulation-shift-vector = 0.1 -0.2 0.3\n");
    //! A shift vector that is lacking an entry
    const std::string mdpTranslationSetWrongValues_ =
            formatString("density-guided-simulation-shift-vector = 0.1 -0.2\n");
    //! A 45 degree rotation around the y axis expressed as matrix transformation
    const std::string mdpTransformationMatrix1degAroundY_ = formatString(
            "density-guided-simulation-transformation-matrix = 0.9998477 0.0000000 0.0174524 "
            "0.0000000 1.0000000 0.0000000 -0.0174524 0.0000000 0.9998477 \n");
    //! The identity matrix as transformation matrix
    const std::string mdpTransformationMatrixIdentity_ = formatString(
            "density-guided-simulation-transformation-matrix = 1 0 0 "
            "0 1 0 0 0 1 \n");
    //! A transformation matrix string where only eight values are given
    const std::string mdpTransformationMatrixWrongValues_ = formatString(
            "density-guided-simulation-transformation-matrix = 0.7071068 0.0000000 0.7071068 "
            "0.0000000 0.0000000 -0.7071068 0.0000000 0.7071068 \n");

    //! The command line to call mdrun
    CommandLine commandLineForMdrun_;
};

/* Fit a subset of three of twelve argon atoms into a reference density
 * whose origin is offset from the simulation box origin.
 *
 * All density fitting mdp parameters are set to defaults
 */
TEST_F(DensityFittingTest, EnergyMinimizationEnergyCorrectInnerProduct)
{
    runner_.useStringAsMdpFile(mdpEminDensfitYesUnsetValues);

    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun(commandLineForMdrun_));

    const real expectedEnergyTermMagnitude = -3378.825928;
    checkMdrun(expectedEnergyTermMagnitude);
}

/* Fit a subset of three of twelve argon atoms into a reference density
 * whose origin is offset from the simulation box origin.
 *
 * All density fitting mdp parameters are set to defaults
 */
TEST_F(DensityFittingTest, EnergyMinimizationEnergyCorrectInnerProductTranslation)
{
    runner_.useStringAsMdpFile(mdpEminDensfitYesUnsetValues + mdpTranslationSet_);

    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun(commandLineForMdrun_));

    const real expectedEnergyTermMagnitude = -8991;
    checkMdrun(expectedEnergyTermMagnitude);
}

/* Fit a subset of three of twelve argon atoms into a reference density
 * whose origin is offset from the simulation box origin.
 *
 * All density fitting mdp parameters are set to defaults
 */
TEST_F(DensityFittingTest, EnergyMinimizationEnergyTranslationParametersOff)
{
    runner_.useStringAsMdpFile(mdpEminDensfitYesUnsetValues + mdpTranslationSetWrongValues_);

    GMX_EXPECT_DEATH_IF_SUPPORTED(runner_.callGrompp(), ".*Reading three real values.*");
}

/* Fit a subset of three of twelve argon atoms into a reference density
 * that are rotated around the simulation box origin by a matrix multiplication.
 *
 * All density fitting mdp parameters are set to defaults
 */
TEST_F(DensityFittingTest, EnergyMinimizationEnergyCorrectInnerProductTranslationAndTransformationMatrix)
{
    runner_.useStringAsMdpFile(mdpEminDensfitYesUnsetValues + mdpTranslationSet_
                               + mdpTransformationMatrix1degAroundY_);

    runner_.nsteps_ = 4;
    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun(commandLineForMdrun_));

    const real expectedEnergyTermMagnitude = -8991;
    checkMdrun(expectedEnergyTermMagnitude);
}

/* Fit a subset of three of twelve argon atoms into a reference density
 * whose origin is offset from the simulation box origin.
 *
 * All density fitting mdp parameters are set to defaults
 */
TEST_F(DensityFittingTest, EnergyMinimizationEnergyMatrixTransfromationOff)
{
    runner_.useStringAsMdpFile(mdpEminDensfitYesUnsetValues + mdpTransformationMatrixWrongValues_);

    GMX_EXPECT_DEATH_IF_SUPPORTED(runner_.callGrompp(), ".*Reading nine real values.*");
}

/* Fit a subset of three of twelve argon atoms into a reference density
 * where the given matrix transformation is the identity transformation.
 *
 * All density fitting mdp parameters are set to defaults
 */
TEST_F(DensityFittingTest, EnergyMinimizationEnergyCorrectInnerProductIdentityMatrix)
{
    runner_.useStringAsMdpFile(mdpEminDensfitYesUnsetValues + mdpTransformationMatrixIdentity_);

    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun(commandLineForMdrun_));

    const real expectedEnergyTermMagnitude = -8991;
    checkMdrun(expectedEnergyTermMagnitude);
}

/* Like above, but with as many parameters reversed as possible
 */
TEST_F(DensityFittingTest, EnergyMinimizationEnergyCorrectForRelativeEntropy)
{
    runner_.useStringAsMdpFile(mdpEminDensfitYesUnsetValues + mdpDensiftAllDefaultsChanged_);

    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun(commandLineForMdrun_));

    const real expectedEnergyTermMagnitude = -27000;
    checkMdrun(expectedEnergyTermMagnitude);
}

/* Test that grompp exits with error message if energy evaluation frequencies
 * do not match.
 */
TEST_F(DensityFittingTest, GromppErrorWhenEnergyEvaluationFrequencyMismatch)
{
    runner_.useStringAsMdpFile(mdpMdDensfitYesUnsetValues + mdpEnergyAndDensityfittingIntervalMismatch_);

    GMX_EXPECT_DEATH_IF_SUPPORTED(runner_.callGrompp(),
                                  ".*is not a multiple of density-guided-simulation-nst.*");
}

/* Fit a subset of three of twelve argon atoms into a reference density
 * whose origin is offset from the simulation box origin. Stop the simulation,
 * then restart.
 *
 * All density fitting mdp parameters are set to defaults
 */
TEST_F(DensityFittingTest, CheckpointWorks)
{
    runner_.useStringAsMdpFile(mdpMdDensfitYesUnsetValues + mdpSkipDensityfittingEveryOtherStep_);

    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun(commandLineForMdrun_));

    // checkMdrun(expectedEnergyTermMagnitude);

    CommandLine commandLineForRestart;
    commandLineForRestart.addOption("-cpi", runner_.cptOutputFileName_);
    commandLineForRestart.addOption("-noappend");
    runner_.nsteps_ = 4;
    ASSERT_EQ(0, runner_.callMdrun(commandLineForRestart));

    const real expectedEnergyTermMagnitude = -3378.825928;
    checkMdrun(expectedEnergyTermMagnitude);
}


} // namespace test
} // namespace gmx
