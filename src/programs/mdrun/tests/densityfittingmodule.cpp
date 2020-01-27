/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Tests utilities for "Densityfitting" setups.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
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
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(".edr");
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
            TestFileManager::getInputFilePath("ellipsoid-density.mrc").c_str());

    //! Mdp values for md integrator with default density fitting parameters.
    const std::string mdpMdDensfitYesUnsetValues = formatString(
            "integrator                       = md\n"
            "nsteps                           = 2\n"
            "cutoff-scheme                    = verlet\n"
            "density-guided-simulation-active = yes\n"
            "density-guided-simulation-group  = FirstThreeOfTwelve\n"
            "density-guided-simulation-reference-density-filename = %s\n",
            TestFileManager::getInputFilePath("ellipsoid-density.mrc").c_str());

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

    EXPECT_DEATH_IF_SUPPORTED(runner_.callGrompp(),
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
    runner_.cptFileName_ = fileManager_.getTemporaryFilePath(".cpt");
    commandLineForMdrun_.addOption("-cpo", runner_.cptFileName_);

    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun(commandLineForMdrun_));

    // checkMdrun(expectedEnergyTermMagnitude);

    CommandLine commandLineForRestart;
    commandLineForRestart.addOption("-cpi", runner_.cptFileName_);
    commandLineForRestart.addOption("-noappend");
    runner_.nsteps_ = 4;
    ASSERT_EQ(0, runner_.callMdrun(commandLineForRestart));

    const real expectedEnergyTermMagnitude = -3378.825928;
    checkMdrun(expectedEnergyTermMagnitude);
}


} // namespace test
} // namespace gmx
