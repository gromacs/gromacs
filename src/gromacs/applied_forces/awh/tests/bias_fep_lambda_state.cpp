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
#include "gmxpre.h"

#include <cmath>

#include <memory>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/applied_forces/awh/bias.h"
#include "gromacs/applied_forces/awh/correlationgrid.h"
#include "gromacs/applied_forces/awh/pointstate.h"
#include "gromacs/applied_forces/awh/tests/awh_setup.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

//! The number of lambda states to use in the tests.
static constexpr int c_numLambdaStates = 16;


//! Convenience typedef: growth type enum, potential type enum, disable update skips
typedef std::tuple<AwhHistogramGrowthType, AwhPotentialType, BiasParams::DisableUpdateSkips> BiasTestParameters;

/*! \brief Test fixture for testing Bias updates
 */
class BiasFepLambdaStateTest : public ::testing::TestWithParam<BiasTestParameters>
{
private:
    //! Storage for test parameters.
    std::unique_ptr<AwhTestParameters> params_;

public:
    //! Random seed for AWH MC sampling
    int64_t seed_;

    //! The awh Bias
    std::unique_ptr<Bias> bias_;

    BiasFepLambdaStateTest()
    {
        /* We test all combinations of:
         *   eawhgrowth:
         *     eawhgrowthLINEAR:     final, normal update phase
         *     ewahgrowthEXP_LINEAR: initial phase, updated size is constant
         *   eawhpotential (test only eawhpotentialUMBRELLA (MC) for FEP lambda dimensions)
         *   disableUpdateSkips (should not affect the results):
         *     BiasParams::DisableUpdateSkips::yes: update the point state for every sample
         *     BiasParams::DisableUpdateSkips::no:  update the point state at an interval > 1 sample
         *
         * Note: It would be nice to explicitly check that eawhpotential
         *       and disableUpdateSkips do not affect the point state.
         *       But the reference data will also ensure this.
         */
        AwhHistogramGrowthType         eawhgrowth;
        AwhPotentialType               eawhpotential;
        BiasParams::DisableUpdateSkips disableUpdateSkips;
        std::tie(eawhgrowth, eawhpotential, disableUpdateSkips) = GetParam();

        /* Set up a basic AWH setup with a single, 1D bias with parameters
         * such that we can measure the effects of different parameters.
         */
        constexpr AwhCoordinateProviderType coordinateProvider = AwhCoordinateProviderType::FreeEnergyLambda;
        constexpr int                       coordIndex = 0;
        constexpr double                    origin     = 0;
        constexpr double                    end        = c_numLambdaStates - 1;
        constexpr double                    period     = 0;
        // Correction for removal of GaussianGeometryFactor/2 in histogram size
        constexpr double diffusion = 1e-4 / (0.12927243028700 * 2);
        const auto       awhDimBuffer =
                awhDimParamSerialized(coordinateProvider, coordIndex, origin, end, period, diffusion);
        const auto awhDimArrayRef = gmx::arrayRefFromArray(&awhDimBuffer, 1);
        params_                   = std::make_unique<AwhTestParameters>(getAwhTestParameters(
                eawhgrowth, eawhpotential, awhDimArrayRef, false, 0.4, true, 1.0, c_numLambdaStates));

        seed_ = params_->awhParams.seed();

        double mdTimeStep = 0.1;

        bias_ = std::make_unique<Bias>(-1,
                                       params_->awhParams,
                                       params_->awhParams.awhBiasParams()[0],
                                       params_->dimParams,
                                       params_->beta,
                                       mdTimeStep,
                                       nullptr,
                                       "",
                                       Bias::ThisRankWillDoIO::No,
                                       disableUpdateSkips);
    }
};

TEST_P(BiasFepLambdaStateTest, ForcesBiasPmf)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    Bias& bias = *bias_;

    /* Make strings with the properties we expect to be different in the tests.
     * These also helps to interpret the reference data.
     */
    std::vector<std::string> props;
    props.push_back(formatString("stage:           %s", bias.state().inInitialStage() ? "initial" : "final"));
    props.push_back(formatString("convolve forces: %s", bias.params().convolveForce ? "yes" : "no"));
    props.push_back(formatString("skip updates:    %s", bias.params().skipUpdates() ? "yes" : "no"));

    SCOPED_TRACE(gmx::formatString("%s, %s, %s", props[0].c_str(), props[1].c_str(), props[2].c_str()));

    std::vector<double> force, pot;

    double potentialJump = 0;
    double mdTimeStep    = 0.1;
    int    nSteps        = 501;

    /* Some energies to use as base values (to which some noise is added later on). */
    std::vector<double> neighborLambdaEnergies(c_numLambdaStates);
    std::vector<double> neighborLambdaDhdl(c_numLambdaStates);
    const double        magnitude = 12.0;
    for (int i = 0; i < c_numLambdaStates; i++)
    {
        neighborLambdaEnergies[i] = magnitude * std::sin(i * 0.1);
        neighborLambdaDhdl[i]     = magnitude * std::cos(i * 0.1);
    }

    for (int step = 0; step < nSteps; step++)
    {
        int      umbrellaGridpointIndex = bias.state().coordState().umbrellaGridpoint();
        awh_dvec coordValue = { bias.getGridCoordValue(umbrellaGridpointIndex)[0], 0, 0, 0 };
        double   potential  = 0;
        gmx::ArrayRef<const double> biasForce = bias.calcForceAndUpdateBias(coordValue,
                                                                            neighborLambdaEnergies,
                                                                            neighborLambdaDhdl,
                                                                            &potential,
                                                                            &potentialJump,
                                                                            step * mdTimeStep,
                                                                            step,
                                                                            seed_,
                                                                            nullptr);

        force.push_back(biasForce[0]);
        pot.push_back(potential);
    }

    /* When skipping updates, ensure all skipped updates are performed here.
     * This should result in the same bias state as at output in a normal run.
     */
    if (bias.params().skipUpdates())
    {
        bias.doSkippedUpdatesForAllPoints();
    }

    std::vector<double> pointBias, logPmfsum;
    for (const auto& point : bias.state().points())
    {
        pointBias.push_back(point.bias());
        logPmfsum.push_back(point.logPmfSum());
    }

    constexpr int ulpTol = 10;

    checker.checkSequence(props.begin(), props.end(), "Properties");
    checker.setDefaultTolerance(absoluteTolerance(magnitude * GMX_DOUBLE_EPS * ulpTol));
    checker.checkSequence(force.begin(), force.end(), "Force");
    checker.checkSequence(pot.begin(), pot.end(), "Potential");
    checker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpTol));
    checker.checkSequence(pointBias.begin(), pointBias.end(), "PointBias");
    checker.checkSequence(logPmfsum.begin(), logPmfsum.end(), "PointLogPmfsum");
}

/* Scan initial/final phase, MC/convolved force and update skip (not) allowed
 * Both the convolving and skipping should not affect the bias and PMF.
 * It would be nice if the test would explicitly check for this.
 * Currently this is tested through identical reference data.
 */
INSTANTIATE_TEST_SUITE_P(WithParameters,
                         BiasFepLambdaStateTest,
                         ::testing::Combine(::testing::Values(AwhHistogramGrowthType::Linear,
                                                              AwhHistogramGrowthType::ExponentialLinear),
                                            ::testing::Values(AwhPotentialType::Umbrella),
                                            ::testing::Values(BiasParams::DisableUpdateSkips::yes,
                                                              BiasParams::DisableUpdateSkips::no)));

// Test that we detect coverings and exit the initial stage at the correct step
TEST(BiasFepLambdaStateTest, DetectsCovering)
{
    constexpr AwhCoordinateProviderType coordinateProvider = AwhCoordinateProviderType::FreeEnergyLambda;
    constexpr int                       coordIndex         = 0;
    constexpr double                    origin             = 0;
    constexpr double                    end                = c_numLambdaStates - 1;
    constexpr double                    period             = 0;
    constexpr double                    diffusion          = 1e-4 / (0.12927243028700 * 2);
    auto                                awhDimBuffer =
            awhDimParamSerialized(coordinateProvider, coordIndex, origin, end, period, diffusion);
    auto                    awhDimArrayRef = gmx::arrayRefFromArray(&awhDimBuffer, 1);
    const AwhTestParameters params(getAwhTestParameters(AwhHistogramGrowthType::ExponentialLinear,
                                                        AwhPotentialType::Umbrella,
                                                        awhDimArrayRef,
                                                        false,
                                                        0.4,
                                                        true,
                                                        1.0,
                                                        c_numLambdaStates));

    const double mdTimeStep = 0.1;

    Bias bias(-1,
              params.awhParams,
              params.awhParams.awhBiasParams()[0],
              params.dimParams,
              params.beta,
              mdTimeStep,
              nullptr,
              "",
              Bias::ThisRankWillDoIO::No);

    const int64_t exitStepRef = 320;

    bool inInitialStage = bias.state().inInitialStage();

    /* Some energies to use as base values (to which some noise is added later on). */
    std::vector<double> neighborLambdaEnergies(c_numLambdaStates);
    std::vector<double> neighborLambdaDhdl(c_numLambdaStates);
    const double        magnitude = 12.0;
    for (int i = 0; i < c_numLambdaStates; i++)
    {
        neighborLambdaEnergies[i] = magnitude * std::sin(i * 0.1);
        neighborLambdaDhdl[i]     = magnitude * std::cos(i * 0.1);
    }

    int64_t step;
    /* Normally this loop exits at exitStepRef, but we extend with failure */
    for (step = 0; step <= 2 * exitStepRef; step++)
    {
        int      umbrellaGridpointIndex = bias.state().coordState().umbrellaGridpoint();
        awh_dvec coordValue = { bias.getGridCoordValue(umbrellaGridpointIndex)[0], 0, 0, 0 };

        double potential     = 0;
        double potentialJump = 0;
        bias.calcForceAndUpdateBias(coordValue,
                                    neighborLambdaEnergies,
                                    neighborLambdaDhdl,
                                    &potential,
                                    &potentialJump,
                                    step,
                                    step,
                                    params.awhParams.seed(),
                                    nullptr);

        inInitialStage = bias.state().inInitialStage();
        if (!inInitialStage)
        {
            break;
        }
    }

    EXPECT_EQ(false, inInitialStage);
    if (!inInitialStage)
    {
        EXPECT_EQ(exitStepRef, step);
    }
}

// Test that we catch too large negative foreign energy differencs
TEST(BiasFepLambdaStateTest, DetectsLargeNegativeForeignEnergy)
{
    constexpr AwhCoordinateProviderType coordinateProvider = AwhCoordinateProviderType::FreeEnergyLambda;
    constexpr int                       coordIndex         = 0;
    constexpr double                    origin             = 0;
    constexpr double                    end                = c_numLambdaStates - 1;
    constexpr double                    period             = 0;
    constexpr double                    diffusion          = 1e-4 / (0.12927243028700 * 2);
    auto                                awhDimBuffer =
            awhDimParamSerialized(coordinateProvider, coordIndex, origin, end, period, diffusion);
    auto                    awhDimArrayRef = gmx::arrayRefFromArray(&awhDimBuffer, 1);
    const AwhTestParameters params(getAwhTestParameters(AwhHistogramGrowthType::ExponentialLinear,
                                                        AwhPotentialType::Umbrella,
                                                        awhDimArrayRef,
                                                        false,
                                                        0.4,
                                                        true,
                                                        1.0,
                                                        c_numLambdaStates));

    const double mdTimeStep = 0.1;

    Bias bias(-1,
              params.awhParams,
              params.awhParams.awhBiasParams()[0],
              params.dimParams,
              params.beta,
              mdTimeStep,
              nullptr,
              "",
              Bias::ThisRankWillDoIO::No);

    /* Some energies to use as base values (to which some noise is added later on). */
    std::vector<double> neighborLambdaEnergies(c_numLambdaStates, 0);
    std::vector<double> neighborLambdaDhdl(c_numLambdaStates, 0);

    // Set last foreign energy to a too large negative value compared to zero
    neighborLambdaEnergies[c_numLambdaStates - 1] =
            -0.51 * gmx::detail::c_largePositiveExponent / params.beta;

    int      umbrellaGridpointIndex = bias.state().coordState().umbrellaGridpoint();
    awh_dvec coordValue = { bias.getGridCoordValue(umbrellaGridpointIndex)[0], 0, 0, 0 };

    double potential     = 0;
    double potentialJump = 0;

    EXPECT_THROW_GMX(bias.calcForceAndUpdateBias(coordValue,
                                                 neighborLambdaEnergies,
                                                 neighborLambdaDhdl,
                                                 &potential,
                                                 &potentialJump,
                                                 0,
                                                 0,
                                                 params.awhParams.seed(),
                                                 nullptr),
                     SimulationInstabilityError);
}

} // namespace test
} // namespace gmx
