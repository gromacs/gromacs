/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include <cmath>

#include <memory>
#include <random>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/applied_forces/awh/bias.h"
#include "gromacs/applied_forces/awh/correlationgrid.h"
#include "gromacs/applied_forces/awh/pointstate.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

//! The number of lambda states to use in the tests.
const int numLambdaStates = 16;

/*! \internal \brief
 * Struct that gathers all input for setting up and using a Bias
 */
struct AwhFepLambdaStateTestParameters
{
    AwhFepLambdaStateTestParameters() = default;
    //! Move constructor
    AwhFepLambdaStateTestParameters(AwhFepLambdaStateTestParameters&& o) noexcept :
        beta(o.beta),
        awhDimParams(o.awhDimParams),
        awhBiasParams(o.awhBiasParams),
        awhParams(o.awhParams),
        dimParams(std::move(o.dimParams))
    {
        awhBiasParams.dimParams = &awhDimParams;
        awhParams.awhBiasParams = &awhBiasParams;
    }
    double beta; //!< 1/(kB*T)

    AwhDimParams  awhDimParams;  //!< Dimension parameters pointed to by \p awhBiasParams
    AwhBiasParams awhBiasParams; //!< Bias parameters pointed to by \[ awhParams
    AwhParams     awhParams;     //!< AWH parameters, this is the struct to actually use

    std::vector<DimParams> dimParams; //!< Dimension parameters for setting up Bias
};

/*! \internal \brief Helper function to fill an array with random values (between lowerBound and
 * upperBound) from randomEngine.
 */
static void randomArrayFill(ArrayRef<double>           array,
                            std::default_random_engine randomEngine,
                            double                     lowerBound,
                            double                     upperBound)
{
    std::uniform_real_distribution<double> unif(lowerBound, upperBound);
    for (size_t i = 0; i < array.size(); i++)
    {
        array[i] = unif(randomEngine);
    }
}

//! Helper function to set up the C-style AWH parameters for the test
static AwhFepLambdaStateTestParameters getAwhTestParameters(int eawhgrowth, int eawhpotential)
{
    AwhFepLambdaStateTestParameters params;

    params.beta = 0.4;

    AwhDimParams& awhDimParams = params.awhDimParams;

    awhDimParams.period         = 0;
    awhDimParams.diffusion      = 1e-4;
    awhDimParams.origin         = 0;
    awhDimParams.end            = numLambdaStates - 1;
    awhDimParams.coordValueInit = awhDimParams.origin;
    awhDimParams.coverDiameter  = 0;
    awhDimParams.eCoordProvider = eawhcoordproviderFREE_ENERGY_LAMBDA;

    AwhBiasParams& awhBiasParams = params.awhBiasParams;

    awhBiasParams.ndim                 = 1;
    awhBiasParams.dimParams            = &awhDimParams;
    awhBiasParams.eTarget              = eawhtargetCONSTANT;
    awhBiasParams.targetBetaScaling    = 0;
    awhBiasParams.targetCutoff         = 0;
    awhBiasParams.eGrowth              = eawhgrowth;
    awhBiasParams.bUserData            = FALSE;
    awhBiasParams.errorInitial         = 1.0 / params.beta;
    awhBiasParams.shareGroup           = 0;
    awhBiasParams.equilibrateHistogram = FALSE;

    double  k    = 1000;
    int64_t seed = 93471803;

    params.dimParams.emplace_back(k, params.beta, numLambdaStates);

    AwhParams& awhParams = params.awhParams;

    awhParams.numBias                    = 1;
    awhParams.awhBiasParams              = &awhBiasParams;
    awhParams.seed                       = seed;
    awhParams.nstOut                     = 0;
    awhParams.nstSampleCoord             = 1;
    awhParams.numSamplesUpdateFreeEnergy = 10;
    awhParams.ePotential                 = eawhpotential;
    awhParams.shareBiasMultisim          = FALSE;

    return params;
}

//! Convenience typedef: growth type enum, potential type enum, disable update skips
typedef std::tuple<int, int, BiasParams::DisableUpdateSkips> BiasTestParameters;

/*! \brief Test fixture for testing Bias updates
 */
class BiasFepLambdaStateTest : public ::testing::TestWithParam<BiasTestParameters>
{
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
         *     ewahgrowthEXP_LINEAR: intial phase, updated size is constant
         *   eawhpotential (test both, but for the FEP lambda state dimension MC will in practice be used,
         *                  except that eawhpotentialCONVOLVED also gives a potential output):
         *     eawhpotentialUMBRELLA:  MC on lambda state
         *     eawhpotentialCONVOLVED: MD on a convolved potential landscape (falling back to MC on lambda state)
         *   disableUpdateSkips (should not affect the results):
         *     BiasParams::DisableUpdateSkips::yes: update the point state for every sample
         *     BiasParams::DisableUpdateSkips::no:  update the point state at an interval > 1 sample
         *
         * Note: It would be nice to explicitly check that eawhpotential
         *       and disableUpdateSkips do not affect the point state.
         *       But the reference data will also ensure this.
         */
        int                            eawhgrowth;
        int                            eawhpotential;
        BiasParams::DisableUpdateSkips disableUpdateSkips;
        std::tie(eawhgrowth, eawhpotential, disableUpdateSkips) = GetParam();

        /* Set up a basic AWH setup with a single, 1D bias with parameters
         * such that we can measure the effects of different parameters.
         */
        const AwhFepLambdaStateTestParameters params = getAwhTestParameters(eawhgrowth, eawhpotential);

        seed_ = params.awhParams.seed;

        double mdTimeStep = 0.1;

        bias_ = std::make_unique<Bias>(-1, params.awhParams, params.awhBiasParams, params.dimParams,
                                       params.beta, mdTimeStep, 1, "", Bias::ThisRankWillDoIO::No,
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

    double                     potentialJump        = 0;
    double                     mdTimeStep           = 0.1;
    double                     energyNoiseMagnitude = 1.0;
    double                     dhdlNoiseMagnitude   = 1.5;
    int                        nSteps               = 501;
    std::default_random_engine randomEngine;
    randomEngine.seed(1234);

    /* Some energies to use as base values (to which some noise is added later on). */
    std::vector<double> lambdaEnergyBase(numLambdaStates);
    std::vector<double> lambdaDhdlBase(numLambdaStates);
    const double        magnitude = 12.0;
    for (int i = 0; i < numLambdaStates; i++)
    {
        lambdaEnergyBase[i] = magnitude * std::sin(i * 0.1);
        lambdaDhdlBase[i]   = magnitude * std::cos(i * 0.1);
    }

    for (int step = 0; step < nSteps; step++)
    {
        /* Create some noise and add it to the base values */
        std::vector<double> neighborLambdaEnergyNoise(numLambdaStates);
        std::vector<double> neighborLambdaDhdlNoise(numLambdaStates);
        randomArrayFill(neighborLambdaEnergyNoise, randomEngine, -energyNoiseMagnitude, energyNoiseMagnitude);
        randomArrayFill(neighborLambdaDhdlNoise, randomEngine, -dhdlNoiseMagnitude, dhdlNoiseMagnitude);
        std::vector<double> neighborLambdaEnergies(numLambdaStates);
        std::vector<double> neighborLambdaDhdl(numLambdaStates);
        for (int i = 0; i < numLambdaStates; i++)
        {
            neighborLambdaEnergies[i] = lambdaEnergyBase[i] + neighborLambdaEnergyNoise[i];
            neighborLambdaDhdl[i]     = lambdaDhdlBase[i] + neighborLambdaDhdlNoise[i];
        }

        int      umbrellaGridpointIndex = bias.state().coordState().umbrellaGridpoint();
        awh_dvec coordValue = { bias.getGridCoordValue(umbrellaGridpointIndex)[0], 0, 0, 0 };
        double   potential  = 0;
        gmx::ArrayRef<const double> biasForce = bias.calcForceAndUpdateBias(
                coordValue, neighborLambdaEnergies, neighborLambdaDhdl, &potential, &potentialJump,
                nullptr, nullptr, step * mdTimeStep, step, seed_, nullptr);

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
    for (auto& point : bias.state().points())
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
INSTANTIATE_TEST_CASE_P(WithParameters,
                        BiasFepLambdaStateTest,
                        ::testing::Combine(::testing::Values(eawhgrowthLINEAR, eawhgrowthEXP_LINEAR),
                                           ::testing::Values(eawhpotentialUMBRELLA, eawhpotentialCONVOLVED),
                                           ::testing::Values(BiasParams::DisableUpdateSkips::yes,
                                                             BiasParams::DisableUpdateSkips::no)));

// Test that we detect coverings and exit the initial stage at the correct step
TEST(BiasFepLambdaStateTest, DetectsCovering)
{
    const AwhFepLambdaStateTestParameters params =
            getAwhTestParameters(eawhgrowthEXP_LINEAR, eawhpotentialCONVOLVED);

    const double mdTimeStep = 0.1;

    Bias bias(-1, params.awhParams, params.awhBiasParams, params.dimParams, params.beta, mdTimeStep,
              1, "", Bias::ThisRankWillDoIO::No);

    const int64_t exitStepRef = 380;

    bool inInitialStage = bias.state().inInitialStage();

    double                     energyNoiseMagnitude = 1.0;
    double                     dhdlNoiseMagnitude   = 1.5;
    std::default_random_engine randomEngine;
    randomEngine.seed(1234);

    /* Some energies to use as base values (to which some noise is added later on). */
    std::vector<double> lambdaEnergyBase(numLambdaStates);
    std::vector<double> lambdaDhdlBase(numLambdaStates);
    const double        magnitude = 12.0;
    for (int i = 0; i < numLambdaStates; i++)
    {
        lambdaEnergyBase[i] = magnitude * std::sin(i * 0.1);
        lambdaDhdlBase[i]   = magnitude * std::cos(i * 0.1);
    }

    int64_t step;
    /* Normally this loop exits at exitStepRef, but we extend with failure */
    for (step = 0; step <= 2 * exitStepRef; step++)
    {
        /* Create some noise and add it to the base values */
        std::vector<double> neighborLambdaEnergyNoise(numLambdaStates);
        std::vector<double> neighborLambdaDhdlNoise(numLambdaStates);
        randomArrayFill(neighborLambdaEnergyNoise, randomEngine, -energyNoiseMagnitude, energyNoiseMagnitude);
        randomArrayFill(neighborLambdaDhdlNoise, randomEngine, -dhdlNoiseMagnitude, dhdlNoiseMagnitude);
        std::vector<double> neighborLambdaEnergies(numLambdaStates);
        std::vector<double> neighborLambdaDhdl(numLambdaStates);
        for (int i = 0; i < numLambdaStates; i++)
        {
            neighborLambdaEnergies[i] = lambdaEnergyBase[i] + neighborLambdaEnergyNoise[i];
            neighborLambdaDhdl[i]     = lambdaDhdlBase[i] + neighborLambdaDhdlNoise[i];
        }

        int      umbrellaGridpointIndex = bias.state().coordState().umbrellaGridpoint();
        awh_dvec coordValue = { bias.getGridCoordValue(umbrellaGridpointIndex)[0], 0, 0, 0 };

        double potential     = 0;
        double potentialJump = 0;
        bias.calcForceAndUpdateBias(coordValue, neighborLambdaEnergies, neighborLambdaDhdl,
                                    &potential, &potentialJump, nullptr, nullptr, step, step,
                                    params.awhParams.seed, nullptr);

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

} // namespace test
} // namespace gmx
