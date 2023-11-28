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

#include "gromacs/applied_forces/awh/bias.h"

#include <cmath>

#include <memory>
#include <tuple>
#include <vector>

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/applied_forces/awh/correlationgrid.h"
#include "gromacs/applied_forces/awh/pointstate.h"
#include "gromacs/applied_forces/awh/tests/awh_setup.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

//! Database of 21 test coordinates that represent a trajectory */
constexpr double g_coords[] = { 0.62, 0.70, 0.68, 0.80, 0.93, 0.87, 1.16, 1.14, 0.95, 0.89, 0.91,
                                0.86, 0.88, 0.79, 0.75, 0.82, 0.74, 0.70, 0.68, 0.71, 0.73 };

//! Convenience typedef: growth type enum, potential type enum, disable update skips
typedef std::tuple<AwhHistogramGrowthType, AwhPotentialType, BiasParams::DisableUpdateSkips> BiasTestParameters;

/*! \brief Test fixture for testing Bias updates
 */
class BiasTest : public ::testing::TestWithParam<BiasTestParameters>
{
private:
    //! Storage for test parameters.
    std::unique_ptr<AwhTestParameters> params_;

public:
    //! Random seed for AWH MC sampling
    int64_t seed_;
    //! Coordinates representing a trajectory in time
    std::vector<double> coordinates_;
    //! The awh Bias
    std::unique_ptr<Bias> bias_;

    //! Return the awhParams
    const AwhParams& awhParams() const { return params_->awhParams; }

    BiasTest() : coordinates_(std::begin(g_coords), std::end(g_coords))
    {
        /* We test all combinations of:
         *   eawhgrowth:
         *     eawhgrowthLINEAR:     final, normal update phase
         *     ewahgrowthEXP_LINEAR: initial phase, updated size is constant
         *   eawhpotential (should only affect the force output):
         *     eawhpotentialUMBRELLA:  MC on lambda (umbrella potential location)
         *     eawhpotentialCONVOLVED: MD on a convolved potential landscape
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
         * The idea is to, among other things, have part of the interval
         * not covered by samples.
         */
        auto awhDimBuffer   = awhDimParamSerialized();
        auto awhDimArrayRef = gmx::arrayRefFromArray(&awhDimBuffer, 1);
        params_             = std::make_unique<AwhTestParameters>(getAwhTestParameters(
                eawhgrowth, eawhpotential, awhDimArrayRef, false, 0.4, false, 0.5, 0));

        seed_ = params_->awhParams.seed();

        constexpr double mdTimeStep = 0.1;

        int numSamples = coordinates_.size() - 1; // No sample taken at step 0
        GMX_RELEASE_ASSERT(numSamples % params_->awhParams.numSamplesUpdateFreeEnergy() == 0,
                           "This test is intended to reproduce the situation when the might need "
                           "to write output during a normal AWH run, therefore the number of "
                           "samples should be a multiple of the free-energy update interval (but "
                           "the test should also runs fine without this condition).");

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

TEST_P(BiasTest, ForcesBiasPmfWeightSum)
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

    std::vector<double> force, pot, potJump;

    double  coordMaxValue = 0;
    double  potentialJump = 0;
    int64_t step          = 0;
    for (auto& coord : coordinates_)
    {
        coordMaxValue = std::max(coordMaxValue, std::abs(coord));

        awh_dvec                    coordValue = { coord, 0, 0, 0 };
        double                      potential  = 0;
        gmx::ArrayRef<const double> biasForce  = bias.calcForceAndUpdateBias(
                coordValue, {}, {}, &potential, &potentialJump, step, step, seed_, nullptr);

        force.push_back(biasForce[0]);
        pot.push_back(potential);
        potJump.push_back(potentialJump);

        step++;
    }

    /* When skipping updates, ensure all skipped updates are performed here.
     * This should result in the same bias state as at output in a normal run.
     */
    if (bias.params().skipUpdates())
    {
        bias.doSkippedUpdatesForAllPoints();
    }

    std::vector<double> pointBias, logPmfsum, pointLocalWeightSum, pointWeightSumTot,
            pointWeightSumIteration;
    for (const auto& point : bias.state().points())
    {
        pointBias.push_back(point.bias());
        logPmfsum.push_back(point.logPmfSum());
        pointLocalWeightSum.push_back(point.localWeightSum());
        pointWeightSumTot.push_back(point.weightSumTot());
        pointWeightSumIteration.push_back(point.weightSumIteration());
    }

    /* The umbrella force is computed from the coordinate deviation.
     * In taking this deviation we lose a lot of precision, so we should
     * compare against k*max(coord) instead of the instantaneous force.
     */
    const double kCoordMax = bias.dimParams()[0].pullDimParams().k * coordMaxValue;

    constexpr int lowUlpTol = 10;
    /* Allow higher tolerance for weight sum calculations - base it on the number of steps. */
    const double highUlpTol = 10 * std::sqrt(step);

    checker.checkSequence(props.begin(), props.end(), "Properties");
    checker.setDefaultTolerance(absoluteTolerance(kCoordMax * GMX_DOUBLE_EPS * lowUlpTol));
    checker.checkSequence(force.begin(), force.end(), "Force");
    checker.checkSequence(pot.begin(), pot.end(), "Potential");
    checker.checkSequence(potJump.begin(), potJump.end(), "PotentialJump");
    checker.setDefaultTolerance(relativeToleranceAsUlp(1.0, lowUlpTol));
    checker.checkSequence(pointBias.begin(), pointBias.end(), "PointBias");
    checker.checkSequence(logPmfsum.begin(), logPmfsum.end(), "PointLogPmfsum");
    checker.setDefaultTolerance(relativeToleranceAsUlp(1.0, highUlpTol));
    checker.checkSequence(
            pointLocalWeightSum.begin(), pointLocalWeightSum.end(), "pointLocalWeightSum");
    checker.checkSequence(pointWeightSumTot.begin(), pointWeightSumTot.end(), "pointWeightSumTot");
    checker.checkSequence(pointWeightSumIteration.begin(),
                          pointWeightSumIteration.end(),
                          "pointWeightSumIteration");
}

/* Scan initial/final phase, MC/convolved force and update skip (not) allowed
 * Both the convolving and skipping should not affect the bias and PMF.
 * It would be nice if the test would explicitly check for this.
 * Currently this is tested through identical reference data.
 */
INSTANTIATE_TEST_SUITE_P(WithParameters,
                         BiasTest,
                         ::testing::Combine(::testing::Values(AwhHistogramGrowthType::Linear,
                                                              AwhHistogramGrowthType::ExponentialLinear),
                                            ::testing::Values(AwhPotentialType::Umbrella,
                                                              AwhPotentialType::Convolved),
                                            ::testing::Values(BiasParams::DisableUpdateSkips::yes,
                                                              BiasParams::DisableUpdateSkips::no)));

// Test that we detect coverings and exit the initial stage at the correct step
TEST(BiasTest, DetectsCovering)
{
    const std::vector<char> serializedAwhParametersPerDim = awhDimParamSerialized();
    auto awhDimArrayRef            = gmx::arrayRefFromArray(&serializedAwhParametersPerDim, 1);
    const AwhTestParameters params = getAwhTestParameters(AwhHistogramGrowthType::ExponentialLinear,
                                                          AwhPotentialType::Convolved,
                                                          awhDimArrayRef,
                                                          false,
                                                          0.4,
                                                          false,
                                                          0.5,
                                                          0);
    const AwhDimParams&     awhDimParams = params.awhParams.awhBiasParams(0).dimParams(0);

    constexpr double mdTimeStep = 0.1;

    Bias bias(-1,
              params.awhParams,
              params.awhParams.awhBiasParams(0),
              params.dimParams,
              params.beta,
              mdTimeStep,
              nullptr,
              "",
              Bias::ThisRankWillDoIO::No);

    /* We use a trajectory of the sum of two sines to cover the reaction
     * coordinate range in a semi-realistic way. The period is 4*pi=12.57.
     * We get out of the initial stage after 4 coverings at step 300.
     */
    constexpr int64_t exitStepRef = 300;
    const double      midPoint    = 0.5 * (awhDimParams.end() + awhDimParams.origin());
    const double      halfWidth   = 0.5 * (awhDimParams.end() - awhDimParams.origin());

    bool inInitialStage = bias.state().inInitialStage();
    /* Normally this loop exits at exitStepRef, but we extend with failure */
    int64_t step;
    for (step = 0; step <= 2 * exitStepRef; step++)
    {
        double t     = step * mdTimeStep;
        double coord = midPoint + halfWidth * (0.5 * std::sin(t) + 0.55 * std::sin(1.5 * t));

        awh_dvec coordValue    = { coord, 0, 0, 0 };
        double   potential     = 0;
        double   potentialJump = 0;
        bias.calcForceAndUpdateBias(
                coordValue, {}, {}, &potential, &potentialJump, step, step, params.awhParams.seed(), nullptr);

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
