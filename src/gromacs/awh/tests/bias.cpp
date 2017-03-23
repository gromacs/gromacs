/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include <memory>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/awh/bias.h"
#include "gromacs/awh/grid.h"
#include "gromacs/awh/pointstate.h"
#include "gromacs/mdtypes/awh-params.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

//! Database of 21 test coordinates that represent a trajectory */
const double g_coords[] = {
    0.62,
    0.70,
    0.68,
    0.80,
    0.93,
    0.87,
    1.16,
    1.14,
    0.95,
    0.89,
    0.91,
    0.86,
    0.88,
    0.79,
    0.75,
    0.82,
    0.74,
    0.70,
    0.68,
    0.71,
    0.73
};

//! Convenience typedef
typedef std::tuple<bool, bool, bool> BiasTestParameters;

/*! \brief Test fixture for testing Bias updates
 */
class BiasTest : public ::testing::TestWithParam<BiasTestParameters>
{
    public:
        //! Random seed for AWH MC sampling
        static const int c_seed = 93471803;

        //! Updated water atom positions to constrain (DIM reals per atom)
        std::vector<real>     coordinates_;
        //! The awh Bias
        std::unique_ptr<Bias> bias_;

        //! Constructor
        BiasTest() :
            coordinates_(std::begin(g_coords), std::end(g_coords))
        {
            // We test all combinations of:
            //   useConstantUpdateSize:
            //     false: final, normal update phase
            //     true:  intial phase, updated size is constant
            //   convolveForce (should only affect the force output):
            //     false: MC on lambda (umbrella potential location)
            //     true:  MD on a convolved potential landscape
            //   allowUpdateSkips (should not affect the results):
            //     false: update the point state for every sample
            //     true:  update the point state at an interval > 1 sample
            //
            // Note: It would be nice to explicitly check that convolveForce
            //       and allowUpdateSkips do not affect the point state.
            //       But the reference data will also ensure this.
            //
            bool useConstantUpdateSize, convolveForce, allowUpdateSkips;
            std::tie(useConstantUpdateSize, convolveForce, allowUpdateSkips) = GetParam();

            // Set up a basic AWH setup with a single, 1D bias with parameters
            // such that we can measure the effects of different parameters.
            // The idea is to, among other things, have part of the interval
            // not covered by samples.
            awh_dim_params_t awhDimParams;

            awhDimParams.period           = 0;
            awhDimParams.diffusion        = 0.1;
            awhDimParams.origin           = 0.5;
            awhDimParams.end              = 1.5;
            awhDimParams.ninterval        = 1;
            awhDimParams.interval_overlap = 1;
            awhDimParams.coverDiameter    = 0;

            awh_bias_params_t awhBiasParams;

            awhBiasParams.ndim                 = 1;
            awhBiasParams.dim_params           = &awhDimParams;
            awhBiasParams.eTarget              = eawhtargetCONSTANT;
            awhBiasParams.targetBetaScaling    = 0;
            awhBiasParams.targetCutoff         = 0;
            awhBiasParams.eGrowth              = useConstantUpdateSize ? eawhgrowthEXP_LINEAR : eawhgrowthLINEAR;
            awhBiasParams.bUser_data           = 0;
            awhBiasParams.error_initial        = 0.5;
            awhBiasParams.bShare               = FALSE;
            awhBiasParams.equilibrateHistogram = TRUE;

            double                 convFactor  = 1;
            double                 k           = 1000;
            double                 beta        = 0.4;
            double                 mdTimeStep  = 0.1;

            std::vector<DimParams> dimParams   = { DimParams(convFactor, k, beta) };

            awh_params_t           awhParams;

            awhParams.nbias                       = 1;
            awhParams.awh_bias_params             = &awhBiasParams;
            awhParams.seed                        = c_seed;
            awhParams.nstout                      = 0;
            awhParams.nstsample_coord             = 1;
            awhParams.nsamples_update_free_energy = 10;
            awhParams.ePotential                  = (convolveForce ? eawhpotentialCONVOLVED : eawhpotentialUMBRELLA);

            GMX_RELEASE_ASSERT((coordinates_.size() - 1) % awhParams.nsamples_update_free_energy == 0, "This test is intended to reproduce the situation when the might need to write output during a normal AWH run, therefore the number of samples should be a multiple of the free-energy update interal (but the test should also runs fine without this condition).");

            bias_ = std::unique_ptr<Bias>(new Bias(nullptr, nullptr, -1, awhParams, awhBiasParams, dimParams, beta, mdTimeStep, allowUpdateSkips ? BiasParams::DisableUpdateSkips::no : BiasParams::DisableUpdateSkips::yes));
        }
};

TEST_P(BiasTest, ForcesBiasPmf)
{
    gmx::test::TestReferenceData            data;
    gmx::test::TestReferenceChecker         checker(data.rootChecker());

    std::vector<double>                     force;

    Bias        &bias          = *bias_.get();

    double       potentialJump = 0;
    gmx_int64_t  step          = 0;
    for (auto &coord : coordinates_)
    {
        bias.setCoordValue(0, coord);

        awh_dvec biasForce;
        double   potential = 0;
        bias.doStep(biasForce, &potential, &potentialJump,
                    nullptr, step, step, c_seed, nullptr);

        force.push_back(biasForce[0]);

        step++;
    }

    // When skipping updates, ensure all skipped updates are performed here.
    // This should result in the same bias state as at output in a normal run.
    if (bias.params().skipUpdates())
    {
        bias.doSkippedUpdatesForAllPoints();
    }

    std::vector<double> pointBias, logPmfsum;
    for (auto &point : bias.pointState())
    {
        pointBias.push_back(point.bias());
        logPmfsum.push_back(point.logPmfsum());
    }
    checker.checkSequence(force.begin(),     force.end(),     "Force");
    checker.checkSequence(pointBias.begin(), pointBias.end(), "PointBias");
    checker.checkSequence(logPmfsum.begin(), logPmfsum.end(), "PointLogPmfsum");
}

// Scan initial/final phase, MC/convolved force and update skip (not) allowed
INSTANTIATE_TEST_CASE_P(WithParameters, BiasTest,
                            ::testing::Combine(::testing::Bool(),
                                                   ::testing::Bool(),
                                                   ::testing::Bool()));
} // namespace
} // namespace
