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

#include "gromacs/awh/biasstate.h"

#include <cmath>

#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/awh/grid.h"
#include "gromacs/awh/pointstate.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

/*! \internal \brief
 * Struct that gathers all input for setting up and using a Bias
 */
struct AwhTestParameters
{
    double        beta;            //!< 1/(kB*T)

    AwhDimParams  awhDimParams[2]; //!< Dimension parameters pointed to by \p awhBiasParams
    AwhBiasParams awhBiasParams;   //!< Bias parameters pointed to by \[ awhParams
    AwhParams     awhParams;       //!< AWH parameters, this is the struct to actually use
};

//! Helper function to set up the C-style AWH parameters for the test
static AwhTestParameters getAwhTestParameters()
{
    AwhTestParameters params;

    params.beta = 1.0;

    AwhParams     &awhParams      = params.awhParams;
    snew(params.awhParams.awhBiasParams, 1);
    AwhBiasParams &awhBiasParams  = params.awhParams.awhBiasParams[0];
    snew(awhBiasParams.dimParams, 2);

    AwhDimParams  &awhDimParams0  = awhBiasParams.dimParams[0];

    awhDimParams0.period          = 0;
    awhDimParams0.diffusion       = 0.1;
    awhDimParams0.origin          = 0.5;
    awhDimParams0.end             = 1.5;
    awhDimParams0.coordValueInit  = awhDimParams0.origin;
    awhDimParams0.coverDiameter   = 0;

    AwhDimParams &awhDimParams1   = awhBiasParams.dimParams[1];

    awhDimParams1.period          = 0;
    awhDimParams1.diffusion       = 0.1;
    awhDimParams1.origin          = 0.8;
    awhDimParams1.end             = 1.3;
    awhDimParams1.coordValueInit  = awhDimParams1.origin;
    awhDimParams1.coverDiameter   = 0;

    awhBiasParams.ndim                 = 2;
    awhBiasParams.eTarget              = eawhtargetCONSTANT;
    awhBiasParams.targetBetaScaling    = 0;
    awhBiasParams.targetCutoff         = 0;
    awhBiasParams.eGrowth              = eawhgrowthLINEAR;
    awhBiasParams.bUserData            = TRUE;
    awhBiasParams.errorInitial         = 0.5;
    awhBiasParams.shareGroup           = 0;
    awhBiasParams.equilibrateHistogram = FALSE;

    awhParams.numBias                    = 1;
    awhParams.seed                       = 93471803;
    awhParams.nstOut                     = 0;
    awhParams.nstSampleCoord             = 1;
    awhParams.numSamplesUpdateFreeEnergy = 10;
    awhParams.ePotential                 = eawhpotentialCONVOLVED;
    awhParams.shareBiasMultisim          = FALSE;

    return params;
}

/*! \brief Test fixture for testing Bias updates
 */
class BiasStateTest : public ::testing::TestWithParam<const char *>
{
    public:
        std::unique_ptr<BiasState> biasState_;  //!< The bias state

        //! Constructor
        BiasStateTest()
        {
            AwhTestParameters       params        = getAwhTestParameters();
            const AwhParams        &awhParams     = params.awhParams;
            const AwhBiasParams    &awhBiasParams = awhParams.awhBiasParams[0];
            std::vector<DimParams>  dimParams;
            dimParams.push_back(DimParams(1.0, 15.0, params.beta));
            dimParams.push_back(DimParams(1.0, 15.0, params.beta));
            Grid                    grid(dimParams, awhBiasParams.dimParams);
            BiasParams              biasParams(awhParams, awhBiasParams, dimParams, 1.0, 1.0, BiasParams::DisableUpdateSkips::no, 1, grid.axis(), 0);
            biasState_ = std::unique_ptr<BiasState>(new BiasState(awhBiasParams, 1.0, dimParams, grid));

            // Here we initialize the grid point state using the input file
            std::string             filename = gmx::test::TestFileManager::getInputFilePath(GetParam());
            biasState_->initGridPointState(awhBiasParams, dimParams, grid, biasParams, filename, params.awhParams.numBias);

            sfree(params.awhParams.awhBiasParams[0].dimParams);
            sfree(params.awhParams.awhBiasParams);
        }
};

TEST_P(BiasStateTest, InitializesFromFile)
{
    gmx::ArrayRef<const PointState> points = biasState_->points();

    /* Compute the mean square deviation from the expected values in the file.
     * The PMF values are spaced by 0.5 per points and logPmfsum has opposite sign.
     * The target is (index + 1)/120.
     */
    double msdPmf    = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
        msdPmf += gmx::square(points[i].logPmfSum() - points[0].logPmfSum() + 0.5*i) / points.size();
        EXPECT_DOUBLE_EQ(points[i].target(), (i + 1)/120.0);
    }

    EXPECT_NEAR(0.0, msdPmf, 1e-31);
}

// Test that Bias initialization open and reads the correct initialization
// files and the the correct PMF and target distribution is set.
INSTANTIATE_TEST_CASE_P(WithParameters, BiasStateTest,
                            ::testing::Values("pmf_target_format0.xvg",
                                              "pmf_target_format1.xvg"));

} // namespace
} // namespace
