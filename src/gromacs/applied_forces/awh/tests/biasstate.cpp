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

#include "gromacs/applied_forces/awh/biasstate.h"

#include <cmath>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/applied_forces/awh/biasgrid.h"
#include "gromacs/applied_forces/awh/biasparams.h"
#include "gromacs/applied_forces/awh/correlationgrid.h"
#include "gromacs/applied_forces/awh/dimparams.h"
#include "gromacs/applied_forces/awh/pointstate.h"
#include "gromacs/applied_forces/awh/tests/awh_setup.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

/*! \brief Test fixture for testing Bias updates
 */
class BiasStateTest : public ::testing::TestWithParam<const char*>
{
private:
    std::unique_ptr<AwhTestParameters> params_;

public:
    std::unique_ptr<BiasState> biasState_; //!< The bias state

    BiasStateTest()
    {
        std::vector<std::vector<char>> awhDimParameters;
        AwhCoordinateProviderType      coordinateProvider = AwhCoordinateProviderType::Pull;
        double                         diffusion          = 0.1;
        {
            int    coordIndex = 0;
            double origin     = 0.5;
            double end        = 1.5;
            double period     = 0;
            awhDimParameters.emplace_back(awhDimParamSerialized(
                    coordinateProvider, coordIndex, origin, end, period, diffusion));
        }
        {
            int    coordIndex = 1;
            double origin     = 0.8;
            double end        = 1.3;
            double period     = 0;
            awhDimParameters.emplace_back(awhDimParamSerialized(
                    coordinateProvider, coordIndex, origin, end, period, diffusion));
        }
        params_                          = std::make_unique<AwhTestParameters>(getAwhTestParameters(
                AwhHistogramGrowthType::Linear, AwhPotentialType::Convolved, awhDimParameters, true, 1.0, false, 0.5, 0));
        const AwhParams&       awhParams = params_->awhParams;
        const AwhBiasParams&   awhBiasParams = awhParams.awhBiasParams(0);
        std::vector<DimParams> dimParams;
        dimParams.push_back(DimParams::pullDimParams(1.0, 15.0, params_->beta));
        dimParams.push_back(DimParams::pullDimParams(1.0, 15.0, params_->beta));
        BiasGrid   grid(dimParams, awhBiasParams.dimParams());
        BiasParams biasParams(
                awhParams, awhBiasParams, dimParams, 1.0, 1.0, BiasParams::DisableUpdateSkips::no, 1, grid.axis(), 0);
        biasState_ = std::make_unique<BiasState>(awhBiasParams, 1.0, dimParams, grid, nullptr);

        /* We let the correlation init function set its parameters
         * to something useful for now.
         */
        double blockLength = 0;
        double mdTimeStep  = 0.002;
        /* Construct the force correlation object. */
        CorrelationGrid forceCorrelationGrid = CorrelationGrid(biasState_->points().size(),
                                                               dimParams.size(),
                                                               blockLength,
                                                               CorrelationGrid::BlockLengthMeasure::Time,
                                                               awhParams.nstSampleCoord() * mdTimeStep);

        // Here we initialize the grid point state using the input file
        std::string filename = gmx::test::TestFileManager::getInputFilePath(GetParam()).string();
        biasState_->initGridPointState(awhBiasParams,
                                       dimParams,
                                       grid,
                                       biasParams,
                                       forceCorrelationGrid,
                                       filename,
                                       params_->awhParams.numBias());
    }
};

/*! \brief Test user input data reading
 */
class UserInputTest : public ::testing::TestWithParam<const char*>
{
private:
    std::unique_ptr<AwhTestParameters> params_;

public:
    std::unique_ptr<BiasGrid>                                      grid_;
    std::vector<int>                                               gridIndexToDataIndex_;
    gmx::MultiDimArray<std::vector<double>, gmx::dynamicExtents2D> data_;
    int                                                            numColumns_;
    int                                                            numRows_;
    std::filesystem::path                                          filename_;

    UserInputTest()
    {
        std::vector<std::vector<char>> awhDimParameters;
        AwhCoordinateProviderType      coordinateProvider = AwhCoordinateProviderType::Pull;
        double                         diffusion          = 0.1;
        {
            int    coordIndex = 0;
            double origin     = 0.5;
            double end        = 1.5;
            double period     = 0;
            awhDimParameters.emplace_back(awhDimParamSerialized(
                    coordinateProvider, coordIndex, origin, end, period, diffusion));
        }
        {
            int    coordIndex = 1;
            double origin     = 0.8;
            double end        = 1.3;
            double period     = 0;
            awhDimParameters.emplace_back(awhDimParamSerialized(
                    coordinateProvider, coordIndex, origin, end, period, diffusion));
        }
        {
            int    coordIndex = 2;
            double origin     = 0.5;
            double end        = 1.0;
            double period     = 0;
            awhDimParameters.emplace_back(awhDimParamSerialized(
                    coordinateProvider, coordIndex, origin, end, period, diffusion));
        }
        params_                          = std::make_unique<AwhTestParameters>(getAwhTestParameters(
                AwhHistogramGrowthType::Linear, AwhPotentialType::Convolved, awhDimParameters, true, 1.0, false, 0.5, 0));
        const AwhParams&       awhParams = params_->awhParams;
        const AwhBiasParams&   awhBiasParams = awhParams.awhBiasParams(0);
        std::vector<DimParams> dimParams;
        dimParams.push_back(DimParams::pullDimParams(1.0, 15.0, params_->beta));
        dimParams.push_back(DimParams::pullDimParams(1.0, 15.0, params_->beta));
        dimParams.push_back(DimParams::pullDimParams(1.0, 15.0, params_->beta));
        grid_                 = std::make_unique<BiasGrid>(dimParams, awhBiasParams.dimParams());
        gridIndexToDataIndex_ = std::vector<int>(grid_->numPoints());

        // Here we read the input file
        filename_ = gmx::test::TestFileManager::getInputFilePath(GetParam());

        data_       = readXvgData(filename_);
        numColumns_ = data_.extent(0);
        numRows_    = data_.extent(1);
    }
};

TEST_P(BiasStateTest, InitializesFromFile)
{
    gmx::ArrayRef<const PointState> points = biasState_->points();

    /* Compute the mean square deviation from the expected values in the file.
     * The PMF values are spaced by 0.5 per points and logPmfsum has opposite sign.
     * The target is (index + 1)/120.
     */
    double msdPmf = 0;
    for (Index i = 0; i < points.ssize(); i++)
    {
        msdPmf += gmx::square(points[i].logPmfSum() - points[0].logPmfSum() + 0.5 * i) / points.size();
        EXPECT_DOUBLE_EQ(points[i].target(), (i + 1) / 120.0);
    }

    EXPECT_NEAR(0.0, msdPmf, 1e-31);
}

TEST_P(UserInputTest, ParsesUser3DInput)
{
    const BiasGrid& grid = *grid_;
    std::string     correctFormatMessage;
    /* Get a data point for each AWH grid point so that they all get data. */
    EXPECT_NO_THROW(mapGridToDataGrid(
            &gridIndexToDataIndex_, data_, numRows_, filename_.string(), grid, correctFormatMessage));
    EXPECT_EQ(numRows_, 30);
    EXPECT_EQ(numColumns_, 8);
}


// Test that Bias initialization open and reads the correct initialization
// files and the correct PMF and target distribution is set.
INSTANTIATE_TEST_SUITE_P(WithParameters,
                         BiasStateTest,
                         ::testing::Values("pmf_target_format0.xvg", "pmf_target_format1.xvg"));

INSTANTIATE_TEST_SUITE_P(WithParameters, UserInputTest, ::testing::Values("pmf_target_format2.xvg"));

} // namespace test
} // namespace gmx
