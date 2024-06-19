/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
#include <cstddef>
#include <cstdint>

#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "gromacs/applied_forces/awh/correlationgrid.h"
#include "gromacs/applied_forces/awh/correlationhistory.h"
#include "gromacs/applied_forces/awh/correlationtensor.h"
#include "gromacs/mdtypes/awh_correlation_history.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

//! Database of 80 test coordinates (fraction of point range) that represent a trajectory.
constexpr double g_coords[] = { 0.0, 0.2, 0.2, 0.3, 0.1, 0.3, 0.5, 0.6, 0.8, 0.6, 0.7, 0.7,
                                0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.8, 0.7, 0.7, 0.7, 0.6, 0.6,
                                0.6, 0.5, 0.4, 0.4, 0.2, 0.1, 0.0, 0.0, 0.1, 0.1, 0.1, 0.3,
                                0.3, 0.3, 0.3, 0.4, 0.5, 0.4, 0.2, 0.2, 0.1, 0.1, 0.9, 0.1,
                                0.3, 0.7, 0.6, 0.5, 0.5, 0.4, 0.1, 0.0, 0.1, 0.0, 0.2, 0.4,
                                0.4, 0.3, 0.3, 0.3, 0.5, 0.6, 0.8, 0.9, 0.9, 0.7, 0.5, 0.4,
                                0.3, 0.3, 0.2, 0.1, 0.0, 0.1, 0.1, 0.2 };
//! The forces, per sample, below just provide a template to base the forces on.
constexpr double g_forces[] = {
    0.10, 0.12, 0.15, 0.18, 0.10, 0.24, 0.26, 0.80, 0.84, 0.85, 0.83, 0.82, 0.65, 0.60, 0.73, 0.65,
    0.45, 0.50, 0.65, 0.55, 0.56, 0.55, 0.50, 0.50, 0.53, 0.50, 0.40, 0.38, 0.18, 0.18, 0.10, 0.09,
    0.13, 0.13, 0.13, 0.15, 0.20, 0.35, 0.50, 0.54, 0.30, 0.25, 0.35, 0.38, 0.23, 0.40, 0.50, 0.15,
    0.25, 0.40, 0.50, 0.55, 0.53, 0.50, 0.32, 0.14, 0.20, 0.12, 0.15, 0.20, 0.33, 0.33, 0.35, 0.38,
    0.62, 0.56, 0.64, 0.58, 0.66, 0.58, 0.58, 0.40, 0.32, 0.38, 0.50, 0.50, 0.30, 0.20, 0.20, 0.20
};

constexpr double g_sampleWeightCentral  = 0.8;
constexpr double g_sampleWeightNeighbor = 0.1;
constexpr double g_sampleTimeStep       = 0.1;
constexpr int    g_numPointsPerDim      = 10;

/*! \brief Test fixture for testing Friction metric, including storing to, and restoring from, history.
 */
class FrictionMetricTest : public ::testing::TestWithParam<int>
{
public:
    //! Coordinates representing a trajectory in time
    std::vector<double> coordinates_;
    //! Forces that will be used for the trajectory (modified in dimensions > 1)
    std::vector<double> forces_;
    //! The correlation grid
    std::unique_ptr<CorrelationGrid> correlationGrid_;
    //! Number of dimensions
    int64_t numDim_;
    //! Number of grid points
    int64_t numPoints_;
    //! Time between samples
    double sampleTimeStep_ = g_sampleTimeStep;

    FrictionMetricTest() :
        coordinates_(std::begin(g_coords), std::end(g_coords)),
        forces_(std::begin(g_forces), std::end(g_forces))
    {
        /* Set up a basic Correlation Grid. */
        constexpr double                    blockLengthInit = 0;
        CorrelationGrid::BlockLengthMeasure blockLengthMeasure = CorrelationGrid::BlockLengthMeasure::Time;
        numDim_                                                = GetParam();
        GMX_RELEASE_ASSERT(numDim_ < 4, "Too high dimensionality.");
        numPoints_ = std::pow(g_numPointsPerDim, numDim_);

        correlationGrid_ = std::make_unique<CorrelationGrid>(
                numPoints_, numDim_, blockLengthInit, blockLengthMeasure, sampleTimeStep_);
    }
};

TEST_P(FrictionMetricTest, FrictionMetric)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    GMX_RELEASE_ASSERT(coordinates_.size() == forces_.size(),
                       "The number of coordinate steps and force steps do not match.");
    size_t           numSteps        = coordinates_.size();
    CorrelationGrid& correlationGrid = *correlationGrid_;

    for (size_t step = 0; step < numSteps; step++)
    {
        const double coord      = coordinates_[step];
        const int    pointIndex = coord * numPoints_;
        /* Make a set of neighbors just based on the grid point indices. Set their sampling weight arbitrarily. */
        std::vector<int>    neighbors;
        std::vector<double> neighborSampleWeights;
        neighbors.push_back(pointIndex);
        neighborSampleWeights.push_back(g_sampleWeightCentral);
        if (pointIndex > 0)
        {
            neighbors.push_back(pointIndex - 1);
            neighborSampleWeights.push_back(g_sampleWeightNeighbor);
        }
        if (pointIndex < numPoints_ - 1)
        {
            neighbors.push_back(pointIndex + 1);
            neighborSampleWeights.push_back(g_sampleWeightNeighbor);
        }
        const size_t numNeighbors = neighbors.size();
        for (size_t neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++)
        {
            std::vector<double> force;
            for (int64_t dim = 0; dim < numDim_; dim++)
            {
                /* Make the forces a little different in dimensions > 1 */
                force.push_back(g_forces[step] * (dim + 1));
            }
            correlationGrid.addData(neighbors[neighborIndex],
                                    neighborSampleWeights[neighborIndex],
                                    force,
                                    step * sampleTimeStep_);
        }
    }
    std::vector<CorrelationTensor> tensors = correlationGrid.tensors();
    std::vector<double>            tensorCorrelationIntegrals;
    int                            correlationIntegralSize;
    switch (numDim_)
    {
        case 1: correlationIntegralSize = 1; break;
        case 2: correlationIntegralSize = 3; break;
        case 3: correlationIntegralSize = 6; break;
        default: GMX_THROW(gmx::InternalError("Too high dimensionality."));
    }
    for (const gmx::CorrelationTensor& tensor : tensors)
    {
        for (int correlationIntegralIndex = 0; correlationIntegralIndex < correlationIntegralSize;
             correlationIntegralIndex++)
        {
            tensorCorrelationIntegrals.push_back(
                    tensor.getTimeIntegral(correlationIntegralIndex, sampleTimeStep_));
        }
    }
    checker.checkSequence(tensorCorrelationIntegrals.begin(),
                          tensorCorrelationIntegrals.end(),
                          "tensorCorrelationIntegrals");

    /* Store to history and restore from history. */
    CorrelationGridHistory correlationHistory = initCorrelationGridHistoryFromState(correlationGrid);
    updateCorrelationGridHistory(&correlationHistory, correlationGrid);
    correlationGrid.restoreStateFromHistory(correlationHistory);
    tensors     = correlationGrid.tensors();
    int counter = 0;
    /* Check that the tensor integrals are the same before and after restoring from history. */
    for (const gmx::CorrelationTensor& tensor : tensors)
    {
        for (int correlationIntegralIndex = 0; correlationIntegralIndex < correlationIntegralSize;
             correlationIntegralIndex++)
        {
            EXPECT_DOUBLE_EQ(tensor.getTimeIntegral(correlationIntegralIndex, sampleTimeStep_),
                             tensorCorrelationIntegrals[counter++]);
        }
    }
}

/* Test correlation grids of the dimensions listed in Values.
 */
INSTANTIATE_TEST_SUITE_P(WithParameters, FrictionMetricTest, ::testing::Values(1, 2, 3));

} // namespace test
} // namespace gmx
