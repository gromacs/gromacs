/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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

#include "gromacs/applied_forces/awh/biassharing.h"

#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "gromacs/mdtypes/commrec.h"
#include "thread_mpi/tmpi.h"

#include "gromacs/applied_forces/awh/tests/awh_setup.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

// This test requires thread-MPI
#if GMX_THREAD_MPI

namespace
{

//! The number of thread-MPI ranks to run this test on
const int c_numRanks = 4;

//! The number of simulations sharing the same bias
const int c_numSharingBiases = 2;

/*! \brief The actual test body, executed by each MPI rank
 *
 * Sets ups sharing and sums over sahring simulations.
 */
void parallelTestFunction(const void gmx_unused* dummy)
{
    int numRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    GMX_RELEASE_ASSERT(numRanks == c_numRanks, "Expect c_numRanks thread-MPI ranks");

    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    const int shareGroup = 1 + (myRank / c_numSharingBiases);

    t_commrec commRecord = { 0 };
    commRecord.nnodes    = 1;

    const std::vector<char> serializedAwhParametersPerDim = awhDimParamSerialized();
    auto              awhDimArrayRef = gmx::arrayRefFromArray(&serializedAwhParametersPerDim, 1);
    AwhTestParameters params = getAwhTestParameters(AwhHistogramGrowthType::ExponentialLinear,
                                                    AwhPotentialType::Convolved,
                                                    awhDimArrayRef,
                                                    false,
                                                    0.4,
                                                    false,
                                                    0.5,
                                                    0,
                                                    shareGroup);

    BiasSharing biasSharing(params.awhParams, commRecord, MPI_COMM_WORLD);

    EXPECT_EQ(biasSharing.numSharingSimulations(0), c_numSharingBiases);
    EXPECT_EQ(biasSharing.sharingSimulationIndex(0), myRank % c_numSharingBiases);

    const std::array<int, c_numRanks> input{ 1, 2, 4, 8 };
    std::array<int, 1>                buffer{ input[myRank] };

    biasSharing.sumOverSharingMasterRanks(buffer, 0);
    int expectedSum = 0;
    for (int i = 0; i < c_numSharingBiases; i++)
    {
        expectedSum += input[(myRank / c_numSharingBiases) * c_numSharingBiases + i];
    }
    EXPECT_EQ(buffer[0], expectedSum);
}

} // namespace

TEST(BiasSharingTest, SharingWorks)
{
    if (tMPI_Init_fn(FALSE, c_numRanks, TMPI_AFFINITY_NONE, parallelTestFunction, static_cast<const void*>(this))
        != TMPI_SUCCESS)
    {
        GMX_THROW(gmx::InternalError("Failed to spawn thread-MPI threads"));
    }
}

#endif

} // namespace test
} // namespace gmx
