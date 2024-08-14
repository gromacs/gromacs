/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief Tests for the halo exchange
 *
 *  The test sets up the rank topology and performs a coordinate halo
 *  exchange (for both CPU and GPU codepaths) for several 1D and 2D
 *  pulse configurations. Each pulse involves a few non-contiguous
 *  indices. The sending rank, atom number and spatial 3D index are
 *  encoded in the x values, to allow correctness checking following
 *  the halo exchange.
 *
 * \todo Add 3D case
 *
 * \author Alan Gray <alang@nvidia.com>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "config.h"

#include <array>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/domdec/atomdistribution.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/device_stream.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#endif
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdtypes/inputrec.h"

#include "testutils/mpitest.h"
#include "testutils/test_hardware_environment.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Get encoded numerical value for sending rank, atom number and spatial 3D index
 *
 * \param [in] sendRank         MPI rank of sender
 * \param [in] atomNumber       Atom number
 * \param [in] spatial3dIndex   Spatial 3D Index
 *
 * \returns                     Encoded value
 */
float encodedValue(const int sendRank, const int atomNumber, const int spatial3dIndex)
{
    return sendRank * 1000 + atomNumber * 100 + spatial3dIndex;
}

/*! \brief Initialize halo array
 *
 * \param [in] x              Atom coordinate data array
 * \param [in] numHomeAtoms   Number of home atoms
 * \param [in] numAtomsTotal  Total number of atoms, including halo
 */
void initHaloData(RVec* x, const int numHomeAtoms, const int numAtomsTotal)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < numAtomsTotal; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            x[i][j] = i < numHomeAtoms ? encodedValue(rank, i, j) : -1;
        }
    }
}

/*! \brief Perform GPU halo exchange, including required setup and data transfers
 *
 * \param [in] dd             Domain decomposition object
 * \param [in] box            Box matrix
 * \param [in] h_x            Atom coordinate data array on host
 * \param [in] numAtomsTotal  Total number of atoms, including halo
 */
void gpuHalo(gmx_domdec_t* dd, matrix box, HostVector<RVec>* h_x, int numAtomsTotal)
{
#if (GMX_GPU_CUDA && GMX_THREAD_MPI)
    // pin memory if possible
    changePinningPolicy(h_x, PinningPolicy::PinnedIfSupported);
    // Set up GPU hardware environment and assign this MPI rank to a device
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int         numDevices = getTestHardwareEnvironment()->getTestDeviceList().size();
    const auto& testDevice = getTestHardwareEnvironment()->getTestDeviceList()[rank % numDevices];
    const auto& deviceContext = testDevice->deviceContext();
    deviceContext.activate();
    DeviceStream deviceStream(deviceContext, DeviceStreamPriority::Normal, false);

    // Set up GPU buffer and copy input data from host
    DeviceBuffer<RVec> d_x;
    int                d_x_size       = -1;
    int                d_x_size_alloc = -1;
    reallocateDeviceBuffer(&d_x, numAtomsTotal, &d_x_size, &d_x_size_alloc, deviceContext);

    copyToDeviceBuffer(&d_x, h_x->data(), 0, numAtomsTotal, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    const int numPulses = std::accumulate(
            dd->comm->cd.begin(), dd->comm->cd.end(), 0, [](const int a, const auto& b) {
                return a + b.numPulses();
            });
    const int numExtraConsumptions = GMX_THREAD_MPI ? dd->ndim : 0;
    // Will be consumed once for each pulse, and, with tMPI, once more for each dim on the first pulse
    GpuEventSynchronizer coordinatesReadyOnDeviceEvent(numPulses + numExtraConsumptions,
                                                       numPulses + numExtraConsumptions);
    coordinatesReadyOnDeviceEvent.markEvent(deviceStream);

    std::array<std::vector<GpuHaloExchange>, DIM> gpuHaloExchange;

    // Create halo exchange objects
    for (int d = 0; d < dd->ndim; d++)
    {
        for (int pulse = 0; pulse < dd->comm->cd[d].numPulses(); pulse++)
        {
            gpuHaloExchange[d].push_back(
                    GpuHaloExchange(dd, d, MPI_COMM_WORLD, deviceContext, pulse, nullptr));
        }
    }

    // Perform GPU halo exchange
    for (int d = 0; d < dd->ndim; d++)
    {
        for (int pulse = 0; pulse < dd->comm->cd[d].numPulses(); pulse++)
        {
            gpuHaloExchange[d][pulse].reinitHalo(d_x, nullptr);
            gpuHaloExchange[d][pulse].communicateHaloCoordinates(box, &coordinatesReadyOnDeviceEvent);
        }
    }
    // Barrier is needed to avoid other threads using events after its owner has exited and destroyed the context.
    MPI_Barrier(MPI_COMM_WORLD);

    deviceStream.synchronize();

    // Copy results back to host
    copyFromDeviceBuffer(
            h_x->data(), &d_x, 0, numAtomsTotal, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    freeDeviceBuffer(d_x);
#else
    GMX_UNUSED_VALUE(dd);
    GMX_UNUSED_VALUE(box);
    GMX_UNUSED_VALUE(h_x);
    GMX_UNUSED_VALUE(numAtomsTotal);
#endif
}

/*! \brief Define 1D rank topology with 4 MPI tasks
 *
 * \param [in] dd  Domain decomposition object
 */
void define1dRankTopology(gmx_domdec_t* dd)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int numRanks = getNumberOfTestMpiRanks();
    dd->neighbor[0][0] = (rank + 1) % numRanks;
    dd->neighbor[0][1] = (rank == 0) ? (numRanks - 1) : rank - 1;
}

/*! \brief Define 2D rank topology with 4 MPI tasks
 *
 *    -----
 *   | 2 3 |
 *   | 0 1 |
 *    -----
 *
 * \param [in] dd  Domain decomposition object
 */
void define2dRankTopology(gmx_domdec_t* dd)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    switch (rank)
    {
        case 0:
            dd->neighbor[0][0] = 1;
            dd->neighbor[0][1] = 1;
            dd->neighbor[1][0] = 2;
            dd->neighbor[1][1] = 2;
            break;
        case 1:
            dd->neighbor[0][0] = 0;
            dd->neighbor[0][1] = 0;
            dd->neighbor[1][0] = 3;
            dd->neighbor[1][1] = 3;
            break;
        case 2:
            dd->neighbor[0][0] = 3;
            dd->neighbor[0][1] = 3;
            dd->neighbor[1][0] = 0;
            dd->neighbor[1][1] = 0;
            break;
        case 3:
            dd->neighbor[0][0] = 2;
            dd->neighbor[0][1] = 2;
            dd->neighbor[1][0] = 1;
            dd->neighbor[1][1] = 1;
            break;
    }
}

/*! \brief Define a 1D halo with 1 pulses
 *
 * \param [in] dd      Domain decomposition object
 * \param [in] indvec  Vector of index vectors
 */
void define1dHaloWith1Pulse(gmx_domdec_t* dd, std::vector<gmx_domdec_ind_t>* indvec)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> indexvec;
    gmx_domdec_ind_t ind;

    dd->ndim     = 1;
    int nzone    = 1;
    int dimIndex = 0;

    // Set up indices involved in halo
    indexvec.clear();
    indvec->clear();

    dd->comm->cd[dimIndex].receiveInPlace = true;
    dd->dim[dimIndex]                     = 0;
    dd->ci[dimIndex]                      = rank;

    // First pulse involves (arbitrary) indices 1 and 3
    indexvec.push_back(1);
    indexvec.push_back(3);

    ind.index            = indexvec;
    ind.nsend[nzone + 1] = 2;
    ind.nrecv[nzone + 1] = 2;
    indvec->push_back(ind);

    dd->comm->cd[dimIndex].ind = *indvec;
}

/*! \brief Define a 1D halo with 2 pulses
 *
 * \param [in] dd      Domain decomposition object
 * \param [in] indvec  Vector of index vectors
 */
void define1dHaloWith2Pulses(gmx_domdec_t* dd, std::vector<gmx_domdec_ind_t>* indvec)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> indexvec;
    gmx_domdec_ind_t ind;

    dd->ndim     = 1;
    int nzone    = 1;
    int dimIndex = 0;

    // Set up indices involved in halo
    indexvec.clear();
    indvec->clear();

    dd->comm->cd[dimIndex].receiveInPlace = true;
    dd->dim[dimIndex]                     = 0;
    dd->ci[dimIndex]                      = rank;

    // First pulse involves (arbitrary) indices 1 and 3
    indexvec.push_back(1);
    indexvec.push_back(3);

    ind.index            = indexvec;
    ind.nsend[nzone + 1] = 2;
    ind.nrecv[nzone + 1] = 2;
    indvec->push_back(ind);

    // Add another pulse with (arbitrary) indices 4,5,7
    indexvec.clear();

    indexvec.push_back(4);
    indexvec.push_back(5);
    indexvec.push_back(7);

    ind.index            = indexvec;
    ind.nsend[nzone + 1] = 3;
    ind.nrecv[nzone + 1] = 3;
    indvec->push_back(ind);

    dd->comm->cd[dimIndex].ind = *indvec;
}

/*! \brief Define a 2D halo with 1 pulse in each dimension
 *
 * \param [in] dd      Domain decomposition object
 * \param [in] indvec  Vector of index vectors
 */
void define2dHaloWith1PulseInEachDim(gmx_domdec_t* dd, std::vector<gmx_domdec_ind_t>* indvec)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> indexvec;
    gmx_domdec_ind_t ind;

    dd->ndim  = 2;
    int nzone = 1;
    for (int dimIndex = 0; dimIndex < dd->ndim; dimIndex++)
    {

        // Set up indices involved in halo
        indexvec.clear();
        indvec->clear();

        dd->comm->cd[dimIndex].receiveInPlace = true;
        dd->dim[dimIndex]                     = 0;
        dd->ci[dimIndex]                      = rank;

        // Single pulse involving (arbitrary) indices 1 and 3
        indexvec.push_back(1);
        indexvec.push_back(3);

        ind.index            = indexvec;
        ind.nsend[nzone + 1] = 2;
        ind.nrecv[nzone + 1] = 2;
        indvec->push_back(ind);

        dd->comm->cd[dimIndex].ind = *indvec;

        nzone += nzone;
    }
}

/*! \brief Define a 2D halo with 2 pulses in the first dimension
 *
 * \param [in] dd      Domain decomposition object
 * \param [in] indvec  Vector of index vectors
 */
void define2dHaloWith2PulsesInDim1(gmx_domdec_t* dd, std::vector<gmx_domdec_ind_t>* indvec)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> indexvec;
    gmx_domdec_ind_t ind;

    dd->ndim  = 2;
    int nzone = 1;
    for (int dimIndex = 0; dimIndex < dd->ndim; dimIndex++)
    {

        // Set up indices involved in halo
        indexvec.clear();
        indvec->clear();

        dd->comm->cd[dimIndex].receiveInPlace = true;
        dd->dim[dimIndex]                     = 0;
        dd->ci[dimIndex]                      = rank;

        // First pulse involves (arbitrary) indices 1 and 3
        indexvec.push_back(1);
        indexvec.push_back(3);

        ind.index            = indexvec;
        ind.nsend[nzone + 1] = 2;
        ind.nrecv[nzone + 1] = 2;
        indvec->push_back(ind);

        if (dimIndex == 0) // Add another pulse with (arbitrary) indices 4,5,7
        {
            indexvec.clear();

            indexvec.push_back(4);
            indexvec.push_back(5);
            indexvec.push_back(7);

            ind.index            = indexvec;
            ind.nsend[nzone + 1] = 3;
            ind.nrecv[nzone + 1] = 3;
            indvec->push_back(ind);
        }

        dd->comm->cd[dimIndex].ind = *indvec;

        nzone += nzone;
    }
}

/*! \brief Check results for above-defined 1D halo with 1 pulse
 *
 * \param [in] x             Atom coordinate data array
 * \param [in] dd            Domain decomposition object
 * \param [in] numHomeAtoms  Number of home atoms
 */
void checkResults1dHaloWith1Pulse(const RVec* x, const gmx_domdec_t* dd, const int numHomeAtoms)
{
    // Check results are expected from values encoded in x data
    for (int j = 0; j < DIM; j++)
    {
        // First Pulse in first dim: atoms 1 and 3 from forward horizontal neighbour
        EXPECT_EQ(x[numHomeAtoms][j], encodedValue(dd->neighbor[0][0], 1, j));
        EXPECT_EQ(x[numHomeAtoms + 1][j], encodedValue(dd->neighbor[0][0], 3, j));
    }
}

/*! \brief Check results for above-defined 1D halo with 2 pulses
 *
 * \param [in] x             Atom coordinate data array
 * \param [in] dd            Domain decomposition object
 * \param [in] numHomeAtoms  Number of home atoms
 */
void checkResults1dHaloWith2Pulses(const RVec* x, const gmx_domdec_t* dd, const int numHomeAtoms)
{
    // Check results are expected from values encoded in x data
    for (int j = 0; j < DIM; j++)
    {
        // First Pulse in first dim: atoms 1 and 3 from forward horizontal neighbour
        EXPECT_EQ(x[numHomeAtoms][j], encodedValue(dd->neighbor[0][0], 1, j));
        EXPECT_EQ(x[numHomeAtoms + 1][j], encodedValue(dd->neighbor[0][0], 3, j));
        // Second Pulse in first dim: atoms 4,5,7 from forward horizontal neighbour
        EXPECT_EQ(x[numHomeAtoms + 2][j], encodedValue(dd->neighbor[0][0], 4, j));
        EXPECT_EQ(x[numHomeAtoms + 3][j], encodedValue(dd->neighbor[0][0], 5, j));
        EXPECT_EQ(x[numHomeAtoms + 4][j], encodedValue(dd->neighbor[0][0], 7, j));
    }
}

/*! \brief Check results for above-defined 2D halo with 1 pulse in each dimension
 *
 * \param [in] x             Atom coordinate data array
 * \param [in] dd            Domain decomposition object
 * \param [in] numHomeAtoms  Number of home atoms
 */
void checkResults2dHaloWith1PulseInEachDim(const RVec* x, const gmx_domdec_t* dd, const int numHomeAtoms)
{
    // Check results are expected from values encoded in x data
    for (int j = 0; j < DIM; j++)
    {
        // First Pulse in first dim: atoms 1 and 3 from forward horizontal neighbour
        EXPECT_EQ(x[numHomeAtoms][j], encodedValue(dd->neighbor[0][0], 1, j));
        EXPECT_EQ(x[numHomeAtoms + 1][j], encodedValue(dd->neighbor[0][0], 3, j));
        // First Pulse in second dim: atoms 1 and 3 from forward vertical neighbour
        EXPECT_EQ(x[numHomeAtoms + 2][j], encodedValue(dd->neighbor[1][0], 1, j));
        EXPECT_EQ(x[numHomeAtoms + 3][j], encodedValue(dd->neighbor[1][0], 3, j));
    }
}

/*! \brief Check results for above-defined 2D halo with 2 pulses in the first dimension
 *
 * \param [in] x             Atom coordinate data array
 * \param [in] dd            Domain decomposition object
 * \param [in] numHomeAtoms  Number of home atoms
 */
void checkResults2dHaloWith2PulsesInDim1(const RVec* x, const gmx_domdec_t* dd, const int numHomeAtoms)
{
    // Check results are expected from values encoded in x data
    for (int j = 0; j < DIM; j++)
    {
        // First Pulse in first dim: atoms 1 and 3 from forward horizontal neighbour
        EXPECT_EQ(x[numHomeAtoms][j], encodedValue(dd->neighbor[0][0], 1, j));
        EXPECT_EQ(x[numHomeAtoms + 1][j], encodedValue(dd->neighbor[0][0], 3, j));
        // Second Pulse in first dim: atoms 4,5,7 from forward horizontal neighbour
        EXPECT_EQ(x[numHomeAtoms + 2][j], encodedValue(dd->neighbor[0][0], 4, j));
        EXPECT_EQ(x[numHomeAtoms + 3][j], encodedValue(dd->neighbor[0][0], 5, j));
        EXPECT_EQ(x[numHomeAtoms + 4][j], encodedValue(dd->neighbor[0][0], 7, j));
        // First Pulse in second dim: atoms 1 and 3 from forward vertical neighbour
        EXPECT_EQ(x[numHomeAtoms + 5][j], encodedValue(dd->neighbor[1][0], 1, j));
        EXPECT_EQ(x[numHomeAtoms + 6][j], encodedValue(dd->neighbor[1][0], 3, j));
    }
}

TEST(HaloExchangeTest, Coordinates1dHaloWith1Pulse)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    // Set up atom data
    const int        numHomeAtoms  = 10;
    const int        numHaloAtoms  = 2;
    const int        numAtomsTotal = numHomeAtoms + numHaloAtoms;
    HostVector<RVec> h_x;
    h_x.resize(numAtomsTotal);

    initHaloData(h_x.data(), numHomeAtoms, numAtomsTotal);

    // Set up dd
    t_inputrec         ir;
    std::array<int, 1> ddDims = { 0 };
    gmx_domdec_t       dd(ir, ddDims);
    dd.mpi_comm_all              = MPI_COMM_WORLD;
    dd.comm                      = std::make_unique<gmx_domdec_comm_t>();
    dd.unitCellInfo.haveScrewPBC = false;

    DDAtomRanges atomRanges;
    atomRanges.setEnd(DDAtomRanges::Type::Home, numHomeAtoms);
    dd.comm->atomRanges = atomRanges;

    define1dRankTopology(&dd);

    std::vector<gmx_domdec_ind_t> indvec;
    define1dHaloWith1Pulse(&dd, &indvec);

    // Perform halo exchange
    matrix box = { { 0., 0., 0. } };
    dd_move_x(&dd, box, static_cast<ArrayRef<RVec>>(h_x), nullptr);

    // Check results
    checkResults1dHaloWith1Pulse(h_x.data(), &dd, numHomeAtoms);

    if (GMX_GPU_CUDA && GMX_THREAD_MPI) // repeat with GPU halo codepath
    {
        // early return if no devices are available.
        if (getTestHardwareEnvironment()->getTestDeviceList().empty())
        {
            return;
        }

        // Re-initialize input
        initHaloData(h_x.data(), numHomeAtoms, numAtomsTotal);

        // Perform GPU halo exchange
        gpuHalo(&dd, box, &h_x, numAtomsTotal);

        // Check results
        checkResults1dHaloWith1Pulse(h_x.data(), &dd, numHomeAtoms);
    }
}

TEST(HaloExchangeTest, Coordinates1dHaloWith2Pulses)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    // Set up atom data
    const int        numHomeAtoms  = 10;
    const int        numHaloAtoms  = 5;
    const int        numAtomsTotal = numHomeAtoms + numHaloAtoms;
    HostVector<RVec> h_x;
    h_x.resize(numAtomsTotal);

    initHaloData(h_x.data(), numHomeAtoms, numAtomsTotal);

    // Set up dd
    t_inputrec         ir;
    std::array<int, 1> ddDims = { 0 };
    gmx_domdec_t       dd(ir, ddDims);
    dd.mpi_comm_all              = MPI_COMM_WORLD;
    dd.comm                      = std::make_unique<gmx_domdec_comm_t>();
    dd.unitCellInfo.haveScrewPBC = false;

    DDAtomRanges atomRanges;
    atomRanges.setEnd(DDAtomRanges::Type::Home, numHomeAtoms);
    dd.comm->atomRanges = atomRanges;

    define1dRankTopology(&dd);

    std::vector<gmx_domdec_ind_t> indvec;
    define1dHaloWith2Pulses(&dd, &indvec);

    // Perform halo exchange
    matrix box = { { 0., 0., 0. } };
    dd_move_x(&dd, box, static_cast<ArrayRef<RVec>>(h_x), nullptr);

    // Check results
    checkResults1dHaloWith2Pulses(h_x.data(), &dd, numHomeAtoms);

    if (GMX_GPU_CUDA && GMX_THREAD_MPI) // repeat with GPU halo codepath
    {
        // early return if no devices are available.
        if (getTestHardwareEnvironment()->getTestDeviceList().empty())
        {
            return;
        }

        // Re-initialize input
        initHaloData(h_x.data(), numHomeAtoms, numAtomsTotal);

        // Perform GPU halo exchange
        gpuHalo(&dd, box, &h_x, numAtomsTotal);

        // Check results
        checkResults1dHaloWith2Pulses(h_x.data(), &dd, numHomeAtoms);
    }
}


TEST(HaloExchangeTest, Coordinates2dHaloWith1PulseInEachDim)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    // Set up atom data
    const int        numHomeAtoms  = 10;
    const int        numHaloAtoms  = 4;
    const int        numAtomsTotal = numHomeAtoms + numHaloAtoms;
    HostVector<RVec> h_x;
    h_x.resize(numAtomsTotal);

    initHaloData(h_x.data(), numHomeAtoms, numAtomsTotal);

    // Set up dd
    t_inputrec         ir;
    std::array<int, 2> ddDims = { 0, 1 };
    gmx_domdec_t       dd(ir, ddDims);
    dd.mpi_comm_all              = MPI_COMM_WORLD;
    dd.comm                      = std::make_unique<gmx_domdec_comm_t>();
    dd.unitCellInfo.haveScrewPBC = false;

    DDAtomRanges atomRanges;
    atomRanges.setEnd(DDAtomRanges::Type::Home, numHomeAtoms);
    dd.comm->atomRanges = atomRanges;

    define2dRankTopology(&dd);

    std::vector<gmx_domdec_ind_t> indvec;
    define2dHaloWith1PulseInEachDim(&dd, &indvec);

    // Perform halo exchange
    matrix box = { { 0., 0., 0. } };
    dd_move_x(&dd, box, static_cast<ArrayRef<RVec>>(h_x), nullptr);

    // Check results
    checkResults2dHaloWith1PulseInEachDim(h_x.data(), &dd, numHomeAtoms);

    if (GMX_GPU_CUDA && GMX_THREAD_MPI) // repeat with GPU halo codepath
    {
        // early return if no devices are available.
        if (getTestHardwareEnvironment()->getTestDeviceList().empty())
        {
            return;
        }

        // Re-initialize input
        initHaloData(h_x.data(), numHomeAtoms, numAtomsTotal);

        // Perform GPU halo exchange
        gpuHalo(&dd, box, &h_x, numAtomsTotal);

        // Check results
        checkResults2dHaloWith1PulseInEachDim(h_x.data(), &dd, numHomeAtoms);
    }
}

TEST(HaloExchangeTest, Coordinates2dHaloWith2PulsesInDim1)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    // Set up atom data
    const int        numHomeAtoms  = 10;
    const int        numHaloAtoms  = 7;
    const int        numAtomsTotal = numHomeAtoms + numHaloAtoms;
    HostVector<RVec> h_x;
    h_x.resize(numAtomsTotal);

    initHaloData(h_x.data(), numHomeAtoms, numAtomsTotal);

    // Set up dd
    t_inputrec         ir;
    std::array<int, 2> ddDims = { 0, 1 };
    gmx_domdec_t       dd(ir, ddDims);
    dd.mpi_comm_all              = MPI_COMM_WORLD;
    dd.comm                      = std::make_unique<gmx_domdec_comm_t>();
    dd.unitCellInfo.haveScrewPBC = false;

    DDAtomRanges atomRanges;
    atomRanges.setEnd(DDAtomRanges::Type::Home, numHomeAtoms);
    dd.comm->atomRanges = atomRanges;

    define2dRankTopology(&dd);

    std::vector<gmx_domdec_ind_t> indvec;
    define2dHaloWith2PulsesInDim1(&dd, &indvec);

    // Perform halo exchange
    matrix box = { { 0., 0., 0. } };
    dd_move_x(&dd, box, static_cast<ArrayRef<RVec>>(h_x), nullptr);

    // Check results
    checkResults2dHaloWith2PulsesInDim1(h_x.data(), &dd, numHomeAtoms);

    if (GMX_GPU_CUDA && GMX_THREAD_MPI) // repeat with GPU halo codepath
    {
        // early return if no devices are available.
        if (getTestHardwareEnvironment()->getTestDeviceList().empty())
        {
            return;
        }

        // Re-initialize input
        initHaloData(h_x.data(), numHomeAtoms, numAtomsTotal);

        // Perform GPU halo exchange
        gpuHalo(&dd, box, &h_x, numAtomsTotal);

        // Check results
        checkResults2dHaloWith2PulsesInDim1(h_x.data(), &dd, numHomeAtoms);
    }
}

} // namespace
} // namespace test
} // namespace gmx
