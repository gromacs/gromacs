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
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/domdec/atomdistribution.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/gpu_utils/capabilities.h"
#if GMX_GPU
#    include "gromacs/gpu_utils/device_stream.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#endif
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/mpiinfo.h"
#include "gromacs/utility/stringutil.h"

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
 * \param [in]  rank           Rank within MPI communicator
 * \param [out] x              Atom coordinate data array
 * \param [in]  numHomeAtoms   Number of home atoms
 * \param [in]  numAtomsTotal  Total number of atoms, including halo
 */
void initHaloData(const int rank, RVec* x, const int numHomeAtoms, const int numAtomsTotal)
{
    for (int i = 0; i < numAtomsTotal; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            x[i][j] = i < numHomeAtoms ? encodedValue(rank, i, j) : -1;
        }
    }
}

#if GMX_GPU
//! Return whether the given value is the same on all ranks of \c comm
bool valueIsCommonOnAllRanks(const GpuAwareMpiStatus status, const MpiComm& comm)
{
#    if GMX_MPI
    static_assert(
            std::is_same_v<int, std::underlying_type_t<decltype(status)>>,
            "MPI operation expects GPU-aware MPI status enum to use int as the underlying type");
    // Use standard trick to have one collective operation effectively
    // provide the minimum and maximum values across all MPI ranks to
    // each rank, by reducing both the value and its negation separately.
    int valuesToReduce[2] = { -static_cast<int>(status), static_cast<int>(status) };
    MPI_Allreduce(MPI_IN_PLACE, valuesToReduce, 2, MPI_INT, MPI_MIN, comm.comm());
    return valuesToReduce[0] == -valuesToReduce[1];
#    else
    GMX_UNUSED_VALUE(status);
    GMX_UNUSED_VALUE(comm);
    return true
#    endif
}
#endif

/*! \brief Perform GPU halo exchange, including required setup and data transfers
 *
 * \param [in] dd             Domain decomposition object
 * \param [in] box            Box matrix
 * \param [in] h_x            Atom coordinate data array on host
 * \param [in] numAtomsTotal  Total number of atoms, including halo
 *
 * \returns A collection of messages explaining why the halo exchange
 * needs to be skipped, which is empty when the halo exchange was run.
 */
MessageStringCollector gpuHalo(gmx_domdec_t* dd, matrix box, HostVector<RVec>* h_x, int numAtomsTotal)
{
    MessageStringCollector errorReasons;
#if GMX_GPU
    // get communicator
    const MpiComm& mpiComm = dd->mpiComm();
    // Set up GPU hardware environment and assign this MPI rank to a device
    const int         numDevices = getTestHardwareEnvironment()->getTestDeviceList().size();
    const TestDevice* testDevice =
            getTestHardwareEnvironment()->getTestDeviceList()[mpiComm.rank() % numDevices].get();
    const DeviceContext&     deviceContext = testDevice->deviceContext();
    const DeviceInformation& deviceInfo    = testDevice->deviceInfo();
    const GpuAwareMpiStatus  status        = deviceInfo.gpuAwareMpiStatus;
    deviceContext.activate();
    DeviceStream deviceStream(deviceContext, DeviceStreamPriority::Normal, false);

    if (GMX_LIB_MPI)
    {
        errorReasons.appendIf(
                !valueIsCommonOnAllRanks(status, mpiComm),
                formatString(
                        "This rank had GPU-aware MPI support level '%s', but some other rank "
                        "did not. Skipping tests unless all ranks have the same support level.",
                        enumValueToString(status)));
        errorReasons.appendIf(status == GpuAwareMpiStatus::NotSupported,
                              "GPU-aware MPI not supported, so GPU halo exchange cannot be tested");
    }
    // Skip the halo exchange if it cannot work
    if (!errorReasons.isEmpty())
    {
        return errorReasons;
    }

    // pin memory if possible
    changePinningPolicy(h_x, PinningPolicy::PinnedIfSupported);
    // Set up GPU buffer and copy input data from host
    DeviceBuffer<RVec> d_x;
    int                d_x_size       = -1;
    int                d_x_size_alloc = -1;
    reallocateDeviceBuffer(&d_x, numAtomsTotal, &d_x_size, &d_x_size_alloc, deviceContext);

    copyToDeviceBuffer(&d_x, h_x->data(), 0, numAtomsTotal, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    const int numPulses =
            std::accumulate(dd->comm->cd.begin(),
                            dd->comm->cd.end(),
                            0,
                            [](const int a, const auto& b) { return a + b.numPulses(); });
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
            gpuHaloExchange[d].push_back(GpuHaloExchange(
                    dd, d, mpiComm.comm(), mpiComm.comm(), deviceContext, pulse, nullptr));
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
    MPI_Barrier(mpiComm.comm());

    deviceStream.synchronize();

    // Copy results back to host
    copyFromDeviceBuffer(
            h_x->data(), &d_x, 0, numAtomsTotal, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    freeDeviceBuffer(&d_x);
#else
    GMX_UNUSED_VALUE(dd);
    GMX_UNUSED_VALUE(box);
    GMX_UNUSED_VALUE(h_x);
    GMX_UNUSED_VALUE(numAtomsTotal);
#endif
    return errorReasons;
}

/*! \brief Define 1D rank topology with 4 MPI tasks
 *
 * \param [in]  rank  Rank within MPI communicator
 * \param [in]  size  Size of MPI communicator
 * \param [out] dd    Domain decomposition object
 */
void define1dRankTopology(const int rank, const int size, gmx_domdec_t* dd)
{
    dd->neighbor[0][0] = (rank + 1) % size;
    dd->neighbor[0][1] = (rank == 0) ? (size - 1) : rank - 1;
}

/*! \brief Define 2D rank topology with 4 MPI tasks
 *
 *    -----
 *   | 2 3 |
 *   | 0 1 |
 *    -----
 *
 * \param [in]  rank  Rank within MPI communicator
 * \param [out] dd    Domain decomposition object
 */
void define2dRankTopology(const int rank, const int /* size */, gmx_domdec_t* dd)
{
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
 * \param [in]  rank    Rank within MPI communicator
 * \param [out] dd      Domain decomposition object
 * \param [out] indvec  Vector of index vectors
 */
void define1dHaloWith1Pulse(const int rank,
                            const int /* size */,
                            gmx_domdec_t*                  dd,
                            std::vector<gmx_domdec_ind_t>* indvec)
{
    gmx::HostVector<int> indexvec;
    gmx_domdec_ind_t     ind;

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
 * \param [in]  rank    Rank within MPI communicator
 * \param [out] dd      Domain decomposition object
 * \param [out] indvec  Vector of index vectors
 */
void define1dHaloWith2Pulses(const int rank,
                             const int /* size */,
                             gmx_domdec_t*                  dd,
                             std::vector<gmx_domdec_ind_t>* indvec)
{
    gmx::HostVector<int> indexvec;
    gmx_domdec_ind_t     ind;

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
 * \param [in]  rank    Rank within MPI communicator
 * \param [out] dd      Domain decomposition object
 * \param [out] indvec  Vector of index vectors
 */
void define2dHaloWith1PulseInEachDim(const int rank,
                                     const int /*size*/,
                                     gmx_domdec_t*                  dd,
                                     std::vector<gmx_domdec_ind_t>* indvec)
{
    gmx::HostVector<int> indexvec;
    gmx_domdec_ind_t     ind;

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
 * \param [in]  rank    Rank within MPI communicator
 * \param [out] dd      Domain decomposition object
 * \param [out] indvec  Vector of index vectors
 */
void define2dHaloWith2PulsesInDim1(const int rank,
                                   const int /* size */,
                                   gmx_domdec_t*                  dd,
                                   std::vector<gmx_domdec_ind_t>* indvec)
{
    HostVector<int>  indexvec;
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

//! Parameters over which halo exchange is tested
struct HaloExchangeTestParameters
{
    //! Human-readable description for scoped traces
    std::string description;
    //! Number of home atoms
    const int numHomeAtoms;
    //! Number of halo atoms
    const int numHaloAtoms;
    //! The dimensions of the domain decomposition
    std::vector<int> ddDims;
    //! The DD topology setup function to use
    void (*domainTopologySetupFunction)(int, int, gmx_domdec_t*);
    //! The DD pulse-structure setup function to use
    void (*pulseSetupFunction)(int, int, gmx_domdec_t*, std::vector<gmx_domdec_ind_t>*);
    //! The results-checking funciton to use
    void (*checkResults)(const RVec*, const gmx_domdec_t*, const int);
};

/*! \brief Data used in test body for halo exchange
 *
 * This data cannot be a member of the test fixture class because with
 * thread-MPI then they would be shared across all ranks, which can't
 * work. */
class HaloExchangeTestData
{
public:
    HaloExchangeTestData(const HaloExchangeTestParameters& parameters) :
        dd_{ mpiComm_, ir_, parameters.ddDims },
        numAtomsTotal_{ parameters.numHomeAtoms + parameters.numHaloAtoms }
    {
        SCOPED_TRACE("Testing " + parameters.description);
        dd_.comm                      = std::make_unique<gmx_domdec_comm_t>(mpiComm_);
        dd_.unitCellInfo.haveScrewPBC = false;

        DDAtomRanges atomRanges;
        atomRanges.setEnd(DDAtomRanges::Type::Home, parameters.numHomeAtoms);
        dd_.comm->atomRanges = atomRanges;

        parameters.domainTopologySetupFunction(mpiComm_.rank(), mpiComm_.size(), &dd_);
        parameters.pulseSetupFunction(mpiComm_.rank(), mpiComm_.size(), &dd_, &indvec_);
    }
    //! MPI communicator
    MpiComm mpiComm_{ MPI_COMM_WORLD };
    //! Input record
    t_inputrec ir_;
    //! DD manager
    gmx_domdec_t dd_;
    //! Describes DD pulse structure
    std::vector<gmx_domdec_ind_t> indvec_;
    //! Total number of atoms known to each domain
    int numAtomsTotal_;
    //! Position coordinates
    HostVector<RVec> h_x_{ static_cast<size_t>(numAtomsTotal_) };
};

//! Test fixture for halo exchange
class HaloExchangeTest : public ::testing::TestWithParam<HaloExchangeTestParameters>
{
public:
    //! Box matrix (unused by halo exchange in practice)
    matrix box_ = { { 0., 0., 0. } };
};

TEST_P(HaloExchangeTest, WithParametersOnCpu)
{
    SCOPED_TRACE("Testing " + GetParam().description);
    GMX_MPI_TEST(RequireRankCount<4>);

    HaloExchangeTestData data(GetParam());
    const int            numHomeAtoms = GetParam().numHomeAtoms;
    initHaloData(data.mpiComm_.rank(), data.h_x_.data(), numHomeAtoms, data.numAtomsTotal_);
    dd_move_x(&data.dd_, box_, static_cast<ArrayRef<RVec>>(data.h_x_), nullptr);
    GetParam().checkResults(data.h_x_.data(), &data.dd_, numHomeAtoms);
}

TEST_P(HaloExchangeTest, WithParametersOnGpu)
{
    SCOPED_TRACE("Testing " + GetParam().description);
    GMX_MPI_TEST(RequireRankCount<4>);

    if (GMX_THREAD_MPI && !GpuConfigurationCapabilities::HaloExchangeDirectComm)
    {
        GTEST_SKIP() << "With thread-MPI, GPU halo exchange is only supported on CUDA";
    }
    if (GMX_LIB_MPI && !GpuConfigurationCapabilities::HaloExchangeDirectComm)
    {
        GTEST_SKIP() << "With library MPI, GPU halo exchange is not supported for this build";
    }
    if (getTestHardwareEnvironment()->getTestDeviceList().empty())
    {
        GTEST_SKIP() << "No GPUs detected";
    }

    HaloExchangeTestData data(GetParam());
    const int            numHomeAtoms = GetParam().numHomeAtoms;
    initHaloData(data.mpiComm_.rank(), data.h_x_.data(), numHomeAtoms, data.numAtomsTotal_);
    MessageStringCollector skipReasons = gpuHalo(&data.dd_, box_, &data.h_x_, data.numAtomsTotal_);
    if (skipReasons.isEmpty())
    {
        // Halo exchange ran, so check the results
        GetParam().checkResults(data.h_x_.data(), &data.dd_, numHomeAtoms);
    }
    else
    {
        GTEST_SKIP() << skipReasons.toString();
    }
}

static const std::vector<HaloExchangeTestParameters> c_testSetups = {
    { "1D halo with 1 pulse", 10, 2, { 0 }, define1dRankTopology, define1dHaloWith1Pulse, checkResults1dHaloWith1Pulse },
    { "1D halo with 2 pulses", 10, 5, { 0 }, define1dRankTopology, define1dHaloWith2Pulses, checkResults1dHaloWith2Pulses },
    { "2D halo with 1 pulse in each dimension",
      10,
      4,
      { 0, 1 },
      define2dRankTopology,
      define2dHaloWith1PulseInEachDim,
      checkResults2dHaloWith1PulseInEachDim },
    { "2D halo with 2 pulses in first dimension",
      10,
      7,
      { 0, 1 },
      define2dRankTopology,
      define2dHaloWith2PulsesInDim1,
      checkResults2dHaloWith2PulsesInDim1 },
};

INSTANTIATE_TEST_SUITE_P(Works, HaloExchangeTest, ::testing::ValuesIn(c_testSetups));

} // namespace
} // namespace test
} // namespace gmx
