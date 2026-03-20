/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 *
 * \brief This file contains function definitions necessary for
 * managing the offload of long-ranged PME work to separate MPI rank,
 * for computing energies and forces (Coulomb and LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "pme_pp.h"

#include "config.h"

#include <cstdio>

#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/vec.h"

namespace gmx
{

/*! \brief Block to wait for communication to PME ranks to complete
 *
 * This should be faster with a real non-blocking MPI implementation
 */
static constexpr bool c_useDelayedWait = false;

PmePpComm::PmePpComm(const MpiComm&               comm,
                     const int                    rankOfPartnerPmeRank,
                     const VanDerWaalsType        vdwType,
                     const CoulombInteractionType coulombType,
                     const bool                   thisRankReceivesVirialAndEnergy,
                     const bool                   useGpuPmePpCommunication,
                     const bool                   useNvshmem,
                     const DeviceStreamManager*   deviceStreamManager) :
    comm_(comm),
    rankOfPartnerPmeRank_(rankOfPartnerPmeRank),
    vdwType_(vdwType),
    coulombType_(coulombType),
    thisRankReceivesVirialAndEnergy_(thisRankReceivesVirialAndEnergy),
    useGpuPmePpCommunication_(useGpuPmePpCommunication),
    cnb_{ thisRankReceivesVirialAndEnergy_ ? std::make_optional<gmx_pme_comm_n_box_t>() : std::nullopt },
    cpuPmeForceReceiveBuffer_{ HostAllocationPolicy{
            useGpuPmePpCommunication ? PinningPolicy::PinnedIfSupported : PinningPolicy::CannotBePinned } }
{
    if (useGpuPmePpCommunication_)
    {
        GMX_RELEASE_ASSERT(
                deviceStreamManager != nullptr,
                "GPU device stream manager should be valid in order to use PME-PP direct "
                "communications.");
        GMX_RELEASE_ASSERT(deviceStreamManager->streamIsValid(DeviceStreamType::PmePpTransfer),
                           "GPU PP-PME stream should be valid in order to use GPU PME-PP direct "
                           "communications.");
        pmePpCommGpu_.emplace(comm_.comm(),
                              rankOfPartnerPmeRank_,
                              deviceStreamManager->context(),
                              deviceStreamManager->stream(DeviceStreamType::PmePpTransfer),
                              useNvshmem);
    }
}

/*! \brief Wait for the pending data send requests to PME ranks to complete */
static void waitForSentData(std::vector<MPI_Request>* requests)
{
    if (!requests->empty())
    {
#if GMX_MPI
        MPI_Waitall(requests->size(), requests->data(), MPI_STATUSES_IGNORE);
#endif
        requests->clear();
    }
}

void PmePpComm::sendParameters(const int            numHomeAtoms,
                               const bool           bFreeEnergy_q,
                               const bool           bFreeEnergy_lj,
                               ArrayRef<const real> chargeA,
                               ArrayRef<const real> chargeB,
                               ArrayRef<const real> sqrt_c6A,
                               ArrayRef<const real> sqrt_c6B,
                               ArrayRef<const real> sigmaA,
                               ArrayRef<const real> sigmaB,
                               const int            maxshift_x,
                               const int            maxshift_y)
{
    numHomeAtoms_ = numHomeAtoms;

    unsigned int flags = 0;

    if (usingPme(coulombType_))
    {
        flags |= PP_PME_CHARGE;
    }
    if (usingLJPme(vdwType_))
    {
        flags |= (PP_PME_SQRTC6 | PP_PME_SIGMA);
    }
    if (bFreeEnergy_q || bFreeEnergy_lj)
    {
        // Require that the B-state flags are in the bits just above
        // the ones for the A state, ie. double the value.
        static_assert(PP_PME_CHARGEB == PP_PME_CHARGE << 1,
                      "PP-PME communication flag assumption violated");
        static_assert(PP_PME_SQRTC6B == PP_PME_SQRTC6 << 1,
                      "PP-PME communication flag assumption violated");
        static_assert(PP_PME_SIGMAB == PP_PME_SIGMA << 1,
                      "PP-PME communication flag assumption violated");
        flags |= (flags << 1);
    }

    if (debug)
    {
        fprintf(debug,
                "PP rank %d sending to PME rank %d: %d %s%s%s\n",
                comm_.rank(),
                rankOfPartnerPmeRank_,
                numHomeAtoms_,
                (flags & PP_PME_CHARGE) ? "charges " : "",
                (flags & PP_PME_SQRTC6) ? "sqrtC6 " : "",
                (flags & PP_PME_SIGMA) ? " sigma " : "");
    }

    if (c_useDelayedWait)
    {
        waitForSentData(&requests_);
    }

#if GMX_MPI
    MPI_Comm comm = comm_.comm();
#endif

    if (thisRankReceivesVirialAndEnergy_)
    {
        // Peer PP node: communicate all per-domain-lifetime data
        cnb_->flags      = flags;
        cnb_->natoms     = numHomeAtoms_;
        cnb_->maxshift_x = maxshift_x;
        cnb_->maxshift_y = maxshift_y;
#if GMX_MPI
        requests_.push_back(MPI_Request{});
        MPI_Isend(&cnb_.value(),
                  sizeof(cnb_.value()),
                  MPI_BYTE,
                  rankOfPartnerPmeRank_,
                  eCommType_CNB,
                  comm,
                  &requests_.back());
#endif
    }
    else
    {
#if GMX_MPI
        // Communicate only the number of atoms
        requests_.push_back(MPI_Request{});
        MPI_Isend(&numHomeAtoms_,
                  sizeof(numHomeAtoms_),
                  MPI_BYTE,
                  rankOfPartnerPmeRank_,
                  eCommType_CNB,
                  comm,
                  &requests_.back());
#endif
    }

#if GMX_MPI
    if (numHomeAtoms_ > 0)
    {
        if (flags & PP_PME_CHARGE)
        {
            requests_.push_back(MPI_Request{});
            GMX_ASSERT(gmx::ssize(chargeA) >= numHomeAtoms_, "A-state charge send buffer too small");
            MPI_Isend(chargeA.data(),
                      numHomeAtoms_ * sizeof(chargeA[0]),
                      MPI_BYTE,
                      rankOfPartnerPmeRank_,
                      eCommType_ChargeA,
                      comm,
                      &requests_.back());
        }
        if (flags & PP_PME_CHARGEB)
        {
            requests_.push_back(MPI_Request{});
            GMX_ASSERT(gmx::ssize(chargeB) >= numHomeAtoms_, "B-state charge send buffer too small");
            MPI_Isend(chargeB.data(),
                      numHomeAtoms_ * sizeof(chargeB[0]),
                      MPI_BYTE,
                      rankOfPartnerPmeRank_,
                      eCommType_ChargeB,
                      comm,
                      &requests_.back());
        }
        if (flags & PP_PME_SQRTC6)
        {
            requests_.push_back(MPI_Request{});
            GMX_ASSERT(gmx::ssize(sqrt_c6A) >= numHomeAtoms_,
                       "A-state sqrt C6 send buffer too small");
            MPI_Isend(sqrt_c6A.data(),
                      numHomeAtoms_ * sizeof(sqrt_c6A[0]),
                      MPI_BYTE,
                      rankOfPartnerPmeRank_,
                      eCommType_SQRTC6A,
                      comm,
                      &requests_.back());
        }
        if (flags & PP_PME_SQRTC6B)
        {
            requests_.push_back(MPI_Request{});
            GMX_ASSERT(gmx::ssize(sqrt_c6B) >= numHomeAtoms_,
                       "B-state sqrt C6 send buffer too small");
            MPI_Isend(sqrt_c6B.data(),
                      numHomeAtoms_ * sizeof(sqrt_c6B[0]),
                      MPI_BYTE,
                      rankOfPartnerPmeRank_,
                      eCommType_SQRTC6B,
                      comm,
                      &requests_.back());
        }
        if (flags & PP_PME_SIGMA)
        {
            requests_.push_back(MPI_Request{});
            GMX_ASSERT(gmx::ssize(sigmaA) >= numHomeAtoms_, "A-state sigma send buffer too small");
            MPI_Isend(sigmaA.data(),
                      numHomeAtoms_ * sizeof(sigmaA[0]),
                      MPI_BYTE,
                      rankOfPartnerPmeRank_,
                      eCommType_SigmaA,
                      comm,
                      &requests_.back());
        }
        if (flags & PP_PME_SIGMAB)
        {
            requests_.push_back(MPI_Request{});
            GMX_ASSERT(gmx::ssize(sigmaB) >= numHomeAtoms_, "B-state sigma send buffer too small");
            MPI_Isend(sigmaB.data(),
                      numHomeAtoms_ * sizeof(sigmaB[0]),
                      MPI_BYTE,
                      rankOfPartnerPmeRank_,
                      eCommType_SigmaB,
                      comm,
                      &requests_.back());
        }
    }
#else
    GMX_UNUSED_VALUE(chargeA);
    GMX_UNUSED_VALUE(chargeB);
    GMX_UNUSED_VALUE(sqrt_c6A);
    GMX_UNUSED_VALUE(sqrt_c6B);
    GMX_UNUSED_VALUE(sigmaA);
    GMX_UNUSED_VALUE(sigmaB);
#endif

    if (!c_useDelayedWait)
    {
        // Wait for pending communication to finish
        waitForSentData(&requests_);
    }
}

void PmePpComm::sendCoordinates(DeviceBuffer<RVec>    coordinates,
                                const matrix          box,
                                ArrayRef<const RVec>  x,
                                const real            lambda_q,
                                const real            lambda_lj,
                                const bool            computeEnergyAndVirial,
                                const int64_t         step,
                                const bool            reinitGpuPmePpComms,
                                const bool            sendCoordinatesFromGpu,
                                const bool            receiveForcesToGpu,
                                GpuEventSynchronizer* coordinatesReadyOnDeviceEvent,
                                const bool            useMdGpuGraph,
                                gmx_wallcycle*        wcycle)
{
    wallcycle_start(wcycle, WallCycleCounter::PpPmeSendX);

    unsigned int flags = PP_PME_COORD;
    if (computeEnergyAndVirial)
    {
        flags |= PP_PME_ENER_VIR;
    }

    if (debug)
    {
        fprintf(debug,
                "PP rank %d sending to PME rank %d: %d %s\n",
                comm_.rank(),
                rankOfPartnerPmeRank_,
                numHomeAtoms_,
                (flags & PP_PME_COORD) ? "coordinates" : "");
    }

    if (receiveForcesToGpu)
    {
        flags |= PP_PME_RECVFTOGPU;
    }

    if (useMdGpuGraph)
    {
        flags |= PP_PME_MDGPUGRAPH;
    }

    if (c_useDelayedWait)
    {
        waitForSentData(&requests_);
    }

#if GMX_MPI
    MPI_Comm comm = comm_.comm();
#endif

    if (thisRankReceivesVirialAndEnergy_)
    {
        // Peer PP node: communicate all per-step data
        cnb_->flags     = flags;
        cnb_->natoms    = numHomeAtoms_;
        cnb_->lambda_q  = lambda_q;
        cnb_->lambda_lj = lambda_lj;
        cnb_->step      = step;
        copy_mat(box, cnb_->box);
#if GMX_MPI
        requests_.push_back(MPI_Request{});
        MPI_Isend(&cnb_.value(),
                  sizeof(cnb_.value()),
                  MPI_BYTE,
                  rankOfPartnerPmeRank_,
                  eCommType_CNB,
                  comm,
                  &requests_.back());
#endif
    }
    // No CNB message from other PP ranks

    // With direct-GPU PME-PP communication, coordinates and
    // forces are always transferred each step, even for empty
    // domains.

    // Resize the PME force-receive buffer on the CPU.
    //
    // TODO If/when this moves, adjust the related comment in
    // receiveForces()
    cpuPmeForceReceiveBuffer_.resize(numHomeAtoms_);
    if (reinitGpuPmePpComms)
    {
        GMX_ASSERT(useGpuPmePpCommunication_, "Flag mismatch, cannot reinitialize PME-PP comms");
        // When using GPU-direct PP-PME communication, notify the
        // receiver where it should store forces when they arrive
        // on the CPU.
        pmePpCommGpu_->reinit(cpuPmeForceReceiveBuffer_);
    }

    // Ensure that with GPU PME-PP comms, even empty domains still
    // have the chance to post a possible non-blocking matching
    // force receive, since the PME rank always transfers
    // forces to each PP rank.
    if (useGpuPmePpCommunication_)
    {
        if (sendCoordinatesFromGpu)
        {
            GMX_ASSERT(coordinatesReadyOnDeviceEvent != nullptr,
                       "When sending coordinates from GPU, a synchronization event should "
                       "be provided");
            pmePpCommGpu_->sendCoordinatesToPmeFromGpu(
                    coordinates, numHomeAtoms_, coordinatesReadyOnDeviceEvent, receiveForcesToGpu);
        }
        else
        {
            pmePpCommGpu_->sendCoordinatesToPmeFromCpu(x.data(), numHomeAtoms_, receiveForcesToGpu);
        }
    }
    else if (numHomeAtoms_ > 0)
    {
#if GMX_MPI
        // With CPU PME-PP comms, coordinate messages are only
        // sent for non-empty domains.
        requests_.push_back(MPI_Request{});
        GMX_ASSERT(gmx::ssize(x) >= numHomeAtoms_, "Position send buffer too small");
        MPI_Isend(x.data(),
                  numHomeAtoms_ * sizeof(x[0]),
                  MPI_BYTE,
                  rankOfPartnerPmeRank_,
                  eCommType_COORD,
                  comm,
                  &requests_.back());
#endif
    }

    if (!c_useDelayedWait)
    {
        // We can delay this wait as we are sure x and q will not be modified
        // before the next call to send coordinates or receive results
        waitForSentData(&requests_);
    }

    wallcycle_stop(wcycle, WallCycleCounter::PpPmeSendX);
}

void PmePpComm::sendFinish() const
{
    // Only let one PP rank signal each PME rank
    if (thisRankReceivesVirialAndEnergy_)
    {
        // Blocking MPI call is used, can use a stack variable
        gmx_pme_comm_n_box_t cnb;
        cnb.flags = PP_PME_FINISH;

#if GMX_MPI
        // We send this, uncommon, message blocking to simplify the code
        MPI_Send(&cnb, sizeof(cnb), MPI_BYTE, rankOfPartnerPmeRank_, eCommType_CNB, comm_.comm());
#endif
    }
}

void PmePpComm::sendSwitchGrid(ivec grid_size, const real ewaldcoeff_q, const real ewaldcoeff_lj) const
{
    // Only let one PP rank signal each PME rank
    if (thisRankReceivesVirialAndEnergy_)
    {
        // Blocking MPI call is used, can use a stack variable
        gmx_pme_comm_n_box_t cnb;
        cnb.flags = PP_PME_SWITCHGRID;
        copy_ivec(grid_size, cnb.grid_size);
        cnb.ewaldcoeff_q  = ewaldcoeff_q;
        cnb.ewaldcoeff_lj = ewaldcoeff_lj;

#if GMX_MPI
        // We send this, uncommon, message blocking to simplify the code
        MPI_Send(&cnb, sizeof(cnb), MPI_BYTE, rankOfPartnerPmeRank_, eCommType_CNB, comm_.comm());
#endif
    }
}

void PmePpComm::sendResetCounters(const int64_t step) const
{
    // Only let one PP rank signal each PME rank
    if (thisRankReceivesVirialAndEnergy_)
    {
        // Blocking MPI call is used, can use a stack variable
        gmx_pme_comm_n_box_t cnb;

        cnb.flags = PP_PME_RESETCOUNTERS;
        cnb.step  = step;

#if GMX_MPI
        // We send this, uncommon, message blocking to simplify the code
        MPI_Send(&cnb, sizeof(cnb), MPI_BYTE, rankOfPartnerPmeRank_, eCommType_CNB, comm_.comm());
#endif
    }
}

void PmePpComm::receiveVirialAndEnergy(ForceWithVirial* forceWithVirial,
                                       real*            energy_q,
                                       real*            energy_lj,
                                       real*            dvdlambda_q,
                                       real*            dvdlambda_lj,
                                       float*           pme_cycles)
{
    gmx_pme_comm_vir_ene_t cve;

    if (thisRankReceivesVirialAndEnergy_)
    {
        if (debug)
        {
            fprintf(debug,
                    "PP rank %d receiving from PME rank %d: virial and energy\n",
                    comm_.rank(),
                    rankOfPartnerPmeRank_);
        }
#if GMX_MPI
        MPI_Recv(&cve,
                 sizeof(cve),
                 MPI_BYTE,
                 rankOfPartnerPmeRank_,
                 eCommType_ENERGY_VIRIAL_DVDL,
                 comm_.comm(),
                 MPI_STATUS_IGNORE);
#else
        std::memset(&cve, 0, sizeof(cve));
#endif

        forceWithVirial->addVirialContribution(cve.vir_q);
        forceWithVirial->addVirialContribution(cve.vir_lj);
        *energy_q  = cve.energy_q;
        *energy_lj = cve.energy_lj;
        *dvdlambda_q += cve.dvdlambda_q;
        *dvdlambda_lj += cve.dvdlambda_lj;
        *pme_cycles = cve.cycles;

        if (cve.stop_cond != StopCondition::None)
        {
            gmx_set_stop_condition(cve.stop_cond);
        }
    }
    else
    {
        *energy_q   = 0;
        *energy_lj  = 0;
        *pme_cycles = 0;
    }
}

/*! \brief Receive force data from PME ranks */
void PmePpComm::receiveForces(const bool receivePmeForceToGpu)
{
    // With all kinds of PME-PP communication, forces are always
    // returned each step, even to empty domains.
    if (useGpuPmePpCommunication_)
    {
        // Receive forces from PME rank to cpuPmeForceReceiveBuffer_
        // (which was set up at repartition time).
        pmePpCommGpu_->receiveForceFromPme(receivePmeForceToGpu);
    }
    else
    {
        // Receive data using MPI
#if GMX_MPI
        GMX_ASSERT(gmx::ssize(cpuPmeForceReceiveBuffer_) >= numHomeAtoms_,
                   "CPU PME force receive buffer too small");
        MPI_Recv(cpuPmeForceReceiveBuffer_.data(),
                 numHomeAtoms_ * sizeof(cpuPmeForceReceiveBuffer_[0]),
                 MPI_BYTE,
                 rankOfPartnerPmeRank_,
                 eCommType_FORCES,
                 comm_.comm(),
                 MPI_STATUS_IGNORE);
#endif
    }
}


void PmePpComm::receiveResults(ForceWithVirial* forceWithVirial,
                               real*            energy_q,
                               real*            energy_lj,
                               real*            dvdlambda_q,
                               real*            dvdlambda_lj,
                               const bool       receivePmeForceToGpu,
                               float*           pme_cycles)
{
    if (c_useDelayedWait)
    {
        /* Wait for the x request to finish */
        waitForSentData(&requests_);
    }

    receiveForces(receivePmeForceToGpu);

    int nt = gmx_omp_nthreads_get_simple_rvec_task(ModuleMultiThread::Default, numHomeAtoms_);

    ArrayRef<RVec> f = forceWithVirial->force_;

    if (!receivePmeForceToGpu)
    {
        /* Note that we would like to avoid this conditional by putting it
         * into the omp pragma instead, but then we still take the full
         * omp parallel for overhead (at least with gcc5).
         */
        if (nt == 1)
        {
            for (int i = 0; i < numHomeAtoms_; i++)
            {
                f[i] += cpuPmeForceReceiveBuffer_[i];
            }
        }
        else
        {
#pragma omp parallel for num_threads(nt) schedule(static)
            for (int i = 0; i < numHomeAtoms_; i++)
            {
                f[i] += cpuPmeForceReceiveBuffer_[i];
            }
        }
    }

    receiveVirialAndEnergy(forceWithVirial, energy_q, energy_lj, dvdlambda_q, dvdlambda_lj, pme_cycles);
}

std::optional<int> PmePpComm::rankOfControlledPmeRank() const
{
    return thisRankReceivesVirialAndEnergy_ ? std::make_optional<int>(rankOfPartnerPmeRank_) : std::nullopt;
}

PmePpCommGpu* PmePpComm::pmePpCommGpu()
{
    if (pmePpCommGpu_.has_value())
    {
        return &pmePpCommGpu_.value();
    }
    return nullptr;
}

} // namespace gmx
