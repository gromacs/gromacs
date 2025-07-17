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
#include <cstring>

#include <memory>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vec.h"

#include "pme_pp_communication.h"

/*! \brief Block to wait for communication to PME ranks to complete
 *
 * This should be faster with a real non-blocking MPI implementation
 */
static constexpr bool c_useDelayedWait = false;

/*! \brief Wait for the pending data send requests to PME ranks to complete */
static void gmx_pme_send_coeffs_coords_wait(gmx_domdec_t* dd)
{
    if (dd->nreq_pme)
    {
#if GMX_MPI
        MPI_Waitall(dd->nreq_pme, dd->req_pme, MPI_STATUSES_IGNORE);
#endif
        dd->nreq_pme = 0;
    }
}

/*! \brief Send data to PME ranks */
static void gmx_pme_send_coeffs_coords(t_forcerec*                    fr,
                                       const gmx::MpiComm&            mpiCommMySim,
                                       gmx_domdec_t*                  dd,
                                       unsigned int                   flags,
                                       gmx::ArrayRef<const real>      chargeA,
                                       gmx::ArrayRef<const real>      chargeB,
                                       gmx::ArrayRef<const real>      c6A,
                                       gmx::ArrayRef<const real>      c6B,
                                       gmx::ArrayRef<const real>      sigmaA,
                                       gmx::ArrayRef<const real>      sigmaB,
                                       const matrix                   box,
                                       gmx::ArrayRef<const gmx::RVec> x,
                                       real                           lambda_q,
                                       real                           lambda_lj,
                                       int                            maxshift_x,
                                       int                            maxshift_y,
                                       int64_t                        step,
                                       bool                           useGpuPmePpComms,
                                       bool                           reinitGpuPmePpComms,
                                       bool                           sendCoordinatesFromGpu,
                                       bool                           receiveForcesToGpu,
                                       bool                           useMdGpuGraph,
                                       GpuEventSynchronizer*          coordinatesReadyOnDeviceEvent)
{
    const int n = dd_numHomeAtoms(*dd);

    if (debug)
    {
        fprintf(debug,
                "PP rank %d sending to PME rank %d: %d%s%s%s%s\n",
                mpiCommMySim.rank(),
                dd->pme_nodeid,
                n,
                (flags & PP_PME_CHARGE) ? " charges" : "",
                (flags & PP_PME_SQRTC6) ? " sqrtC6" : "",
                (flags & PP_PME_SIGMA) ? " sigma" : "",
                (flags & PP_PME_COORD) ? " coordinates" : "");
    }

    if (useGpuPmePpComms)
    {
        flags |= PP_PME_GPUCOMMS;
        if (receiveForcesToGpu)
        {
            flags |= PP_PME_RECVFTOGPU;
        }
    }

    if (useMdGpuGraph)
    {
        flags |= PP_PME_MDGPUGRAPH;
    }

    if (c_useDelayedWait)
    {
        /* We can not use cnb until pending communication has finished */
        gmx_pme_send_coeffs_coords_wait(dd);
    }

    if (dd->pme_receive_vir_ener)
    {
        /* Peer PP node: communicate all data */
        if (dd->cnb == nullptr)
        {
            snew(dd->cnb, 1);
        }
        gmx_pme_comm_n_box_t* cnb = dd->cnb;

        cnb->flags      = flags;
        cnb->natoms     = n;
        cnb->maxshift_x = maxshift_x;
        cnb->maxshift_y = maxshift_y;
        cnb->lambda_q   = lambda_q;
        cnb->lambda_lj  = lambda_lj;
        cnb->step       = step;
        if (flags & PP_PME_COORD)
        {
            copy_mat(box, cnb->box);
        }
#if GMX_MPI
        MPI_Isend(cnb,
                  sizeof(*cnb),
                  MPI_BYTE,
                  dd->pme_nodeid,
                  eCommType_CNB,
                  mpiCommMySim.comm(),
                  &dd->req_pme[dd->nreq_pme++]);
#endif
    }
    else if (flags & (PP_PME_CHARGE | PP_PME_SQRTC6 | PP_PME_SIGMA))
    {
#if GMX_MPI
        /* Communicate only the number of atoms */
        MPI_Isend(&n,
                  sizeof(n),
                  MPI_BYTE,
                  dd->pme_nodeid,
                  eCommType_CNB,
                  mpiCommMySim.comm(),
                  &dd->req_pme[dd->nreq_pme++]);
#endif
    }

#if GMX_MPI
    if (n > 0)
    {
        if (flags & PP_PME_CHARGE)
        {
            MPI_Isend(chargeA.data(),
                      n * sizeof(real),
                      MPI_BYTE,
                      dd->pme_nodeid,
                      eCommType_ChargeA,
                      mpiCommMySim.comm(),
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_CHARGEB)
        {
            MPI_Isend(chargeB.data(),
                      n * sizeof(real),
                      MPI_BYTE,
                      dd->pme_nodeid,
                      eCommType_ChargeB,
                      mpiCommMySim.comm(),
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SQRTC6)
        {
            MPI_Isend(c6A.data(),
                      n * sizeof(real),
                      MPI_BYTE,
                      dd->pme_nodeid,
                      eCommType_SQRTC6A,
                      mpiCommMySim.comm(),
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SQRTC6B)
        {
            MPI_Isend(c6B.data(),
                      n * sizeof(real),
                      MPI_BYTE,
                      dd->pme_nodeid,
                      eCommType_SQRTC6B,
                      mpiCommMySim.comm(),
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SIGMA)
        {
            MPI_Isend(sigmaA.data(),
                      n * sizeof(real),
                      MPI_BYTE,
                      dd->pme_nodeid,
                      eCommType_SigmaA,
                      mpiCommMySim.comm(),
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SIGMAB)
        {
            MPI_Isend(sigmaB.data(),
                      n * sizeof(real),
                      MPI_BYTE,
                      dd->pme_nodeid,
                      eCommType_SigmaB,
                      mpiCommMySim.comm(),
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_COORD)
        {
            if (reinitGpuPmePpComms)
            {
                changePinningPolicy(&dd->pmeForceReceiveBuffer, gmx::PinningPolicy::PinnedIfSupported);
                dd->pmeForceReceiveBuffer.resize(n);
                fr->pmePpCommGpu->reinit(n);
            }

            if (useGpuPmePpComms && (fr != nullptr))
            {
                if (sendCoordinatesFromGpu)
                {
                    GMX_ASSERT(coordinatesReadyOnDeviceEvent != nullptr,
                               "When sending coordinates from GPU, a synchronization event should "
                               "be provided");
                    fr->pmePpCommGpu->sendCoordinatesToPmeFromGpu(
                            fr->stateGpu->getCoordinates(), n, coordinatesReadyOnDeviceEvent);
                }
                else
                {
                    fr->pmePpCommGpu->sendCoordinatesToPmeFromCpu(const_cast<gmx::RVec*>(x.data()), n);
                }
            }
            else
            {
                MPI_Isend(x.data(),
                          n * sizeof(rvec),
                          MPI_BYTE,
                          dd->pme_nodeid,
                          eCommType_COORD,
                          mpiCommMySim.comm(),
                          &dd->req_pme[dd->nreq_pme++]);
            }
        }
    }
#else
    GMX_UNUSED_VALUE(fr);
    GMX_UNUSED_VALUE(chargeA);
    GMX_UNUSED_VALUE(chargeB);
    GMX_UNUSED_VALUE(c6A);
    GMX_UNUSED_VALUE(c6B);
    GMX_UNUSED_VALUE(sigmaA);
    GMX_UNUSED_VALUE(sigmaB);
    GMX_UNUSED_VALUE(x);
    GMX_UNUSED_VALUE(reinitGpuPmePpComms);
    GMX_UNUSED_VALUE(sendCoordinatesFromGpu);
    GMX_UNUSED_VALUE(coordinatesReadyOnDeviceEvent);
#endif
    if (!c_useDelayedWait)
    {
        /* Wait for the data to arrive */
        /* We can skip this wait as we are sure x and q will not be modified
         * before the next call to gmx_pme_send_x_q or gmx_pme_receive_f.
         */
        gmx_pme_send_coeffs_coords_wait(dd);
    }
}

void gmx_pme_send_parameters(const gmx::MpiComm&        mpiCommMySim,
                             gmx_domdec_t*              dd,
                             const interaction_const_t& interactionConst,
                             bool                       bFreeEnergy_q,
                             bool                       bFreeEnergy_lj,
                             gmx::ArrayRef<const real>  chargeA,
                             gmx::ArrayRef<const real>  chargeB,
                             gmx::ArrayRef<const real>  sqrt_c6A,
                             gmx::ArrayRef<const real>  sqrt_c6B,
                             gmx::ArrayRef<const real>  sigmaA,
                             gmx::ArrayRef<const real>  sigmaB,
                             int                        maxshift_x,
                             int                        maxshift_y)
{
    unsigned int flags = 0;

    if (usingPme(interactionConst.coulomb.type))
    {
        flags |= PP_PME_CHARGE;
    }
    if (usingLJPme(interactionConst.vdw.type))
    {
        flags |= (PP_PME_SQRTC6 | PP_PME_SIGMA);
    }
    if (bFreeEnergy_q || bFreeEnergy_lj)
    {
        /* Assumes that the B state flags are in the bits just above
         * the ones for the A state. */
        flags |= (flags << 1);
    }

    gmx_pme_send_coeffs_coords(nullptr,
                               mpiCommMySim,
                               dd,
                               flags,
                               chargeA,
                               chargeB,
                               sqrt_c6A,
                               sqrt_c6B,
                               sigmaA,
                               sigmaB,
                               nullptr,
                               gmx::ArrayRef<gmx::RVec>(),
                               0,
                               0,
                               maxshift_x,
                               maxshift_y,
                               -1,
                               false,
                               false,
                               false,
                               false,
                               false,
                               nullptr);
}

void gmx_pme_send_coordinates(t_forcerec*                    fr,
                              const gmx::MpiComm&            mpiCommMySim,
                              gmx_domdec_t*                  dd,
                              const matrix                   box,
                              gmx::ArrayRef<const gmx::RVec> x,
                              real                           lambda_q,
                              real                           lambda_lj,
                              bool                           computeEnergyAndVirial,
                              int64_t                        step,
                              bool                           useGpuPmePpComms,
                              bool                           receiveCoordinateAddressFromPme,
                              bool                           sendCoordinatesFromGpu,
                              bool                           receiveForcesToGpu,
                              GpuEventSynchronizer*          coordinatesReadyOnDeviceEvent,
                              bool                           useMdGpuGraph,
                              gmx_wallcycle*                 wcycle)
{
    wallcycle_start(wcycle, WallCycleCounter::PpPmeSendX);

    unsigned int flags = PP_PME_COORD;
    if (computeEnergyAndVirial)
    {
        flags |= PP_PME_ENER_VIR;
    }
    gmx_pme_send_coeffs_coords(fr,
                               mpiCommMySim,
                               dd,
                               flags,
                               {},
                               {},
                               {},
                               {},
                               {},
                               {},
                               box,
                               x,
                               lambda_q,
                               lambda_lj,
                               0,
                               0,
                               step,
                               useGpuPmePpComms,
                               receiveCoordinateAddressFromPme,
                               sendCoordinatesFromGpu,
                               receiveForcesToGpu,
                               useMdGpuGraph,
                               coordinatesReadyOnDeviceEvent);

    wallcycle_stop(wcycle, WallCycleCounter::PpPmeSendX);
}

void gmx_pme_send_finish(const gmx::MpiComm& mpiCommMySim, gmx_domdec_t* dd)
{
    unsigned int flags = PP_PME_FINISH;

    gmx_pme_send_coeffs_coords(nullptr,
                               mpiCommMySim,
                               dd,
                               flags,
                               {},
                               {},
                               {},
                               {},
                               {},
                               {},
                               nullptr,
                               gmx::ArrayRef<gmx::RVec>(),
                               0,
                               0,
                               0,
                               0,
                               -1,
                               false,
                               false,
                               false,
                               false,
                               false,
                               nullptr);
}

void gmx_pme_send_switchgrid(const gmx::MpiComm& mpiCommMySim,
                             const gmx_domdec_t& dd,
                             ivec                grid_size,
                             real                ewaldcoeff_q,
                             real                ewaldcoeff_lj)
{
#if GMX_MPI
    gmx_pme_comm_n_box_t cnb;

    /* Only let one PP node signal each PME node */
    if (dd.pme_receive_vir_ener)
    {
        cnb.flags = PP_PME_SWITCHGRID;
        copy_ivec(grid_size, cnb.grid_size);
        cnb.ewaldcoeff_q  = ewaldcoeff_q;
        cnb.ewaldcoeff_lj = ewaldcoeff_lj;

        /* We send this, uncommon, message blocking to simplify the code */
        MPI_Send(&cnb, sizeof(cnb), MPI_BYTE, dd.pme_nodeid, eCommType_CNB, mpiCommMySim.comm());
    }
#else
    GMX_UNUSED_VALUE(mpiCommMySim);
    GMX_UNUSED_VALUE(dd);
    GMX_UNUSED_VALUE(grid_size);
    GMX_UNUSED_VALUE(ewaldcoeff_q);
    GMX_UNUSED_VALUE(ewaldcoeff_lj);
#endif
}

void gmx_pme_send_resetcounters(const gmx::MpiComm& mpiCommMySim, gmx_domdec_t* dd, int64_t gmx_unused step)
{
#if GMX_MPI
    gmx_pme_comm_n_box_t cnb;

    /* Only let one PP node signal each PME node */
    if (dd->pme_receive_vir_ener)
    {
        cnb.flags = PP_PME_RESETCOUNTERS;
        cnb.step  = step;

        /* We send this, uncommon, message blocking to simplify the code */
        MPI_Send(&cnb, sizeof(cnb), MPI_BYTE, dd->pme_nodeid, eCommType_CNB, mpiCommMySim.comm());
    }
#else
    GMX_UNUSED_VALUE(mpiCommMySim);
    GMX_UNUSED_VALUE(dd);
    GMX_UNUSED_VALUE(step);
#endif
}

/*! \brief Receive virial and energy from PME rank */
static void receive_virial_energy(const gmx::MpiComm&   mpiCommMySim,
                                  const gmx_domdec_t*   dd,
                                  gmx::ForceWithVirial* forceWithVirial,
                                  real*                 energy_q,
                                  real*                 energy_lj,
                                  real*                 dvdlambda_q,
                                  real*                 dvdlambda_lj,
                                  float*                pme_cycles)
{
    gmx_pme_comm_vir_ene_t cve;

    if (dd->pme_receive_vir_ener)
    {
        if (debug)
        {
            fprintf(debug,
                    "PP rank %d receiving from PME rank %d: virial and energy\n",
                    mpiCommMySim.rank(),
                    dd->pme_nodeid);
        }
#if GMX_MPI
        MPI_Recv(&cve,
                 sizeof(cve),
                 MPI_BYTE,
                 dd->pme_nodeid,
                 eCommType_ENERGY_VIRIAL_DVDL,
                 mpiCommMySim.comm(),
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

/*! \brief Recieve force data from PME ranks */
static void recvFFromPme(gmx::PmePpCommGpu*  pmePpCommGpu,
                         void*               recvptr,
                         int                 n,
                         const gmx::MpiComm& mpiCommMySim,
                         const gmx_domdec_t& dd,
                         bool                useGpuPmePpComms,
                         bool                receivePmeForceToGpu)
{
    if (useGpuPmePpComms)
    {
        GMX_ASSERT(pmePpCommGpu != nullptr, "Need valid pmePpCommGpu");
        // Receive forces from PME rank
        pmePpCommGpu->receiveForceFromPme(static_cast<gmx::RVec*>(recvptr), n, receivePmeForceToGpu);
    }
    else
    {
        // Receive data using MPI
#if GMX_MPI
        MPI_Recv(recvptr, n * sizeof(rvec), MPI_BYTE, dd.pme_nodeid, eCommType_FORCES, mpiCommMySim.comm(), MPI_STATUS_IGNORE);
#else
        GMX_UNUSED_VALUE(mpiCommMySim);
        GMX_UNUSED_VALUE(dd);
        GMX_UNUSED_VALUE(n);
#endif
    }
}


void gmx_pme_receive_f(gmx::PmePpCommGpu*    pmePpCommGpu,
                       const gmx::MpiComm&   mpiCommMySim,
                       gmx_domdec_t*         dd,
                       gmx::ForceWithVirial* forceWithVirial,
                       real*                 energy_q,
                       real*                 energy_lj,
                       real*                 dvdlambda_q,
                       real*                 dvdlambda_lj,
                       bool                  useGpuPmePpComms,
                       bool                  receivePmeForceToGpu,
                       float*                pme_cycles)
{
    if (c_useDelayedWait)
    {
        /* Wait for the x request to finish */
        gmx_pme_send_coeffs_coords_wait(dd);
    }

    const int                   natoms = dd_numHomeAtoms(*dd);
    gmx::HostVector<gmx::RVec>& buffer = dd->pmeForceReceiveBuffer;
    buffer.resize(natoms);

    void* recvptr = reinterpret_cast<void*>(buffer.data());
    recvFFromPme(pmePpCommGpu, recvptr, natoms, mpiCommMySim, *dd, useGpuPmePpComms, receivePmeForceToGpu);

    int nt = gmx_omp_nthreads_get_simple_rvec_task(ModuleMultiThread::Default, natoms);

    gmx::ArrayRef<gmx::RVec> f = forceWithVirial->force_;

    if (!receivePmeForceToGpu)
    {
        /* Note that we would like to avoid this conditional by putting it
         * into the omp pragma instead, but then we still take the full
         * omp parallel for overhead (at least with gcc5).
         */
        if (nt == 1)
        {
            for (int i = 0; i < natoms; i++)
            {
                f[i] += buffer[i];
            }
        }
        else
        {
#pragma omp parallel for num_threads(nt) schedule(static)
            for (int i = 0; i < natoms; i++)
            {
                f[i] += buffer[i];
            }
        }
    }

    receive_virial_energy(
            mpiCommMySim, dd, forceWithVirial, energy_q, energy_lj, dvdlambda_q, dvdlambda_lj, pme_cycles);
}
