/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
/* IMPORTANT FOR DEVELOPERS:
 *
 * Triclinic pme stuff isn't entirely trivial, and we've experienced
 * some bugs during development (many of them due to me). To avoid
 * this in the future, please check the following things if you make
 * changes in this file:
 *
 * 1. You should obtain identical (at least to the PME precision)
 *    energies, forces, and virial for
 *    a rectangular box and a triclinic one where the z (or y) axis is
 *    tilted a whole box side. For instance you could use these boxes:
 *
 *    rectangular       triclinic
 *     2  0  0           2  0  0
 *     0  2  0           0  2  0
 *     0  0  6           2  2  6
 *
 * 2. You should check the energy conservation in a triclinic box.
 *
 * It might seem an overkill, but better safe than sorry.
 * /Erik 001109
 */

#include "gmxpre.h"

#include "pme_only.h"

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <numeric>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_coordinate_receiver_gpu.h"
#include "gromacs/ewald/pme_force_sender_gpu.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

#include "pme_gpu_internal.h"
#include "pme_output.h"
#include "pme_pp_communication.h"

/*! \brief environment variable to enable GPU P2P communication */
static const bool c_enableGpuPmePpComms =
        (getenv("GMX_GPU_PME_PP_COMMS") != nullptr) && GMX_THREAD_MPI && (GMX_GPU == GMX_GPU_CUDA);

/*! \brief Master PP-PME communication data structure */
struct gmx_pme_pp
{
    MPI_Comm             mpi_comm_mysim; /**< MPI communicator for this simulation */
    std::vector<PpRanks> ppRanks;        /**< The PP partner ranks                 */
    int                  peerRankId;     /**< The peer PP rank id                  */
    //@{
    /**< Vectors of A- and B-state parameters used to transfer vectors to PME ranks  */
    gmx::PaddedHostVector<real> chargeA;
    std::vector<real>           chargeB;
    std::vector<real>           sqrt_c6A;
    std::vector<real>           sqrt_c6B;
    std::vector<real>           sigmaA;
    std::vector<real>           sigmaB;
    //@}
    gmx::HostVector<gmx::RVec> x; /**< Vector of atom coordinates to transfer to PME ranks */
    std::vector<gmx::RVec>     f; /**< Vector of atom forces received from PME ranks */
    //@{
    /**< Vectors of MPI objects used in non-blocking communication between multiple PP ranks per PME rank */
    std::vector<MPI_Request> req;
    std::vector<MPI_Status>  stat;
    //@}

    /*! \brief object for receiving coordinates using communications operating on GPU memory space */
    std::unique_ptr<gmx::PmeCoordinateReceiverGpu> pmeCoordinateReceiverGpu;
    /*! \brief object for sending PME force using communications operating on GPU memory space */
    std::unique_ptr<gmx::PmeForceSenderGpu> pmeForceSenderGpu;

    /*! \brief whether GPU direct communications are active for PME-PP transfers */
    bool useGpuDirectComm = false;
};

/*! \brief Initialize the PME-only side of the PME <-> PP communication */
static std::unique_ptr<gmx_pme_pp> gmx_pme_pp_init(const t_commrec* cr)
{
    auto pme_pp = std::make_unique<gmx_pme_pp>();

#if GMX_MPI
    int rank;

    pme_pp->mpi_comm_mysim = cr->mpi_comm_mysim;
    MPI_Comm_rank(cr->mpi_comm_mygroup, &rank);
    auto ppRanks = get_pme_ddranks(cr, rank);
    pme_pp->ppRanks.reserve(ppRanks.size());
    for (const auto& ppRankId : ppRanks)
    {
        pme_pp->ppRanks.push_back({ ppRankId, 0 });
    }
    // The peer PP rank is the last one.
    pme_pp->peerRankId = pme_pp->ppRanks.back().rankId;
    pme_pp->req.resize(eCommType_NR * pme_pp->ppRanks.size());
    pme_pp->stat.resize(eCommType_NR * pme_pp->ppRanks.size());
#else
    GMX_UNUSED_VALUE(cr);
#endif

    return pme_pp;
}

static void reset_pmeonly_counters(gmx_wallcycle_t           wcycle,
                                   gmx_walltime_accounting_t walltime_accounting,
                                   t_nrnb*                   nrnb,
                                   int64_t                   step,
                                   bool                      useGpuForPme)
{
    /* Reset all the counters related to performance over the run */
    wallcycle_stop(wcycle, ewcRUN);
    wallcycle_reset_all(wcycle);
    *nrnb = { 0 };
    wallcycle_start(wcycle, ewcRUN);
    walltime_accounting_reset_time(walltime_accounting, step);

    if (useGpuForPme)
    {
        resetGpuProfiler();
    }
}

static gmx_pme_t* gmx_pmeonly_switch(std::vector<gmx_pme_t*>* pmedata,
                                     const ivec               grid_size,
                                     real                     ewaldcoeff_q,
                                     real                     ewaldcoeff_lj,
                                     const t_commrec*         cr,
                                     const t_inputrec*        ir)
{
    GMX_ASSERT(pmedata, "Bad PME tuning list pointer");
    for (auto& pme : *pmedata)
    {
        GMX_ASSERT(pme, "Bad PME tuning list element pointer");
        if (gmx_pme_grid_matches(*pme, grid_size))
        {
            /* Here we have found an existing PME data structure that suits us.
             * However, in the GPU case, we have to reinitialize it - there's only one GPU structure.
             * This should not cause actual GPU reallocations, at least (the allocated buffers are never shrunk).
             * So, just some grid size updates in the GPU kernel parameters.
             * TODO: this should be something like gmx_pme_update_split_params()
             */
            gmx_pme_reinit(&pme, cr, pme, ir, grid_size, ewaldcoeff_q, ewaldcoeff_lj);
            return pme;
        }
    }

    const auto& pme          = pmedata->back();
    gmx_pme_t*  newStructure = nullptr;
    // Copy last structure with new grid params
    gmx_pme_reinit(&newStructure, cr, pme, ir, grid_size, ewaldcoeff_q, ewaldcoeff_lj);
    pmedata->push_back(newStructure);
    return newStructure;
}

/*! \brief Called by PME-only ranks to receive coefficients and coordinates
 *
 * \param[in] pme                     PME data structure.
 * \param[in,out] pme_pp              PME-PP communication structure.
 * \param[out] natoms                 Number of received atoms.
 * \param[out] box                    System box, if received.
 * \param[out] maxshift_x             Maximum shift in X direction, if received.
 * \param[out] maxshift_y             Maximum shift in Y direction, if received.
 * \param[out] lambda_q               Free-energy lambda for electrostatics, if received.
 * \param[out] lambda_lj              Free-energy lambda for Lennard-Jones, if received.
 * \param[out] computeEnergyAndVirial Set to true if this is an energy/virial calculation
 *                                    step, otherwise set to false.
 * \param[out] step                   MD integration step number.
 * \param[out] grid_size              PME grid size, if received.
 * \param[out] ewaldcoeff_q           Ewald cut-off parameter for electrostatics, if received.
 * \param[out] ewaldcoeff_lj          Ewald cut-off parameter for Lennard-Jones, if received.
 * \param[in]  useGpuForPme           Flag on whether PME is on GPU.
 * \param[in]  stateGpu               GPU state propagator object.
 * \param[in]  runMode                PME run mode.
 *
 * \retval pmerecvqxX                 All parameters were set, chargeA and chargeB can be NULL.
 * \retval pmerecvqxFINISH            No parameters were set.
 * \retval pmerecvqxSWITCHGRID        Only grid_size and *ewaldcoeff were set.
 * \retval pmerecvqxRESETCOUNTERS     *step was set.
 */
static int gmx_pme_recv_coeffs_coords(struct gmx_pme_t*            pme,
                                      gmx_pme_pp*                  pme_pp,
                                      int*                         natoms,
                                      matrix                       box,
                                      int*                         maxshift_x,
                                      int*                         maxshift_y,
                                      real*                        lambda_q,
                                      real*                        lambda_lj,
                                      gmx_bool*                    computeEnergyAndVirial,
                                      int64_t*                     step,
                                      ivec*                        grid_size,
                                      real*                        ewaldcoeff_q,
                                      real*                        ewaldcoeff_lj,
                                      bool                         useGpuForPme,
                                      gmx::StatePropagatorDataGpu* stateGpu,
                                      PmeRunMode gmx_unused runMode)
{
    int status = -1;
    int nat    = 0;

#if GMX_MPI
    unsigned int flags          = 0;
    int          messages       = 0;
    bool         atomSetChanged = false;

    do
    {
        gmx_pme_comm_n_box_t cnb;
        cnb.flags = 0;

        /* Receive the send count, box and time step from the peer PP node */
        MPI_Recv(&cnb, sizeof(cnb), MPI_BYTE, pme_pp->peerRankId, eCommType_CNB,
                 pme_pp->mpi_comm_mysim, MPI_STATUS_IGNORE);

        /* We accumulate all received flags */
        flags |= cnb.flags;

        *step = cnb.step;

        if (debug)
        {
            fprintf(debug, "PME only rank receiving:%s%s%s%s%s\n",
                    (cnb.flags & PP_PME_CHARGE) ? " charges" : "",
                    (cnb.flags & PP_PME_COORD) ? " coordinates" : "",
                    (cnb.flags & PP_PME_FINISH) ? " finish" : "",
                    (cnb.flags & PP_PME_SWITCHGRID) ? " switch grid" : "",
                    (cnb.flags & PP_PME_RESETCOUNTERS) ? " reset counters" : "");
        }

        pme_pp->useGpuDirectComm = ((cnb.flags & PP_PME_GPUCOMMS) != 0);
        GMX_ASSERT(!pme_pp->useGpuDirectComm || (pme_pp->pmeForceSenderGpu != nullptr),
                   "The use of GPU direct communication for PME-PP is enabled, "
                   "but the PME GPU force reciever object does not exist");

        if (cnb.flags & PP_PME_FINISH)
        {
            status = pmerecvqxFINISH;
        }

        if (cnb.flags & PP_PME_SWITCHGRID)
        {
            /* Special case, receive the new parameters and return */
            copy_ivec(cnb.grid_size, *grid_size);
            *ewaldcoeff_q  = cnb.ewaldcoeff_q;
            *ewaldcoeff_lj = cnb.ewaldcoeff_lj;

            status = pmerecvqxSWITCHGRID;
        }

        if (cnb.flags & PP_PME_RESETCOUNTERS)
        {
            /* Special case, receive the step (set above) and return */
            status = pmerecvqxRESETCOUNTERS;
        }

        if (cnb.flags & (PP_PME_CHARGE | PP_PME_SQRTC6 | PP_PME_SIGMA))
        {
            atomSetChanged = true;

            /* Receive the send counts from the other PP nodes */
            for (auto& sender : pme_pp->ppRanks)
            {
                if (sender.rankId == pme_pp->peerRankId)
                {
                    sender.numAtoms = cnb.natoms;
                }
                else
                {
                    MPI_Irecv(&sender.numAtoms, sizeof(sender.numAtoms), MPI_BYTE, sender.rankId,
                              eCommType_CNB, pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);
                }
            }
            MPI_Waitall(messages, pme_pp->req.data(), pme_pp->stat.data());
            messages = 0;

            nat = 0;
            for (const auto& sender : pme_pp->ppRanks)
            {
                nat += sender.numAtoms;
            }

            if (cnb.flags & PP_PME_CHARGE)
            {
                pme_pp->chargeA.resizeWithPadding(nat);
            }
            if (cnb.flags & PP_PME_CHARGEB)
            {
                pme_pp->chargeB.resize(nat);
            }
            if (cnb.flags & PP_PME_SQRTC6)
            {
                pme_pp->sqrt_c6A.resize(nat);
            }
            if (cnb.flags & PP_PME_SQRTC6B)
            {
                pme_pp->sqrt_c6B.resize(nat);
            }
            if (cnb.flags & PP_PME_SIGMA)
            {
                pme_pp->sigmaA.resize(nat);
            }
            if (cnb.flags & PP_PME_SIGMAB)
            {
                pme_pp->sigmaB.resize(nat);
            }
            pme_pp->x.resize(nat);
            pme_pp->f.resize(nat);

            /* maxshift is sent when the charges are sent */
            *maxshift_x = cnb.maxshift_x;
            *maxshift_y = cnb.maxshift_y;

            /* Receive the charges in place */
            for (int q = 0; q < eCommType_NR; q++)
            {
                real* bufferPtr;

                if (!(cnb.flags & (PP_PME_CHARGE << q)))
                {
                    continue;
                }
                switch (q)
                {
                    case eCommType_ChargeA: bufferPtr = pme_pp->chargeA.data(); break;
                    case eCommType_ChargeB: bufferPtr = pme_pp->chargeB.data(); break;
                    case eCommType_SQRTC6A: bufferPtr = pme_pp->sqrt_c6A.data(); break;
                    case eCommType_SQRTC6B: bufferPtr = pme_pp->sqrt_c6B.data(); break;
                    case eCommType_SigmaA: bufferPtr = pme_pp->sigmaA.data(); break;
                    case eCommType_SigmaB: bufferPtr = pme_pp->sigmaB.data(); break;
                    default: gmx_incons("Wrong eCommType");
                }
                nat = 0;
                for (const auto& sender : pme_pp->ppRanks)
                {
                    if (sender.numAtoms > 0)
                    {
                        MPI_Irecv(bufferPtr + nat, sender.numAtoms * sizeof(real), MPI_BYTE,
                                  sender.rankId, q, pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);
                        nat += sender.numAtoms;
                        if (debug)
                        {
                            fprintf(debug, "Received from PP rank %d: %d %s\n", sender.rankId,
                                    sender.numAtoms,
                                    (q == eCommType_ChargeA || q == eCommType_ChargeB) ? "charges"
                                                                                       : "params");
                        }
                    }
                }
            }
        }

        if (cnb.flags & PP_PME_COORD)
        {
            if (atomSetChanged)
            {
                gmx_pme_reinit_atoms(pme, nat, pme_pp->chargeA.data());
                if (useGpuForPme)
                {
                    stateGpu->reinit(nat, nat);
                    pme_gpu_set_device_x(pme, stateGpu->getCoordinates());
                }
                if (pme_pp->useGpuDirectComm)
                {
                    GMX_ASSERT(runMode == PmeRunMode::GPU,
                               "GPU Direct PME-PP communication has been enabled, "
                               "but PME run mode is not PmeRunMode::GPU\n");

                    // This rank will have its data accessed directly by PP rank, so needs to send the remote addresses.
                    pme_pp->pmeCoordinateReceiverGpu->sendCoordinateBufferAddressToPpRanks(
                            stateGpu->getCoordinates());
                    pme_pp->pmeForceSenderGpu->sendForceBufferAddressToPpRanks(
                            reinterpret_cast<rvec*>(pme_gpu_get_device_f(pme)));
                }
            }


            /* The box, FE flag and lambda are sent along with the coordinates
             *  */
            copy_mat(cnb.box, box);
            *lambda_q               = cnb.lambda_q;
            *lambda_lj              = cnb.lambda_lj;
            *computeEnergyAndVirial = ((cnb.flags & PP_PME_ENER_VIR) != 0U);
            *step                   = cnb.step;

            /* Receive the coordinates in place */
            nat = 0;
            for (const auto& sender : pme_pp->ppRanks)
            {
                if (sender.numAtoms > 0)
                {
                    if (pme_pp->useGpuDirectComm)
                    {
                        pme_pp->pmeCoordinateReceiverGpu->launchReceiveCoordinatesFromPpCudaDirect(
                                sender.rankId);
                    }
                    else
                    {
                        MPI_Irecv(pme_pp->x[nat], sender.numAtoms * sizeof(rvec), MPI_BYTE, sender.rankId,
                                  eCommType_COORD, pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);
                    }
                    nat += sender.numAtoms;
                    if (debug)
                    {
                        fprintf(debug,
                                "Received from PP rank %d: %d "
                                "coordinates\n",
                                sender.rankId, sender.numAtoms);
                    }
                }
            }

            if (pme_pp->useGpuDirectComm)
            {
                pme_pp->pmeCoordinateReceiverGpu->enqueueWaitReceiveCoordinatesFromPpCudaDirect();
            }

            status = pmerecvqxX;
        }

        /* Wait for the coordinates and/or charges to arrive */
        MPI_Waitall(messages, pme_pp->req.data(), pme_pp->stat.data());
        messages = 0;
    } while (status == -1);
#else
    GMX_UNUSED_VALUE(pme);
    GMX_UNUSED_VALUE(pme_pp);
    GMX_UNUSED_VALUE(box);
    GMX_UNUSED_VALUE(maxshift_x);
    GMX_UNUSED_VALUE(maxshift_y);
    GMX_UNUSED_VALUE(lambda_q);
    GMX_UNUSED_VALUE(lambda_lj);
    GMX_UNUSED_VALUE(computeEnergyAndVirial);
    GMX_UNUSED_VALUE(step);
    GMX_UNUSED_VALUE(grid_size);
    GMX_UNUSED_VALUE(ewaldcoeff_q);
    GMX_UNUSED_VALUE(ewaldcoeff_lj);
    GMX_UNUSED_VALUE(useGpuForPme);
    GMX_UNUSED_VALUE(stateGpu);

    status = pmerecvqxX;
#endif

    if (status == pmerecvqxX)
    {
        *natoms = nat;
    }

    return status;
}

#if GMX_MPI
/*! \brief Send force data to PP ranks */
static void sendFToPP(void* sendbuf, PpRanks receiver, gmx_pme_pp* pme_pp, int* messages)
{

    if (pme_pp->useGpuDirectComm)
    {
        GMX_ASSERT((pme_pp->pmeForceSenderGpu != nullptr),
                   "The use of GPU direct communication for PME-PP is enabled, "
                   "but the PME GPU force reciever object does not exist");

        pme_pp->pmeForceSenderGpu->sendFToPpCudaDirect(receiver.rankId);
    }
    else
    {
        // Send using MPI
        MPI_Isend(sendbuf, receiver.numAtoms * sizeof(rvec), MPI_BYTE, receiver.rankId, 0,
                  pme_pp->mpi_comm_mysim, &pme_pp->req[*messages]);
        *messages = *messages + 1;
    }
}
#endif

/*! \brief Send the PME mesh force, virial and energy to the PP-only ranks. */
static void gmx_pme_send_force_vir_ener(const gmx_pme_t& pme,
                                        gmx_pme_pp*      pme_pp,
                                        const PmeOutput& output,
                                        real             dvdlambda_q,
                                        real             dvdlambda_lj,
                                        float            cycles)
{
#if GMX_MPI
    gmx_pme_comm_vir_ene_t cve;
    int                    messages, ind_start, ind_end;
    cve.cycles = cycles;

    /* Now the evaluated forces have to be transferred to the PP nodes */
    messages = 0;
    ind_end  = 0;
    for (const auto& receiver : pme_pp->ppRanks)
    {
        ind_start     = ind_end;
        ind_end       = ind_start + receiver.numAtoms;
        void* sendbuf = const_cast<void*>(static_cast<const void*>(output.forces_[ind_start]));
        if (pme_pp->useGpuDirectComm)
        {
            // Data will be transferred directly from GPU.
            rvec* d_f = reinterpret_cast<rvec*>(pme_gpu_get_device_f(&pme));
            sendbuf   = reinterpret_cast<void*>(&d_f[ind_start]);
        }
        sendFToPP(sendbuf, receiver, pme_pp, &messages);
    }

    /* send virial and energy to our last PP node */
    copy_mat(output.coulombVirial_, cve.vir_q);
    copy_mat(output.lennardJonesVirial_, cve.vir_lj);
    cve.energy_q     = output.coulombEnergy_;
    cve.energy_lj    = output.lennardJonesEnergy_;
    cve.dvdlambda_q  = dvdlambda_q;
    cve.dvdlambda_lj = dvdlambda_lj;
    /* check for the signals to send back to a PP node */
    cve.stop_cond = gmx_get_stop_condition();

    cve.cycles = cycles;

    if (debug)
    {
        fprintf(debug, "PME rank sending to PP rank %d: virial and energy\n", pme_pp->peerRankId);
    }
    MPI_Isend(&cve, sizeof(cve), MPI_BYTE, pme_pp->peerRankId, 1, pme_pp->mpi_comm_mysim,
              &pme_pp->req[messages++]);

    /* Wait for the forces to arrive */
    MPI_Waitall(messages, pme_pp->req.data(), pme_pp->stat.data());
#else
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_pme_send_force_vir_ener");
    GMX_UNUSED_VALUE(pme);
    GMX_UNUSED_VALUE(pme_pp);
    GMX_UNUSED_VALUE(output);
    GMX_UNUSED_VALUE(dvdlambda_q);
    GMX_UNUSED_VALUE(dvdlambda_lj);
    GMX_UNUSED_VALUE(cycles);
#endif
}

int gmx_pmeonly(struct gmx_pme_t*               pme,
                const t_commrec*                cr,
                t_nrnb*                         mynrnb,
                gmx_wallcycle*                  wcycle,
                gmx_walltime_accounting_t       walltime_accounting,
                t_inputrec*                     ir,
                PmeRunMode                      runMode,
                const gmx::DeviceStreamManager* deviceStreamManager)
{
    int     ret;
    int     natoms = 0;
    matrix  box;
    real    lambda_q   = 0;
    real    lambda_lj  = 0;
    int     maxshift_x = 0, maxshift_y = 0;
    real    dvdlambda_q, dvdlambda_lj;
    float   cycles;
    int     count;
    bool    computeEnergyAndVirial = false;
    int64_t step;

    /* This data will only use with PME tuning, i.e. switching PME grids */
    std::vector<gmx_pme_t*> pmedata;
    pmedata.push_back(pme);

    auto pme_pp = gmx_pme_pp_init(cr);

    std::unique_ptr<gmx::StatePropagatorDataGpu> stateGpu;
    // TODO the variable below should be queried from the task assignment info
    const bool useGpuForPme = (runMode == PmeRunMode::GPU) || (runMode == PmeRunMode::Mixed);
    if (useGpuForPme)
    {
        GMX_RELEASE_ASSERT(
                deviceStreamManager != nullptr,
                "Device stream manager can not be nullptr when using GPU in PME-only rank.");
        GMX_RELEASE_ASSERT(deviceStreamManager->streamIsValid(gmx::DeviceStreamType::Pme),
                           "Device stream can not be nullptr when using GPU in PME-only rank");
        changePinningPolicy(&pme_pp->chargeA, pme_get_pinning_policy());
        changePinningPolicy(&pme_pp->x, pme_get_pinning_policy());
        if (c_enableGpuPmePpComms)
        {
            pme_pp->pmeCoordinateReceiverGpu = std::make_unique<gmx::PmeCoordinateReceiverGpu>(
                    deviceStreamManager->stream(gmx::DeviceStreamType::Pme), pme_pp->mpi_comm_mysim,
                    pme_pp->ppRanks);
            pme_pp->pmeForceSenderGpu = std::make_unique<gmx::PmeForceSenderGpu>(
                    deviceStreamManager->stream(gmx::DeviceStreamType::Pme), pme_pp->mpi_comm_mysim,
                    pme_pp->ppRanks);
        }
        // TODO: Special PME-only constructor is used here. There is no mechanism to prevent from using the other constructor here.
        //       This should be made safer.
        stateGpu = std::make_unique<gmx::StatePropagatorDataGpu>(
                &deviceStreamManager->stream(gmx::DeviceStreamType::Pme), deviceStreamManager->context(),
                GpuApiCallBehavior::Async, pme_gpu_get_block_size(pme), wcycle);
    }

    clear_nrnb(mynrnb);

    count = 0;
    do /****** this is a quasi-loop over time steps! */
    {
        /* The reason for having a loop here is PME grid tuning/switching */
        do
        {
            /* Domain decomposition */
            ivec newGridSize;
            real ewaldcoeff_q = 0, ewaldcoeff_lj = 0;
            ret = gmx_pme_recv_coeffs_coords(pme, pme_pp.get(), &natoms, box, &maxshift_x, &maxshift_y,
                                             &lambda_q, &lambda_lj, &computeEnergyAndVirial, &step,
                                             &newGridSize, &ewaldcoeff_q, &ewaldcoeff_lj,
                                             useGpuForPme, stateGpu.get(), runMode);

            if (ret == pmerecvqxSWITCHGRID)
            {
                /* Switch the PME grid to newGridSize */
                pme = gmx_pmeonly_switch(&pmedata, newGridSize, ewaldcoeff_q, ewaldcoeff_lj, cr, ir);
            }

            if (ret == pmerecvqxRESETCOUNTERS)
            {
                /* Reset the cycle and flop counters */
                reset_pmeonly_counters(wcycle, walltime_accounting, mynrnb, step, useGpuForPme);
            }
        } while (ret == pmerecvqxSWITCHGRID || ret == pmerecvqxRESETCOUNTERS);

        if (ret == pmerecvqxFINISH)
        {
            /* We should stop: break out of the loop */
            break;
        }

        if (count == 0)
        {
            wallcycle_start(wcycle, ewcRUN);
            walltime_accounting_start_time(walltime_accounting);
        }

        wallcycle_start(wcycle, ewcPMEMESH);

        dvdlambda_q  = 0;
        dvdlambda_lj = 0;

        // TODO Make a struct of array refs onto these per-atom fields
        // of pme_pp (maybe box, energy and virial, too; and likewise
        // from mdatoms for the other call to gmx_pme_do), so we have
        // fewer lines of code and less parameter passing.
        gmx::StepWorkload stepWork;
        stepWork.computeVirial = computeEnergyAndVirial;
        stepWork.computeEnergy = computeEnergyAndVirial;
        stepWork.computeForces = true;
        PmeOutput output;
        if (useGpuForPme)
        {
            stepWork.haveDynamicBox      = false;
            stepWork.useGpuPmeFReduction = pme_pp->useGpuDirectComm;
            // TODO this should be set properly by gmx_pme_recv_coeffs_coords,
            // or maybe use inputrecDynamicBox(ir), at the very least - change this when this codepath is tested!
            pme_gpu_prepare_computation(pme, box, wcycle, stepWork);
            if (!pme_pp->useGpuDirectComm)
            {
                stateGpu->copyCoordinatesToGpu(gmx::ArrayRef<gmx::RVec>(pme_pp->x), gmx::AtomLocality::All);
            }
            // On the separate PME rank we do not need a synchronizer as we schedule everything in a single stream
            // TODO: with pme on GPU the receive should make a list of synchronizers and pass it here #3157
            auto xReadyOnDevice = nullptr;

            pme_gpu_launch_spread(pme, xReadyOnDevice, wcycle);
            pme_gpu_launch_complex_transforms(pme, wcycle, stepWork);
            pme_gpu_launch_gather(pme, wcycle);
            output = pme_gpu_wait_finish_task(pme, computeEnergyAndVirial, wcycle);
            pme_gpu_reinit_computation(pme, wcycle);
        }
        else
        {
            GMX_ASSERT(pme_pp->x.size() == static_cast<size_t>(natoms),
                       "The coordinate buffer should have size natoms");

            gmx_pme_do(pme, pme_pp->x, pme_pp->f, pme_pp->chargeA.data(), pme_pp->chargeB.data(),
                       pme_pp->sqrt_c6A.data(), pme_pp->sqrt_c6B.data(), pme_pp->sigmaA.data(),
                       pme_pp->sigmaB.data(), box, cr, maxshift_x, maxshift_y, mynrnb, wcycle,
                       output.coulombVirial_, output.lennardJonesVirial_, &output.coulombEnergy_,
                       &output.lennardJonesEnergy_, lambda_q, lambda_lj, &dvdlambda_q,
                       &dvdlambda_lj, stepWork);
            output.forces_ = pme_pp->f;
        }

        cycles = wallcycle_stop(wcycle, ewcPMEMESH);
        gmx_pme_send_force_vir_ener(*pme, pme_pp.get(), output, dvdlambda_q, dvdlambda_lj, cycles);

        count++;
    } /***** end of quasi-loop, we stop with the break above */
    while (TRUE);

    walltime_accounting_end_time(walltime_accounting);

    return 0;
}
