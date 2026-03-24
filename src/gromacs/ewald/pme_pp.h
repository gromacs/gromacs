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
 *
 * \brief This file contains function declarations necessary for
 * mananging the PP side of PME-only ranks.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_PP_H
#define GMX_EWALD_PME_PP_H

#include <cstdint>

#include <optional>
#include <vector>

#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

#include "pme_pp_communication.h"

struct gmx_wallcycle;

class DeviceContext;
class DeviceStream;
class GpuEventSynchronizer;

namespace gmx
{
class DeviceStreamManager;
class ForceWithVirial;
template<typename>
class ArrayRef;

/*! \brief Coordinate communication from PP rank with its PME-only rank
 *
 * When domains are repartitioned, call \c sendParameters() to inform the
 * PME-only rank how this domain will be handled. Every MD step, call \c
 * sendCoordinates() and later \c receiveResults().
 */
class PmePpComm
{
public:
    //! Constructor
    PmePpComm(const MpiComm&             comm,
              int                        rankOfPartnerPmeRank,
              VanDerWaalsType            vdwType,
              CoulombInteractionType     coulombType,
              bool                       thisRankReceivesVirialAndEnergy,
              bool                       useGpuPmePpCommunication,
              bool                       useNvshmem,
              const DeviceStreamManager* deviceStreamManager);
    /*! \brief Send the atomic parameters, settings, and maxshift to
     * the PME-only rank
     *
     * Must be called at the beginning of the lifetime of a
     * domain repartitioning. */
    void sendParameters(int                  numHomeAtoms,
                        bool                 bFreeEnergy_q,
                        bool                 bFreeEnergy_lj,
                        ArrayRef<const real> chargeA,
                        ArrayRef<const real> chargeB,
                        ArrayRef<const real> sqrt_c6A,
                        ArrayRef<const real> sqrt_c6B,
                        ArrayRef<const real> sigmaA,
                        ArrayRef<const real> sigmaB,
                        int                  maxshift_x,
                        int                  maxshift_y);
    //! Send the coordinates to our PME-only rank and request a PME calculation
    void sendCoordinates(DeviceBuffer<RVec>    coordinates,
                         const matrix          box,
                         ArrayRef<const RVec>  x,
                         real                  lambda_q,
                         real                  lambda_lj,
                         bool                  computeEnergyAndVirial,
                         int64_t               step,
                         bool                  sendCoordinatesFromGpu,
                         bool                  receiveForcesToGpu,
                         GpuEventSynchronizer* coordinatesReadyOnDeviceEvent,
                         bool                  useMdGpuGraph,
                         gmx_wallcycle*        wcycle);
    /*! \brief Receive results from PME-only ranks
     *
     * All PP ranks receive forces. Only one PP rank per PME rank
     * receives non-zero energy and virial contributions, and only
     * on the appropriate steps. */
    void receiveResults(ForceWithVirial* forceWithVirial,
                        real*            energy_q,
                        real*            energy_lj,
                        real*            dvdlambda_q,
                        real*            dvdlambda_lj,
                        bool             receivePmeForceToGpu,
                        float*           pme_cycles);
    //! Tell our PME-only rank to finish
    void sendFinish() const;
    //! Tell our PME-only rank to reset all cycle and flop counters
    void sendResetCounters(int64_t step) const;
    //! Tell our PME-only rank to switch to a new grid size
    void sendSwitchGrid(ivec grid_size, real ewaldcoeff_q, real ewaldcoeff_lj) const;
    //! Return the MPI rank of the partner PME rank controlled by this PP rank, if any
    std::optional<int> rankOfControlledPmeRank() const;
    //! Get handle to direct-GPU PME-PP comm object, if active
    PmePpCommGpu* pmePpCommGpu();

    /*! \brief Deleted move and copy operations
     *
     * These would be potentially unsafe during ongoing MPI
     * communication involving \c cnb_, and not useful, so
     * deleted. */
    PmePpComm(const PmePpComm&)            = delete;
    PmePpComm& operator=(const PmePpComm&) = delete;
    PmePpComm(PmePpComm&&)                 = delete;
    PmePpComm& operator=(PmePpComm&&)      = delete;

private:
    //! Receive virial and energy from PME rank
    void receiveVirialAndEnergy(ForceWithVirial* forceWithVirial,
                                real*            energy_q,
                                real*            energy_lj,
                                real*            dvdlambda_q,
                                real*            dvdlambda_lj,
                                float*           pme_cycles);
    //! Receive force data from PME ranks
    void receiveForces(bool receivePmeForceToGpu);

    //! Communicator between PME and PP ranks
    const MpiComm& comm_;
    /*! \brief Number of atoms home to this PP rank's domain
     *
     * Must be updated every domain repartitioning. */
    int numHomeAtoms_;
    //! Rank of partner PME rank
    int rankOfPartnerPmeRank_;
    //! The type of Van der Waals interaction
    VanDerWaalsType vdwType_;
    //! The type of Coulomb interaction
    CoulombInteractionType coulombType_;
    /*! \brief Whether this rank receives virial and energy from
     * the PME rank
     *
     * When multiple PP ranks interact with a PME-only rank, only one
     * receives the virial and energy contributions for all of
     * them. It also sends the control signals to the PME-only
     * rank. */
    bool thisRankReceivesVirialAndEnergy_;
    //! If direct PP-PME communication between GPUs is used.
    bool useGpuPmePpCommunication_;
    /*! \brief Struct used by the controlling PP rank
     *
     * Contains values sent to the PME rank with non-blocking
     * messages, so needs to remain in scope until the matching
     * MPI wait call is made. */
    std::optional<gmx_pme_comm_n_box_t> cnb_;
    //! MPI requests for communication with the PME rank
    std::vector<MPI_Request> requests_;
    //! CPU buffer into which this PP rank receives PME forces, possibly from a GPU
    HostVector<RVec> cpuPmeForceReceiveBuffer_;
    //! Manager of direct GPU PME-PP communication, when active.
    std::optional<PmePpCommGpu> pmePpCommGpu_;
};

} // namespace gmx

#endif
