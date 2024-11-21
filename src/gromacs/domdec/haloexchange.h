/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

/*! \libinternal \file
 *
 * \brief Declares the HaloExchange class for halo communication of coordinates and forces
 *
 * The algorithm this class implements uses direct communication between the home domain
 * and all domains in the, up to 7, zones. This contrasts with the original eighth shell
 * domain decomposition which only communicated with neighbor domains along the three
 * dimensions and aggregrated data over multiple pulses to communicate with more distant
 * domains. This is also referred to as staged communication. The setup here can be
 * referred to as "direct" instead, not to be confused with "GPU direct"
 * communication.
 *
 * The HaloExchange class contains a list of DomainPairComm objects, one for every domain
 * in the halo. The communication setup and calls for each domain pair are completely
 * independent from the other pairs. This make the class setup simple and also allows
 * for overlap of all communication for the different domains. The only coupling between
 * DomainPairComm objects is that they together should provide the atoms necessary
 * to compute all required non-bonded and bonded interactions in the system.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_HALOEXCHANGE_H
#define GMX_DOMDEC_HALOEXCHANGE_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/range.h"

struct gmx_ddbox_t;
struct gmx_domdec_t;
struct nonbonded_verlet_t;
enum class PbcType : int;
struct t_forcerec;
class t_state;

namespace gmx
{
template<typename>
class ArrayRef;
class DomainPairComm;

//! Storage for MPI request for halo MPI receive and send operations
struct HaloMpiRequests
{
    //! Receive requests for all domain pairs in the halo
    std::vector<MPI_Request> receive;
    //! Send requests for all domain pairs in the halo
    std::vector<MPI_Request> send;
};

//! Handles the halo communication of coordinates and forces
class HaloExchange
{
public:
    HaloExchange(PbcType pbcType);

    ~HaloExchange();

    //! Set up the halo communication, should be called after (re)partitioning
    void setup(gmx_domdec_t* dd, t_state* localState, const gmx_ddbox_t& ddbox, t_forcerec* fr, const bool cellsChanged);

    /*! \brief Iniatiate a non-blocking receive of the halo coordinates x
     *
     * Should be called as early as possible for fast communication.
     */
    void initiateReceiveX(ArrayRef<RVec> x);

    /*! \brief Initiate a non-blocking send of the halo coordinates x
     *
     * Should be called as early as possible when x is ready.
     */
    void initiateSendX(const matrix box, ArrayRef<RVec> x);

    /*! \brief Complete a non-blocking receive of coordinates for the complete halo
     *
     * Needs to be called after initiateRecvX(). Call as late as possible.
     */
    void completeReceiveX();

    /*! \brief Ensures that x has been sent
     *
     * Needs to be called after initiateSendX(). Call as late as possible.
     */
    void completeSendX();

    //! Send and receive the halo coordinates, call the 4 methods above
    void moveX(const matrix box, ArrayRef<RVec> x);

    /*! \brief Initiate A non-blocking receive of the halo forces f
     *
     * Should be called as early as possible for fast communication.
     */
    void initiateReceiveF();

    /*! \brief Initiate the non-blocking receive of the halo forces f
     *
     * Should be called as early as possible for fast communication.
     */
    void initiateSendF(ArrayRef<const RVec> f);

    /*! \brief Ensures that f has been received and reduced the received forces
     *
     * If \p shiftForces is not empty, also updates the shift forces.
     * Needs to be called after initiateReceiveF(). Call as late as possible.
     */
    void completeReceiveF(ArrayRef<RVec> forces, ArrayRef<RVec> shiftForces);

    /*! \brief Ensures that f has been sent
     *
     * Needs to be called after initiateSendF(). Call as late as possible.
     */
    void completeSendF();

    //! Send and receive the halo force, accumulates shift forces to \p shiftForces when non-empty
    void moveF(ArrayRef<RVec> f, ArrayRef<RVec> shiftForces);

private:
    //! Checks whether the allocation range is sufficient, if not: re-initializes \p domainPairComm_
    void checkDomainRangeAllocation(const gmx_domdec_t& dd, const IVec& domainRange);

    //! The type of PBC
    PbcType pbcType_;

    //! List of objects for communicating between pairs of domains
    std::vector<DomainPairComm> domainPairComm_;
    //! The grid range of the entries in \p domainPairComm_, index 0,0,0 is not present
    IVec allocatedDomainRange_ = { 0, 0, 0 };

    //! List of MPI requests for coordinate communiation
    HaloMpiRequests mpiCoordinateRequests_;
    //! List of MPI requests for force communiation
    HaloMpiRequests mpiForceRequests_;

    //! Buffer for storing the MPI status for MPI waits
    std::vector<MPI_Status> mpiStatus_;
};

} // namespace gmx

#endif // GMX_DOMDEC_HALOEXCHANGE_H
