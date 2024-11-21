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

/*! \internal \file
 *
 * \brief Declares classes for communication between domain pairs in the halo
 *
 * The DomainPairComm class is a holder for an object of type DomainPairCommBackward
 * and an object of type DomainPairCommForward. The backward object is used to send
 * coordinates backward along the DD grid and receive forces from the same rank/domain
 * the coordinates were sent to. The forward object is used for receiving the coordinates
 * and sending the forces in the other direction. Both objects couple to a ranks/domains
 * with the same displacement vector in DD-grid coordinates but with opposite sign.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMAINPAIRCOMM_H
#define GMX_DOMDEC_DOMAINPAIRCOMM_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/nbnxm/grid.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/defaultinitializationallocator.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"

enum class PbcType : int;
struct gmx_domdec_t;

namespace gmx
{

/*! \brief MPI tags for non-blocking x and f communication.
 *
 * With the current call order we don't need this.
 * But it's safer to have them, in case one would e.g. like to post
 * the force receive init before the coordinate communication is completed.
 */
enum class HaloMpiTag
{
    X,              //!< Coordinates
    F,              //!< Forces
    GridCounts,     //! The number of grid columns
    GridColumns,    //! The contents of the grid column
    GridDimensions, //! The dimensions of the grid
    AtomIndices     //! Global atom indices
};

/*! \brief Setup for selecting halo atoms to be sent and sending coordinates to another domain
 *
 * Also used for receiving forces.
 *
 * This object is used for communicating with a domain that resides in backward direction along
 * the domain decomposition grid. This object holds the NBNxM grid cell indices that needs to
 * be sent and does the packing of coordinates into a buffer in this object. It also does
 * the reverse indexing for the halo force reduction of forces received.
 */
class DomainCommBackward
{
public:
    //! Struct for collecting information on grid columns
    struct ColumnInfo
    {
        //! The index of the grid column
        int index;
        //! The cell range to communicate
        gmx::Range<int> cellRange;
    };

    /*! \brief Constructor
     *
     * \param[in] rank         The MPI rank we communicate with
     * \param[in] zone         The domain decomposition zone this pair of domains belongs to
     * \param[in] domainShift  The shift >=0 in domain indices between the two domains
     * \param[in] pbcType      The type of PBC
     * \param[in] commOverPbc  Whether we are communicating over a periodic boundary
     * \param[in] pbcCoordinateShift  The PBC coordinate shift, 0 or 1 for each dimension
     * \param[in] mpiCommAll   MPI communicator for all PP ranks
     */
    DomainCommBackward(int         rank,
                       int         zone,
                       const IVec& domainShift,
                       PbcType     pbcType,
                       bool        commOverPbc,
                       const IVec& pbcCoordinateShift,
                       MPI_Comm    mpiCommAll);

    //! Clears this communication, set no columns and atoms to send
    void clear();

    /*! \brief Determine which NBNxM grid cells (and atoms) we need to send
     *
     * \param[in] dd                      The domain decomposition struct
     * \param[in] grid                    The local NBNxM pair-search grid
     * \param[in] cutoffSquaredNonbonded  The cutoff^2 for non-bonded interactions
     * \param[in] cutoffSquaredBonded     The cutoff^2 for bonded interactions
     * \param[in] box                     The box
     * \param[in] dimensionIsTriclinic    Tells whether the dimensions require
     *                                    triclinic distance checks
     * \param[in] normal                  The normal vectors to planes separating domains along
     *                                    the triclinic dimensions
     * \param[in] isCellMissingLinks      Tells whether cells are missing bonded interactions,
     *                                    only used when filtering bonded communication
     */
    void selectHaloAtoms(const gmx_domdec_t&      dd,
                         const Grid&              grid,
                         const real               cutoffSquaredNonbonded,
                         const real               cutoffSquaredBonded,
                         const matrix             box,
                         const ivec               dimensionIsTriclinic,
                         ArrayRef<const RVec>     normal,
                         const std::vector<bool>& isCellMissingLinks);

    //! Creates and returns a buffer with column indices and cell counts to be sent
    FastVector<std::pair<int, int>> makeColumnsSendBuffer() const;

    //! Copies the coordinates to commnicate to the send buffer
    void packCoordinateSendBuffer(const matrix box, ArrayRef<const RVec> x, ArrayRef<RVec> sendBuffer) const;

    //! Accumulates the received forces to the force buffer and shift force buffer, when not empty
    void accumulateReceivedForces(ArrayRef<RVec> forces, ArrayRef<RVec> shiftForces) const;

    //! Returns the rank we communicate with
    int rank() const { return rank_; }

    //! Return the zone this domain pair resides in
    int zone() const { return zone_; }

    //! Returns the shift in domains in Cartesian coordinates between the communicating domains
    int domainShift(int dim) const { return domainShift_[dim]; }

    //! Returns whether the domain shift is more than one domain along at least one dimension
    bool shiftMultipleDomains() const { return shiftMultipleDomains_; }

    //! Returns the number of atoms per cell
    int numAtomsPerCell() const { return numAtomsPerCell_; }

    //! Returns the list of all columns to send with column information
    ArrayRef<const ColumnInfo> columnsToSend() const { return columnsToSend_; }

    //! Returns the number of atoms to send
    int numAtoms() const { return numAtomsToSend_; }

    //! Returns the buffer with global atom indices to communicate
    ArrayRef<const int> globalAtomIndices() const { return globalAtomIndices_; }

    //! Returns the buffer for sending coordinates and receiving forces
    ArrayRef<RVec> rvecBuffer() { return rvecBuffer_; }

    //! Returns the buffer for sending coordinates and receiving forces
    ArrayRef<const RVec> rvecBuffer() const { return rvecBuffer_; }

    //! The MPI communication for halo exchange
    MPI_Comm mpiCommAll() const { return mpiCommAll_; }

private:
    //! The rank to communicate with
    int rank_;
    //! The zone this part of the halo belongs to
    int zone_;
    //! The shift in domains in Cartesian coordinates between the communicating domains (always >= 0)
    IVec domainShift_;
    //! The type of PBC
    PbcType pbcType_;
    //! Whether we communicate over PBC
    bool commOverPbc_;
    //! Whether we use screw PBC
    bool usesScrewPbc_;
    //! The number of box vectors to shift coordinate to communicate by
    IVec pbcCoordinateShift_;
    //! The shift vector index to accumulate shift forces to
    int pbcForceShiftIndex_;
    //! Whether the domain shift is more than one domain along at least one dimension
    bool shiftMultipleDomains_;
    //! The number of atoms per cell, note that this is max over i/j, unlike Grid which uses i
    int numAtomsPerCell_;
    //! The cell ranges to commnicate
    FastVector<ColumnInfo> columnsToSend_;
    //! The number of atoms to send (or receive in case of forces)
    int numAtomsToSend_;
    //! Buffer for communicating global atom indices
    FastVector<int> globalAtomIndices_;
    //! Buffer for sending coordinates and receiving forces
    FastVector<RVec> rvecBuffer_;

    //! The MPI communication for all PP ranks
    MPI_Comm mpiCommAll_;
};

/*! \brief Setup for receiving halo coordinates from another domain and sending halo forces
 *
 * This object is used for communicating with a domain that resides in forward direction along
 * the domain decomposition grid. The coordinates are received in place in the coordinate
 * vector on this rank. The forces to send are read directly from the force buffer on this rank.
 */
class DomainCommForward
{
public:
    /*! \brief Constructor
     *
     * \param[in] rank        The MPI rank we communicate with
     * \param[in] zone        The domain decomposition zone this pair of domains belongs to
     * \param[in] mpiCommAll  MPI communicator for all PP ranks
     */
    DomainCommForward(int rank, int zone, MPI_Comm mpiCommAll);

    /*! \brief Sets up the forward communication, receiving the count from a remote \p send
     *
     * \param[in] send                      The send data we need to send foward
     * \param[in] offsetInCoordinateBuffer  The offset in the coordinate buffer for receiving coordinates
     */
    void setup(const DomainCommBackward& send, int offsetInCoordinateBuffer);

    //! Return the remote rank we communicate with
    int rank() const { return rank_; }

    //! Return the zone this domain pair resides in
    int zone() const { return zone_; }

    //! Returns the list of pairs of column indices and cell counts that we receive
    ArrayRef<const std::pair<int, int>> columnsReceived() const { return columnsReceived_; }

    //! The number of atoms to receive
    int numAtoms() const { return numAtomsToReceive_; }

    //! The atom range we store received coordinates in and send forces from
    Range<int> atomRange() const { return atomRange_; }

    //! The MPI communicator for halo exchange
    MPI_Comm mpiCommAll() const { return mpiCommAll_; }

private:
    //! The rank we communicate with
    int rank_;
    //! The zone this part of the halo belongs to
    int zone_;
    //! Pairs of column indices and cell counts (matching the Grid cell size definition)
    FastVector<std::pair<int, int>> columnsReceived_;
    //! The number of atoms to receive
    int numAtomsToReceive_;

    //! The atom range to receive the coordinates in and send forces from
    Range<int> atomRange_;

    //! The MPI communication for all PP ranks
    MPI_Comm mpiCommAll_;
};

/*! \brief Setup for communication between pairs of domains, both backward and forward along the DD grid
 *
 * This object is a holder a DomainCommBackward and DomainCommForward object. Both objects communicate
 * along the same DD-grid displacement vector but with opposite direction.
 */
class DomainPairComm
{
public:
    /*! \brief Constructor
     *
     * \param[in] backwardRank  The MPI rank we communicate with in backward direction
     * \param[in] forwardRank   The MPI rank we communicate with in forward direction
     * \param[in] zone          The domain decomposition zone this pair of domains belongs to
     * \param[in] domainShift   The shift >=0 in domain indices between the two domains
     * \param[in] pbcType       The type of PBC
     * \param[in] commOverPbc   Whether we are communicating backward over a periodic boundary
     * \param[in] pbcCoordinateShift  The PBC coordinate shift, 0 or 1 for each dimension
     * \param[in] mpiCommAll    MPI communicator for all PP ranks
     */
    DomainPairComm(int         backwardRank,
                   int         forwardRank,
                   int         zone,
                   const IVec& domainShift,
                   PbcType     pbcType,
                   bool        commOverPbc,
                   IVec        pbcCoordinateShift,
                   MPI_Comm    mpiCommAll);

    //! Returns the object for communicating backward along the DD grid
    DomainCommBackward& backward() { return backward_; }

    //! Returns the object for communicating backward along the DD grid
    const DomainCommBackward& backward() const { return backward_; }

    //! Returns the object for communicating forward along the DD grid
    DomainCommForward& forward() { return forward_; }

    //! Returns the object for communicating forward along the DD grid
    const DomainCommForward& forward() const { return forward_; }

private:
    //! Communication backward along the DD grid
    DomainCommBackward backward_;
    //! Communication forward along the DD grid
    DomainCommForward forward_;
};

} // namespace gmx

#endif // GMX_DOMDEC_DOMAINPAIRCOMM_H
