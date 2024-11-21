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
 * \brief Implements the setup() method for the HaloCommunication class
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/grid.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"

#include "domainpaircomm.h"
#include "domdec_internal.h"
#include "ga2la.h"
#include "haloexchange.h"
#include "utility.h"

namespace gmx
{

namespace
{

/*! \brief Move data of type \p T forward or backward between zones
 *
 * Moves in the dimension indexed by ddDimensionIndex, either forward
 * (direction=dddirFoward) or backward (direction=dddirBackward).
 */
template<typename T>
void ddSendReceive(const DomainCommBackward& domainCommBackward,
                   const DomainCommForward&  domainCommForward,
                   const int                 direction,
                   const T*                  sendBuffer,
                   const int                 numElementsToSend,
                   T*                        receiveBuffer,
                   const int                 numElementsToReceive,
                   const HaloMpiTag          tag)
{
#if GMX_MPI
    const int sendRank =
            (direction == dddirForward ? domainCommForward.rank() : domainCommBackward.rank());
    const int receiveRank =
            (direction == dddirForward ? domainCommBackward.rank() : domainCommForward.rank());

    const int  mpiTag  = static_cast<int>(tag);
    MPI_Comm   mpiComm = domainCommForward.mpiCommAll();
    MPI_Status mpiStatus;
    if (numElementsToSend > 0 && numElementsToReceive > 0)
    {
        int gmx_unused ret = MPI_Sendrecv(const_cast<T*>(sendBuffer),
                                          numElementsToSend * sizeof(T),
                                          MPI_BYTE,
                                          sendRank,
                                          mpiTag,
                                          receiveBuffer,
                                          numElementsToReceive * sizeof(T),
                                          MPI_BYTE,
                                          receiveRank,
                                          mpiTag,
                                          mpiComm,
                                          &mpiStatus);

        GMX_ASSERT(ret == 0, "Expect success");
    }
    else if (numElementsToSend > 0)
    {
        MPI_Send(sendBuffer, numElementsToSend * sizeof(T), MPI_BYTE, sendRank, mpiTag, mpiComm);
    }
    else if (numElementsToReceive > 0)
    {
        MPI_Recv(receiveBuffer, numElementsToReceive * sizeof(T), MPI_BYTE, receiveRank, mpiTag, mpiComm, &mpiStatus);
    }
#else  // GMX_MPI
    GMX_UNUSED_VALUE(domainCommBackward);
    GMX_UNUSED_VALUE(domainCommForward);
    GMX_UNUSED_VALUE(direction);
    GMX_UNUSED_VALUE(sendBuffer);
    GMX_UNUSED_VALUE(numElementsToSend);
    GMX_UNUSED_VALUE(receiveBuffer);
    GMX_UNUSED_VALUE(numElementsToReceive);
    GMX_UNUSED_VALUE(tag);
#endif // GMX_MPI
}

//! Returns ddCoord with each value within (0, numDDCells[] - 1) by applying modulo
IVec pbcDDCoord(const IVec& ddCoord, const IVec& numDDCells)
{
    IVec v = ddCoord;

    for (int d = 0; d < DIM; d++)
    {
        v[d] = (v[d] + numDDCells[d]) % numDDCells[d];
    }

    return v;
}

/* Check if atom with global index at_gl has any bonded interactions
 * that involve non-local atoms.
 */
inline bool isMissingALink(const ArrayRef<const int>& links, const gmx_ga2la_t& ga2la)
{
    return std::any_of(
            links.begin(), links.end(), [&ga2la](int a) { return ga2la.findHome(a) == nullptr; });
}

//! Return a list of of bools for grid cells, which tell if a cell has non-local bonded interactions
std::vector<bool> flagCellsForBondcomm(const gmx_domdec_t&           dd,
                                       const ArrayRef<const int32_t> atinfo,
                                       const nonbonded_verlet_t&     nbv)
{
    const auto& bondedLinks = *dd.comm->bondedLinks;

    // Note: numAtomsPerCell should match what is used in DomainCommBackward::selectHaloAtoms()
    const int numAtomsPerCell = nbv.localGrid().geometry().numAtomsPerCell_;

    ArrayRef<const int> atomOrder = nbv.getLocalAtomOrder();

    const int numCells = atomOrder.ssize() / numAtomsPerCell;
    GMX_ASSERT(numCells * numAtomsPerCell == atomOrder.ssize(),
               "The number of atoms on the grid and in the sort struct should match");

    std::vector<bool> isCellMissingLinks(numCells);

    /* Loop over all cells of the local grid */
    for (int c = 0; c < numCells; c++)
    {
        /* Loop over all atoms in this cell */
        bool thisCellIsMissingALink = false;
        for (int i = c * numAtomsPerCell; i < (c + 1) * numAtomsPerCell && !thisCellIsMissingALink; i++)
        {
            const int a = atomOrder[i];

            /* Check if this is a real atom (not a filler atom) */
            if (a >= 0 && (atinfo[a] & sc_atomInfo_BondCommunication))
            {
                if (isMissingALink(bondedLinks[dd.globalAtomIndices[a]], *dd.ga2la))
                {
                    thisCellIsMissingALink = true;
                }
            }

            isCellMissingLinks[c] = thisCellIsMissingALink;
        }
    }

    return isCellMissingLinks;
}

//! Sets the atom info structures.
void ddSetAtominfo(ArrayRef<const int> globalAtomIndices, const int atomStart, int atomEnd, t_forcerec* fr)
{
    if (fr != nullptr)
    {
        ArrayRef<AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock =
                fr->atomInfoForEachMoleculeBlock;
        ArrayRef<int32_t> atomInfo = fr->atomInfo;

        for (int a = atomStart; a < atomEnd; a++)
        {
            const int globalAtomIndex = globalAtomIndices[a];
            if (globalAtomIndex >= 0)
            {
                atomInfo[a] = ddGetAtomInfo(atomInfoForEachMoleculeBlock, globalAtomIndex);
            }
            else
            {
                atomInfo[a] = sc_atomInfo_IsFillerParticle;
            }
        }
    }
}

//! Returns the maximum domain communication range, indexed by DD dimension index
IVec getDomainCommunicationRange(const gmx_domdec_t& dd)
{
    IVec range = { 0, 0, 0 };

    for (int d = 0; d < dd.ndim; d++)
    {
        range[d] = dd.numPulses[dd.dim[d]];
    }

    return range;
}

} // namespace

void HaloExchange::checkDomainRangeAllocation(const gmx_domdec_t& dd, const IVec& domainRange)
{
    bool realloc              = false;
    int  totNumDomainsInZones = 1;
    for (int dimIndex = 0; dimIndex < dd.ndim; dimIndex++)
    {
        if (domainRange[dimIndex] > allocatedDomainRange_[dimIndex])
        {
            allocatedDomainRange_[dimIndex] = domainRange[dimIndex];
            realloc                         = true;
        }
        totNumDomainsInZones *= 1 + allocatedDomainRange_[dimIndex];
    }
    if (!realloc)
    {
        return;
    }

    domainPairComm_.clear();

    for (int zone = 1; zone < dd.zones.numZones(); zone++)
    {
        const IVec& zoneShift = dd.zones.shift(zone);

        int numDomains = 1;
        for (int dimIndex = 0; dimIndex < dd.ndim; dimIndex++)
        {
            if (zoneShift[dd.dim[dimIndex]] > 0)
            {
                numDomains *= allocatedDomainRange_[dimIndex];
            }
        }

        for (int domainIndex = 0; domainIndex < numDomains; domainIndex++)
        {
            IVec domainShift = { 0, 0, 0 };
            int  rest        = domainIndex;
            for (int dimIndex = dd.ndim - 1; dimIndex >= 0; dimIndex--)
            {
                if (zoneShift[dd.dim[dimIndex]] > 0)
                {
                    domainShift[dd.dim[dimIndex]] = 1 + (rest % allocatedDomainRange_[dimIndex]);

                    rest /= allocatedDomainRange_[dimIndex];
                }
            }

            const IVec backwardCoord = pbcDDCoord(dd.ci - domainShift, dd.numCells);
            const IVec forwardCoord  = pbcDDCoord(dd.ci + domainShift, dd.numCells);
            const int  backwardRank  = ddRankFromDDCoord(dd, backwardCoord);
            const int  forwardRank   = ddRankFromDDCoord(dd, forwardCoord);

            const bool commOverPbc = (backwardCoord != dd.ci - domainShift);
            ivec       pbcShift;
            for (int d = 0; d < DIM; d++)
            {
                pbcShift[d] = (backwardCoord[d] != dd.ci[d] - domainShift[d]) ? 1 : 0;
            }
            if (debug)
            {
                fprintf(debug,
                        "zone %d domainShift %d %d %d pbc %d\n",
                        zone,
                        domainShift[0],
                        domainShift[1],
                        domainShift[2],
                        static_cast<int>(commOverPbc));
            }
            domainPairComm_.emplace_back(
                    backwardRank, forwardRank, zone, domainShift, pbcType_, commOverPbc, pbcShift, dd.mpi_comm_all);
        }
    }

    GMX_RELEASE_ASSERT(int(domainPairComm_.size()) == totNumDomainsInZones - 1,
                       "We should have as many comm entries as non-local domains");

    mpiStatus_.resize(totNumDomainsInZones - 1);
}

void HaloExchange::setup(gmx_domdec_t*         dd,
                         t_state*              localState,
                         const gmx_ddbox_t&    ddbox,
                         t_forcerec*           fr,
                         const bool gmx_unused cellsChanged)
{
    const IVec domainRange = getDomainCommunicationRange(*dd);

    if (debug)
    {
        fprintf(debug, "Setting up DD communication, range");
        for (int dimIndex = 0; dimIndex < dd->ndim; dimIndex++)
        {
            fprintf(debug, " %d", domainRange[dimIndex]);
        }
        fprintf(debug, "\n");
    }

    checkDomainRangeAllocation(*dd, domainRange);

    gmx_domdec_comm_t& comm = *dd->comm;

    GMX_RELEASE_ASSERT(!isDlbOn(comm.dlbState), "DLB is not supported here yet");

    /* The naming is not fully consistent here (yet).
     * But the 2-body bonded interaction cut-off is max(cutoff, cutoff_mbody),
     * so there is not short and exact naming.
     */
    const real cutoffSquaredNonbonded = gmx::square(comm.systemInfo.cutoff);
    const real cutoffSquaredBonded    = gmx::square(comm.cutoff_mbody);

    /* Determine the normals of the zone planes */
    std::array<RVec, 3> normal;
    for (int i = 0; i < DIM; i++)
    {
        if (ddbox.tric_dir[i])
        {
            /* ddbox->normal has length skew_fac, normalize it */
            svmul(1 / ddbox.skew_fac[i], ddbox.normal[i], normal[i]);
        }
        else
        {
            clear_rvec(normal[i]);
            normal[i][i] = 1;
        }
    }

    DomdecZones& zones = dd->zones;

    auto& atomInfo = fr->atomInfo;
    auto& nbv      = *fr->nbv;

    std::vector<bool> missingLinkInCells;
    if (comm.systemInfo.filterBondedCommunication)
    {
        /* Create a boolean array for cells telling if bondeds linked to atoms
         * are not locally present, so we need to communicate those cells.
         */
        missingLinkInCells = flagCellsForBondcomm(*dd, atomInfo, nbv);
    }

    zones.setAtomRangeEnd(0, dd->numHomeAtoms, true);
    comm.atomRanges.setEnd(DDAtomRanges::Type::Home, dd->numHomeAtoms);

    /* Here we distribute the comm setup calculation over threads.
     * Note that selectHaloAtomsForDomainPair() is actually not very time consuming.
     * Alternatively we could overlap selectHaloAtomsForZone() with comm_zone_comm_setup
     * for the previous zone.
     */
    const int gmx_unused nthread = gmx_omp_nthreads_get(ModuleMultiThread::Domdec);
#pragma omp parallel for num_threads(nthread) schedule(static, 1)
    for (Index dpcIndex = 0; dpcIndex < gmx::ssize(domainPairComm_); dpcIndex++)
    {
        DomainPairComm& dpc = domainPairComm_[dpcIndex];

        DomainCommBackward& send = dpc.backward();

        bool domainIsInRange = true;

        for (int dimIndex = 0; dimIndex < dd->ndim; dimIndex++)
        {
            if (send.domainShift(dd->dim[dimIndex]) > domainRange[dimIndex])
            {
                // This domain is currently out of range, we can skip the selection
                domainIsInRange = false;
            }
        }

        if (domainIsInRange)
        {
            // Determine which atoms we need to send
            send.selectHaloAtoms(*dd,
                                 nbv.localGrid(),
                                 cutoffSquaredNonbonded,
                                 cutoffSquaredBonded,
                                 localState->box,
                                 ddbox.tric_dir,
                                 normal,
                                 missingLinkInCells);
        }
        else
        {
            send.clear();
        }
    }

    // We need a separate, non-threaded, loop for the MPI communication
    int numAtomsTotal = dd->numHomeAtoms;
    for (DomainPairComm& dpc : domainPairComm_)
    {
        DomainCommBackward& send    = dpc.backward();
        DomainCommForward&  receive = dpc.forward();

        const int zone = send.zone();

        // This communication could be overlapped with selectHaloAtomForZone() for the next zone
        receive.setup(send, numAtomsTotal);

        numAtomsTotal += receive.numAtoms();

        dd->globalAtomIndices.resize(numAtomsTotal);

        // Communicate the global atom indices
        ddSendReceive(send,
                      receive,
                      dddirBackward,
                      send.globalAtomIndices().data(),
                      send.numAtoms(),
                      dd->globalAtomIndices.data() + *receive.atomRange().begin(),
                      receive.numAtoms(),
                      HaloMpiTag::AtomIndices);

        // Update the zone atom count, will be updated again with multiple domains in this zone
        zones.setAtomRangeEnd(zone, numAtomsTotal, !send.shiftMultipleDomains());
    }

    comm.atomRanges.setEnd(DDAtomRanges::Type::Zones, numAtomsTotal);

    // Set atominfo, so we can pass it to setNonLocalGrid()
    atomInfo.resize(numAtomsTotal);
    ddSetAtominfo(dd->globalAtomIndices, dd->numHomeAtoms, numAtomsTotal, fr);

    // Communicate the coordinates, needed for calculating the grid cell bounding boxes
    localState->changeNumAtoms(numAtomsTotal);
    moveX(localState->box, localState->x);
    nbv.nbat().resizeCoordinateBuffer(numAtomsTotal, dd->zones.numZones() - 1);

    // Set up the non-local grids using the communicated cells.
    int nbnxmGridIndex = 1;
    for (const DomainPairComm& dpc : domainPairComm_)
    {
        const DomainCommBackward& send    = dpc.backward();
        const DomainCommForward&  receive = dpc.forward();

        /* Communicate the grid geometry */
        GridDimensions gridDimensions;
        ddSendReceive(
                send, receive, dddirBackward, &nbv.localGrid().dimensions(), 1, &gridDimensions, 1, HaloMpiTag::GridDimensions);

        /* Apply PBC on the receiver side */
        for (int dimIndex = 0; dimIndex < dd->ndim; dimIndex++)
        {
            const int dim = dd->dim[dimIndex];

            if (dim < dd->unitCellInfo.npbcdim && dd->zones.shift(receive.zone())[dim] == 1
                && dd->ci[dim] >= dd->numCells[dim] - send.domainShift(dim))
            {
                /* We are communicating over pbc along dim */
                gridDimensions.lowerCorner += localState->box[dim];
                gridDimensions.upperCorner += localState->box[dim];
            }
        }

        // Add 1 to domainIndex as the grid with index 0 is the local grid
        nbv.setNonLocalGrid(nbnxmGridIndex,
                            receive.zone(),
                            gridDimensions,
                            receive.columnsReceived(),
                            fr->atomInfo,
                            localState->x);
        nbnxmGridIndex++;
    }

    if (debug)
    {
        fprintf(debug, "Finished setting up DD communication, domain atom receive counts:");
        int previousZone = 0;
        for (const DomainPairComm& dpc : domainPairComm_)
        {
            const DomainCommForward& receive = dpc.forward();

            if (receive.zone() != previousZone)
            {
                fprintf(debug, "\n");
                fprintf(debug, "zone %d:", receive.zone());
                previousZone = receive.zone();
            }
            fprintf(debug, " %d", receive.numAtoms());
        }
        fprintf(debug, "\n");
    }

    nbv.nbat().resizeForceBuffers();
}

} // namespace gmx
