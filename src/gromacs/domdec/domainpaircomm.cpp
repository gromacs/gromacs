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
 * \brief Implements the DomainPairComm, DomainCommBackward and DomainCommForward classes
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "domainpaircomm.h"

#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/vec.h"
#include "gromacs/nbnxm/grid.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fixedcapacityvector.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/range.h"

#include "domdec_internal.h"
#include "utility.h"

namespace gmx
{

namespace
{

//! The upper corners of a zone, used for computing which halo cell need to be sent
struct ZoneCorners
{
    //! Corner for two-body interations, involves all pair-interacting zones
    RVec twoBody = { 0.0_real, 0.0_real, 0.0_real };
    //! Corner for multi-body interactions, involves all zones
    RVec multiBody = { 0.0_real, 0.0_real, 0.0_real };
    //! The corner of our own zone
    RVec zone = { 0.0_real, 0.0_real, 0.0_real };
    //! Whether \p twoBody and \p multiBody differ
    bool cornersDiffer = false;
};

struct GridColumnInfo
{
    //! The bounding box of the column
    BoundingBox columnBB;
    //! Bounding boxes for the cells, used with CPU grids, otherwise empty
    ArrayRef<const BoundingBox> cellBBs;
    //! Bounding boxes along Z only, used with GPU grids, otherwise empty
    ArrayRef<const BoundingBox1D> cellBBsZonly;
    //! The number of i-cells per bounding box
    int bbToCellFactor;
};

//! Returns information on a column of the local grid
GridColumnInfo nbnxmGetLocalGridColumn(const Grid& grid, const int columnIndex)
{
    const GridDimensions& dims = grid.dimensions();

    GMX_ASSERT(columnIndex >= 0 && columnIndex < grid.numColumns(),
               "columnIndex should be in range");

    // Determine the x-column through division with rounding down
    const int cx = columnIndex / dims.numCells[YY];
    const int cy = columnIndex - cx * dims.numCells[YY];

    GridColumnInfo gci;

    // Atoms can be at most out of the column by a distance of maxAtomGroupRadius
    const real magr = grid.dimensions().maxAtomGroupRadius;

    gci.columnBB.lower.x = dims.lowerCorner[XX] + cx * dims.cellSize[XX] - magr;
    gci.columnBB.upper.x = dims.lowerCorner[XX] + (cx + 1) * dims.cellSize[XX] + magr;
    gci.columnBB.lower.y = dims.lowerCorner[YY] + cy * dims.cellSize[YY] - magr;
    gci.columnBB.upper.y = dims.lowerCorner[YY] + (cy + 1) * dims.cellSize[YY] + magr;
    gci.columnBB.lower.z = dims.lowerCorner[ZZ] - magr;
    gci.columnBB.upper.z = dims.upperCorner[ZZ] + magr;

    const Grid::Geometry& geometry = grid.geometry();

    gci.bbToCellFactor = 1;

    if (geometry.isSimple_)
    {
        int numBBs = grid.cxy_ind()[columnIndex + 1] - grid.cxy_ind()[columnIndex];
        // We return the largest of the x/y bounding boxes
        if (geometry.numAtomsJCluster_ <= geometry.numAtomsICluster_)
        {
            gci.cellBBs = grid.iBoundingBoxes().subArray(grid.cxy_ind()[columnIndex], numBBs);
        }
        else
        {
            GMX_ASSERT(geometry.numAtomsJCluster_ == 2 * geometry.numAtomsICluster_,
                       "Currently only equal i/j cells and j-cell size double i-cell size are "
                       "supported");
            // We have 2 i-cells per bounding box
            gci.bbToCellFactor = 2;
            // j-clusters are twice as large as i, need to divide counts by 2
            numBBs /= 2;
            gci.cellBBs = grid.jBoundingBoxes().subArray(grid.cxy_ind()[columnIndex] / 2, numBBs);
        }
    }
    else
    {
        const int numBBs = grid.cxy_ind()[columnIndex + 1] - grid.cxy_ind()[columnIndex];
        gci.cellBBsZonly = grid.zBoundingBoxes().subArray(grid.cxy_ind()[columnIndex], numBBs);
    }

    return gci;
}

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
        int ret = MPI_Sendrecv(const_cast<T*>(sendBuffer),
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
        GMX_UNUSED_VALUE(ret);
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

//! Information for distance calculations between domains
struct DistanceCalculationInfo
{
    //! Whether this zone is shifted wrt zone 0 for Cartesian dimensions
    IVec zoneShift = { 0, 0, 0 };
    //! Whether the unit cell is triclinic
    bool isTriclinic = false;
    //! Normals for the unit cell faces
    ArrayRef<const RVec> normal;
    //! Whether we should take into account non-orthogonal normals of faces
    IVec sumSquares = { 0, 0, 0 };
    //! Whether we need to check distances for multi-body interactions when selecting cells
    bool checkMultiBodyDistance = false;
    //! Whether we need to check distances for two-body bonded interactions when selecting cells
    bool checkTwoBodyDistance = false;
    //! Cut-off squared for non-bonded interaction distances
    real cutoffSquaredNonbonded;
    //! Cut-off squared for bonded interaction distances
    real cutoffSquaredBonded;

    // Whether to separately filter atoms for communication for long-distance bonded interactions.
    bool filterBondComm = false;
};

//! Struct for returning results from distance calculation of corners to bounding boxes
struct DistancesSquared
{
    //! The maximum distance squared for pair interactions
    real pair = 0;
    //! The maximum distance squared for multi-body interactions
    real multiBody = 0;
};

//! Returns the squared distances for non-bonded and bonded interactions of a bounding box to zone corners
DistancesSquared cornerToBoundingBoxDistanceRectangular(const DistanceCalculationInfo& dci,
                                                        const RVec&        cornerTwoBody,
                                                        const RVec&        cornerMultiBody,
                                                        const BoundingBox& bb)
{
    /* Here we use rounding and calculate the actual distance^2
     * to the corner(s).
     */
    DistancesSquared ds;

    /* Rectangular unit-cell, easy */
    for (int d = 0; d < DIM; d++)
    {
        if (dci.zoneShift[d])
        {
            /* Distance to the corner for pair-wise interactions */
            real r = bb.lower.ptr()[d] - cornerTwoBody[d];
            if (r > 0)
            {
                ds.pair += r * r;
            }
            if (dci.checkMultiBodyDistance)
            {
                /* Distance to the corner for multi-body interactions */
                r = bb.lower.ptr()[d] - cornerMultiBody[d];
                if (r > 0)
                {
                    ds.multiBody += r * r;
                }
            }
        }
    }

    return ds;
}

//! Returns the squared distances for non-bonded and bonded interactions of a bounding box to zone corners
DistancesSquared cornerToBoundingBoxDistanceTrilinic(const DistanceCalculationInfo& dci,
                                                     const RVec&                    cornerTwoBody,
                                                     const RVec&                    cornerMultiBody,
                                                     const BoundingBox&             bb)
{
    /* Here we use partial and approximate rounding */
    DistancesSquared ds;

    for (int d = 0; d < DIM; d++)
    {
        if (dci.zoneShift[d])
        {
            real r = 0;
            for (int d2 = d; d2 < DIM; d2++)
            {
                r += ((dci.normal[d][d2] >= 0 ? bb.lower.ptr()[d2] : bb.upper.ptr()[d2]) - cornerTwoBody[d2])
                     * dci.normal[d][d2];
            }
            if (r > 0)
            {
                if (dci.sumSquares[d])
                {
                    /* The angle(s) between the normal of this zone plane
                     * an the preceding zone plane(s) is/are >= 90 degrees.
                     * Add up the squares of the distance. This underestimates
                     * the distance for angles > 90 degrees.
                     */
                    ds.pair += r * r;
                }
                else
                {
                    /* The/A angle between the normals is < 90 degrees.
                     * We use the maximum distance, which is an underestimate.
                     */
                    ds.pair = std::max(ds.pair, r * r);
                }
            }

            if (dci.checkMultiBodyDistance)
            {
                r = 0;
                for (int d2 = d; d2 < DIM; d2++)
                {
                    r += ((dci.normal[d][d2] >= 0 ? bb.lower.ptr()[d2] : bb.upper.ptr()[d2])
                          - cornerMultiBody[d2])
                         * dci.normal[d][d2];
                }
                if (r > 0)
                {
                    if (dci.sumSquares[d])
                    {
                        ds.multiBody += r * r;
                    }
                    else
                    {
                        ds.multiBody = std::max(ds.multiBody, r * r);
                    }
                }
            }
        }
    }

    return ds;
}

/*! \brief Wrapper function for corner - bounding-box distance calculation.
 *
 * Only splits triclinic vs non-triclinic distance calculations.
 */
DistancesSquared cornerToBoundingBoxDistance(const DistanceCalculationInfo& dci,
                                             const RVec&                    cornerTwoBody,
                                             const RVec&                    cornerMultiBody,
                                             const BoundingBox&             bb)
{
    if (dci.isTriclinic)
    {
        return cornerToBoundingBoxDistanceTrilinic(dci, cornerTwoBody, cornerMultiBody, bb);
    }
    else
    {
        return cornerToBoundingBoxDistanceRectangular(dci, cornerTwoBody, cornerMultiBody, bb);
    }
}

//! Determines the corner for 2-body, corner_2b, and multi-body, corner_mb, communication distances
ZoneCorners getZoneCorners(const gmx_domdec_t& dd, const matrix box, const int zone, const real maxAtomGroupRadius)
{
    const DomdecZones& zones = dd.zones;

    const bool                                canHaveGridJumps = isDlbOn(dd.comm->dlbState);
    FixedCapacityVector<int, sc_maxNumIZones> interactingIZones;
    if (canHaveGridJumps)
    {
        /* Make a list of i-zones that see our zone on the receiving end */
        for (int iZone = 0; iZone < zones.numIZones(); iZone++)
        {
            if (zones.jZoneRange(iZone).isInRange(zone))
            {
                interactingIZones.push_back(iZone);
            }
        }
    }

    GMX_RELEASE_ASSERT(!canHaveGridJumps, "Grid jumps are not supported here yet");

    ZoneCorners zc;

    for (int dim = 0; dim < DIM; dim++)
    {
        /* This is the zone corner (simple, no staggering) */
        const real corner_z_d = dd.comm->cell_x0[dim];

        /* No staggering, all bounds are equal to our local bounds */
        const real corner_2b_d = dd.comm->cell_x0[dim];
        const real corner_mb_d = dd.comm->cell_x0[dim];

        zc.cornersDiffer = (corner_mb_d != corner_2b_d);

        /* Add the off-diagonal couplings to the corners */
        for (int d = 0; d < DIM; d++)
        {
            /* With 1D domain decomposition the cg's are not in
             * a triclinic box, but triclinic x-y and rectangular y/x-z.
             * So we should ignore the coupling for the non
             * domain-decomposed dimension of the pair x and y.
             */
            if (!(dd.ndim == 1 && ((dd.dim[0] == XX && d == YY) || (dd.dim[0] == YY && d == XX))))
            {
                zc.twoBody[d] += corner_2b_d * box[dim][d] / box[dim][dim];
                zc.multiBody[d] += corner_mb_d * box[dim][d] / box[dim][dim];
                zc.zone[d] += corner_z_d * box[dim][d] / box[dim][dim];
            }
        }
    }

    // Atoms can be outside the zone boundary by at most a distance of maxAtomGroupRadius
    for (int d = 0; d < DIM; d++)
    {
        zc.twoBody[d] += maxAtomGroupRadius;
        zc.multiBody[d] += maxAtomGroupRadius;
        zc.zone[d] += maxAtomGroupRadius;
    }

    if (debug)
    {
        fprintf(debug,
                "halo corners home %5.2f %5.2f %5.2f zone %d 2b %5.2f %5.2f %5.2f mb %5.2f %5.2f "
                "%5.2f\n",
                dd.comm->cell_x0[XX],
                dd.comm->cell_x0[YY],
                dd.comm->cell_x0[ZZ],
                zone,
                zc.twoBody[XX],
                zc.twoBody[YY],
                zc.twoBody[ZZ],
                zc.multiBody[XX],
                zc.multiBody[YY],
                zc.multiBody[ZZ]);
    }

    return zc;
}

//! Computes and sets the cell range we will communicatie for grid column \p columnIndex
template<bool doChecksForBondeds>
Range<int> getCellRangeForGridColumn(const Grid&                    grid,
                                     const int                      columnIndex,
                                     const ZoneCorners&             zoneCorners,
                                     const DistanceCalculationInfo& dci,
                                     const std::vector<bool>&       isCellMissingLinks)
{
    const GridColumnInfo gci = nbnxmGetLocalGridColumn(grid, columnIndex);

    const int numCellsInColumn = std::max(gci.cellBBs.ssize(), gci.cellBBsZonly.ssize());
    if (numCellsInColumn == 0)
    {
        /* Empty column */
        return {};
    }

    const auto distancesSquared =
            cornerToBoundingBoxDistance(dci, zoneCorners.twoBody, zoneCorners.multiBody, gci.columnBB);

    const bool columnIsInRange =
            distancesSquared.pair < dci.cutoffSquaredNonbonded
            || (dci.checkMultiBodyDistance && distancesSquared.multiBody < dci.cutoffSquaredBonded)
            || (dci.checkTwoBodyDistance && distancesSquared.pair < dci.cutoffSquaredBonded);
    if (!columnIsInRange)
    {
        /* This whole column is out of range */
        return {};
    }

    const bool         useColumnBBForXY = gci.cellBBs.empty();
    BoundingBox        bbFull;
    const BoundingBox* bbPointer;
    if (useColumnBBForXY)
    {
        /* Take the x & y components from the column */
        bbFull    = gci.columnBB;
        bbPointer = &bbFull;
    }

    // The cells we operate on here are the largest of the i- and j-cells,
    // so we need to convert the size used in grid, which is always of the i-cells
    const int firstCellInColumn = grid.firstCellInColumn(columnIndex) / gci.bbToCellFactor;

    int  firstCell = 0;
    int  lastCell  = numCellsInColumn - 1;
    bool isInRange;
    do
    {
        if (useColumnBBForXY)
        {
            bbFull.lower.z = gci.cellBBsZonly[lastCell].lower;
            bbFull.upper.z = gci.cellBBsZonly[lastCell].upper;
        }
        else
        {
            bbPointer = &gci.cellBBs[lastCell];
        }

        const auto distancesSquared = cornerToBoundingBoxDistance(
                dci, zoneCorners.twoBody, zoneCorners.multiBody, *bbPointer);

        /* Here we check:
         * the 2-atom distance against the non-bonded cut-off,
         * the multi-body distance against the bonded cut-off
         * The 2-atom distance against the bonded cut-off.
         * The bonded check only triggers communication without bBondComm
         * or when the cell has missing bonded interactions.
         */
        isInRange = (distancesSquared.pair < dci.cutoffSquaredNonbonded);

        if constexpr (doChecksForBondeds)
        {
            isInRange =
                    isInRange
                    || (((dci.checkMultiBodyDistance && distancesSquared.multiBody < dci.cutoffSquaredBonded)
                         || (dci.checkTwoBodyDistance && distancesSquared.pair < dci.cutoffSquaredBonded))
                        && (!dci.filterBondComm || isCellMissingLinks[firstCellInColumn + lastCell]));
        }

        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
        if (!isInRange)
        {
            lastCell--;
        }
    } while (lastCell >= firstCell && !isInRange);

    /* For rectangular grids without bondcomm we do not try to eliminate
     * cells from the bottom, since if lastCell is within range,
     * there is a high chance that firstCell=0 is also in range.
     */
    if ((dci.isTriclinic || dci.filterBondComm) && lastCell > firstCell)
    {
        /* This loop is a copy of the one above with lastCell replaced
         * by firstCell. Putting it in a function would be cleaner,
         * but this would require 20 parameters.
         */
        bool isInRange;
        do
        {
            if (!gci.cellBBs.empty())
            {
                bbPointer = &gci.cellBBs[firstCell];
            }
            else
            {
                bbFull.lower.z = gci.cellBBsZonly[firstCell].lower;
                bbFull.upper.z = gci.cellBBsZonly[firstCell].upper;
            }

            const auto distancesSquared = cornerToBoundingBoxDistance(
                    dci, zoneCorners.twoBody, zoneCorners.multiBody, *bbPointer);

            isInRange = (distancesSquared.pair < dci.cutoffSquaredNonbonded);

            if constexpr (doChecksForBondeds)
            {
                isInRange =
                        isInRange
                        || (((dci.checkMultiBodyDistance && distancesSquared.multiBody < dci.cutoffSquaredBonded)
                             || (dci.checkTwoBodyDistance && distancesSquared.pair < dci.cutoffSquaredBonded))
                            && (!dci.filterBondComm || isCellMissingLinks[firstCellInColumn + firstCell]));
            }

            // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
            if (!isInRange)
            {
                firstCell++;
            }
        } while (firstCell < lastCell && !isInRange);
    }

    return { (firstCellInColumn + firstCell) * gci.bbToCellFactor,
             (firstCellInColumn + lastCell + 1) * gci.bbToCellFactor };
}

} // namespace

DomainCommBackward::DomainCommBackward(int         rank,
                                       int         zone,
                                       const IVec& domainShift,
                                       PbcType     pbcType,
                                       bool        commOverPbc,
                                       const IVec& pbcCoordinateShift,
                                       MPI_Comm    mpiCommAll) :
    rank_(rank),
    zone_(zone),
    domainShift_(domainShift),
    pbcType_(pbcType),
    commOverPbc_(commOverPbc),
    usesScrewPbc_(false),
    pbcCoordinateShift_(pbcCoordinateShift),
    pbcForceShiftIndex_(ivecToShiftIndex(pbcCoordinateShift)),
    shiftMultipleDomains_(domainShift[XX] > 1 || domainShift[YY] > 1 || domainShift[ZZ] > 1),
    mpiCommAll_(mpiCommAll)
{
}

void DomainCommBackward::clear()
{
    columnsToSend_.clear();
    numAtomsToSend_ = 0;
}

void DomainCommBackward::selectHaloAtoms(const gmx_domdec_t&      dd,
                                         const Grid&              grid,
                                         const real               cutoffSquaredNonbonded,
                                         const real               cutoffSquaredBonded,
                                         const matrix             box,
                                         const ivec               dimensionIsTriclinic,
                                         ArrayRef<const RVec>     normal,
                                         const std::vector<bool>& isCellMissingLinks)
{
    if (pbcType_ == PbcType::XY && pbcCoordinateShift_[ZZ] != 0)
    {
        clear();

        return;
    }

    const gmx_domdec_comm_t& comm = *dd.comm;

    const DomdecZones& zones = dd.zones;

    const bool filterBondComm = comm.systemInfo.filterBondedCommunication;

    DistanceCalculationInfo dci;

    dci.cutoffSquaredNonbonded = cutoffSquaredNonbonded;
    dci.cutoffSquaredBonded    = cutoffSquaredBonded;

    for (int dimIndex = 0; dimIndex < dd.ndim; dimIndex++)
    {
        const int dim = dd.dim[dimIndex];

        if (zones.shift(zone_)[dim])
        {
            dci.zoneShift[dim] = 1;

            if (dimensionIsTriclinic[dim])
            {
                dci.isTriclinic = true;
            }
        }
    }
    dci.normal = normal;

    /* For triclinic zones with shifts along multiple unit-cell vectors,
     * the exact distance calculation gets very complex, since the normals
     * to the zone planes are not orthogonal. This makes rounding the edges
     * of those zones hard.
     * We apply approximate rounding in case the angles between the normals
     * are >= 90 degrees and no rouding when < 90 degrees. This leads to
     * a few more atoms communicated than strictly necessary, but results
     * in relatively simple, efficient and less bug-prone code.
     */
    int dim_count = 0;
    int dim_prev  = 0;
    for (int dim = 0; dim < DIM; dim++)
    {
        if (dci.zoneShift[dim])
        {
            switch (dim_count)
            {
                case 0:
                    /* First dimension, doesn't matter what we do here */
                    dci.sumSquares[dim] = 1;
                    break;
                case 1:
                    /* Second dimension, determine the angle between normals;
                     * if angle >= 90 degrees: sum squares of distances, rounds
                     *                         the zone, but not optimally, as
                     *                         the distance is underestimated
                     * if angle < 90 degrees:  use the maximum, no rounding,
                     *                         underestimates the distance.
                     */
                    dci.sumSquares[dim] = static_cast<int>(iprod(normal[dim], normal[dim_prev]) <= 0);
                    break;
                case 2:
                    /* As case 1, but check the angels with both planes */
                    GMX_ASSERT(dim == ZZ, "With 2 zone shifts, the dimension has to be z");
                    /* The normal along z is always (0,0,1) */
                    dci.sumSquares[dim] = static_cast<int>(normal[0][ZZ] <= 0 && normal[1][ZZ] <= 0);
                    break;
            }
            dim_prev = dim;
            dim_count++;
        }
    }

    /* Get the corners for non-bonded and bonded distance calculations */
    const ZoneCorners zoneCorners = getZoneCorners(dd, box, zone_, grid.dimensions().maxAtomGroupRadius);

    /* Do we need to determine extra distances for multi-body bondeds?
     * Note that with filterBondComm we might need distances longer than
     * the non-bonded cut-off, but with a grid without staggering (zoneCorners.cornersDiffer)
     * this check is indentical to the one triggered by checkTwoBodyDistance below.
     */
    dci.checkMultiBodyDistance =
            (comm.systemInfo.haveInterDomainMultiBodyBondeds && zoneCorners.cornersDiffer);

    /* Do we need to determine extra distances for only two-body bondeds? */
    dci.checkTwoBodyDistance = (filterBondComm && !dci.checkMultiBodyDistance);

    if (debug)
    {
        fprintf(debug,
                "setup zone %d BondComm %d DistMB %d Dist2B %d\n",
                zone_,
                static_cast<int>(filterBondComm),
                static_cast<int>(dci.checkMultiBodyDistance),
                static_cast<int>(dci.checkTwoBodyDistance));
    }

    numAtomsPerCell_ = grid.geometry().numAtomsPerCell_;

    // Clear the send counts
    clear();

    // We only need to check bonded distances when only single domain shifts are involved
    const bool checkBondedDistances =
            (!shiftMultipleDomains_ && (dci.checkMultiBodyDistance || dci.checkTwoBodyDistance));
    for (int columnIndex = 0; columnIndex < grid.numColumns(); columnIndex++)
    {
        Range<int> cellRange;

        if (checkBondedDistances)
        {
            cellRange = getCellRangeForGridColumn<true>(
                    grid, columnIndex, zoneCorners, dci, isCellMissingLinks);
        }
        else
        {
            cellRange = getCellRangeForGridColumn<false>(
                    grid, columnIndex, zoneCorners, dci, isCellMissingLinks);
        }
        if (!cellRange.empty())
        {
            columnsToSend_.push_back({ columnIndex, cellRange });

            numAtomsToSend_ += cellRange.size() * numAtomsPerCell_;
        }
    }

    if (debug)
    {
        fprintf(debug,
                "rank %d #atoms %d zone %d shift %d %d %d pbc %d %d %d sending %d atoms\n",
                dd.rank,
                dd.numHomeAtoms,
                zone_,
                domainShift_[0],
                domainShift_[1],
                domainShift_[2],
                pbcCoordinateShift_[0],
                pbcCoordinateShift_[1],
                pbcCoordinateShift_[2],
                numAtomsToSend_);
    }

    // Resize the communication buffer
    rvecBuffer_.resize(numAtomsToSend_);

    // Copy the global atom indices to the send buffer
    globalAtomIndices_.clear();
    for (const auto& columnInfo : columnsToSend_)
    {
        const int at_start = *columnInfo.cellRange.begin() * numAtomsPerCell_;
        const int at_end   = *columnInfo.cellRange.end() * numAtomsPerCell_;
        globalAtomIndices_.insert(globalAtomIndices_.end(),
                                  dd.globalAtomIndices.begin() + at_start,
                                  dd.globalAtomIndices.begin() + at_end);
    }
    GMX_ASSERT(int(globalAtomIndices_.size()) == numAtomsToSend_,
               "Here we should have all global atom indices to send");
}

FastVector<std::pair<int, int>> DomainCommBackward::makeColumnsSendBuffer() const
{
    // Store the grid info with only index and cell count per column
    FastVector<std::pair<int, int>> sendBuffer;
    sendBuffer.reserve(columnsToSend_.size());
    for (const auto& columnInfo : columnsToSend_)
    {
        sendBuffer.emplace_back(columnInfo.index, columnInfo.cellRange.size());
    }

    return sendBuffer;
}

DomainCommForward::DomainCommForward(int rank, int zone, MPI_Comm mpiCommAll) :
    rank_(rank), zone_(zone), mpiCommAll_(mpiCommAll)
{
}

void DomainCommForward::setup(const DomainCommBackward& send, const int offsetInCoordinateBuffer)
{
    std::array<int, 2> sendSizes = { int(send.columnsToSend().size()), send.numAtoms() };
    std::array<int, 2> receiveSizes;

    ddSendReceive(send, *this, dddirBackward, sendSizes.data(), 2, receiveSizes.data(), 2, HaloMpiTag::GridCounts);

    numAtomsToReceive_ = receiveSizes[1];

    atomRange_ = { offsetInCoordinateBuffer, offsetInCoordinateBuffer + numAtomsToReceive_ };

    if (debug)
    {
        fprintf(debug, "For zone %d, receiving %d atoms\n", zone_, numAtoms());
    }

    // Store the grid info with only index and cell count per column
    const FastVector<std::pair<int, int>> sendBuffer = send.makeColumnsSendBuffer();

    columnsReceived_.resize(receiveSizes[0]);

    ddSendReceive(send,
                  *this,
                  dddirBackward,
                  sendBuffer.data(),
                  send.columnsToSend().size(),
                  columnsReceived_.data(),
                  numAtomsToReceive_,
                  HaloMpiTag::GridColumns);
}

DomainPairComm::DomainPairComm(int         backwardRank,
                               int         forwardRank,
                               int         zone,
                               const IVec& domainShift,
                               PbcType     pbcType,
                               bool        commOverPbc,
                               IVec        pbcCoordinateShift,
                               MPI_Comm    mpiCommAll) :
    backward_(backwardRank, zone, domainShift, pbcType, commOverPbc, pbcCoordinateShift, mpiCommAll),
    forward_(forwardRank, zone, mpiCommAll)
{
}

} // namespace gmx
