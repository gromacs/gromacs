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
 * \brief Implements the DomdecZones class
 *
 * \inlibraryapi
 * \ingroup module_domdec
 *
 * \author Berk Hess <hess@kth.se>
 */

#include "domdec_zones.h"

#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/utility.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

//! The DD zone order
static const ivec sc_ddZoneOrder[sc_maxNumZones] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 },
                                                     { 0, 1, 0 }, { 0, 1, 1 }, { 0, 0, 1 },
                                                     { 1, 0, 1 }, { 1, 1, 1 } };

/*! \brief The non-bonded zone-pair setup for domain decomposition
 *
 * The first number is the i-zone, the second number the first j-zone seen by
 * this i-zone, the third number the last+1 j-zone seen by this i-zone.
 * As is, this is for 3D decomposition, where there are 4 i-zones.
 * With 2D decomposition use only the first 2 i-zones and a last+1 j-zone of 4.
 * With 1D decomposition use only the first i-zone and a last+1 j-zone of 2.
 */
static const int ddNonbondedZonePairRanges[sc_maxNumIZones][3] = { { 0, 0, 8 },
                                                                   { 1, 3, 6 },
                                                                   { 2, 5, 6 },
                                                                   { 3, 5, 7 } };

gmx::DomdecZones::DomdecZones(gmx::ArrayRef<const int> ddDims) :
    num_(1 << ddDims.size()), numIZones_(std::max(num_ / 2, 1))
{
    for (int iZone = 0; iZone < numIZones_; iZone++)
    {
        GMX_RELEASE_ASSERT(
                ddNonbondedZonePairRanges[iZone][0] == iZone,
                "The first element for each ddNonbondedZonePairRanges should match its index");

        jZoneRanges_[iZone] = gmx::Range<int>(std::min(ddNonbondedZonePairRanges[iZone][1], num_),
                                              std::min(ddNonbondedZonePairRanges[iZone][2], num_));
    }

    for (int i = 0; i < num_; i++)
    {
        int m      = 0;
        shifts_[i] = { 0, 0, 0 };
        for (int dim : ddDims)
        {
            shifts_[i][dim] = sc_ddZoneOrder[i][m++];
        }
    }
}

void gmx::DomdecZones::setSizes(const gmx_domdec_t&   dd,
                                const matrix          box,
                                const gmx_ddbox_t*    ddbox,
                                const gmx::Range<int> zoneRange)
{
    const gmx_domdec_comm_t& comm = *dd.comm;

    /* Do we need to determine extra distances for multi-body bondeds? */
    const bool checkMultiBodyDistances =
            (comm.systemInfo.haveInterDomainMultiBodyBondeds && isDlbOn(comm.dlbState) && dd.ndim > 1);

    for (int z : zoneRange)
    {
        /* Copy cell limits to zone limits.
         * Valid for non-DD dims and non-shifted dims.
         */
        sizes_[z].x0 = comm.cell_x0;
        sizes_[z].x1 = comm.cell_x1;
    }

    for (int d = 0; d < dd.ndim; d++)
    {
        const int dim = dd.dim[d];

        for (int z = 0; z < num_; z++)
        {
            /* With a staggered grid we have different sizes
             * for non-shifted dimensions.
             */
            if (isDlbOn(comm.dlbState) && shifts_[z][dim] == 0)
            {
                if (d == 1)
                {
                    sizes_[z].x0[dim] = comm.zone_d1[shifts_[z][dd.dim[d - 1]]].min0;
                    sizes_[z].x1[dim] = comm.zone_d1[shifts_[z][dd.dim[d - 1]]].max1;
                }
                else if (d == 2)
                {
                    sizes_[z].x0[dim] =
                            comm.zone_d2[shifts_[z][dd.dim[d - 2]]][shifts_[z][dd.dim[d - 1]]].min0;
                    sizes_[z].x1[dim] =
                            comm.zone_d2[shifts_[z][dd.dim[d - 2]]][shifts_[z][dd.dim[d - 1]]].max1;
                }
            }
        }

        real cutoffTwoBody   = comm.systemInfo.cutoff;
        real cutoffMultiBody = comm.cutoff_mbody;
        if (ddbox->tric_dir[dim])
        {
            cutoffTwoBody /= ddbox->skew_fac[dim];
            cutoffMultiBody /= ddbox->skew_fac[dim];
        }

        /* Set the lower limit for the shifted zone dimensions */
        for (int z : zoneRange)
        {
            if (shifts_[z][dim] > 0)
            {
                if (!isDlbOn(comm.dlbState) || d == 0)
                {
                    sizes_[z].x0[dim] = comm.cell_x1[dim];
                    sizes_[z].x1[dim] = comm.cell_x1[dim] + cutoffTwoBody;
                }
                else
                {
                    /* Here we take the lower limit of the zone from
                     * the lowest domain of the zone below.
                     */
                    if (z < 4)
                    {
                        sizes_[z].x0[dim] = comm.zone_d1[shifts_[z][dd.dim[d - 1]]].min1;
                    }
                    else
                    {
                        if (d == 1)
                        {
                            sizes_[z].x0[dim] = sizes_[zone_perm[2][z - 4]].x0[dim];
                        }
                        else
                        {
                            sizes_[z].x0[dim] =
                                    comm.zone_d2[shifts_[z][dd.dim[d - 2]]][shifts_[z][dd.dim[d - 1]]]
                                            .min1;
                        }
                    }
                    /* A temporary limit, is updated below */
                    sizes_[z].x1[dim] = sizes_[z].x0[dim];

                    if (checkMultiBodyDistances)
                    {
                        for (int zi = 0; zi < numIZones(); zi++)
                        {
                            if (shifts_[zi][dim] == 0)
                            {
                                /* This takes the whole zone into account.
                                 * With multiple pulses this will lead
                                 * to a larger zone then strictly necessary.
                                 */
                                sizes_[z].x1[dim] = std::max(sizes_[z].x1[dim],
                                                             sizes_[zi].x1[dim] + cutoffMultiBody);
                            }
                        }
                    }
                }
            }
        }

        /* Loop over the i-zones to set the upper limit of each
         * j-zone they see.
         */
        for (int iZone = 0; iZone < numIZones_; iZone++)
        {
            if (shifts_[iZone][dim] == 0)
            {
                /* We should only use zones up to zone_end */
                const auto& jZoneRangeFull = jZoneRanges_[iZone];
                if (*zoneRange.end() <= *jZoneRangeFull.begin())
                {
                    continue;
                }
                const gmx::Range<int> jZoneRange(*jZoneRangeFull.begin(),
                                                 std::min(*jZoneRangeFull.end(), *zoneRange.end()));
                for (int jZone : jZoneRange)
                {
                    if (shifts_[jZone][dim] > 0)
                    {
                        sizes_[jZone].x1[dim] =
                                std::max(sizes_[jZone].x1[dim], sizes_[iZone].x1[dim] + cutoffTwoBody);
                    }
                }
            }
        }
    }

    for (int z : zoneRange)
    {
        /* Initialization only required to keep the compiler happy */
        RVec corner_min = { 0, 0, 0 };
        RVec corner_max = { 0, 0, 0 };

        /* To determine the bounding box for a zone we need to find
         * the extreme corners of 4, 2 or 1 corners.
         */
        const int nc = 1 << (ddbox->nboundeddim - 1);

        for (int c = 0; c < nc; c++)
        {
            RVec corner;

            /* Set up a zone corner at x=0, ignoring trilinic couplings */
            corner[XX] = 0;
            if ((c & 1) == 0)
            {
                corner[YY] = sizes_[z].x0[YY];
            }
            else
            {
                corner[YY] = sizes_[z].x1[YY];
            }
            if ((c & 2) == 0)
            {
                corner[ZZ] = sizes_[z].x0[ZZ];
            }
            else
            {
                corner[ZZ] = sizes_[z].x1[ZZ];
            }
            if (dd.ndim == 1 && dd.dim[0] < ZZ && ZZ < dd.unitCellInfo.npbcdim
                && box[ZZ][1 - dd.dim[0]] != 0)
            {
                /* With 1D domain decomposition the atom groups are not in
                 * the triclinic box, but triclinic x-y and rectangular y/x-z.
                 * Shift the corner of the z-vector back to along the box
                 * vector of dimension d, so it will later end up at 0 along d.
                 * This can affect the location of this corner along dd.dim[0]
                 * through the matrix operation below if box[d][dd.dim[0]]!=0.
                 */
                int d = 1 - dd.dim[0];

                corner[d] -= corner[ZZ] * box[ZZ][d] / box[ZZ][ZZ];
            }
            /* Apply the triclinic couplings */
            for (int i = YY; i < ddbox->npbcdim && i < DIM; i++)
            {
                for (int j = XX; j < i; j++)
                {
                    corner[j] += corner[i] * box[i][j] / box[i][i];
                }
            }
            if (c == 0)
            {
                corner_min = corner;
                corner_max = corner;
            }
            else
            {
                corner_min = elementWiseMin(corner_min, corner);
                corner_max = elementWiseMax(corner_max, corner);
            }
        }
        /* Copy the extreme cornes without offset along x */
        sizes_[z].bb_x0 = corner_min;
        sizes_[z].bb_x1 = corner_max;

        /* Add the offset along x */
        sizes_[z].bb_x0[XX] += sizes_[z].x0[XX];
        sizes_[z].bb_x1[XX] += sizes_[z].x1[XX];
    }

    if (debug)
    {
        for (int z : zoneRange)
        {
            fprintf(debug,
                    "zone %d    %6.3f - %6.3f  %6.3f - %6.3f  %6.3f - %6.3f\n",
                    z,
                    sizes_[z].x0[XX],
                    sizes_[z].x1[XX],
                    sizes_[z].x0[YY],
                    sizes_[z].x1[YY],
                    sizes_[z].x0[ZZ],
                    sizes_[z].x1[ZZ]);
            fprintf(debug,
                    "zone %d bb %6.3f - %6.3f  %6.3f - %6.3f  %6.3f - %6.3f\n",
                    z,
                    sizes_[z].bb_x0[XX],
                    sizes_[z].bb_x1[XX],
                    sizes_[z].bb_x0[YY],
                    sizes_[z].bb_x1[YY],
                    sizes_[z].bb_x0[ZZ],
                    sizes_[z].bb_x1[ZZ]);
        }
    }
}

} // namespace gmx
