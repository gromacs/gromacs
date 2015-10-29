/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "domdec_halo.h"

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_grid.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/qsort_threadsafe.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_internal.h"
#include "domdec_utility.h"


/* MPI tags for non-blocking x and f communication.
 * We assign values that are unlikely to occur in other MPI calls.
 * This is to avoid incorrect MPI message matching, which could otherwise
 * occur with other MPI communication happening during the non-blocking
 * halo communcation initiate and complete calls.
 */
enum {
    halo_tag_X = 9990, halo_tag_F
};

typedef struct {
    real x0;
    real max_f0;
    real min_f1;
} zone_dim_bounds_t;


template <typename T>
void dd_sendrecv_zone(const gmx_domdec_t gmx_unused *dd,
                      int gmx_unused zone, int gmx_unused direction,
                      T gmx_unused *buf_s, int gmx_unused n_s,
                      T gmx_unused *buf_r, int gmx_unused n_r)
{
#if GMX_MPI
    int rank_s, rank_r;

    if (direction == dddirBackward)
    {
        rank_s = dd->comm->domain_backward[zone].rank;
        rank_r = dd->comm->domain_forward[zone].rank;
    }
    else
    {
        rank_s = dd->comm->domain_forward[zone].rank;
        rank_r = dd->comm->domain_backward[zone].rank;
    }

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(T), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(T), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, MPI_STATUS_IGNORE);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(T), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r, n_r*sizeof(T), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, MPI_STATUS_IGNORE);
    }
#endif // GMX_MPI
}

static void
dd_isend_zone_rvec(const gmx_domdec_t gmx_unused *dd,
                   const domain_comm_t gmx_unused *send, rvec gmx_unused *buf_s,
                   int gmx_unused tag,
                   zones_mpi_dir_t gmx_unused *zones_mpi_dir)
{
#if GMX_MPI
    if (send->dims.natoms > 0)
    {
        MPI_Isend(buf_s[0], send->dims.natoms*sizeof(rvec), MPI_BYTE,
                  send->rank, tag,
                  dd->mpi_comm_all,
                  &zones_mpi_dir->request[zones_mpi_dir->nrequest]);
        zones_mpi_dir->nrequest++;
    }
#endif
}

static void
dd_irecv_zone_recv(const gmx_domdec_t gmx_unused *dd,
                   const domain_comm_t gmx_unused *recv, rvec gmx_unused *buf_r,
                   int gmx_unused tag,
                   zones_mpi_dir_t gmx_unused *zones_mpi_dir)
{
#if GMX_MPI
    if (recv->dims.natoms > 0)
    {
        MPI_Irecv(buf_r[0], recv->dims.natoms*sizeof(rvec), MPI_BYTE,
                  recv->rank, tag,
                  dd->mpi_comm_all,
                  &zones_mpi_dir->request[zones_mpi_dir->nrequest]);
        zones_mpi_dir->nrequest++;
    }
#endif
}

static void
dd_recv_zone_rvec(const gmx_domdec_t gmx_unused *dd,
                  const domain_comm_t gmx_unused *recv, rvec gmx_unused *buf_r)
{
#if GMX_MPI
    if (recv->dims.natoms > 0)
    {
        MPI_Recv(buf_r[0], recv->dims.natoms*sizeof(rvec), MPI_BYTE,
                 recv->rank, 0,
                 dd->mpi_comm_all, MPI_STATUS_IGNORE);
    }
#endif
}

static void dd_halo_send_x_zone(gmx_domdec_t *dd, int zone,
                                const matrix box, rvec x[],
                                zones_mpi_dir_t *zones_mpi_dir)
{
    const domain_comm_t *send      = &dd->comm->domain_backward[zone];
    rvec                 shift_vec = { 0, 0, 0 };

    if (send->bPBC)
    {
        const ivec *shift = &dd->comm->zones.shift[zone];

        for (int d = 0; d < DIM; d++)
        {
            if ((*shift)[d] != 0 && dd->ci[d] == 0)
            {
                for (int e = 0; e < DIM; e++)
                {
                    shift_vec[e] += (*shift)[d]*box[d][e];
                }
            }
        }
    }

    rvec *sbuf = send->at_buf;

    int   j = 0;
    for (int c = 0; c < send->ncolumn; c++)
    {
        int start = send->column_atom_range[2*c];
        int end   = send->column_atom_range[2*c+1];

        if (!send->bPBC)
        {
            for (int i = start; i < end; i++)
            {
                copy_rvec(x[i], sbuf[j++]);
            }
        }
        else if (!send->bScrew)
        {
            for (int i = start; i < end; i++)
            {
                rvec_add(x[i], shift_vec, sbuf[j++]);
            }
        }
        else
        {
            /* Screw PBC */
            for (int i = start; i < end; i++)
            {
                /* Shift x */
                sbuf[j][XX] = x[i][XX] + shift_vec[XX];
                /* Rotate y and z.
                 * This operation requires a special shift force
                 * treatment, which is performed in calc_vir.
                 */
                sbuf[j][YY] = box[YY][YY] - x[i][YY];
                sbuf[j][ZZ] = box[YY][ZZ] - x[i][ZZ];
                j++;
            }
        }
    }

    assert(j == send->dims.natoms);

    dd_isend_zone_rvec(dd,
                       send, sbuf,
                       halo_tag_X,
                       zones_mpi_dir);
}

void dd_halo_initiate_recv_x(gmx_domdec_t *dd, rvec x[])
{
    zones_mpi_t *zmpi = &dd->comm->zones_mpi_x;

    assert(zmpi->recv.nrequest == 0);

    /* Post all the non-blocking receives */
    for (int zone = 1; zone < dd->comm->zones.n; zone++)
    {
        dd_irecv_zone_recv(dd,
                           &dd->comm->domain_forward[zone],
                           x + dd->comm->zones.at_range[zone],
                           halo_tag_X,
                           &zmpi->recv);
    }
}

void dd_halo_initiate_send_x(gmx_domdec_t *dd, matrix box, rvec x[])
{
    zones_mpi_dir_t *zmpid = &dd->comm->zones_mpi_x.send;

    assert(zmpid->nrequest == 0);

    for (int zone = 1; zone < dd->comm->zones.n; zone++)
    {
        /* Copy x to the send buffer and post non-blocking send */
        dd_halo_send_x_zone(dd, zone, box, x, zmpid);
    }
}

void dd_halo_complete_recv_x(gmx_domdec_t *dd)
{
    zones_mpi_t *zmpi = &dd->comm->zones_mpi_x;

#if GMX_MPI
    /* Wait for all non-blocking receives to complete */
    MPI_Waitall(zmpi->recv.nrequest, zmpi->recv.request, MPI_STATUSES_IGNORE);
#endif
    zmpi->recv.nrequest = 0;
}

void dd_halo_complete_send_x(gmx_domdec_t *dd)
{
    zones_mpi_t *zmpi = &dd->comm->zones_mpi_x;

#if GMX_MPI
    /* Wait for all non-blocking sends to complete */
    MPI_Waitall(zmpi->send.nrequest, zmpi->send.request, MPI_STATUSES_IGNORE);
#endif
    zmpi->send.nrequest = 0;
}

void dd_halo_move_x(gmx_domdec_t *dd, matrix box, rvec x[])
{
    dd_halo_initiate_recv_x(dd, x);
    dd_halo_initiate_send_x(dd, box, x);
    dd_halo_complete_recv_x(dd);
    dd_halo_complete_send_x(dd);
}

static void dd_halo_reduce_f_zone(gmx_domdec_t *dd, int zone,
                                  const rvec f_buf[],
                                  rvec f[], rvec fshift[])
{
    domain_comm_t *recv           = &dd->comm->domain_backward[zone];

    bool           bForcesNeedPbc = ((fshift != NULL && recv->bPBC) || recv->bScrew);

    int            j = 0;
    for (int c = 0; c < recv->ncolumn; c++)
    {
        int start, end;

        start = recv->column_atom_range[2*c];
        end   = recv->column_atom_range[2*c+1];

        if (!bForcesNeedPbc)
        {
            for (int i = start; i < end; i++)
            {
                rvec_inc(f[i], f_buf[j++]);
            }
        }
        else if (!recv->bScrew)
        {
            for (int i = start; i < end; i++)
            {
                rvec_inc(f[i], f_buf[j]);
                rvec_inc(fshift[recv->shift_ind], f_buf[j]);
                j++;
            }
        }
        else
        {
            /* Screw PBC */
            for (int i = start; i < end; i++)
            {
                /* Rotate the force */
                f[i][XX] += f_buf[j][XX];
                f[i][YY] -= f_buf[j][YY];
                f[i][ZZ] -= f_buf[j][ZZ];

                if (fshift != NULL)
                {
                    /* Add this force to the shift force */
                    rvec_inc(fshift[recv->shift_ind], f_buf[j]);
                }
                j++;
            }
        }
    }

    assert(j == recv->dims.natoms);
}

void dd_halo_initiate_recv_f(gmx_domdec_t *dd)
{
    zones_mpi_t *zmpi;

    zmpi = &dd->comm->zones_mpi_f;

    assert(zmpi->recv.nrequest == 0);

    /* Post all the non-blocking receives */
    int zone;

    for (zone = 1; zone < dd->comm->zones.n; zone++)
    {
        /* We can reuse the send x buffer as the receive buffer for f,
         * since the received x need to be processed before f can be
         * calculated and communicated.
         */
        dd_irecv_zone_recv(dd,
                           &dd->comm->domain_backward[zone],
                           dd->comm->domain_backward[zone].at_buf,
                           halo_tag_F,
                           &zmpi->recv);
    }
}

void dd_halo_initiate_send_f(gmx_domdec_t *dd, rvec f[])
{
    zones_mpi_t *zmpi = &dd->comm->zones_mpi_f;

    assert(zmpi->send.nrequest == 0);

    /* Non-blocking send using direct force buffer pointers */
    for (int zone = 1; zone < dd->comm->zones.n; zone++)
    {
        domain_comm_t *send = &dd->comm->domain_forward[zone];

        dd_isend_zone_rvec(dd,
                           send,
                           f + dd->comm->zones.at_range[zone],
                           halo_tag_F,
                           &zmpi->send);
    }
}

void dd_halo_complete_recv_f(gmx_domdec_t *dd, rvec f[], rvec fshift[])
{
    zones_mpi_t *zmpi = &dd->comm->zones_mpi_f;

#if GMX_MPI
    if (debug)
    {
        fprintf(debug, "Waiting for force recv: %d messages\n",
                zmpi->recv.nrequest);
    }

    /* Wait for all non-blocking communication to complete */
    MPI_Waitall(zmpi->recv.nrequest, zmpi->recv.request, MPI_STATUSES_IGNORE);
#endif
    zmpi->recv.nrequest = 0;

    /* Reduce the received non-local forces with our local forces */
    for (int zone = 1; zone < dd->comm->zones.n; zone++)
    {
        const rvec *f_buf = dd->comm->domain_backward[zone].at_buf;
        dd_halo_reduce_f_zone(dd, zone, f_buf, f, fshift);
    }
}

void dd_halo_complete_send_f(gmx_domdec_t *dd)
{
    zones_mpi_t *zmpi = &dd->comm->zones_mpi_f;

#if GMX_MPI
    if (debug)
    {
        fprintf(debug, "Waiting for force send: %d messages\n",
                zmpi->send.nrequest);
    }

    /* Wait for all non-blocking sends to complete */
    MPI_Waitall(zmpi->send.nrequest, zmpi->send.request, MPI_STATUSES_IGNORE);
#endif
    zmpi->send.nrequest = 0;
}

void dd_halo_move_f(gmx_domdec_t *dd, rvec f[], rvec fshift[])
{
    dd_halo_initiate_recv_f(dd);
    dd_halo_initiate_send_f(dd, f);
    dd_halo_complete_recv_f(dd, f, fshift);
    dd_halo_complete_send_f(dd);
}

void dd_halo_move_int(gmx_domdec_t *dd, int *buf)
{
    const gmx_domdec_zones_t *zones = &dd->comm->zones;

    for (int zone = 1; zone < zones->n; zone++)
    {
        domain_comm_t *send = &dd->comm->domain_backward[zone];
        domain_comm_t *recv = &dd->comm->domain_forward[zone];

        int            j = 0;
        for (int c = 0; c < send->ncolumn; c++)
        {
            int start = send->column_atom_range[2*c];
            int end   = send->column_atom_range[2*c+1];

            for (int i = start; i < end; i++)
            {
                send->atominfo[j++] = buf[i];
            }
        }

        /* Communicate the atom info flags, cheaper than extracting from mtop */
        dd_sendrecv_zone<int>(dd, zone, dddirBackward,
                              send->atominfo, send->dims.natoms,
                              buf + zones->at_range[zone], recv->dims.natoms);
    }
}

static gmx_inline real
corner_bb_distance2_rect(const ivec        zone_shift,
                         const real        corner[DIM],
                         const nbnxn_bb_t *bb)
{
    /* Here we use rounding and calculate the actual distance^2
     * to the corner(s).
     */
    real r2 = 0;

    /* Rectangular unit-cell, easy */
    for (int d = 0; d < DIM; d++)
    {
        if (zone_shift[d])
        {
            /* Distance to the corner for pair-wise interactions */
            real r = bb->lower[d] - corner[d];
            if (r > 0)
            {
                r2 += r*r;
            }
        }
    }

    return r2;
}

static gmx_inline real
corner_bb_distance2_tric(const ivec        zone_shift,
                         const rvec       *normal,
                         const ivec        sumSquares,
                         const real        corner[DIM],
                         const nbnxn_bb_t *bb)
{
    /* Here we use partial and approximate rounding */
    real r2 = 0;

    for (int d = 0; d < DIM; d++)
    {
        if (zone_shift[d])
        {
            real r = 0;
            for (int d2 = d; d2 < DIM; d2++)
            {
                r += ((normal[d][d2] >= 0 ? bb->lower[d2] : bb->upper[d2]) - corner[d2])*normal[d][d2];
            }
            if (r > 0)
            {
                if (sumSquares[d])
                {
                    /* The angle(s) between the normal of this zone plane
                     * an the preceding zone plane(s) is/are >= 90 degrees.
                     * Add up the squares of the distance. This underestimates
                     * the distance for angles > 90 degrees.
                     */
                    r2 += r*r;
                }
                else
                {
                    /* The/A angle between the normals is < 90 degrees.
                     * We use the maximum distance, which is an underestimate.
                     */
                    r2 = std::max(r2, r*r);
                }
            }
        }
    }

    return r2;
}

/* Wrapper function for corner - bounding-box distance calculation.
 * Only splits triclinic vs non-triclinic distance calculations.
 */
static gmx_inline real
corner_bb_distance2(const ivec        zone_shift,
                    gmx_bool          bTriclinic,
                    const rvec       *normal,
                    const ivec        sumSquares,
                    const real        corner[DIM],
                    const nbnxn_bb_t *bb)
{
    if (bTriclinic)
    {
        return corner_bb_distance2_tric(zone_shift, normal, sumSquares,
                                        corner, bb);
    }
    else
    {
        return corner_bb_distance2_rect(zone_shift, corner, bb);
    }
}

/* Set the cell count along x and y in zone_comm.
 * Ensures enough memory is present for things that depend on the cell count.
 */
static void set_domain_comm_ncolumns(domain_comm_t *domain_comm,
                                     int ncolumn_x, int ncolumn_y)
{
    domain_comm->dims.ncolumn_x = ncolumn_x;
    domain_comm->dims.ncolumn_y = ncolumn_y;
    if (domain_comm->dims.ncolumn_x*domain_comm->dims.ncolumn_y > domain_comm->column_nalloc)
    {
        domain_comm->column_nalloc = over_alloc_dd(domain_comm->dims.ncolumn_x*domain_comm->dims.ncolumn_y);
        srenew(domain_comm->column_natoms, domain_comm->column_nalloc);
        /* We usually need less space, but this is negligible on the total */
        srenew(domain_comm->column_atom_range, domain_comm->column_nalloc*2);
    }
}

/* Determine the corner for 2-body and for multi-body interaction atom
 * communication distance calculations for combination of i-zones for the
 * remote rank that sees our zone as index 'zone' with shift 'zone_shift'.
 * *bCornersDiffers tells if corner_2body and corner_multibody are different.
 *
 * When there is a single i-zone or multiple i-zones without staggering (DLB),
 * both corners are identical to the upper corner (union of) the i-zone(s).
 * When there are multiple i-zones with staggering, the minimum communication
 * volume needs to be calculated using the minimum distane to multiple corners.
 * But since this is complex, we compute a single corner that leads to
 * a communication volume that ecompasses the volumes for all i-zones. This
 * leads to a bit of a lot of extra, useless communication.
 */
static void get_zone_corners(const gmx_domdec_t *dd, const matrix box,
                             int zone, const ivec zone_shift,
                             const zone_dim_bounds_t *zone_bounds_dim1,
                             const zone_dim_bounds_t *zone_bounds_dim2,
                             rvec corner_2body, rvec corner_multibody,
                             gmx_bool *bCornersDiffer)
{
    int nizone = 0;
    int izone[4];

    if (dd->comm->dlbState == edlbsOn)
    {
        /* Make a list of i-zones that see our zone as a j-zone */
        nizone = 0;
        for (int z = 0; z < (1 << (dd->ndim - 1)); z++)
        {
            if (zone >= ddNonbondedZonePairRanges[z][1] &&
                zone <  ddNonbondedZonePairRanges[z][2])
            {
                izone[nizone++] = ddNonbondedZonePairRanges[z][0];
            }
        }
    }

    clear_rvec(corner_2body);
    clear_rvec(corner_multibody);

    *bCornersDiffer = FALSE;

    for (int dim = 0; dim < DIM; dim++)
    {
        real corner_2body_dim, corner_multibody_dim;

        if (zone_shift[dim] == 0 ||
            dim == dd->dim[0] ||
            dd->comm->dlbState != edlbsOn)
        {
            /* No staggering, all bounds are equal to our local bounds */
            corner_2body_dim     = dd->comm->cell_x0[dim];
            corner_multibody_dim = dd->comm->cell_x0[dim];
        }
        else if (dim == dd->dim[1])
        {
            /* Multi-body bonded interactions can involve all zones,
             * so we need to use the maximum over all corners.
             */
            int recv_ind0        = 1 - ddZoneOrder[zone][0];
            corner_multibody_dim = std::max(zone_bounds_dim1[recv_ind0    ].x0,
                                            zone_bounds_dim1[recv_ind0 + 1].x0);

            /* Determine the maximum corner over the i-zones for our zone */
            corner_2body_dim = 0;
            for (int z = 0; z < nizone; z++)
            {
                corner_2body_dim = std::max(corner_2body_dim,
                                            zone_bounds_dim1[recv_ind0 + ddZoneOrder[izone[z]][0]].x0);
            }
        }
        else
        {
            assert(dim == dd->dim[2]);

            int recv_ind0        = 1 - ddZoneOrder[zone][0];
            int recv_ind1        = 1 - ddZoneOrder[zone][1];
            corner_multibody_dim = 0;
            for (int sd0 = 0; sd0 < 2; sd0++)
            {
                for (int sd1 = 0; sd1 < 2; sd1++)
                {
                    corner_multibody_dim = std::max(corner_multibody_dim,
                                                    zone_bounds_dim2[(recv_ind0 + sd0)*3 + recv_ind1 + sd1].x0);
                }
            }

            corner_2body_dim = 0;
            for (int z = 0; z < nizone; z++)
            {
                corner_2body_dim = std::max(corner_2body_dim,
                                            zone_bounds_dim2[(recv_ind0 + ddZoneOrder[izone[z]][0])*3 + recv_ind1 + ddZoneOrder[izone[z]][1]].x0);
            }
        }

        if (corner_multibody_dim != corner_2body_dim)
        {
            *bCornersDiffer = TRUE;
        }

        corner_2body[dim]     = corner_2body_dim;
        corner_multibody[dim] = corner_multibody_dim;
    }

    dd_zone_bounds_to_corner(dd, box, corner_2body);
    dd_zone_bounds_to_corner(dd, box, corner_multibody);

    if (debug)
    {
        fprintf(debug, "rank %3d halo zone %d corner 2b %6.3f %6.3f %6.3f mb %6.3f %6.3f %6.3f\n",
                dd->rank, zone,
                corner_2body[XX], corner_2body[YY], corner_2body[ZZ],
                corner_multibody[XX], corner_multibody[YY], corner_multibody[ZZ] );
    }
}

/* Sets sumSquares that tells if to use rounding for corners of zones */
static void get_sumSquares(const ivec  zone_shift,
                           const rvec *normal,
                           ivec        sumSquares)
{
    /* For triclinic zones with shifts along multiple unit-cell vectors,
     * the exact distance calculation gets very complex, since the normals
     * to the zone planes are not orthogonal. This makes rounding the edges
     * of those zones hard.
     * We apply approximate rounding in case the angles between the normals
     * are >= 90 degrees and no rouding when < 90 degrees. This leads to
     * a few more atoms communicated than strictly necessary, but results
     * in relatively simple, efficient and less bug-prone code.
     */
    int  dim_count = 0;
    int  dim_prev  = 0;
    for (int dim = 0; dim < DIM; dim++)
    {
        if (zone_shift[dim])
        {
            switch (dim_count)
            {
                case 0:
                    /* First dimension, doesn't matter what we do here */
                    sumSquares[dim] = 1;
                    break;
                case 1:
                    /* Second dimension, determine the angle between normals;
                     * if angle >= 90 degrees: sum squares of distances, rounds
                     *                         the zone, but not optimally, as
                     *                         the distance is underestimated
                     * if angle < 90 degrees:  use the maximum, no rounding,
                     *                         underestimates the distance.
                     */
                    sumSquares[dim] = (iprod(normal[dim],
                                             normal[dim_prev]) <= 0);
                    break;
                case 2:
                    /* As case 1, but check the angels with both planes */
                    assert(dim == ZZ);
                    /* The normal along z is always (0,0,1) */
                    sumSquares[dim] = (normal[0][ZZ] <= 0 &&
                                       normal[1][ZZ] <= 0);
                    break;
            }
            dim_prev = dim;
            dim_count++;
        }
    }
}

/* Determine which nbnxn grid cells (and atoms) we need to send for this zone */
static void
setup_domain_comm(gmx_domdec_t            *dd,
                  int                      zone,
                  const nbnxn_search_t     nbs,
                  real                     rc2_nonbonded,
                  real                     rc2_bonded,
                  const matrix             box,
                  const ivec               tric_dir,
                  const rvec              *normal,
                  const rvec               skew_fac,
                  const zone_dim_bounds_t *zone_bounds_dim1,
                  const zone_dim_bounds_t *zone_bounds_dim2,
                  const gmx_bool          *cell_missing_link,
                  const int               *atinfo,
                  domain_comm_t           *send)
{
    gmx_domdec_comm_t *comm       = dd->comm;
    gmx_bool           bTriclinic = FALSE;
    gmx_bool           bBondComm  = comm->bBondComm;

    /* Convert the shift array indexed by dim_ind to indexed by dim */
    ivec zone_shift = { 0, 0, 0 };
    for (int dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
    {
        if (ddZoneOrder[zone][dim_ind])
        {
            int dim = dd->dim[dim_ind];

            zone_shift[dim] = 1;

            if (tric_dir[dim])
            {
                bTriclinic = TRUE;
            }
        }
    }

    /* Get if we should use rounding of corners */
    ivec sumSquares;

    get_sumSquares(zone_shift, normal, sumSquares);

    /* Get the corners for non-bonded and bonded distance calculations */
    rvec     corner_2body, corner_multibody;
    gmx_bool bCornersDiffer;
    get_zone_corners(dd, box,
                     zone, zone_shift, zone_bounds_dim1, zone_bounds_dim2,
                     corner_2body, corner_multibody,
                     &bCornersDiffer);

    /* Do we need to determine extra distances for multi-body bondeds?
     * Note that with bBondComm we might need distances longer than
     * the non-bonded cut-off, but with a grid without staggering (bGridJump)
     * this check is indentical to the one triggered by bDist2Body below.
     */
    gmx_bool bDistMultiBody = (comm->bInterCGMultiBody && bCornersDiffer);

    /* Do we need to determine extra distances for only two-body bondeds? */
    gmx_bool bDist2Body     = (bBondComm && !bDistMultiBody);

    if (debug)
    {
        fprintf(debug, "setup zone %d BondComm %d DistMultiBody %d Dist2Body %d\n",
                zone, bBondComm, bDistMultiBody, bDist2Body);
    }

    nbnxn_get_local_grid_dimensions(nbs, &send->dims);

    int ncolumn_x        = send->dims.ncolumn_x;
    int ncolumn_y        = send->dims.ncolumn_y;
    set_domain_comm_ncolumns(send, ncolumn_x, ncolumn_y);
    /* send->ncolumn_x will be updated when adding cells to send */
    send->dims.ncolumn_x = 0;
    send->ncolumn        = 0;
    send->dims.natoms    = 0;
    int   column_y_max   = 0;
    float bb_z_max       = 0;

    for (int column = 0; column < ncolumn_x*ncolumn_y; column++)
    {
        nbnxn_bb_t        bb_full;
        const nbnxn_bb_t *bb_ptr = NULL;
        gmx_bool          bInRange;

        int               column_x = column/ncolumn_y;
        int               column_y = column - column_x*ncolumn_y;

        /* Get all the information for the grid column */
        nbnxn_grid_column_t col;
        nbnxn_get_local_grid_column(nbs, column_x, column_y, &col);

        /* Set this colum empty, for now (done here because of the continue) */
        send->column_natoms[column] = 0;

        if (col.nbb == 0)
        {
            /* Empty column */
            continue;
        }

        real r2  = corner_bb_distance2(zone_shift, bTriclinic,
                                       normal, sumSquares,
                                       corner_2body, &col.column_bb);
        real rb2 = 0;
        if (bDistMultiBody)
        {
            rb2  = corner_bb_distance2(zone_shift, bTriclinic,
                                       normal, sumSquares,
                                       corner_multibody, &col.column_bb);
        }

        if (!(r2 < rc2_nonbonded ||
              ((bDistMultiBody && rb2 < rc2_bonded) ||
               (bDist2Body     && r2  < rc2_bonded))))
        {
            /* This whole column is out of range */
            continue;
        }

        if (col.bb_z != NULL)
        {
            /* Take the x & y components from the column */
            bb_full = col.column_bb;
            bb_ptr  = &bb_full;
        }

        int bb_first = 0;
        int bb_last  = col.nbb - 1;
        do
        {
            if (col.bb != NULL)
            {
                bb_ptr = &col.bb[bb_last];
            }
            else
            {
                assert(col.bb_z != NULL);
                bb_full.lower[ZZ] = col.bb_z[bb_last*2 + 0];
                bb_full.upper[ZZ] = col.bb_z[bb_last*2 + 1];
            }

            r2      = corner_bb_distance2(zone_shift, bTriclinic,
                                          normal, sumSquares,
                                          corner_2body, bb_ptr);
            if (bDistMultiBody)
            {
                rb2 = corner_bb_distance2(zone_shift, bTriclinic,
                                          normal, sumSquares,
                                          corner_multibody, bb_ptr);
            }

            /* Here we check:
             * the 2-atom distance against the non-bonded cut-off,
             * the multi-body distance against the bonded cut-off
             * The 2-atom distance against the bonded cut-off.
             * The bonded check only triggers communication without bBondComm
             * or when the cell has missing bonded interactions.
             */
            bInRange = (r2 < rc2_nonbonded ||
                        (((bDistMultiBody && rb2 < rc2_bonded) ||
                          (bDist2Body     && r2  < rc2_bonded)) &&
                         (!bBondComm || cell_missing_link[col.bb_start + bb_last])));

            if (!bInRange)
            {
                bb_last--;
            }
        }
        while (bb_last >= bb_first && !bInRange);

        /* For rectangular grids without bondcomm we do not try to eliminate
         * cells from the bottom, since if bb_last is within range,
         * there is a high chance that bb_first=0 is also in range.
         */
        if ((bTriclinic || bBondComm) && bb_last > bb_first)
        {
            /* This loop is a copy of the one above with bb_last replaced
             * by bb_first. Putting it in a function would be cleaner,
             * but this would require 20 parameters.
             */
            do
            {
                if (col.bb != NULL)
                {
                    bb_ptr = &col.bb[bb_first];
                }
                else
                {
                    bb_full.lower[ZZ] = col.bb_z[bb_first*2 + 0];
                    bb_full.upper[ZZ] = col.bb_z[bb_first*2 + 1];
                }

                real r2  = corner_bb_distance2(zone_shift, bTriclinic,
                                               normal, sumSquares,
                                               corner_2body, bb_ptr);
                real rb2 = 0;
                if (bDistMultiBody)
                {
                    rb2  = corner_bb_distance2(zone_shift, bTriclinic,
                                               normal, sumSquares,
                                               corner_multibody, bb_ptr);
                }

                bInRange = (r2 < rc2_nonbonded ||
                            (((bDistMultiBody && rb2 < rc2_bonded) ||
                              (bDist2Body     && r2  < rc2_bonded)) &&
                             (!bBondComm || cell_missing_link[col.bb_start + bb_first])));

                if (!bInRange)
                {
                    bb_first++;
                }
            }
            while (bb_first < bb_last && !bInRange);
        }

        /* Check if we have cells within range */
        if (bb_last >= bb_first)
        {
            /* Determine the range of real atoms to send.
             * The filler atoms are always at the end, so we use min() there.
             */
            int at_start = col.atom_start + bb_first*col.bb_natoms;
            int at_end   = col.atom_start + std::min((bb_last + 1)*col.bb_natoms,
                                                     col.natoms);

            /* Set the grid column atom count in the send buffer */
            send->column_natoms[column] = at_end - at_start;

            /* Store the atom range, needed for preparing data for sending */
            send->column_atom_range[send->ncolumn*2  ] = at_start;
            send->column_atom_range[send->ncolumn*2+1] = at_end;
            send->ncolumn++;

            send->dims.natoms  += at_end - at_start;

            /* We need to send the grid up to and including x-column cx */
            send->dims.ncolumn_x = column_x + 1;
            /* Also keep track of the maximum range in y and z */
            column_y_max         = std::max(column_y_max, column_y);
            if (col.bb != NULL)
            {
                bb_z_max         = std::max(bb_z_max, col.bb[bb_last].upper[ZZ]);
            }
            else
            {
                bb_z_max         = std::max(bb_z_max, col.bb_z[bb_last*2 + 1]);
            }
        }
    }

    clear_ivec(send->rCheckBonded);
    if (comm->bInterCGBondeds)
    {
        /* Check if we need to check distances for bonded assigments
         * on the receiver side.
         */
        for (int dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
        {
            int dim = dd->dim[dim_ind];
            if (dim < dd->npbcdim && zone_shift[dim] && dd->nc[dim] == 2)
            {
                /* With only 2 domains, the same pair (or triplet, ...)
                 * of atoms can potentially be present on multiple ranks.
                 * We need to check if we need to check atom distances when
                 * assigning bonded interactions on the receiving rank.
                 */
                if (tric_dir[dim] || (comm->dlbState == edlbsOn && dim_ind > 0))
                {
                    /* For triclinic or staggered dimensions determining
                     * the range of communicated atoms is tedious and
                     * bug-prone, so we always check.
                     */
                    send->rCheckBonded[dim] = TRUE;
                }
                else
                {
                    real cellThickness =
                        (comm->cell_x1[dim] - comm->cell_x0[dim])*skew_fac[dim];

                    /* To avoid checking distances when assigning bondeds,
                     * a slab of the cell with a thickness of at least
                     * the maximum allowed distance between atoms in bondeds
                     * should not have any atoms communicated.
                     */
                    real rCheckBondedMargin =
                        cellThickness - std::max(comm->cutoff,
                                                 comm->cutoff_mbody);

                    real commDist = 0;
                    switch (dim)
                    {
                        case XX:
                            commDist = send->dims.ncolumn_x*send->dims.size_x;
                            break;
                        case YY:
                            commDist = (column_y_max + 1)*send->dims.size_y;
                            break;
                        case ZZ:
                            /* No staggering here, so all corners are equal */
                            commDist = bb_z_max - corner_2body[ZZ];
                            break;
                    }

                    if (commDist > rCheckBondedMargin)
                    {
                        send->rCheckBonded[dim] = TRUE;
                    }
                }
            }
        }
    }

    /* Ensure the atom buffers have enough space */
    if (send->dims.natoms > send->at_nalloc)
    {
        send->at_nalloc = over_alloc_dd(send->dims.natoms);
        srenew(send->index_gl, send->at_nalloc);
        srenew(send->atominfo, send->at_nalloc);
        /* Reallocate the x/f buffer for MPI_Isend and MPI_Irecv */
        srenew(send->at_buf, send->at_nalloc);
    }

    /* Copy the global atom indices to the send buffer */

    int ind = 0;
    for (int column = 0; column < send->ncolumn; column++)
    {
        int at_start = send->column_atom_range[column*2  ];
        int at_end   = send->column_atom_range[column*2+1];
        /* Copy the global atom indices to the send buffer */
        for (int at = at_start; at < at_end; at++)
        {
            send->index_gl[ind] = dd->index_gl[at];
            send->atominfo[ind] = atinfo[at];
            ind++;
        }
    }
}

/* Communicate the halo communication setup */
static void
comm_domain_comm_setup(gmx_domdec_t  *dd,
                       int            zone,
                       const matrix   box,
                       domain_comm_t *send,
                       domain_comm_t *recv)
{
    /* Communicate the grid dimensions */
    copy_ivec(send->rCheckBonded, send->dims.work);
    dd_sendrecv_zone<nbnxn_grid_dims_t>(dd, zone, dddirBackward,
                                        &send->dims, 1,
                                        &recv->dims, 1);
    copy_ivec(recv->dims.work, recv->rCheckBonded);

    if (debug)
    {
        fprintf(debug, "For zone %d, receiving %d atoms\n",
                zone, recv->dims.natoms);
    }

    /* Apply PBC on the receiver side */
    for (int dim_ind = 0; dim_ind < DIM; dim_ind++)
    {
        int dim = dd->dim[dim_ind];

        if (dim < dd->npbcdim &&
            ddZoneOrder[zone][dim_ind] == 1 &&
            dd->ci[dim] == dd->nc[dim] - 1)
        {
            /* We are communicating over pbc along dim */
            rvec_inc(recv->dims.corner0, box[dim]);
            rvec_inc(recv->dims.corner1, box[dim]);
        }
    }

    set_domain_comm_ncolumns(recv, recv->dims.ncolumn_x, recv->dims.ncolumn_y);

    dd_sendrecv_zone<int>(dd, zone, dddirBackward,
                          send->column_natoms, send->dims.ncolumn_x*send->dims.ncolumn_y,
                          recv->column_natoms, recv->dims.ncolumn_x*recv->dims.ncolumn_y);
}

/* Gathers the zone bounds for backward and forward neighboring cells
 * along decomposition dimension index dim_ind
 * over decomposition dimension index < dim_ind.
 * Stores the bounds for dim_ind=1 in b of size 3.
 * Stores the bounds for dim_ind=2 in b of size 3x3=9, major index is dim_ind=0.
 * Note that the extremes are only valid for the central entry.
 */
static void dd_halo_move_zone_bounds(gmx_domdec_t *dd, int dim_ind,
                                     zone_dim_bounds_t *b)
{
    gmx_domdec_comm_t *comm;
    int                dim, central, cstart, cnr, dic, ind;
    const int          zdb_size = sizeof(zone_dim_bounds_t)/sizeof(real);

    comm = dd->comm;

    dim  = dd->dim[dim_ind];

    /* Determine the index of the central, our home, zone.
     * For dim_ind=1 we have a total of 3 bounds, center=1,
     * for dim_ind=2 we have a total of 3*3=9 bounds, center=4;
     */
    central           = (dim_ind - 1)*3 + 1;

    b[central].x0     = comm->cell_x0[dim];
    b[central].max_f0 = comm->cell_f0[dim_ind];
    b[central].min_f1 = comm->cell_f1[dim_ind];

    /* We start communicating from the central element, sending 1 element */
    cstart = central;
    cnr    = 1;
    /* Loop over the 1 or 2 dimensions for communicating the bounds */
    for (dic = dim_ind - 1; dic >= 0; dic--)
    {
        dd_sendrecv_real(dd, dic, dddirForward,
                         (real *)&b[cstart],       cnr*zdb_size,
                         (real *)&b[cstart - cnr], cnr*zdb_size);

        if (dd->nc[dd->dim[dic]] > 2)
        {
            dd_sendrecv_real(dd, dic, dddirBackward,
                             (real *)&b[cstart],       cnr*zdb_size,
                             (real *)&b[cstart + cnr], cnr*zdb_size);
        }
        else
        {
            /* We have only 2 domain cells along this dimension,
             * we save some communication by copying instead.
             */
            int i;

            for (i = 0; i < cnr; i++)
            {
                b[cstart + cnr + i] = b[cstart - cnr + i];
            }
        }

        /* Loop over the central elements of the two received entries/rows */
        for (ind = central - cnr; ind < central + 2*cnr; ind += 2*cnr)
        {
            b[central].max_f0 = std::max(b[central].max_f0, b[ind].max_f0);
            b[central].min_f1 = std::min(b[central].min_f1, b[ind].min_f1);
        }
        /* Note that the max and min values in entries 0 and 2 are not up to date */

        /* We now have now added cnr elements before as well as after the cnr
         * elements we started from: shift the start and increase the count.
         */
        cstart -= cnr;
        cnr    *= 3;
    }
    if (debug)
    {
        int i;

        fprintf(debug, "Bounds, di=%d:", dim_ind);
        for (i = 0; i < cnr; i++)
        {
            fprintf(debug, " %5.2f", b[i].x0);
        }
        fprintf(debug, "\n");
    }

    /* Store the relative extremes for use during the next DLB step */
    comm->cell_f_max0[dim_ind] = b[central].max_f0;
    comm->cell_f_min1[dim_ind] = b[central].min_f1;
}

/* Check if atom with global index at_gl has any bonded interactions
 * that involve non-local atoms.
 */
static gmx_bool gmx_inline missing_link(const t_blocka *link,
                                        int at_gl,
                                        const gmx_ga2la_t *ga2la)
{
    int i;

    for (i = link->index[at_gl]; i < link->index[at_gl+1]; i++)
    {
        if (!ga2la_is_home(ga2la, link->a[i]))
        {
            return TRUE;
        }
    }

    return FALSE;
}

/* Allocate, generate and return *cell_missing_link_ptr, a list of of booleans
 * for a grid cells, which tell if a cell has non-local bonded interactions.
 */
static gmx_bool *mark_nbnxn_cells_bondcomm(const gmx_domdec_t  *dd,
                                           const int           *atinfo,
                                           const nbnxn_search_t nbs)
{
    int bb_natoms;
    {
        nbnxn_grid_column_t col;

        /* We only need bb_natoms, but we use the same function call as
         * used later in the communication setup, so we don't need to write
         * an extra function and avoid more bug opportunities.
         */
        nbnxn_get_local_grid_column(nbs, 0, 0, &col);
        bb_natoms = col.bb_natoms;
    }
    const int *atom_order;
    int        na;
    nbnxn_get_atomorder(nbs, &atom_order, &na);

    int nc = na/bb_natoms;
    assert(nc*bb_natoms == na);

    const t_blocka *link  = dd->comm->cglink;
    gmx_bool       *cell_missing_link;
    snew(cell_missing_link, nc);

    /* Loop over all cells of the local grid */
    for (int c = 0; c < nc; c++)
    {
        /* Loop over all atoms in this cell */
        gmx_bool bMissingLink = FALSE;
        for (int i = c*bb_natoms; i < (c + 1)*bb_natoms && !bMissingLink; i++)
        {
            int a = atom_order[i];

            /* Check if this is a real atom (not a filler atom) */
            if (a >= 0 && GET_CGINFO_BOND_INTER(atinfo[a]))
            {
                if (missing_link(link, dd->index_gl[a], dd->ga2la))
                {
                    bMissingLink = TRUE;
                }
            }
        }
        cell_missing_link[c] = bMissingLink;
    }

    return cell_missing_link;
}

void setup_halo_communication(gmx_domdec_t *dd,
                              const matrix box, const gmx_ddbox_t *ddbox,
                              t_forcerec *fr,
                              gmx_bool bCellsChanged)
{
    if (debug)
    {
        fprintf(debug, "Setting up DD communication\n");
    }

    assert(fr->cutoff_scheme == ecutsVERLET);

    gmx_domdec_comm_t *comm  = dd->comm;
    zone_dim_bounds_t  zone_bounds_dim1[3], zone_bounds_dim2[9];

    if (dd->comm->dlbState == edlbsOn && dd->ndim > 1 && bCellsChanged)
    {
        dd_halo_move_zone_bounds(dd, 1, zone_bounds_dim1);

        if (dd->ndim == 3)
        {
            /* We could save a few communication calls by integrating
             * the communication for dim_ind 1 and 2, but that leads
             * to far more complex code.
             */
            dd_halo_move_zone_bounds(dd, 2, zone_bounds_dim2);
        }
    }

    /* Determine the normals of the zone planes */
    rvec normal[DIM];
    for (int i = 0; i < DIM; i++)
    {
        if (ddbox->tric_dir[i])
        {
            /* ddbox->normal has length skew_fac, normalize it */
            svmul(1/ddbox->skew_fac[i], ddbox->normal[i], normal[i]);
        }
        else
        {
            clear_rvec(normal[i]);
            normal[i][i] = 1;
        }
    }

    gmx_bool *cell_missing_link = NULL;
    if (comm->bBondComm && dd->nnodes > 1)
    {
        /* Create a boolean array for cells telling if bondeds linked to atoms
         * are not locally present, so we need to communicate those cells.
         */
        cell_missing_link =
            mark_nbnxn_cells_bondcomm(dd, fr->cginfo, fr->nbv->nbs);
    }

    comm->zones.at_range[0] = 0;
    comm->zones.at_range[1] = dd->nat_home;
    dd->nat_tot             = dd->nat_home;
    dd->ncg_tot             = dd->nat_home;

    clear_ivec(comm->rCheckBonded);

    gmx_domdec_zones_t *zones         = &comm->zones;

    /* The naming is not fully consistent here (yet).
     * But the 2-body bonded interaction cut-off is max(cutoff, cutoff_mbody),
     * so there is not short and exact naming.
     */
    real rc2_nonbonded = gmx::square(comm->cutoff);
    real rc2_bonded    = gmx::square(comm->cutoff_mbody);

    /* Here we distribute the comm setup calculation over threads.
     * Note that setup_zone_comm is actually not very time consuming.
     * Alternatively we could overlap setup_zone_comm with comm_zone_comm_setup
     * for the previous zone.
     */
    int nthread = gmx_omp_nthreads_get(emntDomdec);
#pragma omp parallel for num_threads(nthread) schedule(static, 1)
    for (int zone_ind = 1; zone_ind < zones->n; zone_ind++)
    {
        int zone;

        if (nthread == 1)
        {
            zone = zone_ind;
        }
        else
        {
            /* When ordering the zones Cartesian, we get better load balancing
             * with 2 and 4 threads.
             */
            zone = zone_reorder_cartesian[zone_ind];
        }

        domain_comm_t *send = &comm->domain_backward[zone];

        /* Get the cg's for this pulse in this zone/domain */
        setup_domain_comm(dd, zone,
                          fr->nbv->nbs,
                          rc2_nonbonded, rc2_bonded,
                          box, ddbox->tric_dir,
                          normal, ddbox->skew_fac,
                          zone_bounds_dim1,
                          zone_bounds_dim2,
                          cell_missing_link,
                          fr->cginfo,
                          send);
    }

    if (cell_missing_link != NULL)
    {
        sfree(cell_missing_link);
    }

    for (int zone = 1; zone < zones->n; zone++)
    {
        domain_comm_t *send = &comm->domain_backward[zone];
        domain_comm_t *recv = &comm->domain_forward[zone];

        /* This communication could be overlapped with setup_zone_comm
         * for the next zone.
         */
        comm_domain_comm_setup(dd, zone, box, send, recv);

        /* Set the atom ranges for the zones */
        zones->at_range[zone+1] = zones->at_range[zone] + recv->dims.natoms;

        dd->nat_tot += recv->dims.natoms;
        dd->ncg_tot += recv->dims.natoms;
    }

    /* Ensure that we have space for all the atoms in the halo */
    if (dd->nat_tot > dd->cg_nalloc)
    {
        dd->cg_nalloc = over_alloc_dd(dd->nat_tot);
        srenew(dd->index_gl, dd->cg_nalloc);
        srenew(dd->cgindex, dd->cg_nalloc+1);
    }
    if (dd->nat_tot > fr->cg_nalloc)
    {
        fr->cg_nalloc = over_alloc_dd(dd->nat_tot);
        srenew(fr->cginfo, fr->cg_nalloc);
    }

    for (int zone = 1; zone < zones->n; zone++)
    {
        domain_comm_t *send = &comm->domain_backward[zone];
        domain_comm_t *recv = &comm->domain_forward[zone];

        /* NOTE: We should consider replacing the index_gl array by an array
         *       that combines index_gl and atominfo. This not only saves
         *       MPI calls here, but since they are mostly used together,
         *       it might also help caching. We should then copy atominfo
         *       fr->cginfo or only use the new combined array.
         *       This should be done when the group scheme is removed,
         *       since currently some code use both cg and atom indices.
         */

        /* Communicate the global atom indices */
        dd_sendrecv_zone<int>(dd, zone, dddirBackward,
                              send->index_gl, send->dims.natoms,
                              dd->index_gl + zones->at_range[zone], recv->dims.natoms);

        /* Communicate the atom info flags, cheaper than extracting from mtop */
        dd_sendrecv_zone<int>(dd, zone, dddirBackward,
                              send->atominfo, send->dims.natoms,
                              fr->cginfo + zones->at_range[zone], recv->dims.natoms);
        /* Set up the nbnxn grid for this non-local zone/domain
         * using the grid setup of this domain from the home rank.
         */
        nbnxn_set_domain_grid(fr->nbv->nbs, zone,
                              &recv->dims, recv->column_natoms);

        for (int i = 0; i < DIM; i++)
        {
            if (recv->rCheckBonded[i])
            {
                comm->rCheckBonded[i] = 1;
            }
        }
    }

    comm->nat[ddnatHOME] = dd->nat_home;
    for (int i = ddnatZONE; i < ddnatNR; i++)
    {
        comm->nat[i] = dd->nat_tot;
    }

    if (debug)
    {
        fprintf(debug, "Finished setting up DD communication, zones:");
        for (int c = 1; c < zones->n; c++)
        {
            fprintf(debug, " %d", comm->domain_forward[c].dims.natoms);
        }
        fprintf(debug, "\n");
    }

    /* TODO: remove this charge group code */
    for (int z = 0; z < zones->n + 1; z++)
    {
        zones->cg_range[z] = zones->at_range[z];
    }
}
