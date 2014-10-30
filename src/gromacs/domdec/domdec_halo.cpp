/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
#include "gromacs/legacyheaders/gmx_ga2la.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_search.h"
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


static void dd_sendrecv_zone_int(const gmx_domdec_t gmx_unused *dd,
                                 int gmx_unused zone, int gmx_unused direction,
                                 int gmx_unused *buf_s, int gmx_unused n_s,
                                 int gmx_unused *buf_r, int gmx_unused n_r)
{
#ifdef GMX_MPI
    int        rank_s, rank_r;
    MPI_Status stat;

    if (direction == dddirBackward)
    {
        rank_s = dd->comm->zone_bw[zone].rank;
        rank_r = dd->comm->zone_fw[zone].rank;
    }
    else
    {
        rank_s = dd->comm->zone_fw[zone].rank;
        rank_r = dd->comm->zone_bw[zone].rank;
    }

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(int), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(int), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(int), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r, n_r*sizeof(int), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

#endif
}
static void dd_sendrecv_zone_real(const gmx_domdec_t gmx_unused *dd,
                                  int gmx_unused zone, int gmx_unused direction,
                                  real gmx_unused *buf_s, int gmx_unused n_s,
                                  real gmx_unused *buf_r, int gmx_unused n_r)
{
#ifdef GMX_MPI
    int        rank_s, rank_r;
    MPI_Status stat;

    if (direction == dddirBackward)
    {
        rank_s = dd->comm->zone_bw[zone].rank;
        rank_r = dd->comm->zone_fw[zone].rank;
    }
    else
    {
        rank_s = dd->comm->zone_fw[zone].rank;
        rank_r = dd->comm->zone_bw[zone].rank;
    }

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(real), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(real), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(real), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r, n_r*sizeof(real), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

#endif
}

static void
dd_sendrecv_zone_rvec(const gmx_domdec_t gmx_unused *dd,
                      const zone_comm_t gmx_unused *send, rvec gmx_unused *buf_s,
                      const zone_comm_t gmx_unused *recv, rvec gmx_unused *buf_r)
{
#ifdef GMX_MPI
    MPI_Status stat;

    if (send->natoms > 0 && recv->natoms > 0)
    {
        MPI_Sendrecv(buf_s[0], send->natoms*sizeof(rvec), MPI_BYTE, send->rank, 0,
                     buf_r[0], recv->natoms*sizeof(rvec), MPI_BYTE, recv->rank, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (send->natoms > 0)
    {
        MPI_Send(    buf_s[0], send->natoms*sizeof(rvec), MPI_BYTE, send->rank, 0,
                     dd->mpi_comm_all);
    }
    else if (recv->natoms > 0)
    {
        MPI_Recv(    buf_r[0], recv->natoms*sizeof(rvec), MPI_BYTE, recv->rank, 0,
                     dd->mpi_comm_all, &stat);
    }
#endif

}

static void
dd_isend_zone_rvec(const gmx_domdec_t gmx_unused *dd,
                   const zone_comm_t gmx_unused *send, rvec gmx_unused *buf_s,
                   int gmx_unused tag,
                   zones_mpi_dir_t gmx_unused *zones_mpi_dir)
{
#ifdef GMX_MPI
    if (send->natoms > 0)
    {
        MPI_Isend(buf_s[0], send->natoms*sizeof(rvec), MPI_BYTE,
                  send->rank, tag,
                  dd->mpi_comm_all, &zones_mpi_dir->req[zones_mpi_dir->nreq]);
        zones_mpi_dir->nreq++;
    }
#endif
}

static void
dd_irecv_zone_recv(const gmx_domdec_t gmx_unused *dd,
                   const zone_comm_t gmx_unused *recv, rvec gmx_unused *buf_r,
                   int gmx_unused tag,
                   zones_mpi_dir_t gmx_unused *zones_mpi_dir)
{
#ifdef GMX_MPI
    if (recv->natoms > 0)
    {
        MPI_Irecv(buf_r[0], recv->natoms*sizeof(rvec), MPI_BYTE,
                  recv->rank, tag,
                  dd->mpi_comm_all, &zones_mpi_dir->req[zones_mpi_dir->nreq]);
        zones_mpi_dir->nreq++;
    }
#endif
}

static void
dd_recv_zone_rvec(const gmx_domdec_t gmx_unused *dd,
                  const zone_comm_t gmx_unused *recv, rvec gmx_unused *buf_r)
{
#ifdef GMX_MPI
    MPI_Status stat;

    if (recv->natoms > 0)
    {
        MPI_Recv(buf_r[0], recv->natoms*sizeof(rvec), MPI_BYTE, recv->rank, 0,
                 dd->mpi_comm_all, &stat);
    }
#endif
}

static void dd_halo_send_x_zone(gmx_domdec_t *dd, int zone,
                                const matrix box, rvec x[],
                                zones_mpi_dir_t *zones_mpi_dir)
{
    const zone_comm_t *send;
    rvec               shift_vec;
    rvec              *sbuf;
    int                d, e, c, i, j;

    send = &dd->comm->zone_bw[zone];

    clear_rvec(shift_vec);

    if (send->bPBC)
    {
        const ivec *shift;

        shift = &dd->comm->zones.shift[zone];

        for (d = 0; d < DIM; d++)
        {
            if ((*shift)[d] != 0 && dd->ci[d] == 0)
            {
                for (e = 0; e < DIM; e++)
                {
                    shift_vec[e] += (*shift)[d]*box[d][e];
                }
            }
        }
    }

    sbuf = send->at_buf;

    j = 0;
    for (c = 0; c < send->ncolumn; c++)
    {
        int start, end;

        start = send->column_atom_range[2*c];
        end   = send->column_atom_range[2*c+1];

        if (!send->bPBC)
        {
            for (i = start; i < end; i++)
            {
                copy_rvec(x[i], sbuf[j++]);
            }
        }
        else if (!send->bScrew)
        {
            for (i = start; i < end; i++)
            {
                rvec_add(x[i], shift_vec, sbuf[j++]);
            }
        }
        else
        {
            /* Screw PBC */
            for (i = start; i < end; i++)
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

    assert(j == send->natoms);

    dd_isend_zone_rvec(dd,
                       send, sbuf,
                       halo_tag_X,
                       zones_mpi_dir);
}

void dd_halo_initiate_recv_x(gmx_domdec_t *dd, rvec x[])
{
    zones_mpi_t *zmpi;

    zmpi = &dd->comm->zones_mpi_x;

    assert(zmpi->recv.nreq == 0);

    int zone;

    /* Post all the non-blocking receives */
    for (zone = 1; zone < dd->comm->zones.n; zone++)
    {
        dd_irecv_zone_recv(dd,
                           &dd->comm->zone_fw[zone],
                           x + dd->comm->zones.at_range[zone],
                           halo_tag_X,
                           &zmpi->recv);
    }
}

void dd_halo_initiate_send_x(gmx_domdec_t *dd, matrix box, rvec x[])
{
    zones_mpi_dir_t *zmpid;
    int              zone;

    zmpid = &dd->comm->zones_mpi_x.send;

    assert(zmpid->nreq == 0);

    for (zone = 1; zone < dd->comm->zones.n; zone++)
    {
        /* Copy x to the send buffer and post non-blocking send */
        dd_halo_send_x_zone(dd, zone, box, x, zmpid);
    }
}

void dd_halo_complete_recv_x(gmx_domdec_t *dd)
{
    zones_mpi_t *zmpi;

    zmpi = &dd->comm->zones_mpi_x;

#ifdef GMX_MPI
    MPI_Status      stat[DD_MAXZONE];

    /* Wait for all non-blocking receives to complete */
    MPI_Waitall(zmpi->recv.nreq, zmpi->recv.req, stat);
#endif
    zmpi->recv.nreq = 0;
}

void dd_halo_complete_send_x(gmx_domdec_t *dd)
{
    zones_mpi_t *zmpi;

    zmpi = &dd->comm->zones_mpi_x;

#ifdef GMX_MPI
    MPI_Status      stat[DD_MAXZONE];

    /* Wait for all non-blocking sends to complete */
    MPI_Waitall(zmpi->send.nreq, zmpi->send.req, stat);
#endif
    zmpi->send.nreq = 0;
}

void dd_halo_move_x(gmx_domdec_t *dd, matrix box, rvec x[])
{
    dd_halo_initiate_recv_x(dd, x);
    dd_halo_initiate_send_x(dd, box, x);
    dd_halo_complete_recv_x(dd);
    dd_halo_complete_send_x(dd);
}

static void dd_halo_recv_f_zone(gmx_domdec_t *dd, int zone,
                                const rvec f_buf[],
                                rvec f[], rvec fshift[])
{
    zone_comm_t *recv;
    gmx_bool     bForcesNeedPbc;
    int          c, i, j;

    recv = &dd->comm->zone_bw[zone];

    bForcesNeedPbc = ((fshift != NULL && recv->bPBC) || recv->bScrew);

    j = 0;
    for (c = 0; c < recv->ncolumn; c++)
    {
        int start, end;

        start = recv->column_atom_range[2*c];
        end   = recv->column_atom_range[2*c+1];

        if (!bForcesNeedPbc)
        {
            for (i = start; i < end; i++)
            {
                rvec_inc(f[i], f_buf[j++]);
            }
        }
        else if (!recv->bScrew)
        {
            for (i = start; i < end; i++)
            {
                rvec_inc(f[i], f_buf[j]);
                rvec_inc(fshift[recv->shift_ind], f_buf[j]);
                j++;
            }
        }
        else
        {
            /* Screw PBC */
            for (i = start; i < end; i++)
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

    assert(j == recv->natoms);
}

void dd_halo_initiate_recv_f(gmx_domdec_t *dd)
{
    zones_mpi_t *zmpi;

    zmpi = &dd->comm->zones_mpi_f;

    assert(zmpi->recv.nreq == 0);

    /* Post all the non-blocking receives */
    int zone;

    for (zone = 1; zone < dd->comm->zones.n; zone++)
    {
        /* We can reuse the send x buffer as the receive buffer for f,
         * since the received x need to be processed before f can be
         * calculated and communicated.
         */
        dd_irecv_zone_recv(dd,
                           &dd->comm->zone_bw[zone],
                           dd->comm->zone_bw[zone].at_buf,
                           halo_tag_F,
                           &zmpi->recv);
    }
}

void dd_halo_initiate_send_f(gmx_domdec_t *dd, rvec f[])
{
    zones_mpi_t *zmpi;
    int          zone;

    zmpi = &dd->comm->zones_mpi_f;

    assert(zmpi->send.nreq == 0);

    /* Non-blocking send using direct force buffer pointers */
    for (zone = 1; zone < dd->comm->zones.n; zone++)
    {
        zone_comm_t *send;

        send = &dd->comm->zone_fw[zone];

        dd_isend_zone_rvec(dd,
                           send,
                           f + dd->comm->zones.at_range[zone],
                           halo_tag_F,
                           &zmpi->send);
    }
}

void dd_halo_complete_recv_f(gmx_domdec_t *dd, rvec f[], rvec fshift[])
{
    zones_mpi_t *zmpi;
    int          zone;

    zmpi = &dd->comm->zones_mpi_f;

#ifdef GMX_MPI
    MPI_Status      stat[DD_MAXZONE];

    if (debug)
    {
        fprintf(debug, "Waiting for force recv: %d messages\n",
                zmpi->recv.nreq);
    }

    /* Wait for all non-blocking communication to complete */
    MPI_Waitall(zmpi->recv.nreq, zmpi->recv.req, stat);
#endif
    zmpi->recv.nreq = 0;

    /* Reduce the received non-local forces with our local forces */
    for (zone = 1; zone < dd->comm->zones.n; zone++)
    {
        const rvec *f_buf;

        f_buf = dd->comm->zone_bw[zone].at_buf;
        dd_halo_recv_f_zone(dd, zone, f_buf, f, fshift);
    }
}

void dd_halo_complete_send_f(gmx_domdec_t *dd)
{
    zones_mpi_t *zmpi;

    zmpi = &dd->comm->zones_mpi_f;

#ifdef GMX_MPI
    MPI_Status      stat[DD_MAXZONE];

    if (debug)
    {
        fprintf(debug, "Waiting for force send: %d messages\n",
                zmpi->send.nreq);
    }

    /* Wait for all non-blocking sends to complete */
    MPI_Waitall(zmpi->send.nreq, zmpi->send.req, stat);
#endif
    zmpi->send.nreq = 0;
}

void dd_halo_move_f(gmx_domdec_t *dd, rvec f[], rvec fshift[])
{
    dd_halo_initiate_recv_f(dd);
    dd_halo_initiate_send_f(dd, f);
    dd_halo_complete_recv_f(dd, f, fshift);
    dd_halo_complete_send_f(dd);
}

static gmx_inline void
corner_bb_distance_rect(const ivec        zone_shift_dim,
                        gmx_bool          bDistMB,
                        const real        corner_nb[DIM],
                        const real        corner_b[DIM],
                        const nbnxn_bb_t *bb,
                        real             *r2,
                        real             *rb2)
{
    int  d;
    real r;

    /* Here we use rounding and calculate the actual distance^2
     * to the corner(s).
     */
    *r2  = 0;
    *rb2 = 0;

    /* Rectangular unit-cell, easy */
    for (d = 0; d < DIM; d++)
    {
        if (zone_shift_dim[d])
        {
            /* Distance to the corner for pair-wise interactions */
            r = bb->lower[d] - corner_nb[d];
            if (r > 0)
            {
                *r2 += r*r;
            }
            if (bDistMB)
            {
                /* Distance to the corner for multi-body interactions */
                r = bb->lower[d] - corner_b[d];
                if (r > 0)
                {
                    *rb2 += r*r;
                }
            }
        }
    }
}

static gmx_inline void
corner_bb_distance_tric(const ivec        zone_shift_dim,
                        const rvec       *normal,
                        const ivec        sumSquares,
                        gmx_bool          bDistMB,
                        const real        corner_nb[DIM],
                        const real        corner_b[DIM],
                        const nbnxn_bb_t *bb,
                        real             *r2,
                        real             *rb2)
{
    real r;
    int  d, d2;

    /* Here we use partial and approximate rounding */
    *r2  = 0;
    *rb2 = 0;
    for (d = 0; d < DIM; d++)
    {
        if (zone_shift_dim[d])
        {
            r = 0;
            for (d2 = d; d2 < DIM; d2++)
            {
                r += ((normal[d][d2] >= 0 ? bb->lower[d2] : bb->upper[d2]) - corner_nb[d2])*normal[d][d2];
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
                    *r2 += r*r;
                }
                else
                {
                    /* The/A angle between the normals is < 90 degrees.
                     * We use the maximum distance, which is an underestimate.
                     */
                    *r2 = std::max(*r2, r*r);
                }
            }

            if (bDistMB)
            {
                r = 0;
                for (d2 = d; d2 < DIM; d2++)
                {
                    r += ((normal[d][d2] >= 0 ? bb->lower[d2] : bb->upper[d2]) - corner_b[d2])*normal[d][d2];
                }
                if (r > 0)
                {
                    if (sumSquares[d])
                    {
                        *rb2 += r*r;
                    }
                    else
                    {
                        *rb2 = std::max(*rb2, r*r);
                    }
                }
            }
        }
    }
}

/* Wrapper function for corner - bounding-box distance calculation.
 * Only splits triclinic vs non-triclinic distance calculations.
 */
static gmx_inline void
corner_bb_distance(const ivec        zone_shift_dim,
                   gmx_bool          bTriclinic,
                   const rvec       *normal,
                   const ivec        sumSquares,
                   gmx_bool          bDistMB,
                   const real        corner_nb[DIM],
                   const real        corner_b[DIM],
                   const nbnxn_bb_t *bb,
                   real             *r2,
                   real             *rb2)
{
    if (bTriclinic)
    {
        corner_bb_distance_tric(zone_shift_dim, normal, sumSquares,
                                bDistMB, corner_nb, corner_b, bb,
                                r2, rb2);
    }
    else
    {
        corner_bb_distance_rect(zone_shift_dim,
                                bDistMB, corner_nb, corner_b, bb,
                                r2, rb2);
    }
}

/* Set the cell count along x and y in zone_comm.
 * Ensures enough memory is present for things that depend on the cell count.
 */
static void set_zone_comm_nc(zone_comm_t *zone_comm, int ncx, int ncy)
{
    zone_comm->ncx = ncx;
    zone_comm->ncy = ncy;
    if (zone_comm->ncx*zone_comm->ncy > zone_comm->cxy_nalloc)
    {
        zone_comm->cxy_nalloc = over_alloc_dd(zone_comm->ncx*zone_comm->ncy);
        srenew(zone_comm->cxy_natoms, zone_comm->cxy_nalloc);
        /* We usually need less space, but this is negligible on the total */
        srenew(zone_comm->column_atom_range, zone_comm->cxy_nalloc*2);
    }
}

/* Determine the corner for 2-body, corner_2b, and multi-body, corner_mb,
 * communication distance calculations.
 * Also sets the corner of our zone: corner_zone.
 * *bCornersDiffers tells if corner_2b and corner_mb are different.
 */
static void get_zone_corners(const gmx_domdec_t *dd, const matrix box,
                             int zone, const ivec zone_shift_dim,
                             const zone_dim_bounds_t *zone_bounds_dim1,
                             const zone_dim_bounds_t *zone_bounds_dim2,
                             rvec corner_2b, rvec corner_mb, rvec corner_zone,
                             gmx_bool *bCornersDiffer)
{
    int nizone = 0, izone[4], z;
    int dim;

    if (dd->bGridJump)
    {
        /* Make a list of i-zones that see our zone on the receiving end */
        nizone = 0;
        for (z = 0; z < (1 << (dd->ndim - 1)); z++)
        {
            if (zone >= dd_zp3[z][1] && zone < dd_zp3[z][2])
            {
                izone[nizone++] = dd_zp3[z][0];
            }
        }
    }

    clear_rvec(corner_2b);
    clear_rvec(corner_mb);
    clear_rvec(corner_zone);

    *bCornersDiffer = FALSE;

    for (dim = 0; dim < DIM; dim++)
    {
        real corner_2b_d, corner_mb_d, corner_z_d;

        /* This is the zone corner (simple, no staggering) */
        corner_z_d = dd->comm->cell_x0[dim];

        if (zone_shift_dim[dim] == 0 ||
            dim == dd->dim[0] ||
            !dd->bGridJump)
        {
            /* No staggering, all bounds are equal to our local bounds */
            corner_2b_d = dd->comm->cell_x0[dim];
            corner_mb_d = dd->comm->cell_x0[dim];
        }
        else if (dim == dd->dim[1])
        {
            int recv_ind0;

            /* Multi-body bonded interactions can involve all zones,
             * so we need to use the maximum over all corners.
             */
            recv_ind0   = 1 - dd_zo[zone][0];
            corner_mb_d = std::max(zone_bounds_dim1[recv_ind0    ].x0,
                                   zone_bounds_dim1[recv_ind0 + 1].x0);

            /* Determine the maximum corner over the i-zones for our zone */
            corner_2b_d = 0;
            for (z = 0; z < nizone; z++)
            {
                corner_2b_d = std::max(corner_2b_d,
                                       zone_bounds_dim1[recv_ind0 + dd_zo[izone[z]][0]].x0);
            }
        }
        else
        {
            int recv_ind0, recv_ind1;
            int sd0, sd1;

            assert(dim == dd->dim[2]);

            recv_ind0   = 1 - dd_zo[zone][0];
            recv_ind1   = 1 - dd_zo[zone][1];
            corner_mb_d = 0;
            for (sd0 = 0; sd0 < 2; sd0++)
            {
                for (sd1 = 0; sd1 < 2; sd1++)
                {
                    corner_mb_d = std::max(corner_mb_d,
                                           zone_bounds_dim2[(recv_ind0 + sd0)*3 + recv_ind1 + sd1].x0);
                }
            }

            corner_2b_d = 0;
            for (z = 0; z < nizone; z++)
            {
                corner_2b_d = std::max(corner_2b_d,
                                       zone_bounds_dim2[(recv_ind0 + dd_zo[izone[z]][0])*3 + recv_ind1 + dd_zo[izone[z]][1]].x0);
            }
        }

        if (corner_mb_d != corner_2b_d)
        {
            *bCornersDiffer = TRUE;
        }

        corner_2b[dim]   = corner_2b_d;
        corner_mb[dim]   = corner_mb_d;
        corner_zone[dim] = corner_z_d;
    }

    dd_zone_bounds_to_corner(dd, box, corner_2b);
    dd_zone_bounds_to_corner(dd, box, corner_mb);
    dd_zone_bounds_to_corner(dd, box, corner_zone);

    if (debug)
    {
        fprintf(debug, "rank %3d halo zone %d corner 2b %6.3f %6.3f %6.3f mb %6.3f %6.3f %6.3f\n",
                dd->rank, zone,
                corner_2b[XX], corner_2b[YY], corner_2b[ZZ],
                corner_mb[XX], corner_mb[YY], corner_mb[ZZ] );
    }
}

/* Determine which nbnxn grid cells (and atoms) we need to send for this zone */
static void
setup_zone_comm(gmx_domdec_t            *dd,
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
                zone_comm_t             *send)
{
    gmx_domdec_comm_t *comm;
    gmx_bool           bTriclinic;
    gmx_bool           bBondComm;
    ivec               zone_shift_dim;
    gmx_bool           bCheckRCheckBonded;
    ivec               checkRCheckBonded;
    rvec               rCheckBondedMargin2;
    rvec               corner_nb, corner_b, corner_zone;
    gmx_bool           bCornersDiffer, bDistMB, bDist2B;
    int                ncx, ncy, cxy;
    int                dim_ind, dim;

    comm = dd->comm;

    bTriclinic         = FALSE;

    bBondComm          = comm->bBondComm;

    bCheckRCheckBonded = FALSE;
    clear_ivec(checkRCheckBonded);
    clear_ivec(send->rCheckBonded);

    /* Convert the shift array indexed by dim_ind to indexed by dim */
    clear_ivec(zone_shift_dim);
    for (dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
    {
        if (dd_zo[zone][dim_ind])
        {
            dim = dd->dim[dim_ind];

            zone_shift_dim[dim] = 1;

            if (tric_dir[dim])
            {
                bTriclinic = TRUE;
            }

            if (comm->bInterCGBondeds &&
                dim < dd->npbcdim && dd->nc[dim] == 2)
            {
                /* With only 2 domains, the same pair (or triplet, ...)
                 * of atoms can potentially be present on multiple ranks.
                 * We need to check if we need to check atom distances when
                 * assigning bonded interactions on the receiving rank.
                 */
                real cellThickness, rCheckBondedMargin;

                bCheckRCheckBonded = TRUE;

                cellThickness =
                    (comm->cell_x1[dim] - comm->cell_x0[dim])*skew_fac[dim];

                /* To avoid checking distances when assigning bondeds,
                 * a slab of the cell with a thicknes of at least the maximum
                 * allowed distance between atoms in bondeds should not have
                 * any atoms communicated.
                 */
                rCheckBondedMargin       = cellThickness - std::max(comm->cutoff, comm->cutoff_mbody);
                rCheckBondedMargin2[dim] = rCheckBondedMargin*rCheckBondedMargin;

                /* If we the part of our cell that we do not communicate is
                 * less than the communication cut-off, we need to check
                 * bonded distances.
                 */
                if (rCheckBondedMargin < comm->cutoff ||
                    (!bBondComm && rCheckBondedMargin < comm->cutoff_mbody))
                {
                    send->rCheckBonded[dim] = TRUE;
                }
                else
                {
                    checkRCheckBonded[dim]  = TRUE;
                    bCheckRCheckBonded      = TRUE;
                }
            }
        }
    }

    /* For triclinic zones with shifts along multiple unit-cell vectors,
     * the exact distance calculation gets very complex, since the normals
     * to the zone planes are not orthogonal. This makes rounding the edges
     * of those zones hard.
     * We apply approximate rounding in case the angles between the normals
     * are >= 90 degrees and no rouding when < 90 degrees. This leads to
     * a few more atoms communicated than strictly necessary, but results
     * in relatively simple, efficient and less bug-prone code.
     */
    ivec sumSquares;
    int  dim_count = 0;
    int  dim_prev  = 0;
    for (dim = 0; dim < DIM; dim++)
    {
        if (zone_shift_dim[dim])
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

    /* Get the corners for non-bonded and bonded distance calculations */
    get_zone_corners(dd, box,
                     zone, zone_shift_dim, zone_bounds_dim1, zone_bounds_dim2,
                     corner_nb, corner_b, corner_zone,
                     &bCornersDiffer);

    /* Do we need to determine extra distances for multi-body bondeds?
     * Note that with bBondComm we might need distances longer than
     * the non-bonded cut-off, but with a grid without staggering (bGridJump)
     * this check is indentical to the one triggered by bDist2B below.
     */
    bDistMB = (comm->bInterCGMultiBody && bCornersDiffer);

    /* Do we need to determine extra distances for only two-body bondeds? */
    bDist2B = (bBondComm && !bDistMB);

    if (debug)
    {
        fprintf(debug, "setup zone %d BondComm %d DistMB %d Dist2B %d\n",
                zone, bBondComm, bDistMB, bDist2B);
    }

    nbnxn_get_local_grid_sizes(nbs,
                               &ncx, &ncy,
                               &send->corner0, &send->corner1,
                               &send->size_x, &send->size_y);

    set_zone_comm_nc(send, ncx, ncy);
    /* send->ncx will be updated when adding cells to send */
    send->ncx     = 0;
    send->ncolumn = 0;
    send->natoms  = 0;

    for (cxy = 0; cxy < ncx*ncy; cxy++)
    {
        int         cx, cy;
        nbnxn_bb_t  column_bb, *bb, bb_full, *bb_ptr = NULL;
        float      *bb_z_only;
        int         bb_start, nbb, atom_start, bb_natoms, column_natoms;
        real        r2, rb2;
        int         bb_first, bb_last;
        gmx_bool    bInRange;

        cx = cxy/ncy;
        cy = cxy - cx*ncy;

        /* Set this column to empty, for now */
        send->cxy_natoms[cxy] = 0;

        /* Get the bounding boxes of the column either in bb or bb_z_only */
        nbnxn_get_local_grid_column(nbs, cx, cy, &column_bb,
                                    &bb_start,
                                    &nbb, &bb, &bb_z_only,
                                    &atom_start, &bb_natoms, &column_natoms);

        if (nbb == 0)
        {
            /* Empty column */
            continue;
        }

        corner_bb_distance(zone_shift_dim, bTriclinic, normal, sumSquares,
                           bDistMB, corner_nb, corner_b, &column_bb,
                           &r2, &rb2);

        if (!(r2 < rc2_nonbonded ||
              ((bDistMB && rb2 < rc2_bonded) ||
               (bDist2B && r2  < rc2_bonded))))
        {
            /* This whole column is out of range */
            continue;
        }

        if (bb_z_only != NULL)
        {
            /* Take the x & y components from the column */
            bb_full = column_bb;
            bb_ptr  = &bb_full;
        }

        bb_first = 0;
        bb_last  = nbb - 1;
        do
        {
            if (bb != NULL)
            {
                bb_ptr = &bb[bb_last];
            }
            else
            {
                assert(bb_z_only != NULL);
                bb_full.lower[ZZ] = bb_z_only[bb_last*2+0];
                bb_full.upper[ZZ] = bb_z_only[bb_last*2+1];
            }

            corner_bb_distance(zone_shift_dim, bTriclinic, normal, sumSquares,
                               bDistMB, corner_nb, corner_b, bb_ptr,
                               &r2, &rb2);

            /* Here we check:
             * the 2-atom distance against the non-bonded cut-off,
             * the multi-body distance against the bonded cut-off
             * The 2-atom distance against the bonded cut-off.
             * The bonded check only triggers communication without bBondComm
             * or when the cell has missing bonded interactions.
             */
            bInRange = (r2 < rc2_nonbonded ||
                        (((bDistMB && rb2 < rc2_bonded) ||
                          (bDist2B && r2  < rc2_bonded)) &&
                         (!bBondComm || cell_missing_link[bb_start + bb_last])));

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
                if (bb != NULL)
                {
                    bb_ptr = &bb[bb_first];
                }
                else
                {
                    bb_full.lower[ZZ] = bb_z_only[bb_first*2+0];
                    bb_full.upper[ZZ] = bb_z_only[bb_first*2+1];
                }

                corner_bb_distance(zone_shift_dim, bTriclinic, normal, sumSquares,
                                   bDistMB, corner_nb, corner_b, bb_ptr,
                                   &r2, &rb2);

                bInRange = (r2 < rc2_nonbonded ||
                            (((bDistMB && rb2 < rc2_bonded) ||
                              (bDist2B && r2  < rc2_bonded)) &&
                             (!bBondComm || cell_missing_link[bb_start + bb_first])));

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
            int at_start, at_end;

            /* Determine the range of real atoms to send.
             * The filler atoms are always at the end, so we use min() there.
             */
            at_start = atom_start + bb_first*bb_natoms;
            at_end   = atom_start + std::min((bb_last + 1)*bb_natoms, column_natoms);

            /* Set the grid column atom count in the send buffer */
            send->cxy_natoms[cxy] = at_end - at_start;

            /* Store the atom range, needed for preparing data for sending */
            send->column_atom_range[send->ncolumn*2  ] = at_start;
            send->column_atom_range[send->ncolumn*2+1] = at_end;
            send->ncolumn++;

            send->natoms += at_end - at_start;

            /* We need to send the grid up to and including x-column cx */
            send->ncx = cx + 1;

            if (bCheckRCheckBonded)
            {
                /* It is annoying that we need all this code and effort
                 * for checking for bonded interactions, but since there
                 * is no bound on size along z of communicated cells,
                 * we can not avoid this.
                 */
                /* Construct a bounding for all cells sent for this column */
                nbnxn_bb_t send_bb;

                /* We invert upper and lower to be able to call
                 * corner_bb_distance to get the longest iso shortest distance.
                 */
                send_bb.upper[XX]     = send->corner0[XX] + cx*send->size_x;
                send_bb.lower[XX]     = send_bb.upper[XX] + send->size_x;
                send_bb.upper[YY]     = send->corner0[YY] + cy*send->size_y;
                send_bb.lower[YY]     = send_bb.upper[YY] + send->size_y;
                if (bb != NULL)
                {
                    send_bb.upper[ZZ] = bb[bb_first].lower[ZZ];
                    send_bb.lower[ZZ] = bb[bb_last ].upper[ZZ];
                }
                else
                {
                    send_bb.upper[ZZ] = bb_z_only[bb_first*2   ];
                    send_bb.lower[ZZ] = bb_z_only[bb_last*2 + 1];
                }

                for (dim = 0; dim < DIM; dim++)
                {
                    if (checkRCheckBonded[dim])
                    {
                        ivec zone_shift_one = { 0, 0, 0 };

                        zone_shift_one[dim] = 1;
                        corner_bb_distance(zone_shift_one,
                                           bTriclinic, normal, sumSquares,
                                           FALSE, corner_zone, corner_zone,
                                           &send_bb,
                                           &r2, &rb2);

                        if (r2 > rCheckBondedMargin2[dim])
                        {
                            send->rCheckBonded[dim] = 1;
                            checkRCheckBonded[dim]  = 0;
                        }
                    }
                }
            }
        }
    }

    /* Ensure the atom buffers have enough space */
    if (send->natoms > send->at_nalloc)
    {
        send->at_nalloc = over_alloc_dd(send->natoms);
        srenew(send->index_gl, send->at_nalloc);
        srenew(send->atominfo, send->at_nalloc);
        /* Reallocate the x/f buffer for MPI_Isend and MPI_Irecv */
        srenew(send->at_buf, send->at_nalloc);
    }

    /* Copy the global atom indices to the send buffer */
    int column, ind;

    ind = 0;
    for (column = 0; column < send->ncolumn; column++)
    {
        int at_start, at_end, at;

        at_start = send->column_atom_range[column*2  ];
        at_end   = send->column_atom_range[column*2+1];
        /* Copy the global atom indices to the send buffer */
        for (at = at_start; at < at_end; at++)
        {
            send->index_gl[ind] = dd->index_gl[at];
            send->atominfo[ind] = atinfo[at];
            ind++;
        }
    }
}

/* Communicate the halo communication setup */
static void
comm_zone_comm_setup(gmx_domdec_t *dd,
                     int           zone,
                     const matrix  box,
                     zone_comm_t  *send,
                     zone_comm_t  *recv)
{
    int dim_ind, dim;

    /* Communicate ncx, ncy, natoms, rCheckBonded */
    dd_sendrecv_zone_int(dd, zone, dddirBackward,
                         &send->ncx, 6,
                         &recv->ncx, 6);

    if (debug)
    {
        fprintf(debug, "For zone %d, receiving %d atoms\n", zone, recv->natoms);
    }

    /* Communicate corner0, corner1, sx, sy */
    dd_sendrecv_zone_real(dd, zone, dddirBackward,
                          &send->corner0[0], 8,
                          &recv->corner0[0], 8);

    /* Apply PBC on the receiver side */
    for (dim_ind = 0; dim_ind < DIM; dim_ind++)
    {
        dim = dd->dim[dim_ind];

        if (dim < dd->npbcdim &&
            dd_zo[zone][dim_ind] == 1 &&
            dd->ci[dim] == dd->nc[dim] - 1)
        {
            /* We are communicating over pbc along dim */
            rvec_inc(recv->corner0, box[dim]);
            rvec_inc(recv->corner1, box[dim]);
        }
    }

    set_zone_comm_nc(recv, recv->ncx, recv->ncy);

    dd_sendrecv_zone_int(dd, zone, dddirBackward,
                         send->cxy_natoms, send->ncx*send->ncy,
                         recv->cxy_natoms, recv->ncx*recv->ncy);
}

/* Accumulates the zone bounds for backward and forward neighboring cells
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
                                        const gmx_ga2la_t ga2la)
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

/* Allocate and generate *cell_missing_link_ptr, a list of of booleans
 * for a grid cells, which tell if a cell has non-local bonded interactions.
 */
static void mark_nbnxn_cells_bondcomm(const gmx_domdec_t  *dd,
                                      const int           *atinfo,
                                      const nbnxn_search_t nbs,
                                      gmx_bool           **cell_missing_link_ptr)
{
    int             bb_natoms, na, nc, c;
    const int      *atom_order;
    const t_blocka *link;
    gmx_bool       *cell_missing_link;

    link       = dd->comm->cglink;

    {
        nbnxn_bb_t  column_bb, *bb;
        float      *bb_z_only;
        int         bb_start, nbb, atom_start, column_natoms;

        /* We only need bb_natoms, but we use the same function call as
         * used later in the communication setup, so we don't need to write
         * an extra function and avoid more bug opportunities.
         */
        nbnxn_get_local_grid_column(nbs, 0, 0, &column_bb,
                                    &bb_start, &nbb, &bb, &bb_z_only,
                                    &atom_start, &bb_natoms, &column_natoms);
    }
    nbnxn_get_atomorder(nbs, &atom_order, &na);

    nc = na/bb_natoms;
    assert(nc*bb_natoms == na);

    snew(cell_missing_link, nc);

    /* Loop over all cells of the local grid */
    for (c = 0; c < nc; c++)
    {
        gmx_bool bMissingLink;
        int      i;

        /* Loop over all atoms in this cell */
        bMissingLink = FALSE;
        for (i = c*bb_natoms; i < (c + 1)*bb_natoms && !bMissingLink; i++)
        {
            int a;

            a = atom_order[i];

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

    *cell_missing_link_ptr = cell_missing_link;
}

void setup_halo_communication(gmx_domdec_t *dd,
                              const matrix box, const gmx_ddbox_t *ddbox,
                              t_forcerec *fr,
                              gmx_bool bCellsChanged)
{
    int                 zone;
    int                 c, i;
    gmx_domdec_comm_t  *comm;
    zone_dim_bounds_t   zone_bounds_dim1[3], zone_bounds_dim2[9];
    gmx_domdec_zones_t *zones;
    real                rc2_nonbonded, rc2_bonded;
    rvec                normal[DIM];
    gmx_bool           *cell_missing_link = NULL;

    if (debug)
    {
        fprintf(debug, "Setting up DD communication\n");
    }

    assert(fr->cutoff_scheme == ecutsVERLET);

    comm  = dd->comm;

    if (dd->bGridJump && dd->ndim > 1 && bCellsChanged)
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

    /* The naming is not fully consistent here (yet).
     * But the 2-body bonded interaction cut-off is max(cutoff, cutoff_mbody),
     * so there is not short and exact naming.
     */
    rc2_nonbonded = sqr(comm->cutoff);
    rc2_bonded    = sqr(comm->cutoff_mbody);

    /* Determine the normals of the zone planes */
    for (i = 0; i < DIM; i++)
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

    zones = &comm->zones;

    if (comm->bBondComm)
    {
        /* Create a boolean array for cells telling if bondeds linked to atoms
         * are not locally present, so we need to communicate those cells.
         */
        mark_nbnxn_cells_bondcomm(dd, fr->cginfo, fr->nbv->nbs,
                                  &cell_missing_link);
    }

    comm->zones.at_range[0] = 0;
    comm->zones.at_range[1] = dd->nat_home;
    dd->nat_tot             = dd->nat_home;
    dd->ncg_tot             = dd->nat_home;

    clear_ivec(comm->rCheckBonded);

    int nthread;
    int zone_ind;

    /* Here we distribute the comm setup calculation over threads.
     * Note that setup_zone_comm is actually not very time consuming.
     * Alternatively we could overlap setup_zone_comm with comm_zone_comm_setup
     * for the previous zone.
     */
    nthread = gmx_omp_nthreads_get(emntDomdec);
#pragma omp parallel for num_threads(nthread) schedule(static, 1)
    for (zone_ind = 1; zone_ind < comm->zones.n; zone_ind++)
    {
        int          zone;
        zone_comm_t *send;

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

        send = &comm->zone_bw[zone];

        /* Get the cg's for this pulse in this zone */
        setup_zone_comm(dd, zone,
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

    for (zone = 1; zone < comm->zones.n; zone++)
    {
        zone_comm_t *send, *recv;

        send = &comm->zone_bw[zone];
        recv = &comm->zone_fw[zone];

        /* This communication could be overlapped with setup_zone_comm
         * for the next zone.
         */
        comm_zone_comm_setup(dd, zone, box, send, recv);

        /* Set the atom ranges for the zones */
        dd->comm->zones.at_range[zone+1] = dd->comm->zones.at_range[zone] + recv->natoms;

        dd->nat_tot += recv->natoms;
        dd->ncg_tot += recv->natoms;
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

    for (zone = 1; zone < comm->zones.n; zone++)
    {
        zone_comm_t *send, *recv;

        send = &comm->zone_bw[zone];
        recv = &comm->zone_fw[zone];

        /* NOTE: We should consider replacing the index_gl array by an array
         *       that ombines index_gl and atominfo. This not only saves
         *       MPI calls here, but since they are mostly used together,
         *       it might also help caching. We should then copy atominfo
         *       fr->cginfo or only use the new combined array.
         *       This should be done when the group scheme is removed,
         *       since currently some code use both cg and atom indices.
         */

        /* Communicate the global atom indices */
        dd_sendrecv_zone_int(dd, zone, dddirBackward,
                             send->index_gl, send->natoms,
                             dd->index_gl + dd->comm->zones.at_range[zone], recv->natoms);

        /* Communicate the atom info flags, cheaper than extracting from mtop */
        dd_sendrecv_zone_int(dd, zone, dddirBackward,
                             send->atominfo, send->natoms,
                             fr->cginfo + dd->comm->zones.at_range[zone], recv->natoms);
        /* Set up the nbnxn grid for this non-local zone
         * using the grid setup of this zone from the home rank.
         */
        nbnxn_set_zone_grid(fr->nbv->nbs, zone,
                            recv->ncx, recv->ncy, recv->corner0, recv->corner1,
                            recv->size_x, recv->size_y,
                            recv->cxy_natoms);

        for (i = 0; i < DIM; i++)
        {
            if (recv->rCheckBonded[i])
            {
                comm->rCheckBonded[i] = 1;
            }
        }
    }

    comm->nat[ddnatHOME] = dd->nat_home;
    for (i = ddnatZONE; i < ddnatNR; i++)
    {
        comm->nat[i] = dd->nat_tot;
    }

    if (debug)
    {
        fprintf(debug, "Finished setting up DD communication, zones:");
        for (c = 1; c < zones->n; c++)
        {
            fprintf(debug, " %d", comm->zone_fw[c].natoms);
        }
        fprintf(debug, "\n");
    }

    /* TODO: remove this charge group code */
    {
        int z;

        for (z = 0; z < dd->comm->zones.n + 1; z++)
        {
            dd->comm->zones.cg_range[z] = dd->comm->zones.at_range[z];
        }
    }
}
