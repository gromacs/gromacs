/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "pme-spread.h"

#include "config.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/ewald/pme.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"
#include "pme-simd.h"
#include "pme-spline-work.h"

/* TODO consider split of pme-spline from this file */

static void calc_interpolation_idx(struct gmx_pme_t *pme, pme_atomcomm_t *atc,
                                   int start, int grid_index, int end, int thread)
{
    int             i;
    int            *idxptr, tix, tiy, tiz;
    real           *xptr, *fptr, tx, ty, tz;
    real            rxx, ryx, ryy, rzx, rzy, rzz;
    int             nx, ny, nz;
    int            *g2tx, *g2ty, *g2tz;
    gmx_bool        bThreads;
    int            *thread_idx = NULL;
    thread_plist_t *tpl        = NULL;
    int            *tpl_n      = NULL;
    int             thread_i;

    nx  = pme->nkx;
    ny  = pme->nky;
    nz  = pme->nkz;

    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    g2tx = pme->pmegrid[grid_index].g2t[XX];
    g2ty = pme->pmegrid[grid_index].g2t[YY];
    g2tz = pme->pmegrid[grid_index].g2t[ZZ];

    bThreads = (atc->nthread > 1);
    if (bThreads)
    {
        thread_idx = atc->thread_idx;

        tpl   = &atc->thread_plist[thread];
        tpl_n = tpl->n;
        for (i = 0; i < atc->nthread; i++)
        {
            tpl_n[i] = 0;
        }
    }

    for (i = start; i < end; i++)
    {
        xptr   = atc->x[i];
        idxptr = atc->idx[i];
        fptr   = atc->fractx[i];

        /* Fractional coordinates along box vectors, add 2.0 to make 100% sure we are positive for triclinic boxes */
        tx = nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx + 2.0 );
        ty = ny * (                  xptr[YY] * ryy + xptr[ZZ] * rzy + 2.0 );
        tz = nz * (                                   xptr[ZZ] * rzz + 2.0 );

        tix = (int)(tx);
        tiy = (int)(ty);
        tiz = (int)(tz);

        /* Because decomposition only occurs in x and y,
         * we never have a fraction correction in z.
         */
        fptr[XX] = tx - tix + pme->fshx[tix];
        fptr[YY] = ty - tiy + pme->fshy[tiy];
        fptr[ZZ] = tz - tiz;

        idxptr[XX] = pme->nnx[tix];
        idxptr[YY] = pme->nny[tiy];
        idxptr[ZZ] = pme->nnz[tiz];

#ifdef DEBUG
        range_check(idxptr[XX], 0, pme->pmegrid_nx);
        range_check(idxptr[YY], 0, pme->pmegrid_ny);
        range_check(idxptr[ZZ], 0, pme->pmegrid_nz);
#endif

        if (bThreads)
        {
            thread_i      = g2tx[idxptr[XX]] + g2ty[idxptr[YY]] + g2tz[idxptr[ZZ]];
            thread_idx[i] = thread_i;
            tpl_n[thread_i]++;
        }
    }

    if (bThreads)
    {
        /* Make a list of particle indices sorted on thread */

        /* Get the cumulative count */
        for (i = 1; i < atc->nthread; i++)
        {
            tpl_n[i] += tpl_n[i-1];
        }
        /* The current implementation distributes particles equally
         * over the threads, so we could actually allocate for that
         * in pme_realloc_atomcomm_things.
         */
        if (tpl_n[atc->nthread-1] > tpl->nalloc)
        {
            tpl->nalloc = over_alloc_large(tpl_n[atc->nthread-1]);
            srenew(tpl->i, tpl->nalloc);
        }
        /* Set tpl_n to the cumulative start */
        for (i = atc->nthread-1; i >= 1; i--)
        {
            tpl_n[i] = tpl_n[i-1];
        }
        tpl_n[0] = 0;

        /* Fill our thread local array with indices sorted on thread */
        for (i = start; i < end; i++)
        {
            tpl->i[tpl_n[atc->thread_idx[i]]++] = i;
        }
        /* Now tpl_n contains the cummulative count again */
    }
}

static void make_thread_local_ind(pme_atomcomm_t *atc,
                                  int thread, splinedata_t *spline)
{
    int             n, t, i, start, end;
    thread_plist_t *tpl;

    /* Combine the indices made by each thread into one index */

    n     = 0;
    start = 0;
    for (t = 0; t < atc->nthread; t++)
    {
        tpl = &atc->thread_plist[t];
        /* Copy our part (start - end) from the list of thread t */
        if (thread > 0)
        {
            start = tpl->n[thread-1];
        }
        end = tpl->n[thread];
        for (i = start; i < end; i++)
        {
            spline->ind[n++] = tpl->i[i];
        }
    }

    spline->n = n;
}

/* Macro to force loop unrolling by fixing order.
 * This gives a significant performance gain.
 */
#define CALC_SPLINE(order)                     \
    {                                              \
        for (int j = 0; (j < DIM); j++)            \
        {                                          \
            real dr, div;                          \
            real data[PME_ORDER_MAX];              \
                                                   \
            dr  = xptr[j];                         \
                                               \
            /* dr is relative offset from lower cell limit */ \
            data[order-1] = 0;                     \
            data[1]       = dr;                          \
            data[0]       = 1 - dr;                      \
                                               \
            for (int k = 3; (k < order); k++)      \
            {                                      \
                div       = 1.0/(k - 1.0);               \
                data[k-1] = div*dr*data[k-2];      \
                for (int l = 1; (l < (k-1)); l++)  \
                {                                  \
                    data[k-l-1] = div*((dr+l)*data[k-l-2]+(k-l-dr)* \
                                       data[k-l-1]);                \
                }                                  \
                data[0] = div*(1-dr)*data[0];      \
            }                                      \
            /* differentiate */                    \
            dtheta[j][i*order+0] = -data[0];       \
            for (int k = 1; (k < order); k++)      \
            {                                      \
                dtheta[j][i*order+k] = data[k-1] - data[k]; \
            }                                      \
                                               \
            div           = 1.0/(order - 1);                 \
            data[order-1] = div*dr*data[order-2];  \
            for (int l = 1; (l < (order-1)); l++)  \
            {                                      \
                data[order-l-1] = div*((dr+l)*data[order-l-2]+    \
                                       (order-l-dr)*data[order-l-1]); \
            }                                      \
            data[0] = div*(1 - dr)*data[0];        \
                                               \
            for (int k = 0; k < order; k++)        \
            {                                      \
                theta[j][i*order+k]  = data[k];    \
            }                                      \
        }                                          \
    }

static void make_bsplines(splinevec theta, splinevec dtheta, int order,
                          rvec fractx[], int nr, int ind[], real coefficient[],
                          gmx_bool bDoSplines)
{
    /* construct splines for local atoms */
    int   i, ii;
    real *xptr;

    for (i = 0; i < nr; i++)
    {
        /* With free energy we do not use the coefficient check.
         * In most cases this will be more efficient than calling make_bsplines
         * twice, since usually more than half the particles have non-zero coefficients.
         */
        ii = ind[i];
        if (bDoSplines || coefficient[ii] != 0.0)
        {
            xptr = fractx[ii];
            assert(order >= 4 && order <= PME_ORDER_MAX);
            switch (order)
            {
                case 4:  CALC_SPLINE(4);     break;
                case 5:  CALC_SPLINE(5);     break;
                default: CALC_SPLINE(order); break;
            }
        }
    }
}

/* This has to be a macro to enable full compiler optimization with xlC (and probably others too) */
#define DO_BSPLINE(order)                            \
    for (ithx = 0; (ithx < order); ithx++)                    \
    {                                                    \
        index_x = (i0+ithx)*pny*pnz;                     \
        valx    = coefficient*thx[ithx];                          \
                                                     \
        for (ithy = 0; (ithy < order); ithy++)                \
        {                                                \
            valxy    = valx*thy[ithy];                   \
            index_xy = index_x+(j0+ithy)*pnz;            \
                                                     \
            for (ithz = 0; (ithz < order); ithz++)            \
            {                                            \
                index_xyz        = index_xy+(k0+ithz);   \
                grid[index_xyz] += valxy*thz[ithz];      \
            }                                            \
        }                                                \
    }


static void spread_coefficients_bsplines_thread(pmegrid_t                         *pmegrid,
                                                pme_atomcomm_t                    *atc,
                                                splinedata_t                      *spline,
                                                struct pme_spline_work gmx_unused *work)
{

    /* spread coefficients from home atoms to local grid */
    real          *grid;
    int            i, nn, n, ithx, ithy, ithz, i0, j0, k0;
    int       *    idxptr;
    int            order, norder, index_x, index_xy, index_xyz;
    real           valx, valxy, coefficient;
    real          *thx, *thy, *thz;
    int            pnx, pny, pnz, ndatatot;
    int            offx, offy, offz;

#if defined PME_SIMD4_SPREAD_GATHER && !defined PME_SIMD4_UNALIGNED
    real           thz_buffer[GMX_SIMD4_WIDTH*3], *thz_aligned;

    thz_aligned = gmx_simd4_align_r(thz_buffer);
#endif

    pnx = pmegrid->s[XX];
    pny = pmegrid->s[YY];
    pnz = pmegrid->s[ZZ];

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    ndatatot = pnx*pny*pnz;
    grid     = pmegrid->grid;
    for (i = 0; i < ndatatot; i++)
    {
        grid[i] = 0;
    }

    order = pmegrid->order;

    for (nn = 0; nn < spline->n; nn++)
    {
        n           = spline->ind[nn];
        coefficient = atc->coefficient[n];

        if (coefficient != 0)
        {
            idxptr = atc->idx[n];
            norder = nn*order;

            i0   = idxptr[XX] - offx;
            j0   = idxptr[YY] - offy;
            k0   = idxptr[ZZ] - offz;

            thx = spline->theta[XX] + norder;
            thy = spline->theta[YY] + norder;
            thz = spline->theta[ZZ] + norder;

            switch (order)
            {
                case 4:
#ifdef PME_SIMD4_SPREAD_GATHER
#ifdef PME_SIMD4_UNALIGNED
#define PME_SPREAD_SIMD4_ORDER4
#else
#define PME_SPREAD_SIMD4_ALIGNED
#define PME_ORDER 4
#endif
#include "pme-simd4.h"
#else
                    DO_BSPLINE(4);
#endif
                    break;
                case 5:
#ifdef PME_SIMD4_SPREAD_GATHER
#define PME_SPREAD_SIMD4_ALIGNED
#define PME_ORDER 5
#include "pme-simd4.h"
#else
                    DO_BSPLINE(5);
#endif
                    break;
                default:
                    DO_BSPLINE(order);
                    break;
            }
        }
    }
}

static void copy_local_grid(struct gmx_pme_t *pme, pmegrids_t *pmegrids,
                            int grid_index, int thread, real *fftgrid)
{
    ivec local_fft_ndata, local_fft_offset, local_fft_size;
    int  fft_my, fft_mz;
    int  nsy, nsz;
    ivec nf;
    int  offx, offy, offz, x, y, z, i0, i0t;
    int  d;
    pmegrid_t *pmegrid;
    real *grid_th;

    gmx_parallel_3dfft_real_limits(pme->pfft_setup[grid_index],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);
    fft_my = local_fft_size[YY];
    fft_mz = local_fft_size[ZZ];

    pmegrid = &pmegrids->grid_th[thread];

    nsy = pmegrid->s[YY];
    nsz = pmegrid->s[ZZ];

    for (d = 0; d < DIM; d++)
    {
        nf[d] = std::min(pmegrid->n[d] - (pmegrid->order - 1),
                         local_fft_ndata[d] - pmegrid->offset[d]);
    }

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    /* Directly copy the non-overlapping parts of the local grids.
     * This also initializes the full grid.
     */
    grid_th = pmegrid->grid;
    for (x = 0; x < nf[XX]; x++)
    {
        for (y = 0; y < nf[YY]; y++)
        {
            i0  = ((offx + x)*fft_my + (offy + y))*fft_mz + offz;
            i0t = (x*nsy + y)*nsz;
            for (z = 0; z < nf[ZZ]; z++)
            {
                fftgrid[i0+z] = grid_th[i0t+z];
            }
        }
    }
}

static void
reduce_threadgrid_overlap(struct gmx_pme_t *pme,
                          const pmegrids_t *pmegrids, int thread,
                          real *fftgrid, real *commbuf_x, real *commbuf_y,
                          int grid_index)
{
    ivec local_fft_ndata, local_fft_offset, local_fft_size;
    int  fft_nx, fft_ny, fft_nz;
    int  fft_my, fft_mz;
    int  buf_my = -1;
    int  nsy, nsz;
    ivec localcopy_end, commcopy_end;
    int  offx, offy, offz, x, y, z, i0, i0t;
    int  sx, sy, sz, fx, fy, fz, tx1, ty1, tz1, ox, oy, oz;
    gmx_bool bClearBufX, bClearBufY, bClearBufXY, bClearBuf;
    gmx_bool bCommX, bCommY;
    int  d;
    int  thread_f;
    const pmegrid_t *pmegrid, *pmegrid_g, *pmegrid_f;
    const real *grid_th;
    real *commbuf = NULL;

    gmx_parallel_3dfft_real_limits(pme->pfft_setup[grid_index],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);
    fft_nx = local_fft_ndata[XX];
    fft_ny = local_fft_ndata[YY];
    fft_nz = local_fft_ndata[ZZ];

    fft_my = local_fft_size[YY];
    fft_mz = local_fft_size[ZZ];

    /* This routine is called when all thread have finished spreading.
     * Here each thread sums grid contributions calculated by other threads
     * to the thread local grid volume.
     * To minimize the number of grid copying operations,
     * this routines sums immediately from the pmegrid to the fftgrid.
     */

    /* Determine which part of the full node grid we should operate on,
     * this is our thread local part of the full grid.
     */
    pmegrid = &pmegrids->grid_th[thread];

    for (d = 0; d < DIM; d++)
    {
        /* Determine up to where our thread needs to copy from the
         * thread-local charge spreading grid to the rank-local FFT grid.
         * This is up to our spreading grid end minus order-1 and
         * not beyond the local FFT grid.
         */
        localcopy_end[d] =
            std::min(pmegrid->offset[d] + pmegrid->n[d] - (pmegrid->order - 1),
                     local_fft_ndata[d]);

        /* Determine up to where our thread needs to copy from the
         * thread-local charge spreading grid to the communication buffer.
         * Note: only relevant with communication, ignored otherwise.
         */
        commcopy_end[d]  = localcopy_end[d];
        if (pmegrid->ci[d] == pmegrids->nc[d] - 1)
        {
            /* The last thread should copy up to the last pme grid line.
             * When the rank-local FFT grid is narrower than pme-order,
             * we need the max below to ensure copying of all data.
             */
            commcopy_end[d] = std::max(commcopy_end[d], pme->pme_order);
        }
    }

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];


    bClearBufX  = TRUE;
    bClearBufY  = TRUE;
    bClearBufXY = TRUE;

    /* Now loop over all the thread data blocks that contribute
     * to the grid region we (our thread) are operating on.
     */
    /* Note that fft_nx/y is equal to the number of grid points
     * between the first point of our node grid and the one of the next node.
     */
    for (sx = 0; sx >= -pmegrids->nthread_comm[XX]; sx--)
    {
        fx     = pmegrid->ci[XX] + sx;
        ox     = 0;
        bCommX = FALSE;
        if (fx < 0)
        {
            fx    += pmegrids->nc[XX];
            ox    -= fft_nx;
            bCommX = (pme->nnodes_major > 1);
        }
        pmegrid_g = &pmegrids->grid_th[fx*pmegrids->nc[YY]*pmegrids->nc[ZZ]];
        ox       += pmegrid_g->offset[XX];
        /* Determine the end of our part of the source grid.
         * Use our thread local source grid and target grid part
         */
        tx1 = std::min(ox + pmegrid_g->n[XX],
                       !bCommX ? localcopy_end[XX] : commcopy_end[XX]);

        for (sy = 0; sy >= -pmegrids->nthread_comm[YY]; sy--)
        {
            fy     = pmegrid->ci[YY] + sy;
            oy     = 0;
            bCommY = FALSE;
            if (fy < 0)
            {
                fy    += pmegrids->nc[YY];
                oy    -= fft_ny;
                bCommY = (pme->nnodes_minor > 1);
            }
            pmegrid_g = &pmegrids->grid_th[fy*pmegrids->nc[ZZ]];
            oy       += pmegrid_g->offset[YY];
            /* Determine the end of our part of the source grid.
             * Use our thread local source grid and target grid part
             */
            ty1 = std::min(oy + pmegrid_g->n[YY],
                           !bCommY ? localcopy_end[YY] : commcopy_end[YY]);

            for (sz = 0; sz >= -pmegrids->nthread_comm[ZZ]; sz--)
            {
                fz = pmegrid->ci[ZZ] + sz;
                oz = 0;
                if (fz < 0)
                {
                    fz += pmegrids->nc[ZZ];
                    oz -= fft_nz;
                }
                pmegrid_g = &pmegrids->grid_th[fz];
                oz       += pmegrid_g->offset[ZZ];
                tz1       = std::min(oz + pmegrid_g->n[ZZ], localcopy_end[ZZ]);

                if (sx == 0 && sy == 0 && sz == 0)
                {
                    /* We have already added our local contribution
                     * before calling this routine, so skip it here.
                     */
                    continue;
                }

                thread_f = (fx*pmegrids->nc[YY] + fy)*pmegrids->nc[ZZ] + fz;

                pmegrid_f = &pmegrids->grid_th[thread_f];

                grid_th = pmegrid_f->grid;

                nsy = pmegrid_f->s[YY];
                nsz = pmegrid_f->s[ZZ];

#ifdef DEBUG_PME_REDUCE
                printf("n%d t%d add %d  %2d %2d %2d  %2d %2d %2d  %2d-%2d %2d-%2d, %2d-%2d %2d-%2d, %2d-%2d %2d-%2d\n",
                       pme->nodeid, thread, thread_f,
                       pme->pmegrid_start_ix,
                       pme->pmegrid_start_iy,
                       pme->pmegrid_start_iz,
                       sx, sy, sz,
                       offx-ox, tx1-ox, offx, tx1,
                       offy-oy, ty1-oy, offy, ty1,
                       offz-oz, tz1-oz, offz, tz1);
#endif

                if (!(bCommX || bCommY))
                {
                    /* Copy from the thread local grid to the node grid */
                    for (x = offx; x < tx1; x++)
                    {
                        for (y = offy; y < ty1; y++)
                        {
                            i0  = (x*fft_my + y)*fft_mz;
                            i0t = ((x - ox)*nsy + (y - oy))*nsz - oz;
                            for (z = offz; z < tz1; z++)
                            {
                                fftgrid[i0+z] += grid_th[i0t+z];
                            }
                        }
                    }
                }
                else
                {
                    /* The order of this conditional decides
                     * where the corner volume gets stored with x+y decomp.
                     */
                    if (bCommY)
                    {
                        commbuf = commbuf_y;
                        /* The y-size of the communication buffer is set by
                         * the overlap of the grid part of our local slab
                         * with the part starting at the next slab.
                         */
                        buf_my  =
                            pme->overlap[1].s2g1[pme->nodeid_minor] -
                            pme->overlap[1].s2g0[pme->nodeid_minor+1];
                        if (bCommX)
                        {
                            /* We index commbuf modulo the local grid size */
                            commbuf += buf_my*fft_nx*fft_nz;

                            bClearBuf   = bClearBufXY;
                            bClearBufXY = FALSE;
                        }
                        else
                        {
                            bClearBuf  = bClearBufY;
                            bClearBufY = FALSE;
                        }
                    }
                    else
                    {
                        commbuf    = commbuf_x;
                        buf_my     = fft_ny;
                        bClearBuf  = bClearBufX;
                        bClearBufX = FALSE;
                    }

                    /* Copy to the communication buffer */
                    for (x = offx; x < tx1; x++)
                    {
                        for (y = offy; y < ty1; y++)
                        {
                            i0  = (x*buf_my + y)*fft_nz;
                            i0t = ((x - ox)*nsy + (y - oy))*nsz - oz;

                            if (bClearBuf)
                            {
                                /* First access of commbuf, initialize it */
                                for (z = offz; z < tz1; z++)
                                {
                                    commbuf[i0+z]  = grid_th[i0t+z];
                                }
                            }
                            else
                            {
                                for (z = offz; z < tz1; z++)
                                {
                                    commbuf[i0+z] += grid_th[i0t+z];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


static void sum_fftgrid_dd(struct gmx_pme_t *pme, real *fftgrid, int grid_index)
{
    ivec local_fft_ndata, local_fft_offset, local_fft_size;
    pme_overlap_t *overlap;
    int  send_index0, send_nindex;
    int  recv_nindex;
#ifdef GMX_MPI
    MPI_Status stat;
#endif
    int  recv_size_y;
    int  ipulse, size_yx;
    real *sendptr, *recvptr;
    int  x, y, z, indg, indb;

    /* Note that this routine is only used for forward communication.
     * Since the force gathering, unlike the coefficient spreading,
     * can be trivially parallelized over the particles,
     * the backwards process is much simpler and can use the "old"
     * communication setup.
     */

    gmx_parallel_3dfft_real_limits(pme->pfft_setup[grid_index],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);

    if (pme->nnodes_minor > 1)
    {
        /* Major dimension */
        overlap = &pme->overlap[1];

        if (pme->nnodes_major > 1)
        {
            size_yx = pme->overlap[0].comm_data[0].send_nindex;
        }
        else
        {
            size_yx = 0;
        }
#ifdef GMX_MPI
        int datasize = (local_fft_ndata[XX] + size_yx)*local_fft_ndata[ZZ];

        int send_size_y = overlap->send_size;
#endif

        for (ipulse = 0; ipulse < overlap->noverlap_nodes; ipulse++)
        {
            send_index0   =
                overlap->comm_data[ipulse].send_index0 -
                overlap->comm_data[0].send_index0;
            send_nindex   = overlap->comm_data[ipulse].send_nindex;
            /* We don't use recv_index0, as we always receive starting at 0 */
            recv_nindex   = overlap->comm_data[ipulse].recv_nindex;
            recv_size_y   = overlap->comm_data[ipulse].recv_size;

            sendptr = overlap->sendbuf + send_index0*local_fft_ndata[ZZ];
            recvptr = overlap->recvbuf;

            if (debug != NULL)
            {
                fprintf(debug, "PME fftgrid comm y %2d x %2d x %2d\n",
                        local_fft_ndata[XX], send_nindex, local_fft_ndata[ZZ]);
            }

#ifdef GMX_MPI
            int send_id = overlap->send_id[ipulse];
            int recv_id = overlap->recv_id[ipulse];
            MPI_Sendrecv(sendptr, send_size_y*datasize, GMX_MPI_REAL,
                         send_id, ipulse,
                         recvptr, recv_size_y*datasize, GMX_MPI_REAL,
                         recv_id, ipulse,
                         overlap->mpi_comm, &stat);
#endif

            for (x = 0; x < local_fft_ndata[XX]; x++)
            {
                for (y = 0; y < recv_nindex; y++)
                {
                    indg = (x*local_fft_size[YY] + y)*local_fft_size[ZZ];
                    indb = (x*recv_size_y        + y)*local_fft_ndata[ZZ];
                    for (z = 0; z < local_fft_ndata[ZZ]; z++)
                    {
                        fftgrid[indg+z] += recvptr[indb+z];
                    }
                }
            }

            if (pme->nnodes_major > 1)
            {
                /* Copy from the received buffer to the send buffer for dim 0 */
                sendptr = pme->overlap[0].sendbuf;
                for (x = 0; x < size_yx; x++)
                {
                    for (y = 0; y < recv_nindex; y++)
                    {
                        indg = (x*local_fft_ndata[YY] + y)*local_fft_ndata[ZZ];
                        indb = ((local_fft_ndata[XX] + x)*recv_size_y + y)*local_fft_ndata[ZZ];
                        for (z = 0; z < local_fft_ndata[ZZ]; z++)
                        {
                            sendptr[indg+z] += recvptr[indb+z];
                        }
                    }
                }
            }
        }
    }

    /* We only support a single pulse here.
     * This is not a severe limitation, as this code is only used
     * with OpenMP and with OpenMP the (PME) domains can be larger.
     */
    if (pme->nnodes_major > 1)
    {
        /* Major dimension */
        overlap = &pme->overlap[0];

        ipulse = 0;

        send_nindex   = overlap->comm_data[ipulse].send_nindex;
        /* We don't use recv_index0, as we always receive starting at 0 */
        recv_nindex   = overlap->comm_data[ipulse].recv_nindex;

        recvptr = overlap->recvbuf;

        if (debug != NULL)
        {
            fprintf(debug, "PME fftgrid comm x %2d x %2d x %2d\n",
                    send_nindex, local_fft_ndata[YY], local_fft_ndata[ZZ]);
        }

#ifdef GMX_MPI
        int datasize = local_fft_ndata[YY]*local_fft_ndata[ZZ];
        int send_id  = overlap->send_id[ipulse];
        int recv_id  = overlap->recv_id[ipulse];
        sendptr      = overlap->sendbuf;
        MPI_Sendrecv(sendptr, send_nindex*datasize, GMX_MPI_REAL,
                     send_id, ipulse,
                     recvptr, recv_nindex*datasize, GMX_MPI_REAL,
                     recv_id, ipulse,
                     overlap->mpi_comm, &stat);
#endif

        for (x = 0; x < recv_nindex; x++)
        {
            for (y = 0; y < local_fft_ndata[YY]; y++)
            {
                indg = (x*local_fft_size[YY]  + y)*local_fft_size[ZZ];
                indb = (x*local_fft_ndata[YY] + y)*local_fft_ndata[ZZ];
                for (z = 0; z < local_fft_ndata[ZZ]; z++)
                {
                    fftgrid[indg+z] += recvptr[indb+z];
                }
            }
        }
    }
}

void spread_on_grid(struct gmx_pme_t *pme,
                    pme_atomcomm_t *atc, pmegrids_t *grids,
                    gmx_bool bCalcSplines, gmx_bool bSpread,
                    real *fftgrid, gmx_bool bDoSplines, int grid_index)
{
    int nthread, thread;
#ifdef PME_TIME_THREADS
    gmx_cycles_t c1, c2, c3, ct1a, ct1b, ct1c;
    static double cs1     = 0, cs2 = 0, cs3 = 0;
    static double cs1a[6] = {0, 0, 0, 0, 0, 0};
    static int cnt        = 0;
#endif

    nthread = pme->nthread;
    assert(nthread > 0);

#ifdef PME_TIME_THREADS
    c1 = omp_cyc_start();
#endif
    if (bCalcSplines)
    {
#pragma omp parallel for num_threads(nthread) schedule(static)
        for (thread = 0; thread < nthread; thread++)
        {
            int start, end;

            start = atc->n* thread   /nthread;
            end   = atc->n*(thread+1)/nthread;

            /* Compute fftgrid index for all atoms,
             * with help of some extra variables.
             */
            calc_interpolation_idx(pme, atc, start, grid_index, end, thread);
        }
    }
#ifdef PME_TIME_THREADS
    c1   = omp_cyc_end(c1);
    cs1 += (double)c1;
#endif

#ifdef PME_TIME_THREADS
    c2 = omp_cyc_start();
#endif
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        splinedata_t *spline;
        pmegrid_t *grid = NULL;

        /* make local bsplines  */
        if (grids == NULL || !pme->bUseThreads)
        {
            spline = &atc->spline[0];

            spline->n = atc->n;

            if (bSpread)
            {
                grid = &grids->grid;
            }
        }
        else
        {
            spline = &atc->spline[thread];

            if (grids->nthread == 1)
            {
                /* One thread, we operate on all coefficients */
                spline->n = atc->n;
            }
            else
            {
                /* Get the indices our thread should operate on */
                make_thread_local_ind(atc, thread, spline);
            }

            grid = &grids->grid_th[thread];
        }

        if (bCalcSplines)
        {
            make_bsplines(spline->theta, spline->dtheta, pme->pme_order,
                          atc->fractx, spline->n, spline->ind, atc->coefficient, bDoSplines);
        }

        if (bSpread)
        {
            /* put local atoms on grid. */
#ifdef PME_TIME_SPREAD
            ct1a = omp_cyc_start();
#endif
            spread_coefficients_bsplines_thread(grid, atc, spline, pme->spline_work);

            if (pme->bUseThreads)
            {
                copy_local_grid(pme, grids, grid_index, thread, fftgrid);
            }
#ifdef PME_TIME_SPREAD
            ct1a          = omp_cyc_end(ct1a);
            cs1a[thread] += (double)ct1a;
#endif
        }
    }
#ifdef PME_TIME_THREADS
    c2   = omp_cyc_end(c2);
    cs2 += (double)c2;
#endif

    if (bSpread && pme->bUseThreads)
    {
#ifdef PME_TIME_THREADS
        c3 = omp_cyc_start();
#endif
#pragma omp parallel for num_threads(grids->nthread) schedule(static)
        for (thread = 0; thread < grids->nthread; thread++)
        {
            reduce_threadgrid_overlap(pme, grids, thread,
                                      fftgrid,
                                      pme->overlap[0].sendbuf,
                                      pme->overlap[1].sendbuf,
                                      grid_index);
        }
#ifdef PME_TIME_THREADS
        c3   = omp_cyc_end(c3);
        cs3 += (double)c3;
#endif

        if (pme->nnodes > 1)
        {
            /* Communicate the overlapping part of the fftgrid.
             * For this communication call we need to check pme->bUseThreads
             * to have all ranks communicate here, regardless of pme->nthread.
             */
            sum_fftgrid_dd(pme, fftgrid, grid_index);
        }
    }

#ifdef PME_TIME_THREADS
    cnt++;
    if (cnt % 20 == 0)
    {
        printf("idx %.2f spread %.2f red %.2f",
               cs1*1e-9, cs2*1e-9, cs3*1e-9);
#ifdef PME_TIME_SPREAD
        for (thread = 0; thread < nthread; thread++)
        {
            printf(" %.2f", cs1a[thread]*1e-9);
        }
#endif
        printf("\n");
    }
#endif
}
