/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "pme_spread.h"

#include "config.h"

#include <cassert>

#include <algorithm>

#include "gromacs/ewald/pme.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "pme_grid.h"
#include "pme_internal.h"
#include "pme_simd.h"
#include "pme_spline_work.h"

/* TODO consider split of pme-spline from this file */

static void calc_interpolation_idx(const gmx_pme_t*  pme,
                                   PmeAtomComm*      atc,
                                   int               start,
                                   const pmegrids_t& pmeGrids,
                                   int               end,
                                   int               thread)
{
    int         i;
    int *       idxptr, tix, tiy, tiz;
    const real* xptr;
    real *      fptr, tx, ty, tz;
    real        rxx, ryx, ryy, rzx, rzy, rzz;
    int         nx, ny, nz;
    gmx_bool    bThreads;
    int*        thread_idx = nullptr;
    int*        tpl_n      = nullptr;
    int         thread_i;

    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    const int* g2tx = pmeGrids.g2t[XX].data();
    const int* g2ty = pmeGrids.g2t[YY].data();
    const int* g2tz = pmeGrids.g2t[ZZ].data();

    bThreads = (atc->nthread > 1);
    if (bThreads)
    {
        thread_idx = atc->thread_idx.data();

        tpl_n = atc->threadMap[thread].n;
        for (i = 0; i < atc->nthread; i++)
        {
            tpl_n[i] = 0;
        }
    }

    const real shift = c_pmeMaxUnitcellShift;

    for (i = start; i < end; i++)
    {
        xptr   = atc->x[i];
        idxptr = atc->idx[i];
        fptr   = atc->fractx[i];

        /* Fractional coordinates along box vectors, add a positive shift to ensure tx/ty/tz are positive for triclinic boxes */
        tx = nx * (xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx + shift);
        ty = ny * (xptr[YY] * ryy + xptr[ZZ] * rzy + shift);
        tz = nz * (xptr[ZZ] * rzz + shift);

        tix = static_cast<int>(tx);
        tiy = static_cast<int>(ty);
        tiz = static_cast<int>(tz);

#ifdef DEBUG
        range_check(tix, 0, c_pmeNeighborUnitcellCount * nx);
        range_check(tiy, 0, c_pmeNeighborUnitcellCount * ny);
        range_check(tiz, 0, c_pmeNeighborUnitcellCount * nz);
#endif
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
            tpl_n[i] += tpl_n[i - 1];
        }
        /* The current implementation distributes particles equally
         * over the threads, so we could actually allocate for that
         * in pme_realloc_atomcomm_things.
         */
        AtomToThreadMap& threadMap = atc->threadMap[thread];
        threadMap.i.resize(tpl_n[atc->nthread - 1]);
        /* Set tpl_n to the cumulative start */
        for (i = atc->nthread - 1; i >= 1; i--)
        {
            tpl_n[i] = tpl_n[i - 1];
        }
        tpl_n[0] = 0;

        /* Fill our thread local array with indices sorted on thread */
        for (i = start; i < end; i++)
        {
            threadMap.i[tpl_n[atc->thread_idx[i]]++] = i;
        }
        /* Now tpl_n contains the cummulative count again */
    }
}

static void make_thread_local_ind(const PmeAtomComm* atc, int thread, splinedata_t* spline)
{
    int n, t, i, start, end;

    /* Combine the indices made by each thread into one index */

    n     = 0;
    start = 0;
    for (t = 0; t < atc->nthread; t++)
    {
        const AtomToThreadMap& threadMap = atc->threadMap[t];
        /* Copy our part (start - end) from the list of thread t */
        if (thread > 0)
        {
            start = threadMap.n[thread - 1];
        }
        end = threadMap.n[thread];
        for (i = start; i < end; i++)
        {
            spline->ind[n++] = threadMap.i[i];
        }
    }

    spline->n = n;
}

// At run time, the values of order used and asserted upon mean that
// indexing out of bounds does not occur. However compilers don't
// always understand that, so we suppress this warning for this code
// region.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"

/* Macro to force loop unrolling by fixing order.
 * This gives a significant performance gain.
 */
#define CALC_SPLINE(order)                                                                               \
    {                                                                                                    \
        for (int j = 0; (j < DIM); j++)                                                                  \
        {                                                                                                \
            real dr, div;                                                                                \
            real data[PME_ORDER_MAX];                                                                    \
                                                                                                         \
            dr = xptr[j];                                                                                \
                                                                                                         \
            /* dr is relative offset from lower cell limit */                                            \
            data[(order)-1] = 0;                                                                         \
            data[1]         = dr;                                                                        \
            data[0]         = 1 - dr;                                                                    \
                                                                                                         \
            for (int k = 3; (k < (order)); k++)                                                          \
            {                                                                                            \
                div         = 1.0 / (k - 1.0);                                                           \
                data[k - 1] = div * dr * data[k - 2];                                                    \
                for (int l = 1; (l < (k - 1)); l++)                                                      \
                {                                                                                        \
                    data[k - l - 1] =                                                                    \
                            div * ((dr + l) * data[k - l - 2] + (k - l - dr) * data[k - l - 1]);         \
                }                                                                                        \
                data[0] = div * (1 - dr) * data[0];                                                      \
            }                                                                                            \
            /* differentiate */                                                                          \
            dtheta[j][i * (order) + 0] = -data[0];                                                       \
            for (int k = 1; (k < (order)); k++)                                                          \
            {                                                                                            \
                dtheta[j][i * (order) + k] = data[k - 1] - data[k];                                      \
            }                                                                                            \
                                                                                                         \
            div             = 1.0 / ((order)-1);                                                         \
            data[(order)-1] = div * dr * data[(order)-2];                                                \
            for (int l = 1; (l < ((order)-1)); l++)                                                      \
            {                                                                                            \
                data[(order)-l - 1] =                                                                    \
                        div * ((dr + l) * data[(order)-l - 2] + ((order)-l - dr) * data[(order)-l - 1]); \
            }                                                                                            \
            data[0] = div * (1 - dr) * data[0];                                                          \
                                                                                                         \
            for (int k = 0; k < (order); k++)                                                            \
            {                                                                                            \
                theta[j][i * (order) + k] = data[k];                                                     \
            }                                                                                            \
        }                                                                                                \
    }

static void make_bsplines(gmx::ArrayRef<real*> theta,
                          gmx::ArrayRef<real*> dtheta,
                          int                  order,
                          rvec                 fractx[],
                          int                  nr,
                          const int            ind[],
                          const real           coefficient[],
                          const bool           computeAllSplineCoefficients)
{
    /* construct splines for local atoms */
    int   i, ii;
    real* xptr;

    for (i = 0; i < nr; i++)
    {
        /* With free energy we do not use the coefficient check.
         * In most cases this will be more efficient than calling make_bsplines
         * twice, since usually more than half the particles have non-zero coefficients.
         */
        ii = ind[i];
        if (computeAllSplineCoefficients || coefficient[ii] != 0.0)
        {
            xptr = fractx[ii];
            assert(order >= 3 && order <= PME_ORDER_MAX);
            switch (order)
            {
                case 4: CALC_SPLINE(4) break;
                case 5: CALC_SPLINE(5) break;
                default: CALC_SPLINE(order) break;
            }
        }
    }
}

#pragma GCC diagnostic pop

/* This has to be a macro to enable full compiler optimization with xlC (and probably others too) */
#define DO_BSPLINE(order)                             \
    for (ithx = 0; (ithx < (order)); ithx++)          \
    {                                                 \
        index_x = (i0 + ithx) * pny * pnz;            \
        valx    = coefficient * thx[ithx];            \
                                                      \
        for (ithy = 0; (ithy < (order)); ithy++)      \
        {                                             \
            valxy    = valx * thy[ithy];              \
            index_xy = index_x + (j0 + ithy) * pnz;   \
                                                      \
            for (ithz = 0; (ithz < (order)); ithz++)  \
            {                                         \
                index_xyz = index_xy + (k0 + ithz);   \
                grid[index_xyz] += valxy * thz[ithz]; \
            }                                         \
        }                                             \
    }


static void spread_coefficients_bsplines_thread(pmegrid_t*            pmegrid,
                                                const PmeAtomComm*    atc,
                                                splinedata_t*         spline,
                                                const pme_spline_work gmx_unused& work)
{

    /* spread coefficients from home atoms to local grid */
    int        i, nn, n, ithx, ithy, ithz, i0, j0, k0;
    const int* idxptr;
    int        order, norder, index_x, index_xy, index_xyz;
    real       valx, valxy, coefficient;
    int        pnx, pny, pnz, ndatatot;
    int        offx, offy, offz;

#if defined PME_SIMD4_SPREAD_GATHER && !defined PME_SIMD4_UNALIGNED
    alignas(GMX_SIMD_ALIGNMENT) real thz_aligned[GMX_SIMD4_WIDTH * 2];
#endif

    pnx = pmegrid->s[XX];
    pny = pmegrid->s[YY];
    pnz = pmegrid->s[ZZ];

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    ndatatot = pnx * pny * pnz;

    real* gmx_restrict grid = pmegrid->grid.data();

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
            norder = nn * order;

            i0 = idxptr[XX] - offx;
            j0 = idxptr[YY] - offy;
            k0 = idxptr[ZZ] - offz;

            const real* thx = spline->theta.coefficients[XX] + norder;
            const real* thy = spline->theta.coefficients[YY] + norder;
            const real* thz = spline->theta.coefficients[ZZ] + norder;

            switch (order)
            {
                case 4:
#ifdef PME_SIMD4_SPREAD_GATHER
#    ifdef PME_SIMD4_UNALIGNED
#        define PME_SPREAD_SIMD4_ORDER4
#    else
#        define PME_SPREAD_SIMD4_ALIGNED
#        define PME_ORDER 4
#    endif
#    include "pme_simd4.h"
#else
                    DO_BSPLINE(4)
#endif
                    break;
                case 5:
#ifdef PME_SIMD4_SPREAD_GATHER
#    define PME_SPREAD_SIMD4_ALIGNED
#    define PME_ORDER 5
#    include "pme_simd4.h"
#else
                    DO_BSPLINE(5)
#endif
                    break;
                default: DO_BSPLINE(order) break;
            }
        }
    }
}

static void copy_local_grid(PmeAndFftGrids* grids, const int thread)
{
    const pmegrids_t*  pmegrids = &grids->pmeGrids;
    real* gmx_restrict fftgrid  = grids->fftgrid;

    ivec local_fft_ndata, local_fft_offset, local_fft_size;
    int  fft_my, fft_mz;
    int  nsy, nsz;
    ivec nf;
    int  offx, offy, offz, x, y, z, i0, i0t;
    int  d;

    gmx_parallel_3dfft_real_limits(
            grids->pfft_setup.get(), local_fft_ndata, local_fft_offset, local_fft_size);
    fft_my = local_fft_size[YY];
    fft_mz = local_fft_size[ZZ];

    const pmegrid_t* pmegrid = &pmegrids->grid_th[thread];

    nsy = pmegrid->s[YY];
    nsz = pmegrid->s[ZZ];

    for (d = 0; d < DIM; d++)
    {
        nf[d] = std::min(pmegrid->n[d] - (pmegrid->order - 1), local_fft_ndata[d] - pmegrid->offset[d]);
    }

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    /* Directly copy the non-overlapping parts of the local grids.
     * This also initializes the full grid.
     */
    const real* gmx_restrict grid_th = pmegrid->grid.data();
    for (x = 0; x < nf[XX]; x++)
    {
        for (y = 0; y < nf[YY]; y++)
        {
            i0  = ((offx + x) * fft_my + (offy + y)) * fft_mz + offz;
            i0t = (x * nsy + y) * nsz;
            for (z = 0; z < nf[ZZ]; z++)
            {
                fftgrid[i0 + z] = grid_th[i0t + z];
            }
        }
    }
}

static void reduce_threadgrid_overlap(const gmx_pme_t* pme,
                                      PmeAndFftGrids*  grids,
                                      int              thread,
                                      real*            commbuf_x,
                                      real*            commbuf_y)
{
    const pmegrids_t*  pmegrids = &grids->pmeGrids;
    real* gmx_restrict fftgrid  = grids->fftgrid;

    ivec             local_fft_ndata, local_fft_offset, local_fft_size;
    int              fft_nx, fft_ny, fft_nz;
    int              fft_my, fft_mz;
    int              buf_my = -1;
    int              nsy, nsz;
    ivec             localcopy_end, commcopy_end;
    int              offx, offy, offz, x, y, z, i0, i0t;
    int              sx, sy, sz, fx, fy, fz, tx1, ty1, tz1, ox, oy, oz;
    gmx_bool         bClearBufX, bClearBufY, bClearBufXY, bClearBuf;
    gmx_bool         bCommX, bCommY;
    int              d;
    int              thread_f;
    const pmegrid_t *pmegrid, *pmegrid_g, *pmegrid_f;
    const real*      grid_th;
    real*            commbuf = nullptr;

    gmx_parallel_3dfft_real_limits(
            grids->pfft_setup.get(), local_fft_ndata, local_fft_offset, local_fft_size);
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
        localcopy_end[d] = std::min(pmegrid->offset[d] + pmegrid->n[d] - (pmegrid->order - 1),
                                    local_fft_ndata[d]);

        /* Determine up to where our thread needs to copy from the
         * thread-local charge spreading grid to the communication buffer.
         * Note: only relevant with communication, ignored otherwise.
         */
        commcopy_end[d] = localcopy_end[d];
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
            fx += pmegrids->nc[XX];
            ox -= fft_nx;
            bCommX = (pme->nnodes_major > 1);
        }
        pmegrid_g = &pmegrids->grid_th[fx * pmegrids->nc[YY] * pmegrids->nc[ZZ]];
        ox += pmegrid_g->offset[XX];
        /* Determine the end of our part of the source grid.
         * Use our thread local source grid and target grid part
         */
        tx1 = std::min(ox + pmegrid_g->n[XX], !bCommX ? localcopy_end[XX] : commcopy_end[XX]);

        for (sy = 0; sy >= -pmegrids->nthread_comm[YY]; sy--)
        {
            fy     = pmegrid->ci[YY] + sy;
            oy     = 0;
            bCommY = FALSE;
            if (fy < 0)
            {
                fy += pmegrids->nc[YY];
                oy -= fft_ny;
                bCommY = (pme->nnodes_minor > 1);
            }
            pmegrid_g = &pmegrids->grid_th[fy * pmegrids->nc[ZZ]];
            oy += pmegrid_g->offset[YY];
            /* Determine the end of our part of the source grid.
             * Use our thread local source grid and target grid part
             */
            ty1 = std::min(oy + pmegrid_g->n[YY], !bCommY ? localcopy_end[YY] : commcopy_end[YY]);

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
                oz += pmegrid_g->offset[ZZ];
                tz1 = std::min(oz + pmegrid_g->n[ZZ], localcopy_end[ZZ]);

                if (sx == 0 && sy == 0 && sz == 0)
                {
                    /* We have already added our local contribution
                     * before calling this routine, so skip it here.
                     */
                    continue;
                }

                thread_f = (fx * pmegrids->nc[YY] + fy) * pmegrids->nc[ZZ] + fz;

                pmegrid_f = &pmegrids->grid_th[thread_f];

                grid_th = pmegrid_f->grid.data();

                nsy = pmegrid_f->s[YY];
                nsz = pmegrid_f->s[ZZ];

#ifdef DEBUG_PME_REDUCE
                printf("n%d t%d add %d  %2d %2d %2d  %2d %2d %2d  %2d-%2d %2d-%2d, %2d-%2d "
                       "%2d-%2d, %2d-%2d %2d-%2d\n",
                       pme->nodeid,
                       thread,
                       thread_f,
                       pme->pmegrid_start_ix,
                       pme->pmegrid_start_iy,
                       pme->pmegrid_start_iz,
                       sx,
                       sy,
                       sz,
                       offx - ox,
                       tx1 - ox,
                       offx,
                       tx1,
                       offy - oy,
                       ty1 - oy,
                       offy,
                       ty1,
                       offz - oz,
                       tz1 - oz,
                       offz,
                       tz1);
#endif

                if (!(bCommX || bCommY))
                {
                    /* Copy from the thread local grid to the node grid */
                    for (x = offx; x < tx1; x++)
                    {
                        for (y = offy; y < ty1; y++)
                        {
                            i0  = (x * fft_my + y) * fft_mz;
                            i0t = ((x - ox) * nsy + (y - oy)) * nsz - oz;
                            for (z = offz; z < tz1; z++)
                            {
                                fftgrid[i0 + z] += grid_th[i0t + z];
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
                        buf_my = pme->overlap[1].s2g1[pme->nodeid_minor]
                                 - pme->overlap[1].s2g0[pme->nodeid_minor + 1];
                        if (bCommX)
                        {
                            /* We index commbuf modulo the local grid size */
                            commbuf += buf_my * fft_nx * fft_nz;

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
                            i0  = (x * buf_my + y) * fft_nz;
                            i0t = ((x - ox) * nsy + (y - oy)) * nsz - oz;

                            if (bClearBuf)
                            {
                                /* First access of commbuf, initialize it */
                                for (z = offz; z < tz1; z++)
                                {
                                    commbuf[i0 + z] = grid_th[i0t + z];
                                }
                            }
                            else
                            {
                                for (z = offz; z < tz1; z++)
                                {
                                    commbuf[i0 + z] += grid_th[i0t + z];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


static void sum_fftgrid_dd(const gmx_pme_t* pme, PmeAndFftGrids* grids)
{
    real* fftgrid = grids->fftgrid;

    ivec local_fft_ndata, local_fft_offset, local_fft_size;
    int  send_index0, send_nindex;
    int  recv_nindex;
#if GMX_MPI
    MPI_Status stat;
#endif
    int recv_size_y;
    int size_yx;
    int x, y, z, indg, indb;

    /* Note that this routine is only used for forward communication.
     * Since the force gathering, unlike the coefficient spreading,
     * can be trivially parallelized over the particles,
     * the backwards process is much simpler and can use the "old"
     * communication setup.
     */

    gmx_parallel_3dfft_real_limits(
            grids->pfft_setup.get(), local_fft_ndata, local_fft_offset, local_fft_size);

    if (pme->nnodes_minor > 1)
    {
        /* Major dimension */
        const pme_overlap_t* overlap = &pme->overlap[1];

        if (pme->nnodes_major > 1)
        {
            size_yx = pme->overlap[0].comm_data[0].send_nindex;
        }
        else
        {
            size_yx = 0;
        }
#if GMX_MPI
        int datasize = (local_fft_ndata[XX] + size_yx) * local_fft_ndata[ZZ];

        int send_size_y = overlap->send_size;
#endif

        for (size_t ipulse = 0; ipulse < overlap->comm_data.size(); ipulse++)
        {
            send_index0 = overlap->comm_data[ipulse].send_index0 - overlap->comm_data[0].send_index0;
            send_nindex = overlap->comm_data[ipulse].send_nindex;
            /* We don't use recv_index0, as we always receive starting at 0 */
            recv_nindex = overlap->comm_data[ipulse].recv_nindex;
            recv_size_y = overlap->comm_data[ipulse].recv_size;

            auto* sendptr =
                    const_cast<real*>(overlap->sendbuf.data()) + send_index0 * local_fft_ndata[ZZ];
            auto* recvptr = const_cast<real*>(overlap->recvbuf.data());

            if (debug != nullptr)
            {
                fprintf(debug,
                        "PME fftgrid comm y %2d x %2d x %2d\n",
                        local_fft_ndata[XX],
                        send_nindex,
                        local_fft_ndata[ZZ]);
            }

#if GMX_MPI
            int send_id = overlap->comm_data[ipulse].send_id;
            int recv_id = overlap->comm_data[ipulse].recv_id;
            MPI_Sendrecv(sendptr,
                         send_size_y * datasize,
                         GMX_MPI_REAL,
                         send_id,
                         ipulse,
                         recvptr,
                         recv_size_y * datasize,
                         GMX_MPI_REAL,
                         recv_id,
                         ipulse,
                         overlap->mpi_comm,
                         &stat);
#endif

            for (x = 0; x < local_fft_ndata[XX]; x++)
            {
                for (y = 0; y < recv_nindex; y++)
                {
                    indg = (x * local_fft_size[YY] + y) * local_fft_size[ZZ];
                    indb = (x * recv_size_y + y) * local_fft_ndata[ZZ];
                    for (z = 0; z < local_fft_ndata[ZZ]; z++)
                    {
                        fftgrid[indg + z] += recvptr[indb + z];
                    }
                }
            }

            if (pme->nnodes_major > 1)
            {
                /* Copy from the received buffer to the send buffer for dim 0 */
                sendptr = const_cast<real*>(pme->overlap[0].sendbuf.data());
                for (x = 0; x < size_yx; x++)
                {
                    for (y = 0; y < recv_nindex; y++)
                    {
                        indg = (x * local_fft_ndata[YY] + y) * local_fft_ndata[ZZ];
                        indb = ((local_fft_ndata[XX] + x) * recv_size_y + y) * local_fft_ndata[ZZ];
                        for (z = 0; z < local_fft_ndata[ZZ]; z++)
                        {
                            sendptr[indg + z] += recvptr[indb + z];
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
        const pme_overlap_t* overlap = &pme->overlap[0];

        size_t ipulse = 0;

        send_nindex = overlap->comm_data[ipulse].send_nindex;
        /* We don't use recv_index0, as we always receive starting at 0 */
        recv_nindex = overlap->comm_data[ipulse].recv_nindex;

        if (debug != nullptr)
        {
            fprintf(debug,
                    "PME fftgrid comm x %2d x %2d x %2d\n",
                    send_nindex,
                    local_fft_ndata[YY],
                    local_fft_ndata[ZZ]);
        }

#if GMX_MPI
        int   datasize = local_fft_ndata[YY] * local_fft_ndata[ZZ];
        int   send_id  = overlap->comm_data[ipulse].send_id;
        int   recv_id  = overlap->comm_data[ipulse].recv_id;
        auto* sendptr  = const_cast<real*>(overlap->sendbuf.data());
        auto* recvptr  = const_cast<real*>(overlap->recvbuf.data());
        MPI_Sendrecv(sendptr,
                     send_nindex * datasize,
                     GMX_MPI_REAL,
                     send_id,
                     ipulse,
                     recvptr,
                     recv_nindex * datasize,
                     GMX_MPI_REAL,
                     recv_id,
                     ipulse,
                     overlap->mpi_comm,
                     &stat);
#endif

        for (x = 0; x < recv_nindex; x++)
        {
            for (y = 0; y < local_fft_ndata[YY]; y++)
            {
                indg = (x * local_fft_size[YY] + y) * local_fft_size[ZZ];
                indb = (x * local_fft_ndata[YY] + y) * local_fft_ndata[ZZ];
                for (z = 0; z < local_fft_ndata[ZZ]; z++)
                {
                    fftgrid[indg + z] += overlap->recvbuf[indb + z];
                }
            }
        }
    }
}

void spread_on_grid(const gmx_pme_t* pme,
                    PmeAtomComm*     atc,
                    PmeAndFftGrids*  grids,
                    const bool       calculateSplines,
                    const bool       doSpreading,
                    const bool       computeAllSplineCoefficients)
{
#ifdef PME_TIME_THREADS
    gmx_cycles_t  c1, c2, c3, ct1a, ct1b, ct1c;
    static double cs1 = 0, cs2 = 0, cs3 = 0;
    static double cs1a[6] = { 0, 0, 0, 0, 0, 0 };
    static int    cnt     = 0;
#endif

    const int nthread = pme->nthread;
    assert(nthread > 0);
    GMX_ASSERT(grids != nullptr || !doSpreading, "If there's no grid, we cannot be spreading");

#ifdef PME_TIME_THREADS
    c1 = omp_cyc_start();
#endif
    if (calculateSplines)
    {
#pragma omp parallel for num_threads(nthread) schedule(static)
        for (int thread = 0; thread < nthread; thread++)
        {
            try
            {
                int start, end;

                start = atc->numAtoms() * thread / nthread;
                end   = atc->numAtoms() * (thread + 1) / nthread;

                /* Compute fftgrid index for all atoms,
                 * with help of some extra variables.
                 */
                calc_interpolation_idx(pme, atc, start, grids->pmeGrids, end, thread);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
    }
#ifdef PME_TIME_THREADS
    c1 = omp_cyc_end(c1);
    cs1 += (double)c1;
#endif

#ifdef PME_TIME_THREADS
    c2 = omp_cyc_start();
#endif
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        try
        {
            splinedata_t* spline;

            /* make local bsplines  */
            if (grids == nullptr || !pme->bUseThreads)
            {
                spline = &atc->spline[0];

                spline->n = atc->numAtoms();
            }
            else
            {
                spline = &atc->spline[thread];

                if (grids->pmeGrids.nthread == 1)
                {
                    /* One thread, we operate on all coefficients */
                    spline->n = atc->numAtoms();
                }
                else
                {
                    /* Get the indices our thread should operate on */
                    make_thread_local_ind(atc, thread, spline);
                }
            }

            if (calculateSplines)
            {
                make_bsplines(spline->theta.coefficients,
                              spline->dtheta.coefficients,
                              pme->pme_order,
                              as_rvec_array(atc->fractx.data()),
                              spline->n,
                              spline->ind.data(),
                              atc->coefficient.data(),
                              computeAllSplineCoefficients);
            }

            if (doSpreading)
            {
                /* put local atoms on grid. */
                pmegrid_t& grid =
                        pme->bUseThreads ? grids->pmeGrids.grid_th[thread] : grids->pmeGrids.grid;

#ifdef PME_TIME_SPREAD
                ct1a = omp_cyc_start();
#endif
                spread_coefficients_bsplines_thread(&grid, atc, spline, *pme->spline_work);

                if (pme->bUseThreads)
                {
                    copy_local_grid(grids, thread);
                }
#ifdef PME_TIME_SPREAD
                ct1a = omp_cyc_end(ct1a);
                cs1a[thread] += (double)ct1a;
#endif
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
#ifdef PME_TIME_THREADS
    c2 = omp_cyc_end(c2);
    cs2 += (double)c2;
#endif

    if (doSpreading && pme->bUseThreads)
    {
#ifdef PME_TIME_THREADS
        c3 = omp_cyc_start();
#endif
#pragma omp parallel for num_threads(grids->pmeGrids.nthread) schedule(static)
        for (int thread = 0; thread < grids->pmeGrids.nthread; thread++)
        {
            try
            {
                reduce_threadgrid_overlap(pme,
                                          grids,
                                          thread,
                                          const_cast<real*>(pme->overlap[0].sendbuf.data()),
                                          const_cast<real*>(pme->overlap[1].sendbuf.data()));
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
#ifdef PME_TIME_THREADS
        c3 = omp_cyc_end(c3);
        cs3 += (double)c3;
#endif

        if (pme->nnodes > 1)
        {
            /* Communicate the overlapping part of the fftgrid.
             * For this communication call we need to check pme->bUseThreads
             * to have all ranks communicate here, regardless of pme->nthread.
             */
            sum_fftgrid_dd(pme, grids);
        }
    }

#ifdef PME_TIME_THREADS
    cnt++;
    if (cnt % 20 == 0)
    {
        printf("idx %.2f spread %.2f red %.2f", cs1 * 1e-9, cs2 * 1e-9, cs3 * 1e-9);
#    ifdef PME_TIME_SPREAD
        for (int thread = 0; thread < nthread; thread++)
        {
            printf(" %.2f", cs1a[thread] * 1e-9);
        }
#    endif
        printf("\n");
    }
#endif
}
