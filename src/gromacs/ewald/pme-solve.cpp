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

#include "pme-solve.h"

#include <math.h>

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/math/vec.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"

#ifdef GMX_SIMD_HAVE_REAL
/* Turn on arbitrary width SIMD intrinsics for PME solve */
#    define PME_SIMD_SOLVE
#endif

struct pme_solve_work_t
{
    /* work data for solve_pme */
    int      nalloc;
    real *   mhx;
    real *   mhy;
    real *   mhz;
    real *   m2;
    real *   denom;
    real *   tmp1_alloc;
    real *   tmp1;
    real *   tmp2;
    real *   eterm;
    real *   m2inv;

    real     energy_q;
    matrix   vir_q;
    real     energy_lj;
    matrix   vir_lj;
};

static void realloc_work(struct pme_solve_work_t *work, int nkx)
{
    if (nkx > work->nalloc)
    {
        int simd_width, i;

        work->nalloc = nkx;
        srenew(work->mhx, work->nalloc);
        srenew(work->mhy, work->nalloc);
        srenew(work->mhz, work->nalloc);
        srenew(work->m2, work->nalloc);
        /* Allocate an aligned pointer for SIMD operations, including extra
         * elements at the end for padding.
         */
#ifdef PME_SIMD_SOLVE
        simd_width = GMX_SIMD_REAL_WIDTH;
#else
        /* We can use any alignment, apart from 0, so we use 4 */
        simd_width = 4;
#endif
        sfree_aligned(work->denom);
        sfree_aligned(work->tmp1);
        sfree_aligned(work->tmp2);
        sfree_aligned(work->eterm);
        snew_aligned(work->denom, work->nalloc+simd_width, simd_width*sizeof(real));
        snew_aligned(work->tmp1,  work->nalloc+simd_width, simd_width*sizeof(real));
        snew_aligned(work->tmp2,  work->nalloc+simd_width, simd_width*sizeof(real));
        snew_aligned(work->eterm, work->nalloc+simd_width, simd_width*sizeof(real));
        srenew(work->m2inv, work->nalloc);

        /* Init all allocated elements of denom to 1 to avoid 1/0 exceptions
         * of simd padded elements.
         */
        for (i = 0; i < work->nalloc+simd_width; i++)
        {
            work->denom[i] = 1;
        }
    }
}

void pme_init_all_work(struct pme_solve_work_t **work, int nthread, int nkx)
{
    int thread;
    /* Use fft5d, order after FFT is y major, z, x minor */

    snew(*work, nthread);
    /* Allocate the work arrays thread local to optimize memory access */
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        realloc_work(&((*work)[thread]), nkx);
    }
}

static void free_work(struct pme_solve_work_t *work)
{
    sfree(work->mhx);
    sfree(work->mhy);
    sfree(work->mhz);
    sfree(work->m2);
    sfree_aligned(work->denom);
    sfree_aligned(work->tmp1);
    sfree_aligned(work->tmp2);
    sfree_aligned(work->eterm);
    sfree(work->m2inv);
}

void pme_free_all_work(struct pme_solve_work_t **work, int nthread)
{
    int thread;

    for (thread = 0; thread < nthread; thread++)
    {
        free_work(&(*work)[thread]);
    }
    sfree(work);
    *work = NULL;
}

void get_pme_ener_vir_q(struct pme_solve_work_t *work, int nthread,
                        real *mesh_energy, matrix vir)
{
    /* This function sums output over threads and should therefore
     * only be called after thread synchronization.
     */
    int thread;

    *mesh_energy = work[0].energy_q;
    copy_mat(work[0].vir_q, vir);

    for (thread = 1; thread < nthread; thread++)
    {
        *mesh_energy += work[thread].energy_q;
        m_add(vir, work[thread].vir_q, vir);
    }
}

void get_pme_ener_vir_lj(struct pme_solve_work_t *work, int nthread,
                         real *mesh_energy, matrix vir)
{
    /* This function sums output over threads and should therefore
     * only be called after thread synchronization.
     */
    int thread;

    *mesh_energy = work[0].energy_lj;
    copy_mat(work[0].vir_lj, vir);

    for (thread = 1; thread < nthread; thread++)
    {
        *mesh_energy += work[thread].energy_lj;
        m_add(vir, work[thread].vir_lj, vir);
    }
}

#if defined PME_SIMD_SOLVE
/* Calculate exponentials through SIMD */
gmx_inline static void calc_exponentials_q(int gmx_unused start, int end, real f, real *d_aligned, real *r_aligned, real *e_aligned)
{
    {
        gmx_simd_real_t       f_simd;
        gmx_simd_real_t       tmp_d1, d_inv, tmp_r, tmp_e;
        int                   kx;
        f_simd = gmx_simd_set1_r(f);
        /* We only need to calculate from start. But since start is 0 or 1
         * and we want to use aligned loads/stores, we always start from 0.
         */
        for (kx = 0; kx < end; kx += GMX_SIMD_REAL_WIDTH)
        {
            tmp_d1   = gmx_simd_load_r(d_aligned+kx);
            d_inv    = gmx_simd_inv_r(tmp_d1);
            tmp_r    = gmx_simd_load_r(r_aligned+kx);
            tmp_r    = gmx_simd_exp_r(tmp_r);
            tmp_e    = gmx_simd_mul_r(f_simd, d_inv);
            tmp_e    = gmx_simd_mul_r(tmp_e, tmp_r);
            gmx_simd_store_r(e_aligned+kx, tmp_e);
        }
    }
}
#else
gmx_inline static void calc_exponentials_q(int start, int end, real f, real *d, real *r, real *e)
{
    int kx;
    for (kx = start; kx < end; kx++)
    {
        d[kx] = 1.0/d[kx];
    }
    for (kx = start; kx < end; kx++)
    {
        r[kx] = exp(r[kx]);
    }
    for (kx = start; kx < end; kx++)
    {
        e[kx] = f*r[kx]*d[kx];
    }
}
#endif

#if defined PME_SIMD_SOLVE
/* Calculate exponentials through SIMD */
gmx_inline static void calc_exponentials_lj(int gmx_unused start, int end, real *r_aligned, real *factor_aligned, real *d_aligned)
{
    gmx_simd_real_t       tmp_r, tmp_d, tmp_fac, d_inv, tmp_mk;
    const gmx_simd_real_t sqr_PI = gmx_simd_sqrt_r(gmx_simd_set1_r(M_PI));
    int                   kx;
    for (kx = 0; kx < end; kx += GMX_SIMD_REAL_WIDTH)
    {
        /* We only need to calculate from start. But since start is 0 or 1
         * and we want to use aligned loads/stores, we always start from 0.
         */
        tmp_d = gmx_simd_load_r(d_aligned+kx);
        d_inv = gmx_simd_inv_r(tmp_d);
        gmx_simd_store_r(d_aligned+kx, d_inv);
        tmp_r = gmx_simd_load_r(r_aligned+kx);
        tmp_r = gmx_simd_exp_r(tmp_r);
        gmx_simd_store_r(r_aligned+kx, tmp_r);
        tmp_mk  = gmx_simd_load_r(factor_aligned+kx);
        tmp_fac = gmx_simd_mul_r(sqr_PI, gmx_simd_mul_r(tmp_mk, gmx_simd_erfc_r(tmp_mk)));
        gmx_simd_store_r(factor_aligned+kx, tmp_fac);
    }
}
#else
gmx_inline static void calc_exponentials_lj(int start, int end, real *r, real *tmp2, real *d)
{
    int  kx;
    real mk;
    for (kx = start; kx < end; kx++)
    {
        d[kx] = 1.0/d[kx];
    }

    for (kx = start; kx < end; kx++)
    {
        r[kx] = exp(r[kx]);
    }

    for (kx = start; kx < end; kx++)
    {
        mk       = tmp2[kx];
        tmp2[kx] = sqrt(M_PI)*mk*gmx_erfc(mk);
    }
}
#endif

int solve_pme_yzx(struct gmx_pme_t *pme, t_complex *grid,
                  real ewaldcoeff, real vol,
                  gmx_bool bEnerVir,
                  int nthread, int thread)
{
    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    t_complex               *p0;
    int                      kx, ky, kz, maxkx, maxky;
    int                      nx, ny, nz, iyz0, iyz1, iyz, iy, iz, kxstart, kxend;
    real                     mx, my, mz;
    real                     factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    real                     ets2, struct2, vfactor, ets2vf;
    real                     d1, d2, energy = 0;
    real                     by, bz;
    real                     virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;
    real                     rxx, ryx, ryy, rzx, rzy, rzz;
    struct pme_solve_work_t *work;
    real                    *mhx, *mhy, *mhz, *m2, *denom, *tmp1, *eterm, *m2inv;
    real                     mhxk, mhyk, mhzk, m2k;
    real                     corner_fac;
    ivec                     complex_order;
    ivec                     local_ndata, local_offset, local_size;
    real                     elfac;

    elfac = ONE_4PI_EPS0/pme->epsilon_r;

    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_complex_limits(pme->pfft_setup[PME_GRID_QA],
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);

    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    maxkx = (nx+1)/2;
    maxky = (ny+1)/2;

    work  = &pme->solve_work[thread];
    mhx   = work->mhx;
    mhy   = work->mhy;
    mhz   = work->mhz;
    m2    = work->m2;
    denom = work->denom;
    tmp1  = work->tmp1;
    eterm = work->eterm;
    m2inv = work->m2inv;

    iyz0 = local_ndata[YY]*local_ndata[ZZ]* thread   /nthread;
    iyz1 = local_ndata[YY]*local_ndata[ZZ]*(thread+1)/nthread;

    for (iyz = iyz0; iyz < iyz1; iyz++)
    {
        iy = iyz/local_ndata[ZZ];
        iz = iyz - iy*local_ndata[ZZ];

        ky = iy + local_offset[YY];

        if (ky < maxky)
        {
            my = ky;
        }
        else
        {
            my = (ky - ny);
        }

        by = M_PI*vol*pme->bsp_mod[YY][ky];

        kz = iz + local_offset[ZZ];

        mz = kz;

        bz = pme->bsp_mod[ZZ][kz];

        /* 0.5 correction for corner points */
        corner_fac = 1;
        if (kz == 0 || kz == (nz+1)/2)
        {
            corner_fac = 0.5;
        }

        p0 = grid + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];

        /* We should skip the k-space point (0,0,0) */
        /* Note that since here x is the minor index, local_offset[XX]=0 */
        if (local_offset[XX] > 0 || ky > 0 || kz > 0)
        {
            kxstart = local_offset[XX];
        }
        else
        {
            kxstart = local_offset[XX] + 1;
            p0++;
        }
        kxend = local_offset[XX] + local_ndata[XX];

        if (bEnerVir)
        {
            /* More expensive inner loop, especially because of the storage
             * of the mh elements in array's.
             * Because x is the minor grid index, all mh elements
             * depend on kx for triclinic unit cells.
             */

            /* Two explicit loops to avoid a conditional inside the loop */
            for (kx = kxstart; kx < maxkx; kx++)
            {
                mx = kx;

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2inv[kx] = 1.0/m2[kx];
            }

            calc_exponentials_q(kxstart, kxend, elfac, denom, tmp1, eterm);

            for (kx = kxstart; kx < kxend; kx++, p0++)
            {
                d1      = p0->re;
                d2      = p0->im;

                p0->re  = d1*eterm[kx];
                p0->im  = d2*eterm[kx];

                struct2 = 2.0*(d1*d1+d2*d2);

                tmp1[kx] = eterm[kx]*struct2;
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                ets2     = corner_fac*tmp1[kx];
                vfactor  = (factor*m2[kx] + 1.0)*2.0*m2inv[kx];
                energy  += ets2;

                ets2vf   = ets2*vfactor;
                virxx   += ets2vf*mhx[kx]*mhx[kx] - ets2;
                virxy   += ets2vf*mhx[kx]*mhy[kx];
                virxz   += ets2vf*mhx[kx]*mhz[kx];
                viryy   += ets2vf*mhy[kx]*mhy[kx] - ets2;
                viryz   += ets2vf*mhy[kx]*mhz[kx];
                virzz   += ets2vf*mhz[kx]*mhz[kx] - ets2;
            }
        }
        else
        {
            /* We don't need to calculate the energy and the virial.
             * In this case the triclinic overhead is small.
             */

            /* Two explicit loops to avoid a conditional inside the loop */

            for (kx = kxstart; kx < maxkx; kx++)
            {
                mx = kx;

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
            }

            calc_exponentials_q(kxstart, kxend, elfac, denom, tmp1, eterm);

            for (kx = kxstart; kx < kxend; kx++, p0++)
            {
                d1      = p0->re;
                d2      = p0->im;

                p0->re  = d1*eterm[kx];
                p0->im  = d2*eterm[kx];
            }
        }
    }

    if (bEnerVir)
    {
        /* Update virial with local values.
         * The virial is symmetric by definition.
         * this virial seems ok for isotropic scaling, but I'm
         * experiencing problems on semiisotropic membranes.
         * IS THAT COMMENT STILL VALID??? (DvdS, 2001/02/07).
         */
        work->vir_q[XX][XX] = 0.25*virxx;
        work->vir_q[YY][YY] = 0.25*viryy;
        work->vir_q[ZZ][ZZ] = 0.25*virzz;
        work->vir_q[XX][YY] = work->vir_q[YY][XX] = 0.25*virxy;
        work->vir_q[XX][ZZ] = work->vir_q[ZZ][XX] = 0.25*virxz;
        work->vir_q[YY][ZZ] = work->vir_q[ZZ][YY] = 0.25*viryz;

        /* This energy should be corrected for a charged system */
        work->energy_q = 0.5*energy;
    }

    /* Return the loop count */
    return local_ndata[YY]*local_ndata[XX];
}

int solve_pme_lj_yzx(struct gmx_pme_t *pme, t_complex **grid, gmx_bool bLB,
                     real ewaldcoeff, real vol,
                     gmx_bool bEnerVir, int nthread, int thread)
{
    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    int                      ig, gcount;
    int                      kx, ky, kz, maxkx, maxky;
    int                      nx, ny, nz, iy, iyz0, iyz1, iyz, iz, kxstart, kxend;
    real                     mx, my, mz;
    real                     factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    real                     ets2, ets2vf;
    real                     eterm, vterm, d1, d2, energy = 0;
    real                     by, bz;
    real                     virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;
    real                     rxx, ryx, ryy, rzx, rzy, rzz;
    real                    *mhx, *mhy, *mhz, *m2, *denom, *tmp1, *tmp2;
    real                     mhxk, mhyk, mhzk, m2k;
    struct pme_solve_work_t *work;
    real                     corner_fac;
    ivec                     complex_order;
    ivec                     local_ndata, local_offset, local_size;
    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_complex_limits(pme->pfft_setup[PME_GRID_C6A],
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);
    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    maxkx = (nx+1)/2;
    maxky = (ny+1)/2;

    work  = &pme->solve_work[thread];
    mhx   = work->mhx;
    mhy   = work->mhy;
    mhz   = work->mhz;
    m2    = work->m2;
    denom = work->denom;
    tmp1  = work->tmp1;
    tmp2  = work->tmp2;

    iyz0 = local_ndata[YY]*local_ndata[ZZ]* thread   /nthread;
    iyz1 = local_ndata[YY]*local_ndata[ZZ]*(thread+1)/nthread;

    for (iyz = iyz0; iyz < iyz1; iyz++)
    {
        iy = iyz/local_ndata[ZZ];
        iz = iyz - iy*local_ndata[ZZ];

        ky = iy + local_offset[YY];

        if (ky < maxky)
        {
            my = ky;
        }
        else
        {
            my = (ky - ny);
        }

        by = 3.0*vol*pme->bsp_mod[YY][ky]
            / (M_PI*sqrt(M_PI)*ewaldcoeff*ewaldcoeff*ewaldcoeff);

        kz = iz + local_offset[ZZ];

        mz = kz;

        bz = pme->bsp_mod[ZZ][kz];

        /* 0.5 correction for corner points */
        corner_fac = 1;
        if (kz == 0 || kz == (nz+1)/2)
        {
            corner_fac = 0.5;
        }

        kxstart = local_offset[XX];
        kxend   = local_offset[XX] + local_ndata[XX];
        if (bEnerVir)
        {
            /* More expensive inner loop, especially because of the
             * storage of the mh elements in array's.  Because x is the
             * minor grid index, all mh elements depend on kx for
             * triclinic unit cells.
             */

            /* Two explicit loops to avoid a conditional inside the loop */
            for (kx = kxstart; kx < maxkx; kx++)
            {
                mx = kx;

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
                tmp2[kx]  = sqrt(factor*m2k);
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
                tmp2[kx]  = sqrt(factor*m2k);
            }

            calc_exponentials_lj(kxstart, kxend, tmp1, tmp2, denom);

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2k   = factor*m2[kx];
                eterm = -((1.0 - 2.0*m2k)*tmp1[kx]
                          + 2.0*m2k*tmp2[kx]);
                vterm    = 3.0*(-tmp1[kx] + tmp2[kx]);
                tmp1[kx] = eterm*denom[kx];
                tmp2[kx] = vterm*denom[kx];
            }

            if (!bLB)
            {
                t_complex *p0;
                real       struct2;

                p0 = grid[0] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                for (kx = kxstart; kx < kxend; kx++, p0++)
                {
                    d1      = p0->re;
                    d2      = p0->im;

                    eterm   = tmp1[kx];
                    vterm   = tmp2[kx];
                    p0->re  = d1*eterm;
                    p0->im  = d2*eterm;

                    struct2 = 2.0*(d1*d1+d2*d2);

                    tmp1[kx] = eterm*struct2;
                    tmp2[kx] = vterm*struct2;
                }
            }
            else
            {
                real *struct2 = denom;
                real  str2;

                for (kx = kxstart; kx < kxend; kx++)
                {
                    struct2[kx] = 0.0;
                }
                /* Due to symmetry we only need to calculate 4 of the 7 terms */
                for (ig = 0; ig <= 3; ++ig)
                {
                    t_complex *p0, *p1;
                    real       scale;

                    p0    = grid[ig] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                    p1    = grid[6-ig] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                    scale = 2.0*lb_scale_factor_symm[ig];
                    for (kx = kxstart; kx < kxend; ++kx, ++p0, ++p1)
                    {
                        struct2[kx] += scale*(p0->re*p1->re + p0->im*p1->im);
                    }

                }
                for (ig = 0; ig <= 6; ++ig)
                {
                    t_complex *p0;

                    p0 = grid[ig] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                    for (kx = kxstart; kx < kxend; kx++, p0++)
                    {
                        d1     = p0->re;
                        d2     = p0->im;

                        eterm  = tmp1[kx];
                        p0->re = d1*eterm;
                        p0->im = d2*eterm;
                    }
                }
                for (kx = kxstart; kx < kxend; kx++)
                {
                    eterm    = tmp1[kx];
                    vterm    = tmp2[kx];
                    str2     = struct2[kx];
                    tmp1[kx] = eterm*str2;
                    tmp2[kx] = vterm*str2;
                }
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                ets2     = corner_fac*tmp1[kx];
                vterm    = 2.0*factor*tmp2[kx];
                energy  += ets2;
                ets2vf   = corner_fac*vterm;
                virxx   += ets2vf*mhx[kx]*mhx[kx] - ets2;
                virxy   += ets2vf*mhx[kx]*mhy[kx];
                virxz   += ets2vf*mhx[kx]*mhz[kx];
                viryy   += ets2vf*mhy[kx]*mhy[kx] - ets2;
                viryz   += ets2vf*mhy[kx]*mhz[kx];
                virzz   += ets2vf*mhz[kx]*mhz[kx] - ets2;
            }
        }
        else
        {
            /* We don't need to calculate the energy and the virial.
             *  In this case the triclinic overhead is small.
             */

            /* Two explicit loops to avoid a conditional inside the loop */

            for (kx = kxstart; kx < maxkx; kx++)
            {
                mx = kx;

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
                tmp2[kx]  = sqrt(factor*m2k);
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
                tmp2[kx]  = sqrt(factor*m2k);
            }

            calc_exponentials_lj(kxstart, kxend, tmp1, tmp2, denom);

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2k    = factor*m2[kx];
                eterm  = -((1.0 - 2.0*m2k)*tmp1[kx]
                           + 2.0*m2k*tmp2[kx]);
                tmp1[kx] = eterm*denom[kx];
            }
            gcount = (bLB ? 7 : 1);
            for (ig = 0; ig < gcount; ++ig)
            {
                t_complex *p0;

                p0 = grid[ig] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                for (kx = kxstart; kx < kxend; kx++, p0++)
                {
                    d1      = p0->re;
                    d2      = p0->im;

                    eterm   = tmp1[kx];

                    p0->re  = d1*eterm;
                    p0->im  = d2*eterm;
                }
            }
        }
    }
    if (bEnerVir)
    {
        work->vir_lj[XX][XX] = 0.25*virxx;
        work->vir_lj[YY][YY] = 0.25*viryy;
        work->vir_lj[ZZ][ZZ] = 0.25*virzz;
        work->vir_lj[XX][YY] = work->vir_lj[YY][XX] = 0.25*virxy;
        work->vir_lj[XX][ZZ] = work->vir_lj[ZZ][XX] = 0.25*virxz;
        work->vir_lj[YY][ZZ] = work->vir_lj[ZZ][YY] = 0.25*viryz;

        /* This energy should be corrected for a charged system */
        work->energy_lj = 0.5*energy;
    }
    /* Return the loop count */
    return local_ndata[YY]*local_ndata[XX];
}
