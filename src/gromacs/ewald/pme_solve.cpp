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

#include "pme_solve.h"

#include <cmath>

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

#include "pme_internal.h"
#include "pme_output.h"

#if GMX_SIMD_HAVE_REAL
/* Turn on arbitrary width SIMD intrinsics for PME solve */
#    define PME_SIMD_SOLVE
#endif

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

#ifdef PME_SIMD_SOLVE
constexpr int c_simdWidth = GMX_SIMD_REAL_WIDTH;
#else
/* We can use any alignment > 0, so we use 4 */
constexpr int c_simdWidth = 4;
#endif

// Work data for solve for a thread
struct pme_solve_work_t
{
    // Constructor, \p nkx is the number of grid points along X
    pme_solve_work_t(int nkx);

    /* work data for solve_pme */
    std::vector<real>       mhx;
    std::vector<real>       mhy;
    std::vector<real>       mhz;
    std::vector<real>       m2;
    gmx::PaddedVector<real> denom;
    gmx::PaddedVector<real> tmp1;
    gmx::PaddedVector<real> tmp2;
    gmx::PaddedVector<real> eterm;
    std::vector<real>       m2inv;

    real   energy_q;
    matrix vir_q;
    real   energy_lj;
    matrix vir_lj;
};

pme_solve_work_t::pme_solve_work_t(const int nkx)
{
    mhx.resize(nkx);
    mhy.resize(nkx);
    mhz.resize(nkx);
    m2.resize(nkx);
    denom.resizeWithPadding(nkx);
    tmp1.resizeWithPadding(nkx);
    tmp2.resizeWithPadding(nkx);
    eterm.resizeWithPadding(nkx);
    m2inv.resize(nkx);

    /* Init all allocated elements of denom to 1 to avoid 1/0 exceptions
     * of simd padded elements.
     */
    ArrayRef<real> denomPadded = denom.arrayRefWithPadding().paddedArrayRef();
    for (real& d : denomPadded)
    {
        d = 1;
    }
}

PmeSolve::PmeSolve(const int numThreads, const int nkx)
{
    /* Use fft5d, order after FFT is y major, z, x minor */

    workData_.resize(numThreads);
    /* Allocate the work arrays thread local to optimize memory access */
#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int thread = 0; thread < numThreads; thread++)
    {
        try
        {
            workData_[thread] = std::make_unique<pme_solve_work_t>(nkx);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

PmeSolve::~PmeSolve() = default;

/* Returns the smallest number >= \p that is a multiple of \p factor, \p factor must be a power of 2 */
template<unsigned int factor>
static size_t roundUpToMultipleOfFactor(size_t number)
{
    static_assert(gmx::isPowerOfTwo(factor));

    /* We need to add a most factor-1 and because factor is a power of 2,
     * we get the result by masking out the bits corresponding to factor-1.
     */
    return (number + factor - 1) & ~(factor - 1);
}

void PmeSolve::getCoulombEnergyAndVirial(PmeOutput* output) const
{
    GMX_ASSERT(output != nullptr, "Need valid output buffer");
    /* This function sums output over threads and should therefore
     * only be called after thread synchronization.
     */
    output->coulombEnergy_ = workData(0).energy_q;
    copy_mat(workData(0).vir_q, output->coulombVirial_);

    for (int thread = 1; thread < numThreads(); thread++)
    {
        output->coulombEnergy_ += workData(thread).energy_q;
        m_add(output->coulombVirial_, workData(thread).vir_q, output->coulombVirial_);
    }
}

void PmeSolve::getLJEnergyAndVirial(PmeOutput* output) const
{
    GMX_ASSERT(output != nullptr, "Need valid output buffer");
    /* This function sums output over threads and should therefore
     * only be called after thread synchronization.
     */
    output->lennardJonesEnergy_ = workData(0).energy_lj;
    copy_mat(workData(0).vir_lj, output->lennardJonesVirial_);

    for (int thread = 1; thread < numThreads(); thread++)
    {
        output->lennardJonesEnergy_ += workData(thread).energy_lj;
        m_add(output->lennardJonesVirial_, workData(thread).vir_lj, output->lennardJonesVirial_);
    }
}

#if defined PME_SIMD_SOLVE
/* Calculate exponentials through SIMD */
inline static void calc_exponentials_q(int /*unused*/,
                                       int /*unused*/,
                                       real                     f,
                                       ArrayRef<const SimdReal> d_aligned,
                                       ArrayRef<const SimdReal> r_aligned,
                                       ArrayRef<SimdReal>       e_aligned)
{
    {
        SimdReal f_simd(f);
        SimdReal tmp_d1, tmp_r, tmp_e;

        /* We only need to calculate from start. But since start is 0 or 1
         * and we want to use aligned loads/stores, we always start from 0.
         */
        GMX_ASSERT(d_aligned.size() == r_aligned.size(), "d and r must have same size");
        GMX_ASSERT(d_aligned.size() == e_aligned.size(), "d and e must have same size");
        for (size_t kx = 0; kx != d_aligned.size(); ++kx)
        {
            tmp_d1        = d_aligned[kx];
            tmp_r         = r_aligned[kx];
            tmp_r         = gmx::exp(tmp_r);
            tmp_e         = f_simd / tmp_d1;
            tmp_e         = tmp_e * tmp_r;
            e_aligned[kx] = tmp_e;
        }
    }
}
#else
inline static void
calc_exponentials_q(int start, int end, real f, ArrayRef<real> d, ArrayRef<real> r, ArrayRef<real> e)
{
    GMX_ASSERT(d.size() == r.size(), "d and r must have same size");
    GMX_ASSERT(d.size() == e.size(), "d and e must have same size");
    int kx;
    for (kx = start; kx < end; kx++)
    {
        d[kx] = 1.0 / d[kx];
    }
    for (kx = start; kx < end; kx++)
    {
        r[kx] = std::exp(r[kx]);
    }
    for (kx = start; kx < end; kx++)
    {
        e[kx] = f * r[kx] * d[kx];
    }
}
#endif

#if defined PME_SIMD_SOLVE
/* Calculate exponentials through SIMD */
inline static void calc_exponentials_lj(int /*unused*/,
                                        int /*unused*/,
                                        ArrayRef<SimdReal> r_aligned,
                                        ArrayRef<SimdReal> factor_aligned,
                                        ArrayRef<SimdReal> d_aligned)
{
    SimdReal       tmp_r, tmp_d, tmp_fac, d_inv, tmp_mk;
    const SimdReal sqr_PI = sqrt(SimdReal(M_PI));

    GMX_ASSERT(d_aligned.size() == r_aligned.size(), "d and r must have same size");
    GMX_ASSERT(d_aligned.size() == factor_aligned.size(), "d and factor must have same size");
    for (size_t kx = 0; kx != d_aligned.size(); ++kx)
    {
        /* We only need to calculate from start. But since start is 0 or 1
         * and we want to use aligned loads/stores, we always start from 0.
         */
        tmp_d              = d_aligned[kx];
        d_inv              = SimdReal(1.0) / tmp_d;
        d_aligned[kx]      = d_inv;
        tmp_r              = r_aligned[kx];
        tmp_r              = gmx::exp(tmp_r);
        r_aligned[kx]      = tmp_r;
        tmp_mk             = factor_aligned[kx];
        tmp_fac            = sqr_PI * tmp_mk * erfc(tmp_mk);
        factor_aligned[kx] = tmp_fac;
    }
}
#else
inline static void
calc_exponentials_lj(int start, int end, ArrayRef<real> r, ArrayRef<real> tmp2, ArrayRef<real> d)
{
    int  kx;
    real mk;
    GMX_ASSERT(d.size() == r.size(), "d and r must have same size");
    GMX_ASSERT(d.size() == tmp2.size(), "d and tmp2 must have same size");
    for (kx = start; kx < end; kx++)
    {
        d[kx] = 1.0 / d[kx];
    }

    for (kx = start; kx < end; kx++)
    {
        r[kx] = std::exp(r[kx]);
    }

    for (kx = start; kx < end; kx++)
    {
        mk       = tmp2[kx];
        tmp2[kx] = std::sqrt(M_PI) * mk * std::erfc(mk);
    }
}
#endif

#if defined PME_SIMD_SOLVE
using PME_T = SimdReal;
#else
using PME_T = real;
#endif

int PmeSolve::solveCoulombYZX(const gmx_pme_t& pme,
                              t_complex*       grid,
                              const real       vol,
                              const bool       computeEnergyAndVirial,
                              const int        thread)
{
    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    t_complex* p0;
    int        kx, ky, kz, maxkx, maxky;
    int        iyz0, iyz1, iyz, iy, iz, kxstart, kxend;
    real       mx, my, mz;
    real       ewaldcoeff = pme.ewaldcoeff_q;
    real       factor     = M_PI * M_PI / (ewaldcoeff * ewaldcoeff);
    real       ets2, struct2, vfactor, ets2vf;
    real       d1, d2, energy = 0;
    real       by, bz;
    real       virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;
    real       mhxk, mhyk, mhzk, m2k;
    real       corner_fac;
    ivec       complex_order;
    ivec       local_ndata, local_offset, local_size;

    const real elfac = gmx::c_one4PiEps0 / pme.epsilon_r;

    const int nx = pme.nkx;
    const int ny = pme.nky;
    const int nz = pme.nkz;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_complex_limits(
            pme.gridsCoulomb[0].pfft_setup.get(), complex_order, local_ndata, local_offset, local_size);

    const real rxx = pme.recipbox[XX][XX];
    const real ryx = pme.recipbox[YY][XX];
    const real ryy = pme.recipbox[YY][YY];
    const real rzx = pme.recipbox[ZZ][XX];
    const real rzy = pme.recipbox[ZZ][YY];
    const real rzz = pme.recipbox[ZZ][ZZ];

    GMX_ASSERT(rxx != 0.0, "Someone broke the reciprocal box again");

    maxkx = (nx + 1) / 2;
    maxky = (ny + 1) / 2;

    const int nthread = numThreads();

    pme_solve_work_t& work = workData(thread);

    real* gmx_restrict mhx   = work.mhx.data();
    real* gmx_restrict mhy   = work.mhy.data();
    real* gmx_restrict mhz   = work.mhz.data();
    real* gmx_restrict m2    = work.m2.data();
    real* gmx_restrict denom = work.denom.data();
    real* gmx_restrict tmp1  = work.tmp1.data();
    real* gmx_restrict eterm = work.eterm.data();
    real* gmx_restrict m2inv = work.m2inv.data();

    iyz0 = local_ndata[YY] * local_ndata[ZZ] * thread / nthread;
    iyz1 = local_ndata[YY] * local_ndata[ZZ] * (thread + 1) / nthread;

    for (iyz = iyz0; iyz < iyz1; iyz++)
    {
        iy = iyz / local_ndata[ZZ];
        iz = iyz - iy * local_ndata[ZZ];

        ky = iy + local_offset[YY];

        if (ky < maxky)
        {
            my = ky;
        }
        else
        {
            my = (ky - ny);
        }

        by = M_PI * vol * pme.bsp_mod[YY][ky];

        kz = iz + local_offset[ZZ];

        mz = kz;

        bz = pme.bsp_mod[ZZ][kz];

        /* 0.5 correction for corner points */
        corner_fac = 1;
        if (kz == 0 || kz == (nz + 1) / 2)
        {
            corner_fac = 0.5;
        }

        p0 = grid + iy * local_size[ZZ] * local_size[XX] + iz * local_size[XX];

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

        if (computeEnergyAndVirial)
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
                m2k       = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = m2k * bz * by * pme.bsp_mod[XX][kx];
                tmp1[kx]  = -factor * m2k;
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = m2k * bz * by * pme.bsp_mod[XX][kx];
                tmp1[kx]  = -factor * m2k;
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2inv[kx] = 1.0 / m2[kx];
            }

            calc_exponentials_q(
                    kxstart,
                    kxend,
                    elfac,
                    ArrayRef<PME_T>(denom, denom + roundUpToMultipleOfFactor<c_simdWidth>(kxend)),
                    ArrayRef<PME_T>(tmp1, tmp1 + roundUpToMultipleOfFactor<c_simdWidth>(kxend)),
                    ArrayRef<PME_T>(eterm, eterm + roundUpToMultipleOfFactor<c_simdWidth>(kxend)));

            for (kx = kxstart; kx < kxend; kx++, p0++)
            {
                d1 = p0->re;
                d2 = p0->im;

                p0->re = d1 * eterm[kx];
                p0->im = d2 * eterm[kx];

                struct2 = 2.0 * (d1 * d1 + d2 * d2);

                tmp1[kx] = eterm[kx] * struct2;
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                ets2    = corner_fac * tmp1[kx];
                vfactor = (factor * m2[kx] + 1.0) * 2.0 * m2inv[kx];
                energy += ets2;

                ets2vf = ets2 * vfactor;
                virxx += ets2vf * mhx[kx] * mhx[kx] - ets2;
                virxy += ets2vf * mhx[kx] * mhy[kx];
                virxz += ets2vf * mhx[kx] * mhz[kx];
                viryy += ets2vf * mhy[kx] * mhy[kx] - ets2;
                viryz += ets2vf * mhy[kx] * mhz[kx];
                virzz += ets2vf * mhz[kx] * mhz[kx] - ets2;
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
                m2k       = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                denom[kx] = m2k * bz * by * pme.bsp_mod[XX][kx];
                tmp1[kx]  = -factor * m2k;
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                denom[kx] = m2k * bz * by * pme.bsp_mod[XX][kx];
                tmp1[kx]  = -factor * m2k;
            }

            calc_exponentials_q(
                    kxstart,
                    kxend,
                    elfac,
                    ArrayRef<PME_T>(denom, denom + roundUpToMultipleOfFactor<c_simdWidth>(kxend)),
                    ArrayRef<PME_T>(tmp1, tmp1 + roundUpToMultipleOfFactor<c_simdWidth>(kxend)),
                    ArrayRef<PME_T>(eterm, eterm + roundUpToMultipleOfFactor<c_simdWidth>(kxend)));


            for (kx = kxstart; kx < kxend; kx++, p0++)
            {
                d1 = p0->re;
                d2 = p0->im;

                p0->re = d1 * eterm[kx];
                p0->im = d2 * eterm[kx];
            }
        }
    }

    if (computeEnergyAndVirial)
    {
        /* Update virial with local values.
         * The virial is symmetric by definition.
         * this virial seems ok for isotropic scaling, but I'm
         * experiencing problems on semiisotropic membranes.
         * IS THAT COMMENT STILL VALID??? (DvdS, 2001/02/07).
         */
        work.vir_q[XX][XX] = 0.25 * virxx;
        work.vir_q[YY][YY] = 0.25 * viryy;
        work.vir_q[ZZ][ZZ] = 0.25 * virzz;
        work.vir_q[XX][YY] = work.vir_q[YY][XX] = 0.25 * virxy;
        work.vir_q[XX][ZZ] = work.vir_q[ZZ][XX] = 0.25 * virxz;
        work.vir_q[YY][ZZ] = work.vir_q[ZZ][YY] = 0.25 * viryz;

        /* This energy should be corrected for a charged system */
        work.energy_q = 0.5 * energy;
    }

    /* Return the loop count over all threads */
    return local_ndata[YY] * local_ndata[ZZ] * local_ndata[XX];
}

int PmeSolve::solveLJYZX(const gmx_pme_t&              pme,
                         gmx::ArrayRef<PmeAndFftGrids> grids,
                         const bool                    useLBCombinationRule,
                         const real                    vol,
                         const bool                    computeEnergyAndVirial,
                         const int                     thread)
{
    GMX_ASSERT(!useLBCombinationRule || gmx::ssize(grids) == sc_numGridsLJLB,
               "Expect 7 grids for LJ with LB comb.rule.");

    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    int  kx, ky, kz, maxkx, maxky;
    int  iy, iyz0, iyz1, iyz, iz, kxstart, kxend;
    real mx, my, mz;
    real ewaldcoeff = pme.ewaldcoeff_lj;
    real factor     = M_PI * M_PI / (ewaldcoeff * ewaldcoeff);
    real ets2, ets2vf;
    real eterm, vterm, d1, d2, energy = 0;
    real by, bz;
    real virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;
    real mhxk, mhyk, mhzk, m2k;
    real corner_fac;
    ivec complex_order;
    ivec local_ndata, local_offset, local_size;

    const int nx = pme.nkx;
    const int ny = pme.nky;
    const int nz = pme.nkz;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_complex_limits(
            grids[0].pfft_setup.get(), complex_order, local_ndata, local_offset, local_size);
    const real rxx = pme.recipbox[XX][XX];
    const real ryx = pme.recipbox[YY][XX];
    const real ryy = pme.recipbox[YY][YY];
    const real rzx = pme.recipbox[ZZ][XX];
    const real rzy = pme.recipbox[ZZ][YY];
    const real rzz = pme.recipbox[ZZ][ZZ];

    maxkx = (nx + 1) / 2;
    maxky = (ny + 1) / 2;

    const int nthread = numThreads();

    pme_solve_work_t& work = workData(thread);

    real* gmx_restrict mhx   = work.mhx.data();
    real* gmx_restrict mhy   = work.mhy.data();
    real* gmx_restrict mhz   = work.mhz.data();
    real* gmx_restrict m2    = work.m2.data();
    real* gmx_restrict denom = work.denom.data();
    real* gmx_restrict tmp1  = work.tmp1.data();
    real* gmx_restrict tmp2  = work.tmp2.data();

    iyz0 = local_ndata[YY] * local_ndata[ZZ] * thread / nthread;
    iyz1 = local_ndata[YY] * local_ndata[ZZ] * (thread + 1) / nthread;

    for (iyz = iyz0; iyz < iyz1; iyz++)
    {
        iy = iyz / local_ndata[ZZ];
        iz = iyz - iy * local_ndata[ZZ];

        ky = iy + local_offset[YY];

        if (ky < maxky)
        {
            my = ky;
        }
        else
        {
            my = (ky - ny);
        }

        by = 3.0 * vol * pme.bsp_mod[YY][ky]
             / (M_PI * std::sqrt(M_PI) * ewaldcoeff * ewaldcoeff * ewaldcoeff);

        kz = iz + local_offset[ZZ];

        mz = kz;

        bz = pme.bsp_mod[ZZ][kz];

        /* 0.5 correction for corner points */
        corner_fac = 1;
        if (kz == 0 || kz == (nz + 1) / 2)
        {
            corner_fac = 0.5;
        }

        kxstart = local_offset[XX];
        kxend   = local_offset[XX] + local_ndata[XX];
        if (computeEnergyAndVirial)
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
                m2k       = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz * by * pme.bsp_mod[XX][kx];
                tmp1[kx]  = -factor * m2k;
                tmp2[kx]  = std::sqrt(factor * m2k);
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz * by * pme.bsp_mod[XX][kx];
                tmp1[kx]  = -factor * m2k;
                tmp2[kx]  = std::sqrt(factor * m2k);
            }
            /* Clear padding elements to avoid (harmless) fp exceptions */
            const int kxendSimd = roundUpToMultipleOfFactor<c_simdWidth>(kxend);
            for (; kx < kxendSimd; kx++)
            {
                tmp1[kx] = 0;
                tmp2[kx] = 0;
            }

            calc_exponentials_lj(
                    kxstart,
                    kxend,
                    ArrayRef<PME_T>(tmp1, tmp1 + roundUpToMultipleOfFactor<c_simdWidth>(kxend)),
                    ArrayRef<PME_T>(tmp2, tmp2 + roundUpToMultipleOfFactor<c_simdWidth>(kxend)),
                    ArrayRef<PME_T>(denom, denom + roundUpToMultipleOfFactor<c_simdWidth>(kxend)));

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2k      = factor * m2[kx];
                eterm    = -((1.0 - 2.0 * m2k) * tmp1[kx] + 2.0 * m2k * tmp2[kx]);
                vterm    = 3.0 * (-tmp1[kx] + tmp2[kx]);
                tmp1[kx] = eterm * denom[kx];
                tmp2[kx] = vterm * denom[kx];
            }

            if (!useLBCombinationRule)
            {
                t_complex* p0;
                real       struct2;

                p0 = grids[0].cfftgrid + iy * local_size[ZZ] * local_size[XX] + iz * local_size[XX];
                for (kx = kxstart; kx < kxend; kx++, p0++)
                {
                    d1 = p0->re;
                    d2 = p0->im;

                    eterm  = tmp1[kx];
                    vterm  = tmp2[kx];
                    p0->re = d1 * eterm;
                    p0->im = d2 * eterm;

                    struct2 = 2.0 * (d1 * d1 + d2 * d2);

                    tmp1[kx] = eterm * struct2;
                    tmp2[kx] = vterm * struct2;
                }
            }
            else
            {
                real* struct2 = denom;
                real  str2;

                for (kx = kxstart; kx < kxend; kx++)
                {
                    struct2[kx] = 0.0;
                }
                /* Due to symmetry we only need to calculate 4 of the 7 terms */
                for (int ig = 0; ig <= 3; ++ig)
                {
                    t_complex *p0, *p1;
                    real       scale;

                    p0 = grids[ig].cfftgrid + iy * local_size[ZZ] * local_size[XX] + iz * local_size[XX];
                    p1 = grids[6 - ig].cfftgrid + iy * local_size[ZZ] * local_size[XX]
                         + iz * local_size[XX];
                    scale = 2.0 * lb_scale_factor_symm[ig];
                    for (kx = kxstart; kx < kxend; ++kx, ++p0, ++p1)
                    {
                        struct2[kx] += scale * (p0->re * p1->re + p0->im * p1->im);
                    }
                }
                for (int ig = 0; ig <= 6; ++ig)
                {
                    t_complex* p0;

                    p0 = grids[ig].cfftgrid + iy * local_size[ZZ] * local_size[XX] + iz * local_size[XX];
                    for (kx = kxstart; kx < kxend; kx++, p0++)
                    {
                        d1 = p0->re;
                        d2 = p0->im;

                        eterm  = tmp1[kx];
                        p0->re = d1 * eterm;
                        p0->im = d2 * eterm;
                    }
                }
                for (kx = kxstart; kx < kxend; kx++)
                {
                    eterm    = tmp1[kx];
                    vterm    = tmp2[kx];
                    str2     = struct2[kx];
                    tmp1[kx] = eterm * str2;
                    tmp2[kx] = vterm * str2;
                }
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                ets2  = corner_fac * tmp1[kx];
                vterm = 2.0 * factor * tmp2[kx];
                energy += ets2;
                ets2vf = corner_fac * vterm;
                virxx += ets2vf * mhx[kx] * mhx[kx] - ets2;
                virxy += ets2vf * mhx[kx] * mhy[kx];
                virxz += ets2vf * mhx[kx] * mhz[kx];
                viryy += ets2vf * mhy[kx] * mhy[kx] - ets2;
                viryz += ets2vf * mhy[kx] * mhz[kx];
                virzz += ets2vf * mhz[kx] * mhz[kx] - ets2;
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
                m2k       = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz * by * pme.bsp_mod[XX][kx];
                tmp1[kx]  = -factor * m2k;
                tmp2[kx]  = std::sqrt(factor * m2k);
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz * by * pme.bsp_mod[XX][kx];
                tmp1[kx]  = -factor * m2k;
                tmp2[kx]  = std::sqrt(factor * m2k);
            }
            /* Clear padding elements to avoid (harmless) fp exceptions */
            const int kxendSimd = roundUpToMultipleOfFactor<c_simdWidth>(kxend);
            for (; kx < kxendSimd; kx++)
            {
                tmp1[kx] = 0;
                tmp2[kx] = 0;
            }

            calc_exponentials_lj(
                    kxstart,
                    kxend,
                    ArrayRef<PME_T>(tmp1, tmp1 + roundUpToMultipleOfFactor<c_simdWidth>(kxend)),
                    ArrayRef<PME_T>(tmp2, tmp2 + roundUpToMultipleOfFactor<c_simdWidth>(kxend)),
                    ArrayRef<PME_T>(denom, denom + roundUpToMultipleOfFactor<c_simdWidth>(kxend)));

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2k      = factor * m2[kx];
                eterm    = -((1.0 - 2.0 * m2k) * tmp1[kx] + 2.0 * m2k * tmp2[kx]);
                tmp1[kx] = eterm * denom[kx];
            }
            const int gcount = (useLBCombinationRule ? 7 : 1);
            for (int ig = 0; ig < gcount; ++ig)
            {
                t_complex* p0;

                p0 = grids[ig].cfftgrid + iy * local_size[ZZ] * local_size[XX] + iz * local_size[XX];
                for (kx = kxstart; kx < kxend; kx++, p0++)
                {
                    d1 = p0->re;
                    d2 = p0->im;

                    eterm = tmp1[kx];

                    p0->re = d1 * eterm;
                    p0->im = d2 * eterm;
                }
            }
        }
    }
    if (computeEnergyAndVirial)
    {
        work.vir_lj[XX][XX] = 0.25 * virxx;
        work.vir_lj[YY][YY] = 0.25 * viryy;
        work.vir_lj[ZZ][ZZ] = 0.25 * virzz;
        work.vir_lj[XX][YY] = work.vir_lj[YY][XX] = 0.25 * virxy;
        work.vir_lj[XX][ZZ] = work.vir_lj[ZZ][XX] = 0.25 * virxz;
        work.vir_lj[YY][ZZ] = work.vir_lj[ZZ][YY] = 0.25 * viryz;

        /* This energy should be corrected for a charged system */
        work.energy_lj = 0.5 * energy;
    }
    /* Return the loop count over all threads */
    return local_ndata[YY] * local_ndata[ZZ] * local_ndata[XX];
}
