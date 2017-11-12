/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "pme-gather.h"

#include "gromacs/math/vec.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/typetraits.h"

#include "pme-internal.h"
#include "pme-simd.h"
#include "pme-spline-work.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

/* Spline function. Goals: 1) Force compiler to instantiate function separately
   for each compile-time value of order and once more for any possible (runtime)
   value. 2) Allow overloading for specific compile-time values.
   The Int template argument can be either int (runtime) or an object of type
   integral_constant<int, N> (compile-time). Common runtime values can be
   converted to compile-time values with a switch statement. For the compile
   value the compiler is required to instantiate the function separately for
   each value. The function can be overloaded for specific compile-time values
   using integral_constant<int, N> where N is either a specific value or an
   enable_if constrained non-type template parameter. The most specific overload
   (specific value > int template parameter > general function) is called. Inside
   the function the order argument can be used as regular int because
   integral_constant has a proper conversion.

   SIMD do_fspline() template funtions will be used for PME order 4 and 5
   when the SIMD module has support for SIMD4 for the architecture used.
   For SIMD4 without unaligned load/store support:
     order 4 and 5 use the order 4+5 aligned SIMD template
   For SIMD4 with unaligned load/store support:
     order 4 uses the order 4 unaligned SIMD template
     order 5 uses the order 4+5 aligned SIMD template
 */
struct do_fspline
{
    do_fspline (
            const gmx_pme_t *                   pme,
            const real * gmx_restrict           grid,
            const pme_atomcomm_t * gmx_restrict atc,
            const splinedata_t * gmx_restrict   spline,
            int                                 nn)
        : pme(pme), grid(grid), atc(atc), spline(spline), nn(nn) {}

    template <typename Int>
    RVec operator()(Int order) const
    {
        static_assert(isIntegralConstant<Int, int>::value || std::is_same<Int, int>::value,
                      "'order' needs to be either of type integral_constant<int,N> or int.");

        const int  norder = nn*order;

        /* Pointer arithmetic alert, next six statements */
        const real *const gmx_restrict thx  = spline->theta[XX] + norder;
        const real *const gmx_restrict thy  = spline->theta[YY] + norder;
        const real *const gmx_restrict thz  = spline->theta[ZZ] + norder;
        const real *const gmx_restrict dthx = spline->dtheta[XX] + norder;
        const real *const gmx_restrict dthy = spline->dtheta[YY] + norder;
        const real *const gmx_restrict dthz = spline->dtheta[ZZ] + norder;

        RVec                           f(0, 0, 0);

        for (int ithx = 0; (ithx < order); ithx++)
        {
            const int  index_x = (idxX + ithx)*gridNY*gridNZ;
            const real tx      = thx[ithx];
            const real dx      = dthx[ithx];

            for (int ithy = 0; (ithy < order); ithy++)
            {
                const int  index_xy = index_x + (idxY + ithy)*gridNZ;
                const real ty       = thy[ithy];
                const real dy       = dthy[ithy];
                real       fxy1     = 0, fz1 = 0;

                for (int ithz = 0; (ithz < order); ithz++)
                {
                    const real gval = grid[index_xy + (idxZ + ithz)];
                    fxy1 += thz[ithz]*gval;
                    fz1  += dthz[ithz]*gval;
                }
                f[XX] += dx*ty*fxy1;
                f[YY] += tx*dy*fxy1;
                f[ZZ] += tx*ty*fz1;
            }
        }

        return f;
    }

#ifdef PME_SIMD4_UNALIGNED //TODO: Consider always have at least a dummy implementation of Simd (enough for first phase of two-phase lookup) and then use enable_if instead of #ifdef
/* Gather for one charge with pme_order=4 with unaligned SIMD4 load+store.
 * This code does not assume any memory alignment for the grid.
 */
    RVec
    operator()(std::integral_constant<int, 4>) const
    {
        const int                      norder = nn*4;
        /* Pointer arithmetic alert, next six statements */
        const real *const gmx_restrict thx  = spline->theta[XX] + norder;
        const real *const gmx_restrict thy  = spline->theta[YY] + norder;
        const real *const gmx_restrict thz  = spline->theta[ZZ] + norder;
        const real *const gmx_restrict dthx = spline->dtheta[XX] + norder;
        const real *const gmx_restrict dthy = spline->dtheta[YY] + norder;
        const real *const gmx_restrict dthz = spline->dtheta[ZZ] + norder;

        Simd4Real                      fx_S = setZero();
        Simd4Real                      fy_S = setZero();
        Simd4Real                      fz_S = setZero();

        /* With order 4 the z-spline is actually aligned */
        const Simd4Real tz_S = load4(thz);
        const Simd4Real dz_S = load4(dthz);

        for (int ithx = 0; (ithx < 4); ithx++)
        {
            const int       index_x = (idxX + ithx)*gridNY*gridNZ;
            const Simd4Real tx_S    = Simd4Real(thx[ithx]);
            const Simd4Real dx_S    = Simd4Real(dthx[ithx]);

            for (int ithy = 0; (ithy < 4); ithy++)
            {
                const int       index_xy = index_x + (idxY + ithy)*gridNZ;
                const Simd4Real ty_S     = Simd4Real(thy[ithy]);
                const Simd4Real dy_S     = Simd4Real(dthy[ithy]);

                const Simd4Real gval_S = load4U(grid + index_xy + idxZ);

                const Simd4Real fxy1_S = tz_S * gval_S;
                const Simd4Real fz1_S  = dz_S * gval_S;

                fx_S = fma(dx_S * ty_S, fxy1_S, fx_S);
                fy_S = fma(tx_S * dy_S, fxy1_S, fy_S);
                fz_S = fma(tx_S * ty_S, fz1_S, fz_S);
            }
        }

        return {
                   reduce(fx_S), reduce(fy_S), reduce(fz_S)
        };
    }
#endif

#ifdef PME_SIMD4_SPREAD_GATHER
/* Load order elements from unaligned memory into two 4-wide SIMD */
    template<int order>
    static inline void loadOrderU(const real* data, std::integral_constant<int, order>,
                                  int offset, Simd4Real* S0, Simd4Real* S1)
    {
#ifdef PME_SIMD4_UNALIGNED //TODO: Extract into helper function
        *S0 = load4U(data-offset);
        *S1 = load4U(data-offset+4);
#else
        GMX_ALIGNED(real, GMX_SIMD4_WIDTH)  buf_aligned[GMX_SIMD4_WIDTH*2];
        /* Copy data to an aligned buffer */
        for (int i = 0; i < order; i++)
        {
            buf_aligned[offset+i]  = data[i];
        }
        *S0 = load4(buf_aligned);
        *S1 = load4(buf_aligned+4);
#endif
    }
#endif

#ifdef PME_SIMD4_SPREAD_GATHER
/* This code assumes that the grid is allocated 4-real aligned
 * and that pme->pmegrid_nz is a multiple of 4.
 * This code supports pme_order <= 5.
 */
    template <int Order>
    typename std::enable_if<Order == 4 || Order == 5, RVec>::type
    operator()(std::integral_constant<int, Order>  order) const
    {
        const int                     norder = nn*order;
        GMX_ASSERT(gridNZ % 4 == 0, "For aligned SIMD4 operations the grid size has to be padded up to a multiple of 4");
        /* Pointer arithmetic alert, next six statements */
        const real *const gmx_restrict thx  = spline->theta[XX] + norder;
        const real *const gmx_restrict thy  = spline->theta[YY] + norder;
        const real *const gmx_restrict thz  = spline->theta[ZZ] + norder;
        const real *const gmx_restrict dthx = spline->dtheta[XX] + norder;
        const real *const gmx_restrict dthy = spline->dtheta[YY] + norder;
        const real *const gmx_restrict dthz = spline->dtheta[ZZ] + norder;

        struct pme_spline_work *const  work = pme->spline_work;

        const int                      offset = idxZ & 3;

        Simd4Real                      fx_S = setZero();
        Simd4Real                      fy_S = setZero();
        Simd4Real                      fz_S = setZero();

        Simd4Real                      tz_S0, tz_S1, dz_S0, dz_S1;
        loadOrderU(thz,  order, offset, &tz_S0, &tz_S1);
        loadOrderU(dthz, order, offset, &dz_S0, &dz_S1);

        tz_S0 = selectByMask(tz_S0, work->mask_S0[offset]);
        dz_S0 = selectByMask(dz_S0, work->mask_S0[offset]);
        tz_S1 = selectByMask(tz_S1, work->mask_S1[offset]);
        dz_S1 = selectByMask(dz_S1, work->mask_S1[offset]);

        for (int ithx = 0; (ithx < order); ithx++)
        {
            const int       index_x  = (idxX + ithx)*gridNY*gridNZ;
            const Simd4Real tx_S     = Simd4Real(thx[ithx]);
            const Simd4Real dx_S     = Simd4Real(dthx[ithx]);

            for (int ithy = 0; (ithy < order); ithy++)
            {
                const int       index_xy = index_x + (idxY + ithy)*gridNZ;
                const Simd4Real ty_S     = Simd4Real(thy[ithy]);
                const Simd4Real dy_S     = Simd4Real(dthy[ithy]);

                const Simd4Real gval_S0 = load4(grid + index_xy + idxZ - offset);
                const Simd4Real gval_S1 = load4(grid + index_xy + idxZ - offset + 4);

                const Simd4Real fxy1_S0 = tz_S0 * gval_S0;
                const Simd4Real fz1_S0  = dz_S0 * gval_S0;
                const Simd4Real fxy1_S1 = tz_S1 * gval_S1;
                const Simd4Real fz1_S1  = dz_S1 * gval_S1;

                const Simd4Real fxy1_S = fxy1_S0 + fxy1_S1;
                const Simd4Real fz1_S  = fz1_S0 + fz1_S1;

                fx_S = fma(dx_S * ty_S, fxy1_S, fx_S);
                fy_S = fma(tx_S * dy_S, fxy1_S, fy_S);
                fz_S = fma(tx_S * ty_S, fz1_S, fz_S);
            }
        }

        return {
                   reduce(fx_S), reduce(fy_S), reduce(fz_S)
        };
    }
#endif
    private:
        const gmx_pme_t *const                   pme;
        const real *const gmx_restrict           grid;
        const pme_atomcomm_t *const gmx_restrict atc;
        const splinedata_t *const gmx_restrict   spline;
        const int                                nn;

        const int                                gridNY = pme->pmegrid_ny;
        const int                                gridNZ = pme->pmegrid_nz;

        const int *const                         idxptr = atc->idx[spline->ind[nn]];
        const int                                idxX   = idxptr[XX];
        const int                                idxY   = idxptr[YY];
        const int                                idxZ   = idxptr[ZZ];
};


void gather_f_bsplines(const gmx_pme_t *pme, const real *grid,
                       gmx_bool bClearF, const pme_atomcomm_t *atc,
                       const splinedata_t *spline,
                       real scale)
{
    /* sum forces for local particles */

    const int  order = pme->pme_order;
    const int  nx    = pme->nkx;
    const int  ny    = pme->nky;
    const int  nz    = pme->nkz;

    const real rxx   = pme->recipbox[XX][XX];
    const real ryx   = pme->recipbox[YY][XX];
    const real ryy   = pme->recipbox[YY][YY];
    const real rzx   = pme->recipbox[ZZ][XX];
    const real rzy   = pme->recipbox[ZZ][YY];
    const real rzz   = pme->recipbox[ZZ][ZZ];

    /* Extract the buffer for force output */
    rvec * gmx_restrict force = atc->f;

    /* Note that unrolling this loop by templating this function on order
     * deteriorates performance significantly with gcc5/6/7.
     */
    for (int nn = 0; nn < spline->n; nn++)
    {
        const int  n           = spline->ind[nn];
        const real coefficient = scale*atc->coefficient[n];

        if (bClearF)
        {
            force[n][XX] = 0;
            force[n][YY] = 0;
            force[n][ZZ] = 0;
        }
        if (coefficient != 0)
        {
            RVec       f;
            const auto spline_func = do_fspline(pme, grid, atc, spline, nn);

            switch (order)
            {
                case 4:
                    f = spline_func(std::integral_constant<int, 4>());
                    break;
                case 5:
                    f = spline_func(std::integral_constant<int, 5>());
                    break;
                default:
                    f = spline_func(order);
                    break;
            }

            force[n][XX] += -coefficient*( f[XX]*nx*rxx );
            force[n][YY] += -coefficient*( f[XX]*nx*ryx + f[YY]*ny*ryy );
            force[n][ZZ] += -coefficient*( f[XX]*nx*rzx + f[YY]*ny*rzy + f[ZZ]*nz*rzz );
        }
    }
    /* Since the energy and not forces are interpolated
     * the net force might not be exactly zero.
     * This can be solved by also interpolating F, but
     * that comes at a cost.
     * A better hack is to remove the net force every
     * step, but that must be done at a higher level
     * since this routine doesn't see all atoms if running
     * in parallel. Don't know how important it is?  EL 990726
     */
}


real gather_energy_bsplines(gmx_pme_t *pme, real *grid,
                            pme_atomcomm_t *atc)
{
    splinedata_t *spline;
    int           n, ithx, ithy, ithz, i0, j0, k0;
    int           index_x, index_xy;
    int          *idxptr;
    real          energy, pot, tx, ty, coefficient, gval;
    real         *thx, *thy, *thz;
    int           norder;
    int           order;

    spline = &atc->spline[0];

    order = pme->pme_order;

    energy = 0;
    for (n = 0; (n < atc->n); n++)
    {
        coefficient      = atc->coefficient[n];

        if (coefficient != 0)
        {
            idxptr = atc->idx[n];
            norder = n*order;

            i0   = idxptr[XX];
            j0   = idxptr[YY];
            k0   = idxptr[ZZ];

            /* Pointer arithmetic alert, next three statements */
            thx  = spline->theta[XX] + norder;
            thy  = spline->theta[YY] + norder;
            thz  = spline->theta[ZZ] + norder;

            pot = 0;
            for (ithx = 0; (ithx < order); ithx++)
            {
                index_x = (i0+ithx)*pme->pmegrid_ny*pme->pmegrid_nz;
                tx      = thx[ithx];

                for (ithy = 0; (ithy < order); ithy++)
                {
                    index_xy = index_x+(j0+ithy)*pme->pmegrid_nz;
                    ty       = thy[ithy];

                    for (ithz = 0; (ithz < order); ithz++)
                    {
                        gval  = grid[index_xy+(k0+ithz)];
                        pot  += tx*ty*thz[ithz]*gval;
                    }

                }
            }

            energy += pot*coefficient;
        }
    }

    return energy;
}
