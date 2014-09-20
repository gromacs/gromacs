/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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

#include <math.h>

#include "gromacs/legacyheaders/genborn.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/genborn_allvsall.h"
#include "gromacs/utility/smalloc.h"

#if 0 && defined (GMX_SIMD_X86_SSE2_OR_HIGHER)

#include <gmx_sse2_single.h>


#define SIMD_WIDTH 4
#define UNROLLI    4
#define UNROLLJ    4









typedef struct
{
    int *      jindex_gb;
    int **     prologue_mask_gb;
    int **     epilogue_mask;
    int *      imask;
    real *     gb_radius;
    real *     workparam;
    real *     work;
    real *     x_align;
    real *     y_align;
    real *     z_align;
    real *     fx_align;
    real *     fy_align;
    real *     fz_align;
}
gmx_allvsallgb2_data_t;


static int
calc_maxoffset(int i, int natoms)
{
    int maxoffset;

    if ((natoms % 2) == 1)
    {
        /* Odd number of atoms, easy */
        maxoffset = natoms/2;
    }
    else if ((natoms % 4) == 0)
    {
        /* Multiple of four is hard */
        if (i < natoms/2)
        {
            if ((i % 2) == 0)
            {
                maxoffset = natoms/2;
            }
            else
            {
                maxoffset = natoms/2-1;
            }
        }
        else
        {
            if ((i % 2) == 1)
            {
                maxoffset = natoms/2;
            }
            else
            {
                maxoffset = natoms/2-1;
            }
        }
    }
    else
    {
        /* natoms/2 = odd */
        if ((i % 2) == 0)
        {
            maxoffset = natoms/2;
        }
        else
        {
            maxoffset = natoms/2-1;
        }
    }

    return maxoffset;
}

static void
setup_gb_exclusions_and_indices(gmx_allvsallgb2_data_t     *   aadata,
                                t_ilist     *                  ilist,
                                int                            start,
                                int                            end,
                                int                            natoms,
                                gmx_bool                       bInclude12,
                                gmx_bool                       bInclude13,
                                gmx_bool                       bInclude14)
{
    int   i, j, k, tp;
    int   a1, a2;
    int   ni0, ni1, nj0, nj1, nj;
    int   imin, imax, iexcl;
    int   max_offset;
    int   max_excl_offset;
    int   firstinteraction;
    int   ibase;
    int  *pi;

    /* This routine can appear to be a bit complex, but it is mostly book-keeping.
     * To enable the fast all-vs-all kernel we need to be able to stream through all coordinates
     * whether they should interact or not.
     *
     * To avoid looping over the exclusions, we create a simple mask that is 1 if the interaction
     * should be present, otherwise 0. Since exclusions typically only occur when i & j are close,
     * we create a jindex array with three elements per i atom: the starting point, the point to
     * which we need to check exclusions, and the end point.
     * This way we only have to allocate a short exclusion mask per i atom.
     */

    ni0 = (start/UNROLLI)*UNROLLI;
    ni1 = ((end+UNROLLI-1)/UNROLLI)*UNROLLI;

    /* Set the interaction mask to only enable the i atoms we want to include */
    snew(pi, natoms+UNROLLI+2*SIMD_WIDTH);
    aadata->imask = (int *) (((size_t) pi + 16) & (~((size_t) 15)));
    for (i = 0; i < natoms+UNROLLI; i++)
    {
        aadata->imask[i] = (i >= start && i < end) ? 0xFFFFFFFF : 0;
    }

    /* Allocate memory for our modified jindex array */
    snew(aadata->jindex_gb, 4*(natoms+UNROLLI));
    for (i = 0; i < 4*(natoms+UNROLLI); i++)
    {
        aadata->jindex_gb[i] = 0;
    }

    /* Create the exclusion masks for the prologue part */
    snew(aadata->prologue_mask_gb, natoms+UNROLLI); /* list of pointers */

    /* First zero everything to avoid uninitialized data */
    for (i = 0; i < natoms+UNROLLI; i++)
    {
        aadata->prologue_mask_gb[i] = NULL;
    }

    /* Calculate the largest exclusion range we need for each UNROLLI-tuplet of i atoms. */
    for (ibase = ni0; ibase < ni1; ibase += UNROLLI)
    {
        max_excl_offset = -1;

        /* First find maxoffset for the next 4 atoms (or fewer if we are close to end) */
        imax = ((ibase+UNROLLI) < end) ? (ibase+UNROLLI) : end;

        /* Which atom is the first we (might) interact with? */
        imin = natoms; /* Guaranteed to be overwritten by one of 'firstinteraction' */
        for (i = ibase; i < imax; i++)
        {
            /* Before exclusions, which atom is the first we (might) interact with? */
            firstinteraction = i+1;
            max_offset       = calc_maxoffset(i, natoms);

            if (!bInclude12)
            {
                for (j = 0; j < ilist[F_GB12].nr; j += 3)
                {
                    a1 = ilist[F_GB12].iatoms[j+1];
                    a2 = ilist[F_GB12].iatoms[j+2];

                    if (a1 == i)
                    {
                        k = a2;
                    }
                    else if (a2 == i)
                    {
                        k = a1;
                    }
                    else
                    {
                        continue;
                    }

                    if (k == firstinteraction)
                    {
                        firstinteraction++;
                    }
                }
            }
            if (!bInclude13)
            {
                for (j = 0; j < ilist[F_GB13].nr; j += 3)
                {
                    a1 = ilist[F_GB13].iatoms[j+1];
                    a2 = ilist[F_GB13].iatoms[j+2];

                    if (a1 == i)
                    {
                        k = a2;
                    }
                    else if (a2 == i)
                    {
                        k = a1;
                    }
                    else
                    {
                        continue;
                    }

                    if (k == firstinteraction)
                    {
                        firstinteraction++;
                    }
                }
            }
            if (!bInclude14)
            {
                for (j = 0; j < ilist[F_GB14].nr; j += 3)
                {
                    a1 = ilist[F_GB14].iatoms[j+1];
                    a2 = ilist[F_GB14].iatoms[j+2];
                    if (a1 == i)
                    {
                        k = a2;
                    }
                    else if (a2 == i)
                    {
                        k = a1;
                    }
                    else
                    {
                        continue;
                    }

                    if (k == firstinteraction)
                    {
                        firstinteraction++;
                    }
                }
            }
            imin = (firstinteraction < imin) ? firstinteraction : imin;
        }
        /* round down to j unrolling factor */
        imin = (imin/UNROLLJ)*UNROLLJ;

        for (i = ibase; i < imax; i++)
        {
            max_offset = calc_maxoffset(i, natoms);

            if (!bInclude12)
            {
                for (j = 0; j < ilist[F_GB12].nr; j += 3)
                {
                    a1 = ilist[F_GB12].iatoms[j+1];
                    a2 = ilist[F_GB12].iatoms[j+2];

                    if (a1 == i)
                    {
                        k = a2;
                    }
                    else if (a2 == i)
                    {
                        k = a1;
                    }
                    else
                    {
                        continue;
                    }

                    if (k < imin)
                    {
                        k += natoms;
                    }

                    if (k > i+max_offset)
                    {
                        continue;
                    }

                    k = k - imin;

                    if (k+natoms <= max_offset)
                    {
                        k += natoms;
                    }
                    max_excl_offset = (k > max_excl_offset) ? k : max_excl_offset;
                }
            }
            if (!bInclude13)
            {
                for (j = 0; j < ilist[F_GB13].nr; j += 3)
                {
                    a1 = ilist[F_GB13].iatoms[j+1];
                    a2 = ilist[F_GB13].iatoms[j+2];

                    if (a1 == i)
                    {
                        k = a2;
                    }
                    else if (a2 == i)
                    {
                        k = a1;
                    }
                    else
                    {
                        continue;
                    }

                    if (k < imin)
                    {
                        k += natoms;
                    }

                    if (k > i+max_offset)
                    {
                        continue;
                    }

                    k = k - imin;

                    if (k+natoms <= max_offset)
                    {
                        k += natoms;
                    }
                    max_excl_offset = (k > max_excl_offset) ? k : max_excl_offset;
                }
            }
            if (!bInclude14)
            {
                for (j = 0; j < ilist[F_GB14].nr; j += 3)
                {
                    a1 = ilist[F_GB14].iatoms[j+1];
                    a2 = ilist[F_GB14].iatoms[j+2];

                    if (a1 == i)
                    {
                        k = a2;
                    }
                    else if (a2 == i)
                    {
                        k = a1;
                    }
                    else
                    {
                        continue;
                    }

                    if (k < imin)
                    {
                        k += natoms;
                    }

                    if (k > i+max_offset)
                    {
                        continue;
                    }

                    k = k - imin;

                    if (k+natoms <= max_offset)
                    {
                        k += natoms;
                    }
                    max_excl_offset = (k > max_excl_offset) ? k : max_excl_offset;
                }
            }
        }

        /* The offset specifies the last atom to be excluded, so add one unit to get an upper loop limit */
        max_excl_offset++;
        /* round up to j unrolling factor */
        max_excl_offset = (max_excl_offset/UNROLLJ+1)*UNROLLJ;

        /* Set all the prologue masks length to this value (even for i>end) */
        for (i = ibase; i < ibase+UNROLLI; i++)
        {
            aadata->jindex_gb[4*i]   = imin;
            aadata->jindex_gb[4*i+1] = imin+max_excl_offset;
        }
    }

    /* Now the hard part, loop over it all again to calculate the actual contents of the prologue masks */
    for (ibase = ni0; ibase < ni1; ibase += UNROLLI)
    {
        for (i = ibase; i < ibase+UNROLLI; i++)
        {
            nj   = aadata->jindex_gb[4*i+1] - aadata->jindex_gb[4*i];
            imin = aadata->jindex_gb[4*i];

            /* Allocate aligned memory */
            snew(pi, nj+2*SIMD_WIDTH);
            aadata->prologue_mask_gb[i] = (int *) (((size_t) pi + 16) & (~((size_t) 15)));

            max_offset = calc_maxoffset(i, natoms);

            /* Include interactions i+1 <= j < i+maxoffset */
            for (k = 0; k < nj; k++)
            {
                j = imin + k;

                if ( (j > i) && (j <= i+max_offset) )
                {
                    aadata->prologue_mask_gb[i][k] = 0xFFFFFFFF;
                }
                else
                {
                    aadata->prologue_mask_gb[i][k] = 0;
                }
            }

            /* Clear out the explicit exclusions */
            if (i < end)
            {
                if (!bInclude12)
                {
                    for (j = 0; j < ilist[F_GB12].nr; j += 3)
                    {
                        a1 = ilist[F_GB12].iatoms[j+1];
                        a2 = ilist[F_GB12].iatoms[j+2];

                        if (a1 == i)
                        {
                            k = a2;
                        }
                        else if (a2 == i)
                        {
                            k = a1;
                        }
                        else
                        {
                            continue;
                        }

                        if (k > i+max_offset)
                        {
                            continue;
                        }
                        k = k-i;

                        if (k+natoms <= max_offset)
                        {
                            k += natoms;
                        }

                        k = k+i-imin;
                        if (k >= 0)
                        {
                            aadata->prologue_mask_gb[i][k] = 0;
                        }
                    }
                }
                if (!bInclude13)
                {
                    for (j = 0; j < ilist[F_GB13].nr; j += 3)
                    {
                        a1 = ilist[F_GB13].iatoms[j+1];
                        a2 = ilist[F_GB13].iatoms[j+2];

                        if (a1 == i)
                        {
                            k = a2;
                        }
                        else if (a2 == i)
                        {
                            k = a1;
                        }
                        else
                        {
                            continue;
                        }

                        if (k > i+max_offset)
                        {
                            continue;
                        }
                        k = k-i;

                        if (k+natoms <= max_offset)
                        {
                            k += natoms;
                        }

                        k = k+i-imin;
                        if (k >= 0)
                        {
                            aadata->prologue_mask_gb[i][k] = 0;
                        }
                    }
                }
                if (!bInclude14)
                {
                    for (j = 0; j < ilist[F_GB14].nr; j += 3)
                    {
                        a1 = ilist[F_GB14].iatoms[j+1];
                        a2 = ilist[F_GB14].iatoms[j+2];

                        if (a1 == i)
                        {
                            k = a2;
                        }
                        else if (a2 == i)
                        {
                            k = a1;
                        }
                        else
                        {
                            continue;
                        }

                        if (k > i+max_offset)
                        {
                            continue;
                        }
                        k = k-i;

                        if (k+natoms <= max_offset)
                        {
                            k += natoms;
                        }

                        k = k+i-imin;
                        if (k >= 0)
                        {
                            aadata->prologue_mask_gb[i][k] = 0;
                        }
                    }
                }
            }
        }
    }

    /* Construct the epilogue mask - this just contains the check for maxoffset */
    snew(aadata->epilogue_mask, natoms+UNROLLI);

    /* First zero everything to avoid uninitialized data */
    for (i = 0; i < natoms+UNROLLI; i++)
    {
        aadata->jindex_gb[4*i+2]    = aadata->jindex_gb[4*i+1];
        aadata->jindex_gb[4*i+3]    = aadata->jindex_gb[4*i+1];
        aadata->epilogue_mask[i]    = NULL;
    }

    for (ibase = ni0; ibase < ni1; ibase += UNROLLI)
    {
        /* Find the lowest index for which we need to use the epilogue */
        imin       = ibase;
        max_offset = calc_maxoffset(imin, natoms);

        imin = imin + 1 + max_offset;

        /* Find largest index for which we need to use the epilogue */
        imax = ibase + UNROLLI-1;
        imax = (imax < end) ? imax : end;

        max_offset = calc_maxoffset(imax, natoms);
        imax       = imax + 1 + max_offset + UNROLLJ - 1;

        for (i = ibase; i < ibase+UNROLLI; i++)
        {
            /* Start of epilogue - round down to j tile limit */
            aadata->jindex_gb[4*i+2] = (imin/UNROLLJ)*UNROLLJ;
            /* Make sure we dont overlap - for small systems everything is done in the prologue */
            aadata->jindex_gb[4*i+2] = (aadata->jindex_gb[4*i+1] > aadata->jindex_gb[4*i+2]) ? aadata->jindex_gb[4*i+1] : aadata->jindex_gb[4*i+2];
            /* Round upwards to j tile limit */
            aadata->jindex_gb[4*i+3] = (imax/UNROLLJ)*UNROLLJ;
            /* Make sure we dont have a negative range for the epilogue */
            aadata->jindex_gb[4*i+3] = (aadata->jindex_gb[4*i+2] > aadata->jindex_gb[4*i+3]) ? aadata->jindex_gb[4*i+2] : aadata->jindex_gb[4*i+3];
        }
    }

    /* And fill it with data... */

    for (ibase = ni0; ibase < ni1; ibase += UNROLLI)
    {
        for (i = ibase; i < ibase+UNROLLI; i++)
        {

            nj = aadata->jindex_gb[4*i+3] - aadata->jindex_gb[4*i+2];

            /* Allocate aligned memory */
            snew(pi, nj+2*SIMD_WIDTH);
            aadata->epilogue_mask[i] = (int *) (((size_t) pi + 16) & (~((size_t) 15)));

            max_offset = calc_maxoffset(i, natoms);

            for (k = 0; k < nj; k++)
            {
                j = aadata->jindex_gb[4*i+2] + k;
                aadata->epilogue_mask[i][k] = (j <= i+max_offset) ? 0xFFFFFFFF : 0;
            }
        }
    }
}


static void
genborn_allvsall_setup(gmx_allvsallgb2_data_t     **  p_aadata,
                       gmx_localtop_t     *           top,
                       gmx_genborn_t     *            born,
                       t_mdatoms     *                mdatoms,
                       real                           radius_offset,
                       int                            gb_algorithm,
                       gmx_bool                       bInclude12,
                       gmx_bool                       bInclude13,
                       gmx_bool                       bInclude14)
{
    int                     i, j, idx;
    int                     natoms;
    gmx_allvsallgb2_data_t *aadata;
    real                   *p;

    natoms = mdatoms->nr;

    snew(aadata, 1);
    *p_aadata = aadata;

    snew(p, 2*natoms+2*SIMD_WIDTH);
    aadata->x_align = (real *) (((size_t) p + 16) & (~((size_t) 15)));
    snew(p, 2*natoms+2*SIMD_WIDTH);
    aadata->y_align = (real *) (((size_t) p + 16) & (~((size_t) 15)));
    snew(p, 2*natoms+2*SIMD_WIDTH);
    aadata->z_align = (real *) (((size_t) p + 16) & (~((size_t) 15)));
    snew(p, 2*natoms+2*SIMD_WIDTH);
    aadata->fx_align = (real *) (((size_t) p + 16) & (~((size_t) 15)));
    snew(p, 2*natoms+2*SIMD_WIDTH);
    aadata->fy_align = (real *) (((size_t) p + 16) & (~((size_t) 15)));
    snew(p, 2*natoms+2*SIMD_WIDTH);
    aadata->fz_align = (real *) (((size_t) p + 16) & (~((size_t) 15)));

    snew(p, 2*natoms+UNROLLJ+SIMD_WIDTH);
    aadata->gb_radius = (real *) (((size_t) p + 16) & (~((size_t) 15)));

    snew(p, 2*natoms+UNROLLJ+SIMD_WIDTH);
    aadata->workparam = (real *) (((size_t) p + 16) & (~((size_t) 15)));

    snew(p, 2*natoms+UNROLLJ+SIMD_WIDTH);
    aadata->work = (real *) (((size_t) p + 16) & (~((size_t) 15)));

    for (i = 0; i < mdatoms->nr; i++)
    {
        aadata->gb_radius[i] = top->atomtypes.gb_radius[mdatoms->typeA[i]] - radius_offset;
        if (gb_algorithm == egbSTILL)
        {
            aadata->workparam[i] = born->vsolv[i];
        }
        else if (gb_algorithm == egbOBC)
        {
            aadata->workparam[i] = born->param[i];
        }
        aadata->work[i]      = 0.0;
    }
    for (i = 0; i < mdatoms->nr; i++)
    {
        aadata->gb_radius[natoms+i] = aadata->gb_radius[i];
        aadata->workparam[natoms+i] = aadata->workparam[i];
        aadata->work[natoms+i]      = aadata->work[i];
    }

    for (i = 0; i < 2*natoms+SIMD_WIDTH; i++)
    {
        aadata->x_align[i]  = 0.0;
        aadata->y_align[i]  = 0.0;
        aadata->z_align[i]  = 0.0;
        aadata->fx_align[i] = 0.0;
        aadata->fy_align[i] = 0.0;
        aadata->fz_align[i] = 0.0;
    }

    setup_gb_exclusions_and_indices(aadata, top->idef.il, 0, mdatoms->homenr, mdatoms->nr,
                                    bInclude12, bInclude13, bInclude14);
}


int
genborn_allvsall_calc_still_radii_sse2_single(t_forcerec *           fr,
                                              t_mdatoms *            mdatoms,
                                              gmx_genborn_t *        born,
                                              gmx_localtop_t *       top,
                                              real *                 x,
                                              t_commrec *            cr,
                                              void *                 paadata)
{
    gmx_allvsallgb2_data_t *aadata;
    int                     natoms;
    int                     ni0, ni1;
    int                     nj0, nj1, nj2, nj3;
    int                     i, j, k, n;
    int              *      mask;
    int              *      pmask0;
    int              *      pmask1;
    int              *      pmask2;
    int              *      pmask3;
    int              *      emask0;
    int              *      emask1;
    int              *      emask2;
    int              *      emask3;
    real                    ix, iy, iz;
    real                    jx, jy, jz;
    real                    dx, dy, dz;
    real                    rsq, rinv;
    real                    gpi, rai, vai;
    real                    prod_ai;
    real                    irsq, idr4, idr6;
    real                    raj, rvdw, ratio;
    real                    vaj, ccf, dccf, theta, cosq;
    real                    term, prod, icf4, icf6, gpi2, factor, sinq;
    real              *     gb_radius;
    real              *     vsolv;
    real              *     work;
    real                    tmpsum[4];
    real              *     x_align;
    real              *     y_align;
    real              *     z_align;
    int              *      jindex;
    real              *     dadx;

    __m128                  ix_SSE0, iy_SSE0, iz_SSE0;
    __m128                  ix_SSE1, iy_SSE1, iz_SSE1;
    __m128                  ix_SSE2, iy_SSE2, iz_SSE2;
    __m128                  ix_SSE3, iy_SSE3, iz_SSE3;
    __m128                  gpi_SSE0, rai_SSE0, prod_ai_SSE0;
    __m128                  gpi_SSE1, rai_SSE1, prod_ai_SSE1;
    __m128                  gpi_SSE2, rai_SSE2, prod_ai_SSE2;
    __m128                  gpi_SSE3, rai_SSE3, prod_ai_SSE3;
    __m128                  imask_SSE0, jmask_SSE0;
    __m128                  imask_SSE1, jmask_SSE1;
    __m128                  imask_SSE2, jmask_SSE2;
    __m128                  imask_SSE3, jmask_SSE3;
    __m128                  jx_SSE, jy_SSE, jz_SSE;
    __m128                  dx_SSE0, dy_SSE0, dz_SSE0;
    __m128                  dx_SSE1, dy_SSE1, dz_SSE1;
    __m128                  dx_SSE2, dy_SSE2, dz_SSE2;
    __m128                  dx_SSE3, dy_SSE3, dz_SSE3;
    __m128                  rsq_SSE0, rinv_SSE0, irsq_SSE0, idr4_SSE0, idr6_SSE0;
    __m128                  rsq_SSE1, rinv_SSE1, irsq_SSE1, idr4_SSE1, idr6_SSE1;
    __m128                  rsq_SSE2, rinv_SSE2, irsq_SSE2, idr4_SSE2, idr6_SSE2;
    __m128                  rsq_SSE3, rinv_SSE3, irsq_SSE3, idr4_SSE3, idr6_SSE3;
    __m128                  raj_SSE, vaj_SSE, prod_SSE;
    __m128                  rvdw_SSE0, ratio_SSE0;
    __m128                  rvdw_SSE1, ratio_SSE1;
    __m128                  rvdw_SSE2, ratio_SSE2;
    __m128                  rvdw_SSE3, ratio_SSE3;
    __m128                  theta_SSE0, sinq_SSE0, cosq_SSE0, term_SSE0;
    __m128                  theta_SSE1, sinq_SSE1, cosq_SSE1, term_SSE1;
    __m128                  theta_SSE2, sinq_SSE2, cosq_SSE2, term_SSE2;
    __m128                  theta_SSE3, sinq_SSE3, cosq_SSE3, term_SSE3;
    __m128                  ccf_SSE0, dccf_SSE0;
    __m128                  ccf_SSE1, dccf_SSE1;
    __m128                  ccf_SSE2, dccf_SSE2;
    __m128                  ccf_SSE3, dccf_SSE3;
    __m128                  icf4_SSE0, icf6_SSE0;
    __m128                  icf4_SSE1, icf6_SSE1;
    __m128                  icf4_SSE2, icf6_SSE2;
    __m128                  icf4_SSE3, icf6_SSE3;
    __m128                  half_SSE, one_SSE, two_SSE, four_SSE;
    __m128                  still_p4_SSE, still_p5inv_SSE, still_pip5_SSE;

    natoms              = mdatoms->nr;
    ni0                 = 0;
    ni1                 = mdatoms->homenr;

    n = 0;

    aadata = *((gmx_allvsallgb2_data_t **)paadata);


    if (aadata == NULL)
    {
        genborn_allvsall_setup(&aadata, top, born, mdatoms, 0.0,
                               egbSTILL, FALSE, FALSE, TRUE);
        *((gmx_allvsallgb2_data_t **)paadata) = aadata;
    }

    x_align = aadata->x_align;
    y_align = aadata->y_align;
    z_align = aadata->z_align;

    gb_radius = aadata->gb_radius;
    vsolv     = aadata->workparam;
    work      = aadata->work;
    jindex    = aadata->jindex_gb;
    dadx      = fr->dadx;

    still_p4_SSE    = _mm_set1_ps(STILL_P4);
    still_p5inv_SSE = _mm_set1_ps(STILL_P5INV);
    still_pip5_SSE  = _mm_set1_ps(STILL_PIP5);
    half_SSE        = _mm_set1_ps(0.5);
    one_SSE         = _mm_set1_ps(1.0);
    two_SSE         = _mm_set1_ps(2.0);
    four_SSE        = _mm_set1_ps(4.0);

    /* This will be summed, so it has to extend to natoms + buffer */
    for (i = 0; i < natoms+1+natoms/2; i++)
    {
        work[i] = 0;
    }

    for (i = ni0; i < ni1+1+natoms/2; i++)
    {
        k           = i%natoms;
        x_align[i]  = x[3*k];
        y_align[i]  = x[3*k+1];
        z_align[i]  = x[3*k+2];
        work[i]     = 0;
    }


    for (i = ni0; i < ni1; i += UNROLLI)
    {
        /* We assume shifts are NOT used for all-vs-all interactions */

        /* Load i atom data */
        ix_SSE0          = _mm_load1_ps(x_align+i);
        iy_SSE0          = _mm_load1_ps(y_align+i);
        iz_SSE0          = _mm_load1_ps(z_align+i);
        ix_SSE1          = _mm_load1_ps(x_align+i+1);
        iy_SSE1          = _mm_load1_ps(y_align+i+1);
        iz_SSE1          = _mm_load1_ps(z_align+i+1);
        ix_SSE2          = _mm_load1_ps(x_align+i+2);
        iy_SSE2          = _mm_load1_ps(y_align+i+2);
        iz_SSE2          = _mm_load1_ps(z_align+i+2);
        ix_SSE3          = _mm_load1_ps(x_align+i+3);
        iy_SSE3          = _mm_load1_ps(y_align+i+3);
        iz_SSE3          = _mm_load1_ps(z_align+i+3);

        gpi_SSE0         = _mm_setzero_ps();
        gpi_SSE1         = _mm_setzero_ps();
        gpi_SSE2         = _mm_setzero_ps();
        gpi_SSE3         = _mm_setzero_ps();

        rai_SSE0         = _mm_load1_ps(gb_radius+i);
        rai_SSE1         = _mm_load1_ps(gb_radius+i+1);
        rai_SSE2         = _mm_load1_ps(gb_radius+i+2);
        rai_SSE3         = _mm_load1_ps(gb_radius+i+3);

        prod_ai_SSE0     = _mm_set1_ps(STILL_P4*vsolv[i]);
        prod_ai_SSE1     = _mm_set1_ps(STILL_P4*vsolv[i+1]);
        prod_ai_SSE2     = _mm_set1_ps(STILL_P4*vsolv[i+2]);
        prod_ai_SSE3     = _mm_set1_ps(STILL_P4*vsolv[i+3]);

        /* Load limits for loop over neighbors */
        nj0              = jindex[4*i];
        nj1              = jindex[4*i+1];
        nj2              = jindex[4*i+2];
        nj3              = jindex[4*i+3];

        pmask0           = aadata->prologue_mask_gb[i];
        pmask1           = aadata->prologue_mask_gb[i+1];
        pmask2           = aadata->prologue_mask_gb[i+2];
        pmask3           = aadata->prologue_mask_gb[i+3];
        emask0           = aadata->epilogue_mask[i];
        emask1           = aadata->epilogue_mask[i+1];
        emask2           = aadata->epilogue_mask[i+2];
        emask3           = aadata->epilogue_mask[i+3];

        imask_SSE0        = _mm_load1_ps((real *)(aadata->imask+i));
        imask_SSE1        = _mm_load1_ps((real *)(aadata->imask+i+1));
        imask_SSE2        = _mm_load1_ps((real *)(aadata->imask+i+2));
        imask_SSE3        = _mm_load1_ps((real *)(aadata->imask+i+3));

        /* Prologue part, including exclusion mask */
        for (j = nj0; j < nj1; j += UNROLLJ)
        {
            jmask_SSE0 = _mm_load_ps((real *)pmask0);
            jmask_SSE1 = _mm_load_ps((real *)pmask1);
            jmask_SSE2 = _mm_load_ps((real *)pmask2);
            jmask_SSE3 = _mm_load_ps((real *)pmask3);
            pmask0    += UNROLLJ;
            pmask1    += UNROLLJ;
            pmask2    += UNROLLJ;
            pmask3    += UNROLLJ;

            /* load j atom coordinates */
            jx_SSE            = _mm_load_ps(x_align+j);
            jy_SSE            = _mm_load_ps(y_align+j);
            jz_SSE            = _mm_load_ps(z_align+j);

            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0, jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0, jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0, jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1, jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1, jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1, jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2, jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2, jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2, jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3, jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3, jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1, dy_SSE1, dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2, dy_SSE2, dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3, dy_SSE3, dz_SSE3);

            /* Combine masks */
            jmask_SSE0         = _mm_and_ps(jmask_SSE0, imask_SSE0);
            jmask_SSE1         = _mm_and_ps(jmask_SSE1, imask_SSE1);
            jmask_SSE2         = _mm_and_ps(jmask_SSE2, imask_SSE2);
            jmask_SSE3         = _mm_and_ps(jmask_SSE3, imask_SSE3);

            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);

            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0, jmask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1, jmask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2, jmask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3, jmask_SSE3);

            irsq_SSE0          = _mm_mul_ps(rinv_SSE0, rinv_SSE0);
            irsq_SSE1          = _mm_mul_ps(rinv_SSE1, rinv_SSE1);
            irsq_SSE2          = _mm_mul_ps(rinv_SSE2, rinv_SSE2);
            irsq_SSE3          = _mm_mul_ps(rinv_SSE3, rinv_SSE3);
            idr4_SSE0          = _mm_mul_ps(irsq_SSE0, irsq_SSE0);
            idr4_SSE1          = _mm_mul_ps(irsq_SSE1, irsq_SSE1);
            idr4_SSE2          = _mm_mul_ps(irsq_SSE2, irsq_SSE2);
            idr4_SSE3          = _mm_mul_ps(irsq_SSE3, irsq_SSE3);
            idr6_SSE0          = _mm_mul_ps(idr4_SSE0, irsq_SSE0);
            idr6_SSE1          = _mm_mul_ps(idr4_SSE1, irsq_SSE1);
            idr6_SSE2          = _mm_mul_ps(idr4_SSE2, irsq_SSE2);
            idr6_SSE3          = _mm_mul_ps(idr4_SSE3, irsq_SSE3);

            raj_SSE            = _mm_load_ps(gb_radius+j);
            vaj_SSE            = _mm_load_ps(vsolv+j);

            rvdw_SSE0          = _mm_add_ps(rai_SSE0, raj_SSE);
            rvdw_SSE1          = _mm_add_ps(rai_SSE1, raj_SSE);
            rvdw_SSE2          = _mm_add_ps(rai_SSE2, raj_SSE);
            rvdw_SSE3          = _mm_add_ps(rai_SSE3, raj_SSE);

            ratio_SSE0         = _mm_mul_ps(rsq_SSE0, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE0, rvdw_SSE0)));
            ratio_SSE1         = _mm_mul_ps(rsq_SSE1, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE1, rvdw_SSE1)));
            ratio_SSE2         = _mm_mul_ps(rsq_SSE2, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE2, rvdw_SSE2)));
            ratio_SSE3         = _mm_mul_ps(rsq_SSE3, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE3, rvdw_SSE3)));

            ratio_SSE0         = _mm_min_ps(ratio_SSE0, still_p5inv_SSE);
            ratio_SSE1         = _mm_min_ps(ratio_SSE1, still_p5inv_SSE);
            ratio_SSE2         = _mm_min_ps(ratio_SSE2, still_p5inv_SSE);
            ratio_SSE3         = _mm_min_ps(ratio_SSE3, still_p5inv_SSE);
            theta_SSE0         = _mm_mul_ps(ratio_SSE0, still_pip5_SSE);
            theta_SSE1         = _mm_mul_ps(ratio_SSE1, still_pip5_SSE);
            theta_SSE2         = _mm_mul_ps(ratio_SSE2, still_pip5_SSE);
            theta_SSE3         = _mm_mul_ps(ratio_SSE3, still_pip5_SSE);
            gmx_mm_sincos_ps(theta_SSE0, &sinq_SSE0, &cosq_SSE0);
            gmx_mm_sincos_ps(theta_SSE1, &sinq_SSE1, &cosq_SSE1);
            gmx_mm_sincos_ps(theta_SSE2, &sinq_SSE2, &cosq_SSE2);
            gmx_mm_sincos_ps(theta_SSE3, &sinq_SSE3, &cosq_SSE3);
            term_SSE0          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE0));
            term_SSE1          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE1));
            term_SSE2          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE2));
            term_SSE3          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE3));
            ccf_SSE0           = _mm_mul_ps(term_SSE0, term_SSE0);
            ccf_SSE1           = _mm_mul_ps(term_SSE1, term_SSE1);
            ccf_SSE2           = _mm_mul_ps(term_SSE2, term_SSE2);
            ccf_SSE3           = _mm_mul_ps(term_SSE3, term_SSE3);
            dccf_SSE0          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE0),
                                            _mm_mul_ps(sinq_SSE0, theta_SSE0));
            dccf_SSE1          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE1),
                                            _mm_mul_ps(sinq_SSE1, theta_SSE1));
            dccf_SSE2          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE2),
                                            _mm_mul_ps(sinq_SSE2, theta_SSE2));
            dccf_SSE3          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE3),
                                            _mm_mul_ps(sinq_SSE3, theta_SSE3));

            prod_SSE           = _mm_mul_ps(still_p4_SSE, vaj_SSE);
            icf4_SSE0          = _mm_mul_ps(ccf_SSE0, idr4_SSE0);
            icf4_SSE1          = _mm_mul_ps(ccf_SSE1, idr4_SSE1);
            icf4_SSE2          = _mm_mul_ps(ccf_SSE2, idr4_SSE2);
            icf4_SSE3          = _mm_mul_ps(ccf_SSE3, idr4_SSE3);
            icf6_SSE0          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE0), dccf_SSE0), idr6_SSE0);
            icf6_SSE1          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE1), dccf_SSE1), idr6_SSE1);
            icf6_SSE2          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE2), dccf_SSE2), idr6_SSE2);
            icf6_SSE3          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE3), dccf_SSE3), idr6_SSE3);

            _mm_store_ps(work+j, _mm_add_ps(_mm_load_ps(work+j),
                                            gmx_mm_sum4_ps(_mm_mul_ps(prod_ai_SSE0, icf4_SSE0),
                                                           _mm_mul_ps(prod_ai_SSE1, icf4_SSE1),
                                                           _mm_mul_ps(prod_ai_SSE2, icf4_SSE2),
                                                           _mm_mul_ps(prod_ai_SSE3, icf4_SSE3))));

            gpi_SSE0           = _mm_add_ps(gpi_SSE0, _mm_mul_ps(prod_SSE, icf4_SSE0));
            gpi_SSE1           = _mm_add_ps(gpi_SSE1, _mm_mul_ps(prod_SSE, icf4_SSE1));
            gpi_SSE2           = _mm_add_ps(gpi_SSE2, _mm_mul_ps(prod_SSE, icf4_SSE2));
            gpi_SSE3           = _mm_add_ps(gpi_SSE3, _mm_mul_ps(prod_SSE, icf4_SSE3));

            /* Save ai->aj and aj->ai chain rule terms */
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE3));
            dadx += 4;

            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE0, icf6_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE1, icf6_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE2, icf6_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE3, icf6_SSE3));
            dadx += 4;
        }

        /* Main part, no exclusions */
        for (j = nj1; j < nj2; j += UNROLLJ)
        {
            /* load j atom coordinates */
            jx_SSE            = _mm_load_ps(x_align+j);
            jy_SSE            = _mm_load_ps(y_align+j);
            jz_SSE            = _mm_load_ps(z_align+j);

            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0, jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0, jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0, jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1, jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1, jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1, jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2, jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2, jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2, jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3, jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3, jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1, dy_SSE1, dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2, dy_SSE2, dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3, dy_SSE3, dz_SSE3);

            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);

            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0, imask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1, imask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2, imask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3, imask_SSE3);

            irsq_SSE0          = _mm_mul_ps(rinv_SSE0, rinv_SSE0);
            irsq_SSE1          = _mm_mul_ps(rinv_SSE1, rinv_SSE1);
            irsq_SSE2          = _mm_mul_ps(rinv_SSE2, rinv_SSE2);
            irsq_SSE3          = _mm_mul_ps(rinv_SSE3, rinv_SSE3);
            idr4_SSE0          = _mm_mul_ps(irsq_SSE0, irsq_SSE0);
            idr4_SSE1          = _mm_mul_ps(irsq_SSE1, irsq_SSE1);
            idr4_SSE2          = _mm_mul_ps(irsq_SSE2, irsq_SSE2);
            idr4_SSE3          = _mm_mul_ps(irsq_SSE3, irsq_SSE3);
            idr6_SSE0          = _mm_mul_ps(idr4_SSE0, irsq_SSE0);
            idr6_SSE1          = _mm_mul_ps(idr4_SSE1, irsq_SSE1);
            idr6_SSE2          = _mm_mul_ps(idr4_SSE2, irsq_SSE2);
            idr6_SSE3          = _mm_mul_ps(idr4_SSE3, irsq_SSE3);

            raj_SSE            = _mm_load_ps(gb_radius+j);

            rvdw_SSE0          = _mm_add_ps(rai_SSE0, raj_SSE);
            rvdw_SSE1          = _mm_add_ps(rai_SSE1, raj_SSE);
            rvdw_SSE2          = _mm_add_ps(rai_SSE2, raj_SSE);
            rvdw_SSE3          = _mm_add_ps(rai_SSE3, raj_SSE);
            vaj_SSE            = _mm_load_ps(vsolv+j);

            ratio_SSE0         = _mm_mul_ps(rsq_SSE0, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE0, rvdw_SSE0)));
            ratio_SSE1         = _mm_mul_ps(rsq_SSE1, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE1, rvdw_SSE1)));
            ratio_SSE2         = _mm_mul_ps(rsq_SSE2, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE2, rvdw_SSE2)));
            ratio_SSE3         = _mm_mul_ps(rsq_SSE3, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE3, rvdw_SSE3)));

            ratio_SSE0         = _mm_min_ps(ratio_SSE0, still_p5inv_SSE);
            ratio_SSE1         = _mm_min_ps(ratio_SSE1, still_p5inv_SSE);
            ratio_SSE2         = _mm_min_ps(ratio_SSE2, still_p5inv_SSE);
            ratio_SSE3         = _mm_min_ps(ratio_SSE3, still_p5inv_SSE);
            theta_SSE0         = _mm_mul_ps(ratio_SSE0, still_pip5_SSE);
            theta_SSE1         = _mm_mul_ps(ratio_SSE1, still_pip5_SSE);
            theta_SSE2         = _mm_mul_ps(ratio_SSE2, still_pip5_SSE);
            theta_SSE3         = _mm_mul_ps(ratio_SSE3, still_pip5_SSE);
            gmx_mm_sincos_ps(theta_SSE0, &sinq_SSE0, &cosq_SSE0);
            gmx_mm_sincos_ps(theta_SSE1, &sinq_SSE1, &cosq_SSE1);
            gmx_mm_sincos_ps(theta_SSE2, &sinq_SSE2, &cosq_SSE2);
            gmx_mm_sincos_ps(theta_SSE3, &sinq_SSE3, &cosq_SSE3);
            term_SSE0          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE0));
            term_SSE1          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE1));
            term_SSE2          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE2));
            term_SSE3          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE3));
            ccf_SSE0           = _mm_mul_ps(term_SSE0, term_SSE0);
            ccf_SSE1           = _mm_mul_ps(term_SSE1, term_SSE1);
            ccf_SSE2           = _mm_mul_ps(term_SSE2, term_SSE2);
            ccf_SSE3           = _mm_mul_ps(term_SSE3, term_SSE3);
            dccf_SSE0          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE0),
                                            _mm_mul_ps(sinq_SSE0, theta_SSE0));
            dccf_SSE1          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE1),
                                            _mm_mul_ps(sinq_SSE1, theta_SSE1));
            dccf_SSE2          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE2),
                                            _mm_mul_ps(sinq_SSE2, theta_SSE2));
            dccf_SSE3          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE3),
                                            _mm_mul_ps(sinq_SSE3, theta_SSE3));

            prod_SSE           = _mm_mul_ps(still_p4_SSE, vaj_SSE );
            icf4_SSE0          = _mm_mul_ps(ccf_SSE0, idr4_SSE0);
            icf4_SSE1          = _mm_mul_ps(ccf_SSE1, idr4_SSE1);
            icf4_SSE2          = _mm_mul_ps(ccf_SSE2, idr4_SSE2);
            icf4_SSE3          = _mm_mul_ps(ccf_SSE3, idr4_SSE3);
            icf6_SSE0          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE0), dccf_SSE0), idr6_SSE0);
            icf6_SSE1          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE1), dccf_SSE1), idr6_SSE1);
            icf6_SSE2          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE2), dccf_SSE2), idr6_SSE2);
            icf6_SSE3          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE3), dccf_SSE3), idr6_SSE3);

            _mm_store_ps(work+j, _mm_add_ps(_mm_load_ps(work+j),
                                            gmx_mm_sum4_ps(_mm_mul_ps(prod_ai_SSE0, icf4_SSE0),
                                                           _mm_mul_ps(prod_ai_SSE1, icf4_SSE1),
                                                           _mm_mul_ps(prod_ai_SSE2, icf4_SSE2),
                                                           _mm_mul_ps(prod_ai_SSE3, icf4_SSE3))));

            gpi_SSE0           = _mm_add_ps(gpi_SSE0, _mm_mul_ps(prod_SSE, icf4_SSE0));
            gpi_SSE1           = _mm_add_ps(gpi_SSE1, _mm_mul_ps(prod_SSE, icf4_SSE1));
            gpi_SSE2           = _mm_add_ps(gpi_SSE2, _mm_mul_ps(prod_SSE, icf4_SSE2));
            gpi_SSE3           = _mm_add_ps(gpi_SSE3, _mm_mul_ps(prod_SSE, icf4_SSE3));

            /* Save ai->aj and aj->ai chain rule terms */
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE3));
            dadx += 4;

            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE0, icf6_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE1, icf6_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE2, icf6_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE3, icf6_SSE3));
            dadx += 4;
        }
        /* Epilogue part, including exclusion mask */
        for (j = nj2; j < nj3; j += UNROLLJ)
        {
            jmask_SSE0 = _mm_load_ps((real *)emask0);
            jmask_SSE1 = _mm_load_ps((real *)emask1);
            jmask_SSE2 = _mm_load_ps((real *)emask2);
            jmask_SSE3 = _mm_load_ps((real *)emask3);
            emask0    += UNROLLJ;
            emask1    += UNROLLJ;
            emask2    += UNROLLJ;
            emask3    += UNROLLJ;

            /* load j atom coordinates */
            jx_SSE            = _mm_load_ps(x_align+j);
            jy_SSE            = _mm_load_ps(y_align+j);
            jz_SSE            = _mm_load_ps(z_align+j);

            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0, jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0, jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0, jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1, jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1, jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1, jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2, jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2, jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2, jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3, jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3, jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1, dy_SSE1, dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2, dy_SSE2, dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3, dy_SSE3, dz_SSE3);

            /* Combine masks */
            jmask_SSE0         = _mm_and_ps(jmask_SSE0, imask_SSE0);
            jmask_SSE1         = _mm_and_ps(jmask_SSE1, imask_SSE1);
            jmask_SSE2         = _mm_and_ps(jmask_SSE2, imask_SSE2);
            jmask_SSE3         = _mm_and_ps(jmask_SSE3, imask_SSE3);

            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);

            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0, jmask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1, jmask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2, jmask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3, jmask_SSE3);

            irsq_SSE0          = _mm_mul_ps(rinv_SSE0, rinv_SSE0);
            irsq_SSE1          = _mm_mul_ps(rinv_SSE1, rinv_SSE1);
            irsq_SSE2          = _mm_mul_ps(rinv_SSE2, rinv_SSE2);
            irsq_SSE3          = _mm_mul_ps(rinv_SSE3, rinv_SSE3);
            idr4_SSE0          = _mm_mul_ps(irsq_SSE0, irsq_SSE0);
            idr4_SSE1          = _mm_mul_ps(irsq_SSE1, irsq_SSE1);
            idr4_SSE2          = _mm_mul_ps(irsq_SSE2, irsq_SSE2);
            idr4_SSE3          = _mm_mul_ps(irsq_SSE3, irsq_SSE3);
            idr6_SSE0          = _mm_mul_ps(idr4_SSE0, irsq_SSE0);
            idr6_SSE1          = _mm_mul_ps(idr4_SSE1, irsq_SSE1);
            idr6_SSE2          = _mm_mul_ps(idr4_SSE2, irsq_SSE2);
            idr6_SSE3          = _mm_mul_ps(idr4_SSE3, irsq_SSE3);

            raj_SSE            = _mm_load_ps(gb_radius+j);
            vaj_SSE            = _mm_load_ps(vsolv+j);

            rvdw_SSE0          = _mm_add_ps(rai_SSE0, raj_SSE);
            rvdw_SSE1          = _mm_add_ps(rai_SSE1, raj_SSE);
            rvdw_SSE2          = _mm_add_ps(rai_SSE2, raj_SSE);
            rvdw_SSE3          = _mm_add_ps(rai_SSE3, raj_SSE);

            ratio_SSE0         = _mm_mul_ps(rsq_SSE0, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE0, rvdw_SSE0)));
            ratio_SSE1         = _mm_mul_ps(rsq_SSE1, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE1, rvdw_SSE1)));
            ratio_SSE2         = _mm_mul_ps(rsq_SSE2, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE2, rvdw_SSE2)));
            ratio_SSE3         = _mm_mul_ps(rsq_SSE3, gmx_mm_inv_ps( _mm_mul_ps(rvdw_SSE3, rvdw_SSE3)));

            ratio_SSE0         = _mm_min_ps(ratio_SSE0, still_p5inv_SSE);
            ratio_SSE1         = _mm_min_ps(ratio_SSE1, still_p5inv_SSE);
            ratio_SSE2         = _mm_min_ps(ratio_SSE2, still_p5inv_SSE);
            ratio_SSE3         = _mm_min_ps(ratio_SSE3, still_p5inv_SSE);
            theta_SSE0         = _mm_mul_ps(ratio_SSE0, still_pip5_SSE);
            theta_SSE1         = _mm_mul_ps(ratio_SSE1, still_pip5_SSE);
            theta_SSE2         = _mm_mul_ps(ratio_SSE2, still_pip5_SSE);
            theta_SSE3         = _mm_mul_ps(ratio_SSE3, still_pip5_SSE);
            gmx_mm_sincos_ps(theta_SSE0, &sinq_SSE0, &cosq_SSE0);
            gmx_mm_sincos_ps(theta_SSE1, &sinq_SSE1, &cosq_SSE1);
            gmx_mm_sincos_ps(theta_SSE2, &sinq_SSE2, &cosq_SSE2);
            gmx_mm_sincos_ps(theta_SSE3, &sinq_SSE3, &cosq_SSE3);
            term_SSE0          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE0));
            term_SSE1          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE1));
            term_SSE2          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE2));
            term_SSE3          = _mm_mul_ps(half_SSE, _mm_sub_ps(one_SSE, cosq_SSE3));
            ccf_SSE0           = _mm_mul_ps(term_SSE0, term_SSE0);
            ccf_SSE1           = _mm_mul_ps(term_SSE1, term_SSE1);
            ccf_SSE2           = _mm_mul_ps(term_SSE2, term_SSE2);
            ccf_SSE3           = _mm_mul_ps(term_SSE3, term_SSE3);
            dccf_SSE0          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE0),
                                            _mm_mul_ps(sinq_SSE0, theta_SSE0));
            dccf_SSE1          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE1),
                                            _mm_mul_ps(sinq_SSE1, theta_SSE1));
            dccf_SSE2          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE2),
                                            _mm_mul_ps(sinq_SSE2, theta_SSE2));
            dccf_SSE3          = _mm_mul_ps(_mm_mul_ps(two_SSE, term_SSE3),
                                            _mm_mul_ps(sinq_SSE3, theta_SSE3));

            prod_SSE           = _mm_mul_ps(still_p4_SSE, vaj_SSE);
            icf4_SSE0          = _mm_mul_ps(ccf_SSE0, idr4_SSE0);
            icf4_SSE1          = _mm_mul_ps(ccf_SSE1, idr4_SSE1);
            icf4_SSE2          = _mm_mul_ps(ccf_SSE2, idr4_SSE2);
            icf4_SSE3          = _mm_mul_ps(ccf_SSE3, idr4_SSE3);
            icf6_SSE0          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE0), dccf_SSE0), idr6_SSE0);
            icf6_SSE1          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE1), dccf_SSE1), idr6_SSE1);
            icf6_SSE2          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE2), dccf_SSE2), idr6_SSE2);
            icf6_SSE3          = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four_SSE, ccf_SSE3), dccf_SSE3), idr6_SSE3);

            _mm_store_ps(work+j, _mm_add_ps(_mm_load_ps(work+j),
                                            gmx_mm_sum4_ps(_mm_mul_ps(prod_ai_SSE0, icf4_SSE0),
                                                           _mm_mul_ps(prod_ai_SSE1, icf4_SSE1),
                                                           _mm_mul_ps(prod_ai_SSE2, icf4_SSE2),
                                                           _mm_mul_ps(prod_ai_SSE3, icf4_SSE3))));

            gpi_SSE0           = _mm_add_ps(gpi_SSE0, _mm_mul_ps(prod_SSE, icf4_SSE0));
            gpi_SSE1           = _mm_add_ps(gpi_SSE1, _mm_mul_ps(prod_SSE, icf4_SSE1));
            gpi_SSE2           = _mm_add_ps(gpi_SSE2, _mm_mul_ps(prod_SSE, icf4_SSE2));
            gpi_SSE3           = _mm_add_ps(gpi_SSE3, _mm_mul_ps(prod_SSE, icf4_SSE3));

            /* Save ai->aj and aj->ai chain rule terms */
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_SSE, icf6_SSE3));
            dadx += 4;

            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE0, icf6_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE1, icf6_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE2, icf6_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_mul_ps(prod_ai_SSE3, icf6_SSE3));
            dadx += 4;
        }
        _MM_TRANSPOSE4_PS(gpi_SSE0, gpi_SSE1, gpi_SSE2, gpi_SSE3);
        gpi_SSE0 = _mm_add_ps(gpi_SSE0, gpi_SSE1);
        gpi_SSE2 = _mm_add_ps(gpi_SSE2, gpi_SSE3);
        gpi_SSE0 = _mm_add_ps(gpi_SSE0, gpi_SSE2);
        _mm_store_ps(work+i, _mm_add_ps(gpi_SSE0, _mm_load_ps(work+i)));
    }

    /* In case we have written anything beyond natoms, move it back.
     * Never mind that we leave stuff above natoms; that will not
     * be accessed later in the routine.
     * In principle this should be a move rather than sum, but this
     * way we dont have to worry about even/odd offsets...
     */
    for (i = natoms; i < ni1+1+natoms/2; i++)
    {
        work[i-natoms] += work[i];
    }

    /* Parallel summations would go here if ever implemented with DD */

    factor  = 0.5 * ONE_4PI_EPS0;
    /* Calculate the radii - should we do all atoms, or just our local ones? */
    for (i = 0; i < natoms; i++)
    {
        if (born->use[i] != 0)
        {
            gpi             = born->gpol[i]+work[i];
            gpi2            = gpi * gpi;
            born->bRad[i]   = factor*gmx_invsqrt(gpi2);
            fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);
        }
    }

    return 0;
}



int
genborn_allvsall_calc_hct_obc_radii_sse2_single(t_forcerec *           fr,
                                                t_mdatoms *            mdatoms,
                                                gmx_genborn_t *        born,
                                                int                    gb_algorithm,
                                                gmx_localtop_t *       top,
                                                real *                 x,
                                                t_commrec *            cr,
                                                void *                 paadata)
{
    gmx_allvsallgb2_data_t *aadata;
    int                     natoms;
    int                     ni0, ni1;
    int                     nj0, nj1, nj2, nj3;
    int                     i, j, k, n;
    int              *      mask;
    int              *      pmask0;
    int              *      pmask1;
    int              *      pmask2;
    int              *      pmask3;
    int              *      emask0;
    int              *      emask1;
    int              *      emask2;
    int              *      emask3;
    real              *     gb_radius;
    real              *     vsolv;
    real              *     work;
    real                    tmpsum[4];
    real              *     x_align;
    real              *     y_align;
    real              *     z_align;
    int              *      jindex;
    real              *     dadx;
    real              *     obc_param;
    real                    rad, min_rad;
    real                    rai, rai_inv, rai_inv2, sum_ai, sum_ai2, sum_ai3, tsum, tchain;

    __m128                  ix_SSE0, iy_SSE0, iz_SSE0;
    __m128                  ix_SSE1, iy_SSE1, iz_SSE1;
    __m128                  ix_SSE2, iy_SSE2, iz_SSE2;
    __m128                  ix_SSE3, iy_SSE3, iz_SSE3;
    __m128                  gpi_SSE0, rai_SSE0, prod_ai_SSE0;
    __m128                  gpi_SSE1, rai_SSE1, prod_ai_SSE1;
    __m128                  gpi_SSE2, rai_SSE2, prod_ai_SSE2;
    __m128                  gpi_SSE3, rai_SSE3, prod_ai_SSE3;
    __m128                  imask_SSE0, jmask_SSE0;
    __m128                  imask_SSE1, jmask_SSE1;
    __m128                  imask_SSE2, jmask_SSE2;
    __m128                  imask_SSE3, jmask_SSE3;
    __m128                  jx_SSE, jy_SSE, jz_SSE;
    __m128                  dx_SSE0, dy_SSE0, dz_SSE0;
    __m128                  dx_SSE1, dy_SSE1, dz_SSE1;
    __m128                  dx_SSE2, dy_SSE2, dz_SSE2;
    __m128                  dx_SSE3, dy_SSE3, dz_SSE3;
    __m128                  rsq_SSE0, rinv_SSE0, irsq_SSE0, idr4_SSE0, idr6_SSE0;
    __m128                  rsq_SSE1, rinv_SSE1, irsq_SSE1, idr4_SSE1, idr6_SSE1;
    __m128                  rsq_SSE2, rinv_SSE2, irsq_SSE2, idr4_SSE2, idr6_SSE2;
    __m128                  rsq_SSE3, rinv_SSE3, irsq_SSE3, idr4_SSE3, idr6_SSE3;
    __m128                  raj_SSE, raj_inv_SSE, sk_aj_SSE, sk2_aj_SSE;
    __m128                  ccf_SSE0, dccf_SSE0, prod_SSE0;
    __m128                  ccf_SSE1, dccf_SSE1, prod_SSE1;
    __m128                  ccf_SSE2, dccf_SSE2, prod_SSE2;
    __m128                  ccf_SSE3, dccf_SSE3, prod_SSE3;
    __m128                  icf4_SSE0, icf6_SSE0;
    __m128                  icf4_SSE1, icf6_SSE1;
    __m128                  icf4_SSE2, icf6_SSE2;
    __m128                  icf4_SSE3, icf6_SSE3;
    __m128                  oneeighth_SSE, onefourth_SSE, half_SSE, one_SSE, two_SSE, four_SSE;
    __m128                  still_p4_SSE, still_p5inv_SSE, still_pip5_SSE;
    __m128                  rai_inv_SSE0;
    __m128                  rai_inv_SSE1;
    __m128                  rai_inv_SSE2;
    __m128                  rai_inv_SSE3;
    __m128                  sk_ai_SSE0, sk2_ai_SSE0, sum_ai_SSE0;
    __m128                  sk_ai_SSE1, sk2_ai_SSE1, sum_ai_SSE1;
    __m128                  sk_ai_SSE2, sk2_ai_SSE2, sum_ai_SSE2;
    __m128                  sk_ai_SSE3, sk2_ai_SSE3, sum_ai_SSE3;
    __m128                  lij_inv_SSE0, sk2_rinv_SSE0;
    __m128                  lij_inv_SSE1, sk2_rinv_SSE1;
    __m128                  lij_inv_SSE2, sk2_rinv_SSE2;
    __m128                  lij_inv_SSE3, sk2_rinv_SSE3;
    __m128                  dr_SSE0;
    __m128                  dr_SSE1;
    __m128                  dr_SSE2;
    __m128                  dr_SSE3;
    __m128                  t1_SSE0, t2_SSE0, t3_SSE0, t4_SSE0;
    __m128                  t1_SSE1, t2_SSE1, t3_SSE1, t4_SSE1;
    __m128                  t1_SSE2, t2_SSE2, t3_SSE2, t4_SSE2;
    __m128                  t1_SSE3, t2_SSE3, t3_SSE3, t4_SSE3;
    __m128                  obc_mask1_SSE0, obc_mask2_SSE0, obc_mask3_SSE0;
    __m128                  obc_mask1_SSE1, obc_mask2_SSE1, obc_mask3_SSE1;
    __m128                  obc_mask1_SSE2, obc_mask2_SSE2, obc_mask3_SSE2;
    __m128                  obc_mask1_SSE3, obc_mask2_SSE3, obc_mask3_SSE3;
    __m128                  uij_SSE0, uij2_SSE0, uij3_SSE0;
    __m128                  uij_SSE1, uij2_SSE1, uij3_SSE1;
    __m128                  uij_SSE2, uij2_SSE2, uij3_SSE2;
    __m128                  uij_SSE3, uij2_SSE3, uij3_SSE3;
    __m128                  lij_SSE0, lij2_SSE0, lij3_SSE0;
    __m128                  lij_SSE1, lij2_SSE1, lij3_SSE1;
    __m128                  lij_SSE2, lij2_SSE2, lij3_SSE2;
    __m128                  lij_SSE3, lij2_SSE3, lij3_SSE3;
    __m128                  dlij_SSE0, diff2_SSE0, logterm_SSE0;
    __m128                  dlij_SSE1, diff2_SSE1, logterm_SSE1;
    __m128                  dlij_SSE2, diff2_SSE2, logterm_SSE2;
    __m128                  dlij_SSE3, diff2_SSE3, logterm_SSE3;
    __m128                  doffset_SSE;

    natoms              = mdatoms->nr;
    ni0                 = 0;
    ni1                 = mdatoms->homenr;

    n = 0;

    aadata = *((gmx_allvsallgb2_data_t **)paadata);


    if (aadata == NULL)
    {
        genborn_allvsall_setup(&aadata, top, born, mdatoms, born->gb_doffset,
                               egbOBC, TRUE, TRUE, TRUE);
        *((gmx_allvsallgb2_data_t **)paadata) = aadata;
    }

    x_align = aadata->x_align;
    y_align = aadata->y_align;
    z_align = aadata->z_align;

    gb_radius = aadata->gb_radius;
    work      = aadata->work;
    jindex    = aadata->jindex_gb;
    dadx      = fr->dadx;
    obc_param = aadata->workparam;

    oneeighth_SSE   = _mm_set1_ps(0.125);
    onefourth_SSE   = _mm_set1_ps(0.25);
    half_SSE        = _mm_set1_ps(0.5);
    one_SSE         = _mm_set1_ps(1.0);
    two_SSE         = _mm_set1_ps(2.0);
    four_SSE        = _mm_set1_ps(4.0);
    doffset_SSE     = _mm_set1_ps(born->gb_doffset);

    for (i = 0; i < natoms; i++)
    {
        x_align[i]  = x[3*i];
        y_align[i]  = x[3*i+1];
        z_align[i]  = x[3*i+2];
    }

    /* Copy again */
    for (i = 0; i < natoms/2+1; i++)
    {
        x_align[natoms+i]  = x_align[i];
        y_align[natoms+i]  = y_align[i];
        z_align[natoms+i]  = z_align[i];
    }

    for (i = 0; i < natoms+natoms/2+1; i++)
    {
        work[i] = 0;
    }

    for (i = ni0; i < ni1; i += UNROLLI)
    {
        /* We assume shifts are NOT used for all-vs-all interactions */

        /* Load i atom data */
        ix_SSE0          = _mm_load1_ps(x_align+i);
        iy_SSE0          = _mm_load1_ps(y_align+i);
        iz_SSE0          = _mm_load1_ps(z_align+i);
        ix_SSE1          = _mm_load1_ps(x_align+i+1);
        iy_SSE1          = _mm_load1_ps(y_align+i+1);
        iz_SSE1          = _mm_load1_ps(z_align+i+1);
        ix_SSE2          = _mm_load1_ps(x_align+i+2);
        iy_SSE2          = _mm_load1_ps(y_align+i+2);
        iz_SSE2          = _mm_load1_ps(z_align+i+2);
        ix_SSE3          = _mm_load1_ps(x_align+i+3);
        iy_SSE3          = _mm_load1_ps(y_align+i+3);
        iz_SSE3          = _mm_load1_ps(z_align+i+3);

        rai_SSE0         = _mm_load1_ps(gb_radius+i);
        rai_SSE1         = _mm_load1_ps(gb_radius+i+1);
        rai_SSE2         = _mm_load1_ps(gb_radius+i+2);
        rai_SSE3         = _mm_load1_ps(gb_radius+i+3);
        rai_inv_SSE0     = gmx_mm_inv_ps(rai_SSE0);
        rai_inv_SSE1     = gmx_mm_inv_ps(rai_SSE1);
        rai_inv_SSE2     = gmx_mm_inv_ps(rai_SSE2);
        rai_inv_SSE3     = gmx_mm_inv_ps(rai_SSE3);

        sk_ai_SSE0       = _mm_load1_ps(obc_param+i);
        sk_ai_SSE1       = _mm_load1_ps(obc_param+i+1);
        sk_ai_SSE2       = _mm_load1_ps(obc_param+i+2);
        sk_ai_SSE3       = _mm_load1_ps(obc_param+i+3);
        sk2_ai_SSE0      = _mm_mul_ps(sk_ai_SSE0, sk_ai_SSE0);
        sk2_ai_SSE1      = _mm_mul_ps(sk_ai_SSE1, sk_ai_SSE1);
        sk2_ai_SSE2      = _mm_mul_ps(sk_ai_SSE2, sk_ai_SSE2);
        sk2_ai_SSE3      = _mm_mul_ps(sk_ai_SSE3, sk_ai_SSE3);

        sum_ai_SSE0      = _mm_setzero_ps();
        sum_ai_SSE1      = _mm_setzero_ps();
        sum_ai_SSE2      = _mm_setzero_ps();
        sum_ai_SSE3      = _mm_setzero_ps();

        /* Load limits for loop over neighbors */
        nj0              = jindex[4*i];
        nj1              = jindex[4*i+1];
        nj2              = jindex[4*i+2];
        nj3              = jindex[4*i+3];

        pmask0           = aadata->prologue_mask_gb[i];
        pmask1           = aadata->prologue_mask_gb[i+1];
        pmask2           = aadata->prologue_mask_gb[i+2];
        pmask3           = aadata->prologue_mask_gb[i+3];
        emask0           = aadata->epilogue_mask[i];
        emask1           = aadata->epilogue_mask[i+1];
        emask2           = aadata->epilogue_mask[i+2];
        emask3           = aadata->epilogue_mask[i+3];

        imask_SSE0        = _mm_load1_ps((real *)(aadata->imask+i));
        imask_SSE1        = _mm_load1_ps((real *)(aadata->imask+i+1));
        imask_SSE2        = _mm_load1_ps((real *)(aadata->imask+i+2));
        imask_SSE3        = _mm_load1_ps((real *)(aadata->imask+i+3));

        /* Prologue part, including exclusion mask */
        for (j = nj0; j < nj1; j += UNROLLJ)
        {
            jmask_SSE0 = _mm_load_ps((real *)pmask0);
            jmask_SSE1 = _mm_load_ps((real *)pmask1);
            jmask_SSE2 = _mm_load_ps((real *)pmask2);
            jmask_SSE3 = _mm_load_ps((real *)pmask3);
            pmask0    += UNROLLJ;
            pmask1    += UNROLLJ;
            pmask2    += UNROLLJ;
            pmask3    += UNROLLJ;

            /* load j atom coordinates */
            jx_SSE            = _mm_load_ps(x_align+j);
            jy_SSE            = _mm_load_ps(y_align+j);
            jz_SSE            = _mm_load_ps(z_align+j);

            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0, jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0, jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0, jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1, jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1, jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1, jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2, jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2, jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2, jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3, jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3, jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1, dy_SSE1, dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2, dy_SSE2, dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3, dy_SSE3, dz_SSE3);

            /* Combine masks */
            jmask_SSE0         = _mm_and_ps(jmask_SSE0, imask_SSE0);
            jmask_SSE1         = _mm_and_ps(jmask_SSE1, imask_SSE1);
            jmask_SSE2         = _mm_and_ps(jmask_SSE2, imask_SSE2);
            jmask_SSE3         = _mm_and_ps(jmask_SSE3, imask_SSE3);

            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);

            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0, jmask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1, jmask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2, jmask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3, jmask_SSE3);

            dr_SSE0            = _mm_mul_ps(rsq_SSE0, rinv_SSE0);
            dr_SSE1            = _mm_mul_ps(rsq_SSE1, rinv_SSE1);
            dr_SSE2            = _mm_mul_ps(rsq_SSE2, rinv_SSE2);
            dr_SSE3            = _mm_mul_ps(rsq_SSE3, rinv_SSE3);

            sk_aj_SSE          = _mm_load_ps(obc_param+j);
            raj_SSE            = _mm_load_ps(gb_radius+j);
            raj_inv_SSE        = gmx_mm_inv_ps(raj_SSE);

            /* Evaluate influence of atom aj -> ai */
            t1_SSE0            = _mm_add_ps(dr_SSE0, sk_aj_SSE);
            t1_SSE1            = _mm_add_ps(dr_SSE1, sk_aj_SSE);
            t1_SSE2            = _mm_add_ps(dr_SSE2, sk_aj_SSE);
            t1_SSE3            = _mm_add_ps(dr_SSE3, sk_aj_SSE);
            t2_SSE0            = _mm_sub_ps(dr_SSE0, sk_aj_SSE);
            t2_SSE1            = _mm_sub_ps(dr_SSE1, sk_aj_SSE);
            t2_SSE2            = _mm_sub_ps(dr_SSE2, sk_aj_SSE);
            t2_SSE3            = _mm_sub_ps(dr_SSE3, sk_aj_SSE);
            t3_SSE0            = _mm_sub_ps(sk_aj_SSE, dr_SSE0);
            t3_SSE1            = _mm_sub_ps(sk_aj_SSE, dr_SSE1);
            t3_SSE2            = _mm_sub_ps(sk_aj_SSE, dr_SSE2);
            t3_SSE3            = _mm_sub_ps(sk_aj_SSE, dr_SSE3);

            obc_mask1_SSE0     = _mm_cmplt_ps(rai_SSE0, t1_SSE0);
            obc_mask1_SSE1     = _mm_cmplt_ps(rai_SSE1, t1_SSE1);
            obc_mask1_SSE2     = _mm_cmplt_ps(rai_SSE2, t1_SSE2);
            obc_mask1_SSE3     = _mm_cmplt_ps(rai_SSE3, t1_SSE3);
            obc_mask2_SSE0     = _mm_cmplt_ps(rai_SSE0, t2_SSE0);
            obc_mask2_SSE1     = _mm_cmplt_ps(rai_SSE1, t2_SSE1);
            obc_mask2_SSE2     = _mm_cmplt_ps(rai_SSE2, t2_SSE2);
            obc_mask2_SSE3     = _mm_cmplt_ps(rai_SSE3, t2_SSE3);
            obc_mask3_SSE0     = _mm_cmplt_ps(rai_SSE0, t3_SSE0);
            obc_mask3_SSE1     = _mm_cmplt_ps(rai_SSE1, t3_SSE1);
            obc_mask3_SSE2     = _mm_cmplt_ps(rai_SSE2, t3_SSE2);
            obc_mask3_SSE3     = _mm_cmplt_ps(rai_SSE3, t3_SSE3);
            obc_mask1_SSE0     = _mm_and_ps(obc_mask1_SSE0, jmask_SSE0);
            obc_mask1_SSE1     = _mm_and_ps(obc_mask1_SSE1, jmask_SSE1);
            obc_mask1_SSE2     = _mm_and_ps(obc_mask1_SSE2, jmask_SSE2);
            obc_mask1_SSE3     = _mm_and_ps(obc_mask1_SSE3, jmask_SSE3);

            uij_SSE0           = gmx_mm_inv_ps(t1_SSE0);
            uij_SSE1           = gmx_mm_inv_ps(t1_SSE1);
            uij_SSE2           = gmx_mm_inv_ps(t1_SSE2);
            uij_SSE3           = gmx_mm_inv_ps(t1_SSE3);
            lij_SSE0           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE0, gmx_mm_inv_ps(t2_SSE0)),
                                              _mm_andnot_ps(obc_mask2_SSE0, rai_inv_SSE0));
            lij_SSE1           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE1, gmx_mm_inv_ps(t2_SSE1)),
                                              _mm_andnot_ps(obc_mask2_SSE1, rai_inv_SSE1));
            lij_SSE2           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE2, gmx_mm_inv_ps(t2_SSE2)),
                                              _mm_andnot_ps(obc_mask2_SSE2, rai_inv_SSE2));
            lij_SSE3           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE3, gmx_mm_inv_ps(t2_SSE3)),
                                              _mm_andnot_ps(obc_mask2_SSE3, rai_inv_SSE3));
            dlij_SSE0          = _mm_and_ps(one_SSE, obc_mask2_SSE0);
            dlij_SSE1          = _mm_and_ps(one_SSE, obc_mask2_SSE1);
            dlij_SSE2          = _mm_and_ps(one_SSE, obc_mask2_SSE2);
            dlij_SSE3          = _mm_and_ps(one_SSE, obc_mask2_SSE3);

            uij2_SSE0          = _mm_mul_ps(uij_SSE0, uij_SSE0);
            uij2_SSE1          = _mm_mul_ps(uij_SSE1, uij_SSE1);
            uij2_SSE2          = _mm_mul_ps(uij_SSE2, uij_SSE2);
            uij2_SSE3          = _mm_mul_ps(uij_SSE3, uij_SSE3);
            uij3_SSE0          = _mm_mul_ps(uij2_SSE0, uij_SSE0);
            uij3_SSE1          = _mm_mul_ps(uij2_SSE1, uij_SSE1);
            uij3_SSE2          = _mm_mul_ps(uij2_SSE2, uij_SSE2);
            uij3_SSE3          = _mm_mul_ps(uij2_SSE3, uij_SSE3);
            lij2_SSE0          = _mm_mul_ps(lij_SSE0, lij_SSE0);
            lij2_SSE1          = _mm_mul_ps(lij_SSE1, lij_SSE1);
            lij2_SSE2          = _mm_mul_ps(lij_SSE2, lij_SSE2);
            lij2_SSE3          = _mm_mul_ps(lij_SSE3, lij_SSE3);
            lij3_SSE0          = _mm_mul_ps(lij2_SSE0, lij_SSE0);
            lij3_SSE1          = _mm_mul_ps(lij2_SSE1, lij_SSE1);
            lij3_SSE2          = _mm_mul_ps(lij2_SSE2, lij_SSE2);
            lij3_SSE3          = _mm_mul_ps(lij2_SSE3, lij_SSE3);

            diff2_SSE0         = _mm_sub_ps(uij2_SSE0, lij2_SSE0);
            diff2_SSE1         = _mm_sub_ps(uij2_SSE1, lij2_SSE1);
            diff2_SSE2         = _mm_sub_ps(uij2_SSE2, lij2_SSE2);
            diff2_SSE3         = _mm_sub_ps(uij2_SSE3, lij2_SSE3);
            lij_inv_SSE0       = gmx_mm_invsqrt_ps(lij2_SSE0);
            lij_inv_SSE1       = gmx_mm_invsqrt_ps(lij2_SSE1);
            lij_inv_SSE2       = gmx_mm_invsqrt_ps(lij2_SSE2);
            lij_inv_SSE3       = gmx_mm_invsqrt_ps(lij2_SSE3);
            sk2_aj_SSE         = _mm_mul_ps(sk_aj_SSE, sk_aj_SSE);
            sk2_rinv_SSE0      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE0);
            sk2_rinv_SSE1      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE1);
            sk2_rinv_SSE2      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE2);
            sk2_rinv_SSE3      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE3);
            prod_SSE0          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE0);
            prod_SSE1          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE1);
            prod_SSE2          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE2);
            prod_SSE3          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE3);

            logterm_SSE0       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE0, lij_inv_SSE0));
            logterm_SSE1       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE1, lij_inv_SSE1));
            logterm_SSE2       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE2, lij_inv_SSE2));
            logterm_SSE3       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE3, lij_inv_SSE3));

            t1_SSE0            = _mm_sub_ps(lij_SSE0, uij_SSE0);
            t1_SSE1            = _mm_sub_ps(lij_SSE1, uij_SSE1);
            t1_SSE2            = _mm_sub_ps(lij_SSE2, uij_SSE2);
            t1_SSE3            = _mm_sub_ps(lij_SSE3, uij_SSE3);
            t2_SSE0            = _mm_mul_ps(diff2_SSE0,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE0),
                                                       prod_SSE0));
            t2_SSE1            = _mm_mul_ps(diff2_SSE1,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE1),
                                                       prod_SSE1));
            t2_SSE2            = _mm_mul_ps(diff2_SSE2,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE2),
                                                       prod_SSE2));
            t2_SSE3            = _mm_mul_ps(diff2_SSE3,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE3),
                                                       prod_SSE3));

            t3_SSE0            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE0, logterm_SSE0));
            t3_SSE1            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE1, logterm_SSE1));
            t3_SSE2            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE2, logterm_SSE2));
            t3_SSE3            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE3, logterm_SSE3));
            t1_SSE0            = _mm_add_ps(t1_SSE0, _mm_add_ps(t2_SSE0, t3_SSE0));
            t1_SSE1            = _mm_add_ps(t1_SSE1, _mm_add_ps(t2_SSE1, t3_SSE1));
            t1_SSE2            = _mm_add_ps(t1_SSE2, _mm_add_ps(t2_SSE2, t3_SSE2));
            t1_SSE3            = _mm_add_ps(t1_SSE3, _mm_add_ps(t2_SSE3, t3_SSE3));
            t4_SSE0            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE0, lij_SSE0));
            t4_SSE1            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE1, lij_SSE1));
            t4_SSE2            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE2, lij_SSE2));
            t4_SSE3            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE3, lij_SSE3));
            t4_SSE0            = _mm_and_ps(t4_SSE0, obc_mask3_SSE0);
            t4_SSE1            = _mm_and_ps(t4_SSE1, obc_mask3_SSE1);
            t4_SSE2            = _mm_and_ps(t4_SSE2, obc_mask3_SSE2);
            t4_SSE3            = _mm_and_ps(t4_SSE3, obc_mask3_SSE3);
            t1_SSE0            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE0, t4_SSE0));
            t1_SSE1            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE1, t4_SSE1));
            t1_SSE2            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE2, t4_SSE2));
            t1_SSE3            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE3, t4_SSE3));

            sum_ai_SSE0        = _mm_add_ps(sum_ai_SSE0, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            sum_ai_SSE1        = _mm_add_ps(sum_ai_SSE1, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            sum_ai_SSE2        = _mm_add_ps(sum_ai_SSE2, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            sum_ai_SSE3        = _mm_add_ps(sum_ai_SSE3, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));

            t1_SSE0            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE0),
                                            _mm_mul_ps(prod_SSE0, lij3_SSE0));
            t1_SSE1            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE1),
                                            _mm_mul_ps(prod_SSE1, lij3_SSE1));
            t1_SSE2            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE2),
                                            _mm_mul_ps(prod_SSE2, lij3_SSE2));
            t1_SSE3            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE3),
                                            _mm_mul_ps(prod_SSE3, lij3_SSE3));
            t1_SSE0            = _mm_sub_ps(t1_SSE0,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE0, rinv_SSE0),
                                                                  _mm_mul_ps(lij3_SSE0, dr_SSE0))));
            t1_SSE1            = _mm_sub_ps(t1_SSE1,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE1, rinv_SSE1),
                                                                  _mm_mul_ps(lij3_SSE1, dr_SSE1))));
            t1_SSE2            = _mm_sub_ps(t1_SSE2,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE2, rinv_SSE2),
                                                                  _mm_mul_ps(lij3_SSE2, dr_SSE2))));
            t1_SSE3            = _mm_sub_ps(t1_SSE3,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE3, rinv_SSE3),
                                                                  _mm_mul_ps(lij3_SSE3, dr_SSE3))));

            t2_SSE0            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE0, rinv_SSE0),
                                                       _mm_mul_ps(uij3_SSE0, dr_SSE0)));
            t2_SSE1            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE1, rinv_SSE1),
                                                       _mm_mul_ps(uij3_SSE1, dr_SSE1)));
            t2_SSE2            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE2, rinv_SSE2),
                                                       _mm_mul_ps(uij3_SSE2, dr_SSE2)));
            t2_SSE3            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE3, rinv_SSE3),
                                                       _mm_mul_ps(uij3_SSE3, dr_SSE3)));
            t2_SSE0            = _mm_sub_ps(t2_SSE0,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE0),
                                                       _mm_mul_ps(prod_SSE0, uij3_SSE0)));
            t2_SSE1            = _mm_sub_ps(t2_SSE1,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE1),
                                                       _mm_mul_ps(prod_SSE1, uij3_SSE1)));
            t2_SSE2            = _mm_sub_ps(t2_SSE2,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE2),
                                                       _mm_mul_ps(prod_SSE2, uij3_SSE2)));
            t2_SSE3            = _mm_sub_ps(t2_SSE3,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE3),
                                                       _mm_mul_ps(prod_SSE3, uij3_SSE3)));
            t3_SSE0            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE0),
                                            _mm_mul_ps(rinv_SSE0, rinv_SSE0));
            t3_SSE1            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE1),
                                            _mm_mul_ps(rinv_SSE1, rinv_SSE1));
            t3_SSE2            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE2),
                                            _mm_mul_ps(rinv_SSE2, rinv_SSE2));
            t3_SSE3            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE3),
                                            _mm_mul_ps(rinv_SSE3, rinv_SSE3));
            t3_SSE0            = _mm_sub_ps(t3_SSE0,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE0, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE0, rinv_SSE0))));
            t3_SSE1            = _mm_sub_ps(t3_SSE1,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE1, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE1, rinv_SSE1))));
            t3_SSE2            = _mm_sub_ps(t3_SSE2,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE2, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE2, rinv_SSE2))));
            t3_SSE3            = _mm_sub_ps(t3_SSE3,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE3, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE3, rinv_SSE3))));

            t1_SSE0            = _mm_mul_ps(rinv_SSE0,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE0, t1_SSE0),
                                                       _mm_add_ps(t2_SSE0, t3_SSE0)));
            t1_SSE1            = _mm_mul_ps(rinv_SSE1,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE1, t1_SSE1),
                                                       _mm_add_ps(t2_SSE1, t3_SSE1)));
            t1_SSE2            = _mm_mul_ps(rinv_SSE2,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE2, t1_SSE2),
                                                       _mm_add_ps(t2_SSE2, t3_SSE2)));
            t1_SSE3            = _mm_mul_ps(rinv_SSE3,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE3, t1_SSE3),
                                                       _mm_add_ps(t2_SSE3, t3_SSE3)));

            _mm_store_ps(dadx, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));
            dadx += 4;

            /* Evaluate influence of atom ai -> aj */
            t1_SSE0            = _mm_add_ps(dr_SSE0, sk_ai_SSE0);
            t1_SSE1            = _mm_add_ps(dr_SSE1, sk_ai_SSE1);
            t1_SSE2            = _mm_add_ps(dr_SSE2, sk_ai_SSE2);
            t1_SSE3            = _mm_add_ps(dr_SSE3, sk_ai_SSE3);
            t2_SSE0            = _mm_sub_ps(dr_SSE0, sk_ai_SSE0);
            t2_SSE1            = _mm_sub_ps(dr_SSE1, sk_ai_SSE1);
            t2_SSE2            = _mm_sub_ps(dr_SSE2, sk_ai_SSE2);
            t2_SSE3            = _mm_sub_ps(dr_SSE3, sk_ai_SSE3);
            t3_SSE0            = _mm_sub_ps(sk_ai_SSE0, dr_SSE0);
            t3_SSE1            = _mm_sub_ps(sk_ai_SSE1, dr_SSE1);
            t3_SSE2            = _mm_sub_ps(sk_ai_SSE2, dr_SSE2);
            t3_SSE3            = _mm_sub_ps(sk_ai_SSE3, dr_SSE3);

            obc_mask1_SSE0     = _mm_cmplt_ps(raj_SSE, t1_SSE0);
            obc_mask1_SSE1     = _mm_cmplt_ps(raj_SSE, t1_SSE1);
            obc_mask1_SSE2     = _mm_cmplt_ps(raj_SSE, t1_SSE2);
            obc_mask1_SSE3     = _mm_cmplt_ps(raj_SSE, t1_SSE3);
            obc_mask2_SSE0     = _mm_cmplt_ps(raj_SSE, t2_SSE0);
            obc_mask2_SSE1     = _mm_cmplt_ps(raj_SSE, t2_SSE1);
            obc_mask2_SSE2     = _mm_cmplt_ps(raj_SSE, t2_SSE2);
            obc_mask2_SSE3     = _mm_cmplt_ps(raj_SSE, t2_SSE3);
            obc_mask3_SSE0     = _mm_cmplt_ps(raj_SSE, t3_SSE0);
            obc_mask3_SSE1     = _mm_cmplt_ps(raj_SSE, t3_SSE1);
            obc_mask3_SSE2     = _mm_cmplt_ps(raj_SSE, t3_SSE2);
            obc_mask3_SSE3     = _mm_cmplt_ps(raj_SSE, t3_SSE3);
            obc_mask1_SSE0     = _mm_and_ps(obc_mask1_SSE0, jmask_SSE0);
            obc_mask1_SSE1     = _mm_and_ps(obc_mask1_SSE1, jmask_SSE1);
            obc_mask1_SSE2     = _mm_and_ps(obc_mask1_SSE2, jmask_SSE2);
            obc_mask1_SSE3     = _mm_and_ps(obc_mask1_SSE3, jmask_SSE3);

            uij_SSE0           = gmx_mm_inv_ps(t1_SSE0);
            uij_SSE1           = gmx_mm_inv_ps(t1_SSE1);
            uij_SSE2           = gmx_mm_inv_ps(t1_SSE2);
            uij_SSE3           = gmx_mm_inv_ps(t1_SSE3);
            lij_SSE0           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE0, gmx_mm_inv_ps(t2_SSE0)),
                                              _mm_andnot_ps(obc_mask2_SSE0, raj_inv_SSE));
            lij_SSE1           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE1, gmx_mm_inv_ps(t2_SSE1)),
                                              _mm_andnot_ps(obc_mask2_SSE1, raj_inv_SSE));
            lij_SSE2           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE2, gmx_mm_inv_ps(t2_SSE2)),
                                              _mm_andnot_ps(obc_mask2_SSE2, raj_inv_SSE));
            lij_SSE3           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE3, gmx_mm_inv_ps(t2_SSE3)),
                                              _mm_andnot_ps(obc_mask2_SSE3, raj_inv_SSE));
            dlij_SSE0          = _mm_and_ps(one_SSE, obc_mask2_SSE0);
            dlij_SSE1          = _mm_and_ps(one_SSE, obc_mask2_SSE1);
            dlij_SSE2          = _mm_and_ps(one_SSE, obc_mask2_SSE2);
            dlij_SSE3          = _mm_and_ps(one_SSE, obc_mask2_SSE3);

            uij2_SSE0          = _mm_mul_ps(uij_SSE0, uij_SSE0);
            uij2_SSE1          = _mm_mul_ps(uij_SSE1, uij_SSE1);
            uij2_SSE2          = _mm_mul_ps(uij_SSE2, uij_SSE2);
            uij2_SSE3          = _mm_mul_ps(uij_SSE3, uij_SSE3);
            uij3_SSE0          = _mm_mul_ps(uij2_SSE0, uij_SSE0);
            uij3_SSE1          = _mm_mul_ps(uij2_SSE1, uij_SSE1);
            uij3_SSE2          = _mm_mul_ps(uij2_SSE2, uij_SSE2);
            uij3_SSE3          = _mm_mul_ps(uij2_SSE3, uij_SSE3);
            lij2_SSE0          = _mm_mul_ps(lij_SSE0, lij_SSE0);
            lij2_SSE1          = _mm_mul_ps(lij_SSE1, lij_SSE1);
            lij2_SSE2          = _mm_mul_ps(lij_SSE2, lij_SSE2);
            lij2_SSE3          = _mm_mul_ps(lij_SSE3, lij_SSE3);
            lij3_SSE0          = _mm_mul_ps(lij2_SSE0, lij_SSE0);
            lij3_SSE1          = _mm_mul_ps(lij2_SSE1, lij_SSE1);
            lij3_SSE2          = _mm_mul_ps(lij2_SSE2, lij_SSE2);
            lij3_SSE3          = _mm_mul_ps(lij2_SSE3, lij_SSE3);

            diff2_SSE0         = _mm_sub_ps(uij2_SSE0, lij2_SSE0);
            diff2_SSE1         = _mm_sub_ps(uij2_SSE1, lij2_SSE1);
            diff2_SSE2         = _mm_sub_ps(uij2_SSE2, lij2_SSE2);
            diff2_SSE3         = _mm_sub_ps(uij2_SSE3, lij2_SSE3);
            lij_inv_SSE0       = gmx_mm_invsqrt_ps(lij2_SSE0);
            lij_inv_SSE1       = gmx_mm_invsqrt_ps(lij2_SSE1);
            lij_inv_SSE2       = gmx_mm_invsqrt_ps(lij2_SSE2);
            lij_inv_SSE3       = gmx_mm_invsqrt_ps(lij2_SSE3);
            sk2_rinv_SSE0      = _mm_mul_ps(sk2_ai_SSE0, rinv_SSE0);
            sk2_rinv_SSE1      = _mm_mul_ps(sk2_ai_SSE1, rinv_SSE1);
            sk2_rinv_SSE2      = _mm_mul_ps(sk2_ai_SSE2, rinv_SSE2);
            sk2_rinv_SSE3      = _mm_mul_ps(sk2_ai_SSE3, rinv_SSE3);
            prod_SSE0          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE0);
            prod_SSE1          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE1);
            prod_SSE2          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE2);
            prod_SSE3          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE3);

            logterm_SSE0       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE0, lij_inv_SSE0));
            logterm_SSE1       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE1, lij_inv_SSE1));
            logterm_SSE2       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE2, lij_inv_SSE2));
            logterm_SSE3       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE3, lij_inv_SSE3));
            t1_SSE0            = _mm_sub_ps(lij_SSE0, uij_SSE0);
            t1_SSE1            = _mm_sub_ps(lij_SSE1, uij_SSE1);
            t1_SSE2            = _mm_sub_ps(lij_SSE2, uij_SSE2);
            t1_SSE3            = _mm_sub_ps(lij_SSE3, uij_SSE3);
            t2_SSE0            = _mm_mul_ps(diff2_SSE0,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE0),
                                                       prod_SSE0));
            t2_SSE1            = _mm_mul_ps(diff2_SSE1,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE1),
                                                       prod_SSE1));
            t2_SSE2            = _mm_mul_ps(diff2_SSE2,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE2),
                                                       prod_SSE2));
            t2_SSE3            = _mm_mul_ps(diff2_SSE3,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE3),
                                                       prod_SSE3));
            t3_SSE0            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE0, logterm_SSE0));
            t3_SSE1            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE1, logterm_SSE1));
            t3_SSE2            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE2, logterm_SSE2));
            t3_SSE3            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE3, logterm_SSE3));
            t1_SSE0            = _mm_add_ps(t1_SSE0, _mm_add_ps(t2_SSE0, t3_SSE0));
            t1_SSE1            = _mm_add_ps(t1_SSE1, _mm_add_ps(t2_SSE1, t3_SSE1));
            t1_SSE2            = _mm_add_ps(t1_SSE2, _mm_add_ps(t2_SSE2, t3_SSE2));
            t1_SSE3            = _mm_add_ps(t1_SSE3, _mm_add_ps(t2_SSE3, t3_SSE3));
            t4_SSE0            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE0));
            t4_SSE1            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE1));
            t4_SSE2            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE2));
            t4_SSE3            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE3));
            t4_SSE0            = _mm_and_ps(t4_SSE0, obc_mask3_SSE0);
            t4_SSE1            = _mm_and_ps(t4_SSE1, obc_mask3_SSE1);
            t4_SSE2            = _mm_and_ps(t4_SSE2, obc_mask3_SSE2);
            t4_SSE3            = _mm_and_ps(t4_SSE3, obc_mask3_SSE3);
            t1_SSE0            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE0, t4_SSE0));
            t1_SSE1            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE1, t4_SSE1));
            t1_SSE2            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE2, t4_SSE2));
            t1_SSE3            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE3, t4_SSE3));

            _mm_store_ps(work+j, _mm_add_ps(_mm_load_ps(work+j),
                                            gmx_mm_sum4_ps(_mm_and_ps(t1_SSE0, obc_mask1_SSE0),
                                                           _mm_and_ps(t1_SSE1, obc_mask1_SSE1),
                                                           _mm_and_ps(t1_SSE2, obc_mask1_SSE2),
                                                           _mm_and_ps(t1_SSE3, obc_mask1_SSE3))));

            t1_SSE0            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE0),
                                            _mm_mul_ps(prod_SSE0, lij3_SSE0));
            t1_SSE1            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE1),
                                            _mm_mul_ps(prod_SSE1, lij3_SSE1));
            t1_SSE2            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE2),
                                            _mm_mul_ps(prod_SSE2, lij3_SSE2));
            t1_SSE3            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE3),
                                            _mm_mul_ps(prod_SSE3, lij3_SSE3));
            t1_SSE0            = _mm_sub_ps(t1_SSE0,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE0, rinv_SSE0),
                                                                  _mm_mul_ps(lij3_SSE0, dr_SSE0))));
            t1_SSE1            = _mm_sub_ps(t1_SSE1,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE1, rinv_SSE1),
                                                                  _mm_mul_ps(lij3_SSE1, dr_SSE1))));
            t1_SSE2            = _mm_sub_ps(t1_SSE2,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE2, rinv_SSE2),
                                                                  _mm_mul_ps(lij3_SSE2, dr_SSE2))));
            t1_SSE3            = _mm_sub_ps(t1_SSE3,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE3, rinv_SSE3),
                                                                  _mm_mul_ps(lij3_SSE3, dr_SSE3))));
            t2_SSE0            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE0, rinv_SSE0),
                                                       _mm_mul_ps(uij3_SSE0, dr_SSE0)));
            t2_SSE1            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE1, rinv_SSE1),
                                                       _mm_mul_ps(uij3_SSE1, dr_SSE1)));
            t2_SSE2            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE2, rinv_SSE2),
                                                       _mm_mul_ps(uij3_SSE2, dr_SSE2)));
            t2_SSE3            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE3, rinv_SSE3),
                                                       _mm_mul_ps(uij3_SSE3, dr_SSE3)));
            t2_SSE0            = _mm_sub_ps(t2_SSE0,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE0),
                                                       _mm_mul_ps(prod_SSE0, uij3_SSE0)));
            t2_SSE1            = _mm_sub_ps(t2_SSE1,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE1),
                                                       _mm_mul_ps(prod_SSE1, uij3_SSE1)));
            t2_SSE2            = _mm_sub_ps(t2_SSE2,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE2),
                                                       _mm_mul_ps(prod_SSE2, uij3_SSE2)));
            t2_SSE3            = _mm_sub_ps(t2_SSE3,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE3),
                                                       _mm_mul_ps(prod_SSE3, uij3_SSE3)));

            t3_SSE0            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE0),
                                            _mm_mul_ps(rinv_SSE0, rinv_SSE0));
            t3_SSE1            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE1),
                                            _mm_mul_ps(rinv_SSE1, rinv_SSE1));
            t3_SSE2            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE2),
                                            _mm_mul_ps(rinv_SSE2, rinv_SSE2));
            t3_SSE3            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE3),
                                            _mm_mul_ps(rinv_SSE3, rinv_SSE3));

            t3_SSE0            = _mm_sub_ps(t3_SSE0,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE0, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE0, rinv_SSE0))));
            t3_SSE1            = _mm_sub_ps(t3_SSE1,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE1, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE1, rinv_SSE1))));
            t3_SSE2            = _mm_sub_ps(t3_SSE2,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE2, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE2, rinv_SSE2))));
            t3_SSE3            = _mm_sub_ps(t3_SSE3,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE3, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE3, rinv_SSE3))));


            t1_SSE0            = _mm_mul_ps(rinv_SSE0,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE0, t1_SSE0),
                                                       _mm_add_ps(t2_SSE0, t3_SSE0)));
            t1_SSE1            = _mm_mul_ps(rinv_SSE1,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE1, t1_SSE1),
                                                       _mm_add_ps(t2_SSE1, t3_SSE1)));
            t1_SSE2            = _mm_mul_ps(rinv_SSE2,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE2, t1_SSE2),
                                                       _mm_add_ps(t2_SSE2, t3_SSE2)));
            t1_SSE3            = _mm_mul_ps(rinv_SSE3,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE3, t1_SSE3),
                                                       _mm_add_ps(t2_SSE3, t3_SSE3)));

            _mm_store_ps(dadx, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));
            dadx += 4;

        }

        /* Main part, no exclusions */
        for (j = nj1; j < nj2; j += UNROLLJ)
        {
            /* load j atom coordinates */
            jx_SSE            = _mm_load_ps(x_align+j);
            jy_SSE            = _mm_load_ps(y_align+j);
            jz_SSE            = _mm_load_ps(z_align+j);

            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0, jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0, jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0, jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1, jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1, jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1, jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2, jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2, jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2, jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3, jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3, jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1, dy_SSE1, dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2, dy_SSE2, dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3, dy_SSE3, dz_SSE3);

            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);

            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0, imask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1, imask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2, imask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3, imask_SSE3);

            dr_SSE0            = _mm_mul_ps(rsq_SSE0, rinv_SSE0);
            dr_SSE1            = _mm_mul_ps(rsq_SSE1, rinv_SSE1);
            dr_SSE2            = _mm_mul_ps(rsq_SSE2, rinv_SSE2);
            dr_SSE3            = _mm_mul_ps(rsq_SSE3, rinv_SSE3);

            sk_aj_SSE          = _mm_load_ps(obc_param+j);
            raj_SSE            = _mm_load_ps(gb_radius+j);

            raj_inv_SSE        = gmx_mm_inv_ps(raj_SSE);

            /* Evaluate influence of atom aj -> ai */
            t1_SSE0            = _mm_add_ps(dr_SSE0, sk_aj_SSE);
            t1_SSE1            = _mm_add_ps(dr_SSE1, sk_aj_SSE);
            t1_SSE2            = _mm_add_ps(dr_SSE2, sk_aj_SSE);
            t1_SSE3            = _mm_add_ps(dr_SSE3, sk_aj_SSE);
            t2_SSE0            = _mm_sub_ps(dr_SSE0, sk_aj_SSE);
            t2_SSE1            = _mm_sub_ps(dr_SSE1, sk_aj_SSE);
            t2_SSE2            = _mm_sub_ps(dr_SSE2, sk_aj_SSE);
            t2_SSE3            = _mm_sub_ps(dr_SSE3, sk_aj_SSE);
            t3_SSE0            = _mm_sub_ps(sk_aj_SSE, dr_SSE0);
            t3_SSE1            = _mm_sub_ps(sk_aj_SSE, dr_SSE1);
            t3_SSE2            = _mm_sub_ps(sk_aj_SSE, dr_SSE2);
            t3_SSE3            = _mm_sub_ps(sk_aj_SSE, dr_SSE3);

            obc_mask1_SSE0     = _mm_cmplt_ps(rai_SSE0, t1_SSE0);
            obc_mask1_SSE1     = _mm_cmplt_ps(rai_SSE1, t1_SSE1);
            obc_mask1_SSE2     = _mm_cmplt_ps(rai_SSE2, t1_SSE2);
            obc_mask1_SSE3     = _mm_cmplt_ps(rai_SSE3, t1_SSE3);
            obc_mask2_SSE0     = _mm_cmplt_ps(rai_SSE0, t2_SSE0);
            obc_mask2_SSE1     = _mm_cmplt_ps(rai_SSE1, t2_SSE1);
            obc_mask2_SSE2     = _mm_cmplt_ps(rai_SSE2, t2_SSE2);
            obc_mask2_SSE3     = _mm_cmplt_ps(rai_SSE3, t2_SSE3);
            obc_mask3_SSE0     = _mm_cmplt_ps(rai_SSE0, t3_SSE0);
            obc_mask3_SSE1     = _mm_cmplt_ps(rai_SSE1, t3_SSE1);
            obc_mask3_SSE2     = _mm_cmplt_ps(rai_SSE2, t3_SSE2);
            obc_mask3_SSE3     = _mm_cmplt_ps(rai_SSE3, t3_SSE3);
            obc_mask1_SSE0     = _mm_and_ps(obc_mask1_SSE0, imask_SSE0);
            obc_mask1_SSE1     = _mm_and_ps(obc_mask1_SSE1, imask_SSE1);
            obc_mask1_SSE2     = _mm_and_ps(obc_mask1_SSE2, imask_SSE2);
            obc_mask1_SSE3     = _mm_and_ps(obc_mask1_SSE3, imask_SSE3);

            uij_SSE0           = gmx_mm_inv_ps(t1_SSE0);
            uij_SSE1           = gmx_mm_inv_ps(t1_SSE1);
            uij_SSE2           = gmx_mm_inv_ps(t1_SSE2);
            uij_SSE3           = gmx_mm_inv_ps(t1_SSE3);
            lij_SSE0           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE0, gmx_mm_inv_ps(t2_SSE0)),
                                              _mm_andnot_ps(obc_mask2_SSE0, rai_inv_SSE0));
            lij_SSE1           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE1, gmx_mm_inv_ps(t2_SSE1)),
                                              _mm_andnot_ps(obc_mask2_SSE1, rai_inv_SSE1));
            lij_SSE2           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE2, gmx_mm_inv_ps(t2_SSE2)),
                                              _mm_andnot_ps(obc_mask2_SSE2, rai_inv_SSE2));
            lij_SSE3           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE3, gmx_mm_inv_ps(t2_SSE3)),
                                              _mm_andnot_ps(obc_mask2_SSE3, rai_inv_SSE3));
            dlij_SSE0          = _mm_and_ps(one_SSE, obc_mask2_SSE0);
            dlij_SSE1          = _mm_and_ps(one_SSE, obc_mask2_SSE1);
            dlij_SSE2          = _mm_and_ps(one_SSE, obc_mask2_SSE2);
            dlij_SSE3          = _mm_and_ps(one_SSE, obc_mask2_SSE3);

            uij2_SSE0          = _mm_mul_ps(uij_SSE0, uij_SSE0);
            uij2_SSE1          = _mm_mul_ps(uij_SSE1, uij_SSE1);
            uij2_SSE2          = _mm_mul_ps(uij_SSE2, uij_SSE2);
            uij2_SSE3          = _mm_mul_ps(uij_SSE3, uij_SSE3);
            uij3_SSE0          = _mm_mul_ps(uij2_SSE0, uij_SSE0);
            uij3_SSE1          = _mm_mul_ps(uij2_SSE1, uij_SSE1);
            uij3_SSE2          = _mm_mul_ps(uij2_SSE2, uij_SSE2);
            uij3_SSE3          = _mm_mul_ps(uij2_SSE3, uij_SSE3);
            lij2_SSE0          = _mm_mul_ps(lij_SSE0, lij_SSE0);
            lij2_SSE1          = _mm_mul_ps(lij_SSE1, lij_SSE1);
            lij2_SSE2          = _mm_mul_ps(lij_SSE2, lij_SSE2);
            lij2_SSE3          = _mm_mul_ps(lij_SSE3, lij_SSE3);
            lij3_SSE0          = _mm_mul_ps(lij2_SSE0, lij_SSE0);
            lij3_SSE1          = _mm_mul_ps(lij2_SSE1, lij_SSE1);
            lij3_SSE2          = _mm_mul_ps(lij2_SSE2, lij_SSE2);
            lij3_SSE3          = _mm_mul_ps(lij2_SSE3, lij_SSE3);

            diff2_SSE0         = _mm_sub_ps(uij2_SSE0, lij2_SSE0);
            diff2_SSE1         = _mm_sub_ps(uij2_SSE1, lij2_SSE1);
            diff2_SSE2         = _mm_sub_ps(uij2_SSE2, lij2_SSE2);
            diff2_SSE3         = _mm_sub_ps(uij2_SSE3, lij2_SSE3);
            lij_inv_SSE0       = gmx_mm_invsqrt_ps(lij2_SSE0);
            lij_inv_SSE1       = gmx_mm_invsqrt_ps(lij2_SSE1);
            lij_inv_SSE2       = gmx_mm_invsqrt_ps(lij2_SSE2);
            lij_inv_SSE3       = gmx_mm_invsqrt_ps(lij2_SSE3);
            sk2_aj_SSE         = _mm_mul_ps(sk_aj_SSE, sk_aj_SSE);
            sk2_rinv_SSE0      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE0);
            sk2_rinv_SSE1      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE1);
            sk2_rinv_SSE2      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE2);
            sk2_rinv_SSE3      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE3);
            prod_SSE0          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE0);
            prod_SSE1          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE1);
            prod_SSE2          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE2);
            prod_SSE3          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE3);

            logterm_SSE0       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE0, lij_inv_SSE0));
            logterm_SSE1       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE1, lij_inv_SSE1));
            logterm_SSE2       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE2, lij_inv_SSE2));
            logterm_SSE3       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE3, lij_inv_SSE3));

            t1_SSE0            = _mm_sub_ps(lij_SSE0, uij_SSE0);
            t1_SSE1            = _mm_sub_ps(lij_SSE1, uij_SSE1);
            t1_SSE2            = _mm_sub_ps(lij_SSE2, uij_SSE2);
            t1_SSE3            = _mm_sub_ps(lij_SSE3, uij_SSE3);
            t2_SSE0            = _mm_mul_ps(diff2_SSE0,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE0),
                                                       prod_SSE0));
            t2_SSE1            = _mm_mul_ps(diff2_SSE1,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE1),
                                                       prod_SSE1));
            t2_SSE2            = _mm_mul_ps(diff2_SSE2,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE2),
                                                       prod_SSE2));
            t2_SSE3            = _mm_mul_ps(diff2_SSE3,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE3),
                                                       prod_SSE3));

            t3_SSE0            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE0, logterm_SSE0));
            t3_SSE1            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE1, logterm_SSE1));
            t3_SSE2            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE2, logterm_SSE2));
            t3_SSE3            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE3, logterm_SSE3));
            t1_SSE0            = _mm_add_ps(t1_SSE0, _mm_add_ps(t2_SSE0, t3_SSE0));
            t1_SSE1            = _mm_add_ps(t1_SSE1, _mm_add_ps(t2_SSE1, t3_SSE1));
            t1_SSE2            = _mm_add_ps(t1_SSE2, _mm_add_ps(t2_SSE2, t3_SSE2));
            t1_SSE3            = _mm_add_ps(t1_SSE3, _mm_add_ps(t2_SSE3, t3_SSE3));
            t4_SSE0            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE0, lij_SSE0));
            t4_SSE1            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE1, lij_SSE1));
            t4_SSE2            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE2, lij_SSE2));
            t4_SSE3            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE3, lij_SSE3));
            t4_SSE0            = _mm_and_ps(t4_SSE0, obc_mask3_SSE0);
            t4_SSE1            = _mm_and_ps(t4_SSE1, obc_mask3_SSE1);
            t4_SSE2            = _mm_and_ps(t4_SSE2, obc_mask3_SSE2);
            t4_SSE3            = _mm_and_ps(t4_SSE3, obc_mask3_SSE3);
            t1_SSE0            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE0, t4_SSE0));
            t1_SSE1            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE1, t4_SSE1));
            t1_SSE2            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE2, t4_SSE2));
            t1_SSE3            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE3, t4_SSE3));

            sum_ai_SSE0        = _mm_add_ps(sum_ai_SSE0, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            sum_ai_SSE1        = _mm_add_ps(sum_ai_SSE1, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            sum_ai_SSE2        = _mm_add_ps(sum_ai_SSE2, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            sum_ai_SSE3        = _mm_add_ps(sum_ai_SSE3, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));

            t1_SSE0            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE0),
                                            _mm_mul_ps(prod_SSE0, lij3_SSE0));
            t1_SSE1            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE1),
                                            _mm_mul_ps(prod_SSE1, lij3_SSE1));
            t1_SSE2            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE2),
                                            _mm_mul_ps(prod_SSE2, lij3_SSE2));
            t1_SSE3            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE3),
                                            _mm_mul_ps(prod_SSE3, lij3_SSE3));
            t1_SSE0            = _mm_sub_ps(t1_SSE0,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE0, rinv_SSE0),
                                                                  _mm_mul_ps(lij3_SSE0, dr_SSE0))));
            t1_SSE1            = _mm_sub_ps(t1_SSE1,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE1, rinv_SSE1),
                                                                  _mm_mul_ps(lij3_SSE1, dr_SSE1))));
            t1_SSE2            = _mm_sub_ps(t1_SSE2,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE2, rinv_SSE2),
                                                                  _mm_mul_ps(lij3_SSE2, dr_SSE2))));
            t1_SSE3            = _mm_sub_ps(t1_SSE3,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE3, rinv_SSE3),
                                                                  _mm_mul_ps(lij3_SSE3, dr_SSE3))));

            t2_SSE0            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE0, rinv_SSE0),
                                                       _mm_mul_ps(uij3_SSE0, dr_SSE0)));
            t2_SSE1            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE1, rinv_SSE1),
                                                       _mm_mul_ps(uij3_SSE1, dr_SSE1)));
            t2_SSE2            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE2, rinv_SSE2),
                                                       _mm_mul_ps(uij3_SSE2, dr_SSE2)));
            t2_SSE3            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE3, rinv_SSE3),
                                                       _mm_mul_ps(uij3_SSE3, dr_SSE3)));
            t2_SSE0            = _mm_sub_ps(t2_SSE0,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE0),
                                                       _mm_mul_ps(prod_SSE0, uij3_SSE0)));
            t2_SSE1            = _mm_sub_ps(t2_SSE1,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE1),
                                                       _mm_mul_ps(prod_SSE1, uij3_SSE1)));
            t2_SSE2            = _mm_sub_ps(t2_SSE2,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE2),
                                                       _mm_mul_ps(prod_SSE2, uij3_SSE2)));
            t2_SSE3            = _mm_sub_ps(t2_SSE3,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE3),
                                                       _mm_mul_ps(prod_SSE3, uij3_SSE3)));
            t3_SSE0            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE0),
                                            _mm_mul_ps(rinv_SSE0, rinv_SSE0));
            t3_SSE1            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE1),
                                            _mm_mul_ps(rinv_SSE1, rinv_SSE1));
            t3_SSE2            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE2),
                                            _mm_mul_ps(rinv_SSE2, rinv_SSE2));
            t3_SSE3            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE3),
                                            _mm_mul_ps(rinv_SSE3, rinv_SSE3));
            t3_SSE0            = _mm_sub_ps(t3_SSE0,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE0, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE0, rinv_SSE0))));
            t3_SSE1            = _mm_sub_ps(t3_SSE1,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE1, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE1, rinv_SSE1))));
            t3_SSE2            = _mm_sub_ps(t3_SSE2,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE2, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE2, rinv_SSE2))));
            t3_SSE3            = _mm_sub_ps(t3_SSE3,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE3, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE3, rinv_SSE3))));

            t1_SSE0            = _mm_mul_ps(rinv_SSE0,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE0, t1_SSE0),
                                                       _mm_add_ps(t2_SSE0, t3_SSE0)));
            t1_SSE1            = _mm_mul_ps(rinv_SSE1,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE1, t1_SSE1),
                                                       _mm_add_ps(t2_SSE1, t3_SSE1)));
            t1_SSE2            = _mm_mul_ps(rinv_SSE2,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE2, t1_SSE2),
                                                       _mm_add_ps(t2_SSE2, t3_SSE2)));
            t1_SSE3            = _mm_mul_ps(rinv_SSE3,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE3, t1_SSE3),
                                                       _mm_add_ps(t2_SSE3, t3_SSE3)));

            _mm_store_ps(dadx, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));
            dadx += 4;

            /* Evaluate influence of atom ai -> aj */
            t1_SSE0            = _mm_add_ps(dr_SSE0, sk_ai_SSE0);
            t1_SSE1            = _mm_add_ps(dr_SSE1, sk_ai_SSE1);
            t1_SSE2            = _mm_add_ps(dr_SSE2, sk_ai_SSE2);
            t1_SSE3            = _mm_add_ps(dr_SSE3, sk_ai_SSE3);
            t2_SSE0            = _mm_sub_ps(dr_SSE0, sk_ai_SSE0);
            t2_SSE1            = _mm_sub_ps(dr_SSE1, sk_ai_SSE1);
            t2_SSE2            = _mm_sub_ps(dr_SSE2, sk_ai_SSE2);
            t2_SSE3            = _mm_sub_ps(dr_SSE3, sk_ai_SSE3);
            t3_SSE0            = _mm_sub_ps(sk_ai_SSE0, dr_SSE0);
            t3_SSE1            = _mm_sub_ps(sk_ai_SSE1, dr_SSE1);
            t3_SSE2            = _mm_sub_ps(sk_ai_SSE2, dr_SSE2);
            t3_SSE3            = _mm_sub_ps(sk_ai_SSE3, dr_SSE3);

            obc_mask1_SSE0     = _mm_cmplt_ps(raj_SSE, t1_SSE0);
            obc_mask1_SSE1     = _mm_cmplt_ps(raj_SSE, t1_SSE1);
            obc_mask1_SSE2     = _mm_cmplt_ps(raj_SSE, t1_SSE2);
            obc_mask1_SSE3     = _mm_cmplt_ps(raj_SSE, t1_SSE3);
            obc_mask2_SSE0     = _mm_cmplt_ps(raj_SSE, t2_SSE0);
            obc_mask2_SSE1     = _mm_cmplt_ps(raj_SSE, t2_SSE1);
            obc_mask2_SSE2     = _mm_cmplt_ps(raj_SSE, t2_SSE2);
            obc_mask2_SSE3     = _mm_cmplt_ps(raj_SSE, t2_SSE3);
            obc_mask3_SSE0     = _mm_cmplt_ps(raj_SSE, t3_SSE0);
            obc_mask3_SSE1     = _mm_cmplt_ps(raj_SSE, t3_SSE1);
            obc_mask3_SSE2     = _mm_cmplt_ps(raj_SSE, t3_SSE2);
            obc_mask3_SSE3     = _mm_cmplt_ps(raj_SSE, t3_SSE3);
            obc_mask1_SSE0     = _mm_and_ps(obc_mask1_SSE0, imask_SSE0);
            obc_mask1_SSE1     = _mm_and_ps(obc_mask1_SSE1, imask_SSE1);
            obc_mask1_SSE2     = _mm_and_ps(obc_mask1_SSE2, imask_SSE2);
            obc_mask1_SSE3     = _mm_and_ps(obc_mask1_SSE3, imask_SSE3);

            uij_SSE0           = gmx_mm_inv_ps(t1_SSE0);
            uij_SSE1           = gmx_mm_inv_ps(t1_SSE1);
            uij_SSE2           = gmx_mm_inv_ps(t1_SSE2);
            uij_SSE3           = gmx_mm_inv_ps(t1_SSE3);
            lij_SSE0           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE0, gmx_mm_inv_ps(t2_SSE0)),
                                              _mm_andnot_ps(obc_mask2_SSE0, raj_inv_SSE));
            lij_SSE1           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE1, gmx_mm_inv_ps(t2_SSE1)),
                                              _mm_andnot_ps(obc_mask2_SSE1, raj_inv_SSE));
            lij_SSE2           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE2, gmx_mm_inv_ps(t2_SSE2)),
                                              _mm_andnot_ps(obc_mask2_SSE2, raj_inv_SSE));
            lij_SSE3           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE3, gmx_mm_inv_ps(t2_SSE3)),
                                              _mm_andnot_ps(obc_mask2_SSE3, raj_inv_SSE));
            dlij_SSE0          = _mm_and_ps(one_SSE, obc_mask2_SSE0);
            dlij_SSE1          = _mm_and_ps(one_SSE, obc_mask2_SSE1);
            dlij_SSE2          = _mm_and_ps(one_SSE, obc_mask2_SSE2);
            dlij_SSE3          = _mm_and_ps(one_SSE, obc_mask2_SSE3);

            uij2_SSE0          = _mm_mul_ps(uij_SSE0, uij_SSE0);
            uij2_SSE1          = _mm_mul_ps(uij_SSE1, uij_SSE1);
            uij2_SSE2          = _mm_mul_ps(uij_SSE2, uij_SSE2);
            uij2_SSE3          = _mm_mul_ps(uij_SSE3, uij_SSE3);
            uij3_SSE0          = _mm_mul_ps(uij2_SSE0, uij_SSE0);
            uij3_SSE1          = _mm_mul_ps(uij2_SSE1, uij_SSE1);
            uij3_SSE2          = _mm_mul_ps(uij2_SSE2, uij_SSE2);
            uij3_SSE3          = _mm_mul_ps(uij2_SSE3, uij_SSE3);
            lij2_SSE0          = _mm_mul_ps(lij_SSE0, lij_SSE0);
            lij2_SSE1          = _mm_mul_ps(lij_SSE1, lij_SSE1);
            lij2_SSE2          = _mm_mul_ps(lij_SSE2, lij_SSE2);
            lij2_SSE3          = _mm_mul_ps(lij_SSE3, lij_SSE3);
            lij3_SSE0          = _mm_mul_ps(lij2_SSE0, lij_SSE0);
            lij3_SSE1          = _mm_mul_ps(lij2_SSE1, lij_SSE1);
            lij3_SSE2          = _mm_mul_ps(lij2_SSE2, lij_SSE2);
            lij3_SSE3          = _mm_mul_ps(lij2_SSE3, lij_SSE3);

            diff2_SSE0         = _mm_sub_ps(uij2_SSE0, lij2_SSE0);
            diff2_SSE1         = _mm_sub_ps(uij2_SSE1, lij2_SSE1);
            diff2_SSE2         = _mm_sub_ps(uij2_SSE2, lij2_SSE2);
            diff2_SSE3         = _mm_sub_ps(uij2_SSE3, lij2_SSE3);
            lij_inv_SSE0       = gmx_mm_invsqrt_ps(lij2_SSE0);
            lij_inv_SSE1       = gmx_mm_invsqrt_ps(lij2_SSE1);
            lij_inv_SSE2       = gmx_mm_invsqrt_ps(lij2_SSE2);
            lij_inv_SSE3       = gmx_mm_invsqrt_ps(lij2_SSE3);
            sk2_rinv_SSE0      = _mm_mul_ps(sk2_ai_SSE0, rinv_SSE0);
            sk2_rinv_SSE1      = _mm_mul_ps(sk2_ai_SSE1, rinv_SSE1);
            sk2_rinv_SSE2      = _mm_mul_ps(sk2_ai_SSE2, rinv_SSE2);
            sk2_rinv_SSE3      = _mm_mul_ps(sk2_ai_SSE3, rinv_SSE3);
            prod_SSE0          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE0);
            prod_SSE1          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE1);
            prod_SSE2          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE2);
            prod_SSE3          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE3);

            logterm_SSE0       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE0, lij_inv_SSE0));
            logterm_SSE1       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE1, lij_inv_SSE1));
            logterm_SSE2       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE2, lij_inv_SSE2));
            logterm_SSE3       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE3, lij_inv_SSE3));
            t1_SSE0            = _mm_sub_ps(lij_SSE0, uij_SSE0);
            t1_SSE1            = _mm_sub_ps(lij_SSE1, uij_SSE1);
            t1_SSE2            = _mm_sub_ps(lij_SSE2, uij_SSE2);
            t1_SSE3            = _mm_sub_ps(lij_SSE3, uij_SSE3);
            t2_SSE0            = _mm_mul_ps(diff2_SSE0,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE0),
                                                       prod_SSE0));
            t2_SSE1            = _mm_mul_ps(diff2_SSE1,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE1),
                                                       prod_SSE1));
            t2_SSE2            = _mm_mul_ps(diff2_SSE2,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE2),
                                                       prod_SSE2));
            t2_SSE3            = _mm_mul_ps(diff2_SSE3,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE3),
                                                       prod_SSE3));
            t3_SSE0            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE0, logterm_SSE0));
            t3_SSE1            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE1, logterm_SSE1));
            t3_SSE2            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE2, logterm_SSE2));
            t3_SSE3            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE3, logterm_SSE3));
            t1_SSE0            = _mm_add_ps(t1_SSE0, _mm_add_ps(t2_SSE0, t3_SSE0));
            t1_SSE1            = _mm_add_ps(t1_SSE1, _mm_add_ps(t2_SSE1, t3_SSE1));
            t1_SSE2            = _mm_add_ps(t1_SSE2, _mm_add_ps(t2_SSE2, t3_SSE2));
            t1_SSE3            = _mm_add_ps(t1_SSE3, _mm_add_ps(t2_SSE3, t3_SSE3));
            t4_SSE0            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE0));
            t4_SSE1            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE1));
            t4_SSE2            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE2));
            t4_SSE3            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE3));
            t4_SSE0            = _mm_and_ps(t4_SSE0, obc_mask3_SSE0);
            t4_SSE1            = _mm_and_ps(t4_SSE1, obc_mask3_SSE1);
            t4_SSE2            = _mm_and_ps(t4_SSE2, obc_mask3_SSE2);
            t4_SSE3            = _mm_and_ps(t4_SSE3, obc_mask3_SSE3);
            t1_SSE0            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE0, t4_SSE0));
            t1_SSE1            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE1, t4_SSE1));
            t1_SSE2            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE2, t4_SSE2));
            t1_SSE3            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE3, t4_SSE3));

            _mm_store_ps(work+j, _mm_add_ps(_mm_load_ps(work+j),
                                            gmx_mm_sum4_ps(_mm_and_ps(t1_SSE0, obc_mask1_SSE0),
                                                           _mm_and_ps(t1_SSE1, obc_mask1_SSE1),
                                                           _mm_and_ps(t1_SSE2, obc_mask1_SSE2),
                                                           _mm_and_ps(t1_SSE3, obc_mask1_SSE3))));

            t1_SSE0            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE0),
                                            _mm_mul_ps(prod_SSE0, lij3_SSE0));
            t1_SSE1            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE1),
                                            _mm_mul_ps(prod_SSE1, lij3_SSE1));
            t1_SSE2            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE2),
                                            _mm_mul_ps(prod_SSE2, lij3_SSE2));
            t1_SSE3            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE3),
                                            _mm_mul_ps(prod_SSE3, lij3_SSE3));
            t1_SSE0            = _mm_sub_ps(t1_SSE0,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE0, rinv_SSE0),
                                                                  _mm_mul_ps(lij3_SSE0, dr_SSE0))));
            t1_SSE1            = _mm_sub_ps(t1_SSE1,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE1, rinv_SSE1),
                                                                  _mm_mul_ps(lij3_SSE1, dr_SSE1))));
            t1_SSE2            = _mm_sub_ps(t1_SSE2,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE2, rinv_SSE2),
                                                                  _mm_mul_ps(lij3_SSE2, dr_SSE2))));
            t1_SSE3            = _mm_sub_ps(t1_SSE3,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE3, rinv_SSE3),
                                                                  _mm_mul_ps(lij3_SSE3, dr_SSE3))));
            t2_SSE0            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE0, rinv_SSE0),
                                                       _mm_mul_ps(uij3_SSE0, dr_SSE0)));
            t2_SSE1            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE1, rinv_SSE1),
                                                       _mm_mul_ps(uij3_SSE1, dr_SSE1)));
            t2_SSE2            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE2, rinv_SSE2),
                                                       _mm_mul_ps(uij3_SSE2, dr_SSE2)));
            t2_SSE3            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE3, rinv_SSE3),
                                                       _mm_mul_ps(uij3_SSE3, dr_SSE3)));
            t2_SSE0            = _mm_sub_ps(t2_SSE0,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE0),
                                                       _mm_mul_ps(prod_SSE0, uij3_SSE0)));
            t2_SSE1            = _mm_sub_ps(t2_SSE1,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE1),
                                                       _mm_mul_ps(prod_SSE1, uij3_SSE1)));
            t2_SSE2            = _mm_sub_ps(t2_SSE2,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE2),
                                                       _mm_mul_ps(prod_SSE2, uij3_SSE2)));
            t2_SSE3            = _mm_sub_ps(t2_SSE3,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE3),
                                                       _mm_mul_ps(prod_SSE3, uij3_SSE3)));

            t3_SSE0            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE0),
                                            _mm_mul_ps(rinv_SSE0, rinv_SSE0));
            t3_SSE1            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE1),
                                            _mm_mul_ps(rinv_SSE1, rinv_SSE1));
            t3_SSE2            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE2),
                                            _mm_mul_ps(rinv_SSE2, rinv_SSE2));
            t3_SSE3            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE3),
                                            _mm_mul_ps(rinv_SSE3, rinv_SSE3));

            t3_SSE0            = _mm_sub_ps(t3_SSE0,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE0, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE0, rinv_SSE0))));
            t3_SSE1            = _mm_sub_ps(t3_SSE1,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE1, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE1, rinv_SSE1))));
            t3_SSE2            = _mm_sub_ps(t3_SSE2,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE2, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE2, rinv_SSE2))));
            t3_SSE3            = _mm_sub_ps(t3_SSE3,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE3, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE3, rinv_SSE3))));

            t1_SSE0            = _mm_mul_ps(rinv_SSE0,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE0, t1_SSE0),
                                                       _mm_add_ps(t2_SSE0, t3_SSE0)));
            t1_SSE1            = _mm_mul_ps(rinv_SSE1,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE1, t1_SSE1),
                                                       _mm_add_ps(t2_SSE1, t3_SSE1)));
            t1_SSE2            = _mm_mul_ps(rinv_SSE2,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE2, t1_SSE2),
                                                       _mm_add_ps(t2_SSE2, t3_SSE2)));
            t1_SSE3            = _mm_mul_ps(rinv_SSE3,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE3, t1_SSE3),
                                                       _mm_add_ps(t2_SSE3, t3_SSE3)));

            _mm_store_ps(dadx, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));
            dadx += 4;
        }

        /* Epilogue part, including exclusion mask */
        for (j = nj2; j < nj3; j += UNROLLJ)
        {
            jmask_SSE0 = _mm_load_ps((real *)emask0);
            jmask_SSE1 = _mm_load_ps((real *)emask1);
            jmask_SSE2 = _mm_load_ps((real *)emask2);
            jmask_SSE3 = _mm_load_ps((real *)emask3);
            emask0    += UNROLLJ;
            emask1    += UNROLLJ;
            emask2    += UNROLLJ;
            emask3    += UNROLLJ;

            /* load j atom coordinates */
            jx_SSE            = _mm_load_ps(x_align+j);
            jy_SSE            = _mm_load_ps(y_align+j);
            jz_SSE            = _mm_load_ps(z_align+j);

            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0, jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0, jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0, jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1, jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1, jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1, jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2, jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2, jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2, jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3, jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3, jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1, dy_SSE1, dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2, dy_SSE2, dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3, dy_SSE3, dz_SSE3);

            /* Combine masks */
            jmask_SSE0         = _mm_and_ps(jmask_SSE0, imask_SSE0);
            jmask_SSE1         = _mm_and_ps(jmask_SSE1, imask_SSE1);
            jmask_SSE2         = _mm_and_ps(jmask_SSE2, imask_SSE2);
            jmask_SSE3         = _mm_and_ps(jmask_SSE3, imask_SSE3);

            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);

            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0, jmask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1, jmask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2, jmask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3, jmask_SSE3);

            dr_SSE0            = _mm_mul_ps(rsq_SSE0, rinv_SSE0);
            dr_SSE1            = _mm_mul_ps(rsq_SSE1, rinv_SSE1);
            dr_SSE2            = _mm_mul_ps(rsq_SSE2, rinv_SSE2);
            dr_SSE3            = _mm_mul_ps(rsq_SSE3, rinv_SSE3);

            sk_aj_SSE          = _mm_load_ps(obc_param+j);
            raj_SSE            = _mm_load_ps(gb_radius+j);

            raj_inv_SSE        = gmx_mm_inv_ps(raj_SSE);

            /* Evaluate influence of atom aj -> ai */
            t1_SSE0            = _mm_add_ps(dr_SSE0, sk_aj_SSE);
            t1_SSE1            = _mm_add_ps(dr_SSE1, sk_aj_SSE);
            t1_SSE2            = _mm_add_ps(dr_SSE2, sk_aj_SSE);
            t1_SSE3            = _mm_add_ps(dr_SSE3, sk_aj_SSE);
            t2_SSE0            = _mm_sub_ps(dr_SSE0, sk_aj_SSE);
            t2_SSE1            = _mm_sub_ps(dr_SSE1, sk_aj_SSE);
            t2_SSE2            = _mm_sub_ps(dr_SSE2, sk_aj_SSE);
            t2_SSE3            = _mm_sub_ps(dr_SSE3, sk_aj_SSE);
            t3_SSE0            = _mm_sub_ps(sk_aj_SSE, dr_SSE0);
            t3_SSE1            = _mm_sub_ps(sk_aj_SSE, dr_SSE1);
            t3_SSE2            = _mm_sub_ps(sk_aj_SSE, dr_SSE2);
            t3_SSE3            = _mm_sub_ps(sk_aj_SSE, dr_SSE3);

            obc_mask1_SSE0     = _mm_cmplt_ps(rai_SSE0, t1_SSE0);
            obc_mask1_SSE1     = _mm_cmplt_ps(rai_SSE1, t1_SSE1);
            obc_mask1_SSE2     = _mm_cmplt_ps(rai_SSE2, t1_SSE2);
            obc_mask1_SSE3     = _mm_cmplt_ps(rai_SSE3, t1_SSE3);
            obc_mask2_SSE0     = _mm_cmplt_ps(rai_SSE0, t2_SSE0);
            obc_mask2_SSE1     = _mm_cmplt_ps(rai_SSE1, t2_SSE1);
            obc_mask2_SSE2     = _mm_cmplt_ps(rai_SSE2, t2_SSE2);
            obc_mask2_SSE3     = _mm_cmplt_ps(rai_SSE3, t2_SSE3);
            obc_mask3_SSE0     = _mm_cmplt_ps(rai_SSE0, t3_SSE0);
            obc_mask3_SSE1     = _mm_cmplt_ps(rai_SSE1, t3_SSE1);
            obc_mask3_SSE2     = _mm_cmplt_ps(rai_SSE2, t3_SSE2);
            obc_mask3_SSE3     = _mm_cmplt_ps(rai_SSE3, t3_SSE3);
            obc_mask1_SSE0     = _mm_and_ps(obc_mask1_SSE0, jmask_SSE0);
            obc_mask1_SSE1     = _mm_and_ps(obc_mask1_SSE1, jmask_SSE1);
            obc_mask1_SSE2     = _mm_and_ps(obc_mask1_SSE2, jmask_SSE2);
            obc_mask1_SSE3     = _mm_and_ps(obc_mask1_SSE3, jmask_SSE3);

            uij_SSE0           = gmx_mm_inv_ps(t1_SSE0);
            uij_SSE1           = gmx_mm_inv_ps(t1_SSE1);
            uij_SSE2           = gmx_mm_inv_ps(t1_SSE2);
            uij_SSE3           = gmx_mm_inv_ps(t1_SSE3);
            lij_SSE0           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE0, gmx_mm_inv_ps(t2_SSE0)),
                                              _mm_andnot_ps(obc_mask2_SSE0, rai_inv_SSE0));
            lij_SSE1           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE1, gmx_mm_inv_ps(t2_SSE1)),
                                              _mm_andnot_ps(obc_mask2_SSE1, rai_inv_SSE1));
            lij_SSE2           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE2, gmx_mm_inv_ps(t2_SSE2)),
                                              _mm_andnot_ps(obc_mask2_SSE2, rai_inv_SSE2));
            lij_SSE3           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE3, gmx_mm_inv_ps(t2_SSE3)),
                                              _mm_andnot_ps(obc_mask2_SSE3, rai_inv_SSE3));
            dlij_SSE0          = _mm_and_ps(one_SSE, obc_mask2_SSE0);
            dlij_SSE1          = _mm_and_ps(one_SSE, obc_mask2_SSE1);
            dlij_SSE2          = _mm_and_ps(one_SSE, obc_mask2_SSE2);
            dlij_SSE3          = _mm_and_ps(one_SSE, obc_mask2_SSE3);

            uij2_SSE0          = _mm_mul_ps(uij_SSE0, uij_SSE0);
            uij2_SSE1          = _mm_mul_ps(uij_SSE1, uij_SSE1);
            uij2_SSE2          = _mm_mul_ps(uij_SSE2, uij_SSE2);
            uij2_SSE3          = _mm_mul_ps(uij_SSE3, uij_SSE3);
            uij3_SSE0          = _mm_mul_ps(uij2_SSE0, uij_SSE0);
            uij3_SSE1          = _mm_mul_ps(uij2_SSE1, uij_SSE1);
            uij3_SSE2          = _mm_mul_ps(uij2_SSE2, uij_SSE2);
            uij3_SSE3          = _mm_mul_ps(uij2_SSE3, uij_SSE3);
            lij2_SSE0          = _mm_mul_ps(lij_SSE0, lij_SSE0);
            lij2_SSE1          = _mm_mul_ps(lij_SSE1, lij_SSE1);
            lij2_SSE2          = _mm_mul_ps(lij_SSE2, lij_SSE2);
            lij2_SSE3          = _mm_mul_ps(lij_SSE3, lij_SSE3);
            lij3_SSE0          = _mm_mul_ps(lij2_SSE0, lij_SSE0);
            lij3_SSE1          = _mm_mul_ps(lij2_SSE1, lij_SSE1);
            lij3_SSE2          = _mm_mul_ps(lij2_SSE2, lij_SSE2);
            lij3_SSE3          = _mm_mul_ps(lij2_SSE3, lij_SSE3);

            diff2_SSE0         = _mm_sub_ps(uij2_SSE0, lij2_SSE0);
            diff2_SSE1         = _mm_sub_ps(uij2_SSE1, lij2_SSE1);
            diff2_SSE2         = _mm_sub_ps(uij2_SSE2, lij2_SSE2);
            diff2_SSE3         = _mm_sub_ps(uij2_SSE3, lij2_SSE3);
            lij_inv_SSE0       = gmx_mm_invsqrt_ps(lij2_SSE0);
            lij_inv_SSE1       = gmx_mm_invsqrt_ps(lij2_SSE1);
            lij_inv_SSE2       = gmx_mm_invsqrt_ps(lij2_SSE2);
            lij_inv_SSE3       = gmx_mm_invsqrt_ps(lij2_SSE3);
            sk2_aj_SSE         = _mm_mul_ps(sk_aj_SSE, sk_aj_SSE);
            sk2_rinv_SSE0      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE0);
            sk2_rinv_SSE1      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE1);
            sk2_rinv_SSE2      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE2);
            sk2_rinv_SSE3      = _mm_mul_ps(sk2_aj_SSE, rinv_SSE3);
            prod_SSE0          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE0);
            prod_SSE1          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE1);
            prod_SSE2          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE2);
            prod_SSE3          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE3);

            logterm_SSE0       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE0, lij_inv_SSE0));
            logterm_SSE1       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE1, lij_inv_SSE1));
            logterm_SSE2       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE2, lij_inv_SSE2));
            logterm_SSE3       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE3, lij_inv_SSE3));

            t1_SSE0            = _mm_sub_ps(lij_SSE0, uij_SSE0);
            t1_SSE1            = _mm_sub_ps(lij_SSE1, uij_SSE1);
            t1_SSE2            = _mm_sub_ps(lij_SSE2, uij_SSE2);
            t1_SSE3            = _mm_sub_ps(lij_SSE3, uij_SSE3);
            t2_SSE0            = _mm_mul_ps(diff2_SSE0,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE0),
                                                       prod_SSE0));
            t2_SSE1            = _mm_mul_ps(diff2_SSE1,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE1),
                                                       prod_SSE1));
            t2_SSE2            = _mm_mul_ps(diff2_SSE2,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE2),
                                                       prod_SSE2));
            t2_SSE3            = _mm_mul_ps(diff2_SSE3,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE3),
                                                       prod_SSE3));

            t3_SSE0            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE0, logterm_SSE0));
            t3_SSE1            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE1, logterm_SSE1));
            t3_SSE2            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE2, logterm_SSE2));
            t3_SSE3            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE3, logterm_SSE3));
            t1_SSE0            = _mm_add_ps(t1_SSE0, _mm_add_ps(t2_SSE0, t3_SSE0));
            t1_SSE1            = _mm_add_ps(t1_SSE1, _mm_add_ps(t2_SSE1, t3_SSE1));
            t1_SSE2            = _mm_add_ps(t1_SSE2, _mm_add_ps(t2_SSE2, t3_SSE2));
            t1_SSE3            = _mm_add_ps(t1_SSE3, _mm_add_ps(t2_SSE3, t3_SSE3));
            t4_SSE0            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE0, lij_SSE0));
            t4_SSE1            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE1, lij_SSE1));
            t4_SSE2            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE2, lij_SSE2));
            t4_SSE3            = _mm_mul_ps(two_SSE, _mm_sub_ps(rai_inv_SSE3, lij_SSE3));
            t4_SSE0            = _mm_and_ps(t4_SSE0, obc_mask3_SSE0);
            t4_SSE1            = _mm_and_ps(t4_SSE1, obc_mask3_SSE1);
            t4_SSE2            = _mm_and_ps(t4_SSE2, obc_mask3_SSE2);
            t4_SSE3            = _mm_and_ps(t4_SSE3, obc_mask3_SSE3);
            t1_SSE0            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE0, t4_SSE0));
            t1_SSE1            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE1, t4_SSE1));
            t1_SSE2            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE2, t4_SSE2));
            t1_SSE3            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE3, t4_SSE3));

            sum_ai_SSE0        = _mm_add_ps(sum_ai_SSE0, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            sum_ai_SSE1        = _mm_add_ps(sum_ai_SSE1, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            sum_ai_SSE2        = _mm_add_ps(sum_ai_SSE2, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            sum_ai_SSE3        = _mm_add_ps(sum_ai_SSE3, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));

            t1_SSE0            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE0),
                                            _mm_mul_ps(prod_SSE0, lij3_SSE0));
            t1_SSE1            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE1),
                                            _mm_mul_ps(prod_SSE1, lij3_SSE1));
            t1_SSE2            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE2),
                                            _mm_mul_ps(prod_SSE2, lij3_SSE2));
            t1_SSE3            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE3),
                                            _mm_mul_ps(prod_SSE3, lij3_SSE3));
            t1_SSE0            = _mm_sub_ps(t1_SSE0,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE0, rinv_SSE0),
                                                                  _mm_mul_ps(lij3_SSE0, dr_SSE0))));
            t1_SSE1            = _mm_sub_ps(t1_SSE1,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE1, rinv_SSE1),
                                                                  _mm_mul_ps(lij3_SSE1, dr_SSE1))));
            t1_SSE2            = _mm_sub_ps(t1_SSE2,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE2, rinv_SSE2),
                                                                  _mm_mul_ps(lij3_SSE2, dr_SSE2))));
            t1_SSE3            = _mm_sub_ps(t1_SSE3,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE3, rinv_SSE3),
                                                                  _mm_mul_ps(lij3_SSE3, dr_SSE3))));

            t2_SSE0            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE0, rinv_SSE0),
                                                       _mm_mul_ps(uij3_SSE0, dr_SSE0)));
            t2_SSE1            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE1, rinv_SSE1),
                                                       _mm_mul_ps(uij3_SSE1, dr_SSE1)));
            t2_SSE2            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE2, rinv_SSE2),
                                                       _mm_mul_ps(uij3_SSE2, dr_SSE2)));
            t2_SSE3            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE3, rinv_SSE3),
                                                       _mm_mul_ps(uij3_SSE3, dr_SSE3)));
            t2_SSE0            = _mm_sub_ps(t2_SSE0,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE0),
                                                       _mm_mul_ps(prod_SSE0, uij3_SSE0)));
            t2_SSE1            = _mm_sub_ps(t2_SSE1,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE1),
                                                       _mm_mul_ps(prod_SSE1, uij3_SSE1)));
            t2_SSE2            = _mm_sub_ps(t2_SSE2,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE2),
                                                       _mm_mul_ps(prod_SSE2, uij3_SSE2)));
            t2_SSE3            = _mm_sub_ps(t2_SSE3,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE3),
                                                       _mm_mul_ps(prod_SSE3, uij3_SSE3)));
            t3_SSE0            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE0),
                                            _mm_mul_ps(rinv_SSE0, rinv_SSE0));
            t3_SSE1            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE1),
                                            _mm_mul_ps(rinv_SSE1, rinv_SSE1));
            t3_SSE2            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE2),
                                            _mm_mul_ps(rinv_SSE2, rinv_SSE2));
            t3_SSE3            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE3),
                                            _mm_mul_ps(rinv_SSE3, rinv_SSE3));
            t3_SSE0            = _mm_sub_ps(t3_SSE0,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE0, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE0, rinv_SSE0))));
            t3_SSE1            = _mm_sub_ps(t3_SSE1,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE1, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE1, rinv_SSE1))));
            t3_SSE2            = _mm_sub_ps(t3_SSE2,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE2, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE2, rinv_SSE2))));
            t3_SSE3            = _mm_sub_ps(t3_SSE3,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE3, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE3, rinv_SSE3))));

            t1_SSE0            = _mm_mul_ps(rinv_SSE0,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE0, t1_SSE0),
                                                       _mm_add_ps(t2_SSE0, t3_SSE0)));
            t1_SSE1            = _mm_mul_ps(rinv_SSE1,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE1, t1_SSE1),
                                                       _mm_add_ps(t2_SSE1, t3_SSE1)));
            t1_SSE2            = _mm_mul_ps(rinv_SSE2,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE2, t1_SSE2),
                                                       _mm_add_ps(t2_SSE2, t3_SSE2)));
            t1_SSE3            = _mm_mul_ps(rinv_SSE3,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE3, t1_SSE3),
                                                       _mm_add_ps(t2_SSE3, t3_SSE3)));

            _mm_store_ps(dadx, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));
            dadx += 4;

            /* Evaluate influence of atom ai -> aj */
            t1_SSE0            = _mm_add_ps(dr_SSE0, sk_ai_SSE0);
            t1_SSE1            = _mm_add_ps(dr_SSE1, sk_ai_SSE1);
            t1_SSE2            = _mm_add_ps(dr_SSE2, sk_ai_SSE2);
            t1_SSE3            = _mm_add_ps(dr_SSE3, sk_ai_SSE3);
            t2_SSE0            = _mm_sub_ps(dr_SSE0, sk_ai_SSE0);
            t2_SSE1            = _mm_sub_ps(dr_SSE1, sk_ai_SSE1);
            t2_SSE2            = _mm_sub_ps(dr_SSE2, sk_ai_SSE2);
            t2_SSE3            = _mm_sub_ps(dr_SSE3, sk_ai_SSE3);
            t3_SSE0            = _mm_sub_ps(sk_ai_SSE0, dr_SSE0);
            t3_SSE1            = _mm_sub_ps(sk_ai_SSE1, dr_SSE1);
            t3_SSE2            = _mm_sub_ps(sk_ai_SSE2, dr_SSE2);
            t3_SSE3            = _mm_sub_ps(sk_ai_SSE3, dr_SSE3);

            obc_mask1_SSE0     = _mm_cmplt_ps(raj_SSE, t1_SSE0);
            obc_mask1_SSE1     = _mm_cmplt_ps(raj_SSE, t1_SSE1);
            obc_mask1_SSE2     = _mm_cmplt_ps(raj_SSE, t1_SSE2);
            obc_mask1_SSE3     = _mm_cmplt_ps(raj_SSE, t1_SSE3);
            obc_mask2_SSE0     = _mm_cmplt_ps(raj_SSE, t2_SSE0);
            obc_mask2_SSE1     = _mm_cmplt_ps(raj_SSE, t2_SSE1);
            obc_mask2_SSE2     = _mm_cmplt_ps(raj_SSE, t2_SSE2);
            obc_mask2_SSE3     = _mm_cmplt_ps(raj_SSE, t2_SSE3);
            obc_mask3_SSE0     = _mm_cmplt_ps(raj_SSE, t3_SSE0);
            obc_mask3_SSE1     = _mm_cmplt_ps(raj_SSE, t3_SSE1);
            obc_mask3_SSE2     = _mm_cmplt_ps(raj_SSE, t3_SSE2);
            obc_mask3_SSE3     = _mm_cmplt_ps(raj_SSE, t3_SSE3);
            obc_mask1_SSE0     = _mm_and_ps(obc_mask1_SSE0, jmask_SSE0);
            obc_mask1_SSE1     = _mm_and_ps(obc_mask1_SSE1, jmask_SSE1);
            obc_mask1_SSE2     = _mm_and_ps(obc_mask1_SSE2, jmask_SSE2);
            obc_mask1_SSE3     = _mm_and_ps(obc_mask1_SSE3, jmask_SSE3);

            uij_SSE0           = gmx_mm_inv_ps(t1_SSE0);
            uij_SSE1           = gmx_mm_inv_ps(t1_SSE1);
            uij_SSE2           = gmx_mm_inv_ps(t1_SSE2);
            uij_SSE3           = gmx_mm_inv_ps(t1_SSE3);
            lij_SSE0           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE0, gmx_mm_inv_ps(t2_SSE0)),
                                              _mm_andnot_ps(obc_mask2_SSE0, raj_inv_SSE));
            lij_SSE1           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE1, gmx_mm_inv_ps(t2_SSE1)),
                                              _mm_andnot_ps(obc_mask2_SSE1, raj_inv_SSE));
            lij_SSE2           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE2, gmx_mm_inv_ps(t2_SSE2)),
                                              _mm_andnot_ps(obc_mask2_SSE2, raj_inv_SSE));
            lij_SSE3           = _mm_or_ps(   _mm_and_ps(obc_mask2_SSE3, gmx_mm_inv_ps(t2_SSE3)),
                                              _mm_andnot_ps(obc_mask2_SSE3, raj_inv_SSE));
            dlij_SSE0          = _mm_and_ps(one_SSE, obc_mask2_SSE0);
            dlij_SSE1          = _mm_and_ps(one_SSE, obc_mask2_SSE1);
            dlij_SSE2          = _mm_and_ps(one_SSE, obc_mask2_SSE2);
            dlij_SSE3          = _mm_and_ps(one_SSE, obc_mask2_SSE3);

            uij2_SSE0          = _mm_mul_ps(uij_SSE0, uij_SSE0);
            uij2_SSE1          = _mm_mul_ps(uij_SSE1, uij_SSE1);
            uij2_SSE2          = _mm_mul_ps(uij_SSE2, uij_SSE2);
            uij2_SSE3          = _mm_mul_ps(uij_SSE3, uij_SSE3);
            uij3_SSE0          = _mm_mul_ps(uij2_SSE0, uij_SSE0);
            uij3_SSE1          = _mm_mul_ps(uij2_SSE1, uij_SSE1);
            uij3_SSE2          = _mm_mul_ps(uij2_SSE2, uij_SSE2);
            uij3_SSE3          = _mm_mul_ps(uij2_SSE3, uij_SSE3);
            lij2_SSE0          = _mm_mul_ps(lij_SSE0, lij_SSE0);
            lij2_SSE1          = _mm_mul_ps(lij_SSE1, lij_SSE1);
            lij2_SSE2          = _mm_mul_ps(lij_SSE2, lij_SSE2);
            lij2_SSE3          = _mm_mul_ps(lij_SSE3, lij_SSE3);
            lij3_SSE0          = _mm_mul_ps(lij2_SSE0, lij_SSE0);
            lij3_SSE1          = _mm_mul_ps(lij2_SSE1, lij_SSE1);
            lij3_SSE2          = _mm_mul_ps(lij2_SSE2, lij_SSE2);
            lij3_SSE3          = _mm_mul_ps(lij2_SSE3, lij_SSE3);

            diff2_SSE0         = _mm_sub_ps(uij2_SSE0, lij2_SSE0);
            diff2_SSE1         = _mm_sub_ps(uij2_SSE1, lij2_SSE1);
            diff2_SSE2         = _mm_sub_ps(uij2_SSE2, lij2_SSE2);
            diff2_SSE3         = _mm_sub_ps(uij2_SSE3, lij2_SSE3);
            lij_inv_SSE0       = gmx_mm_invsqrt_ps(lij2_SSE0);
            lij_inv_SSE1       = gmx_mm_invsqrt_ps(lij2_SSE1);
            lij_inv_SSE2       = gmx_mm_invsqrt_ps(lij2_SSE2);
            lij_inv_SSE3       = gmx_mm_invsqrt_ps(lij2_SSE3);
            sk2_rinv_SSE0      = _mm_mul_ps(sk2_ai_SSE0, rinv_SSE0);
            sk2_rinv_SSE1      = _mm_mul_ps(sk2_ai_SSE1, rinv_SSE1);
            sk2_rinv_SSE2      = _mm_mul_ps(sk2_ai_SSE2, rinv_SSE2);
            sk2_rinv_SSE3      = _mm_mul_ps(sk2_ai_SSE3, rinv_SSE3);
            prod_SSE0          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE0);
            prod_SSE1          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE1);
            prod_SSE2          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE2);
            prod_SSE3          = _mm_mul_ps(onefourth_SSE, sk2_rinv_SSE3);

            logterm_SSE0       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE0, lij_inv_SSE0));
            logterm_SSE1       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE1, lij_inv_SSE1));
            logterm_SSE2       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE2, lij_inv_SSE2));
            logterm_SSE3       = gmx_mm_log_ps(_mm_mul_ps(uij_SSE3, lij_inv_SSE3));
            t1_SSE0            = _mm_sub_ps(lij_SSE0, uij_SSE0);
            t1_SSE1            = _mm_sub_ps(lij_SSE1, uij_SSE1);
            t1_SSE2            = _mm_sub_ps(lij_SSE2, uij_SSE2);
            t1_SSE3            = _mm_sub_ps(lij_SSE3, uij_SSE3);
            t2_SSE0            = _mm_mul_ps(diff2_SSE0,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE0),
                                                       prod_SSE0));
            t2_SSE1            = _mm_mul_ps(diff2_SSE1,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE1),
                                                       prod_SSE1));
            t2_SSE2            = _mm_mul_ps(diff2_SSE2,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE2),
                                                       prod_SSE2));
            t2_SSE3            = _mm_mul_ps(diff2_SSE3,
                                            _mm_sub_ps(_mm_mul_ps(onefourth_SSE, dr_SSE3),
                                                       prod_SSE3));
            t3_SSE0            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE0, logterm_SSE0));
            t3_SSE1            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE1, logterm_SSE1));
            t3_SSE2            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE2, logterm_SSE2));
            t3_SSE3            = _mm_mul_ps(half_SSE, _mm_mul_ps(rinv_SSE3, logterm_SSE3));
            t1_SSE0            = _mm_add_ps(t1_SSE0, _mm_add_ps(t2_SSE0, t3_SSE0));
            t1_SSE1            = _mm_add_ps(t1_SSE1, _mm_add_ps(t2_SSE1, t3_SSE1));
            t1_SSE2            = _mm_add_ps(t1_SSE2, _mm_add_ps(t2_SSE2, t3_SSE2));
            t1_SSE3            = _mm_add_ps(t1_SSE3, _mm_add_ps(t2_SSE3, t3_SSE3));
            t4_SSE0            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE0));
            t4_SSE1            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE1));
            t4_SSE2            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE2));
            t4_SSE3            = _mm_mul_ps(two_SSE, _mm_sub_ps(raj_inv_SSE, lij_SSE3));
            t4_SSE0            = _mm_and_ps(t4_SSE0, obc_mask3_SSE0);
            t4_SSE1            = _mm_and_ps(t4_SSE1, obc_mask3_SSE1);
            t4_SSE2            = _mm_and_ps(t4_SSE2, obc_mask3_SSE2);
            t4_SSE3            = _mm_and_ps(t4_SSE3, obc_mask3_SSE3);
            t1_SSE0            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE0, t4_SSE0));
            t1_SSE1            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE1, t4_SSE1));
            t1_SSE2            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE2, t4_SSE2));
            t1_SSE3            = _mm_mul_ps(half_SSE, _mm_add_ps(t1_SSE3, t4_SSE3));

            _mm_store_ps(work+j, _mm_add_ps(_mm_load_ps(work+j),
                                            gmx_mm_sum4_ps(_mm_and_ps(t1_SSE0, obc_mask1_SSE0),
                                                           _mm_and_ps(t1_SSE1, obc_mask1_SSE1),
                                                           _mm_and_ps(t1_SSE2, obc_mask1_SSE2),
                                                           _mm_and_ps(t1_SSE3, obc_mask1_SSE3))));

            t1_SSE0            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE0),
                                            _mm_mul_ps(prod_SSE0, lij3_SSE0));
            t1_SSE1            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE1),
                                            _mm_mul_ps(prod_SSE1, lij3_SSE1));
            t1_SSE2            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE2),
                                            _mm_mul_ps(prod_SSE2, lij3_SSE2));
            t1_SSE3            = _mm_add_ps(_mm_mul_ps(half_SSE, lij2_SSE3),
                                            _mm_mul_ps(prod_SSE3, lij3_SSE3));
            t1_SSE0            = _mm_sub_ps(t1_SSE0,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE0, rinv_SSE0),
                                                                  _mm_mul_ps(lij3_SSE0, dr_SSE0))));
            t1_SSE1            = _mm_sub_ps(t1_SSE1,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE1, rinv_SSE1),
                                                                  _mm_mul_ps(lij3_SSE1, dr_SSE1))));
            t1_SSE2            = _mm_sub_ps(t1_SSE2,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE2, rinv_SSE2),
                                                                  _mm_mul_ps(lij3_SSE2, dr_SSE2))));
            t1_SSE3            = _mm_sub_ps(t1_SSE3,
                                            _mm_mul_ps(onefourth_SSE,
                                                       _mm_add_ps(_mm_mul_ps(lij_SSE3, rinv_SSE3),
                                                                  _mm_mul_ps(lij3_SSE3, dr_SSE3))));
            t2_SSE0            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE0, rinv_SSE0),
                                                       _mm_mul_ps(uij3_SSE0, dr_SSE0)));
            t2_SSE1            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE1, rinv_SSE1),
                                                       _mm_mul_ps(uij3_SSE1, dr_SSE1)));
            t2_SSE2            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE2, rinv_SSE2),
                                                       _mm_mul_ps(uij3_SSE2, dr_SSE2)));
            t2_SSE3            = _mm_mul_ps(onefourth_SSE,
                                            _mm_add_ps(_mm_mul_ps(uij_SSE3, rinv_SSE3),
                                                       _mm_mul_ps(uij3_SSE3, dr_SSE3)));
            t2_SSE0            = _mm_sub_ps(t2_SSE0,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE0),
                                                       _mm_mul_ps(prod_SSE0, uij3_SSE0)));
            t2_SSE1            = _mm_sub_ps(t2_SSE1,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE1),
                                                       _mm_mul_ps(prod_SSE1, uij3_SSE1)));
            t2_SSE2            = _mm_sub_ps(t2_SSE2,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE2),
                                                       _mm_mul_ps(prod_SSE2, uij3_SSE2)));
            t2_SSE3            = _mm_sub_ps(t2_SSE3,
                                            _mm_add_ps(_mm_mul_ps(half_SSE, uij2_SSE3),
                                                       _mm_mul_ps(prod_SSE3, uij3_SSE3)));

            t3_SSE0            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE0),
                                            _mm_mul_ps(rinv_SSE0, rinv_SSE0));
            t3_SSE1            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE1),
                                            _mm_mul_ps(rinv_SSE1, rinv_SSE1));
            t3_SSE2            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE2),
                                            _mm_mul_ps(rinv_SSE2, rinv_SSE2));
            t3_SSE3            = _mm_mul_ps(_mm_mul_ps(onefourth_SSE, logterm_SSE3),
                                            _mm_mul_ps(rinv_SSE3, rinv_SSE3));

            t3_SSE0            = _mm_sub_ps(t3_SSE0,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE0, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE0, rinv_SSE0))));
            t3_SSE1            = _mm_sub_ps(t3_SSE1,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE1, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE1, rinv_SSE1))));
            t3_SSE2            = _mm_sub_ps(t3_SSE2,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE2, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE2, rinv_SSE2))));
            t3_SSE3            = _mm_sub_ps(t3_SSE3,
                                            _mm_mul_ps(_mm_mul_ps(diff2_SSE3, oneeighth_SSE),
                                                       _mm_add_ps(one_SSE,
                                                                  _mm_mul_ps(sk2_rinv_SSE3, rinv_SSE3))));


            t1_SSE0            = _mm_mul_ps(rinv_SSE0,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE0, t1_SSE0),
                                                       _mm_add_ps(t2_SSE0, t3_SSE0)));
            t1_SSE1            = _mm_mul_ps(rinv_SSE1,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE1, t1_SSE1),
                                                       _mm_add_ps(t2_SSE1, t3_SSE1)));
            t1_SSE2            = _mm_mul_ps(rinv_SSE2,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE2, t1_SSE2),
                                                       _mm_add_ps(t2_SSE2, t3_SSE2)));
            t1_SSE3            = _mm_mul_ps(rinv_SSE3,
                                            _mm_add_ps(_mm_mul_ps(dlij_SSE3, t1_SSE3),
                                                       _mm_add_ps(t2_SSE3, t3_SSE3)));

            _mm_store_ps(dadx, _mm_and_ps(t1_SSE0, obc_mask1_SSE0));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE1, obc_mask1_SSE1));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE2, obc_mask1_SSE2));
            dadx += 4;
            _mm_store_ps(dadx, _mm_and_ps(t1_SSE3, obc_mask1_SSE3));
            dadx += 4;
        }
        _MM_TRANSPOSE4_PS(sum_ai_SSE0, sum_ai_SSE1, sum_ai_SSE2, sum_ai_SSE3);
        sum_ai_SSE0 = _mm_add_ps(sum_ai_SSE0, sum_ai_SSE1);
        sum_ai_SSE2 = _mm_add_ps(sum_ai_SSE2, sum_ai_SSE3);
        sum_ai_SSE0 = _mm_add_ps(sum_ai_SSE0, sum_ai_SSE2);
        _mm_store_ps(work+i, _mm_add_ps(sum_ai_SSE0, _mm_load_ps(work+i)));
    }


    for (i = 0; i < natoms/2+1; i++)
    {
        work[i] += work[natoms+i];
    }

    /* Parallel summations would go here if ever implemented with DD */

    if (gb_algorithm == egbHCT)
    {
        /* HCT */
        for (i = 0; i < natoms; i++)
        {
            if (born->use[i] != 0)
            {
                rai     = top->atomtypes.gb_radius[mdatoms->typeA[i]]-born->gb_doffset;
                sum_ai  = 1.0/rai - work[i];
                min_rad = rai + born->gb_doffset;
                rad     = 1.0/sum_ai;

                born->bRad[i]   = rad > min_rad ? rad : min_rad;
                fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);
            }
        }

    }
    else
    {
        /* OBC */

        /* Calculate the radii */
        for (i = 0; i < natoms; i++)
        {

            if (born->use[i] != 0)
            {
                rai        = top->atomtypes.gb_radius[mdatoms->typeA[i]];
                rai_inv2   = 1.0/rai;
                rai        = rai-born->gb_doffset;
                rai_inv    = 1.0/rai;
                sum_ai     = rai * work[i];
                sum_ai2    = sum_ai  * sum_ai;
                sum_ai3    = sum_ai2 * sum_ai;

                tsum          = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
                born->bRad[i] = rai_inv - tsum*rai_inv2;
                born->bRad[i] = 1.0 / born->bRad[i];

                fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);

                tchain         = rai * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
                born->drobc[i] = (1.0-tsum*tsum)*tchain*rai_inv2;
            }
        }
    }

    return 0;
}








int
genborn_allvsall_calc_chainrule_sse2_single(t_forcerec *           fr,
                                            t_mdatoms *            mdatoms,
                                            gmx_genborn_t *        born,
                                            real *                 x,
                                            real *                 f,
                                            int                    gb_algorithm,
                                            void *                 paadata)
{
    gmx_allvsallgb2_data_t *aadata;
    int                     natoms;
    int                     ni0, ni1;
    int                     nj0, nj1, nj2, nj3;
    int                     i, j, k, n;
    int                     idx;
    int              *      mask;
    int              *      pmask0;
    int              *      emask0;
    int              *      jindex;

    real                    ix, iy, iz;
    real                    fix, fiy, fiz;
    real                    jx, jy, jz;
    real                    dx, dy, dz;
    real                    tx, ty, tz;
    real                    rbai, rbaj, fgb, fgb_ai, rbi;
    real              *     rb;
    real              *     dadx;
    real              *     x_align;
    real              *     y_align;
    real              *     z_align;
    real              *     fx_align;
    real              *     fy_align;
    real              *     fz_align;
    real                    tmpsum[4];

    __m128                  jmask_SSE0, jmask_SSE1, jmask_SSE2, jmask_SSE3;
    __m128                  ix_SSE0, iy_SSE0, iz_SSE0;
    __m128                  ix_SSE1, iy_SSE1, iz_SSE1;
    __m128                  ix_SSE2, iy_SSE2, iz_SSE2;
    __m128                  ix_SSE3, iy_SSE3, iz_SSE3;
    __m128                  fix_SSE0, fiy_SSE0, fiz_SSE0;
    __m128                  fix_SSE1, fiy_SSE1, fiz_SSE1;
    __m128                  fix_SSE2, fiy_SSE2, fiz_SSE2;
    __m128                  fix_SSE3, fiy_SSE3, fiz_SSE3;
    __m128                  rbai_SSE0, rbai_SSE1, rbai_SSE2, rbai_SSE3;
    __m128                  imask_SSE0, imask_SSE1, imask_SSE2, imask_SSE3;
    __m128                  jx_SSE, jy_SSE, jz_SSE, rbaj_SSE;
    __m128                  dx_SSE0, dy_SSE0, dz_SSE0;
    __m128                  dx_SSE1, dy_SSE1, dz_SSE1;
    __m128                  dx_SSE2, dy_SSE2, dz_SSE2;
    __m128                  dx_SSE3, dy_SSE3, dz_SSE3;
    __m128                  fgb_SSE0, fgb_ai_SSE0;
    __m128                  fgb_SSE1, fgb_ai_SSE1;
    __m128                  fgb_SSE2, fgb_ai_SSE2;
    __m128                  fgb_SSE3, fgb_ai_SSE3;
    __m128                  tx_SSE0, ty_SSE0, tz_SSE0;
    __m128                  tx_SSE1, ty_SSE1, tz_SSE1;
    __m128                  tx_SSE2, ty_SSE2, tz_SSE2;
    __m128                  tx_SSE3, ty_SSE3, tz_SSE3;
    __m128                  t1, t2;

    natoms              = mdatoms->nr;
    ni0                 = 0;
    ni1                 = mdatoms->homenr;
    dadx                = fr->dadx;

    aadata = (gmx_allvsallgb2_data_t *)paadata;

    x_align  = aadata->x_align;
    y_align  = aadata->y_align;
    z_align  = aadata->z_align;
    fx_align = aadata->fx_align;
    fy_align = aadata->fy_align;
    fz_align = aadata->fz_align;

    jindex    = aadata->jindex_gb;
    dadx      = fr->dadx;

    n  = 0;
    rb = aadata->work;

    /* Loop to get the proper form for the Born radius term */
    if (gb_algorithm == egbSTILL)
    {
        for (i = 0; i < natoms; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = (2 * rbi * rbi * fr->dvda[i])/ONE_4PI_EPS0;
        }
    }
    else if (gb_algorithm == egbHCT)
    {
        for (i = 0; i < natoms; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = rbi * rbi * fr->dvda[i];
        }
    }
    else if (gb_algorithm == egbOBC)
    {
        for (idx = 0; idx < natoms; idx++)
        {
            rbi     = born->bRad[idx];
            rb[idx] = rbi * rbi * born->drobc[idx] * fr->dvda[idx];
        }
    }

    for (i = 0; i < 2*natoms; i++)
    {
        fx_align[i]       = 0;
        fy_align[i]       = 0;
        fz_align[i]       = 0;
    }


    for (i = 0; i < natoms; i++)
    {
        rb[i+natoms] = rb[i];
    }

    for (i = ni0; i < ni1; i += UNROLLI)
    {
        /* We assume shifts are NOT used for all-vs-all interactions */

        /* Load i atom data */
        ix_SSE0          = _mm_load1_ps(x_align+i);
        iy_SSE0          = _mm_load1_ps(y_align+i);
        iz_SSE0          = _mm_load1_ps(z_align+i);
        ix_SSE1          = _mm_load1_ps(x_align+i+1);
        iy_SSE1          = _mm_load1_ps(y_align+i+1);
        iz_SSE1          = _mm_load1_ps(z_align+i+1);
        ix_SSE2          = _mm_load1_ps(x_align+i+2);
        iy_SSE2          = _mm_load1_ps(y_align+i+2);
        iz_SSE2          = _mm_load1_ps(z_align+i+2);
        ix_SSE3          = _mm_load1_ps(x_align+i+3);
        iy_SSE3          = _mm_load1_ps(y_align+i+3);
        iz_SSE3          = _mm_load1_ps(z_align+i+3);

        fix_SSE0         = _mm_setzero_ps();
        fiy_SSE0         = _mm_setzero_ps();
        fiz_SSE0         = _mm_setzero_ps();
        fix_SSE1         = _mm_setzero_ps();
        fiy_SSE1         = _mm_setzero_ps();
        fiz_SSE1         = _mm_setzero_ps();
        fix_SSE2         = _mm_setzero_ps();
        fiy_SSE2         = _mm_setzero_ps();
        fiz_SSE2         = _mm_setzero_ps();
        fix_SSE3         = _mm_setzero_ps();
        fiy_SSE3         = _mm_setzero_ps();
        fiz_SSE3         = _mm_setzero_ps();

        rbai_SSE0        = _mm_load1_ps(rb+i);
        rbai_SSE1        = _mm_load1_ps(rb+i+1);
        rbai_SSE2        = _mm_load1_ps(rb+i+2);
        rbai_SSE3        = _mm_load1_ps(rb+i+3);

        /* Load limits for loop over neighbors */
        nj0              = jindex[4*i];
        nj3              = jindex[4*i+3];

        /* No masks necessary, since the stored chain rule derivatives will be zero in those cases! */
        for (j = nj0; j < nj3; j += UNROLLJ)
        {
            /* load j atom coordinates */
            jx_SSE           = _mm_load_ps(x_align+j);
            jy_SSE           = _mm_load_ps(y_align+j);
            jz_SSE           = _mm_load_ps(z_align+j);

            /* Calculate distance */
            dx_SSE0          = _mm_sub_ps(ix_SSE0, jx_SSE);
            dy_SSE0          = _mm_sub_ps(iy_SSE0, jy_SSE);
            dz_SSE0          = _mm_sub_ps(iz_SSE0, jz_SSE);
            dx_SSE1          = _mm_sub_ps(ix_SSE1, jx_SSE);
            dy_SSE1          = _mm_sub_ps(iy_SSE1, jy_SSE);
            dz_SSE1          = _mm_sub_ps(iz_SSE1, jz_SSE);
            dx_SSE2          = _mm_sub_ps(ix_SSE2, jx_SSE);
            dy_SSE2          = _mm_sub_ps(iy_SSE2, jy_SSE);
            dz_SSE2          = _mm_sub_ps(iz_SSE2, jz_SSE);
            dx_SSE3          = _mm_sub_ps(ix_SSE3, jx_SSE);
            dy_SSE3          = _mm_sub_ps(iy_SSE3, jy_SSE);
            dz_SSE3          = _mm_sub_ps(iz_SSE3, jz_SSE);

            rbaj_SSE         = _mm_load_ps(rb+j);

            fgb_SSE0         = _mm_mul_ps(rbai_SSE0, _mm_load_ps(dadx));
            dadx            += 4;
            fgb_SSE1         = _mm_mul_ps(rbai_SSE1, _mm_load_ps(dadx));
            dadx            += 4;
            fgb_SSE2         = _mm_mul_ps(rbai_SSE2, _mm_load_ps(dadx));
            dadx            += 4;
            fgb_SSE3         = _mm_mul_ps(rbai_SSE3, _mm_load_ps(dadx));
            dadx            += 4;

            fgb_ai_SSE0      = _mm_mul_ps(rbaj_SSE, _mm_load_ps(dadx));
            dadx            += 4;
            fgb_ai_SSE1      = _mm_mul_ps(rbaj_SSE, _mm_load_ps(dadx));
            dadx            += 4;
            fgb_ai_SSE2      = _mm_mul_ps(rbaj_SSE, _mm_load_ps(dadx));
            dadx            += 4;
            fgb_ai_SSE3      = _mm_mul_ps(rbaj_SSE, _mm_load_ps(dadx));
            dadx            += 4;

            /* Total force between ai and aj is the sum of ai->aj and aj->ai */
            fgb_SSE0         = _mm_add_ps(fgb_SSE0, fgb_ai_SSE0);
            fgb_SSE1         = _mm_add_ps(fgb_SSE1, fgb_ai_SSE1);
            fgb_SSE2         = _mm_add_ps(fgb_SSE2, fgb_ai_SSE2);
            fgb_SSE3         = _mm_add_ps(fgb_SSE3, fgb_ai_SSE3);

            /* Calculate temporary vectorial force */
            tx_SSE0            = _mm_mul_ps(fgb_SSE0, dx_SSE0);
            ty_SSE0            = _mm_mul_ps(fgb_SSE0, dy_SSE0);
            tz_SSE0            = _mm_mul_ps(fgb_SSE0, dz_SSE0);
            tx_SSE1            = _mm_mul_ps(fgb_SSE1, dx_SSE1);
            ty_SSE1            = _mm_mul_ps(fgb_SSE1, dy_SSE1);
            tz_SSE1            = _mm_mul_ps(fgb_SSE1, dz_SSE1);
            tx_SSE2            = _mm_mul_ps(fgb_SSE2, dx_SSE2);
            ty_SSE2            = _mm_mul_ps(fgb_SSE2, dy_SSE2);
            tz_SSE2            = _mm_mul_ps(fgb_SSE2, dz_SSE2);
            tx_SSE3            = _mm_mul_ps(fgb_SSE3, dx_SSE3);
            ty_SSE3            = _mm_mul_ps(fgb_SSE3, dy_SSE3);
            tz_SSE3            = _mm_mul_ps(fgb_SSE3, dz_SSE3);

            /* Increment i atom force */
            fix_SSE0          = _mm_add_ps(fix_SSE0, tx_SSE0);
            fiy_SSE0          = _mm_add_ps(fiy_SSE0, ty_SSE0);
            fiz_SSE0          = _mm_add_ps(fiz_SSE0, tz_SSE0);
            fix_SSE1          = _mm_add_ps(fix_SSE1, tx_SSE1);
            fiy_SSE1          = _mm_add_ps(fiy_SSE1, ty_SSE1);
            fiz_SSE1          = _mm_add_ps(fiz_SSE1, tz_SSE1);
            fix_SSE2          = _mm_add_ps(fix_SSE2, tx_SSE2);
            fiy_SSE2          = _mm_add_ps(fiy_SSE2, ty_SSE2);
            fiz_SSE2          = _mm_add_ps(fiz_SSE2, tz_SSE2);
            fix_SSE3          = _mm_add_ps(fix_SSE3, tx_SSE3);
            fiy_SSE3          = _mm_add_ps(fiy_SSE3, ty_SSE3);
            fiz_SSE3          = _mm_add_ps(fiz_SSE3, tz_SSE3);

            /* Decrement j atom force */
            _mm_store_ps(fx_align+j,
                         _mm_sub_ps( _mm_load_ps(fx_align+j), gmx_mm_sum4_ps(tx_SSE0, tx_SSE1, tx_SSE2, tx_SSE3) ));
            _mm_store_ps(fy_align+j,
                         _mm_sub_ps( _mm_load_ps(fy_align+j), gmx_mm_sum4_ps(ty_SSE0, ty_SSE1, ty_SSE2, ty_SSE3) ));
            _mm_store_ps(fz_align+j,
                         _mm_sub_ps( _mm_load_ps(fz_align+j), gmx_mm_sum4_ps(tz_SSE0, tz_SSE1, tz_SSE2, tz_SSE3) ));
        }
        /* Add i forces to mem and shifted force list */
        _MM_TRANSPOSE4_PS(fix_SSE0, fix_SSE1, fix_SSE2, fix_SSE3);
        fix_SSE0 = _mm_add_ps(fix_SSE0, fix_SSE1);
        fix_SSE2 = _mm_add_ps(fix_SSE2, fix_SSE3);
        fix_SSE0 = _mm_add_ps(fix_SSE0, fix_SSE2);
        _mm_store_ps(fx_align+i, _mm_add_ps(fix_SSE0, _mm_load_ps(fx_align+i)));

        _MM_TRANSPOSE4_PS(fiy_SSE0, fiy_SSE1, fiy_SSE2, fiy_SSE3);
        fiy_SSE0 = _mm_add_ps(fiy_SSE0, fiy_SSE1);
        fiy_SSE2 = _mm_add_ps(fiy_SSE2, fiy_SSE3);
        fiy_SSE0 = _mm_add_ps(fiy_SSE0, fiy_SSE2);
        _mm_store_ps(fy_align+i, _mm_add_ps(fiy_SSE0, _mm_load_ps(fy_align+i)));

        _MM_TRANSPOSE4_PS(fiz_SSE0, fiz_SSE1, fiz_SSE2, fiz_SSE3);
        fiz_SSE0 = _mm_add_ps(fiz_SSE0, fiz_SSE1);
        fiz_SSE2 = _mm_add_ps(fiz_SSE2, fiz_SSE3);
        fiz_SSE0 = _mm_add_ps(fiz_SSE0, fiz_SSE2);
        _mm_store_ps(fz_align+i, _mm_add_ps(fiz_SSE0, _mm_load_ps(fz_align+i)));
    }

    for (i = 0; i < natoms; i++)
    {
        f[3*i]       += fx_align[i] + fx_align[natoms+i];
        f[3*i+1]     += fy_align[i] + fy_align[natoms+i];
        f[3*i+2]     += fz_align[i] + fz_align[natoms+i];
    }

    return 0;
}

#else
/* dummy variable when not using SSE */
int genborn_allvsall_sse2_single_dummy;


#endif
