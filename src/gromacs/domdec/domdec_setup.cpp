/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,2011,2012,2013,2014,2015,2017, by the GROMACS development team, led by
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

/*! \internal \file
 *
 * \brief This file defines functions used by the domdec module
 * in its initial setup phase.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/perf_est.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/*! \brief Margin for setting up the DD grid */
#define DD_GRID_MARGIN_PRES_SCALE 1.05

/*! \brief Factorize \p n.
 *
 * \param[in]    n     Value to factorize
 * \param[out]   fac   Pointer to array of factors (to be allocated in this function)
 * \param[out]   mfac  Pointer to array of the number of times each factor repeats in the factorization (to be allocated in this function)
 * \return             The number of unique factors
 */
static int factorize(int n, int **fac, int **mfac)
{
    int d, ndiv;

    if (n <= 0)
    {
        gmx_fatal(FARGS, "Can only factorize positive integers.");
    }

    /* Decompose n in factors */
    snew(*fac, n/2);
    snew(*mfac, n/2);
    d    = 2;
    ndiv = 0;
    while (n > 1)
    {
        while (n % d == 0)
        {
            if (ndiv == 0 || (*fac)[ndiv-1] != d)
            {
                ndiv++;
                (*fac)[ndiv-1] = d;
            }
            (*mfac)[ndiv-1]++;
            n /= d;
        }
        d++;
    }

    return ndiv;
}

/*! \brief Find largest divisor of \p n smaller than \p n*/
static gmx_bool largest_divisor(int n)
{
    int ndiv, *div, *mdiv, ldiv;

    ndiv = factorize(n, &div, &mdiv);
    ldiv = div[ndiv-1];
    sfree(div);
    sfree(mdiv);

    return ldiv;
}

/*! \brief Compute largest common divisor of \p n1 and \b n2 */
static int lcd(int n1, int n2)
{
    int d, i;

    d = 1;
    for (i = 2; (i <= n1 && i <= n2); i++)
    {
        if (n1 % i == 0 && n2 % i == 0)
        {
            d = i;
        }
    }

    return d;
}

/*! \brief Returns TRUE when there are enough PME ranks for the ratio */
static gmx_bool fits_pme_ratio(int nrank_tot, int nrank_pme, float ratio)
{
    return ((double)nrank_pme/(double)nrank_tot > 0.95*ratio);
}

/*! \brief Returns TRUE when npme out of ntot ranks doing PME is expected to give reasonable performance */
static gmx_bool fits_pp_pme_perf(int ntot, int npme, float ratio)
{
    int ndiv, *div, *mdiv, ldiv;
    int npp_root3, npme_root2;

    ndiv = factorize(ntot - npme, &div, &mdiv);
    ldiv = div[ndiv-1];
    sfree(div);
    sfree(mdiv);

    npp_root3  = static_cast<int>(std::cbrt(ntot - npme) + 0.5);
    npme_root2 = static_cast<int>(std::sqrt(static_cast<double>(npme)) + 0.5);

    /* The check below gives a reasonable division:
     * factor 5 allowed at 5 or more PP ranks,
     * factor 7 allowed at 49 or more PP ranks.
     */
    if (ldiv > 3 + npp_root3)
    {
        return FALSE;
    }

    /* Check if the number of PP and PME ranks have a reasonable sized
     * denominator in common, such that we can use 2D PME decomposition
     * when required (which requires nx_pp == nx_pme).
     * The factor of 2 allows for a maximum ratio of 2^2=4
     * between nx_pme and ny_pme.
     */
    if (lcd(ntot - npme, npme)*2 < npme_root2)
    {
        return FALSE;
    }

    /* Does this division gives a reasonable PME load? */
    return fits_pme_ratio(ntot, npme, ratio);
}

/*! \brief Make a guess for the number of PME ranks to use. */
static int guess_npme(FILE *fplog, const gmx_mtop_t *mtop, const t_inputrec *ir,
                      const matrix box,
                      int nrank_tot)
{
    float      ratio;
    int        npme;

    ratio = pme_load_estimate(mtop, ir, box);

    if (fplog)
    {
        fprintf(fplog, "Guess for relative PME load: %.2f\n", ratio);
    }

    /* We assume the optimal rank ratio is close to the load ratio.
     * The communication load is neglected,
     * but (hopefully) this will balance out between PP and PME.
     */

    if (!fits_pme_ratio(nrank_tot, nrank_tot/2, ratio))
    {
        /* We would need more than nrank_tot/2 PME only nodes,
         * which is not possible. Since the PME load is very high,
         * we will not loose much performance when all ranks do PME.
         */

        return 0;
    }

    /* First try to find npme as a factor of nrank_tot up to nrank_tot/3.
     * We start with a minimum PME node fraction of 1/16
     * and avoid ratios which lead to large prime factors in nnodes-npme.
     */
    npme = (nrank_tot + 15)/16;
    while (npme <= nrank_tot/3)
    {
        if (nrank_tot % npme == 0)
        {
            /* Note that fits_perf might change the PME grid,
             * in the current implementation it does not.
             */
            if (fits_pp_pme_perf(nrank_tot, npme, ratio))
            {
                break;
            }
        }
        npme++;
    }
    if (npme > nrank_tot/3)
    {
        /* Try any possible number for npme */
        npme = 1;
        while (npme <= nrank_tot/2)
        {
            /* Note that fits_perf may change the PME grid */
            if (fits_pp_pme_perf(nrank_tot, npme, ratio))
            {
                break;
            }
            npme++;
        }
    }
    if (npme > nrank_tot/2)
    {
        gmx_fatal(FARGS, "Could not find an appropriate number of separate PME ranks. i.e. >= %5f*#ranks (%d) and <= #ranks/2 (%d) and reasonable performance wise (grid_x=%d, grid_y=%d).\n"
                  "Use the -npme option of mdrun or change the number of ranks or the PME grid dimensions, see the manual for details.",
                  ratio, (int)(0.95*ratio*nrank_tot + 0.5), nrank_tot/2, ir->nkx, ir->nky);
        /* Keep the compiler happy */
        npme = 0;
    }
    else
    {
        if (fplog)
        {
            fprintf(fplog,
                    "Will use %d particle-particle and %d PME only ranks\n"
                    "This is a guess, check the performance at the end of the log file\n",
                    nrank_tot - npme, npme);
        }
        fprintf(stderr, "\n"
                "Will use %d particle-particle and %d PME only ranks\n"
                "This is a guess, check the performance at the end of the log file\n",
                nrank_tot - npme, npme);
    }

    return npme;
}

/*! \brief Return \p n divided by \p f rounded up to the next integer. */
static int div_up(int n, int f)
{
    return (n + f - 1)/f;
}

real comm_box_frac(const ivec dd_nc, real cutoff, const gmx_ddbox_t *ddbox)
{
    int  i, j, k;
    rvec nw;
    real comm_vol;

    for (i = 0; i < DIM; i++)
    {
        real bt = ddbox->box_size[i]*ddbox->skew_fac[i];
        nw[i] = dd_nc[i]*cutoff/bt;
    }

    comm_vol = 0;
    for (i = 0; i < DIM; i++)
    {
        if (dd_nc[i] > 1)
        {
            comm_vol += nw[i];
            for (j = i+1; j < DIM; j++)
            {
                if (dd_nc[j] > 1)
                {
                    comm_vol += nw[i]*nw[j]*M_PI/4;
                    for (k = j+1; k < DIM; k++)
                    {
                        if (dd_nc[k] > 1)
                        {
                            comm_vol += nw[i]*nw[j]*nw[k]*M_PI/6;
                        }
                    }
                }
            }
        }
    }

    return comm_vol;
}

/*! \brief Return whether the DD inhomogeneous in the z direction */
static gmx_bool inhomogeneous_z(const t_inputrec *ir)
{
    return ((EEL_PME(ir->coulombtype) || ir->coulombtype == eelEWALD) &&
            ir->ePBC == epbcXYZ && ir->ewald_geometry == eewg3DC);
}

/*! \brief Estimate cost of PME FFT communication
 *
 * This only takes the communication into account and not imbalance
 * in the calculation. But the imbalance in communication and calculation
 * are similar and therefore these formulas also prefer load balance
 * in the FFT and pme_solve calculation.
 */
static float comm_pme_cost_vol(int npme, int a, int b, int c)
{
    /* We use a float here, since an integer might overflow */
    float comm_vol;

    comm_vol  = npme - 1;
    comm_vol *= npme;
    comm_vol *= div_up(a, npme);
    comm_vol *= div_up(b, npme);
    comm_vol *= c;

    return comm_vol;
}

/*! \brief Estimate cost of communication for a possible domain decomposition. */
static float comm_cost_est(real limit, real cutoff,
                           const matrix box, const gmx_ddbox_t *ddbox,
                           int natoms, const t_inputrec *ir,
                           float pbcdxr,
                           int npme_tot, ivec nc)
{
    ivec  npme = {1, 1, 1};
    int   i, j, nk, overlap;
    rvec  bt;
    float comm_vol, comm_vol_xf, comm_pme, cost_pbcdx;
    /* This is the cost of a pbc_dx call relative to the cost
     * of communicating the coordinate and force of an atom.
     * This will be machine dependent.
     * These factors are for x86 with SMP or Infiniband.
     */
    float pbcdx_rect_fac = 0.1;
    float pbcdx_tric_fac = 0.2;
    float temp;

    /* Check the DD algorithm restrictions */
    if ((ir->ePBC == epbcXY && ir->nwall < 2 && nc[ZZ] > 1) ||
        (ir->ePBC == epbcSCREW && (nc[XX] == 1 || nc[YY] > 1 || nc[ZZ] > 1)))
    {
        return -1;
    }

    if (inhomogeneous_z(ir) && nc[ZZ] > 1)
    {
        return -1;
    }

    assert(ddbox->npbcdim <= DIM);

    /* Check if the triclinic requirements are met */
    for (i = 0; i < DIM; i++)
    {
        for (j = i+1; j < ddbox->npbcdim; j++)
        {
            if (box[j][i] != 0 || ir->deform[j][i] != 0 ||
                (ir->epc != epcNO && ir->compress[j][i] != 0))
            {
                if (nc[j] > 1 && nc[i] == 1)
                {
                    return -1;
                }
            }
        }
    }

    for (i = 0; i < DIM; i++)
    {
        bt[i] = ddbox->box_size[i]*ddbox->skew_fac[i];

        /* Without PBC and with 2 cells, there are no lower limits on the cell size */
        if (!(i >= ddbox->npbcdim && nc[i] <= 2) && bt[i] < nc[i]*limit)
        {
            return -1;
        }
        /* With PBC, check if the cut-off fits in nc[i]-1 cells */
        if (i < ddbox->npbcdim && nc[i] > 1 && (nc[i] - 1)*bt[i] < nc[i]*cutoff)
        {
            return -1;
        }
    }

    if (npme_tot > 1)
    {
        /* The following choices should match those
         * in init_domain_decomposition in domdec.c.
         */
        if (nc[XX] == 1 && nc[YY] > 1)
        {
            npme[XX] = 1;
            npme[YY] = npme_tot;
        }
        else if (nc[YY] == 1)
        {
            npme[XX] = npme_tot;
            npme[YY] = 1;
        }
        else
        {
            /* Will we use 1D or 2D PME decomposition? */
            npme[XX] = (npme_tot % nc[XX] == 0) ? nc[XX] : npme_tot;
            npme[YY] = npme_tot/npme[XX];
        }
    }

    if (EEL_PME(ir->coulombtype) || EVDW_PME(ir->vdwtype))
    {
        /* Check the PME grid restrictions.
         * Currently these can only be invalid here with too few grid lines
         * along the x dimension per rank doing PME.
         */
        int npme_x = (npme_tot > 1 ? npme[XX] : nc[XX]);

        /* Currently we don't have the OpenMP thread count available here.
         * But with threads we have only tighter restrictions and it's
         * probably better anyhow to avoid settings where we need to reduce
         * grid lines over multiple ranks, as the thread check will do.
         */
        bool useThreads     = true;
        bool errorsAreFatal = false;
        if (!gmx_pme_check_restrictions(ir->pme_order, ir->nkx, ir->nky, ir->nkz,
                                        npme_x, useThreads, errorsAreFatal))
        {
            return -1;
        }
    }

    /* When two dimensions are (nearly) equal, use more cells
     * for the smallest index, so the decomposition does not
     * depend sensitively on the rounding of the box elements.
     */
    for (i = 0; i < DIM; i++)
    {
        for (j = i+1; j < DIM; j++)
        {
            /* Check if the box size is nearly identical,
             * in that case we prefer nx > ny  and ny > nz.
             */
            if (fabs(bt[j] - bt[i]) < 0.01*bt[i] && nc[j] > nc[i])
            {
                /* The XX/YY check is a bit compact. If nc[YY]==npme[YY]
                 * this means the swapped nc has nc[XX]==npme[XX],
                 * and we can also swap X and Y for PME.
                 */
                /* Check if dimension i and j are equivalent for PME.
                 * For x/y: if nc[YY]!=npme[YY], we can not swap x/y
                 * For y/z: we can not have PME decomposition in z
                 */
                if (npme_tot <= 1 ||
                    !((i == XX && j == YY && nc[YY] != npme[YY]) ||
                      (i == YY && j == ZZ && npme[YY] > 1)))
                {
                    return -1;
                }
            }
        }
    }

    /* This function determines only half of the communication cost.
     * All PP, PME and PP-PME communication is symmetric
     * and the "back"-communication cost is identical to the forward cost.
     */

    comm_vol = comm_box_frac(nc, cutoff, ddbox);

    comm_pme = 0;
    for (i = 0; i < 2; i++)
    {
        /* Determine the largest volume for PME x/f redistribution */
        if (nc[i] % npme[i] != 0)
        {
            if (nc[i] > npme[i])
            {
                comm_vol_xf = (npme[i] == 2 ? 1.0/3.0 : 0.5);
            }
            else
            {
                comm_vol_xf = 1.0 - lcd(nc[i], npme[i])/(double)npme[i];
            }
            comm_pme += 3*natoms*comm_vol_xf;
        }

        /* Grid overlap communication */
        if (npme[i] > 1)
        {
            nk        = (i == 0 ? ir->nkx : ir->nky);
            overlap   = (nk % npme[i] == 0 ? ir->pme_order-1 : ir->pme_order);
            temp      = npme[i];
            temp     *= overlap;
            temp     *= ir->nkx;
            temp     *= ir->nky;
            temp     *= ir->nkz;
            temp     /= nk;
            comm_pme += temp;
/* Old line comm_pme += npme[i]*overlap*ir->nkx*ir->nky*ir->nkz/nk; */
        }
    }

    comm_pme += comm_pme_cost_vol(npme[YY], ir->nky, ir->nkz, ir->nkx);
    comm_pme += comm_pme_cost_vol(npme[XX], ir->nkx, ir->nky, ir->nkz);

    /* Add cost of pbc_dx for bondeds */
    cost_pbcdx = 0;
    if ((nc[XX] == 1 || nc[YY] == 1) || (nc[ZZ] == 1 && ir->ePBC != epbcXY))
    {
        if ((ddbox->tric_dir[XX] && nc[XX] == 1) ||
            (ddbox->tric_dir[YY] && nc[YY] == 1))
        {
            cost_pbcdx = pbcdxr*pbcdx_tric_fac;
        }
        else
        {
            cost_pbcdx = pbcdxr*pbcdx_rect_fac;
        }
    }

    if (debug)
    {
        fprintf(debug,
                "nc %2d %2d %2d %2d %2d vol pp %6.4f pbcdx %6.4f pme %9.3e tot %9.3e\n",
                nc[XX], nc[YY], nc[ZZ], npme[XX], npme[YY],
                comm_vol, cost_pbcdx, comm_pme/(3*natoms),
                comm_vol + cost_pbcdx + comm_pme/(3*natoms));
    }

    return 3*natoms*(comm_vol + cost_pbcdx) + comm_pme;
}

/*! \brief Assign penalty factors to possible domain decompositions, based on the estimated communication costs. */
static void assign_factors(const gmx_domdec_t *dd,
                           real limit, real cutoff,
                           const matrix box, const gmx_ddbox_t *ddbox,
                           int natoms, const t_inputrec *ir,
                           float pbcdxr, int npme,
                           int ndiv, int *div, int *mdiv, ivec ir_try, ivec opt)
{
    int   x, y, i;
    float ce;

    if (ndiv == 0)
    {
        ce = comm_cost_est(limit, cutoff, box, ddbox,
                           natoms, ir, pbcdxr, npme, ir_try);
        if (ce >= 0 && (opt[XX] == 0 ||
                        ce < comm_cost_est(limit, cutoff, box, ddbox,
                                           natoms, ir, pbcdxr,
                                           npme, opt)))
        {
            copy_ivec(ir_try, opt);
        }

        return;
    }

    for (x = mdiv[0]; x >= 0; x--)
    {
        for (i = 0; i < x; i++)
        {
            ir_try[XX] *= div[0];
        }
        for (y = mdiv[0]-x; y >= 0; y--)
        {
            for (i = 0; i < y; i++)
            {
                ir_try[YY] *= div[0];
            }
            for (i = 0; i < mdiv[0]-x-y; i++)
            {
                ir_try[ZZ] *= div[0];
            }

            /* recurse */
            assign_factors(dd, limit, cutoff, box, ddbox, natoms, ir, pbcdxr, npme,
                           ndiv-1, div+1, mdiv+1, ir_try, opt);

            for (i = 0; i < mdiv[0]-x-y; i++)
            {
                ir_try[ZZ] /= div[0];
            }
            for (i = 0; i < y; i++)
            {
                ir_try[YY] /= div[0];
            }
        }
        for (i = 0; i < x; i++)
        {
            ir_try[XX] /= div[0];
        }
    }
}

/*! \brief Determine the optimal distribution of DD cells for the simulation system and number of MPI ranks */
static real optimize_ncells(FILE *fplog,
                            int nnodes_tot, int npme_only,
                            gmx_bool bDynLoadBal, real dlb_scale,
                            const gmx_mtop_t *mtop,
                            const matrix box, const gmx_ddbox_t *ddbox,
                            const t_inputrec *ir,
                            gmx_domdec_t *dd,
                            real cellsize_limit, real cutoff,
                            gmx_bool bInterCGBondeds,
                            ivec nc)
{
    int      npp, npme, ndiv, *div, *mdiv, d, nmax;
    double   pbcdxr;
    real     limit;
    ivec     itry;

    limit  = cellsize_limit;

    dd->nc[XX] = 1;
    dd->nc[YY] = 1;
    dd->nc[ZZ] = 1;

    npp = nnodes_tot - npme_only;
    if (EEL_PME(ir->coulombtype))
    {
        npme = (npme_only > 0 ? npme_only : npp);
    }
    else
    {
        npme = 0;
    }

    if (bInterCGBondeds)
    {
        /* If we can skip PBC for distance calculations in plain-C bondeds,
         * we can save some time (e.g. 3D DD with pbc=xyz).
         * Here we ignore SIMD bondeds as they always do (fast) PBC.
         */
        count_bonded_distances(mtop, ir, &pbcdxr, nullptr);
        pbcdxr /= (double)mtop->natoms;
    }
    else
    {
        /* Every molecule is a single charge group: no pbc required */
        pbcdxr = 0;
    }
    /* Add a margin for DLB and/or pressure scaling */
    if (bDynLoadBal)
    {
        if (dlb_scale >= 1.0)
        {
            gmx_fatal(FARGS, "The value for option -dds should be smaller than 1");
        }
        if (fplog)
        {
            fprintf(fplog, "Scaling the initial minimum size with 1/%g (option -dds) = %g\n", dlb_scale, 1/dlb_scale);
        }
        limit /= dlb_scale;
    }
    else if (ir->epc != epcNO)
    {
        if (fplog)
        {
            fprintf(fplog, "To account for pressure scaling, scaling the initial minimum size with %g\n", DD_GRID_MARGIN_PRES_SCALE);
            limit *= DD_GRID_MARGIN_PRES_SCALE;
        }
    }

    if (fplog)
    {
        fprintf(fplog, "Optimizing the DD grid for %d cells with a minimum initial size of %.3f nm\n", npp, limit);

        if (inhomogeneous_z(ir))
        {
            fprintf(fplog, "Ewald_geometry=%s: assuming inhomogeneous particle distribution in z, will not decompose in z.\n", eewg_names[ir->ewald_geometry]);
        }

        if (limit > 0)
        {
            fprintf(fplog, "The maximum allowed number of cells is:");
            for (d = 0; d < DIM; d++)
            {
                nmax = (int)(ddbox->box_size[d]*ddbox->skew_fac[d]/limit);
                if (d >= ddbox->npbcdim && nmax < 2)
                {
                    nmax = 2;
                }
                if (d == ZZ && inhomogeneous_z(ir))
                {
                    nmax = 1;
                }
                fprintf(fplog, " %c %d", 'X' + d, nmax);
            }
            fprintf(fplog, "\n");
        }
    }

    if (debug)
    {
        fprintf(debug, "Average nr of pbc_dx calls per atom %.2f\n", pbcdxr);
    }

    /* Decompose npp in factors */
    ndiv = factorize(npp, &div, &mdiv);

    itry[XX] = 1;
    itry[YY] = 1;
    itry[ZZ] = 1;
    clear_ivec(nc);
    assign_factors(dd, limit, cutoff, box, ddbox, mtop->natoms, ir, pbcdxr,
                   npme, ndiv, div, mdiv, itry, nc);

    sfree(div);
    sfree(mdiv);

    return limit;
}

real dd_choose_grid(FILE *fplog,
                    t_commrec *cr, gmx_domdec_t *dd,
                    const t_inputrec *ir,
                    const gmx_mtop_t *mtop,
                    const matrix box, const gmx_ddbox_t *ddbox,
                    int nPmeRanks,
                    gmx_bool bDynLoadBal, real dlb_scale,
                    real cellsize_limit, real cutoff_dd,
                    gmx_bool bInterCGBondeds)
{
    gmx_int64_t     nnodes_div, ldiv;
    real            limit;

    if (MASTER(cr))
    {
        nnodes_div = cr->nnodes;
        if (EEL_PME(ir->coulombtype))
        {
            if (nPmeRanks > 0)
            {
                if (nPmeRanks >= cr->nnodes)
                {
                    gmx_fatal(FARGS,
                              "Cannot have %d separate PME ranks with just %d total ranks",
                              nPmeRanks, cr->nnodes);
                }

                /* If the user purposely selected the number of PME nodes,
                 * only check for large primes in the PP node count.
                 */
                nnodes_div -= nPmeRanks;
            }
        }
        else
        {
            cr->npmenodes = 0;
        }

        if (nnodes_div > 12)
        {
            ldiv = largest_divisor(nnodes_div);
            /* Check if the largest divisor is more than nnodes^2/3 */
            if (ldiv*ldiv*ldiv > nnodes_div*nnodes_div)
            {
                gmx_fatal(FARGS, "The number of ranks you selected (%d) contains a large prime factor %d. In most cases this will lead to bad performance. Choose a number with smaller prime factors or set the decomposition (option -dd) manually.",
                          nnodes_div, ldiv);
            }
        }

        if (EEL_PME(ir->coulombtype))
        {
            if (nPmeRanks < 0)
            {
                /* Use PME nodes when the number of nodes is more than 16 */
                if (cr->nnodes <= 18)
                {
                    cr->npmenodes = 0;
                    if (fplog)
                    {
                        fprintf(fplog, "Using %d separate PME ranks, as there are too few total\n ranks for efficient splitting\n", cr->npmenodes);
                    }
                }
                else
                {
                    cr->npmenodes = guess_npme(fplog, mtop, ir, box, cr->nnodes);
                    if (fplog)
                    {
                        fprintf(fplog, "Using %d separate PME ranks, as guessed by mdrun\n", cr->npmenodes);
                    }
                }
            }
            else
            {
                /* We checked above that nPmeRanks is a valid number */
                cr->npmenodes = nPmeRanks;
                if (fplog)
                {
                    fprintf(fplog, "Using %d separate PME ranks, per user request\n", cr->npmenodes);
                }
            }
        }

        limit = optimize_ncells(fplog, cr->nnodes, cr->npmenodes,
                                bDynLoadBal, dlb_scale,
                                mtop, box, ddbox, ir, dd,
                                cellsize_limit, cutoff_dd,
                                bInterCGBondeds,
                                dd->nc);
    }
    else
    {
        limit = 0;
    }
    /* Communicate the information set by the master to all nodes */
    gmx_bcast(sizeof(dd->nc), dd->nc, cr);
    if (EEL_PME(ir->coulombtype))
    {
        gmx_bcast(sizeof(cr->npmenodes), &cr->npmenodes, cr);
    }
    else
    {
        cr->npmenodes = 0;
    }

    return limit;
}
