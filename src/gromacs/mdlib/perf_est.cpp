/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "perf_est.h"

#include <cmath>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/simd/simd.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"

/* Computational cost of bonded, non-bonded and PME calculations.
 * This will be machine dependent.
 * The numbers are only used for estimating the relative cost of PME vs PP,
 * so only relative numbers matter.
 * The numbers here are accurate cycle counts for Haswell in single precision
 * compiled with gcc5.2. A correction factor for other architectures is given
 * by simd_cycle_factor().
 * In double precision PME mesh is slightly cheaper, although not so much
 * that the numbers need to be adjusted.
 */

/* Cost of a pair interaction in the "group" cut-off scheme */
static const double c_group_fq        = 18.0;
static const double c_group_qlj_cut   = 18.0;
static const double c_group_qlj_tab   = 24.0;
static const double c_group_lj_cut    = 12.0;
static const double c_group_lj_tab    = 21.0;
/* Cost of 1 water with one Q/LJ atom */
static const double c_group_qljw_cut  = 24.0;
static const double c_group_qljw_tab  = 27.0;
/* Cost of 1 water with one Q atom or with 1/3 water (LJ negligible) */
static const double c_group_qw        = 21.0;

/* Cost of a pair interaction in the "Verlet" cut-off scheme, QEXP is Ewald */
static const double c_nbnxn_lj        =  2.5;
static const double c_nbnxn_qrf_lj    =  2.9;
static const double c_nbnxn_qrf       =  2.4;
static const double c_nbnxn_qexp_lj   =  4.2;
static const double c_nbnxn_qexp      =  3.8;
/* Extra cost for expensive LJ interaction, e.g. pot-switch or LJ-PME */
static const double c_nbnxn_ljexp_add =  1.0;

/* Cost of the different components of PME. */
/* Cost of particle reordering and redistribution (no SIMD correction).
 * This will be zero without MPI and can be very high with load imbalance.
 * Thus we use an approximate value for medium parallelization.
 */
static const double c_pme_redist = 100.0;
/* Cost of q spreading and force interpolation per charge. This part almost
 * doesn't accelerate with SIMD, so we don't use SIMD correction.
 */
static const double c_pme_spread =   5.0;
/* Cost of fft's, will be multiplied with 2 N log2(N) (no SIMD correction)
 * Without MPI the number is 2-3, depending on grid factors and thread count.
 * We take the high limit to be on the safe side and account for some MPI
 * communication cost, which will dominate at high parallelization.
 */
static const double c_pme_fft    =   3.0;
/* Cost of pme_solve, will be multiplied with N */
static const double c_pme_solve  =   9.0;

/* Cost of a bonded interaction divided by the number of distances calculations
 * required in one interaction. The actual cost is nearly propotional to this.
 */
static const double c_bond       =  25.0;


#if GMX_SIMD_HAVE_REAL
static const gmx_bool bHaveSIMD = TRUE;
#else
static const gmx_bool bHaveSIMD = FALSE;
#endif

/* Gives a correction factor for the currently compiled SIMD implementations
 * versus the reference used for the coefficients above (8-wide SIMD with FMA).
 * bUseSIMD sets if we asking for plain-C (FALSE) or SIMD (TRUE) code.
 */
static double simd_cycle_factor(gmx_bool bUseSIMD)
{
    /* The (average) ratio of the time taken by plain-C force calculations
     * relative to SIMD versions, for the reference platform Haswell:
     * 8-wide SIMD with FMA, factor: sqrt(2*8)*1.25 = 5.
     * This factor is used for normalization in simd_cycle_factor().
     */
    const double simd_cycle_no_simd = 5.0;
    double       speedup;

#if GMX_SIMD_HAVE_REAL
    if (bUseSIMD)
    {
        /* We never get full speed-up of a factor GMX_SIMD_REAL_WIDTH.
         * The actual speed-up depends very much on gather+scatter overhead,
         * which is different for different bonded and non-bonded kernels.
         * As a rough, but actually not bad, approximation we use a sqrt
         * dependence on the width which gives a factor 4 for width=8.
         */
        speedup = std::sqrt(2.0*GMX_SIMD_REAL_WIDTH);
#if GMX_SIMD_HAVE_FMA
        /* FMA tends to give a bit more speedup */
        speedup *= 1.25;
#endif
    }
    else
    {
        speedup  = 1.0;
    }
#else
    if (bUseSIMD)
    {
        gmx_incons("gmx_cycle_factor() compiled without SIMD called with bUseSIMD=TRUE");
    }
    /* No SIMD, no speedup */
    speedup      = 1.0;
#endif

    /* Return speed compared to the reference (Haswell).
     * For x86 SIMD, the nbnxn kernels are relatively much slower on
     * Sandy/Ivy Bridge than Haswell, but that only leads to a too high
     * PME load estimate on SB/IB, which is erring on the safe side.
     */
    return simd_cycle_no_simd/speedup;
}

void count_bonded_distances(const gmx_mtop_t *mtop, const t_inputrec *ir,
                            double *ndistance_c, double *ndistance_simd)
{
    gmx_bool       bExcl;
    double         nonsimd_step_frac;
    int            mb, nmol, ftype;
    gmx_moltype_t *molt;
    double         ndtot_c, ndtot_simd;
#if GMX_SIMD_HAVE_REAL
    gmx_bool       bSimdBondeds = TRUE;
#else
    gmx_bool       bSimdBondeds = FALSE;
#endif

    bExcl = (ir->cutoff_scheme == ecutsGROUP && inputrecExclForces(ir)
             && !EEL_FULL(ir->coulombtype));

    if (bSimdBondeds)
    {
        /* We only have SIMD versions of these bondeds without energy and
         * without shift-forces, we take that into account here.
         */
        if (ir->nstcalcenergy > 0)
        {
            nonsimd_step_frac = 1.0/ir->nstcalcenergy;
        }
        else
        {
            nonsimd_step_frac = 0;
        }
        if (ir->epc != epcNO && 1.0/ir->nstpcouple > nonsimd_step_frac)
        {
            nonsimd_step_frac = 1.0/ir->nstpcouple;
        }
    }
    else
    {
        nonsimd_step_frac = 1;
    }

    /* Count the number of pbc_rvec_sub calls required for bonded interactions.
     * This number is also roughly proportional to the computational cost.
     */
    ndtot_c    = 0;
    ndtot_simd = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        nmol = mtop->molblock[mb].nmol;
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            int nbonds;

            if (interaction_function[ftype].flags & IF_BOND)
            {
                double nd_c, nd_simd;

                nd_c    = 0;
                nd_simd = 0;
                /* For all interactions, except for the three exceptions
                 * in the switch below, #distances = #atoms - 1.
                 */
                switch (ftype)
                {
                    case F_POSRES:
                    case F_FBPOSRES:
                        nd_c    = 1;
                        break;
                    case F_CONNBONDS:
                        break;
                    /* These bonded potentially use SIMD */
                    case F_ANGLES:
                    case F_PDIHS:
                    case F_RBDIHS:
                    case F_LJ14:
                        nd_c    =      nonsimd_step_frac *(NRAL(ftype) - 1);
                        nd_simd = (1 - nonsimd_step_frac)*(NRAL(ftype) - 1);
                        break;
                    default:
                        nd_c    = NRAL(ftype) - 1;
                        break;
                }
                nbonds      = nmol*molt->ilist[ftype].nr/(1 + NRAL(ftype));
                ndtot_c    += nbonds*nd_c;
                ndtot_simd += nbonds*nd_simd;
            }
        }
        if (bExcl)
        {
            ndtot_c += nmol*(molt->excls.nra - molt->atoms.nr)/2;
        }
    }

    if (debug)
    {
        fprintf(debug, "nr. of distance calculations in bondeds: C %.1f SIMD %.1f\n", ndtot_c, ndtot_simd);
    }

    if (ndistance_c    != nullptr)
    {
        *ndistance_c    = ndtot_c;
    }
    if (ndistance_simd != nullptr)
    {
        *ndistance_simd = ndtot_simd;
    }
}

static void pp_group_load(const gmx_mtop_t *mtop, const t_inputrec *ir,
                          const matrix box,
                          int *nq_tot, int *nlj_tot,
                          double *cost_pp,
                          gmx_bool *bChargePerturbed, gmx_bool *bTypePerturbed)
{
    t_atom        *atom;
    int            mb, nmol, atnr, cg, a, a0, ncqlj, ncq, nclj;
    gmx_bool       bBHAM, bLJcut, bWater, bQ, bLJ;
    int            nw, nqlj, nq, nlj;
    double         fq, fqlj, flj, fqljw, fqw;
    t_iparams     *iparams;
    gmx_moltype_t *molt;

    bBHAM = (mtop->ffparams.functype[0] == F_BHAM);

    bLJcut = ((ir->vdwtype == evdwCUT) && !bBHAM);

    /* Computational cost of bonded, non-bonded and PME calculations.
     * This will be machine dependent.
     * The numbers here are accurate for Intel Core2 and AMD Athlon 64
     * in single precision. In double precision PME mesh is slightly cheaper,
     * although not so much that the numbers need to be adjusted.
     */
    fq    = c_group_fq;
    fqlj  = (bLJcut ? c_group_qlj_cut : c_group_qlj_tab);
    flj   = (bLJcut ? c_group_lj_cut  : c_group_lj_tab);
    /* Cost of 1 water with one Q/LJ atom */
    fqljw = (bLJcut ? c_group_qljw_cut : c_group_qljw_tab);
    /* Cost of 1 water with one Q atom or with 1/3 water (LJ negligible) */
    fqw   = c_group_qw;

    iparams           = mtop->ffparams.iparams;
    atnr              = mtop->ffparams.atnr;
    nw                = 0;
    nqlj              = 0;
    nq                = 0;
    nlj               = 0;
    *bChargePerturbed = FALSE;
    *bTypePerturbed   = FALSE;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        atom = molt->atoms.atom;
        nmol = mtop->molblock[mb].nmol;
        a    = 0;
        for (cg = 0; cg < molt->cgs.nr; cg++)
        {
            bWater = !bBHAM;
            ncqlj  = 0;
            ncq    = 0;
            nclj   = 0;
            a0     = a;
            while (a < molt->cgs.index[cg+1])
            {
                bQ  = (atom[a].q != 0 || atom[a].qB != 0);
                bLJ = (iparams[(atnr+1)*atom[a].type].lj.c6  != 0 ||
                       iparams[(atnr+1)*atom[a].type].lj.c12 != 0);
                if (atom[a].q != atom[a].qB)
                {
                    *bChargePerturbed = TRUE;
                }
                if (atom[a].type != atom[a].typeB)
                {
                    *bTypePerturbed = TRUE;
                }
                /* This if this atom fits into water optimization */
                if (!((a == a0   &&  bQ &&  bLJ) ||
                      (a == a0+1 &&  bQ && !bLJ) ||
                      (a == a0+2 &&  bQ && !bLJ && atom[a].q == atom[a-1].q) ||
                      (a == a0+3 && !bQ &&  bLJ)))
                {
                    bWater = FALSE;
                }
                if (bQ && bLJ)
                {
                    ncqlj++;
                }
                else
                {
                    if (bQ)
                    {
                        ncq++;
                    }
                    if (bLJ)
                    {
                        nclj++;
                    }
                }
                a++;
            }
            if (bWater)
            {
                nw   += nmol;
            }
            else
            {
                nqlj += nmol*ncqlj;
                nq   += nmol*ncq;
                nlj  += nmol*nclj;
            }
        }
    }

    *nq_tot  = nq  + nqlj + nw*3;
    *nlj_tot = nlj + nqlj + nw;

    if (debug)
    {
        fprintf(debug, "nw %d nqlj %d nq %d nlj %d\n", nw, nqlj, nq, nlj);
    }

    /* For the PP non-bonded cost it is (unrealistically) assumed
     * that all atoms are distributed homogeneously in space.
     * Factor 3 is used because a water molecule has 3 atoms
     * (and TIP4P effectively has 3 interactions with (water) atoms)).
     */
    *cost_pp = 0.5*(fqljw*nw*nqlj +
                    fqw  *nw*(3*nw + nq) +
                    fqlj *nqlj*nqlj +
                    fq   *nq*(3*nw + nqlj + nq) +
                    flj  *nlj*(nw + nqlj + nlj))
        *4/3*M_PI*ir->rlist*ir->rlist*ir->rlist/det(box);

    *cost_pp *= simd_cycle_factor(bHaveSIMD);
}

static void pp_verlet_load(const gmx_mtop_t *mtop, const t_inputrec *ir,
                           const matrix box,
                           int *nq_tot, int *nlj_tot,
                           double *cost_pp,
                           gmx_bool *bChargePerturbed, gmx_bool *bTypePerturbed)
{
    t_atom        *atom;
    int            mb, nmol, atnr, a, nqlj, nq, nlj;
    gmx_bool       bQRF;
    t_iparams     *iparams;
    gmx_moltype_t *molt;
    real           r_eff;
    double         c_qlj, c_q, c_lj;
    double         nppa;
    int            j_cluster_size;
    /* Conversion factor for reference vs SIMD kernel performance.
     * The factor is about right for SSE2/4, but should be 2 higher for AVX256.
     */
#if GMX_DOUBLE
    const real     nbnxn_refkernel_fac = 4.0;
#else
    const real     nbnxn_refkernel_fac = 8.0;
#endif

    bQRF = (EEL_RF(ir->coulombtype) || ir->coulombtype == eelCUT);

    iparams           = mtop->ffparams.iparams;
    atnr              = mtop->ffparams.atnr;
    nqlj              = 0;
    nq                = 0;
    *bChargePerturbed = FALSE;
    *bTypePerturbed   = FALSE;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        atom = molt->atoms.atom;
        nmol = mtop->molblock[mb].nmol;
        for (a = 0; a < molt->atoms.nr; a++)
        {
            if (atom[a].q != 0 || atom[a].qB != 0)
            {
                if (iparams[(atnr+1)*atom[a].type].lj.c6  != 0 ||
                    iparams[(atnr+1)*atom[a].type].lj.c12 != 0)
                {
                    nqlj += nmol;
                }
                else
                {
                    nq += nmol;
                }
            }
            if (atom[a].q != atom[a].qB)
            {
                *bChargePerturbed = TRUE;
            }
            if (atom[a].type != atom[a].typeB)
            {
                *bTypePerturbed = TRUE;
            }
        }
    }

    nlj = mtop->natoms - nqlj - nq;

    *nq_tot  = nqlj + nq;
    *nlj_tot = nqlj + nlj;

    /* Effective cut-off for cluster pair list of 4x4 or 4x8 atoms.
     * This choice should match the one of pick_nbnxn_kernel_cpu().
     * TODO: Make this function use pick_nbnxn_kernel_cpu().
     */
#if GMX_SIMD_HAVE_REAL && ((GMX_SIMD_REAL_WIDTH == 8 && defined GMX_SIMD_HAVE_FMA) || GMX_SIMD_REAL_WIDTH > 8)
    j_cluster_size = 8;
#else
    j_cluster_size = 4;
#endif
    r_eff = ir->rlist + nbnxn_get_rlist_effective_inc(j_cluster_size, mtop->natoms/det(box));

    /* The average number of pairs per atom */
    nppa  = 0.5*4/3*M_PI*r_eff*r_eff*r_eff*mtop->natoms/det(box);

    if (debug)
    {
        fprintf(debug, "nqlj %d nq %d nlj %d rlist %.3f r_eff %.3f pairs per atom %.1f\n",
                nqlj, nq, nlj, ir->rlist, r_eff, nppa);
    }

    /* Determine the cost per pair interaction */
    c_qlj = (bQRF ? c_nbnxn_qrf_lj : c_nbnxn_qexp_lj);
    c_q   = (bQRF ? c_nbnxn_qrf    : c_nbnxn_qexp);
    c_lj  = c_nbnxn_lj;
    if (ir->vdw_modifier == eintmodPOTSWITCH || EVDW_PME(ir->vdwtype))
    {
        c_qlj += c_nbnxn_ljexp_add;
        c_lj  += c_nbnxn_ljexp_add;
    }
    if (EVDW_PME(ir->vdwtype) && ir->ljpme_combination_rule == eljpmeLB)
    {
        /* We don't have LJ-PME LB comb. rule kernels, we use slow kernels */
        c_qlj *= nbnxn_refkernel_fac;
        c_q   *= nbnxn_refkernel_fac;
        c_lj  *= nbnxn_refkernel_fac;
    }

    /* For the PP non-bonded cost it is (unrealistically) assumed
     * that all atoms are distributed homogeneously in space.
     */
    *cost_pp = (nqlj*c_qlj + nq*c_q + nlj*c_lj)*nppa;

    *cost_pp *= simd_cycle_factor(bHaveSIMD);
}

float pme_load_estimate(const gmx_mtop_t *mtop, const t_inputrec *ir,
                        const matrix box)
{
    int            nq_tot, nlj_tot;
    gmx_bool       bChargePerturbed, bTypePerturbed;
    double         ndistance_c, ndistance_simd;
    double         cost_bond, cost_pp, cost_redist, cost_spread, cost_fft, cost_solve, cost_pme;
    float          ratio;

    /* Computational cost of bonded, non-bonded and PME calculations.
     * This will be machine dependent.
     * The numbers here are accurate for Intel Core2 and AMD Athlon 64
     * in single precision. In double precision PME mesh is slightly cheaper,
     * although not so much that the numbers need to be adjusted.
     */

    count_bonded_distances(mtop, ir, &ndistance_c, &ndistance_simd);
    /* C_BOND is the cost for bonded interactions with SIMD implementations,
     * so we need to scale the number of bonded interactions for which there
     * are only C implementations to the number of SIMD equivalents.
     */
    cost_bond = c_bond*(ndistance_c   *simd_cycle_factor(FALSE) +
                        ndistance_simd*simd_cycle_factor(bHaveSIMD));

    if (ir->cutoff_scheme == ecutsGROUP)
    {
        pp_group_load(mtop, ir, box,
                      &nq_tot, &nlj_tot, &cost_pp,
                      &bChargePerturbed, &bTypePerturbed);
    }
    else
    {
        pp_verlet_load(mtop, ir, box,
                       &nq_tot, &nlj_tot, &cost_pp,
                       &bChargePerturbed, &bTypePerturbed);
    }

    cost_redist = 0;
    cost_spread = 0;
    cost_fft    = 0;
    cost_solve  = 0;

    if (EEL_PME(ir->coulombtype))
    {
        double grid = ir->nkx*ir->nky*((ir->nkz + 1)/2);

        int    f     = ((ir->efep != efepNO && bChargePerturbed) ? 2 : 1);
        cost_redist +=   c_pme_redist*nq_tot;
        cost_spread += f*c_pme_spread*nq_tot*gmx::power3(ir->pme_order);
        cost_fft    += f*c_pme_fft*grid*std::log(grid)/std::log(2.0);
        cost_solve  += f*c_pme_solve*grid*simd_cycle_factor(bHaveSIMD);
    }

    if (EVDW_PME(ir->vdwtype))
    {
        double grid = ir->nkx*ir->nky*((ir->nkz + 1)/2);

        int    f     = ((ir->efep != efepNO && bTypePerturbed) ? 2 : 1);
        if (ir->ljpme_combination_rule == eljpmeLB)
        {
            /* LB combination rule: we have 7 mesh terms */
            f       *= 7;
        }
        cost_redist +=   c_pme_redist*nlj_tot;
        cost_spread += f*c_pme_spread*nlj_tot*gmx::power3(ir->pme_order);
        cost_fft    += f*c_pme_fft*2*grid*std::log(grid)/std::log(2.0);
        cost_solve  += f*c_pme_solve*grid*simd_cycle_factor(bHaveSIMD);
    }

    cost_pme = cost_redist + cost_spread + cost_fft + cost_solve;

    ratio = cost_pme/(cost_bond + cost_pp + cost_pme);

    if (debug)
    {
        fprintf(debug,
                "cost_bond   %f\n"
                "cost_pp     %f\n"
                "cost_redist %f\n"
                "cost_spread %f\n"
                "cost_fft    %f\n"
                "cost_solve  %f\n",
                cost_bond, cost_pp, cost_redist, cost_spread, cost_fft, cost_solve);

        fprintf(debug, "Estimate for relative PME load: %.3f\n", ratio);
    }

    return ratio;
}
