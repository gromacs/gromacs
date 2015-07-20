/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/perf_est.h"

#include <math.h>

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"

/* Computational cost of bonded, non-bonded and PME calculations.
 * This will be machine dependent.
 * The numbers here are accurate for Intel Core2 and AMD Athlon 64
 * in single precision. In double precision PME mesh is slightly cheaper,
 * although not so much that the numbers need to be adjusted.
 */

/* Cost of a pair interaction in the "group" cut-off scheme */
#define C_GR_FQ        1.5
#define C_GR_QLJ_CUT   1.5
#define C_GR_QLJ_TAB   2.0
#define C_GR_LJ_CUT    1.0
#define C_GR_LJ_TAB    1.75
/* Cost of 1 water with one Q/LJ atom */
#define C_GR_QLJW_CUT  2.0
#define C_GR_QLJW_TAB  2.25
/* Cost of 1 water with one Q atom or with 1/3 water (LJ negligible) */
#define C_GR_QW        1.75

/* Cost of a pair interaction in the "Verlet" cut-off scheme, QEXP is Ewald */
#define C_VT_LJ        0.30
#define C_VT_QRF_LJ    0.40
#define C_VT_QRF       0.30
#define C_VT_QEXP_LJ   0.55
#define C_VT_QEXP      0.50
/* Extra cost for expensive LJ interaction, e.g. pot-switch or LJ-PME */
#define C_VT_LJEXP_ADD 0.20

/* Cost of PME, with all components running with SSE instructions */
/* Cost of particle reordering and redistribution */
#define C_PME_REDIST  12.0
/* Cost of q spreading and force interpolation per charge (mainly memory) */
#define C_PME_SPREAD  0.30
/* Cost of fft's, will be multiplied with N log(N) */
#define C_PME_FFT     0.20
/* Cost of pme_solve, will be multiplied with N */
#define C_PME_SOLVE   0.50

/* Cost of a bonded interaction divided by the number of (pbc_)dx nrequired */
#define C_BOND        5.0

int n_bonded_dx(gmx_mtop_t *mtop, gmx_bool bExcl)
{
    int            mb, nmol, ftype, ndxb, ndx_excl;
    int            ndx;
    gmx_moltype_t *molt;

    /* Count the number of pbc_rvec_sub calls required for bonded interactions.
     * This number is also roughly proportional to the computational cost.
     */
    ndx      = 0;
    ndx_excl = 0;
#if defined _ICC && __ICC == 1400 || defined __ICL && __ICL == 1400
#pragma novector /* Work-around for incorrect vectorization */
#endif
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        nmol = mtop->molblock[mb].nmol;
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (interaction_function[ftype].flags & IF_BOND)
            {
                switch (ftype)
                {
                    case F_POSRES:
                    case F_FBPOSRES:  ndxb = 1; break;
                    case F_CONNBONDS: ndxb = 0; break;
                    default:     ndxb      = NRAL(ftype) - 1; break;
                }
                ndx += nmol*ndxb*molt->ilist[ftype].nr/(1 + NRAL(ftype));
            }
        }
        if (bExcl)
        {
            ndx_excl += nmol*(molt->excls.nra - molt->atoms.nr)/2;
        }
        else
        {
            ndx_excl = 0;
        }
    }

    if (debug)
    {
        fprintf(debug, "ndx bonded %d exclusions %d\n", ndx, ndx_excl);
    }

    ndx += ndx_excl;

    return ndx;
}

static void pp_group_load(gmx_mtop_t *mtop, t_inputrec *ir, matrix box,
                          int *nq_tot, int *nlj_tot,
                          double *cost_pp,
                          gmx_bool *bChargePerturbed, gmx_bool *bTypePerturbed)
{
    t_atom        *atom;
    int            mb, nmol, atnr, cg, a, a0, ncqlj, ncq, nclj;
    gmx_bool       bBHAM, bLJcut, bWater, bQ, bLJ;
    int            nw, nqlj, nq, nlj;
    float          fq, fqlj, flj, fljtab, fqljw, fqw;
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
    fq    = C_GR_FQ;
    fqlj  = (bLJcut ? C_GR_QLJ_CUT : C_GR_QLJ_TAB);
    flj   = (bLJcut ? C_GR_LJ_CUT  : C_GR_LJ_TAB);
    /* Cost of 1 water with one Q/LJ atom */
    fqljw = (bLJcut ? C_GR_QLJW_CUT : C_GR_QLJW_TAB);
    /* Cost of 1 water with one Q atom or with 1/3 water (LJ negligible) */
    fqw   = C_GR_QW;

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
}

static void pp_verlet_load(gmx_mtop_t *mtop, t_inputrec *ir, matrix box,
                           int *nq_tot, int *nlj_tot,
                           double *cost_pp,
                           gmx_bool *bChargePerturbed, gmx_bool *bTypePerturbed)
{
    t_atom        *atom;
    int            mb, nmol, atnr, cg, a, a0, nqlj, nq, nlj;
    gmx_bool       bQRF;
    t_iparams     *iparams;
    gmx_moltype_t *molt;
    real           r_eff;
    double         c_qlj, c_q, c_lj;
    double         nat;
    /* Conversion factor for reference vs SIMD kernel performance.
     * The factor is about right for SSE2/4, but should be 2 higher for AVX256.
     */
#ifdef GMX_DOUBLE
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
        a    = 0;
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

    /* Effective cut-off for cluster pair list of 4x4 atoms */
    r_eff = ir->rlist + nbnxn_get_rlist_effective_inc(NBNXN_CPU_CLUSTER_I_SIZE, mtop->natoms/det(box));

    if (debug)
    {
        fprintf(debug, "nqlj %d nq %d nlj %d rlist %.3f r_eff %.3f\n",
                nqlj, nq, nlj, ir->rlist, r_eff);
    }

    /* Determine the cost per pair interaction */
    c_qlj = (bQRF ? C_VT_QRF_LJ : C_VT_QEXP_LJ);
    c_q   = (bQRF ? C_VT_QRF    : C_VT_QEXP);
    c_lj  = C_VT_LJ;
    if (ir->vdw_modifier == eintmodPOTSWITCH || EVDW_PME(ir->vdwtype))
    {
        c_qlj += C_VT_LJEXP_ADD;
        c_lj  += C_VT_LJEXP_ADD;
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
    /* Convert mtop->natoms to double to avoid int overflow */
    nat      = mtop->natoms;
    *cost_pp = 0.5*nat*(nqlj*c_qlj + nq*c_q + nlj*c_lj)
        *4/3*M_PI*r_eff*r_eff*r_eff/det(box);
}

float pme_load_estimate(gmx_mtop_t *mtop, t_inputrec *ir, matrix box)
{
    t_atom        *atom;
    int            mb, nmol, atnr, cg, a, a0, nq_tot, nlj_tot, f;
    gmx_bool       bBHAM, bLJcut, bChargePerturbed, bTypePerturbed;
    gmx_bool       bWater, bQ, bLJ;
    double         cost_bond, cost_pp, cost_redist, cost_spread, cost_fft, cost_solve, cost_pme;
    float          ratio;
    t_iparams     *iparams;
    gmx_moltype_t *molt;

    /* Computational cost of bonded, non-bonded and PME calculations.
     * This will be machine dependent.
     * The numbers here are accurate for Intel Core2 and AMD Athlon 64
     * in single precision. In double precision PME mesh is slightly cheaper,
     * although not so much that the numbers need to be adjusted.
     */

    iparams = mtop->ffparams.iparams;
    atnr    = mtop->ffparams.atnr;

    cost_bond = C_BOND*n_bonded_dx(mtop, TRUE);

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
        f            = ((ir->efep != efepNO && bChargePerturbed) ? 2 : 1);
        cost_redist +=   C_PME_REDIST*nq_tot;
        cost_spread += f*C_PME_SPREAD*nq_tot*pow(ir->pme_order, 3);
        cost_fft    += f*C_PME_FFT*ir->nkx*ir->nky*ir->nkz*log(ir->nkx*ir->nky*ir->nkz);
        cost_solve  += f*C_PME_SOLVE*ir->nkx*ir->nky*ir->nkz;
    }

    if (EVDW_PME(ir->vdwtype))
    {
        f            = ((ir->efep != efepNO && bTypePerturbed) ? 2 : 1);
        if (ir->ljpme_combination_rule == eljpmeLB)
        {
            /* LB combination rule: we have 7 mesh terms */
            f       *= 7;
        }
        cost_redist +=   C_PME_REDIST*nlj_tot;
        cost_spread += f*C_PME_SPREAD*nlj_tot*pow(ir->pme_order, 3);
        cost_fft    += f*C_PME_FFT*ir->nkx*ir->nky*ir->nkz*log(ir->nkx*ir->nky*ir->nkz);
        cost_solve  += f*C_PME_SOLVE*ir->nkx*ir->nky*ir->nkz;
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
