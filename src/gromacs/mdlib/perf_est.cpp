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

#include "perf_est.h"

#include <cmath>
#include <cstdio>

#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/nbnxm_geometry.h"
#include "gromacs/simd/simd.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/real.h"

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

/* Cost of a pair interaction in the "Verlet" cut-off scheme, QEXP is Ewald */
static const double c_nbnxn_lj      = 2.5;
static const double c_nbnxn_qrf_lj  = 2.9;
static const double c_nbnxn_qrf     = 2.4;
static const double c_nbnxn_qexp_lj = 4.2;
static const double c_nbnxn_qexp    = 3.8;
/* Extra cost for expensive LJ interaction, e.g. pot-switch or LJ-PME */
static const double c_nbnxn_ljexp_add = 1.0;

/* Cost of the different components of PME. */
/* Cost of particle reordering and redistribution (no SIMD correction).
 * This will be zero without MPI and can be very high with load imbalance.
 * Thus we use an approximate value for medium parallelization.
 */
static const double c_pme_redist = 100.0;
/* Cost of q spreading and force interpolation per charge. This part almost
 * doesn't accelerate with SIMD, so we don't use SIMD correction.
 */
static const double c_pme_spread = 5.0;
/* Cost of fft's, will be multiplied with 2 N log2(N) (no SIMD correction)
 * Without MPI the number is 2-3, depending on grid factors and thread count.
 * We take the high limit to be on the safe side and account for some MPI
 * communication cost, which will dominate at high parallelization.
 */
static const double c_pme_fft = 3.0;
/* Cost of pme_solve, will be multiplied with N */
static const double c_pme_solve = 9.0;

/* Cost of a bonded interaction divided by the number of distances calculations
 * required in one interaction. The actual cost is nearly propotional to this.
 */
static const double c_bond = 25.0;


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
        speedup = std::sqrt(2.0 * GMX_SIMD_REAL_WIDTH);
#    if GMX_SIMD_HAVE_FMA
        /* FMA tends to give a bit more speedup */
        speedup *= 1.25;
#    endif
    }
    else
    {
        speedup = 1.0;
    }
#else
    if (bUseSIMD)
    {
        gmx_incons("gmx_cycle_factor() compiled without SIMD called with bUseSIMD=TRUE");
    }
    /* No SIMD, no speedup */
    speedup                        = 1.0;
#endif

    /* Return speed compared to the reference (Haswell).
     * For x86 SIMD, the nbnxn kernels are relatively much slower on
     * Sandy/Ivy Bridge than Haswell, but that only leads to a too high
     * PME load estimate on SB/IB, which is erring on the safe side.
     */
    return simd_cycle_no_simd / speedup;
}

void count_bonded_distances(const gmx_mtop_t& mtop, const t_inputrec& ir, double* ndistance_c, double* ndistance_simd)
{
    gmx_bool bExcl;
    double   nonsimd_step_frac;
    int      ftype;
    double   ndtot_c, ndtot_simd;
#if GMX_SIMD_HAVE_REAL
    gmx_bool bSimdBondeds = TRUE;
#else
    gmx_bool   bSimdBondeds        = FALSE;
#endif

    bExcl = (ir.cutoff_scheme == CutoffScheme::Group && inputrecExclForces(&ir)
             && !usingFullElectrostatics(ir.coulombtype));

    if (bSimdBondeds)
    {
        /* We only have SIMD versions of these bondeds without energy and
         * without shift-forces, we take that into account here.
         */
        if (ir.nstcalcenergy > 0)
        {
            nonsimd_step_frac = 1.0 / ir.nstcalcenergy;
        }
        else
        {
            nonsimd_step_frac = 0;
        }
        if (ir.pressureCouplingOptions.epc != PressureCoupling::No
            && 1.0 / ir.pressureCouplingOptions.nstpcouple > nonsimd_step_frac)
        {
            nonsimd_step_frac = 1.0 / ir.pressureCouplingOptions.nstpcouple;
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
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t* molt = &mtop.moltype[molb.type];
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
                    case F_FBPOSRES: nd_c = 1; break;
                    case F_CONNBONDS: break;
                    /* These bonded potentially use SIMD */
                    case F_ANGLES:
                    case F_PDIHS:
                    case F_RBDIHS:
                    case F_LJ14:
                        nd_c    = nonsimd_step_frac * (NRAL(ftype) - 1);
                        nd_simd = (1 - nonsimd_step_frac) * (NRAL(ftype) - 1);
                        break;
                    default: nd_c = NRAL(ftype) - 1; break;
                }
                nbonds = molb.nmol * molt->ilist[ftype].size() / (1 + NRAL(ftype));
                ndtot_c += nbonds * nd_c;
                ndtot_simd += nbonds * nd_simd;
            }
        }
        if (bExcl)
        {
            ndtot_c += molb.nmol * (molt->excls.numElements() - molt->atoms.nr) / 2.;
        }
    }

    if (debug)
    {
        fprintf(debug, "nr. of distance calculations in bondeds: C %.1f SIMD %.1f\n", ndtot_c, ndtot_simd);
    }

    if (ndistance_c != nullptr)
    {
        *ndistance_c = ndtot_c;
    }
    if (ndistance_simd != nullptr)
    {
        *ndistance_simd = ndtot_simd;
    }
}

static void pp_verlet_load(const gmx_mtop_t& mtop,
                           const t_inputrec& ir,
                           const matrix      box,
                           int*              nq_tot,
                           int*              nlj_tot,
                           double*           cost_pp,
                           gmx_bool*         bChargePerturbed,
                           gmx_bool*         bTypePerturbed)
{
    int      atnr, a, nqlj, nq, nlj;
    gmx_bool bQRF;
    real     r_eff;
    double   c_qlj, c_q, c_lj;
    double   nppa;
    /* Conversion factor for reference vs SIMD kernel performance.
     * The factor is about right for SSE2/4, but should be 2 higher for AVX256.
     */
#if GMX_DOUBLE
    const real nbnxn_refkernel_fac = 4.0;
#else
    const real nbnxn_refkernel_fac = 8.0;
#endif

    bQRF = (usingRF(ir.coulombtype) || ir.coulombtype == CoulombInteractionType::Cut);

    gmx::ArrayRef<const t_iparams> iparams = mtop.ffparams.iparams;
    atnr                                   = mtop.ffparams.atnr;
    nqlj                                   = 0;
    nq                                     = 0;
    *bChargePerturbed                      = FALSE;
    *bTypePerturbed                        = FALSE;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t* molt = &mtop.moltype[molb.type];
        const t_atom*        atom = molt->atoms.atom;
        for (a = 0; a < molt->atoms.nr; a++)
        {
            if (atom[a].q != 0 || atom[a].qB != 0)
            {
                if (iparams[(atnr + 1) * atom[a].type].lj.c6 != 0
                    || iparams[(atnr + 1) * atom[a].type].lj.c12 != 0)
                {
                    nqlj += molb.nmol;
                }
                else
                {
                    nq += molb.nmol;
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

    nlj = mtop.natoms - nqlj - nq;

    *nq_tot  = nqlj + nq;
    *nlj_tot = nqlj + nlj;

    /* Effective radius of a CPU pairlist including the pairs beyond rlist */
    r_eff = ir.rlist + gmx::nbnxmPairlistVolumeRadiusIncrease(false, mtop.natoms / det(box));

    /* The average number of pairs per atom */
    nppa = 0.5 * 4 / 3 * M_PI * r_eff * r_eff * r_eff * mtop.natoms / det(box);

    if (debug)
    {
        fprintf(debug,
                "nqlj %d nq %d nlj %d rlist %.3f r_eff %.3f pairs per atom %.1f\n",
                nqlj,
                nq,
                nlj,
                ir.rlist,
                r_eff,
                nppa);
    }

    /* Determine the cost per pair interaction */
    c_qlj = (bQRF ? c_nbnxn_qrf_lj : c_nbnxn_qexp_lj);
    c_q   = (bQRF ? c_nbnxn_qrf : c_nbnxn_qexp);
    c_lj  = c_nbnxn_lj;
    if (ir.vdw_modifier == InteractionModifiers::PotSwitch || usingLJPme(ir.vdwtype))
    {
        c_qlj += c_nbnxn_ljexp_add;
        c_lj += c_nbnxn_ljexp_add;
    }
    if (usingLJPme(ir.vdwtype) && ir.ljpme_combination_rule == LongRangeVdW::LB)
    {
        /* We don't have LJ-PME LB comb. rule kernels, we use slow kernels */
        c_qlj *= nbnxn_refkernel_fac;
        c_q *= nbnxn_refkernel_fac;
        c_lj *= nbnxn_refkernel_fac;
    }

    /* For the PP non-bonded cost it is (unrealistically) assumed
     * that all atoms are distributed homogeneously in space.
     */
    *cost_pp = (nqlj * c_qlj + nq * c_q + nlj * c_lj) * nppa;

    *cost_pp *= simd_cycle_factor(bHaveSIMD);
}

float pme_load_estimate(const gmx_mtop_t& mtop, const t_inputrec& ir, const matrix box)
{
    int      nq_tot, nlj_tot;
    gmx_bool bChargePerturbed, bTypePerturbed;
    double   ndistance_c, ndistance_simd;
    double   cost_bond, cost_pp, cost_redist, cost_spread, cost_fft, cost_solve, cost_pme;
    float    ratio;

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
    cost_bond = c_bond
                * (ndistance_c * simd_cycle_factor(FALSE) + ndistance_simd * simd_cycle_factor(bHaveSIMD));

    pp_verlet_load(mtop, ir, box, &nq_tot, &nlj_tot, &cost_pp, &bChargePerturbed, &bTypePerturbed);

    cost_redist = 0;
    cost_spread = 0;
    cost_fft    = 0;
    cost_solve  = 0;

    int gridNkzFactor = int{ (ir.nkz + 1) / 2 };
    if (usingPme(ir.coulombtype))
    {
        double grid = ir.nkx * ir.nky * gridNkzFactor;

        int f = ((ir.efep != FreeEnergyPerturbationType::No && bChargePerturbed) ? 2 : 1);
        cost_redist += c_pme_redist * nq_tot;
        cost_spread += f * c_pme_spread * nq_tot * gmx::power3(ir.pme_order);
        cost_fft += f * c_pme_fft * grid * std::log(grid) / std::log(2.0);
        cost_solve += f * c_pme_solve * grid * simd_cycle_factor(bHaveSIMD);
    }

    if (usingLJPme(ir.vdwtype))
    {
        double grid = ir.nkx * ir.nky * gridNkzFactor;

        int f = ((ir.efep != FreeEnergyPerturbationType::No && bTypePerturbed) ? 2 : 1);
        if (ir.ljpme_combination_rule == LongRangeVdW::LB)
        {
            /* LB combination rule: we have 7 mesh terms */
            f *= 7;
        }
        cost_redist += c_pme_redist * nlj_tot;
        cost_spread += f * c_pme_spread * nlj_tot * gmx::power3(ir.pme_order);
        cost_fft += f * c_pme_fft * 2 * grid * std::log(grid) / std::log(2.0);
        cost_solve += f * c_pme_solve * grid * simd_cycle_factor(bHaveSIMD);
    }

    cost_pme = cost_redist + cost_spread + cost_fft + cost_solve;

    ratio = cost_pme / (cost_bond + cost_pp + cost_pme);

    if (debug)
    {
        fprintf(debug,
                "cost_bond   %f\n"
                "cost_pp     %f\n"
                "cost_redist %f\n"
                "cost_spread %f\n"
                "cost_fft    %f\n"
                "cost_solve  %f\n",
                cost_bond,
                cost_pp,
                cost_redist,
                cost_spread,
                cost_fft,
                cost_solve);

        fprintf(debug, "Estimate for relative PME load: %.3f\n", ratio);
    }

    return ratio;
}
