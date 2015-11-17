/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * \brief This file defines functions for "pair" interactions
 * (i.e. listed non-bonded interactions, e.g. 1-4 interactions)
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed-forces
 */
#include "gmxpre.h"

#include "pairs.h"

#include <cmath>

#include "gromacs/legacyheaders/types/group.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc-simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "listed-internal.h"

namespace
{

/*! \brief Issue a warning if a listed interaction is beyond a table limit */
void
warning_rlimit(const rvec *x, int ai, int aj, int * global_atom_index, real r, real rlimit)
{
    gmx_warning("Listed nonbonded interaction between particles %d and %d\n"
                "at distance %.3f which is larger than the table limit %.3f nm.\n\n"
                "This is likely either a 1,4 interaction, or a listed interaction inside\n"
                "a smaller molecule you are decoupling during a free energy calculation.\n"
                "Since interactions at distances beyond the table cannot be computed,\n"
                "they are skipped until they are inside the table limit again. You will\n"
                "only see this message once, even if it occurs for several interactions.\n\n"
                "IMPORTANT: This should not happen in a stable simulation, so there is\n"
                "probably something wrong with your system. Only change the table-extension\n"
                "distance in the mdp file if you are really sure that is the reason.\n",
                glatnr(global_atom_index, ai), glatnr(global_atom_index, aj), r, rlimit);

    if (debug)
    {
        fprintf(debug,
                "%8f %8f %8f\n%8f %8f %8f\n1-4 (%d,%d) interaction not within cut-off! r=%g. Ignored\n",
                x[ai][XX], x[ai][YY], x[ai][ZZ], x[aj][XX], x[aj][YY], x[aj][ZZ],
                glatnr(global_atom_index, ai), glatnr(global_atom_index, aj), r);
    }
}

/*! \brief Compute the energy and force for a single pair interaction */
real
evaluate_single(real r2, real tabscale, real *vftab, real tableStride,
                real qq, real c6, real c12, real *velec, real *vvdw)
{
    real       rinv, r, rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VVe, FFe, VVd, FFd, VVr, FFr, fscal;
    int        ntab;

    /* Do the tabulated interactions - first table lookup */
    rinv             = gmx_invsqrt(r2);
    r                = r2*rinv;
    rtab             = r*tabscale;
    ntab             = static_cast<int>(rtab);
    eps              = rtab-ntab;
    eps2             = eps*eps;
    ntab             = tableStride*ntab;
    /* Electrostatics */
    Y                = vftab[ntab];
    F                = vftab[ntab+1];
    Geps             = eps*vftab[ntab+2];
    Heps2            = eps2*vftab[ntab+3];
    Fp               = F+Geps+Heps2;
    VVe              = Y+eps*Fp;
    FFe              = Fp+Geps+2.0*Heps2;
    /* Dispersion */
    Y                = vftab[ntab+4];
    F                = vftab[ntab+5];
    Geps             = eps*vftab[ntab+6];
    Heps2            = eps2*vftab[ntab+7];
    Fp               = F+Geps+Heps2;
    VVd              = Y+eps*Fp;
    FFd              = Fp+Geps+2.0*Heps2;
    /* Repulsion */
    Y                = vftab[ntab+8];
    F                = vftab[ntab+9];
    Geps             = eps*vftab[ntab+10];
    Heps2            = eps2*vftab[ntab+11];
    Fp               = F+Geps+Heps2;
    VVr              = Y+eps*Fp;
    FFr              = Fp+Geps+2.0*Heps2;

    *velec           = qq*VVe;
    *vvdw            = c6*VVd+c12*VVr;

    fscal            = -(qq*FFe+c6*FFd+c12*FFr)*tabscale*rinv;

    return fscal;
}

/*! \brief Compute the energy and force for a single pair interaction under FEP */
real
free_energy_evaluate_single(real r2, real sc_r_power, real alpha_coul,
                            real alpha_vdw, real tabscale, real *vftab, real tableStride,
                            real qqA, real c6A, real c12A, real qqB,
                            real c6B, real c12B, real LFC[2], real LFV[2], real DLF[2],
                            real lfac_coul[2], real lfac_vdw[2], real dlfac_coul[2],
                            real dlfac_vdw[2], real sigma6_def, real sigma6_min,
                            real sigma2_def, real sigma2_min,
                            real *velectot, real *vvdwtot, real *dvdl)
{
    real       rp, rpm2, rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VV, FF, fscal;
    real       qq[2], c6[2], c12[2], sigma6[2], sigma2[2], sigma_pow[2];
    real       alpha_coul_eff, alpha_vdw_eff, dvdl_coul, dvdl_vdw;
    real       rpinv, r_coul, r_vdw, velecsum, vvdwsum;
    real       fscal_vdw[2], fscal_elec[2];
    real       velec[2], vvdw[2];
    int        i, ntab;
    const real half        = 0.5;
    const real minusOne    = -1.0;
    const real oneThird    = 1.0 / 3.0;
    const real one         = 1.0;
    const real two         = 2.0;
    const real six         = 6.0;
    const real fourtyeight = 48.0;

    qq[0]    = qqA;
    qq[1]    = qqB;
    c6[0]    = c6A;
    c6[1]    = c6B;
    c12[0]   = c12A;
    c12[1]   = c12B;

    if (sc_r_power == six)
    {
        rpm2             = r2*r2;   /* r4 */
        rp               = rpm2*r2; /* r6 */
    }
    else if (sc_r_power == fourtyeight)
    {
        rp               = r2*r2*r2; /* r6 */
        rp               = rp*rp;    /* r12 */
        rp               = rp*rp;    /* r24 */
        rp               = rp*rp;    /* r48 */
        rpm2             = rp/r2;    /* r46 */
    }
    else
    {
        rp             = std::pow(r2, half * sc_r_power);  /* not currently supported as input, but can handle it */
        rpm2           = rp/r2;
    }

    /* Loop over state A(0) and B(1) */
    for (i = 0; i < 2; i++)
    {
        if ((c6[i] > 0) && (c12[i] > 0))
        {
            /* The c6 & c12 coefficients now contain the constants 6.0 and 12.0, respectively.
             * Correct for this by multiplying with (1/12.0)/(1/6.0)=6.0/12.0=0.5.
             */
            sigma6[i]       = half*c12[i]/c6[i];
            sigma2[i]       = std::pow(half*c12[i]/c6[i], oneThird);
            /* should be able to get rid of this ^^^ internal pow call eventually.  Will require agreement on
               what data to store externally.  Can't be fixed without larger scale changes, so not 5.0 */
            if (sigma6[i] < sigma6_min)   /* for disappearing coul and vdw with soft core at the same time */
            {
                sigma6[i] = sigma6_min;
                sigma2[i] = sigma2_min;
            }
        }
        else
        {
            sigma6[i]       = sigma6_def;
            sigma2[i]       = sigma2_def;
        }
        if (sc_r_power == six)
        {
            sigma_pow[i]    = sigma6[i];
        }
        else if (sc_r_power == fourtyeight)
        {
            sigma_pow[i]    = sigma6[i]*sigma6[i];       /* sigma^12 */
            sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^24 */
            sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^48 */
        }
        else
        {       /* not really supported as input, but in here for testing the general case*/
            sigma_pow[i]    = std::pow(sigma2[i], sc_r_power/2);
        }
    }

    /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
    if ((c12[0] > 0) && (c12[1] > 0))
    {
        alpha_vdw_eff    = 0;
        alpha_coul_eff   = 0;
    }
    else
    {
        alpha_vdw_eff    = alpha_vdw;
        alpha_coul_eff   = alpha_coul;
    }

    /* Loop over A and B states again */
    for (i = 0; i < 2; i++)
    {
        fscal_elec[i] = 0;
        fscal_vdw[i]  = 0;
        velec[i]      = 0;
        vvdw[i]       = 0;

        /* Only spend time on A or B state if it is non-zero */
        if ( (qq[i] != 0) || (c6[i] != 0) || (c12[i] != 0) )
        {
            /* Coulomb */
            rpinv            = one/(alpha_coul_eff*lfac_coul[i]*sigma_pow[i]+rp);
            r_coul           = std::pow(rpinv, minusOne / sc_r_power);

            /* Electrostatics table lookup data */
            rtab             = r_coul*tabscale;
            ntab             = static_cast<int>(rtab);
            eps              = rtab-ntab;
            eps2             = eps*eps;
            ntab             = tableStride*ntab;
            /* Electrostatics */
            Y                = vftab[ntab];
            F                = vftab[ntab+1];
            Geps             = eps*vftab[ntab+2];
            Heps2            = eps2*vftab[ntab+3];
            Fp               = F+Geps+Heps2;
            VV               = Y+eps*Fp;
            FF               = Fp+Geps+two*Heps2;
            velec[i]         = qq[i]*VV;
            fscal_elec[i]    = -qq[i]*FF*r_coul*rpinv*tabscale;

            /* Vdw */
            rpinv            = one/(alpha_vdw_eff*lfac_vdw[i]*sigma_pow[i]+rp);
            r_vdw            = std::pow(rpinv, minusOne / sc_r_power);
            /* Vdw table lookup data */
            rtab             = r_vdw*tabscale;
            ntab             = static_cast<int>(rtab);
            eps              = rtab-ntab;
            eps2             = eps*eps;
            ntab             = 12*ntab;
            /* Dispersion */
            Y                = vftab[ntab+4];
            F                = vftab[ntab+5];
            Geps             = eps*vftab[ntab+6];
            Heps2            = eps2*vftab[ntab+7];
            Fp               = F+Geps+Heps2;
            VV               = Y+eps*Fp;
            FF               = Fp+Geps+two*Heps2;
            vvdw[i]          = c6[i]*VV;
            fscal_vdw[i]     = -c6[i]*FF;

            /* Repulsion */
            Y                = vftab[ntab+8];
            F                = vftab[ntab+9];
            Geps             = eps*vftab[ntab+10];
            Heps2            = eps2*vftab[ntab+11];
            Fp               = F+Geps+Heps2;
            VV               = Y+eps*Fp;
            FF               = Fp+Geps+two*Heps2;
            vvdw[i]         += c12[i]*VV;
            fscal_vdw[i]    -= c12[i]*FF;
            fscal_vdw[i]    *= r_vdw*rpinv*tabscale;
        }
    }
    /* Now we have velec[i], vvdw[i], and fscal[i] for both states */
    /* Assemble A and B states */
    velecsum  = 0;
    vvdwsum   = 0;
    dvdl_coul = 0;
    dvdl_vdw  = 0;
    fscal     = 0;
    for (i = 0; i < 2; i++)
    {
        velecsum      += LFC[i]*velec[i];
        vvdwsum       += LFV[i]*vvdw[i];

        fscal         += (LFC[i]*fscal_elec[i]+LFV[i]*fscal_vdw[i])*rpm2;

        dvdl_coul     += velec[i]*DLF[i] + LFC[i]*alpha_coul_eff*dlfac_coul[i]*fscal_elec[i]*sigma_pow[i];
        dvdl_vdw      += vvdw[i]*DLF[i] + LFV[i]*alpha_vdw_eff*dlfac_vdw[i]*fscal_vdw[i]*sigma_pow[i];
    }

    dvdl[efptCOUL]     += dvdl_coul;
    dvdl[efptVDW]      += dvdl_vdw;

    *velectot           = velecsum;
    *vvdwtot            = vvdwsum;

    return fscal;
}

} // namespace

real
do_pairs(int ftype, int nbonds,
         const t_iatom iatoms[], const t_iparams iparams[],
         const rvec x[], rvec f[], rvec fshift[],
         const struct t_pbc *pbc, const struct t_graph *g,
         real *lambda, real *dvdl,
         const t_mdatoms *md,
         const t_forcerec *fr, gmx_grppairener_t *grppener,
         int *global_atom_index)
{
    real             qq, c6, c12;
    rvec             dx;
    ivec             dt;
    int              i, itype, ai, aj, gid;
    int              fshift_index;
    real             r2;
    real             fscal, velec, vvdw;
    real *           energygrp_elec;
    real *           energygrp_vdw;
    static gmx_bool  warned_rlimit = FALSE;
    /* Free energy stuff */
    gmx_bool         bFreeEnergy;
    real             LFC[2], LFV[2], DLF[2], lfac_coul[2], lfac_vdw[2], dlfac_coul[2], dlfac_vdw[2];
    real             qqB, c6B, c12B, sigma2_def, sigma2_min;
    const real       oneThird = 1.0 / 3.0;

    switch (ftype)
    {
        case F_LJ14:
        case F_LJC14_Q:
            energygrp_elec = grppener->ener[egCOUL14];
            energygrp_vdw  = grppener->ener[egLJ14];
            break;
        case F_LJC_PAIRS_NB:
            energygrp_elec = grppener->ener[egCOULSR];
            energygrp_vdw  = grppener->ener[egLJSR];
            break;
        default:
            energygrp_elec = NULL; /* Keep compiler happy */
            energygrp_vdw  = NULL; /* Keep compiler happy */
            gmx_fatal(FARGS, "Unknown function type %d in do_nonbonded14", ftype);
            break;
    }

    if (fr->efep != efepNO)
    {
        /* Lambda factor for state A=1-lambda and B=lambda */
        LFC[0] = 1.0 - lambda[efptCOUL];
        LFV[0] = 1.0 - lambda[efptVDW];
        LFC[1] = lambda[efptCOUL];
        LFV[1] = lambda[efptVDW];

        /*derivative of the lambda factor for state A and B */
        DLF[0] = -1;
        DLF[1] = 1;

        /* precalculate */
        sigma2_def = pow(fr->sc_sigma6_def, oneThird);
        sigma2_min = pow(fr->sc_sigma6_min, oneThird);

        for (i = 0; i < 2; i++)
        {
            lfac_coul[i]  = (fr->sc_power == 2 ? (1-LFC[i])*(1-LFC[i]) : (1-LFC[i]));
            dlfac_coul[i] = DLF[i]*fr->sc_power/fr->sc_r_power*(fr->sc_power == 2 ? (1-LFC[i]) : 1);
            lfac_vdw[i]   = (fr->sc_power == 2 ? (1-LFV[i])*(1-LFV[i]) : (1-LFV[i]));
            dlfac_vdw[i]  = DLF[i]*fr->sc_power/fr->sc_r_power*(fr->sc_power == 2 ? (1-LFV[i]) : 1);
        }
    }
    else
    {
        sigma2_min = sigma2_def = 0;
    }

    /* TODO This code depends on the logic in tables.c that constructs
       the table layout, which should be made explicit in future
       cleanup. */
    GMX_ASSERT(etiNR == 3, "Pair-interaction code that uses GROMACS interaction tables supports exactly 3 tables");
    GMX_ASSERT(fr->tab14.interaction == GMX_TABLE_INTERACTION_ELEC_VDWREP_VDWDISP,
               "Pair interaction kernels need a table with Coulomb, repulsion and dispersion entries");

    bFreeEnergy = FALSE;
    for (i = 0; (i < nbonds); )
    {
        itype = iatoms[i++];
        ai    = iatoms[i++];
        aj    = iatoms[i++];
        gid   = GID(md->cENER[ai], md->cENER[aj], md->nenergrp);

        /* Get parameters */
        switch (ftype)
        {
            case F_LJ14:
                bFreeEnergy =
                    (fr->efep != efepNO &&
                     ((md->nPerturbed && (md->bPerturbed[ai] || md->bPerturbed[aj])) ||
                      iparams[itype].lj14.c6A != iparams[itype].lj14.c6B ||
                      iparams[itype].lj14.c12A != iparams[itype].lj14.c12B));
                qq               = md->chargeA[ai]*md->chargeA[aj]*fr->epsfac*fr->fudgeQQ;
                c6               = iparams[itype].lj14.c6A;
                c12              = iparams[itype].lj14.c12A;
                break;
            case F_LJC14_Q:
                qq               = iparams[itype].ljc14.qi*iparams[itype].ljc14.qj*fr->epsfac*iparams[itype].ljc14.fqq;
                c6               = iparams[itype].ljc14.c6;
                c12              = iparams[itype].ljc14.c12;
                break;
            case F_LJC_PAIRS_NB:
                qq               = iparams[itype].ljcnb.qi*iparams[itype].ljcnb.qj*fr->epsfac;
                c6               = iparams[itype].ljcnb.c6;
                c12              = iparams[itype].ljcnb.c12;
                break;
            default:
                /* Cannot happen since we called gmx_fatal() above in this case */
                qq = c6 = c12 = 0; /* Keep compiler happy */
                break;
        }

        /* To save flops in the optimized kernels, c6/c12 have 6.0/12.0 derivative prefactors
         * included in the general nfbp array now. This means the tables are scaled down by the
         * same factor, so when we use the original c6/c12 parameters from iparams[] they must
         * be scaled up.
         */
        c6  *= 6.0;
        c12 *= 12.0;

        /* Do we need to apply full periodic boundary conditions? */
        if (fr->bMolPBC == TRUE)
        {
            fshift_index = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            fshift_index = CENTRAL;
            rvec_sub(x[ai], x[aj], dx);
        }
        r2           = norm2(dx);

        if (r2 >= fr->tab14.r*fr->tab14.r)
        {
            /* This check isn't race free. But it doesn't matter because if a race occurs the only
             * disadvantage is that the warning is printed twice */
            if (warned_rlimit == FALSE)
            {
                warning_rlimit(x, ai, aj, global_atom_index, sqrt(r2), fr->tab14.r);
                warned_rlimit = TRUE;
            }
            continue;
        }

        if (bFreeEnergy)
        {
            /* Currently free energy is only supported for F_LJ14, so no need to check for that if we got here */
            qqB              = md->chargeB[ai]*md->chargeB[aj]*fr->epsfac*fr->fudgeQQ;
            c6B              = iparams[itype].lj14.c6B*6.0;
            c12B             = iparams[itype].lj14.c12B*12.0;

            fscal            = free_energy_evaluate_single(r2, fr->sc_r_power, fr->sc_alphacoul, fr->sc_alphavdw,
                                                           fr->tab14.scale, fr->tab14.data, fr->tab14.stride,
                                                           qq, c6, c12, qqB, c6B, c12B,
                                                           LFC, LFV, DLF, lfac_coul, lfac_vdw, dlfac_coul, dlfac_vdw,
                                                           fr->sc_sigma6_def, fr->sc_sigma6_min, sigma2_def, sigma2_min, &velec, &vvdw, dvdl);
        }
        else
        {
            /* Evaluate tabulated interaction without free energy */
            fscal            = evaluate_single(r2, fr->tab14.scale, fr->tab14.data, fr->tab14.stride, qq, c6, c12, &velec, &vvdw);
        }

        energygrp_elec[gid]  += velec;
        energygrp_vdw[gid]   += vvdw;
        svmul(fscal, dx, dx);

        /* Add the forces */
        rvec_inc(f[ai], dx);
        rvec_dec(f[aj], dx);

        if (g)
        {
            /* Correct the shift forces using the graph */
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            fshift_index = IVEC2IS(dt);
        }
        if (fshift_index != CENTRAL)
        {
            rvec_inc(fshift[fshift_index], dx);
            rvec_dec(fshift[CENTRAL], dx);
        }
    }
    return 0.0;
}

#if GMX_SIMD_HAVE_REAL

#if GMX_SIMD_X86_AVX_256 || GMX_SIMD_X86_AVX2_256

// This was originally work-in-progress for augmenting the SIMD module with
// masked load/store operations. Instead, that turned into and extended SIMD
// interface that supports gather/scatter in all platforms, which will be
// part of a future Gromacs version. However, since the code for bonded
// interactions and LINCS was already written it would be a pity not to get
// the performance gains in Gromacs-5.1. For this reason we have added it as
// a bit of a hack in the two files that use it. It will be replaced with the
// new generic functionality after version 5.1

#    ifdef GMX_DOUBLE
static gmx_inline void gmx_simdcall
gmx_hack_simd_transpose4_r(gmx_simd_double_t *row0,
                           gmx_simd_double_t *row1,
                           gmx_simd_double_t *row2,
                           gmx_simd_double_t *row3)
{
    __m256d tmp0, tmp1, tmp2, tmp3;

    tmp0  = _mm256_unpacklo_pd(*row0, *row1);
    tmp2  = _mm256_unpacklo_pd(*row2, *row3);
    tmp1  = _mm256_unpackhi_pd(*row0, *row1);
    tmp3  = _mm256_unpackhi_pd(*row2, *row3);
    *row0 = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
    *row1 = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
    *row2 = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);
    *row3 = _mm256_permute2f128_pd(tmp1, tmp3, 0x31);
}

static gmx_inline void gmx_simdcall
gmx_hack_simd4_transpose_to_simd_r(const gmx_simd4_double_t *a,
                                   gmx_simd_double_t        *row0,
                                   gmx_simd_double_t        *row1,
                                   gmx_simd_double_t        *row2,
                                   gmx_simd_double_t        *row3)
{
    *row0 = a[0];
    *row1 = a[1];
    *row2 = a[2];
    *row3 = a[3];

    gmx_hack_simd_transpose4_r(row0, row1, row2, row3);
}

#    if GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
#        define gmx_hack_simd4_load3_r(mem)    _mm256_maskload_pd((mem), _mm_castsi128_ps(_mm256_set_epi32(0, 0, -1, -1, -1, -1, -1, -1)))
#    else
#        define gmx_hack_simd4_load3_r(mem)    _mm256_maskload_pd((mem), _mm256_set_epi32(0, 0, -1, -1, -1, -1, -1, -1))
#    endif

#    else /* single instead of double */
static gmx_inline void gmx_simdcall
gmx_hack_simd_transpose4_r(gmx_simd_float_t *row0,
                           gmx_simd_float_t *row1,
                           gmx_simd_float_t *row2,
                           gmx_simd_float_t *row3)
{
    __m256 tmp0, tmp1, tmp2, tmp3;

    tmp0  = _mm256_unpacklo_ps(*row0, *row1);
    tmp2  = _mm256_unpacklo_ps(*row2, *row3);
    tmp1  = _mm256_unpackhi_ps(*row0, *row1);
    tmp3  = _mm256_unpackhi_ps(*row2, *row3);
    *row0 = _mm256_shuffle_ps(tmp0, tmp2, 0x44);
    *row1 = _mm256_shuffle_ps(tmp0, tmp2, 0xEE);
    *row2 = _mm256_shuffle_ps(tmp1, tmp3, 0x44);
    *row3 = _mm256_shuffle_ps(tmp1, tmp3, 0xEE);
}

static gmx_inline void gmx_simdcall
gmx_hack_simd4_transpose_to_simd_r(const gmx_simd4_float_t *a,
                                   gmx_simd_float_t        *row0,
                                   gmx_simd_float_t        *row1,
                                   gmx_simd_float_t        *row2,
                                   gmx_simd_float_t        *row3)
{
    *row0 = _mm256_insertf128_ps(_mm256_castps128_ps256(a[0]), a[4], 1);
    *row1 = _mm256_insertf128_ps(_mm256_castps128_ps256(a[1]), a[5], 1);
    *row2 = _mm256_insertf128_ps(_mm256_castps128_ps256(a[2]), a[6], 1);
    *row3 = _mm256_insertf128_ps(_mm256_castps128_ps256(a[3]), a[7], 1);

    gmx_hack_simd_transpose4_r(row0, row1, row2, row3);
}
#if GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
#        define gmx_hack_simd4_load3_r(mem)    _mm_maskload_ps((mem), _mm_castsi256_pd(_mm_set_epi32(0, -1, -1, -1)))
#else
#        define gmx_hack_simd4_load3_r(mem)    _mm_maskload_ps((mem), _mm_set_epi32(0, -1, -1, -1))
#endif

#endif

#endif /* AVX */


/*! \brief Store differences between indexed rvecs in SIMD registers.
 *
 * Returns SIMD register with the difference vectors:
 *     v[index0[i]] - v[index1[i]]
 *
 * \param[in]     v           Array of rvecs
 * \param[in]     index0      Index into the vector array
 * \param[in]     index1      Index into the vector array
 * \param[in,out] buf         Aligned tmp buffer of size 3*GMX_SIMD_REAL_WIDTH
 * \param[out]    dx          SIMD register with x difference
 * \param[out]    dy          SIMD register with y difference
 * \param[out]    dz          SIMD register with z difference
 */
static gmx_inline void gmx_simdcall
gmx_hack_simd_gather_rvec_dist_two_index(const rvec      *v,
                                         const int       *index0,
                                         const int       *index1,
                                         real gmx_unused *buf,
                                         gmx_simd_real_t *dx,
                                         gmx_simd_real_t *dy,
                                         gmx_simd_real_t *dz)
{
#if GMX_SIMD_X86_AVX_256 || GMX_SIMD_X86_AVX2_256
    int              i;
    gmx_simd4_real_t d[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  tmp;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        d[i] = gmx_simd4_sub_r(gmx_hack_simd4_load3_r(&(v[index0[i]][0])),
                               gmx_hack_simd4_load3_r(&(v[index1[i]][0])));

    }
    gmx_hack_simd4_transpose_to_simd_r(d, dx, dy, dz, &tmp);
#else /* generic SIMD */
#if GMX_ALIGNMENT
    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) buf_aligned[3*GMX_SIMD_REAL_WIDTH];
#else
    real* buf_aligned = buf;
#endif

    int i, m;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        /* Store the distances packed and aligned */
        for (m = 0; m < DIM; m++)
        {
            buf_aligned[m*GMX_SIMD_REAL_WIDTH + i] =
                v[index0[i]][m] - v[index1[i]][m];
        }
    }
    *dx = gmx_simd_load_r(buf_aligned + 0*GMX_SIMD_REAL_WIDTH);
    *dy = gmx_simd_load_r(buf_aligned + 1*GMX_SIMD_REAL_WIDTH);
    *dz = gmx_simd_load_r(buf_aligned + 2*GMX_SIMD_REAL_WIDTH);
#endif
}


/* As do_pairs(), but using SIMD to calculate many pairs at once.
 * This routine only supports ftype=F_LJ14.
 * This routine does not calculate energies and shift forces.
 * This routine does not support user tables or FEP.
 */
void
do_pairs_noener_simd(int nbonds,
                     const t_iatom iatoms[], const t_iparams iparams[],
                     const rvec x[], rvec f[],
                     const struct t_pbc *pbc,
                     const t_mdatoms *md,
                     const t_forcerec *fr)
{
    const int  nfa1 = 3;
    real       coeff_array[3*GMX_SIMD_REAL_WIDTH+GMX_SIMD_REAL_WIDTH], *coeff;
    real       dr_array[DIM*GMX_SIMD_REAL_WIDTH+GMX_SIMD_REAL_WIDTH], *dr;
    real       f_buf_array[3*GMX_SIMD_REAL_WIDTH+GMX_SIMD_REAL_WIDTH], *f_buf;
    pbc_simd_t pbc_simd;

    /* Ensure register memory alignment */
    coeff = gmx_simd_align_r(coeff_array);
    dr    = gmx_simd_align_r(dr_array);
    f_buf = gmx_simd_align_r(f_buf_array);

    set_pbc_simd(pbc, &pbc_simd);

    gmx_simd_real_t six    = gmx_simd_set1_r(6.0);
    gmx_simd_real_t twelve = gmx_simd_set1_r(12.0);
    gmx_simd_real_t ef     = gmx_simd_set1_r(fr->epsfac*fr->fudgeQQ);

    /* nbonds is #pairs*nfa1, here we step GMX_SIMD_REAL_WIDTH pairs */
    for (int i = 0; (i < nbonds); i += GMX_SIMD_REAL_WIDTH*nfa1)
    {
        /* Collect atoms for GMX_SIMD_REAL_WIDTH pairs.
         * iu indexes into iatoms, we should not let iu go beyond nbonds.
         */
        int ai[GMX_SIMD_REAL_WIDTH], aj[GMX_SIMD_REAL_WIDTH];
        int iu = i;
        for (int s = 0; s < GMX_SIMD_REAL_WIDTH; s++)
        {
            int itype = iatoms[iu];
            ai[s]     = iatoms[iu + 1];
            aj[s]     = iatoms[iu + 2];

            coeff[0*GMX_SIMD_REAL_WIDTH + s] = iparams[itype].lj14.c6A;
            coeff[1*GMX_SIMD_REAL_WIDTH + s] = iparams[itype].lj14.c12A;
            coeff[2*GMX_SIMD_REAL_WIDTH + s] = md->chargeA[ai[s]]*md->chargeA[aj[s]];

            /* At the end fill the arrays with identical entries */
            if (iu + nfa1 < nbonds)
            {
                iu += nfa1;
            }
        }

        /* Store the non PBC corrected distances packed and aligned */
        gmx_simd_real_t dx, dy, dz;
        gmx_hack_simd_gather_rvec_dist_two_index(x, ai, aj, dr, &dx, &dy, &dz);

        gmx_simd_real_t c6    = gmx_simd_load_r(coeff + 0*GMX_SIMD_REAL_WIDTH);
        gmx_simd_real_t c12   = gmx_simd_load_r(coeff + 1*GMX_SIMD_REAL_WIDTH);
        gmx_simd_real_t qq    = gmx_simd_load_r(coeff + 2*GMX_SIMD_REAL_WIDTH);

        /* We could save these operations by storing 6*C6,12*C12 */
        c6                    = gmx_simd_mul_r(six, c6);
        c12                   = gmx_simd_mul_r(twelve, c12);

        pbc_correct_dx_simd(&dx, &dy, &dz, &pbc_simd);

        gmx_simd_real_t rsq   = gmx_simd_calc_rsq_r(dx, dy, dz);
        gmx_simd_real_t rinv  = gmx_simd_invsqrt_r(rsq);
        gmx_simd_real_t rinv2 = gmx_simd_mul_r(rinv, rinv);
        gmx_simd_real_t rinv6 = gmx_simd_mul_r(rinv2, gmx_simd_mul_r(rinv2, rinv2));

        /* Calculate the Coulomb force * r */
        gmx_simd_real_t cfr   = gmx_simd_mul_r(gmx_simd_mul_r(ef, qq), rinv);

        /* Calculate the LJ force * r and add it to the Coulomb part */
        gmx_simd_real_t fr    = gmx_simd_fmadd_r(gmx_simd_fmsub_r(c12, rinv6, c6), rinv6, cfr);

        gmx_simd_real_t finvr = gmx_simd_mul_r(fr, rinv2);
        gmx_simd_real_t fx    = gmx_simd_mul_r(finvr, dx);
        gmx_simd_real_t fy    = gmx_simd_mul_r(finvr, dy);
        gmx_simd_real_t fz    = gmx_simd_mul_r(finvr, dz);

        gmx_simd_store_r(f_buf + 0*GMX_SIMD_REAL_WIDTH, fx);
        gmx_simd_store_r(f_buf + 1*GMX_SIMD_REAL_WIDTH, fy);
        gmx_simd_store_r(f_buf + 2*GMX_SIMD_REAL_WIDTH, fz);

        iu    = i;
        int s = 0;
        do
        {
            for (int m = 0; m < DIM; m++)
            {
                f[ai[s]][m] += f_buf[s + m*GMX_SIMD_REAL_WIDTH];
                f[aj[s]][m] -= f_buf[s + m*GMX_SIMD_REAL_WIDTH];
            }
            s++;
            iu += nfa1;
        }
        while (s < GMX_SIMD_REAL_WIDTH && iu < nbonds);
    }
}

#endif /* GMX_SIMD_HAVE_REAL */
