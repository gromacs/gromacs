/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#include "listed-nonbonded.h"
#include "bonded-impl.h"

#include "gromacs/math/vec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"

static void
nb_listed_warning_rlimit(const rvec *x, int ai, int aj, int * global_atom_index, real r, real rlimit)
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

static real
nb_evaluate_single(real r2, real tabscale, real *vftab,
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
    ntab             = 12*ntab;
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

real
do_nonbonded_listed(int ftype, int nbonds,
                    const t_iatom iatoms[], const t_iparams iparams[],
                    const rvec x[], rvec f[], rvec fshift[],
                    const t_pbc *pbc, const t_graph *g,
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
                nb_listed_warning_rlimit(x, ai, aj, global_atom_index, sqrt(r2), fr->tab14.r);
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

            fscal            = nb_free_energy_evaluate_single(r2, fr->sc_r_power, fr->sc_alphacoul, fr->sc_alphavdw,
                                                              fr->tab14.scale, fr->tab14.data, qq, c6, c12, qqB, c6B, c12B,
                                                              LFC, LFV, DLF, lfac_coul, lfac_vdw, dlfac_coul, dlfac_vdw,
                                                              fr->sc_sigma6_def, fr->sc_sigma6_min, sigma2_def, sigma2_min, &velec, &vvdw, dvdl);
        }
        else
        {
            /* Evaluate tabulated interaction without free energy */
            fscal            = nb_evaluate_single(r2, fr->tab14.scale, fr->tab14.data, qq, c6, c12, &velec, &vvdw);
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
