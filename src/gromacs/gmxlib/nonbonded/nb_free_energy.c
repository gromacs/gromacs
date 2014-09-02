/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "vec.h"
#include "typedefs.h"
#include "nonbonded.h"
#include "nb_kernel.h"
#include "nrnb.h"
#include "macros.h"
#include "nb_free_energy.h"

#include "gmx_fatal.h"

void
gmx_nb_free_energy_kernel(const t_nblist * gmx_restrict    nlist,
                          rvec * gmx_restrict              xx,
                          rvec * gmx_restrict              ff,
                          t_forcerec * gmx_restrict        fr,
                          const t_mdatoms * gmx_restrict   mdatoms,
                          nb_kernel_data_t * gmx_restrict  kernel_data,
                          t_nrnb * gmx_restrict            nrnb)
{

#define  STATE_A  0
#define  STATE_B  1
#define  NSTATES  2
    int           i, j, n, ii, is3, ii3, k, nj0, nj1, jnr, j3, ggid;
    real          shX, shY, shZ;
    real          tx, ty, tz, Fscal;
    double        FscalC[NSTATES], FscalV[NSTATES];  /* Needs double for sc_power==48 */
    double        Vcoul[NSTATES], Vvdw[NSTATES];     /* Needs double for sc_power==48 */
    real          rinv6, r, rt, rtC, rtV;
    real          iqA, iqB;
    real          qq[NSTATES], vctot, krsq;
    int           ntiA, ntiB, tj[NSTATES];
    real          Vvdw6, Vvdw12, vvtot;
    real          ix, iy, iz, fix, fiy, fiz;
    real          dx, dy, dz, rsq, rinv;
    real          c6[NSTATES], c12[NSTATES], c6grid;
    real          LFC[NSTATES], LFV[NSTATES], DLF[NSTATES];
    double        dvdl_coul, dvdl_vdw;
    real          lfac_coul[NSTATES], dlfac_coul[NSTATES], lfac_vdw[NSTATES], dlfac_vdw[NSTATES];
    real          sigma6[NSTATES], alpha_vdw_eff, alpha_coul_eff, sigma2_def, sigma2_min;
    double        rp, rpm2, rC, rV, rinvC, rpinvC, rinvV, rpinvV; /* Needs double for sc_power==48 */
    real          sigma2[NSTATES], sigma_pow[NSTATES], sigma_powm2[NSTATES], rs, rs2;
    int           do_tab, tab_elemsize;
    int           n0, n1C, n1V, nnn;
    real          Y, F, G, H, Fp, Geps, Heps2, epsC, eps2C, epsV, eps2V, VV, FF;
    int           icoul, ivdw;
    int           nri;
    const int *   iinr;
    const int *   jindex;
    const int *   jjnr;
    const int *   shift;
    const int *   gid;
    const int *   typeA;
    const int *   typeB;
    int           ntype;
    const real *  shiftvec;
    real          dvdl_part;
    real *        fshift;
    real          tabscale = 0;
    const real *  VFtab    = NULL;
    const real *  x;
    real *        f;
    real          facel, krf, crf;
    const real *  chargeA;
    const real *  chargeB;
    real          sigma6_min, sigma6_def, lam_power, sc_power, sc_r_power;
    real          alpha_coul, alpha_vdw, lambda_coul, lambda_vdw, ewc_lj;
    real          ewcljrsq, ewclj, ewclj2, exponent, poly, vvdw_disp, vvdw_rep, sh_lj_ewald;
    real          ewclj6;
    const real *  nbfp, *nbfp_grid;
    real *        dvdl;
    real *        Vv;
    real *        Vc;
    gmx_bool      bDoForces, bDoShiftForces, bDoPotential;
    real          rcoulomb, rvdw, sh_invrc6;
    gmx_bool      bExactElecCutoff, bExactVdwCutoff, bExactCutoffAll;
    gmx_bool      bEwald, bEwaldLJ;
    real          rcutoff_max2;
    const real *  tab_ewald_F_lj;
    const real *  tab_ewald_V_lj;
    real          d, d2, sw, dsw, rinvcorr;
    real          elec_swV3, elec_swV4, elec_swV5, elec_swF2, elec_swF3, elec_swF4;
    real          vdw_swV3, vdw_swV4, vdw_swV5, vdw_swF2, vdw_swF3, vdw_swF4;
    gmx_bool      bConvertEwaldToCoulomb, bConvertLJEwaldToLJ6;
    gmx_bool      bComputeVdwInteraction, bComputeElecInteraction;
    const real *  ewtab;
    int           ewitab;
    real          ewrt, eweps, ewtabscale, ewtabhalfspace, sh_ewald;

    sh_ewald            = fr->ic->sh_ewald;
    ewtab               = fr->ic->tabq_coul_FDV0;
    ewtabscale          = fr->ic->tabq_scale;
    ewtabhalfspace      = 0.5/ewtabscale;
    tab_ewald_F_lj      = fr->ic->tabq_vdw_F;
    tab_ewald_V_lj      = fr->ic->tabq_vdw_V;

    x                   = xx[0];
    f                   = ff[0];

    fshift              = fr->fshift[0];

    nri                 = nlist->nri;
    iinr                = nlist->iinr;
    jindex              = nlist->jindex;
    jjnr                = nlist->jjnr;
    icoul               = nlist->ielec;
    ivdw                = nlist->ivdw;
    shift               = nlist->shift;
    gid                 = nlist->gid;

    shiftvec            = fr->shift_vec[0];
    chargeA             = mdatoms->chargeA;
    chargeB             = mdatoms->chargeB;
    facel               = fr->epsfac;
    krf                 = fr->k_rf;
    crf                 = fr->c_rf;
    ewc_lj              = fr->ewaldcoeff_lj;
    Vc                  = kernel_data->energygrp_elec;
    typeA               = mdatoms->typeA;
    typeB               = mdatoms->typeB;
    ntype               = fr->ntype;
    nbfp                = fr->nbfp;
    nbfp_grid           = fr->ljpme_c6grid;
    Vv                  = kernel_data->energygrp_vdw;
    lambda_coul         = kernel_data->lambda[efptCOUL];
    lambda_vdw          = kernel_data->lambda[efptVDW];
    dvdl                = kernel_data->dvdl;
    alpha_coul          = fr->sc_alphacoul;
    alpha_vdw           = fr->sc_alphavdw;
    lam_power           = fr->sc_power;
    sc_r_power          = fr->sc_r_power;
    sigma6_def          = fr->sc_sigma6_def;
    sigma6_min          = fr->sc_sigma6_min;
    bDoForces           = kernel_data->flags & GMX_NONBONDED_DO_FORCE;
    bDoShiftForces      = kernel_data->flags & GMX_NONBONDED_DO_SHIFTFORCE;
    bDoPotential        = kernel_data->flags & GMX_NONBONDED_DO_POTENTIAL;

    rcoulomb            = fr->rcoulomb;
    sh_ewald            = fr->ic->sh_ewald;
    rvdw                = fr->rvdw;
    sh_invrc6           = fr->ic->sh_invrc6;
    sh_lj_ewald         = fr->ic->sh_lj_ewald;
    ewclj               = fr->ewaldcoeff_lj;
    ewclj2              = ewclj*ewclj;
    ewclj6              = ewclj2*ewclj2*ewclj2;

    if (fr->coulomb_modifier == eintmodPOTSWITCH)
    {
        d               = fr->rcoulomb-fr->rcoulomb_switch;
        elec_swV3       = -10.0/(d*d*d);
        elec_swV4       =  15.0/(d*d*d*d);
        elec_swV5       =  -6.0/(d*d*d*d*d);
        elec_swF2       = -30.0/(d*d*d);
        elec_swF3       =  60.0/(d*d*d*d);
        elec_swF4       = -30.0/(d*d*d*d*d);
    }
    else
    {
        /* Avoid warnings from stupid compilers (looking at you, Clang!) */
        elec_swV3 = elec_swV4 = elec_swV5 = elec_swF2 = elec_swF3 = elec_swF4 = 0.0;
    }

    if (fr->vdw_modifier == eintmodPOTSWITCH)
    {
        d               = fr->rvdw-fr->rvdw_switch;
        vdw_swV3        = -10.0/(d*d*d);
        vdw_swV4        =  15.0/(d*d*d*d);
        vdw_swV5        =  -6.0/(d*d*d*d*d);
        vdw_swF2        = -30.0/(d*d*d);
        vdw_swF3        =  60.0/(d*d*d*d);
        vdw_swF4        = -30.0/(d*d*d*d*d);
    }
    else
    {
        /* Avoid warnings from stupid compilers (looking at you, Clang!) */
        vdw_swV3 = vdw_swV4 = vdw_swV5 = vdw_swF2 = vdw_swF3 = vdw_swF4 = 0.0;
    }

    if (fr->cutoff_scheme == ecutsVERLET)
    {
        const interaction_const_t *ic;

        ic = fr->ic;
        if (EVDW_PME(ic->vdwtype))
        {
            ivdw         = GMX_NBKERNEL_VDW_LJEWALD;
        }
        else
        {
            ivdw         = GMX_NBKERNEL_VDW_LENNARDJONES;
        }

        if (ic->eeltype == eelCUT || EEL_RF(ic->eeltype))
        {
            icoul        = GMX_NBKERNEL_ELEC_REACTIONFIELD;
        }
        else if (EEL_PME_EWALD(ic->eeltype))
        {
            icoul        = GMX_NBKERNEL_ELEC_EWALD;
        }
        else
        {
            gmx_incons("Unsupported eeltype with Verlet and free-energy");
        }

        bExactElecCutoff = TRUE;
        bExactVdwCutoff  = TRUE;
    }
    else
    {
        bExactElecCutoff = (fr->coulomb_modifier != eintmodNONE) || fr->eeltype == eelRF_ZERO;
        bExactVdwCutoff  = (fr->vdw_modifier != eintmodNONE);
    }

    bExactCutoffAll = (bExactElecCutoff && bExactVdwCutoff);
    rcutoff_max2    = max(fr->rcoulomb, fr->rvdw);
    rcutoff_max2    = rcutoff_max2*rcutoff_max2;

    bEwald          = (icoul == GMX_NBKERNEL_ELEC_EWALD);
    bEwaldLJ        = (ivdw == GMX_NBKERNEL_VDW_LJEWALD);

    /* For Ewald/PME interactions we cannot easily apply the soft-core component to
     * reciprocal space. When we use vanilla (not switch/shift) Ewald interactions, we
     * can apply the small trick of subtracting the _reciprocal_ space contribution
     * in this kernel, and instead apply the free energy interaction to the 1/r
     * (standard coulomb) interaction.
     *
     * However, we cannot use this approach for switch-modified since we would then
     * effectively end up evaluating a significantly different interaction here compared to the
     * normal (non-free-energy) kernels, either by applying a cutoff at a different
     * position than what the user requested, or by switching different
     * things (1/r rather than short-range Ewald). For these settings, we just
     * use the traditional short-range Ewald interaction in that case.
     */
    bConvertEwaldToCoulomb = (bEwald && (fr->coulomb_modifier != eintmodPOTSWITCH));
    /* For now the below will always be true (since LJ-PME only works with Shift in Gromacs-5.0),
     * but writing it this way means we stay in sync with coulomb, and it avoids future bugs.
     */
    bConvertLJEwaldToLJ6   = (bEwaldLJ && (fr->vdw_modifier   != eintmodPOTSWITCH));

    /* We currently don't implement exclusion correction, needed with the Verlet cut-off scheme, without conversion */
    if (fr->cutoff_scheme == ecutsVERLET &&
        ((bEwald   && !bConvertEwaldToCoulomb) ||
         (bEwaldLJ && !bConvertLJEwaldToLJ6)))
    {
        gmx_incons("Unimplemented non-bonded setup");
    }

    /* fix compiler warnings */
    nj1   = 0;
    n1C   = n1V   = 0;
    epsC  = epsV  = 0;
    eps2C = eps2V = 0;

    dvdl_coul  = 0;
    dvdl_vdw   = 0;

    /* Lambda factor for state A, 1-lambda*/
    LFC[STATE_A] = 1.0 - lambda_coul;
    LFV[STATE_A] = 1.0 - lambda_vdw;

    /* Lambda factor for state B, lambda*/
    LFC[STATE_B] = lambda_coul;
    LFV[STATE_B] = lambda_vdw;

    /*derivative of the lambda factor for state A and B */
    DLF[STATE_A] = -1;
    DLF[STATE_B] = 1;

    for (i = 0; i < NSTATES; i++)
    {
        lfac_coul[i]  = (lam_power == 2 ? (1-LFC[i])*(1-LFC[i]) : (1-LFC[i]));
        dlfac_coul[i] = DLF[i]*lam_power/sc_r_power*(lam_power == 2 ? (1-LFC[i]) : 1);
        lfac_vdw[i]   = (lam_power == 2 ? (1-LFV[i])*(1-LFV[i]) : (1-LFV[i]));
        dlfac_vdw[i]  = DLF[i]*lam_power/sc_r_power*(lam_power == 2 ? (1-LFV[i]) : 1);
    }
    /* precalculate */
    sigma2_def = pow(sigma6_def, 1.0/3.0);
    sigma2_min = pow(sigma6_min, 1.0/3.0);

    /* Ewald (not PME) table is special (icoul==enbcoulFEWALD) */

    do_tab = (icoul == GMX_NBKERNEL_ELEC_CUBICSPLINETABLE ||
              ivdw == GMX_NBKERNEL_VDW_CUBICSPLINETABLE);
    if (do_tab)
    {
        tabscale         = kernel_data->table_elec_vdw->scale;
        VFtab            = kernel_data->table_elec_vdw->data;
        /* we always use the combined table here */
        tab_elemsize     = 12;
    }

    for (n = 0; (n < nri); n++)
    {
        int npair_within_cutoff;

        npair_within_cutoff = 0;

        is3              = 3*shift[n];
        shX              = shiftvec[is3];
        shY              = shiftvec[is3+1];
        shZ              = shiftvec[is3+2];
        nj0              = jindex[n];
        nj1              = jindex[n+1];
        ii               = iinr[n];
        ii3              = 3*ii;
        ix               = shX + x[ii3+0];
        iy               = shY + x[ii3+1];
        iz               = shZ + x[ii3+2];
        iqA              = facel*chargeA[ii];
        iqB              = facel*chargeB[ii];
        ntiA             = 2*ntype*typeA[ii];
        ntiB             = 2*ntype*typeB[ii];
        vctot            = 0;
        vvtot            = 0;
        fix              = 0;
        fiy              = 0;
        fiz              = 0;

        for (k = nj0; (k < nj1); k++)
        {
            jnr              = jjnr[k];
            j3               = 3*jnr;
            dx               = ix - x[j3];
            dy               = iy - x[j3+1];
            dz               = iz - x[j3+2];
            rsq              = dx*dx + dy*dy + dz*dz;

            if (bExactCutoffAll && rsq >= rcutoff_max2)
            {
                /* We save significant time by skipping all code below.
                 * Note that with soft-core interactions, the actual cut-off
                 * check might be different. But since the soft-core distance
                 * is always larger than r, checking on r here is safe.
                 */
                continue;
            }
            npair_within_cutoff++;

            if (rsq > 0)
            {
                rinv         = gmx_invsqrt(rsq);
                r            = rsq*rinv;
            }
            else
            {
                /* The force at r=0 is zero, because of symmetry.
                 * But note that the potential is in general non-zero,
                 * since the soft-cored r will be non-zero.
                 */
                rinv         = 0;
                r            = 0;
            }

            if (sc_r_power == 6.0)
            {
                rpm2             = rsq*rsq;  /* r4 */
                rp               = rpm2*rsq; /* r6 */
            }
            else if (sc_r_power == 48.0)
            {
                rp               = rsq*rsq*rsq; /* r6 */
                rp               = rp*rp;       /* r12 */
                rp               = rp*rp;       /* r24 */
                rp               = rp*rp;       /* r48 */
                rpm2             = rp/rsq;      /* r46 */
            }
            else
            {
                rp             = pow(r, sc_r_power);  /* not currently supported as input, but can handle it */
                rpm2           = rp/rsq;
            }

            Fscal = 0;

            qq[STATE_A]      = iqA*chargeA[jnr];
            qq[STATE_B]      = iqB*chargeB[jnr];

            tj[STATE_A]      = ntiA+2*typeA[jnr];
            tj[STATE_B]      = ntiB+2*typeB[jnr];

            if (nlist->excl_fep == NULL || nlist->excl_fep[k])
            {
                c6[STATE_A]      = nbfp[tj[STATE_A]];
                c6[STATE_B]      = nbfp[tj[STATE_B]];

                for (i = 0; i < NSTATES; i++)
                {
                    c12[i]             = nbfp[tj[i]+1];
                    if ((c6[i] > 0) && (c12[i] > 0))
                    {
                        /* c12 is stored scaled with 12.0 and c6 is scaled with 6.0 - correct for this */
                        sigma6[i]       = 0.5*c12[i]/c6[i];
                        sigma2[i]       = pow(sigma6[i], 1.0/3.0);
                        /* should be able to get rid of this ^^^ internal pow call eventually.  Will require agreement on
                           what data to store externally.  Can't be fixed without larger scale changes, so not 4.6 */
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
                    if (sc_r_power == 6.0)
                    {
                        sigma_pow[i]    = sigma6[i];
                        sigma_powm2[i]  = sigma6[i]/sigma2[i];
                    }
                    else if (sc_r_power == 48.0)
                    {
                        sigma_pow[i]    = sigma6[i]*sigma6[i];       /* sigma^12 */
                        sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^24 */
                        sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^48 */
                        sigma_powm2[i]  = sigma_pow[i]/sigma2[i];
                    }
                    else
                    {    /* not really supported as input, but in here for testing the general case*/
                        sigma_pow[i]    = pow(sigma2[i], sc_r_power/2);
                        sigma_powm2[i]  = sigma_pow[i]/(sigma2[i]);
                    }
                }

                /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
                if ((c12[STATE_A] > 0) && (c12[STATE_B] > 0))
                {
                    alpha_vdw_eff    = 0;
                    alpha_coul_eff   = 0;
                }
                else
                {
                    alpha_vdw_eff    = alpha_vdw;
                    alpha_coul_eff   = alpha_coul;
                }

                for (i = 0; i < NSTATES; i++)
                {
                    FscalC[i]    = 0;
                    FscalV[i]    = 0;
                    Vcoul[i]     = 0;
                    Vvdw[i]      = 0;

                    /* Only spend time on A or B state if it is non-zero */
                    if ( (qq[i] != 0) || (c6[i] != 0) || (c12[i] != 0) )
                    {
                        /* this section has to be inside the loop because of the dependence on sigma_pow */
                        rpinvC         = 1.0/(alpha_coul_eff*lfac_coul[i]*sigma_pow[i]+rp);
                        rinvC          = pow(rpinvC, 1.0/sc_r_power);
                        rC             = 1.0/rinvC;

                        rpinvV         = 1.0/(alpha_vdw_eff*lfac_vdw[i]*sigma_pow[i]+rp);
                        rinvV          = pow(rpinvV, 1.0/sc_r_power);
                        rV             = 1.0/rinvV;

                        if (do_tab)
                        {
                            rtC        = rC*tabscale;
                            n0         = rtC;
                            epsC       = rtC-n0;
                            eps2C      = epsC*epsC;
                            n1C        = tab_elemsize*n0;

                            rtV        = rV*tabscale;
                            n0         = rtV;
                            epsV       = rtV-n0;
                            eps2V      = epsV*epsV;
                            n1V        = tab_elemsize*n0;
                        }

                        /* Only process the coulomb interactions if we have charges,
                         * and if we either include all entries in the list (no cutoff
                         * used in the kernel), or if we are within the cutoff.
                         */
                        bComputeElecInteraction = !bExactElecCutoff ||
                            ( bConvertEwaldToCoulomb && r < rcoulomb) ||
                            (!bConvertEwaldToCoulomb && rC < rcoulomb);

                        if ( (qq[i] != 0) && bComputeElecInteraction)
                        {
                            switch (icoul)
                            {
                                case GMX_NBKERNEL_ELEC_COULOMB:
                                    /* simple cutoff */
                                    Vcoul[i]   = qq[i]*rinvC;
                                    FscalC[i]  = Vcoul[i];
                                    /* The shift for the Coulomb potential is stored in
                                     * the RF parameter c_rf, which is 0 without shift.
                                     */
                                    Vcoul[i]  -= qq[i]*fr->ic->c_rf;
                                    break;

                                case GMX_NBKERNEL_ELEC_REACTIONFIELD:
                                    /* reaction-field */
                                    Vcoul[i]   = qq[i]*(rinvC + krf*rC*rC-crf);
                                    FscalC[i]  = qq[i]*(rinvC - 2.0*krf*rC*rC);
                                    break;

                                case GMX_NBKERNEL_ELEC_CUBICSPLINETABLE:
                                    /* non-Ewald tabulated coulomb */
                                    nnn        = n1C;
                                    Y          = VFtab[nnn];
                                    F          = VFtab[nnn+1];
                                    Geps       = epsC*VFtab[nnn+2];
                                    Heps2      = eps2C*VFtab[nnn+3];
                                    Fp         = F+Geps+Heps2;
                                    VV         = Y+epsC*Fp;
                                    FF         = Fp+Geps+2.0*Heps2;
                                    Vcoul[i]   = qq[i]*VV;
                                    FscalC[i]  = -qq[i]*tabscale*FF*rC;
                                    break;

                                case GMX_NBKERNEL_ELEC_GENERALIZEDBORN:
                                    gmx_fatal(FARGS, "Free energy and GB not implemented.\n");
                                    break;

                                case GMX_NBKERNEL_ELEC_EWALD:
                                    if (bConvertEwaldToCoulomb)
                                    {
                                        /* Ewald FEP is done only on the 1/r part */
                                        Vcoul[i]   = qq[i]*(rinvC-sh_ewald);
                                        FscalC[i]  = qq[i]*rinvC;
                                    }
                                    else
                                    {
                                        ewrt      = rC*ewtabscale;
                                        ewitab    = (int) ewrt;
                                        eweps     = ewrt-ewitab;
                                        ewitab    = 4*ewitab;
                                        FscalC[i] = ewtab[ewitab]+eweps*ewtab[ewitab+1];
                                        rinvcorr  = rinvC-sh_ewald;
                                        Vcoul[i]  = qq[i]*(rinvcorr-(ewtab[ewitab+2]-ewtabhalfspace*eweps*(ewtab[ewitab]+FscalC[i])));
                                        FscalC[i] = qq[i]*(rinvC-rC*FscalC[i]);
                                    }
                                    break;

                                case GMX_NBKERNEL_ELEC_NONE:
                                    FscalC[i]  = 0.0;
                                    Vcoul[i]   = 0.0;
                                    break;

                                default:
                                    gmx_incons("Invalid icoul in free energy kernel");
                                    break;
                            }

                            if (fr->coulomb_modifier == eintmodPOTSWITCH)
                            {
                                d                = rC-fr->rcoulomb_switch;
                                d                = (d > 0.0) ? d : 0.0;
                                d2               = d*d;
                                sw               = 1.0+d2*d*(elec_swV3+d*(elec_swV4+d*elec_swV5));
                                dsw              = d2*(elec_swF2+d*(elec_swF3+d*elec_swF4));

                                FscalC[i]        = FscalC[i]*sw - rC*Vcoul[i]*dsw;
                                Vcoul[i]        *= sw;

                                FscalC[i]        = (rC < rcoulomb) ? FscalC[i] : 0.0;
                                Vcoul[i]         = (rC < rcoulomb) ? Vcoul[i] : 0.0;
                            }
                        }

                        /* Only process the VDW interactions if we have
                         * some non-zero parameters, and if we either
                         * include all entries in the list (no cutoff used
                         * in the kernel), or if we are within the cutoff.
                         */
                        bComputeVdwInteraction = !bExactVdwCutoff ||
                            ( bConvertLJEwaldToLJ6 && r < rvdw) ||
                            (!bConvertLJEwaldToLJ6 && rV < rvdw);
                        if ((c6[i] != 0 || c12[i] != 0) && bComputeVdwInteraction)
                        {
                            switch (ivdw)
                            {
                                case GMX_NBKERNEL_VDW_LENNARDJONES:
                                    /* cutoff LJ */
                                    if (sc_r_power == 6.0)
                                    {
                                        rinv6            = rpinvV;
                                    }
                                    else
                                    {
                                        rinv6            = rinvV*rinvV;
                                        rinv6            = rinv6*rinv6*rinv6;
                                    }
                                    Vvdw6            = c6[i]*rinv6;
                                    Vvdw12           = c12[i]*rinv6*rinv6;

                                    Vvdw[i]          = ( (Vvdw12 - c12[i]*sh_invrc6*sh_invrc6)*(1.0/12.0)
                                                         - (Vvdw6 - c6[i]*sh_invrc6)*(1.0/6.0));
                                    FscalV[i]        = Vvdw12 - Vvdw6;
                                    break;

                                case GMX_NBKERNEL_VDW_BUCKINGHAM:
                                    gmx_fatal(FARGS, "Buckingham free energy not supported.");
                                    break;

                                case GMX_NBKERNEL_VDW_CUBICSPLINETABLE:
                                    /* Table LJ */
                                    nnn = n1V+4;
                                    /* dispersion */
                                    Y          = VFtab[nnn];
                                    F          = VFtab[nnn+1];
                                    Geps       = epsV*VFtab[nnn+2];
                                    Heps2      = eps2V*VFtab[nnn+3];
                                    Fp         = F+Geps+Heps2;
                                    VV         = Y+epsV*Fp;
                                    FF         = Fp+Geps+2.0*Heps2;
                                    Vvdw[i]   += c6[i]*VV;
                                    FscalV[i] -= c6[i]*tabscale*FF*rV;

                                    /* repulsion */
                                    Y          = VFtab[nnn+4];
                                    F          = VFtab[nnn+5];
                                    Geps       = epsV*VFtab[nnn+6];
                                    Heps2      = eps2V*VFtab[nnn+7];
                                    Fp         = F+Geps+Heps2;
                                    VV         = Y+epsV*Fp;
                                    FF         = Fp+Geps+2.0*Heps2;
                                    Vvdw[i]   += c12[i]*VV;
                                    FscalV[i] -= c12[i]*tabscale*FF*rV;
                                    break;

                                case GMX_NBKERNEL_VDW_LJEWALD:
                                    if (sc_r_power == 6.0)
                                    {
                                        rinv6            = rpinvV;
                                    }
                                    else
                                    {
                                        rinv6            = rinvV*rinvV;
                                        rinv6            = rinv6*rinv6*rinv6;
                                    }
                                    c6grid           = nbfp_grid[tj[i]];

                                    if (bConvertLJEwaldToLJ6)
                                    {
                                        /* cutoff LJ */
                                        Vvdw6            = c6[i]*rinv6;
                                        Vvdw12           = c12[i]*rinv6*rinv6;

                                        Vvdw[i]          = ( (Vvdw12 - c12[i]*sh_invrc6*sh_invrc6)*(1.0/12.0)
                                                             - (Vvdw6 - c6[i]*sh_invrc6 - c6grid*sh_lj_ewald)*(1.0/6.0));
                                        FscalV[i]        = Vvdw12 - Vvdw6;
                                    }
                                    else
                                    {
                                        /* Normal LJ-PME */
                                        ewcljrsq         = ewclj2*rV*rV;
                                        exponent         = exp(-ewcljrsq);
                                        poly             = exponent*(1.0 + ewcljrsq + ewcljrsq*ewcljrsq*0.5);
                                        vvdw_disp        = (c6[i]-c6grid*(1.0-poly))*rinv6;
                                        vvdw_rep         = c12[i]*rinv6*rinv6;
                                        FscalV[i]        = vvdw_rep - vvdw_disp - c6grid*(1.0/6.0)*exponent*ewclj6;
                                        Vvdw[i]          = (vvdw_rep - c12[i]*sh_invrc6*sh_invrc6)/12.0 - (vvdw_disp - c6[i]*sh_invrc6 - c6grid*sh_lj_ewald)/6.0;
                                    }
                                    break;

                                case GMX_NBKERNEL_VDW_NONE:
                                    Vvdw[i]    = 0.0;
                                    FscalV[i]  = 0.0;
                                    break;

                                default:
                                    gmx_incons("Invalid ivdw in free energy kernel");
                                    break;
                            }

                            if (fr->vdw_modifier == eintmodPOTSWITCH)
                            {
                                d                = rV-fr->rvdw_switch;
                                d                = (d > 0.0) ? d : 0.0;
                                d2               = d*d;
                                sw               = 1.0+d2*d*(vdw_swV3+d*(vdw_swV4+d*vdw_swV5));
                                dsw              = d2*(vdw_swF2+d*(vdw_swF3+d*vdw_swF4));

                                FscalV[i]        = FscalV[i]*sw - rV*Vvdw[i]*dsw;
                                Vvdw[i]         *= sw;

                                FscalV[i]  = (rV < rvdw) ? FscalV[i] : 0.0;
                                Vvdw[i]    = (rV < rvdw) ? Vvdw[i] : 0.0;
                            }
                        }

                        /* FscalC (and FscalV) now contain: dV/drC * rC
                         * Now we multiply by rC^-p, so it will be: dV/drC * rC^1-p
                         * Further down we first multiply by r^p-2 and then by
                         * the vector r, which in total gives: dV/drC * (r/rC)^1-p
                         */
                        FscalC[i] *= rpinvC;
                        FscalV[i] *= rpinvV;
                    }
                }

                /* Assemble A and B states */
                for (i = 0; i < NSTATES; i++)
                {
                    vctot         += LFC[i]*Vcoul[i];
                    vvtot         += LFV[i]*Vvdw[i];

                    Fscal         += LFC[i]*FscalC[i]*rpm2;
                    Fscal         += LFV[i]*FscalV[i]*rpm2;

                    dvdl_coul     += Vcoul[i]*DLF[i] + LFC[i]*alpha_coul_eff*dlfac_coul[i]*FscalC[i]*sigma_pow[i];
                    dvdl_vdw      += Vvdw[i]*DLF[i] + LFV[i]*alpha_vdw_eff*dlfac_vdw[i]*FscalV[i]*sigma_pow[i];
                }
            }
            else if (icoul == GMX_NBKERNEL_ELEC_REACTIONFIELD)
            {
                /* For excluded pairs, which are only in this pair list when
                 * using the Verlet scheme, we don't use soft-core.
                 * The group scheme also doesn't soft-core for these.
                 * As there is no singularity, there is no need for soft-core.
                 */
                VV = krf*rsq - crf;
                FF = -2.0*krf;

                if (ii == jnr)
                {
                    VV *= 0.5;
                }

                for (i = 0; i < NSTATES; i++)
                {
                    vctot      += LFC[i]*qq[i]*VV;
                    Fscal      += LFC[i]*qq[i]*FF;
                    dvdl_coul  += DLF[i]*qq[i]*VV;
                }
            }

            if (bConvertEwaldToCoulomb && ( !bExactElecCutoff || r < rcoulomb ) )
            {
                /* See comment in the preamble. When using Ewald interactions
                 * (unless we use a switch modifier) we subtract the reciprocal-space
                 * Ewald component here which made it possible to apply the free
                 * energy interaction to 1/r (vanilla coulomb short-range part)
                 * above. This gets us closer to the ideal case of applying
                 * the softcore to the entire electrostatic interaction,
                 * including the reciprocal-space component.
                 */
                real v_lr, f_lr;

                ewrt      = r*ewtabscale;
                ewitab    = (int) ewrt;
                eweps     = ewrt-ewitab;
                ewitab    = 4*ewitab;
                f_lr      = ewtab[ewitab]+eweps*ewtab[ewitab+1];
                v_lr      = (ewtab[ewitab+2]-ewtabhalfspace*eweps*(ewtab[ewitab]+f_lr));
                f_lr     *= rinv;

                /* Note that any possible Ewald shift has already been applied in
                 * the normal interaction part above.
                 */

                if (ii == jnr)
                {
                    /* If we get here, the i particle (ii) has itself (jnr)
                     * in its neighborlist. This can only happen with the Verlet
                     * scheme, and corresponds to a self-interaction that will
                     * occur twice. Scale it down by 50% to only include it once.
                     */
                    v_lr *= 0.5;
                }

                for (i = 0; i < NSTATES; i++)
                {
                    vctot      -= LFC[i]*qq[i]*v_lr;
                    Fscal      -= LFC[i]*qq[i]*f_lr;
                    dvdl_coul  -= (DLF[i]*qq[i])*v_lr;
                }
            }

            if (bConvertLJEwaldToLJ6 && (!bExactVdwCutoff || r < rvdw))
            {
                /* See comment in the preamble. When using LJ-Ewald interactions
                 * (unless we use a switch modifier) we subtract the reciprocal-space
                 * Ewald component here which made it possible to apply the free
                 * energy interaction to r^-6 (vanilla LJ6 short-range part)
                 * above. This gets us closer to the ideal case of applying
                 * the softcore to the entire VdW interaction,
                 * including the reciprocal-space component.
                 */
                /* We could also use the analytical form here
                 * iso a table, but that can cause issues for
                 * r close to 0 for non-interacting pairs.
                 */
                real rs, frac, f_lr;
                int  ri;

                rs     = rsq*rinv*ewtabscale;
                ri     = (int)rs;
                frac   = rs - ri;
                f_lr   = (1 - frac)*tab_ewald_F_lj[ri] + frac*tab_ewald_F_lj[ri+1];
                /* TODO: Currently the Ewald LJ table does not contain
                 * the factor 1/6, we should add this.
                 */
                FF     = f_lr*rinv/6.0;
                VV     = (tab_ewald_V_lj[ri] - ewtabhalfspace*frac*(tab_ewald_F_lj[ri] + f_lr))/6.0;

                if (ii == jnr)
                {
                    /* If we get here, the i particle (ii) has itself (jnr)
                     * in its neighborlist. This can only happen with the Verlet
                     * scheme, and corresponds to a self-interaction that will
                     * occur twice. Scale it down by 50% to only include it once.
                     */
                    VV *= 0.5;
                }

                for (i = 0; i < NSTATES; i++)
                {
                    c6grid      = nbfp_grid[tj[i]];
                    vvtot      += LFV[i]*c6grid*VV;
                    Fscal      += LFV[i]*c6grid*FF;
                    dvdl_vdw   += (DLF[i]*c6grid)*VV;
                }
            }

            if (bDoForces)
            {
                tx         = Fscal*dx;
                ty         = Fscal*dy;
                tz         = Fscal*dz;
                fix        = fix + tx;
                fiy        = fiy + ty;
                fiz        = fiz + tz;
                /* OpenMP atomics are expensive, but this kernels is also
                 * expensive, so we can take this hit, instead of using
                 * thread-local output buffers and extra reduction.
                 */
#pragma omp atomic
                f[j3]     -= tx;
#pragma omp atomic
                f[j3+1]   -= ty;
#pragma omp atomic
                f[j3+2]   -= tz;
            }
        }

        /* The atomics below are expensive with many OpenMP threads.
         * Here unperturbed i-particles will usually only have a few
         * (perturbed) j-particles in the list. Thus with a buffered list
         * we can skip a significant number of i-reductions with a check.
         */
        if (npair_within_cutoff > 0)
        {
            if (bDoForces)
            {
#pragma omp atomic
                f[ii3]        += fix;
#pragma omp atomic
                f[ii3+1]      += fiy;
#pragma omp atomic
                f[ii3+2]      += fiz;
            }
            if (bDoShiftForces)
            {
#pragma omp atomic
                fshift[is3]   += fix;
#pragma omp atomic
                fshift[is3+1] += fiy;
#pragma omp atomic
                fshift[is3+2] += fiz;
            }
            if (bDoPotential)
            {
                ggid               = gid[n];
#pragma omp atomic
                Vc[ggid]          += vctot;
#pragma omp atomic
                Vv[ggid]          += vvtot;
            }
        }
    }

#pragma omp atomic
    dvdl[efptCOUL]     += dvdl_coul;
 #pragma omp atomic
    dvdl[efptVDW]      += dvdl_vdw;

    /* Estimate flops, average for free energy stuff:
     * 12  flops per outer iteration
     * 150 flops per inner iteration
     */
#pragma omp atomic
    inc_nrnb(nrnb, eNR_NBKERNEL_FREE_ENERGY, nlist->nri*12 + nlist->jindex[n]*150);
}

real
nb_free_energy_evaluate_single(real r2, real sc_r_power, real alpha_coul, real alpha_vdw,
                               real tabscale, real *vftab,
                               real qqA, real c6A, real c12A, real qqB, real c6B, real c12B,
                               real LFC[2], real LFV[2], real DLF[2],
                               real lfac_coul[2], real lfac_vdw[2], real dlfac_coul[2], real dlfac_vdw[2],
                               real sigma6_def, real sigma6_min, real sigma2_def, real sigma2_min,
                               real *velectot, real *vvdwtot, real *dvdl)
{
    real       r, rp, rpm2, rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VV, FF, fscal;
    real       qq[2], c6[2], c12[2], sigma6[2], sigma2[2], sigma_pow[2], sigma_powm2[2];
    real       alpha_coul_eff, alpha_vdw_eff, dvdl_coul, dvdl_vdw;
    real       rpinv, r_coul, r_vdw, velecsum, vvdwsum;
    real       fscal_vdw[2], fscal_elec[2];
    real       velec[2], vvdw[2];
    int        i, ntab;

    qq[0]    = qqA;
    qq[1]    = qqB;
    c6[0]    = c6A;
    c6[1]    = c6B;
    c12[0]   = c12A;
    c12[1]   = c12B;

    if (sc_r_power == 6.0)
    {
        rpm2             = r2*r2;   /* r4 */
        rp               = rpm2*r2; /* r6 */
    }
    else if (sc_r_power == 48.0)
    {
        rp               = r2*r2*r2; /* r6 */
        rp               = rp*rp;    /* r12 */
        rp               = rp*rp;    /* r24 */
        rp               = rp*rp;    /* r48 */
        rpm2             = rp/r2;    /* r46 */
    }
    else
    {
        rp             = pow(r2, 0.5*sc_r_power);  /* not currently supported as input, but can handle it */
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
            sigma6[i]       = 0.5*c12[i]/c6[i];
            sigma2[i]       = pow(0.5*c12[i]/c6[i], 1.0/3.0);
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
        if (sc_r_power == 6.0)
        {
            sigma_pow[i]    = sigma6[i];
            sigma_powm2[i]  = sigma6[i]/sigma2[i];
        }
        else if (sc_r_power == 48.0)
        {
            sigma_pow[i]    = sigma6[i]*sigma6[i];       /* sigma^12 */
            sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^24 */
            sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^48 */
            sigma_powm2[i]  = sigma_pow[i]/sigma2[i];
        }
        else
        {    /* not really supported as input, but in here for testing the general case*/
            sigma_pow[i]    = pow(sigma2[i], sc_r_power/2);
            sigma_powm2[i]  = sigma_pow[i]/(sigma2[i]);
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
            rpinv            = 1.0/(alpha_coul_eff*lfac_coul[i]*sigma_pow[i]+rp);
            r_coul           = pow(rpinv, -1.0/sc_r_power);

            /* Electrostatics table lookup data */
            rtab             = r_coul*tabscale;
            ntab             = rtab;
            eps              = rtab-ntab;
            eps2             = eps*eps;
            ntab             = 12*ntab;
            /* Electrostatics */
            Y                = vftab[ntab];
            F                = vftab[ntab+1];
            Geps             = eps*vftab[ntab+2];
            Heps2            = eps2*vftab[ntab+3];
            Fp               = F+Geps+Heps2;
            VV               = Y+eps*Fp;
            FF               = Fp+Geps+2.0*Heps2;
            velec[i]         = qq[i]*VV;
            fscal_elec[i]    = -qq[i]*FF*r_coul*rpinv*tabscale;

            /* Vdw */
            rpinv            = 1.0/(alpha_vdw_eff*lfac_vdw[i]*sigma_pow[i]+rp);
            r_vdw            = pow(rpinv, -1.0/sc_r_power);
            /* Vdw table lookup data */
            rtab             = r_vdw*tabscale;
            ntab             = rtab;
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
            FF               = Fp+Geps+2.0*Heps2;
            vvdw[i]          = c6[i]*VV;
            fscal_vdw[i]     = -c6[i]*FF;

            /* Repulsion */
            Y                = vftab[ntab+8];
            F                = vftab[ntab+9];
            Geps             = eps*vftab[ntab+10];
            Heps2            = eps2*vftab[ntab+11];
            Fp               = F+Geps+Heps2;
            VV               = Y+eps*Fp;
            FF               = Fp+Geps+2.0*Heps2;
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
