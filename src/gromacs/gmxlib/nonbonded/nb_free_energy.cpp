/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "nb_free_energy.h"

#include <cmath>

#include <algorithm>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"


//! Enum for templating the soft-core treatment in the kernel
enum class SoftCoreTreatment
{
    None,    //!< No soft-core
    RPower6, //!< Soft-core with r-power = 6
    RPower48 //!< Soft-core with r-power = 48
};

//! Most treatments are fine with float in mixed-precision mode.
template <SoftCoreTreatment softCoreTreatment>
struct SoftCoreReal
{
    //! Real type for soft-core calculations
    using Real = real;
};

//! This treatment requires double precision for some computations.
template <>
struct SoftCoreReal<SoftCoreTreatment::RPower48>
{
    //! Real type for soft-core calculations
    using Real = double;
};

//! Computes r^(1/p) and 1/r^(1/p) for the standard p=6
template <SoftCoreTreatment softCoreTreatment>
static inline void pthRoot(const real  r,
                           real       *pthRoot,
                           real       *invPthRoot)
{
    *invPthRoot = gmx::invsqrt(std::cbrt(r));
    *pthRoot    = 1/(*invPthRoot);
}

// We need a double version to make the specialization below work
#if !GMX_DOUBLE
//! Computes r^(1/p) and 1/r^(1/p) for the standard p=6
template <SoftCoreTreatment softCoreTreatment>
static inline void pthRoot(const double  r,
                           real         *pthRoot,
                           double       *invPthRoot)
{
    *invPthRoot = gmx::invsqrt(std::cbrt(r));
    *pthRoot    = 1/(*invPthRoot);
}
#endif

//! Computes r^(1/p) and 1/r^(1/p) for p=48
template <>
inline void pthRoot<SoftCoreTreatment::RPower48>(const double  r,
                                                 real         *pthRoot,
                                                 double       *invPthRoot)
{
    *pthRoot    = std::pow(r, 1.0/48.0);
    *invPthRoot = 1/(*pthRoot);
}

//! Templated free-energy non-bonded kernel
template<SoftCoreTreatment softCoreTreatment, bool scLambdasOrAlphasDiffer>
static void
nb_free_energy_kernel(const t_nblist * gmx_restrict    nlist,
                      rvec * gmx_restrict              xx,
                      gmx::ForceWithShiftForces *      forceWithShiftForces,
                      const t_forcerec * gmx_restrict  fr,
                      const t_mdatoms * gmx_restrict   mdatoms,
                      nb_kernel_data_t * gmx_restrict  kernel_data,
                      t_nrnb * gmx_restrict            nrnb)
{
    using SCReal = typename SoftCoreReal<softCoreTreatment>::Real;

    constexpr bool useSoftCore = (softCoreTreatment != SoftCoreTreatment::None);

#define  STATE_A  0
#define  STATE_B  1
#define  NSTATES  2
    int           i, n, ii, is3, ii3, k, nj0, nj1, jnr, j3, ggid;
    real          shX, shY, shZ;
    real          tx, ty, tz, Fscal;
    SCReal        FscalC[NSTATES], FscalV[NSTATES];  /* Needs double for sc_power==48 */
    real          Vcoul[NSTATES], Vvdw[NSTATES];
    real          rinv6, r;
    real          iqA, iqB;
    real          qq[NSTATES], vctot;
    int           ntiA, ntiB, tj[NSTATES];
    real          Vvdw6, Vvdw12, vvtot;
    real          ix, iy, iz, fix, fiy, fiz;
    real          dx, dy, dz, rsq, rinv;
    real          c6[NSTATES], c12[NSTATES], c6grid;
    real          LFC[NSTATES], LFV[NSTATES], DLF[NSTATES];
    SCReal        dvdl_coul, dvdl_vdw;
    real          lfac_coul[NSTATES], dlfac_coul[NSTATES], lfac_vdw[NSTATES], dlfac_vdw[NSTATES];
    real          sigma6[NSTATES], alpha_vdw_eff, alpha_coul_eff;
    real          rinvC, rinvV;
    SCReal        rp, rpm2, rC, rV, rpinvC, rpinvV; /* Needs double for sc_power==48 */
    real          sigma_pow[NSTATES];
    real          VV, FF;
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
    const real *  chargeA;
    const real *  chargeB;
    real          sigma6_min, sigma6_def, lam_power;
    real          alpha_coul, alpha_vdw, lambda_coul, lambda_vdw;
    const real *  nbfp, *nbfp_grid;
    real *        dvdl;
    real *        Vv;
    real *        Vc;
    gmx_bool      bDoForces, bDoShiftForces, bDoPotential;
    gmx_bool      bEwald, bEwaldLJ;
    real          rcutoff_max2;
    const real *  tab_ewald_F_lj = nullptr;
    const real *  tab_ewald_V_lj = nullptr;
    real          d, d2, sw, dsw;
    real          vdw_swV3, vdw_swV4, vdw_swV5, vdw_swF2, vdw_swF3, vdw_swF4;
    gmx_bool      bComputeVdwInteraction, bComputeElecInteraction;
    const real *  ewtab = nullptr;
    int           ewitab;
    real          ewrt, eweps, ewtabscale = 0, ewtabhalfspace = 0, sh_ewald = 0;

    const real    onetwelfth  = 1.0/12.0;
    const real    onesixth    = 1.0/6.0;
    const real    zero        = 0.0;
    const real    half        = 0.5;
    const real    one         = 1.0;
    const real    two         = 2.0;
    const real    six         = 6.0;

    /* Extract pointer to non-bonded interaction constants */
    const interaction_const_t *ic = fr->ic;

    // TODO: We should get rid of using pointers to real
    const real          *x      = xx[0];
    real * gmx_restrict  f      = &(forceWithShiftForces->force()[0][0]);

    real * gmx_restrict  fshift = &(forceWithShiftForces->shiftForces()[0][0]);

    // Extract pair list data
    nri                 = nlist->nri;
    iinr                = nlist->iinr;
    jindex              = nlist->jindex;
    jjnr                = nlist->jjnr;
    shift               = nlist->shift;
    gid                 = nlist->gid;

    shiftvec            = fr->shift_vec[0];
    chargeA             = mdatoms->chargeA;
    chargeB             = mdatoms->chargeB;
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
    sigma6_def          = fr->sc_sigma6_def;
    sigma6_min          = fr->sc_sigma6_min;
    bDoForces           = ((kernel_data->flags & GMX_NONBONDED_DO_FORCE) != 0);
    bDoShiftForces      = ((kernel_data->flags & GMX_NONBONDED_DO_SHIFTFORCE) != 0);
    bDoPotential        = ((kernel_data->flags & GMX_NONBONDED_DO_POTENTIAL) != 0);

    // Extract data from interaction_const_t
    const real facel           = ic->epsfac;
    const real rcoulomb        = ic->rcoulomb;
    const real krf             = ic->k_rf;
    const real crf             = ic->c_rf;
    const real sh_lj_ewald     = ic->sh_lj_ewald;
    const real rvdw            = ic->rvdw;
    const real dispersionShift = ic->dispersion_shift.cpot;
    const real repulsionShift  = ic->repulsion_shift.cpot;

    // Note that the nbnxm kernels do not support Coulomb potential switching at all
    GMX_ASSERT(ic->coulomb_modifier != eintmodPOTSWITCH,
               "Potential switching is not supported for Coulomb with FEP");

    if (ic->vdw_modifier == eintmodPOTSWITCH)
    {
        d               = ic->rvdw - ic->rvdw_switch;
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

    rcutoff_max2    = std::max(ic->rcoulomb, ic->rvdw);
    rcutoff_max2    = rcutoff_max2*rcutoff_max2;

    bEwald          = (icoul == GMX_NBKERNEL_ELEC_EWALD);
    bEwaldLJ        = (ivdw == GMX_NBKERNEL_VDW_LJEWALD);

    if (bEwald || bEwaldLJ)
    {
        const auto &tables = *ic->coulombEwaldTables;
        sh_ewald       = ic->sh_ewald;
        ewtab          = tables.tableFDV0.data();
        ewtabscale     = tables.scale;
        ewtabhalfspace = half/ewtabscale;
        tab_ewald_F_lj = tables.tableF.data();
        tab_ewald_V_lj = tables.tableV.data();
    }

    /* For Ewald/PME interactions we cannot easily apply the soft-core component to
     * reciprocal space. When we use non-switched Ewald interactions, we
     * assume the soft-coring does not significantly affect the grid contribution
     * and apply the soft-core only to the full 1/r (- shift) pair contribution.
     *
     * However, we cannot use this approach for switch-modified since we would then
     * effectively end up evaluating a significantly different interaction here compared to the
     * normal (non-free-energy) kernels, either by applying a cutoff at a different
     * position than what the user requested, or by switching different
     * things (1/r rather than short-range Ewald). For these settings, we just
     * use the traditional short-range Ewald interaction in that case.
     */
    GMX_RELEASE_ASSERT(!(bEwald   && ic->coulomb_modifier == eintmodPOTSWITCH) &&
                       !(bEwaldLJ && ic->vdw_modifier   == eintmodPOTSWITCH),
                       "Can not apply soft-core to switched Ewald potentials");

    dvdl_coul  = 0;
    dvdl_vdw   = 0;

    /* Lambda factor for state A, 1-lambda*/
    LFC[STATE_A] = one - lambda_coul;
    LFV[STATE_A] = one - lambda_vdw;

    /* Lambda factor for state B, lambda*/
    LFC[STATE_B] = lambda_coul;
    LFV[STATE_B] = lambda_vdw;

    /*derivative of the lambda factor for state A and B */
    DLF[STATE_A] = -1;
    DLF[STATE_B] = 1;

    constexpr real sc_r_power = (softCoreTreatment == SoftCoreTreatment::RPower48 ? 48.0_real : 6.0_real);
    for (i = 0; i < NSTATES; i++)
    {
        lfac_coul[i]  = (lam_power == 2 ? (1-LFC[i])*(1-LFC[i]) : (1-LFC[i]));
        dlfac_coul[i] = DLF[i]*lam_power/sc_r_power*(lam_power == 2 ? (1-LFC[i]) : 1);
        lfac_vdw[i]   = (lam_power == 2 ? (1-LFV[i])*(1-LFV[i]) : (1-LFV[i]));
        dlfac_vdw[i]  = DLF[i]*lam_power/sc_r_power*(lam_power == 2 ? (1-LFV[i]) : 1);
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

            if (rsq >= rcutoff_max2)
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
                /* Note that unlike in the nbnxn kernels, we do not need
                 * to clamp the value of rsq before taking the invsqrt
                 * to avoid NaN in the LJ calculation, since here we do
                 * not calculate LJ interactions when C6 and C12 are zero.
                 */

                rinv         = gmx::invsqrt(rsq);
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

            if (softCoreTreatment == SoftCoreTreatment::None)
            {
                /* The soft-core power p will not affect the results
                 * with not using soft-core, so we use power of 0 which gives
                 * the simplest math and cheapest code.
                 */
                rpm2         = rinv*rinv;
                rp           = 1;
            }
            if (softCoreTreatment == SoftCoreTreatment::RPower6)
            {
                rpm2         = rsq*rsq;  /* r4 */
                rp           = rpm2*rsq; /* r6 */
            }
            if (softCoreTreatment == SoftCoreTreatment::RPower48)
            {
                rp           = rsq*rsq*rsq; /* r6 */
                rp           = rp*rp;       /* r12 */
                rp           = rp*rp;       /* r24 */
                rp           = rp*rp;       /* r48 */
                rpm2         = rp/rsq;      /* r46 */
            }

            Fscal = 0;

            qq[STATE_A]      = iqA*chargeA[jnr];
            qq[STATE_B]      = iqB*chargeB[jnr];

            tj[STATE_A]      = ntiA+2*typeA[jnr];
            tj[STATE_B]      = ntiB+2*typeB[jnr];

            if (nlist->excl_fep == nullptr || nlist->excl_fep[k])
            {
                c6[STATE_A]      = nbfp[tj[STATE_A]];
                c6[STATE_B]      = nbfp[tj[STATE_B]];

                for (i = 0; i < NSTATES; i++)
                {
                    c12[i]             = nbfp[tj[i]+1];
                    if (useSoftCore)
                    {
                        if ((c6[i] > 0) && (c12[i] > 0))
                        {
                            /* c12 is stored scaled with 12.0 and c6 is scaled with 6.0 - correct for this */
                            sigma6[i]       = half*c12[i]/c6[i];
                            if (sigma6[i] < sigma6_min)   /* for disappearing coul and vdw with soft core at the same time */
                            {
                                sigma6[i] = sigma6_min;
                            }
                        }
                        else
                        {
                            sigma6[i]       = sigma6_def;
                        }
                        if (softCoreTreatment == SoftCoreTreatment::RPower6)
                        {
                            sigma_pow[i]    = sigma6[i];
                        }
                        else
                        {
                            sigma_pow[i]    = sigma6[i]*sigma6[i];       /* sigma^12 */
                            sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^24 */
                            sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^48 */
                        }
                    }
                }

                if (useSoftCore)
                {
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
                        if (useSoftCore)
                        {
                            rpinvC     = one/(alpha_coul_eff*lfac_coul[i]*sigma_pow[i]+rp);
                            pthRoot<softCoreTreatment>(rpinvC, &rinvC, &rC);

                            if (scLambdasOrAlphasDiffer)
                            {
                                rpinvV = one/(alpha_vdw_eff*lfac_vdw[i]*sigma_pow[i]+rp);
                                pthRoot<softCoreTreatment>(rpinvV, &rinvV, &rV);
                            }
                            else
                            {
                                /* We can avoid one expensive pow and one / operation */
                                rpinvV = rpinvC;
                                rinvV  = rinvC;
                                rV     = rC;
                            }
                        }
                        else
                        {
                            rpinvC     = 1;
                            rinvC      = rinv;
                            rC         = r;

                            rpinvV     = 1;
                            rinvV      = rinv;
                            rV         = r;
                        }

                        /* Only process the coulomb interactions if we have charges,
                         * and if we either include all entries in the list (no cutoff
                         * used in the kernel), or if we are within the cutoff.
                         */
                        bComputeElecInteraction =
                            ( bEwald && r < rcoulomb) ||
                            (!bEwald && rC < rcoulomb);

                        if ( (qq[i] != 0) && bComputeElecInteraction)
                        {
                            if (bEwald)
                            {
                                /* Ewald FEP is done only on the 1/r part */
                                Vcoul[i]   = qq[i]*(rinvC-sh_ewald);
                                FscalC[i]  = qq[i]*rinvC;
                            }
                            else
                            {
                                /* reaction-field */
                                Vcoul[i]   = qq[i]*(rinvC + krf*rC*rC-crf);
                                FscalC[i]  = qq[i]*(rinvC - two*krf*rC*rC);
                            }
                        }

                        /* Only process the VDW interactions if we have
                         * some non-zero parameters, and if we either
                         * include all entries in the list (no cutoff used
                         * in the kernel), or if we are within the cutoff.
                         */
                        bComputeVdwInteraction =
                            ( bEwaldLJ && r < rvdw) ||
                            (!bEwaldLJ && rV < rvdw);
                        if ((c6[i] != 0 || c12[i] != 0) && bComputeVdwInteraction)
                        {
                            /* cutoff LJ, also handles part of Ewald LJ */
                            if (softCoreTreatment == SoftCoreTreatment::RPower6)
                            {
                                rinv6        = rpinvV;
                            }
                            else
                            {
                                rinv6        = rinvV*rinvV;
                                rinv6        = rinv6*rinv6*rinv6;
                            }
                            Vvdw6            = c6[i]*rinv6;
                            Vvdw12           = c12[i]*rinv6*rinv6;

                            Vvdw[i]          = ( (Vvdw12 + c12[i]*repulsionShift)*onetwelfth
                                                 - (Vvdw6 + c6[i]*dispersionShift)*onesixth);
                            FscalV[i]        = Vvdw12 - Vvdw6;

                            if (bEwaldLJ)
                            {
                                /* Subtract the grid potential at the cut-off */
                                c6grid      = nbfp_grid[tj[i]];
                                Vvdw[i]    += c6grid*sh_lj_ewald*onesixth;
                            }

                            if (ic->vdw_modifier == eintmodPOTSWITCH)
                            {
                                d                = rV - ic->rvdw_switch;
                                d                = (d > zero) ? d : zero;
                                d2               = d*d;
                                sw               = one+d2*d*(vdw_swV3+d*(vdw_swV4+d*vdw_swV5));
                                dsw              = d2*(vdw_swF2+d*(vdw_swF3+d*vdw_swF4));

                                FscalV[i]        = FscalV[i]*sw - rV*Vvdw[i]*dsw;
                                Vvdw[i]         *= sw;

                                FscalV[i]  = (rV < rvdw) ? FscalV[i] : zero;
                                Vvdw[i]    = (rV < rvdw) ? Vvdw[i] : zero;
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

                    Fscal     += LFC[i]*FscalC[i]*rpm2;
                    Fscal     += LFV[i]*FscalV[i]*rpm2;

                    if (useSoftCore)
                    {
                        dvdl_coul += Vcoul[i]*DLF[i] + LFC[i]*alpha_coul_eff*dlfac_coul[i]*FscalC[i]*sigma_pow[i];
                        dvdl_vdw  += Vvdw[i]*DLF[i] + LFV[i]*alpha_vdw_eff*dlfac_vdw[i]*FscalV[i]*sigma_pow[i];
                    }
                    else
                    {
                        dvdl_coul += Vcoul[i]*DLF[i];
                        dvdl_vdw  += Vvdw[i]*DLF[i];
                    }
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
                FF = -two*krf;

                if (ii == jnr)
                {
                    VV *= half;
                }

                for (i = 0; i < NSTATES; i++)
                {
                    vctot      += LFC[i]*qq[i]*VV;
                    Fscal      += LFC[i]*qq[i]*FF;
                    dvdl_coul  += DLF[i]*qq[i]*VV;
                }
            }

            if (bEwald && r < rcoulomb)
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
                ewitab    = static_cast<int>(ewrt);
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
                    v_lr *= half;
                }

                for (i = 0; i < NSTATES; i++)
                {
                    vctot      -= LFC[i]*qq[i]*v_lr;
                    Fscal      -= LFC[i]*qq[i]*f_lr;
                    dvdl_coul  -= (DLF[i]*qq[i])*v_lr;
                }
            }

            if (bEwaldLJ && r < rvdw)
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
                ri     = static_cast<int>(rs);
                frac   = rs - ri;
                f_lr   = (1 - frac)*tab_ewald_F_lj[ri] + frac*tab_ewald_F_lj[ri+1];
                /* TODO: Currently the Ewald LJ table does not contain
                 * the factor 1/6, we should add this.
                 */
                FF     = f_lr*rinv/six;
                VV     = (tab_ewald_V_lj[ri] - ewtabhalfspace*frac*(tab_ewald_F_lj[ri] + f_lr))/six;

                if (ii == jnr)
                {
                    /* If we get here, the i particle (ii) has itself (jnr)
                     * in its neighborlist. This can only happen with the Verlet
                     * scheme, and corresponds to a self-interaction that will
                     * occur twice. Scale it down by 50% to only include it once.
                     */
                    VV *= half;
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
                 *
                 * All the OpenMP regions in this file are trivial and should
                 * not throw, so no need for try/catch.
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

void gmx_nb_free_energy_kernel(const t_nblist            *nlist,
                               rvec                      *xx,
                               gmx::ForceWithShiftForces *ff,
                               const t_forcerec          *fr,
                               const t_mdatoms           *mdatoms,
                               nb_kernel_data_t          *kernel_data,
                               t_nrnb                    *nrnb)
{
    if (fr->sc_alphacoul == 0 && fr->sc_alphavdw == 0)
    {
        nb_free_energy_kernel<SoftCoreTreatment::None, false>(nlist, xx, ff, fr, mdatoms, kernel_data, nrnb);
    }
    else if (fr->sc_r_power == 6.0_real)
    {
        if (kernel_data->lambda[efptCOUL] == kernel_data->lambda[efptVDW] &&
            fr->sc_alphacoul == fr->sc_alphavdw)
        {
            nb_free_energy_kernel<SoftCoreTreatment::RPower6, false>(nlist, xx, ff, fr, mdatoms, kernel_data, nrnb);
        }
        else
        {
            nb_free_energy_kernel<SoftCoreTreatment::RPower6, true>(nlist, xx, ff, fr, mdatoms, kernel_data, nrnb);
        }
    }
    else if (fr->sc_r_power == 48.0_real)
    {
        nb_free_energy_kernel<SoftCoreTreatment::RPower48, true>(nlist, xx, ff, fr, mdatoms, kernel_data, nrnb);
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Unsupported soft-core r-power");
    }
}
