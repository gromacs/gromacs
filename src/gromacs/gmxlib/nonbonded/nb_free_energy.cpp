/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"


//! Scalar (non-SIMD) data types.
struct ScalarDataTypes
{
    using RealType                     = real; //!< The data type to use as real.
    using IntType                      = int;  //!< The data type to use as int.
    static constexpr int simdRealWidth = 1;    //!< The width of the RealType.
    static constexpr int simdIntWidth  = 1;    //!< The width of the IntType.
};

#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32_ARITHMETICS
//! SIMD data types.
struct SimdDataTypes
{
    using RealType                     = gmx::SimdReal;         //!< The data type to use as real.
    using IntType                      = gmx::SimdInt32;        //!< The data type to use as int.
    static constexpr int simdRealWidth = GMX_SIMD_REAL_WIDTH;   //!< The width of the RealType.
    static constexpr int simdIntWidth  = GMX_SIMD_FINT32_WIDTH; //!< The width of the IntType.
};
#endif

//! Computes r^(1/p) and 1/r^(1/p) for the standard p=6
template<class RealType>
static inline void pthRoot(const RealType r, RealType* pthRoot, RealType* invPthRoot)
{
    *invPthRoot = gmx::invsqrt(std::cbrt(r));
    *pthRoot    = 1 / (*invPthRoot);
}

template<class RealType>
static inline RealType calculateRinv6(const RealType rinvV)
{
    RealType rinv6 = rinvV * rinvV;
    return (rinv6 * rinv6 * rinv6);
}

template<class RealType>
static inline RealType calculateVdw6(const RealType c6, const RealType rinv6)
{
    return (c6 * rinv6);
}

template<class RealType>
static inline RealType calculateVdw12(const RealType c12, const RealType rinv6)
{
    return (c12 * rinv6 * rinv6);
}

/* reaction-field electrostatics */
template<class RealType>
static inline RealType reactionFieldScalarForce(const RealType qq,
                                                const RealType rinv,
                                                const RealType r,
                                                const real     krf,
                                                const real     two)
{
    return (qq * (rinv - two * krf * r * r));
}
template<class RealType>
static inline RealType reactionFieldPotential(const RealType qq,
                                              const RealType rinv,
                                              const RealType r,
                                              const real     krf,
                                              const real     potentialShift)
{
    return (qq * (rinv + krf * r * r - potentialShift));
}

/* Ewald electrostatics */
template<class RealType>
static inline RealType ewaldScalarForce(const RealType coulomb, const RealType rinv)
{
    return (coulomb * rinv);
}
template<class RealType>
static inline RealType ewaldPotential(const RealType coulomb, const RealType rinv, const real potentialShift)
{
    return (coulomb * (rinv - potentialShift));
}

/* cutoff LJ */
template<class RealType>
static inline RealType lennardJonesScalarForce(const RealType v6, const RealType v12)
{
    return (v12 - v6);
}
template<class RealType>
static inline RealType lennardJonesPotential(const RealType v6,
                                             const RealType v12,
                                             const RealType c6,
                                             const RealType c12,
                                             const real     repulsionShift,
                                             const real     dispersionShift,
                                             const real     onesixth,
                                             const real     onetwelfth)
{
    return ((v12 + c12 * repulsionShift) * onetwelfth - (v6 + c6 * dispersionShift) * onesixth);
}

/* Ewald LJ */
static inline real ewaldLennardJonesGridSubtract(const real c6grid, const real potentialShift, const real onesixth)
{
    return (c6grid * potentialShift * onesixth);
}

/* LJ Potential switch */
template<class RealType>
static inline RealType potSwitchScalarForceMod(const RealType fScalarInp,
                                               const RealType potential,
                                               const RealType sw,
                                               const RealType r,
                                               const RealType rVdw,
                                               const RealType dsw,
                                               const real     zero)
{
    if (r < rVdw)
    {
        real fScalar = fScalarInp * sw - r * potential * dsw;
        return (fScalar);
    }
    return (zero);
}
template<class RealType>
static inline RealType potSwitchPotentialMod(const RealType potentialInp,
                                             const RealType sw,
                                             const RealType r,
                                             const RealType rVdw,
                                             const real     zero)
{
    if (r < rVdw)
    {
        real potential = potentialInp * sw;
        return (potential);
    }
    return (zero);
}


//! Templated free-energy non-bonded kernel
template<typename DataTypes, bool useSoftCore, bool scLambdasOrAlphasDiffer, bool vdwInteractionTypeIsEwald, bool elecInteractionTypeIsEwald, bool vdwModifierIsPotSwitch>
static void nb_free_energy_kernel(const t_nblist* gmx_restrict nlist,
                                  rvec* gmx_restrict         xx,
                                  gmx::ForceWithShiftForces* forceWithShiftForces,
                                  const t_forcerec* gmx_restrict fr,
                                  const t_mdatoms* gmx_restrict mdatoms,
                                  nb_kernel_data_t* gmx_restrict kernel_data,
                                  t_nrnb* gmx_restrict nrnb)
{
#define STATE_A 0
#define STATE_B 1
#define NSTATES 2

    using RealType = typename DataTypes::RealType;
    using IntType  = typename DataTypes::IntType;

    /* FIXME: How should these be handled with SIMD? */
    constexpr real onetwelfth = 1.0 / 12.0;
    constexpr real onesixth   = 1.0 / 6.0;
    constexpr real zero       = 0.0;
    constexpr real half       = 0.5;
    constexpr real one        = 1.0;
    constexpr real two        = 2.0;
    constexpr real six        = 6.0;

    /* Extract pointer to non-bonded interaction constants */
    const interaction_const_t* ic = fr->ic;

    // Extract pair list data
    const int  nri    = nlist->nri;
    const int* iinr   = nlist->iinr;
    const int* jindex = nlist->jindex;
    const int* jjnr   = nlist->jjnr;
    const int* shift  = nlist->shift;
    const int* gid    = nlist->gid;

    const real* shiftvec      = fr->shift_vec[0];
    const real* chargeA       = mdatoms->chargeA;
    const real* chargeB       = mdatoms->chargeB;
    real*       Vc            = kernel_data->energygrp_elec;
    const int*  typeA         = mdatoms->typeA;
    const int*  typeB         = mdatoms->typeB;
    const int   ntype         = fr->ntype;
    const real* nbfp          = fr->nbfp.data();
    const real* nbfp_grid     = fr->ljpme_c6grid;
    real*       Vv            = kernel_data->energygrp_vdw;
    const real  lambda_coul   = kernel_data->lambda[efptCOUL];
    const real  lambda_vdw    = kernel_data->lambda[efptVDW];
    real*       dvdl          = kernel_data->dvdl;
    const real  alpha_coul    = fr->sc_alphacoul;
    const real  alpha_vdw     = fr->sc_alphavdw;
    const real  lam_power     = fr->sc_power;
    const real  sigma6_def    = fr->sc_sigma6_def;
    const real  sigma6_min    = fr->sc_sigma6_min;
    const bool  doForces      = ((kernel_data->flags & GMX_NONBONDED_DO_FORCE) != 0);
    const bool  doShiftForces = ((kernel_data->flags & GMX_NONBONDED_DO_SHIFTFORCE) != 0);
    const bool  doPotential   = ((kernel_data->flags & GMX_NONBONDED_DO_POTENTIAL) != 0);

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

    real vdw_swV3, vdw_swV4, vdw_swV5, vdw_swF2, vdw_swF3, vdw_swF4;
    if (vdwModifierIsPotSwitch)
    {
        const real d = ic->rvdw - ic->rvdw_switch;
        vdw_swV3     = -10.0 / (d * d * d);
        vdw_swV4     = 15.0 / (d * d * d * d);
        vdw_swV5     = -6.0 / (d * d * d * d * d);
        vdw_swF2     = -30.0 / (d * d * d);
        vdw_swF3     = 60.0 / (d * d * d * d);
        vdw_swF4     = -30.0 / (d * d * d * d * d);
    }
    else
    {
        /* Avoid warnings from stupid compilers (looking at you, Clang!) */
        vdw_swV3 = vdw_swV4 = vdw_swV5 = vdw_swF2 = vdw_swF3 = vdw_swF4 = 0.0;
    }

    int icoul;
    if (ic->eeltype == eelCUT || EEL_RF(ic->eeltype))
    {
        icoul = GMX_NBKERNEL_ELEC_REACTIONFIELD;
    }
    else
    {
        icoul = GMX_NBKERNEL_ELEC_NONE;
    }

    real rcutoff_max2 = std::max(ic->rcoulomb, ic->rvdw);
    rcutoff_max2      = rcutoff_max2 * rcutoff_max2;

    const real* tab_ewald_F_lj = nullptr;
    const real* tab_ewald_V_lj = nullptr;
    const real* ewtab          = nullptr;
    real        ewtabscale     = 0;
    real        ewtabhalfspace = 0;
    real        sh_ewald       = 0;
    if (elecInteractionTypeIsEwald || vdwInteractionTypeIsEwald)
    {
        const auto& tables = *ic->coulombEwaldTables;
        sh_ewald           = ic->sh_ewald;
        ewtab              = tables.tableFDV0.data();
        ewtabscale         = tables.scale;
        ewtabhalfspace     = half / ewtabscale;
        tab_ewald_F_lj     = tables.tableF.data();
        tab_ewald_V_lj     = tables.tableV.data();
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
    GMX_RELEASE_ASSERT(!(vdwInteractionTypeIsEwald && vdwModifierIsPotSwitch),
                       "Can not apply soft-core to switched Ewald potentials");

    real dvdl_coul = 0;
    real dvdl_vdw  = 0;

    /* Lambda factor for state A, 1-lambda*/
    real LFC[NSTATES], LFV[NSTATES];
    LFC[STATE_A] = one - lambda_coul;
    LFV[STATE_A] = one - lambda_vdw;

    /* Lambda factor for state B, lambda*/
    LFC[STATE_B] = lambda_coul;
    LFV[STATE_B] = lambda_vdw;

    /*derivative of the lambda factor for state A and B */
    real DLF[NSTATES];
    DLF[STATE_A] = -1;
    DLF[STATE_B] = 1;

    real           lfac_coul[NSTATES], dlfac_coul[NSTATES], lfac_vdw[NSTATES], dlfac_vdw[NSTATES];
    constexpr real sc_r_power = 6.0_real;
    for (int i = 0; i < NSTATES; i++)
    {
        lfac_coul[i]  = (lam_power == 2 ? (1 - LFC[i]) * (1 - LFC[i]) : (1 - LFC[i]));
        dlfac_coul[i] = DLF[i] * lam_power / sc_r_power * (lam_power == 2 ? (1 - LFC[i]) : 1);
        lfac_vdw[i]   = (lam_power == 2 ? (1 - LFV[i]) * (1 - LFV[i]) : (1 - LFV[i]));
        dlfac_vdw[i]  = DLF[i] * lam_power / sc_r_power * (lam_power == 2 ? (1 - LFV[i]) : 1);
    }

    // TODO: We should get rid of using pointers to real
    const real* x             = xx[0];
    real* gmx_restrict f      = &(forceWithShiftForces->force()[0][0]);
    real* gmx_restrict fshift = &(forceWithShiftForces->shiftForces()[0][0]);

    for (int n = 0; n < nri; n++)
    {
        int npair_within_cutoff = 0;

        const int  is3   = 3 * shift[n];
        const real shX   = shiftvec[is3];
        const real shY   = shiftvec[is3 + 1];
        const real shZ   = shiftvec[is3 + 2];
        const int  nj0   = jindex[n];
        const int  nj1   = jindex[n + 1];
        const int  ii    = iinr[n];
        const int  ii3   = 3 * ii;
        const real ix    = shX + x[ii3 + 0];
        const real iy    = shY + x[ii3 + 1];
        const real iz    = shZ + x[ii3 + 2];
        const real iqA   = facel * chargeA[ii];
        const real iqB   = facel * chargeB[ii];
        const int  ntiA  = 2 * ntype * typeA[ii];
        const int  ntiB  = 2 * ntype * typeB[ii];
        real       vctot = 0;
        real       vvtot = 0;
        real       fix   = 0;
        real       fiy   = 0;
        real       fiz   = 0;

        for (int k = nj0; k < nj1; k++)
        {
            int            tj[NSTATES];
            const int      jnr = jjnr[k];
            const int      j3  = 3 * jnr;
            RealType       c6[NSTATES], c12[NSTATES], qq[NSTATES], Vcoul[NSTATES], Vvdw[NSTATES];
            RealType       r, rinv, rp, rpm2;
            RealType       alpha_vdw_eff, alpha_coul_eff, sigma6[NSTATES];
            const RealType dx  = ix - x[j3];
            const RealType dy  = iy - x[j3 + 1];
            const RealType dz  = iz - x[j3 + 2];
            const RealType rsq = dx * dx + dy * dy + dz * dz;
            RealType       FscalC[NSTATES], FscalV[NSTATES];

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

                rinv = gmx::invsqrt(rsq);
                r    = rsq * rinv;
            }
            else
            {
                /* The force at r=0 is zero, because of symmetry.
                 * But note that the potential is in general non-zero,
                 * since the soft-cored r will be non-zero.
                 */
                rinv = 0;
                r    = 0;
            }

            if (useSoftCore)
            {
                rpm2 = rsq * rsq;  /* r4 */
                rp   = rpm2 * rsq; /* r6 */
            }
            else
            {
                /* The soft-core power p will not affect the results
                 * with not using soft-core, so we use power of 0 which gives
                 * the simplest math and cheapest code.
                 */
                rpm2 = rinv * rinv;
                rp   = 1;
            }

            RealType Fscal = 0;

            qq[STATE_A] = iqA * chargeA[jnr];
            qq[STATE_B] = iqB * chargeB[jnr];

            tj[STATE_A] = ntiA + 2 * typeA[jnr];
            tj[STATE_B] = ntiB + 2 * typeB[jnr];

            if (nlist->excl_fep == nullptr || nlist->excl_fep[k])
            {
                c6[STATE_A] = nbfp[tj[STATE_A]];
                c6[STATE_B] = nbfp[tj[STATE_B]];

                for (int i = 0; i < NSTATES; i++)
                {
                    c12[i] = nbfp[tj[i] + 1];
                    if (useSoftCore)
                    {
                        if ((c6[i] > 0) && (c12[i] > 0))
                        {
                            /* c12 is stored scaled with 12.0 and c6 is scaled with 6.0 - correct for this */
                            sigma6[i] = half * c12[i] / c6[i];
                            if (sigma6[i] < sigma6_min) /* for disappearing coul and vdw with soft core at the same time */
                            {
                                sigma6[i] = sigma6_min;
                            }
                        }
                        else
                        {
                            sigma6[i] = sigma6_def;
                        }
                    }
                }

                if (useSoftCore)
                {
                    /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
                    if ((c12[STATE_A] > 0) && (c12[STATE_B] > 0))
                    {
                        alpha_vdw_eff  = 0;
                        alpha_coul_eff = 0;
                    }
                    else
                    {
                        alpha_vdw_eff  = alpha_vdw;
                        alpha_coul_eff = alpha_coul;
                    }
                }

                for (int i = 0; i < NSTATES; i++)
                {
                    FscalC[i] = 0;
                    FscalV[i] = 0;
                    Vcoul[i]  = 0;
                    Vvdw[i]   = 0;

                    RealType rinvC, rinvV, rC, rV, rpinvC, rpinvV;

                    /* Only spend time on A or B state if it is non-zero */
                    if ((qq[i] != 0) || (c6[i] != 0) || (c12[i] != 0))
                    {
                        /* this section has to be inside the loop because of the dependence on sigma6 */
                        if (useSoftCore)
                        {
                            rpinvC = one / (alpha_coul_eff * lfac_coul[i] * sigma6[i] + rp);
                            pthRoot(rpinvC, &rinvC, &rC);
                            if (scLambdasOrAlphasDiffer)
                            {
                                rpinvV = one / (alpha_vdw_eff * lfac_vdw[i] * sigma6[i] + rp);
                                pthRoot(rpinvV, &rinvV, &rV);
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
                            rpinvC = 1;
                            rinvC  = rinv;
                            rC     = r;

                            rpinvV = 1;
                            rinvV  = rinv;
                            rV     = r;
                        }

                        /* Only process the coulomb interactions if we have charges,
                         * and if we either include all entries in the list (no cutoff
                         * used in the kernel), or if we are within the cutoff.
                         */
                        bool computeElecInteraction = (elecInteractionTypeIsEwald && r < rcoulomb)
                                                      || (!elecInteractionTypeIsEwald && rC < rcoulomb);

                        if ((qq[i] != 0) && computeElecInteraction)
                        {
                            if (elecInteractionTypeIsEwald)
                            {
                                Vcoul[i]  = ewaldPotential(qq[i], rinvC, sh_ewald);
                                FscalC[i] = ewaldScalarForce(qq[i], rinvC);
                            }
                            else
                            {
                                Vcoul[i]  = reactionFieldPotential(qq[i], rinvC, rC, krf, crf);
                                FscalC[i] = reactionFieldScalarForce(qq[i], rinvC, rC, krf, two);
                            }
                        }

                        /* Only process the VDW interactions if we have
                         * some non-zero parameters, and if we either
                         * include all entries in the list (no cutoff used
                         * in the kernel), or if we are within the cutoff.
                         */
                        bool computeVdwInteraction = (vdwInteractionTypeIsEwald && r < rvdw)
                                                     || (!vdwInteractionTypeIsEwald && rV < rvdw);
                        if ((c6[i] != 0 || c12[i] != 0) && computeVdwInteraction)
                        {
                            RealType rinv6;
                            if (useSoftCore)
                            {
                                rinv6 = rpinvV;
                            }
                            else
                            {
                                rinv6 = calculateRinv6(rinvV);
                            }
                            RealType Vvdw6  = calculateVdw6(c6[i], rinv6);
                            RealType Vvdw12 = calculateVdw12(c12[i], rinv6);

                            Vvdw[i] = lennardJonesPotential(Vvdw6, Vvdw12, c6[i], c12[i], repulsionShift,
                                                            dispersionShift, onesixth, onetwelfth);
                            FscalV[i] = lennardJonesScalarForce(Vvdw6, Vvdw12);

                            if (vdwInteractionTypeIsEwald)
                            {
                                /* Subtract the grid potential at the cut-off */
                                Vvdw[i] += ewaldLennardJonesGridSubtract(nbfp_grid[tj[i]],
                                                                         sh_lj_ewald, onesixth);
                            }

                            if (vdwModifierIsPotSwitch)
                            {
                                RealType d        = rV - ic->rvdw_switch;
                                d                 = (d > zero) ? d : zero;
                                const RealType d2 = d * d;
                                const RealType sw =
                                        one + d2 * d * (vdw_swV3 + d * (vdw_swV4 + d * vdw_swV5));
                                const RealType dsw = d2 * (vdw_swF2 + d * (vdw_swF3 + d * vdw_swF4));

                                FscalV[i] = potSwitchScalarForceMod(FscalV[i], Vvdw[i], sw, rV,
                                                                    rvdw, dsw, zero);
                                Vvdw[i]   = potSwitchPotentialMod(Vvdw[i], sw, rV, rvdw, zero);
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
                for (int i = 0; i < NSTATES; i++)
                {
                    vctot += LFC[i] * Vcoul[i];
                    vvtot += LFV[i] * Vvdw[i];

                    Fscal += LFC[i] * FscalC[i] * rpm2;
                    Fscal += LFV[i] * FscalV[i] * rpm2;

                    if (useSoftCore)
                    {
                        dvdl_coul += Vcoul[i] * DLF[i]
                                     + LFC[i] * alpha_coul_eff * dlfac_coul[i] * FscalC[i] * sigma6[i];
                        dvdl_vdw += Vvdw[i] * DLF[i]
                                    + LFV[i] * alpha_vdw_eff * dlfac_vdw[i] * FscalV[i] * sigma6[i];
                    }
                    else
                    {
                        dvdl_coul += Vcoul[i] * DLF[i];
                        dvdl_vdw += Vvdw[i] * DLF[i];
                    }
                }
            }
            else if (icoul == GMX_NBKERNEL_ELEC_REACTIONFIELD)
            {
                /* For excluded pairs, which are only in this pair list when
                 * using the Verlet scheme, we don't use soft-core.
                 * As there is no singularity, there is no need for soft-core.
                 */
                const real FF = -two * krf;
                RealType   VV = krf * rsq - crf;

                if (ii == jnr)
                {
                    VV *= half;
                }

                for (int i = 0; i < NSTATES; i++)
                {
                    vctot += LFC[i] * qq[i] * VV;
                    Fscal += LFC[i] * qq[i] * FF;
                    dvdl_coul += DLF[i] * qq[i] * VV;
                }
            }

            if (elecInteractionTypeIsEwald && r < rcoulomb)
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

                const RealType ewrt   = r * ewtabscale;
                IntType        ewitab = static_cast<IntType>(ewrt);
                const RealType eweps  = ewrt - ewitab;
                ewitab                = 4 * ewitab;
                f_lr                  = ewtab[ewitab] + eweps * ewtab[ewitab + 1];
                v_lr = (ewtab[ewitab + 2] - ewtabhalfspace * eweps * (ewtab[ewitab] + f_lr));
                f_lr *= rinv;

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

                for (int i = 0; i < NSTATES; i++)
                {
                    vctot -= LFC[i] * qq[i] * v_lr;
                    Fscal -= LFC[i] * qq[i] * f_lr;
                    dvdl_coul -= (DLF[i] * qq[i]) * v_lr;
                }
            }

            if (vdwInteractionTypeIsEwald && r < rvdw)
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

                const RealType rs   = rsq * rinv * ewtabscale;
                const IntType  ri   = static_cast<IntType>(rs);
                const RealType frac = rs - ri;
                const RealType f_lr = (1 - frac) * tab_ewald_F_lj[ri] + frac * tab_ewald_F_lj[ri + 1];
                /* TODO: Currently the Ewald LJ table does not contain
                 * the factor 1/6, we should add this.
                 */
                const RealType FF = f_lr * rinv / six;
                RealType       VV =
                        (tab_ewald_V_lj[ri] - ewtabhalfspace * frac * (tab_ewald_F_lj[ri] + f_lr)) / six;

                if (ii == jnr)
                {
                    /* If we get here, the i particle (ii) has itself (jnr)
                     * in its neighborlist. This can only happen with the Verlet
                     * scheme, and corresponds to a self-interaction that will
                     * occur twice. Scale it down by 50% to only include it once.
                     */
                    VV *= half;
                }

                for (int i = 0; i < NSTATES; i++)
                {
                    const real c6grid = nbfp_grid[tj[i]];
                    vvtot += LFV[i] * c6grid * VV;
                    Fscal += LFV[i] * c6grid * FF;
                    dvdl_vdw += (DLF[i] * c6grid) * VV;
                }
            }

            if (doForces)
            {
                const real tx = Fscal * dx;
                const real ty = Fscal * dy;
                const real tz = Fscal * dz;
                fix           = fix + tx;
                fiy           = fiy + ty;
                fiz           = fiz + tz;
                /* OpenMP atomics are expensive, but this kernels is also
                 * expensive, so we can take this hit, instead of using
                 * thread-local output buffers and extra reduction.
                 *
                 * All the OpenMP regions in this file are trivial and should
                 * not throw, so no need for try/catch.
                 */
#pragma omp atomic
                f[j3] -= tx;
#pragma omp atomic
                f[j3 + 1] -= ty;
#pragma omp atomic
                f[j3 + 2] -= tz;
            }
        }

        /* The atomics below are expensive with many OpenMP threads.
         * Here unperturbed i-particles will usually only have a few
         * (perturbed) j-particles in the list. Thus with a buffered list
         * we can skip a significant number of i-reductions with a check.
         */
        if (npair_within_cutoff > 0)
        {
            if (doForces)
            {
#pragma omp atomic
                f[ii3] += fix;
#pragma omp atomic
                f[ii3 + 1] += fiy;
#pragma omp atomic
                f[ii3 + 2] += fiz;
            }
            if (doShiftForces)
            {
#pragma omp atomic
                fshift[is3] += fix;
#pragma omp atomic
                fshift[is3 + 1] += fiy;
#pragma omp atomic
                fshift[is3 + 2] += fiz;
            }
            if (doPotential)
            {
                int ggid = gid[n];
#pragma omp atomic
                Vc[ggid] += vctot;
#pragma omp atomic
                Vv[ggid] += vvtot;
            }
        }
    }

#pragma omp atomic
    dvdl[efptCOUL] += dvdl_coul;
#pragma omp atomic
    dvdl[efptVDW] += dvdl_vdw;

    /* Estimate flops, average for free energy stuff:
     * 12  flops per outer iteration
     * 150 flops per inner iteration
     */
#pragma omp atomic
    inc_nrnb(nrnb, eNR_NBKERNEL_FREE_ENERGY, nlist->nri * 12 + nlist->jindex[nri] * 150);
}

typedef void (*KernelFunction)(const t_nblist* gmx_restrict nlist,
                               rvec* gmx_restrict         xx,
                               gmx::ForceWithShiftForces* forceWithShiftForces,
                               const t_forcerec* gmx_restrict fr,
                               const t_mdatoms* gmx_restrict mdatoms,
                               nb_kernel_data_t* gmx_restrict kernel_data,
                               t_nrnb* gmx_restrict nrnb);

template<bool useSoftCore, bool scLambdasOrAlphasDiffer, bool vdwInteractionTypeIsEwald, bool elecInteractionTypeIsEwald, bool vdwModifierIsPotSwitch>
static KernelFunction dispatchKernelOnUseSimd(const bool useSimd)
{
    if (useSimd)
    {
#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32_ARITHMETICS
        /* FIXME: Here SimdDataTypes should be used to enable SIMD. So far, the code in nb_free_energy_kernel is not adapted to SIMD */
        return (nb_free_energy_kernel<ScalarDataTypes, useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald,
                                      elecInteractionTypeIsEwald, vdwModifierIsPotSwitch>);
#else
        return (nb_free_energy_kernel<ScalarDataTypes, useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald,
                                      elecInteractionTypeIsEwald, vdwModifierIsPotSwitch>);
#endif
    }
    else
    {
        return (nb_free_energy_kernel<ScalarDataTypes, useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald,
                                      elecInteractionTypeIsEwald, vdwModifierIsPotSwitch>);
    }
}

template<bool useSoftCore, bool scLambdasOrAlphasDiffer, bool vdwInteractionTypeIsEwald, bool elecInteractionTypeIsEwald>
static KernelFunction dispatchKernelOnVdwModifier(const bool vdwModifierIsPotSwitch, const bool useSimd)
{
    if (vdwModifierIsPotSwitch)
    {
        return (dispatchKernelOnUseSimd<useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald,
                                        elecInteractionTypeIsEwald, true>(useSimd));
    }
    else
    {
        return (dispatchKernelOnUseSimd<useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald,
                                        elecInteractionTypeIsEwald, false>(useSimd));
    }
}

template<bool useSoftCore, bool scLambdasOrAlphasDiffer, bool vdwInteractionTypeIsEwald>
static KernelFunction dispatchKernelOnElecInteractionType(const bool elecInteractionTypeIsEwald,
                                                          const bool vdwModifierIsPotSwitch,
                                                          const bool useSimd)
{
    if (elecInteractionTypeIsEwald)
    {
        return (dispatchKernelOnVdwModifier<useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, true>(
                vdwModifierIsPotSwitch, useSimd));
    }
    else
    {
        return (dispatchKernelOnVdwModifier<useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, false>(
                vdwModifierIsPotSwitch, useSimd));
    }
}

template<bool useSoftCore, bool scLambdasOrAlphasDiffer>
static KernelFunction dispatchKernelOnVdwInteractionType(const bool vdwInteractionTypeIsEwald,
                                                         const bool elecInteractionTypeIsEwald,
                                                         const bool vdwModifierIsPotSwitch,
                                                         const bool useSimd)
{
    if (vdwInteractionTypeIsEwald)
    {
        return (dispatchKernelOnElecInteractionType<useSoftCore, scLambdasOrAlphasDiffer, true>(
                elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd));
    }
    else
    {
        return (dispatchKernelOnElecInteractionType<useSoftCore, scLambdasOrAlphasDiffer, false>(
                elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd));
    }
}

template<bool useSoftCore>
static KernelFunction dispatchKernelOnScLambdasOrAlphasDifference(const bool scLambdasOrAlphasDiffer,
                                                                  const bool vdwInteractionTypeIsEwald,
                                                                  const bool elecInteractionTypeIsEwald,
                                                                  const bool vdwModifierIsPotSwitch,
                                                                  const bool useSimd)
{
    if (scLambdasOrAlphasDiffer)
    {
        return (dispatchKernelOnVdwInteractionType<useSoftCore, true>(
                vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd));
    }
    else
    {
        return (dispatchKernelOnVdwInteractionType<useSoftCore, false>(
                vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd));
    }
}

static KernelFunction dispatchKernel(const bool        scLambdasOrAlphasDiffer,
                                     const bool        vdwInteractionTypeIsEwald,
                                     const bool        elecInteractionTypeIsEwald,
                                     const bool        vdwModifierIsPotSwitch,
                                     const bool        useSimd,
                                     const t_forcerec* fr)
{
    if (fr->sc_alphacoul == 0 && fr->sc_alphavdw == 0)
    {
        return (dispatchKernelOnScLambdasOrAlphasDifference<false>(
                scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald,
                vdwModifierIsPotSwitch, useSimd));
    }
    else
    {
        return (dispatchKernelOnScLambdasOrAlphasDifference<true>(
                scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald,
                vdwModifierIsPotSwitch, useSimd));
    }
}


void gmx_nb_free_energy_kernel(const t_nblist*            nlist,
                               rvec*                      xx,
                               gmx::ForceWithShiftForces* ff,
                               const t_forcerec*          fr,
                               const t_mdatoms*           mdatoms,
                               nb_kernel_data_t*          kernel_data,
                               t_nrnb*                    nrnb)
{
    GMX_ASSERT(EEL_PME_EWALD(fr->ic->eeltype) || fr->ic->eeltype == eelCUT || EEL_RF(fr->ic->eeltype),
               "Unsupported eeltype with free energy");

    const bool vdwInteractionTypeIsEwald  = (EVDW_PME(fr->ic->vdwtype));
    const bool elecInteractionTypeIsEwald = (EEL_PME_EWALD(fr->ic->eeltype));
    const bool vdwModifierIsPotSwitch     = (fr->ic->vdw_modifier == eintmodPOTSWITCH);
    bool       scLambdasOrAlphasDiffer    = true;
    const bool useSimd                    = fr->use_simd_kernels;

    if (fr->sc_alphacoul == 0 && fr->sc_alphavdw == 0)
    {
        scLambdasOrAlphasDiffer = false;
    }
    else if (fr->sc_r_power == 6.0_real)
    {
        if (kernel_data->lambda[efptCOUL] == kernel_data->lambda[efptVDW] && fr->sc_alphacoul == fr->sc_alphavdw)
        {
            scLambdasOrAlphasDiffer = false;
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Unsupported soft-core r-power");
    }

    KernelFunction kernelFunc;
    kernelFunc = dispatchKernel(scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald,
                                elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd, fr);
    kernelFunc(nlist, xx, ff, fr, mdatoms, kernel_data, nrnb);
}
