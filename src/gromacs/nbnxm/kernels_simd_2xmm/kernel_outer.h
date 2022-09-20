/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

{
    using namespace gmx;

#if GMX_SIMD_J_UNROLL_SIZE == 1
    constexpr KernelLayout kernelLayout = KernelLayout::r4xM;
#else
    constexpr KernelLayout            kernelLayout         = KernelLayout::r2xMM;
#endif

    /* The number of 'i' SIMD registers */
    static_assert(UNROLLI % GMX_SIMD_J_UNROLL_SIZE == 0);
    constexpr int nR = UNROLLI / GMX_SIMD_J_UNROLL_SIZE;

    /* Interaction type and output choice specific constexpr variables */
#ifdef CALC_COUL_RF
    constexpr KernelCoulombType coulombType = KernelCoulombType::RF;
#endif
#ifdef CALC_COUL_EWALD
    constexpr KernelCoulombType coulombType = KernelCoulombType::EwaldAnalytical;
#endif
#ifdef CALC_COUL_TAB
    constexpr KernelCoulombType coulombType = KernelCoulombType::EwaldTabulated;
#endif
#ifdef LJ_COMB_GEOM
    constexpr LJCombinationRule ljCombinationRule = LJCombinationRule::Geometric;
#else
#    ifdef LJ_COMB_LB
    constexpr LJCombinationRule       ljCombinationRule    = LJCombinationRule::LorentzBerthelot;
#    else
    constexpr LJCombinationRule ljCombinationRule = LJCombinationRule::None;
#    endif
#endif
#ifdef LJ_POT_SWITCH
    constexpr InteractionModifiers vdwModifier = InteractionModifiers::PotSwitch;
#else
#    ifdef LJ_FORCE_SWITCH
    constexpr InteractionModifiers    vdwModifier          = InteractionModifiers::ForceSwitch;
#    else
    /* Note the we also use the potential-shift kernel for LJ without shift */
    constexpr InteractionModifiers vdwModifier = InteractionModifiers::PotShift;
#    endif
#endif
#ifdef LJ_EWALD_GEOM
    constexpr bool haveLJEwaldGeometric = true;
#else
    constexpr bool                    haveLJEwaldGeometric = false;
#endif
#ifdef CALC_ENERGIES
    constexpr bool calculateEnergies = true;
#else
    constexpr bool                    calculateEnergies    = false;
#endif
#if defined CALC_ENERGIES && defined ENERGY_GROUPS
    constexpr bool useEnergyGroups = true;
#else
    constexpr bool                    useEnergyGroups      = false;
#endif
#ifdef VDW_CUTOFF_CHECK
    constexpr bool haveVdwCutoffCheck = true;
#else
    constexpr bool                    haveVdwCutoffCheck   = false;
#endif

    /* Unpack pointers for output */
    real* f      = out->f.data();
    real* fshift = out->fshift.data();
    real* Vvdw;
    real* Vc;
    if constexpr (calculateEnergies)
    {
        if constexpr (useEnergyGroups)
        {
            Vvdw = out->VSvdw.data();
            Vc   = out->VSc.data();
        }
        else
        {
            Vvdw = out->Vvdw.data();
            Vc   = out->Vc.data();
        }
    }

    SimdBitMask filter_S0, filter_S1, filter_S2, filter_S3;

    SimdReal zero_S(0.0);

#ifdef COUNT_PAIRS
    int npair = 0;
#endif

    const nbnxn_atomdata_t::Params& nbatParams = nbat->params();

    const real* gmx_restrict ljc;
    if constexpr (ljCombinationRule != LJCombinationRule::None || haveLJEwaldGeometric)
    {
        ljc = nbatParams.lj_comb.data();
    }
    const real* gmx_restrict nbfp_ptr;
    const int* gmx_restrict  type;
    if constexpr (ljCombinationRule == LJCombinationRule::None)
    {
        /* No combination rule used */
        nbfp_ptr = nbatParams.nbfp_aligned.data();
        type     = nbatParams.type.data();
    }

    /* Set up the diagonal exclusion masks */
    const DiagonalMasker<nR, kernelLayout, getDiagonalMaskType<UNROLLI, UNROLLJ>()> diagonalMasker(
            nbat->simdMasks);

#if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
    const std::uint64_t* gmx_restrict exclusion_filter = nbat->simdMasks.exclusion_filter64.data();
#else
    const std::uint32_t* gmx_restrict exclusion_filter = nbat->simdMasks.exclusion_filter.data();
#endif

    /* Here we cast the exclusion filters from unsigned * to int * or real *.
     * Since we only check bits, the actual value they represent does not
     * matter, as long as both filter and mask data are treated the same way.
     */
    SimdBitMask exclusionFilterV[nR];
    for (int i = 0; i < nR; i++)
    {
#if GMX_SIMD_HAVE_INT32_LOGICAL
        exclusionFilterV[i] = load<SimdBitMask>(
                reinterpret_cast<const int*>(exclusion_filter + i * GMX_SIMD_J_UNROLL_SIZE * UNROLLJ));
#else
        exclusionFilterV[i] = load<SimdBitMask>(
                reinterpret_cast<const real*>(exclusion_filter + i * GMX_SIMD_J_UNROLL_SIZE * UNROLLJ));
#endif
    }

    CoulombCalculator<coulombType> coulombCalculator(*ic);

    gmx_unused SimdReal ewaldShift;
    if constexpr (coulombType != KernelCoulombType::RF && calculateEnergies)
    {
        ewaldShift = SimdReal(ic->sh_ewald);
    }

    /* LJ function constants, only actually needed with energies or potential switching */
    SimdReal sixth_S(1.0_real / 6.0_real);
    SimdReal twelveth_S(1.0_real / 12.0_real);

    LennardJonesCalculator<calculateEnergies, vdwModifier> ljCalculator(*ic);

    std::array<SimdReal, haveLJEwaldGeometric ? 5 : 0> ljEwaldParams;
    real                                               lj_ewaldcoeff6_6;
    if constexpr (haveLJEwaldGeometric)
    {
        ljEwaldParams[0]          = SimdReal(1.0_real);
        ljEwaldParams[1]          = SimdReal(0.5_real);
        const real lj_ewaldcoeff2 = ic->ewaldcoeff_lj * ic->ewaldcoeff_lj;
        lj_ewaldcoeff6_6          = lj_ewaldcoeff2 * lj_ewaldcoeff2 * lj_ewaldcoeff2 / 6;
        ljEwaldParams[2]          = SimdReal(lj_ewaldcoeff2);
        ljEwaldParams[3]          = SimdReal(lj_ewaldcoeff6_6);
        /* Determine the grid potential at the cut-off */
        ljEwaldParams[4] = ic->sh_lj_ewald;
    }

    /* The kernel either supports rcoulomb = rvdw or rcoulomb >= rvdw */
    const SimdReal cutoffSquared(ic->rcoulomb * ic->rcoulomb);
    SimdReal       vdwCutoffSquared;
    if constexpr (haveVdwCutoffCheck)
    {
        vdwCutoffSquared = SimdReal(ic->rvdw * ic->rvdw);
    }

    const SimdReal minDistanceSquared(c_nbnxnMinDistanceSquared);

    const real* gmx_restrict q        = nbatParams.q.data();
    const real               facel    = ic->epsfac;
    const real* gmx_restrict shiftvec = shift_vec[0];
    const real* gmx_restrict x        = nbat->x().data();


    // These constants are only used when useEnergyGroups==true
    const int egps_ishift  = nbatParams.neg_2log;
    const int egps_imask   = (1 << egps_ishift) - 1;
    const int egps_jshift  = 2 * nbatParams.neg_2log;
    const int egps_jmask   = (1 << egps_jshift) - 1;
    const int egps_jstride = (UNROLLJ >> 1) * UNROLLJ;
    /* Major division is over i-particle energy groups, determine the stride */
    const int Vstride_i = nbatParams.nenergrp * (1 << nbatParams.neg_2log) * egps_jstride;

    const nbnxn_cj_t* l_cj = nbl->cj.list_.data();

    for (const nbnxn_ci_t& ciEntry : nbl->ci)
    {
        const int ish    = (ciEntry.shift & NBNXN_CI_SHIFT);
        const int ish3   = ish * 3;
        const int cjind0 = ciEntry.cj_ind_start;
        const int cjind1 = ciEntry.cj_ind_end;
        const int ci     = ciEntry.ci;
        const int ci_sh  = (ish == gmx::c_centralShiftIndex ? ci : -1);

        // Load the periodic shift vector for the i-atoms
        const SimdReal iShiftX(shiftvec[ish3]);
        const SimdReal iShiftY(shiftvec[ish3 + 1]);
        const SimdReal iShiftZ(shiftvec[ish3 + 2]);

#if UNROLLJ <= 4
        int sci  = ci * STRIDE;
        int scix = sci * DIM;
        //#    if defined LJ_COMB_LB || defined LJ_COMB_GEOM || defined LJ_EWALD_GEOM
        int sci2 = sci * 2;
        //#    endif
#else
        int sci  = (ci >> 1) * STRIDE;
        int scix = sci * DIM + (ci & 1) * (STRIDE >> 1);
        //#    if defined LJ_COMB_LB || defined LJ_COMB_GEOM || defined LJ_EWALD_GEOM
        int sci2 = sci * 2 + (ci & 1) * (STRIDE >> 1);
        //#    endif
        sci += (ci & 1) * (STRIDE >> 1);
#endif

        /* We have 5 LJ/C combinations, but use only three inner loops,
         * as the other combinations are unlikely and/or not much faster:
         * inner half-LJ + C for half-LJ + C / no-LJ + C
         * inner LJ + C      for full-LJ + C
         * inner LJ          for full-LJ + no-C / half-LJ + no-C
         */
        const bool do_LJ   = ((ciEntry.shift & NBNXN_CI_DO_LJ(0)) != 0);
        const bool do_coul = ((ciEntry.shift & NBNXN_CI_DO_COUL(0)) != 0);
        const bool half_LJ = (((ciEntry.shift & NBNXN_CI_HALF_LJ(0)) != 0) || !do_LJ) && do_coul;

        std::array<real*, useEnergyGroups ? UNROLLI : 0> vvdwtp;
        std::array<real*, useEnergyGroups ? UNROLLI : 0> vctp;
        int                                              egps_i;
        if constexpr (useEnergyGroups)
        {
            egps_i = nbatParams.energrp[ci];
            for (int ia = 0; ia < UNROLLI; ia++)
            {
                int egp_ia = (egps_i >> (ia * egps_ishift)) & egps_imask;
                vvdwtp[ia] = Vvdw + egp_ia * Vstride_i;
                vctp[ia]   = Vc + egp_ia * Vstride_i;
            }
        }

        if constexpr (calculateEnergies)
        {
            // Compute self interaction energies, when present
            const bool do_self = haveLJEwaldGeometric || do_coul;

#if UNROLLJ == 4
            if (do_self && l_cj[ciEntry.cj_ind_start].cj == ci_sh)
#endif
#if UNROLLJ == 2
                if (do_self && l_cj[ciEntry.cj_ind_start].cj == (ci_sh << 1))
#endif
#if UNROLLJ == 8
                    if (do_self && l_cj[ciEntry.cj_ind_start].cj == (ci_sh >> 1))
#endif
                    {
                        if (do_coul)
                        {
                            const real Vc_sub_self = coulombCalculator.selfEnergy();

                            for (int ia = 0; ia < UNROLLI; ia++)
                            {
                                const real qi = q[sci + ia];
                                if constexpr (useEnergyGroups)
                                {
                                    vctp[ia][((egps_i >> (ia * egps_ishift)) & egps_imask) * egps_jstride] -=
                                            facel * qi * qi * Vc_sub_self;
                                }
                                else
                                {
                                    Vc[0] -= facel * qi * qi * Vc_sub_self;
                                }
                            }
                        }

                        if constexpr (haveLJEwaldGeometric)
                        {
                            for (int ia = 0; ia < UNROLLI; ia++)
                            {
                                real c6_i =
                                        nbatParams.nbfp[nbatParams.type[sci + ia] * (nbatParams.numTypes + 1) * 2]
                                        / 6;
                                if constexpr (useEnergyGroups)
                                {
                                    vvdwtp[ia][((egps_i >> (ia * egps_ishift)) & egps_imask) * egps_jstride] +=
                                            0.5 * c6_i * lj_ewaldcoeff6_6;
                                }
                                else
                                {
                                    Vvdw[0] += 0.5 * c6_i * lj_ewaldcoeff6_6;
                                }
                            }
                        }
                    }

        } // calulateEnergies

        /* Load i atom data */
        const int sciy = scix + STRIDE;
        const int sciz = sciy + STRIDE;

        const auto ixV =
                genArr<nR>([&](int i) { return loadIAtomData<kernelLayout>(x, scix, i) + iShiftX; });
        const auto iyV =
                genArr<nR>([&](int i) { return loadIAtomData<kernelLayout>(x, sciy, i) + iShiftY; });
        const auto izV =
                genArr<nR>([&](int i) { return loadIAtomData<kernelLayout>(x, sciz, i) + iShiftZ; });

        std::array<SimdReal, nR> chargeIV;
        if (do_coul)
        {
            chargeIV =
                    genArr<nR>([&](int i) { return facel * loadIAtomData<kernelLayout>(q, sci, i); });
        }

        constexpr bool c_ljCombLB = (ljCombinationRule == LJCombinationRule::LorentzBerthelot);
        // Note that when half_lj==true we actually only need nR/2 LJ parameters
        std::array<SimdReal, c_ljCombLB ? nR : 0> halfSigmaIV;
        std::array<SimdReal, c_ljCombLB ? nR : 0> sqrtEpsilonIV;
        std::array<SimdReal, (ljCombinationRule == LJCombinationRule::Geometric || haveLJEwaldGeometric) ? nR : 0> c6GeomV;
        std::array<SimdReal, (ljCombinationRule == LJCombinationRule::Geometric) ? nR : 0> c12GeomV;
        std::array<const real*, (ljCombinationRule == LJCombinationRule::None) ? UNROLLI : 0> nbfpI;
        if constexpr (c_ljCombLB)
        {
            for (int i = 0; i < nR; i++)
            {
                halfSigmaIV[i]   = loadIAtomData<kernelLayout>(ljc, sci2, i);
                sqrtEpsilonIV[i] = loadIAtomData<kernelLayout>(ljc, sci2 + STRIDE, i);
            }
        }
        else if constexpr (ljCombinationRule == LJCombinationRule::Geometric)
        {
            for (int i = 0; i < nR; i++)
            {
                c6GeomV[i]  = loadIAtomData<kernelLayout>(ljc, sci2, i);
                c12GeomV[i] = loadIAtomData<kernelLayout>(ljc, sci2 + STRIDE, i);
            }
        }
        else
        {
            const int numTypes = nbatParams.numTypes;
            for (int i = 0; i < UNROLLI; i++)
            {
                nbfpI[i] = nbfp_ptr + type[sci + i] * numTypes * c_simdBestPairAlignment;
            }
        }
        if constexpr (haveLJEwaldGeometric && ljCombinationRule != LJCombinationRule::Geometric)
        {
            /* We need the geometrically combined C6 for the PME grid correction */
            for (int i = 0; i < nR; i++)
            {
                c6GeomV[i] = loadIAtomData<kernelLayout>(ljc, sci2, i);
            }
        }

        SimdReal Vvdwtot_S;
        SimdReal vctot_S;
        if constexpr (calculateEnergies)
        {
            /* Zero the potential energy for this list */
            Vvdwtot_S = setZero();
            vctot_S   = setZero();
        }

        /* Declare and clear i atom forces */
        auto forceIXV = genArr<nR>([&](int gmx_unused i) {
            SimdReal tmp = setZero();
            return tmp;
        });
        auto forceIYV = genArr<nR>([&](int gmx_unused i) {
            SimdReal tmp = setZero();
            return tmp;
        });
        auto forceIZV = genArr<nR>([&](int gmx_unused i) {
            SimdReal tmp = setZero();
            return tmp;
        });


        int cjind = cjind0;

        /* Currently all kernels use (at least half) LJ */
        if (half_LJ)
        {
            /* Coulomb: all i-atoms, LJ: first half i-atoms */
            constexpr bool            c_calculateCoulombInteractions = true;
            constexpr ILJInteractions c_iLJInteractions              = ILJInteractions::Half;
            {
                constexpr bool c_needToCheckExclusions = true;
                while (cjind < cjind1 && nbl->cj.excl(cjind) != NBNXN_INTERACTION_MASK_ALL)
                {
#include "kernel_inner.h"
                    cjind++;
                }
            }
            {
                constexpr bool c_needToCheckExclusions = false;
                for (; (cjind < cjind1); cjind++)
                {
#include "kernel_inner.h"
                }
            }
        }
        else if (do_coul)
        {
            /* Coulomb: all i-atoms, LJ: all i-atoms */
            constexpr bool            c_calculateCoulombInteractions = true;
            constexpr ILJInteractions c_iLJInteractions              = ILJInteractions::All;
            {
                constexpr bool c_needToCheckExclusions = true;
                while (cjind < cjind1 && nbl->cj.excl(cjind) != NBNXN_INTERACTION_MASK_ALL)
                {
#include "kernel_inner.h"
                    cjind++;
                }
            }
            {
                constexpr bool c_needToCheckExclusions = false;
                for (; (cjind < cjind1); cjind++)
                {
#include "kernel_inner.h"
                }
            }
        }
        else
        {
            /* Coulomb: none, LJ: all i-atoms */
            constexpr bool            c_calculateCoulombInteractions = false;
            constexpr ILJInteractions c_iLJInteractions              = ILJInteractions::All;
            {
                constexpr bool c_needToCheckExclusions = true;
                while (cjind < cjind1 && nbl->cj.excl(cjind) != NBNXN_INTERACTION_MASK_ALL)
                {
#include "kernel_inner.h"
                    cjind++;
                }
            }
            {
                constexpr bool c_needToCheckExclusions = false;
                for (; (cjind < cjind1); cjind++)
                {
#include "kernel_inner.h"
                }
            }
        }
        /* Add accumulated i-forces to the force array */
#if GMX_SIMD_J_UNROLL_SIZE == 1
        real fShiftX =
                reduceIncr4ReturnSum(f + scix, forceIXV[0], forceIXV[1], forceIXV[2], forceIXV[3]);
        real fShiftY =
                reduceIncr4ReturnSum(f + sciy, forceIYV[0], forceIYV[1], forceIYV[2], forceIYV[3]);
        real fShiftZ =
                reduceIncr4ReturnSum(f + sciz, forceIZV[0], forceIZV[1], forceIZV[2], forceIZV[3]);
#else
        real fShiftX = reduceIncr4ReturnSumHsimd(f + scix, forceIXV[0], forceIXV[1]);
        real fShiftY = reduceIncr4ReturnSumHsimd(f + sciy, forceIYV[0], forceIYV[1]);
        real fShiftZ = reduceIncr4ReturnSumHsimd(f + sciz, forceIZV[0], forceIZV[1]);
#endif

#ifdef CALC_SHIFTFORCES
        fshift[ish3 + 0] += fShiftX;
        fshift[ish3 + 1] += fShiftY;
        fshift[ish3 + 2] += fShiftZ;
#endif

        if constexpr (calculateEnergies)
        {
            if (do_coul)
            {
                *Vc += reduce(vctot_S);
            }

            *Vvdw += reduce(Vvdwtot_S);
        }

        /* Outer loop uses 6 flops/iteration */
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}
