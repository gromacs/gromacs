/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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


{
    const nbnxn_ci_t   *nbln;
    const nbnxn_cj_t   *l_cj;
    const real         *q;
    const real         *shiftvec;
    const real         *x;
    real                facel;
    int                 n, ci, ci_sh;
    int                 ish, ish3;
    gmx_bool            do_LJ, half_LJ, do_coul;
    int                 cjind0, cjind1, cjind;

#ifdef ENERGY_GROUPS
    int         Vstride_i;
    int         egps_ishift, egps_imask;
    int         egps_jshift, egps_jmask, egps_jstride;
    int         egps_i;
    real       *vvdwtp[UNROLLI];
    real       *vctp[UNROLLI];
#endif

    SimdReal  shX_S;
    SimdReal  shY_S;
    SimdReal  shZ_S;
    SimdReal  ix_S0, iy_S0, iz_S0;
    SimdReal  ix_S2, iy_S2, iz_S2;
    SimdReal  fix_S0, fiy_S0, fiz_S0;
    SimdReal  fix_S2, fiy_S2, fiz_S2;

    SimdReal  diagonal_jmi_S;
#if UNROLLI == UNROLLJ
    SimdBool  diagonal_mask_S0, diagonal_mask_S2;
#else
    SimdBool  diagonal_mask0_S0, diagonal_mask0_S2;
    SimdBool  diagonal_mask1_S0, diagonal_mask1_S2;
#endif

    unsigned            *exclusion_filter;
    SimdBitMask          filter_S0, filter_S2;

    SimdReal             zero_S(0.0);

    SimdReal             one_S(1.0);
    SimdReal             iq_S0 = setZero();
    SimdReal             iq_S2 = setZero();

#ifdef CALC_COUL_RF
    SimdReal      mrc_3_S;
#ifdef CALC_ENERGIES
    SimdReal      hrc_3_S, moh_rc_S;
#endif
#endif

#ifdef CALC_COUL_TAB
    /* Coulomb table variables */
    SimdReal          invtsp_S;
    const real       *tab_coul_F;
#if defined CALC_ENERGIES && !defined TAB_FDV0
    const real       *tab_coul_V;
#endif

#ifdef CALC_ENERGIES
    SimdReal   mhalfsp_S;
#endif
#endif

#ifdef CALC_COUL_EWALD
    SimdReal beta2_S, beta_S;
#endif

#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
    SimdReal  sh_ewald_S;
#endif

#if defined LJ_CUT && defined CALC_ENERGIES
    SimdReal   p6_cpot_S, p12_cpot_S;
#endif
#ifdef LJ_POT_SWITCH
    SimdReal   rswitch_S;
    SimdReal   swV3_S, swV4_S, swV5_S;
    SimdReal   swF2_S, swF3_S, swF4_S;
#endif
#ifdef LJ_FORCE_SWITCH
    SimdReal   rswitch_S;
    SimdReal   p6_fc2_S, p6_fc3_S;
    SimdReal   p12_fc2_S, p12_fc3_S;
#ifdef CALC_ENERGIES
    SimdReal   p6_vc3_S, p6_vc4_S;
    SimdReal   p12_vc3_S, p12_vc4_S;
    SimdReal   p6_6cpot_S, p12_12cpot_S;
#endif
#endif
#ifdef LJ_EWALD_GEOM
    real              lj_ewaldcoeff2, lj_ewaldcoeff6_6;
    SimdReal          mone_S, half_S, lje_c2_S, lje_c6_6_S;
#endif

#ifdef LJ_COMB_LB
    const real       *ljc;

    SimdReal          hsig_i_S0, seps_i_S0;
    SimdReal          hsig_i_S2, seps_i_S2;
#else
#ifdef FIX_LJ_C
    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH)  pvdw_c6[2*UNROLLI*UNROLLJ];
    real  *pvdw_c12 = pvdw_c6 + UNROLLI*UNROLLJ;
#endif

#if defined LJ_COMB_GEOM || defined LJ_EWALD_GEOM
    const real       *ljc;
#endif
#endif /* LJ_COMB_LB */

    SimdReal  minRsq_S;
    SimdReal  rc2_S;
#ifdef VDW_CUTOFF_CHECK
    SimdReal  rcvdw2_S;
#endif

    int ninner;

#ifdef COUNT_PAIRS
    int npair = 0;
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB || defined LJ_EWALD_GEOM
    ljc = nbat->lj_comb;
#endif
#if !(defined LJ_COMB_GEOM || defined LJ_COMB_LB || defined FIX_LJ_C)
    /* No combination rule used */
    real      *nbfp_ptr = nbat->nbfp_aligned;
    const int *type     = nbat->type;
#endif

    /* Load j-i for the first i */
    diagonal_jmi_S    = load<SimdReal>(nbat->simd_2xnn_diagonal_j_minus_i);
    /* Generate all the diagonal masks as comparison results */
#if UNROLLI == UNROLLJ
    diagonal_mask_S0  = (zero_S < diagonal_jmi_S);
    diagonal_jmi_S    = diagonal_jmi_S - one_S;
    diagonal_jmi_S    = diagonal_jmi_S - one_S;
    diagonal_mask_S2  = (zero_S < diagonal_jmi_S);
#else
#if 2*UNROLLI == UNROLLJ
    diagonal_mask0_S0 = (zero_S < diagonal_jmi_S);
    diagonal_jmi_S    = diagonal_jmi_S - one_S;
    diagonal_jmi_S    = diagonal_jmi_S - one_S;
    diagonal_mask0_S2 = (zero_S < diagonal_jmi_S);
    diagonal_jmi_S    = diagonal_jmi_S - one_S;
    diagonal_jmi_S    = diagonal_jmi_S - one_S;
    diagonal_mask1_S0 = (zero_S < diagonal_jmi_S);
    diagonal_jmi_S    = diagonal_jmi_S - one_S;
    diagonal_jmi_S    = diagonal_jmi_S - one_S;
    diagonal_mask1_S2 = (zero_S < diagonal_jmi_S);
#endif
#endif

    /* Load masks for topology exclusion masking. filter_stride is
       static const, so the conditional will be optimized away. */
#if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
    exclusion_filter = nbat->simd_exclusion_filter64;
#else
    exclusion_filter = nbat->simd_exclusion_filter;
#endif

    /* Here we cast the exclusion filters from unsigned * to int * or real *.
     * Since we only check bits, the actual value they represent does not
     * matter, as long as both filter and mask data are treated the same way.
     */
#if GMX_SIMD_HAVE_INT32_LOGICAL
    filter_S0 = load<SimdBitMask>(reinterpret_cast<const int *>(exclusion_filter + 0*UNROLLJ));
    filter_S2 = load<SimdBitMask>(reinterpret_cast<const int *>(exclusion_filter + 2*UNROLLJ));
#else
    filter_S0 = load<SimdBitMask>(reinterpret_cast<const real *>(exclusion_filter + 0*UNROLLJ));
    filter_S2 = load<SimdBitMask>(reinterpret_cast<const real *>(exclusion_filter + 2*UNROLLJ));
#endif

#ifdef CALC_COUL_RF
    /* Reaction-field constants */
    mrc_3_S  = SimdReal(-2*ic->k_rf);
#ifdef CALC_ENERGIES
    hrc_3_S  = SimdReal(ic->k_rf);
    moh_rc_S = SimdReal(-ic->c_rf);
#endif
#endif

#ifdef CALC_COUL_TAB

    invtsp_S  = SimdReal(ic->tabq_scale);
#ifdef CALC_ENERGIES
    mhalfsp_S = SimdReal(-0.5/ic->tabq_scale);
#endif

#ifdef TAB_FDV0
    tab_coul_F = ic->tabq_coul_FDV0;
#else
    tab_coul_F = ic->tabq_coul_F;
#ifdef CALC_ENERGIES
    tab_coul_V = ic->tabq_coul_V;
#endif
#endif
#endif /* CALC_COUL_TAB */

#ifdef CALC_COUL_EWALD
    beta2_S = SimdReal(ic->ewaldcoeff_q*ic->ewaldcoeff_q);
    beta_S  = SimdReal(ic->ewaldcoeff_q);
#endif

#if (defined CALC_COUL_TAB || defined CALC_COUL_EWALD) && defined CALC_ENERGIES
    sh_ewald_S = SimdReal(ic->sh_ewald);
#endif

    /* LJ function constants */
#if defined CALC_ENERGIES || defined LJ_POT_SWITCH
    SimdReal sixth_S      = SimdReal(1.0/6.0);
    SimdReal twelveth_S   = SimdReal(1.0/12.0);
#endif

#if defined LJ_CUT && defined CALC_ENERGIES
    /* We shift the potential by cpot, which can be zero */
    p6_cpot_S    = SimdReal(ic->dispersion_shift.cpot);
    p12_cpot_S   = SimdReal(ic->repulsion_shift.cpot);
#endif
#ifdef LJ_POT_SWITCH
    rswitch_S = SimdReal(ic->rvdw_switch);
    swV3_S    = SimdReal(ic->vdw_switch.c3);
    swV4_S    = SimdReal(ic->vdw_switch.c4);
    swV5_S    = SimdReal(ic->vdw_switch.c5);
    swF2_S    = SimdReal(3*ic->vdw_switch.c3);
    swF3_S    = SimdReal(4*ic->vdw_switch.c4);
    swF4_S    = SimdReal(5*ic->vdw_switch.c5);
#endif
#ifdef LJ_FORCE_SWITCH
    rswitch_S = SimdReal(ic->rvdw_switch);
    p6_fc2_S  = SimdReal(ic->dispersion_shift.c2);
    p6_fc3_S  = SimdReal(ic->dispersion_shift.c3);
    p12_fc2_S = SimdReal(ic->repulsion_shift.c2);
    p12_fc3_S = SimdReal(ic->repulsion_shift.c3);
#ifdef CALC_ENERGIES
    {
        SimdReal mthird_S  = SimdReal(-1.0/3.0);
        SimdReal mfourth_S = SimdReal(-1.0/4.0);

        p6_vc3_S     = mthird_S * p6_fc2_S;
        p6_vc4_S     = mfourth_S * p6_fc3_S;
        p6_6cpot_S   = SimdReal(ic->dispersion_shift.cpot/6);
        p12_vc3_S    = mthird_S * p12_fc2_S;
        p12_vc4_S    = mfourth_S * p12_fc3_S;
        p12_12cpot_S = SimdReal(ic->repulsion_shift.cpot/12);
    }
#endif
#endif
#ifdef LJ_EWALD_GEOM
    mone_S           = SimdReal(-1.0);
    half_S           = SimdReal(0.5);
    lj_ewaldcoeff2   = ic->ewaldcoeff_lj*ic->ewaldcoeff_lj;
    lj_ewaldcoeff6_6 = lj_ewaldcoeff2*lj_ewaldcoeff2*lj_ewaldcoeff2/6;
    lje_c2_S         = SimdReal(lj_ewaldcoeff2);
    lje_c6_6_S       = SimdReal(lj_ewaldcoeff6_6);
#ifdef CALC_ENERGIES
    /* Determine the grid potential at the cut-off */
    SimdReal lje_vc_S = SimdReal(ic->sh_lj_ewald);
#endif
#endif

    /* The kernel either supports rcoulomb = rvdw or rcoulomb >= rvdw */
    rc2_S    = SimdReal(ic->rcoulomb*ic->rcoulomb);
#ifdef VDW_CUTOFF_CHECK
    rcvdw2_S = SimdReal(ic->rvdw*ic->rvdw);
#endif

    minRsq_S            = SimdReal(NBNXN_MIN_RSQ);

    q                   = nbat->q;
    facel               = ic->epsfac;
    shiftvec            = shift_vec[0];
    x                   = nbat->x;

#ifdef FIX_LJ_C

    for (jp = 0; jp < UNROLLJ; jp++)
    {
        pvdw_c6 [0*UNROLLJ+jp] = nbat->nbfp[0*2];
        pvdw_c6 [1*UNROLLJ+jp] = nbat->nbfp[0*2];
        pvdw_c6 [2*UNROLLJ+jp] = nbat->nbfp[0*2];
        pvdw_c6 [3*UNROLLJ+jp] = nbat->nbfp[0*2];

        pvdw_c12[0*UNROLLJ+jp] = nbat->nbfp[0*2+1];
        pvdw_c12[1*UNROLLJ+jp] = nbat->nbfp[0*2+1];
        pvdw_c12[2*UNROLLJ+jp] = nbat->nbfp[0*2+1];
        pvdw_c12[3*UNROLLJ+jp] = nbat->nbfp[0*2+1];
    }
    SimdReal c6_S0  = load<SimdReal>(pvdw_c6 +0*UNROLLJ);
    SimdReal c6_S1  = load<SimdReal>(pvdw_c6 +1*UNROLLJ);
    SimdReal c6_S2  = load<SimdReal>(pvdw_c6 +2*UNROLLJ);
    SimdReal c6_S3  = load<SimdReal>(pvdw_c6 +3*UNROLLJ);

    SimdReal c12_S0 = load<SimdReal>(pvdw_c12+0*UNROLLJ);
    SimdReal c12_S1 = load<SimdReal>(pvdw_c12+1*UNROLLJ);
    SimdReal c12_S2 = load<SimdReal>(pvdw_c12+2*UNROLLJ);
    SimdReal c12_S3 = load<SimdReal>(pvdw_c12+3*UNROLLJ);
#endif /* FIX_LJ_C */

#ifdef ENERGY_GROUPS
    egps_ishift  = nbat->neg_2log;
    egps_imask   = (1<<egps_ishift) - 1;
    egps_jshift  = 2*nbat->neg_2log;
    egps_jmask   = (1<<egps_jshift) - 1;
    egps_jstride = (UNROLLJ>>1)*UNROLLJ;
    /* Major division is over i-particle energy groups, determine the stride */
    Vstride_i    = nbat->nenergrp*(1<<nbat->neg_2log)*egps_jstride;
#endif

    l_cj = nbl->cj;

    ninner = 0;
    for (n = 0; n < nbl->nci; n++)
    {
        nbln = &nbl->ci[n];

        ish              = (nbln->shift & NBNXN_CI_SHIFT);
        ish3             = ish*3;
        cjind0           = nbln->cj_ind_start;
        cjind1           = nbln->cj_ind_end;
        ci               = nbln->ci;
        ci_sh            = (ish == CENTRAL ? ci : -1);

        shX_S = SimdReal(shiftvec[ish3]);
        shY_S = SimdReal(shiftvec[ish3+1]);
        shZ_S = SimdReal(shiftvec[ish3+2]);

#if UNROLLJ <= 4
        int sci              = ci*STRIDE;
        int scix             = sci*DIM;
#if defined LJ_COMB_LB || defined LJ_COMB_GEOM || defined LJ_EWALD_GEOM
        int sci2             = sci*2;
#endif
#else
        int sci              = (ci>>1)*STRIDE;
        int scix             = sci*DIM + (ci & 1)*(STRIDE>>1);
#if defined LJ_COMB_LB || defined LJ_COMB_GEOM || defined LJ_EWALD_GEOM
        int sci2             = sci*2 + (ci & 1)*(STRIDE>>1);
#endif
        sci             += (ci & 1)*(STRIDE>>1);
#endif

        /* We have 5 LJ/C combinations, but use only three inner loops,
         * as the other combinations are unlikely and/or not much faster:
         * inner half-LJ + C for half-LJ + C / no-LJ + C
         * inner LJ + C      for full-LJ + C
         * inner LJ          for full-LJ + no-C / half-LJ + no-C
         */
        do_LJ   = (nbln->shift & NBNXN_CI_DO_LJ(0));
        do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));
        half_LJ = ((nbln->shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;

#ifdef ENERGY_GROUPS
        egps_i = nbat->energrp[ci];
        {
            int ia, egp_ia;

            for (ia = 0; ia < UNROLLI; ia++)
            {
                egp_ia     = (egps_i >> (ia*egps_ishift)) & egps_imask;
                vvdwtp[ia] = Vvdw + egp_ia*Vstride_i;
                vctp[ia]   = Vc   + egp_ia*Vstride_i;
            }
        }
#endif

#ifdef CALC_ENERGIES
#ifdef LJ_EWALD_GEOM
        gmx_bool do_self = TRUE;
#else
        gmx_bool do_self = do_coul;
#endif
#if UNROLLJ == 4
        if (do_self && l_cj[nbln->cj_ind_start].cj == ci_sh)
#endif
#if UNROLLJ == 8
        if (do_self && l_cj[nbln->cj_ind_start].cj == (ci_sh>>1))
#endif
        {
            if (do_coul)
            {
                real Vc_sub_self;
                int  ia;

#ifdef CALC_COUL_RF
                Vc_sub_self = 0.5*ic->c_rf;
#endif
#ifdef CALC_COUL_TAB
#ifdef TAB_FDV0
                Vc_sub_self = 0.5*tab_coul_F[2];
#else
                Vc_sub_self = 0.5*tab_coul_V[0];
#endif
#endif
#ifdef CALC_COUL_EWALD
                /* beta/sqrt(pi) */
                Vc_sub_self = 0.5*ic->ewaldcoeff_q*M_2_SQRTPI;
#endif

                for (ia = 0; ia < UNROLLI; ia++)
                {
                    real qi;

                    qi = q[sci+ia];
#ifdef ENERGY_GROUPS
                    vctp[ia][((egps_i>>(ia*egps_ishift)) & egps_imask)*egps_jstride]
#else
                    Vc[0]
#endif
                        -= facel*qi*qi*Vc_sub_self;
                }
            }

#ifdef LJ_EWALD_GEOM
            {
                int  ia;

                for (ia = 0; ia < UNROLLI; ia++)
                {
                    real c6_i;

                    c6_i = nbat->nbfp[nbat->type[sci+ia]*(nbat->ntype + 1)*2]/6;
#ifdef ENERGY_GROUPS
                    vvdwtp[ia][((egps_i>>(ia*egps_ishift)) & egps_imask)*egps_jstride]
#else
                    Vvdw[0]
#endif
                        += 0.5*c6_i*lj_ewaldcoeff6_6;
                }
            }
#endif      /* LJ_EWALD */
        }
#endif

        /* Load i atom data */
        int sciy             = scix + STRIDE;
        int sciz             = sciy + STRIDE;
        ix_S0          = loadU1DualHsimd(x+scix);
        ix_S2          = loadU1DualHsimd(x+scix+2);
        iy_S0          = loadU1DualHsimd(x+sciy);
        iy_S2          = loadU1DualHsimd(x+sciy+2);
        iz_S0          = loadU1DualHsimd(x+sciz);
        iz_S2          = loadU1DualHsimd(x+sciz+2);
        ix_S0          = ix_S0 + shX_S;
        ix_S2          = ix_S2 + shX_S;
        iy_S0          = iy_S0 + shY_S;
        iy_S2          = iy_S2 + shY_S;
        iz_S0          = iz_S0 + shZ_S;
        iz_S2          = iz_S2 + shZ_S;

        if (do_coul)
        {
            SimdReal facel_S;

            facel_S    = SimdReal(facel);

            iq_S0      = loadU1DualHsimd(q+sci);
            iq_S2      = loadU1DualHsimd(q+sci+2);
            iq_S0      = facel_S * iq_S0;
            iq_S2      = facel_S * iq_S2;
        }

#ifdef LJ_COMB_LB
        hsig_i_S0 = loadU1DualHsimd(ljc+sci2);
        hsig_i_S2 = loadU1DualHsimd(ljc+sci2+2);
        seps_i_S0 = loadU1DualHsimd(ljc+sci2+STRIDE);
        seps_i_S2 = loadU1DualHsimd(ljc+sci2+STRIDE+2);
#else
#ifdef LJ_COMB_GEOM
        SimdReal   c6s_S0, c12s_S0;
        SimdReal   c6s_S2, c12s_S2;

        c6s_S0 = loadU1DualHsimd(ljc+sci2);

        if (!half_LJ)
        {
            c6s_S2 = loadU1DualHsimd(ljc+sci2+2);
        }
        c12s_S0 = loadU1DualHsimd(ljc+sci2+STRIDE);
        if (!half_LJ)
        {
            c12s_S2 = loadU1DualHsimd(ljc+sci2+STRIDE+2);
        }
#elif !defined LJ_COMB_LB && !defined FIX_LJ_C
        const real *nbfp0     = nbfp_ptr + type[sci  ]*nbat->ntype*c_simdBestPairAlignment;
        const real *nbfp1     = nbfp_ptr + type[sci+1]*nbat->ntype*c_simdBestPairAlignment;
        const real *nbfp2     = NULL, *nbfp3 = NULL;
        if (!half_LJ)
        {
            nbfp2 = nbfp_ptr + type[sci+2]*nbat->ntype*c_simdBestPairAlignment;
            nbfp3 = nbfp_ptr + type[sci+3]*nbat->ntype*c_simdBestPairAlignment;
        }
#endif
#endif
#ifdef LJ_EWALD_GEOM
        /* We need the geometrically combined C6 for the PME grid correction */
        SimdReal c6s_S0, c6s_S2;
        c6s_S0 = loadU1DualHsimd(ljc+sci2);
        if (!half_LJ)
        {
            c6s_S2 = loadU1DualHsimd(ljc+sci2+2);
        }
#endif

        /* Zero the potential energy for this list */
#ifdef CALC_ENERGIES
        SimdReal Vvdwtot_S = setZero();
        SimdReal vctot_S   = setZero();
#endif

        /* Clear i atom forces */
        fix_S0           = setZero();
        fix_S2           = setZero();
        fiy_S0           = setZero();
        fiy_S2           = setZero();
        fiz_S0           = setZero();
        fiz_S2           = setZero();

        cjind = cjind0;

        /* Currently all kernels use (at least half) LJ */
#define CALC_LJ
        if (half_LJ)
        {
            /* Coulomb: all i-atoms, LJ: first half i-atoms */
#define CALC_COULOMB
#define HALF_LJ
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != NBNXN_INTERACTION_MASK_ALL)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn_inner.h"
            }
#undef HALF_LJ
#undef CALC_COULOMB
        }
        else if (do_coul)
        {
            /* Coulomb: all i-atoms, LJ: all i-atoms */
#define CALC_COULOMB
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != NBNXN_INTERACTION_MASK_ALL)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn_inner.h"
            }
#undef CALC_COULOMB
        }
        else
        {
            /* Coulomb: none, LJ: all i-atoms */
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != NBNXN_INTERACTION_MASK_ALL)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn_inner.h"
            }
        }
#undef CALC_LJ
        ninner += cjind1 - cjind0;

        /* Add accumulated i-forces to the force array */
        real fShiftX = reduceIncr4ReturnSumHsimd(f+scix, fix_S0, fix_S2);
        real fShiftY = reduceIncr4ReturnSumHsimd(f+sciy, fiy_S0, fiy_S2);
        real fShiftZ = reduceIncr4ReturnSumHsimd(f+sciz, fiz_S0, fiz_S2);

#ifdef CALC_SHIFTFORCES
        fshift[ish3+0] += fShiftX;
        fshift[ish3+1] += fShiftY;
        fshift[ish3+2] += fShiftZ;
#endif

#ifdef CALC_ENERGIES
        if (do_coul)
        {
            *Vc += reduce(vctot_S);
        }
        *Vvdw += reduce(Vvdwtot_S);
#endif

        /* Outer loop uses 6 flops/iteration */
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}
