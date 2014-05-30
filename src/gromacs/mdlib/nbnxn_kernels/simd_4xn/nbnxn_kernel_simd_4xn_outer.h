/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
    const int *         type;
    const real *        q;
    const real         *shiftvec;
    const real         *x;
    const real         *nbfp0, *nbfp1, *nbfp2 = NULL, *nbfp3 = NULL;
    real                facel;
    real               *nbfp_ptr;
    int                 n, ci, ci_sh;
    int                 ish, ish3;
    gmx_bool            do_LJ, half_LJ, do_coul, do_self;
    int                 sci, scix, sciy, sciz, sci2;
    int                 cjind0, cjind1, cjind;
    int                 ip, jp;

#ifdef ENERGY_GROUPS
    int         Vstride_i;
    int         egps_ishift, egps_imask;
    int         egps_jshift, egps_jmask, egps_jstride;
    int         egps_i;
    real       *vvdwtp[UNROLLI];
    real       *vctp[UNROLLI];
#endif

    gmx_simd_real_t  shX_S;
    gmx_simd_real_t  shY_S;
    gmx_simd_real_t  shZ_S;
    gmx_simd_real_t  ix_S0, iy_S0, iz_S0;
    gmx_simd_real_t  ix_S1, iy_S1, iz_S1;
    gmx_simd_real_t  ix_S2, iy_S2, iz_S2;
    gmx_simd_real_t  ix_S3, iy_S3, iz_S3;
    gmx_simd_real_t  fix_S0, fiy_S0, fiz_S0;
    gmx_simd_real_t  fix_S1, fiy_S1, fiz_S1;
    gmx_simd_real_t  fix_S2, fiy_S2, fiz_S2;
    gmx_simd_real_t  fix_S3, fiy_S3, fiz_S3;
#if UNROLLJ >= 4
    /* We use an i-force SIMD register width of 4 */
    gmx_simd4_real_t fix_S, fiy_S, fiz_S;
#else
    /* We use an i-force SIMD register width of 2 */
    gmx_simd_real_t  fix0_S, fiy0_S, fiz0_S;
    gmx_simd_real_t  fix2_S, fiy2_S, fiz2_S;
#endif

    gmx_simd_real_t  diagonal_jmi_S;
#if UNROLLI == UNROLLJ
    gmx_simd_bool_t  diagonal_mask_S0, diagonal_mask_S1, diagonal_mask_S2, diagonal_mask_S3;
#else
    gmx_simd_bool_t  diagonal_mask0_S0, diagonal_mask0_S1, diagonal_mask0_S2, diagonal_mask0_S3;
    gmx_simd_bool_t  diagonal_mask1_S0, diagonal_mask1_S1, diagonal_mask1_S2, diagonal_mask1_S3;
#endif

    unsigned            *exclusion_filter;
    gmx_exclfilter       filter_S0, filter_S1, filter_S2, filter_S3;

    gmx_simd_real_t      zero_S = gmx_simd_set1_r(0.0);

    gmx_simd_real_t      one_S  = gmx_simd_set1_r(1.0);
    gmx_simd_real_t      iq_S0  = gmx_simd_setzero_r();
    gmx_simd_real_t      iq_S1  = gmx_simd_setzero_r();
    gmx_simd_real_t      iq_S2  = gmx_simd_setzero_r();
    gmx_simd_real_t      iq_S3  = gmx_simd_setzero_r();

#ifdef CALC_COUL_RF
    gmx_simd_real_t      mrc_3_S;
#ifdef CALC_ENERGIES
    gmx_simd_real_t      hrc_3_S, moh_rc_S;
#endif
#endif

#ifdef CALC_COUL_TAB
    /* Coulomb table variables */
    gmx_simd_real_t   invtsp_S;
    const real *      tab_coul_F;
#ifndef TAB_FDV0
    const real *      tab_coul_V;
#endif
    /* Thread-local working buffers for force and potential lookups */
    int               ti0_array[2*GMX_SIMD_REAL_WIDTH], *ti0 = NULL;
    int               ti1_array[2*GMX_SIMD_REAL_WIDTH], *ti1 = NULL;
    int               ti2_array[2*GMX_SIMD_REAL_WIDTH], *ti2 = NULL;
    int               ti3_array[2*GMX_SIMD_REAL_WIDTH], *ti3 = NULL;
#ifdef CALC_ENERGIES
    gmx_simd_real_t   mhalfsp_S;
#endif
#endif

#ifdef CALC_COUL_EWALD
    gmx_simd_real_t beta2_S, beta_S;
#endif

#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
    gmx_simd_real_t  sh_ewald_S;
#endif

#if defined LJ_CUT && defined CALC_ENERGIES
    gmx_simd_real_t   p6_cpot_S, p12_cpot_S;
#endif
#ifdef LJ_POT_SWITCH
    gmx_simd_real_t   rswitch_S;
    gmx_simd_real_t   swV3_S, swV4_S, swV5_S;
    gmx_simd_real_t   swF2_S, swF3_S, swF4_S;
#endif
#ifdef LJ_FORCE_SWITCH
    gmx_simd_real_t   rswitch_S;
    gmx_simd_real_t   p6_fc2_S, p6_fc3_S;
    gmx_simd_real_t   p12_fc2_S, p12_fc3_S;
#ifdef CALC_ENERGIES
    gmx_simd_real_t   p6_vc3_S, p6_vc4_S;
    gmx_simd_real_t   p12_vc3_S, p12_vc4_S;
    gmx_simd_real_t   p6_6cpot_S, p12_12cpot_S;
#endif
#endif
#ifdef LJ_EWALD_GEOM
    real              lj_ewaldcoeff2, lj_ewaldcoeff6_6;
    gmx_simd_real_t   mone_S, half_S, lje_c2_S, lje_c6_6_S, lje_vc_S;
#endif

#ifdef LJ_COMB_LB
    const real       *ljc;

    gmx_simd_real_t   hsig_i_S0, seps_i_S0;
    gmx_simd_real_t   hsig_i_S1, seps_i_S1;
    gmx_simd_real_t   hsig_i_S2, seps_i_S2;
    gmx_simd_real_t   hsig_i_S3, seps_i_S3;
#else
#ifdef FIX_LJ_C
    real              pvdw_array[2*UNROLLI*UNROLLJ+3];
    real             *pvdw_c6, *pvdw_c12;
    gmx_simd_real_t   c6_S0, c12_S0;
    gmx_simd_real_t   c6_S1, c12_S1;
    gmx_simd_real_t   c6_S2, c12_S2;
    gmx_simd_real_t   c6_S3, c12_S3;
#endif

#if defined LJ_COMB_GEOM || defined LJ_EWALD_GEOM
    const real       *ljc;

    gmx_simd_real_t   c6s_S0, c12s_S0;
    gmx_simd_real_t   c6s_S1, c12s_S1;
    gmx_simd_real_t   c6s_S2  = gmx_simd_setzero_r();
    gmx_simd_real_t   c12s_S2 = gmx_simd_setzero_r();
    gmx_simd_real_t   c6s_S3  = gmx_simd_setzero_r();
    gmx_simd_real_t   c12s_S3 = gmx_simd_setzero_r();
#endif
#endif /* LJ_COMB_LB */

    gmx_simd_real_t  vctot_S, Vvdwtot_S;
    gmx_simd_real_t  sixth_S, twelveth_S;

    gmx_simd_real_t  avoid_sing_S;
    gmx_simd_real_t  rc2_S;
#ifdef VDW_CUTOFF_CHECK
    gmx_simd_real_t  rcvdw2_S;
#endif

    int ninner;

#ifdef COUNT_PAIRS
    int npair = 0;
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB || defined LJ_EWALD_GEOM
    ljc = nbat->lj_comb;
#endif
#if !(defined LJ_COMB_GEOM || defined LJ_COMB_LB)
    /* No combination rule used */
    nbfp_ptr    = (4 == nbfp_stride) ? nbat->nbfp_s4 : nbat->nbfp;
#endif

    /* Load j-i for the first i */
    diagonal_jmi_S    = gmx_simd_load_r(nbat->simd_4xn_diagonal_j_minus_i);
    /* Generate all the diagonal masks as comparison results */
#if UNROLLI == UNROLLJ
    diagonal_mask_S0  = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask_S1  = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask_S2  = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask_S3  = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
#else
#if UNROLLI == 2*UNROLLJ || 2*UNROLLI == UNROLLJ
    diagonal_mask0_S0 = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask0_S1 = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask0_S2 = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask0_S3 = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);

#if UNROLLI == 2*UNROLLJ
    /* Load j-i for the second half of the j-cluster */
    diagonal_jmi_S    = gmx_simd_load_r(nbat->simd_4xn_diagonal_j_minus_i + UNROLLJ);
#endif

    diagonal_mask1_S0 = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask1_S1 = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask1_S2 = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
    diagonal_jmi_S    = gmx_simd_sub_r(diagonal_jmi_S, one_S);
    diagonal_mask1_S3 = gmx_simd_cmplt_r(zero_S, diagonal_jmi_S);
#endif
#endif

    /* Load masks for topology exclusion masking. filter_stride is
       static const, so the conditional will be optimized away. */
    if (1 == filter_stride)
    {
        exclusion_filter = nbat->simd_exclusion_filter1;
    }
    else /* (2 == filter_stride) */
    {
        exclusion_filter = nbat->simd_exclusion_filter2;
    }

    /* Here we cast the exclusion filters from unsigned * to int * or real *.
     * Since we only check bits, the actual value they represent does not
     * matter, as long as both filter and mask data are treated the same way.
     */
    filter_S0 = gmx_load_exclusion_filter(exclusion_filter + 0*UNROLLJ*filter_stride);
    filter_S1 = gmx_load_exclusion_filter(exclusion_filter + 1*UNROLLJ*filter_stride);
    filter_S2 = gmx_load_exclusion_filter(exclusion_filter + 2*UNROLLJ*filter_stride);
    filter_S3 = gmx_load_exclusion_filter(exclusion_filter + 3*UNROLLJ*filter_stride);

#ifdef CALC_COUL_RF
    /* Reaction-field constants */
    mrc_3_S  = gmx_simd_set1_r(-2*ic->k_rf);
#ifdef CALC_ENERGIES
    hrc_3_S  = gmx_simd_set1_r(ic->k_rf);
    moh_rc_S = gmx_simd_set1_r(-ic->c_rf);
#endif
#endif

#ifdef CALC_COUL_TAB
    /* Generate aligned table index pointers */
    ti0 = prepare_table_load_buffer(ti0_array);
    ti1 = prepare_table_load_buffer(ti1_array);
    ti2 = prepare_table_load_buffer(ti2_array);
    ti3 = prepare_table_load_buffer(ti3_array);

    invtsp_S  = gmx_simd_set1_r(ic->tabq_scale);
#ifdef CALC_ENERGIES
    mhalfsp_S = gmx_simd_set1_r(-0.5/ic->tabq_scale);
#endif

#ifdef TAB_FDV0
    tab_coul_F = ic->tabq_coul_FDV0;
#else
    tab_coul_F = ic->tabq_coul_F;
    tab_coul_V = ic->tabq_coul_V;
#endif
#endif /* CALC_COUL_TAB */

#ifdef CALC_COUL_EWALD
    beta2_S = gmx_simd_set1_r(ic->ewaldcoeff_q*ic->ewaldcoeff_q);
    beta_S  = gmx_simd_set1_r(ic->ewaldcoeff_q);
#endif

#if (defined CALC_COUL_TAB || defined CALC_COUL_EWALD) && defined CALC_ENERGIES
    sh_ewald_S = gmx_simd_set1_r(ic->sh_ewald);
#endif

    /* LJ function constants */
#if defined CALC_ENERGIES || defined LJ_POT_SWITCH
    sixth_S      = gmx_simd_set1_r(1.0/6.0);
    twelveth_S   = gmx_simd_set1_r(1.0/12.0);
#endif

#if defined LJ_CUT && defined CALC_ENERGIES
    /* We shift the potential by cpot, which can be zero */
    p6_cpot_S    = gmx_simd_set1_r(ic->dispersion_shift.cpot);
    p12_cpot_S   = gmx_simd_set1_r(ic->repulsion_shift.cpot);
#endif
#ifdef LJ_POT_SWITCH
    rswitch_S = gmx_simd_set1_r(ic->rvdw_switch);
    swV3_S    = gmx_simd_set1_r(ic->vdw_switch.c3);
    swV4_S    = gmx_simd_set1_r(ic->vdw_switch.c4);
    swV5_S    = gmx_simd_set1_r(ic->vdw_switch.c5);
    swF2_S    = gmx_simd_set1_r(3*ic->vdw_switch.c3);
    swF3_S    = gmx_simd_set1_r(4*ic->vdw_switch.c4);
    swF4_S    = gmx_simd_set1_r(5*ic->vdw_switch.c5);
#endif
#ifdef LJ_FORCE_SWITCH
    rswitch_S = gmx_simd_set1_r(ic->rvdw_switch);
    p6_fc2_S  = gmx_simd_set1_r(ic->dispersion_shift.c2);
    p6_fc3_S  = gmx_simd_set1_r(ic->dispersion_shift.c3);
    p12_fc2_S = gmx_simd_set1_r(ic->repulsion_shift.c2);
    p12_fc3_S = gmx_simd_set1_r(ic->repulsion_shift.c3);
#ifdef CALC_ENERGIES
    {
        gmx_simd_real_t mthird_S  = gmx_simd_set1_r(-1.0/3.0);
        gmx_simd_real_t mfourth_S = gmx_simd_set1_r(-1.0/4.0);

        p6_vc3_S     = gmx_simd_mul_r(mthird_S,  p6_fc2_S);
        p6_vc4_S     = gmx_simd_mul_r(mfourth_S, p6_fc3_S);
        p6_6cpot_S   = gmx_simd_set1_r(ic->dispersion_shift.cpot/6);
        p12_vc3_S    = gmx_simd_mul_r(mthird_S,  p12_fc2_S);
        p12_vc4_S    = gmx_simd_mul_r(mfourth_S, p12_fc3_S);
        p12_12cpot_S = gmx_simd_set1_r(ic->repulsion_shift.cpot/12);
    }
#endif
#endif
#ifdef LJ_EWALD_GEOM
    mone_S           = gmx_simd_set1_r(-1.0);
    half_S           = gmx_simd_set1_r(0.5);
    lj_ewaldcoeff2   = ic->ewaldcoeff_lj*ic->ewaldcoeff_lj;
    lj_ewaldcoeff6_6 = lj_ewaldcoeff2*lj_ewaldcoeff2*lj_ewaldcoeff2/6;
    lje_c2_S         = gmx_simd_set1_r(lj_ewaldcoeff2);
    lje_c6_6_S       = gmx_simd_set1_r(lj_ewaldcoeff6_6);
    /* Determine the grid potential at the cut-off */
    lje_vc_S         = gmx_simd_set1_r(ic->sh_lj_ewald);
#endif

    /* The kernel either supports rcoulomb = rvdw or rcoulomb >= rvdw */
    rc2_S    = gmx_simd_set1_r(ic->rcoulomb*ic->rcoulomb);
#ifdef VDW_CUTOFF_CHECK
    rcvdw2_S = gmx_simd_set1_r(ic->rvdw*ic->rvdw);
#endif

    avoid_sing_S = gmx_simd_set1_r(NBNXN_AVOID_SING_R2_INC);

    q                   = nbat->q;
    type                = nbat->type;
    facel               = ic->epsfac;
    shiftvec            = shift_vec[0];
    x                   = nbat->x;

#ifdef FIX_LJ_C
    pvdw_c6  = gmx_simd_align_real(pvdw_array);
    pvdw_c12 = pvdw_c6 + UNROLLI*UNROLLJ;

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
    c6_S0            = gmx_simd_load_r(pvdw_c6 +0*UNROLLJ);
    c6_S1            = gmx_simd_load_r(pvdw_c6 +1*UNROLLJ);
    c6_S2            = gmx_simd_load_r(pvdw_c6 +2*UNROLLJ);
    c6_S3            = gmx_simd_load_r(pvdw_c6 +3*UNROLLJ);

    c12_S0           = gmx_simd_load_r(pvdw_c12+0*UNROLLJ);
    c12_S1           = gmx_simd_load_r(pvdw_c12+1*UNROLLJ);
    c12_S2           = gmx_simd_load_r(pvdw_c12+2*UNROLLJ);
    c12_S3           = gmx_simd_load_r(pvdw_c12+3*UNROLLJ);
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

        shX_S = gmx_simd_load1_r(shiftvec+ish3);
        shY_S = gmx_simd_load1_r(shiftvec+ish3+1);
        shZ_S = gmx_simd_load1_r(shiftvec+ish3+2);

#if UNROLLJ <= 4
        sci              = ci*STRIDE;
        scix             = sci*DIM;
        sci2             = sci*2;
#else
        sci              = (ci>>1)*STRIDE;
        scix             = sci*DIM + (ci & 1)*(STRIDE>>1);
        sci2             = sci*2 + (ci & 1)*(STRIDE>>1);
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
#ifdef LJ_EWALD_GEOM
        do_self = TRUE;
#else
        do_self = do_coul;
#endif

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
#if UNROLLJ == 4
        if (do_self && l_cj[nbln->cj_ind_start].cj == ci_sh)
#endif
#if UNROLLJ == 2
        if (do_self && l_cj[nbln->cj_ind_start].cj == (ci_sh<<1))
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
#endif      /* LJ_EWALD_GEOM */
        }
#endif

        /* Load i atom data */
        sciy             = scix + STRIDE;
        sciz             = sciy + STRIDE;
        ix_S0            = gmx_simd_add_r(gmx_simd_load1_r(x+scix), shX_S);
        ix_S1            = gmx_simd_add_r(gmx_simd_load1_r(x+scix+1), shX_S);
        ix_S2            = gmx_simd_add_r(gmx_simd_load1_r(x+scix+2), shX_S);
        ix_S3            = gmx_simd_add_r(gmx_simd_load1_r(x+scix+3), shX_S);
        iy_S0            = gmx_simd_add_r(gmx_simd_load1_r(x+sciy), shY_S);
        iy_S1            = gmx_simd_add_r(gmx_simd_load1_r(x+sciy+1), shY_S);
        iy_S2            = gmx_simd_add_r(gmx_simd_load1_r(x+sciy+2), shY_S);
        iy_S3            = gmx_simd_add_r(gmx_simd_load1_r(x+sciy+3), shY_S);
        iz_S0            = gmx_simd_add_r(gmx_simd_load1_r(x+sciz), shZ_S);
        iz_S1            = gmx_simd_add_r(gmx_simd_load1_r(x+sciz+1), shZ_S);
        iz_S2            = gmx_simd_add_r(gmx_simd_load1_r(x+sciz+2), shZ_S);
        iz_S3            = gmx_simd_add_r(gmx_simd_load1_r(x+sciz+3), shZ_S);

        if (do_coul)
        {
            iq_S0      = gmx_simd_set1_r(facel*q[sci]);
            iq_S1      = gmx_simd_set1_r(facel*q[sci+1]);
            iq_S2      = gmx_simd_set1_r(facel*q[sci+2]);
            iq_S3      = gmx_simd_set1_r(facel*q[sci+3]);
        }

#ifdef LJ_COMB_LB
        hsig_i_S0      = gmx_simd_load1_r(ljc+sci2+0);
        hsig_i_S1      = gmx_simd_load1_r(ljc+sci2+1);
        hsig_i_S2      = gmx_simd_load1_r(ljc+sci2+2);
        hsig_i_S3      = gmx_simd_load1_r(ljc+sci2+3);
        seps_i_S0      = gmx_simd_load1_r(ljc+sci2+STRIDE+0);
        seps_i_S1      = gmx_simd_load1_r(ljc+sci2+STRIDE+1);
        seps_i_S2      = gmx_simd_load1_r(ljc+sci2+STRIDE+2);
        seps_i_S3      = gmx_simd_load1_r(ljc+sci2+STRIDE+3);
#else
#ifdef LJ_COMB_GEOM
        c6s_S0         = gmx_simd_load1_r(ljc+sci2+0);
        c6s_S1         = gmx_simd_load1_r(ljc+sci2+1);
        if (!half_LJ)
        {
            c6s_S2     = gmx_simd_load1_r(ljc+sci2+2);
            c6s_S3     = gmx_simd_load1_r(ljc+sci2+3);
        }
        c12s_S0        = gmx_simd_load1_r(ljc+sci2+STRIDE+0);
        c12s_S1        = gmx_simd_load1_r(ljc+sci2+STRIDE+1);
        if (!half_LJ)
        {
            c12s_S2    = gmx_simd_load1_r(ljc+sci2+STRIDE+2);
            c12s_S3    = gmx_simd_load1_r(ljc+sci2+STRIDE+3);
        }
#else
        nbfp0     = nbfp_ptr + type[sci  ]*nbat->ntype*nbfp_stride;
        nbfp1     = nbfp_ptr + type[sci+1]*nbat->ntype*nbfp_stride;
        if (!half_LJ)
        {
            nbfp2 = nbfp_ptr + type[sci+2]*nbat->ntype*nbfp_stride;
            nbfp3 = nbfp_ptr + type[sci+3]*nbat->ntype*nbfp_stride;
        }
#endif
#endif
#ifdef LJ_EWALD_GEOM
        /* We need the geometrically combined C6 for the PME grid correction */
        c6s_S0 = gmx_simd_load1_r(ljc+sci2+0);
        c6s_S1 = gmx_simd_load1_r(ljc+sci2+1);
        if (!half_LJ)
        {
            c6s_S2 = gmx_simd_load1_r(ljc+sci2+2);
            c6s_S3 = gmx_simd_load1_r(ljc+sci2+3);
        }
#endif

        /* Zero the potential energy for this list */
        Vvdwtot_S        = gmx_simd_setzero_r();
        vctot_S          = gmx_simd_setzero_r();

        /* Clear i atom forces */
        fix_S0           = gmx_simd_setzero_r();
        fix_S1           = gmx_simd_setzero_r();
        fix_S2           = gmx_simd_setzero_r();
        fix_S3           = gmx_simd_setzero_r();
        fiy_S0           = gmx_simd_setzero_r();
        fiy_S1           = gmx_simd_setzero_r();
        fiy_S2           = gmx_simd_setzero_r();
        fiy_S3           = gmx_simd_setzero_r();
        fiz_S0           = gmx_simd_setzero_r();
        fiz_S1           = gmx_simd_setzero_r();
        fiz_S2           = gmx_simd_setzero_r();
        fiz_S3           = gmx_simd_setzero_r();

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
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn_inner.h"
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
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn_inner.h"
            }
#undef CALC_COULOMB
        }
        else
        {
            /* Coulomb: none, LJ: all i-atoms */
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != NBNXN_INTERACTION_MASK_ALL)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn_inner.h"
            }
        }
#undef CALC_LJ
        ninner += cjind1 - cjind0;

        /* Add accumulated i-forces to the force array */
#if UNROLLJ >= 4
        fix_S = gmx_mm_transpose_sum4_pr(fix_S0, fix_S1, fix_S2, fix_S3);
        gmx_simd4_store_r(f+scix, gmx_simd4_add_r(fix_S, gmx_simd4_load_r(f+scix)));

        fiy_S = gmx_mm_transpose_sum4_pr(fiy_S0, fiy_S1, fiy_S2, fiy_S3);
        gmx_simd4_store_r(f+sciy, gmx_simd4_add_r(fiy_S, gmx_simd4_load_r(f+sciy)));

        fiz_S = gmx_mm_transpose_sum4_pr(fiz_S0, fiz_S1, fiz_S2, fiz_S3);
        gmx_simd4_store_r(f+sciz, gmx_simd4_add_r(fiz_S, gmx_simd4_load_r(f+sciz)));

#ifdef CALC_SHIFTFORCES
        fshift[ish3+0] += gmx_simd4_reduce_r(fix_S);
        fshift[ish3+1] += gmx_simd4_reduce_r(fiy_S);
        fshift[ish3+2] += gmx_simd4_reduce_r(fiz_S);
#endif
#else
        fix0_S = gmx_mm_transpose_sum2_pr(fix_S0, fix_S1);
        gmx_simd_store_r(f+scix, gmx_simd_add_r(fix0_S, gmx_simd_load_r(f+scix)));
        fix2_S = gmx_mm_transpose_sum2_pr(fix_S2, fix_S3);
        gmx_simd_store_r(f+scix+2, gmx_simd_add_r(fix2_S, gmx_simd_load_r(f+scix+2)));

        fiy0_S = gmx_mm_transpose_sum2_pr(fiy_S0, fiy_S1);
        gmx_simd_store_r(f+sciy, gmx_simd_add_r(fiy0_S, gmx_simd_load_r(f+sciy)));
        fiy2_S = gmx_mm_transpose_sum2_pr(fiy_S2, fiy_S3);
        gmx_simd_store_r(f+sciy+2, gmx_simd_add_r(fiy2_S, gmx_simd_load_r(f+sciy+2)));

        fiz0_S = gmx_mm_transpose_sum2_pr(fiz_S0, fiz_S1);
        gmx_simd_store_r(f+sciz, gmx_simd_add_r(fiz0_S, gmx_simd_load_r(f+sciz)));
        fiz2_S = gmx_mm_transpose_sum2_pr(fiz_S2, fiz_S3);
        gmx_simd_store_r(f+sciz+2, gmx_simd_add_r(fiz2_S, gmx_simd_load_r(f+sciz+2)));

#ifdef CALC_SHIFTFORCES
        fshift[ish3+0] += gmx_simd_reduce_r(gmx_simd_add_r(fix0_S, fix2_S));
        fshift[ish3+1] += gmx_simd_reduce_r(gmx_simd_add_r(fiy0_S, fiy2_S));
        fshift[ish3+2] += gmx_simd_reduce_r(gmx_simd_add_r(fiz0_S, fiz2_S));
#endif
#endif

#ifdef CALC_ENERGIES
        if (do_coul)
        {
            *Vc += gmx_simd_reduce_r(vctot_S);
        }

        *Vvdw += gmx_simd_reduce_r(Vvdwtot_S);
#endif

        /* Outer loop uses 6 flops/iteration */
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}
