/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

/* This is the innermost loop contents for the 4 x N atom SIMD kernel.
 * This flavor of the kernel calculates interactions of 4 i-atoms
 * with N j-atoms stored in N wide SIMD registers.
 */


/* When calculating RF or Ewald interactions we calculate the electrostatic
 * forces on excluded atom pairs here in the non-bonded loops.
 * But when energies and/or virial is required we calculate them
 * separately to as then it is easier to separate the energy and virial
 * contributions.
 */
#if defined CHECK_EXCLS && defined CALC_COULOMB
#define EXCL_FORCES
#endif

/* Without exclusions and energies we only need to mask the cut-off,
 * this can be faster when we have defined gmx_blendv_pr, i.e. an instruction
 * that selects from two SIMD registers based on the contents of a third.
 */
#if !(defined CHECK_EXCLS || defined CALC_ENERGIES) && defined GMX_SIMD_HAVE_BLENDV
/* With RF and tabulated Coulomb we replace cmp+and with sub+blendv.
 * With gcc this is slower, except for RF on Sandy Bridge.
 * Tested with gcc 4.6.2, 4.6.3 and 4.7.1.
 */
#if (defined CALC_COUL_RF || defined CALC_COUL_TAB) && (!defined __GNUC__ || (defined CALC_COUL_RF && defined GMX_X86_AVX_256))
#define NBNXN_CUTOFF_USE_BLENDV
#endif
/* With analytical Ewald we replace cmp+and+and with sub+blendv+blendv.
 * This is only faster with icc on Sandy Bridge (PS kernel slower than gcc 4.7).
 * Tested with icc 13.
 */
#if defined CALC_COUL_EWALD && defined __INTEL_COMPILER && defined GMX_X86_AVX_256
#define NBNXN_CUTOFF_USE_BLENDV
#endif
#endif

{
    int        cj, aj, ajx, ajy, ajz;

#ifdef ENERGY_GROUPS
    /* Energy group indices for two atoms packed into one int */
    int        egp_jj[UNROLLJ/2];
#endif

#ifdef CHECK_EXCLS
    /* Interaction (non-exclusion) mask of all 1's or 0's */
    gmx_mm_pb  interact_S0;
    gmx_mm_pb  interact_S1;
    gmx_mm_pb  interact_S2;
    gmx_mm_pb  interact_S3;
#endif

    gmx_mm_pr  jx_S, jy_S, jz_S;
    gmx_mm_pr  dx_S0, dy_S0, dz_S0;
    gmx_mm_pr  dx_S1, dy_S1, dz_S1;
    gmx_mm_pr  dx_S2, dy_S2, dz_S2;
    gmx_mm_pr  dx_S3, dy_S3, dz_S3;
    gmx_mm_pr  tx_S0, ty_S0, tz_S0;
    gmx_mm_pr  tx_S1, ty_S1, tz_S1;
    gmx_mm_pr  tx_S2, ty_S2, tz_S2;
    gmx_mm_pr  tx_S3, ty_S3, tz_S3;
    gmx_mm_pr  rsq_S0, rinv_S0, rinvsq_S0;
    gmx_mm_pr  rsq_S1, rinv_S1, rinvsq_S1;
    gmx_mm_pr  rsq_S2, rinv_S2, rinvsq_S2;
    gmx_mm_pr  rsq_S3, rinv_S3, rinvsq_S3;
#ifndef NBNXN_CUTOFF_USE_BLENDV
    /* wco: within cut-off, mask of all 1's or 0's */
    gmx_mm_pb  wco_S0;
    gmx_mm_pb  wco_S1;
    gmx_mm_pb  wco_S2;
    gmx_mm_pb  wco_S3;
#endif
#ifdef VDW_CUTOFF_CHECK
    gmx_mm_pb  wco_vdw_S0;
    gmx_mm_pb  wco_vdw_S1;
#ifndef HALF_LJ
    gmx_mm_pb  wco_vdw_S2;
    gmx_mm_pb  wco_vdw_S3;
#endif
#endif
#ifdef CALC_COULOMB
#ifdef CHECK_EXCLS
    /* 1/r masked with the interaction mask */
    gmx_mm_pr  rinv_ex_S0;
    gmx_mm_pr  rinv_ex_S1;
    gmx_mm_pr  rinv_ex_S2;
    gmx_mm_pr  rinv_ex_S3;
#endif
    gmx_mm_pr  jq_S;
    gmx_mm_pr  qq_S0;
    gmx_mm_pr  qq_S1;
    gmx_mm_pr  qq_S2;
    gmx_mm_pr  qq_S3;
#ifdef CALC_COUL_TAB
    /* The force (PME mesh force) we need to subtract from 1/r^2 */
    gmx_mm_pr  fsub_S0;
    gmx_mm_pr  fsub_S1;
    gmx_mm_pr  fsub_S2;
    gmx_mm_pr  fsub_S3;
#endif
#ifdef CALC_COUL_EWALD
    gmx_mm_pr  brsq_S0, brsq_S1, brsq_S2, brsq_S3;
    gmx_mm_pr  ewcorr_S0, ewcorr_S1, ewcorr_S2, ewcorr_S3;
#endif

    /* frcoul = (1/r - fsub)*r */
    gmx_mm_pr  frcoul_S0;
    gmx_mm_pr  frcoul_S1;
    gmx_mm_pr  frcoul_S2;
    gmx_mm_pr  frcoul_S3;
#ifdef CALC_COUL_TAB
    /* For tables: r, rs=r/sp, rf=floor(rs), frac=rs-rf */
    gmx_mm_pr  r_S0, rs_S0, rf_S0, frac_S0;
    gmx_mm_pr  r_S1, rs_S1, rf_S1, frac_S1;
    gmx_mm_pr  r_S2, rs_S2, rf_S2, frac_S2;
    gmx_mm_pr  r_S3, rs_S3, rf_S3, frac_S3;
    /* Table index: rs truncated to an int */
    gmx_epi32  ti_S0, ti_S1, ti_S2, ti_S3;
    /* Linear force table values */
    gmx_mm_pr  ctab0_S0, ctab1_S0;
    gmx_mm_pr  ctab0_S1, ctab1_S1;
    gmx_mm_pr  ctab0_S2, ctab1_S2;
    gmx_mm_pr  ctab0_S3, ctab1_S3;
#ifdef CALC_ENERGIES
    /* Quadratic energy table value */
    gmx_mm_pr  ctabv_S0;
    gmx_mm_pr  ctabv_S1;
    gmx_mm_pr  ctabv_S2;
    gmx_mm_pr  ctabv_S3;
#endif
#endif
#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
    /* The potential (PME mesh) we need to subtract from 1/r */
    gmx_mm_pr  vc_sub_S0;
    gmx_mm_pr  vc_sub_S1;
    gmx_mm_pr  vc_sub_S2;
    gmx_mm_pr  vc_sub_S3;
#endif
#ifdef CALC_ENERGIES
    /* Electrostatic potential */
    gmx_mm_pr  vcoul_S0;
    gmx_mm_pr  vcoul_S1;
    gmx_mm_pr  vcoul_S2;
    gmx_mm_pr  vcoul_S3;
#endif
#endif
    /* The force times 1/r */
    gmx_mm_pr  fscal_S0;
    gmx_mm_pr  fscal_S1;
    gmx_mm_pr  fscal_S2;
    gmx_mm_pr  fscal_S3;

#ifdef CALC_LJ
#ifdef LJ_COMB_LB
    /* LJ sigma_j/2 and sqrt(epsilon_j) */
    gmx_mm_pr  hsig_j_S, seps_j_S;
    /* LJ sigma_ij and epsilon_ij */
    gmx_mm_pr  sig_S0, eps_S0;
    gmx_mm_pr  sig_S1, eps_S1;
#ifndef HALF_LJ
    gmx_mm_pr  sig_S2, eps_S2;
    gmx_mm_pr  sig_S3, eps_S3;
#endif
#ifdef CALC_ENERGIES
    gmx_mm_pr  sig2_S0, sig6_S0;
    gmx_mm_pr  sig2_S1, sig6_S1;
#ifndef HALF_LJ
    gmx_mm_pr  sig2_S2, sig6_S2;
    gmx_mm_pr  sig2_S3, sig6_S3;
#endif
#endif /* LJ_COMB_LB */
#endif /* CALC_LJ */

#ifdef LJ_COMB_GEOM
    gmx_mm_pr  c6s_j_S, c12s_j_S;
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB
    /* Index for loading LJ parameters, complicated when interleaving */
    int         aj2;
#endif

#ifndef FIX_LJ_C
    /* LJ C6 and C12 parameters, used with geometric comb. rule */
    gmx_mm_pr  c6_S0, c12_S0;
    gmx_mm_pr  c6_S1, c12_S1;
#ifndef HALF_LJ
    gmx_mm_pr  c6_S2, c12_S2;
    gmx_mm_pr  c6_S3, c12_S3;
#endif
#endif

    /* Intermediate variables for LJ calculation */
#ifndef LJ_COMB_LB
    gmx_mm_pr  rinvsix_S0;
    gmx_mm_pr  rinvsix_S1;
#ifndef HALF_LJ
    gmx_mm_pr  rinvsix_S2;
    gmx_mm_pr  rinvsix_S3;
#endif
#endif
#ifdef LJ_COMB_LB
    gmx_mm_pr  sir_S0, sir2_S0, sir6_S0;
    gmx_mm_pr  sir_S1, sir2_S1, sir6_S1;
#ifndef HALF_LJ
    gmx_mm_pr  sir_S2, sir2_S2, sir6_S2;
    gmx_mm_pr  sir_S3, sir2_S3, sir6_S3;
#endif
#endif

    gmx_mm_pr  FrLJ6_S0, FrLJ12_S0;
    gmx_mm_pr  FrLJ6_S1, FrLJ12_S1;
#ifndef HALF_LJ
    gmx_mm_pr  FrLJ6_S2, FrLJ12_S2;
    gmx_mm_pr  FrLJ6_S3, FrLJ12_S3;
#endif
#ifdef CALC_ENERGIES
    gmx_mm_pr  VLJ6_S0, VLJ12_S0, VLJ_S0;
    gmx_mm_pr  VLJ6_S1, VLJ12_S1, VLJ_S1;
#ifndef HALF_LJ
    gmx_mm_pr  VLJ6_S2, VLJ12_S2, VLJ_S2;
    gmx_mm_pr  VLJ6_S3, VLJ12_S3, VLJ_S3;
#endif
#endif
#endif /* CALC_LJ */

    /* j-cluster index */
    cj            = l_cj[cjind].cj;

    /* Atom indices (of the first atom in the cluster) */
    aj            = cj*UNROLLJ;
#if defined CALC_LJ && (defined LJ_COMB_GEOM || defined LJ_COMB_LB)
#if UNROLLJ == STRIDE
    aj2           = aj*2;
#else
    aj2           = (cj>>1)*2*STRIDE + (cj & 1)*UNROLLJ;
#endif
#endif
#if UNROLLJ == STRIDE
    ajx           = aj*DIM;
#else
    ajx           = (cj>>1)*DIM*STRIDE + (cj & 1)*UNROLLJ;
#endif
    ajy           = ajx + STRIDE;
    ajz           = ajy + STRIDE;

#ifdef CHECK_EXCLS
    gmx_load_simd_4xn_interactions(l_cj[cjind].excl, filter_S0, filter_S1, filter_S2, filter_S3, &interact_S0, &interact_S1, &interact_S2, &interact_S3);
#endif /* CHECK_EXCLS */

    /* load j atom coordinates */
    jx_S        = gmx_load_pr(x+ajx);
    jy_S        = gmx_load_pr(x+ajy);
    jz_S        = gmx_load_pr(x+ajz);

    /* Calculate distance */
    dx_S0       = gmx_sub_pr(ix_S0, jx_S);
    dy_S0       = gmx_sub_pr(iy_S0, jy_S);
    dz_S0       = gmx_sub_pr(iz_S0, jz_S);
    dx_S1       = gmx_sub_pr(ix_S1, jx_S);
    dy_S1       = gmx_sub_pr(iy_S1, jy_S);
    dz_S1       = gmx_sub_pr(iz_S1, jz_S);
    dx_S2       = gmx_sub_pr(ix_S2, jx_S);
    dy_S2       = gmx_sub_pr(iy_S2, jy_S);
    dz_S2       = gmx_sub_pr(iz_S2, jz_S);
    dx_S3       = gmx_sub_pr(ix_S3, jx_S);
    dy_S3       = gmx_sub_pr(iy_S3, jy_S);
    dz_S3       = gmx_sub_pr(iz_S3, jz_S);

    /* rsq = dx*dx+dy*dy+dz*dz */
    rsq_S0      = gmx_calc_rsq_pr(dx_S0, dy_S0, dz_S0);
    rsq_S1      = gmx_calc_rsq_pr(dx_S1, dy_S1, dz_S1);
    rsq_S2      = gmx_calc_rsq_pr(dx_S2, dy_S2, dz_S2);
    rsq_S3      = gmx_calc_rsq_pr(dx_S3, dy_S3, dz_S3);

#ifndef NBNXN_CUTOFF_USE_BLENDV
    wco_S0      = gmx_cmplt_pr(rsq_S0, rc2_S);
    wco_S1      = gmx_cmplt_pr(rsq_S1, rc2_S);
    wco_S2      = gmx_cmplt_pr(rsq_S2, rc2_S);
    wco_S3      = gmx_cmplt_pr(rsq_S3, rc2_S);
#endif

#ifdef CHECK_EXCLS
#ifdef EXCL_FORCES
    /* Only remove the (sub-)diagonal to avoid double counting */
#if UNROLLJ == UNROLLI
    if (cj == ci_sh)
    {
        wco_S0  = gmx_and_pb(wco_S0, diagonal_mask_S0);
        wco_S1  = gmx_and_pb(wco_S1, diagonal_mask_S1);
        wco_S2  = gmx_and_pb(wco_S2, diagonal_mask_S2);
        wco_S3  = gmx_and_pb(wco_S3, diagonal_mask_S3);
    }
#else
#if UNROLLJ < UNROLLI
    if (cj == ci_sh*2)
    {
        wco_S0  = gmx_and_pb(wco_S0, diagonal_mask0_S0);
        wco_S1  = gmx_and_pb(wco_S1, diagonal_mask0_S1);
        wco_S2  = gmx_and_pb(wco_S2, diagonal_mask0_S2);
        wco_S3  = gmx_and_pb(wco_S3, diagonal_mask0_S3);
    }
    if (cj == ci_sh*2 + 1)
    {
        wco_S0  = gmx_and_pb(wco_S0, diagonal_mask1_S0);
        wco_S1  = gmx_and_pb(wco_S1, diagonal_mask1_S1);
        wco_S2  = gmx_and_pb(wco_S2, diagonal_mask1_S2);
        wco_S3  = gmx_and_pb(wco_S3, diagonal_mask1_S3);
    }
#else
    if (cj*2 == ci_sh)
    {
        wco_S0  = gmx_and_pb(wco_S0, diagonal_mask0_S0);
        wco_S1  = gmx_and_pb(wco_S1, diagonal_mask0_S1);
        wco_S2  = gmx_and_pb(wco_S2, diagonal_mask0_S2);
        wco_S3  = gmx_and_pb(wco_S3, diagonal_mask0_S3);
    }
    else if (cj*2 + 1 == ci_sh)
    {
        wco_S0  = gmx_and_pb(wco_S0, diagonal_mask1_S0);
        wco_S1  = gmx_and_pb(wco_S1, diagonal_mask1_S1);
        wco_S2  = gmx_and_pb(wco_S2, diagonal_mask1_S2);
        wco_S3  = gmx_and_pb(wco_S3, diagonal_mask1_S3);
    }
#endif
#endif
#else /* EXCL_FORCES */
    /* No exclusion forces: remove all excluded atom pairs from the list */
    wco_S0      = gmx_and_pb(wco_S0, interact_S0);
    wco_S1      = gmx_and_pb(wco_S1, interact_S1);
    wco_S2      = gmx_and_pb(wco_S2, interact_S2);
    wco_S3      = gmx_and_pb(wco_S3, interact_S3);
#endif
#endif

#ifdef COUNT_PAIRS
    {
        int  i, j;
        real tmpa[2*GMX_SIMD_WIDTH_HERE], *tmp;
        tmp = gmx_simd_align_real(tmpa);
        for (i = 0; i < UNROLLI; i++)
        {
            gmx_store_pr(tmp, gmx_sub_pr(rc2_S, i == 0 ? rsq_S0 : (i == 1 ? rsq_S1 : (i == 2 ? rsq_S2 : rsq_S3))));
            for (j = 0; j < UNROLLJ; j++)
            {
                if (tmp[j] >= 0)
                {
                    npair++;
                }
            }
        }
    }
#endif

#ifdef CHECK_EXCLS
    /* For excluded pairs add a small number to avoid r^-6 = NaN */
    rsq_S0      = gmx_masknot_add_pr(interact_S0, rsq_S0, avoid_sing_S);
    rsq_S1      = gmx_masknot_add_pr(interact_S1, rsq_S1, avoid_sing_S);
    rsq_S2      = gmx_masknot_add_pr(interact_S2, rsq_S2, avoid_sing_S);
    rsq_S3      = gmx_masknot_add_pr(interact_S3, rsq_S3, avoid_sing_S);
#endif

    /* Calculate 1/r */
#ifndef GMX_DOUBLE
    rinv_S0     = gmx_invsqrt_pr(rsq_S0);
    rinv_S1     = gmx_invsqrt_pr(rsq_S1);
    rinv_S2     = gmx_invsqrt_pr(rsq_S2);
    rinv_S3     = gmx_invsqrt_pr(rsq_S3);
#else
    gmx_mm_invsqrt2_pd(rsq_S0, rsq_S1, &rinv_S0, &rinv_S1);
    gmx_mm_invsqrt2_pd(rsq_S2, rsq_S3, &rinv_S2, &rinv_S3);
#endif

#ifdef CALC_COULOMB
    /* Load parameters for j atom */
    jq_S        = gmx_load_pr(q+aj);
    qq_S0       = gmx_mul_pr(iq_S0, jq_S);
    qq_S1       = gmx_mul_pr(iq_S1, jq_S);
    qq_S2       = gmx_mul_pr(iq_S2, jq_S);
    qq_S3       = gmx_mul_pr(iq_S3, jq_S);
#endif

#ifdef CALC_LJ

#if !defined LJ_COMB_GEOM && !defined LJ_COMB_LB && !defined FIX_LJ_C
    load_lj_pair_params(nbfp0, type, aj, &c6_S0, &c12_S0);
    load_lj_pair_params(nbfp1, type, aj, &c6_S1, &c12_S1);
#ifndef HALF_LJ
    load_lj_pair_params(nbfp2, type, aj, &c6_S2, &c12_S2);
    load_lj_pair_params(nbfp3, type, aj, &c6_S3, &c12_S3);
#endif
#endif /* not defined any LJ rule */

#ifdef LJ_COMB_GEOM
    c6s_j_S     = gmx_load_pr(ljc+aj2+0);
    c12s_j_S    = gmx_load_pr(ljc+aj2+STRIDE);
    c6_S0       = gmx_mul_pr(c6s_S0, c6s_j_S );
    c6_S1       = gmx_mul_pr(c6s_S1, c6s_j_S );
#ifndef HALF_LJ
    c6_S2       = gmx_mul_pr(c6s_S2, c6s_j_S );
    c6_S3       = gmx_mul_pr(c6s_S3, c6s_j_S );
#endif
    c12_S0      = gmx_mul_pr(c12s_S0, c12s_j_S);
    c12_S1      = gmx_mul_pr(c12s_S1, c12s_j_S);
#ifndef HALF_LJ
    c12_S2      = gmx_mul_pr(c12s_S2, c12s_j_S);
    c12_S3      = gmx_mul_pr(c12s_S3, c12s_j_S);
#endif
#endif /* LJ_COMB_GEOM */

#ifdef LJ_COMB_LB
    hsig_j_S    = gmx_load_pr(ljc+aj2+0);
    seps_j_S    = gmx_load_pr(ljc+aj2+STRIDE);

    sig_S0      = gmx_add_pr(hsig_i_S0, hsig_j_S);
    sig_S1      = gmx_add_pr(hsig_i_S1, hsig_j_S);
    eps_S0      = gmx_mul_pr(seps_i_S0, seps_j_S);
    eps_S1      = gmx_mul_pr(seps_i_S1, seps_j_S);
#ifndef HALF_LJ
    sig_S2      = gmx_add_pr(hsig_i_S2, hsig_j_S);
    sig_S3      = gmx_add_pr(hsig_i_S3, hsig_j_S);
    eps_S2      = gmx_mul_pr(seps_i_S2, seps_j_S);
    eps_S3      = gmx_mul_pr(seps_i_S3, seps_j_S);
#endif
#endif /* LJ_COMB_LB */

#endif /* CALC_LJ */

#ifndef NBNXN_CUTOFF_USE_BLENDV
    rinv_S0     = gmx_blendzero_pr(rinv_S0, wco_S0);
    rinv_S1     = gmx_blendzero_pr(rinv_S1, wco_S1);
    rinv_S2     = gmx_blendzero_pr(rinv_S2, wco_S2);
    rinv_S3     = gmx_blendzero_pr(rinv_S3, wco_S3);
#else
    /* We only need to mask for the cut-off: blendv is faster */
    rinv_S0     = gmx_blendv_pr(rinv_S0, zero_S, gmx_sub_pr(rc2_S, rsq_S0));
    rinv_S1     = gmx_blendv_pr(rinv_S1, zero_S, gmx_sub_pr(rc2_S, rsq_S1));
    rinv_S2     = gmx_blendv_pr(rinv_S2, zero_S, gmx_sub_pr(rc2_S, rsq_S2));
    rinv_S3     = gmx_blendv_pr(rinv_S3, zero_S, gmx_sub_pr(rc2_S, rsq_S3));
#endif

    rinvsq_S0   = gmx_mul_pr(rinv_S0, rinv_S0);
    rinvsq_S1   = gmx_mul_pr(rinv_S1, rinv_S1);
    rinvsq_S2   = gmx_mul_pr(rinv_S2, rinv_S2);
    rinvsq_S3   = gmx_mul_pr(rinv_S3, rinv_S3);

#ifdef CALC_COULOMB
    /* Note that here we calculate force*r, not the usual force/r.
     * This allows avoiding masking the reaction-field contribution,
     * as frcoul is later multiplied by rinvsq which has been
     * masked with the cut-off check.
     */

#ifdef EXCL_FORCES
    /* Only add 1/r for non-excluded atom pairs */
    rinv_ex_S0  = gmx_blendzero_pr(rinv_S0, interact_S0);
    rinv_ex_S1  = gmx_blendzero_pr(rinv_S1, interact_S1);
    rinv_ex_S2  = gmx_blendzero_pr(rinv_S2, interact_S2);
    rinv_ex_S3  = gmx_blendzero_pr(rinv_S3, interact_S3);
#else
    /* No exclusion forces, we always need 1/r */
#define     rinv_ex_S0    rinv_S0
#define     rinv_ex_S1    rinv_S1
#define     rinv_ex_S2    rinv_S2
#define     rinv_ex_S3    rinv_S3
#endif

#ifdef CALC_COUL_RF
    /* Electrostatic interactions */
    frcoul_S0   = gmx_mul_pr(qq_S0, gmx_madd_pr(rsq_S0, mrc_3_S, rinv_ex_S0));
    frcoul_S1   = gmx_mul_pr(qq_S1, gmx_madd_pr(rsq_S1, mrc_3_S, rinv_ex_S1));
    frcoul_S2   = gmx_mul_pr(qq_S2, gmx_madd_pr(rsq_S2, mrc_3_S, rinv_ex_S2));
    frcoul_S3   = gmx_mul_pr(qq_S3, gmx_madd_pr(rsq_S3, mrc_3_S, rinv_ex_S3));

#ifdef CALC_ENERGIES
    vcoul_S0    = gmx_mul_pr(qq_S0, gmx_add_pr(rinv_ex_S0, gmx_add_pr(gmx_mul_pr(rsq_S0, hrc_3_S), moh_rc_S)));
    vcoul_S1    = gmx_mul_pr(qq_S1, gmx_add_pr(rinv_ex_S1, gmx_add_pr(gmx_mul_pr(rsq_S1, hrc_3_S), moh_rc_S)));
    vcoul_S2    = gmx_mul_pr(qq_S2, gmx_add_pr(rinv_ex_S2, gmx_add_pr(gmx_mul_pr(rsq_S2, hrc_3_S), moh_rc_S)));
    vcoul_S3    = gmx_mul_pr(qq_S3, gmx_add_pr(rinv_ex_S3, gmx_add_pr(gmx_mul_pr(rsq_S3, hrc_3_S), moh_rc_S)));
#endif
#endif

#ifdef CALC_COUL_EWALD
    /* We need to mask (or limit) rsq for the cut-off,
     * as large distances can cause an overflow in gmx_pmecorrF/V.
     */
#ifndef NBNXN_CUTOFF_USE_BLENDV
    brsq_S0     = gmx_mul_pr(beta2_S, gmx_blendzero_pr(rsq_S0, wco_S0));
    brsq_S1     = gmx_mul_pr(beta2_S, gmx_blendzero_pr(rsq_S1, wco_S1));
    brsq_S2     = gmx_mul_pr(beta2_S, gmx_blendzero_pr(rsq_S2, wco_S2));
    brsq_S3     = gmx_mul_pr(beta2_S, gmx_blendzero_pr(rsq_S3, wco_S3));
#else
    /* Strangely, putting mul on a separate line is slower (icc 13) */
    brsq_S0     = gmx_mul_pr(beta2_S, gmx_blendv_pr(rsq_S0, zero_S, gmx_sub_pr(rc2_S, rsq_S0)));
    brsq_S1     = gmx_mul_pr(beta2_S, gmx_blendv_pr(rsq_S1, zero_S, gmx_sub_pr(rc2_S, rsq_S1)));
    brsq_S2     = gmx_mul_pr(beta2_S, gmx_blendv_pr(rsq_S2, zero_S, gmx_sub_pr(rc2_S, rsq_S2)));
    brsq_S3     = gmx_mul_pr(beta2_S, gmx_blendv_pr(rsq_S3, zero_S, gmx_sub_pr(rc2_S, rsq_S3)));
#endif
    ewcorr_S0   = gmx_mul_pr(gmx_pmecorrF_pr(brsq_S0), beta_S);
    ewcorr_S1   = gmx_mul_pr(gmx_pmecorrF_pr(brsq_S1), beta_S);
    ewcorr_S2   = gmx_mul_pr(gmx_pmecorrF_pr(brsq_S2), beta_S);
    ewcorr_S3   = gmx_mul_pr(gmx_pmecorrF_pr(brsq_S3), beta_S);
    frcoul_S0   = gmx_mul_pr(qq_S0, gmx_madd_pr(ewcorr_S0, brsq_S0, rinv_ex_S0));
    frcoul_S1   = gmx_mul_pr(qq_S1, gmx_madd_pr(ewcorr_S1, brsq_S1, rinv_ex_S1));
    frcoul_S2   = gmx_mul_pr(qq_S2, gmx_madd_pr(ewcorr_S2, brsq_S2, rinv_ex_S2));
    frcoul_S3   = gmx_mul_pr(qq_S3, gmx_madd_pr(ewcorr_S3, brsq_S3, rinv_ex_S3));

#ifdef CALC_ENERGIES
    vc_sub_S0   = gmx_mul_pr(gmx_pmecorrV_pr(brsq_S0), beta_S);
    vc_sub_S1   = gmx_mul_pr(gmx_pmecorrV_pr(brsq_S1), beta_S);
    vc_sub_S2   = gmx_mul_pr(gmx_pmecorrV_pr(brsq_S2), beta_S);
    vc_sub_S3   = gmx_mul_pr(gmx_pmecorrV_pr(brsq_S3), beta_S);
#endif

#endif /* CALC_COUL_EWALD */

#ifdef CALC_COUL_TAB
    /* Electrostatic interactions */
    r_S0        = gmx_mul_pr(rsq_S0, rinv_S0);
    r_S1        = gmx_mul_pr(rsq_S1, rinv_S1);
    r_S2        = gmx_mul_pr(rsq_S2, rinv_S2);
    r_S3        = gmx_mul_pr(rsq_S3, rinv_S3);
    /* Convert r to scaled table units */
    rs_S0       = gmx_mul_pr(r_S0, invtsp_S);
    rs_S1       = gmx_mul_pr(r_S1, invtsp_S);
    rs_S2       = gmx_mul_pr(r_S2, invtsp_S);
    rs_S3       = gmx_mul_pr(r_S3, invtsp_S);
    /* Truncate scaled r to an int */
    ti_S0       = gmx_cvttpr_epi32(rs_S0);
    ti_S1       = gmx_cvttpr_epi32(rs_S1);
    ti_S2       = gmx_cvttpr_epi32(rs_S2);
    ti_S3       = gmx_cvttpr_epi32(rs_S3);
#ifdef GMX_SIMD_HAVE_FLOOR
    /* SSE4.1 floor is faster than gmx_cvtepi32_ps int->float cast */
    rf_S0       = gmx_floor_pr(rs_S0);
    rf_S1       = gmx_floor_pr(rs_S1);
    rf_S2       = gmx_floor_pr(rs_S2);
    rf_S3       = gmx_floor_pr(rs_S3);
#else
    rf_S0       = gmx_cvtepi32_pr(ti_S0);
    rf_S1       = gmx_cvtepi32_pr(ti_S1);
    rf_S2       = gmx_cvtepi32_pr(ti_S2);
    rf_S3       = gmx_cvtepi32_pr(ti_S3);
#endif
    frac_S0     = gmx_sub_pr(rs_S0, rf_S0);
    frac_S1     = gmx_sub_pr(rs_S1, rf_S1);
    frac_S2     = gmx_sub_pr(rs_S2, rf_S2);
    frac_S3     = gmx_sub_pr(rs_S3, rf_S3);

    /* Load and interpolate table forces and possibly energies.
     * Force and energy can be combined in one table, stride 4: FDV0
     * or in two separate tables with stride 1: F and V
     * Currently single precision uses FDV0, double F and V.
     */
#ifndef CALC_ENERGIES
    load_table_f(tab_coul_F, ti_S0, ti0, &ctab0_S0, &ctab1_S0);
    load_table_f(tab_coul_F, ti_S1, ti1, &ctab0_S1, &ctab1_S1);
    load_table_f(tab_coul_F, ti_S2, ti2, &ctab0_S2, &ctab1_S2);
    load_table_f(tab_coul_F, ti_S3, ti3, &ctab0_S3, &ctab1_S3);
#else
#ifdef TAB_FDV0
    load_table_f_v(tab_coul_F, ti_S0, ti0, &ctab0_S0, &ctab1_S0, &ctabv_S0);
    load_table_f_v(tab_coul_F, ti_S1, ti1, &ctab0_S1, &ctab1_S1, &ctabv_S1);
    load_table_f_v(tab_coul_F, ti_S2, ti2, &ctab0_S2, &ctab1_S2, &ctabv_S2);
    load_table_f_v(tab_coul_F, ti_S3, ti3, &ctab0_S3, &ctab1_S3, &ctabv_S3);
#else
    load_table_f_v(tab_coul_F, tab_coul_V, ti_S0, ti0, &ctab0_S0, &ctab1_S0, &ctabv_S0);
    load_table_f_v(tab_coul_F, tab_coul_V, ti_S1, ti1, &ctab0_S1, &ctab1_S1, &ctabv_S1);
    load_table_f_v(tab_coul_F, tab_coul_V, ti_S2, ti2, &ctab0_S2, &ctab1_S2, &ctabv_S2);
    load_table_f_v(tab_coul_F, tab_coul_V, ti_S3, ti3, &ctab0_S3, &ctab1_S3, &ctabv_S3);
#endif
#endif
    fsub_S0     = gmx_add_pr(ctab0_S0, gmx_mul_pr(frac_S0, ctab1_S0));
    fsub_S1     = gmx_add_pr(ctab0_S1, gmx_mul_pr(frac_S1, ctab1_S1));
    fsub_S2     = gmx_add_pr(ctab0_S2, gmx_mul_pr(frac_S2, ctab1_S2));
    fsub_S3     = gmx_add_pr(ctab0_S3, gmx_mul_pr(frac_S3, ctab1_S3));
    frcoul_S0   = gmx_mul_pr(qq_S0, gmx_sub_pr(rinv_ex_S0, gmx_mul_pr(fsub_S0, r_S0)));
    frcoul_S1   = gmx_mul_pr(qq_S1, gmx_sub_pr(rinv_ex_S1, gmx_mul_pr(fsub_S1, r_S1)));
    frcoul_S2   = gmx_mul_pr(qq_S2, gmx_sub_pr(rinv_ex_S2, gmx_mul_pr(fsub_S2, r_S2)));
    frcoul_S3   = gmx_mul_pr(qq_S3, gmx_sub_pr(rinv_ex_S3, gmx_mul_pr(fsub_S3, r_S3)));

#ifdef CALC_ENERGIES
    vc_sub_S0   = gmx_add_pr(ctabv_S0, gmx_mul_pr(gmx_mul_pr(mhalfsp_S, frac_S0), gmx_add_pr(ctab0_S0, fsub_S0)));
    vc_sub_S1   = gmx_add_pr(ctabv_S1, gmx_mul_pr(gmx_mul_pr(mhalfsp_S, frac_S1), gmx_add_pr(ctab0_S1, fsub_S1)));
    vc_sub_S2   = gmx_add_pr(ctabv_S2, gmx_mul_pr(gmx_mul_pr(mhalfsp_S, frac_S2), gmx_add_pr(ctab0_S2, fsub_S2)));
    vc_sub_S3   = gmx_add_pr(ctabv_S3, gmx_mul_pr(gmx_mul_pr(mhalfsp_S, frac_S3), gmx_add_pr(ctab0_S3, fsub_S3)));
#endif
#endif /* CALC_COUL_TAB */

#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
#ifndef NO_SHIFT_EWALD
    /* Add Ewald potential shift to vc_sub for convenience */
#ifdef CHECK_EXCLS
    vc_sub_S0   = gmx_add_pr(vc_sub_S0, gmx_blendzero_pr(sh_ewald_S, interact_S0));
    vc_sub_S1   = gmx_add_pr(vc_sub_S1, gmx_blendzero_pr(sh_ewald_S, interact_S1));
    vc_sub_S2   = gmx_add_pr(vc_sub_S2, gmx_blendzero_pr(sh_ewald_S, interact_S2));
    vc_sub_S3   = gmx_add_pr(vc_sub_S3, gmx_blendzero_pr(sh_ewald_S, interact_S3));
#else
    vc_sub_S0   = gmx_add_pr(vc_sub_S0, sh_ewald_S);
    vc_sub_S1   = gmx_add_pr(vc_sub_S1, sh_ewald_S);
    vc_sub_S2   = gmx_add_pr(vc_sub_S2, sh_ewald_S);
    vc_sub_S3   = gmx_add_pr(vc_sub_S3, sh_ewald_S);
#endif
#endif

    vcoul_S0    = gmx_mul_pr(qq_S0, gmx_sub_pr(rinv_ex_S0, vc_sub_S0));
    vcoul_S1    = gmx_mul_pr(qq_S1, gmx_sub_pr(rinv_ex_S1, vc_sub_S1));
    vcoul_S2    = gmx_mul_pr(qq_S2, gmx_sub_pr(rinv_ex_S2, vc_sub_S2));
    vcoul_S3    = gmx_mul_pr(qq_S3, gmx_sub_pr(rinv_ex_S3, vc_sub_S3));

#endif

#ifdef CALC_ENERGIES
    /* Mask energy for cut-off and diagonal */
    vcoul_S0    = gmx_blendzero_pr(vcoul_S0, wco_S0);
    vcoul_S1    = gmx_blendzero_pr(vcoul_S1, wco_S1);
    vcoul_S2    = gmx_blendzero_pr(vcoul_S2, wco_S2);
    vcoul_S3    = gmx_blendzero_pr(vcoul_S3, wco_S3);
#endif

#endif /* CALC_COULOMB */

#ifdef CALC_LJ
    /* Lennard-Jones interaction */

#ifdef VDW_CUTOFF_CHECK
    wco_vdw_S0  = gmx_cmplt_pr(rsq_S0, rcvdw2_S);
    wco_vdw_S1  = gmx_cmplt_pr(rsq_S1, rcvdw2_S);
#ifndef HALF_LJ
    wco_vdw_S2  = gmx_cmplt_pr(rsq_S2, rcvdw2_S);
    wco_vdw_S3  = gmx_cmplt_pr(rsq_S3, rcvdw2_S);
#endif
#else
    /* Same cut-off for Coulomb and VdW, reuse the registers */
#define     wco_vdw_S0    wco_S0
#define     wco_vdw_S1    wco_S1
#define     wco_vdw_S2    wco_S2
#define     wco_vdw_S3    wco_S3
#endif

#ifndef LJ_COMB_LB
    rinvsix_S0  = gmx_mul_pr(rinvsq_S0, gmx_mul_pr(rinvsq_S0, rinvsq_S0));
    rinvsix_S1  = gmx_mul_pr(rinvsq_S1, gmx_mul_pr(rinvsq_S1, rinvsq_S1));
#ifdef EXCL_FORCES
    rinvsix_S0  = gmx_blendzero_pr(rinvsix_S0, interact_S0);
    rinvsix_S1  = gmx_blendzero_pr(rinvsix_S1, interact_S1);
#endif
#ifndef HALF_LJ
    rinvsix_S2  = gmx_mul_pr(rinvsq_S2, gmx_mul_pr(rinvsq_S2, rinvsq_S2));
    rinvsix_S3  = gmx_mul_pr(rinvsq_S3, gmx_mul_pr(rinvsq_S3, rinvsq_S3));
#ifdef EXCL_FORCES
    rinvsix_S2  = gmx_blendzero_pr(rinvsix_S2, interact_S2);
    rinvsix_S3  = gmx_blendzero_pr(rinvsix_S3, interact_S3);
#endif
#endif
#ifdef VDW_CUTOFF_CHECK
    rinvsix_S0  = gmx_blendzero_pr(rinvsix_S0, wco_vdw_S0);
    rinvsix_S1  = gmx_blendzero_pr(rinvsix_S1, wco_vdw_S1);
#ifndef HALF_LJ
    rinvsix_S2  = gmx_blendzero_pr(rinvsix_S2, wco_vdw_S2);
    rinvsix_S3  = gmx_blendzero_pr(rinvsix_S3, wco_vdw_S3);
#endif
#endif
    FrLJ6_S0    = gmx_mul_pr(c6_S0, rinvsix_S0);
    FrLJ6_S1    = gmx_mul_pr(c6_S1, rinvsix_S1);
#ifndef HALF_LJ
    FrLJ6_S2    = gmx_mul_pr(c6_S2, rinvsix_S2);
    FrLJ6_S3    = gmx_mul_pr(c6_S3, rinvsix_S3);
#endif
    FrLJ12_S0   = gmx_mul_pr(c12_S0, gmx_mul_pr(rinvsix_S0, rinvsix_S0));
    FrLJ12_S1   = gmx_mul_pr(c12_S1, gmx_mul_pr(rinvsix_S1, rinvsix_S1));
#ifndef HALF_LJ
    FrLJ12_S2   = gmx_mul_pr(c12_S2, gmx_mul_pr(rinvsix_S2, rinvsix_S2));
    FrLJ12_S3   = gmx_mul_pr(c12_S3, gmx_mul_pr(rinvsix_S3, rinvsix_S3));
#endif
#endif /* not LJ_COMB_LB */

#ifdef LJ_COMB_LB
    sir_S0      = gmx_mul_pr(sig_S0, rinv_S0);
    sir_S1      = gmx_mul_pr(sig_S1, rinv_S1);
#ifndef HALF_LJ
    sir_S2      = gmx_mul_pr(sig_S2, rinv_S2);
    sir_S3      = gmx_mul_pr(sig_S3, rinv_S3);
#endif
    sir2_S0     = gmx_mul_pr(sir_S0, sir_S0);
    sir2_S1     = gmx_mul_pr(sir_S1, sir_S1);
#ifndef HALF_LJ
    sir2_S2     = gmx_mul_pr(sir_S2, sir_S2);
    sir2_S3     = gmx_mul_pr(sir_S3, sir_S3);
#endif
    sir6_S0     = gmx_mul_pr(sir2_S0, gmx_mul_pr(sir2_S0, sir2_S0));
    sir6_S1     = gmx_mul_pr(sir2_S1, gmx_mul_pr(sir2_S1, sir2_S1));
#ifdef EXCL_FORCES
    sir6_S0     = gmx_blendzero_pr(sir6_S0, interact_S0);
    sir6_S1     = gmx_blendzero_pr(sir6_S1, interact_S1);
#endif
#ifndef HALF_LJ
    sir6_S2     = gmx_mul_pr(sir2_S2, gmx_mul_pr(sir2_S2, sir2_S2));
    sir6_S3     = gmx_mul_pr(sir2_S3, gmx_mul_pr(sir2_S3, sir2_S3));
#ifdef EXCL_FORCES
    sir6_S2     = gmx_blendzero_pr(sir6_S2, interact_S2);
    sir6_S3     = gmx_blendzero_pr(sir6_S3, interact_S3);
#endif
#endif
#ifdef VDW_CUTOFF_CHECK
    sir6_S0     = gmx_blendzero_pr(sir6_S0, wco_vdw_S0);
    sir6_S1     = gmx_blendzero_pr(sir6_S1, wco_vdw_S1);
#ifndef HALF_LJ
    sir6_S2     = gmx_blendzero_pr(sir6_S2, wco_vdw_S2);
    sir6_S3     = gmx_blendzero_pr(sir6_S3, wco_vdw_S3);
#endif
#endif
    FrLJ6_S0    = gmx_mul_pr(eps_S0, sir6_S0);
    FrLJ6_S1    = gmx_mul_pr(eps_S1, sir6_S1);
#ifndef HALF_LJ
    FrLJ6_S2    = gmx_mul_pr(eps_S2, sir6_S2);
    FrLJ6_S3    = gmx_mul_pr(eps_S3, sir6_S3);
#endif
    FrLJ12_S0   = gmx_mul_pr(FrLJ6_S0, sir6_S0);
    FrLJ12_S1   = gmx_mul_pr(FrLJ6_S1, sir6_S1);
#ifndef HALF_LJ
    FrLJ12_S2   = gmx_mul_pr(FrLJ6_S2, sir6_S2);
    FrLJ12_S3   = gmx_mul_pr(FrLJ6_S3, sir6_S3);
#endif
#if defined CALC_ENERGIES
    /* We need C6 and C12 to calculate the LJ potential shift */
    sig2_S0     = gmx_mul_pr(sig_S0, sig_S0);
    sig2_S1     = gmx_mul_pr(sig_S1, sig_S1);
#ifndef HALF_LJ
    sig2_S2     = gmx_mul_pr(sig_S2, sig_S2);
    sig2_S3     = gmx_mul_pr(sig_S3, sig_S3);
#endif
    sig6_S0     = gmx_mul_pr(sig2_S0, gmx_mul_pr(sig2_S0, sig2_S0));
    sig6_S1     = gmx_mul_pr(sig2_S1, gmx_mul_pr(sig2_S1, sig2_S1));
#ifndef HALF_LJ
    sig6_S2     = gmx_mul_pr(sig2_S2, gmx_mul_pr(sig2_S2, sig2_S2));
    sig6_S3     = gmx_mul_pr(sig2_S3, gmx_mul_pr(sig2_S3, sig2_S3));
#endif
    c6_S0       = gmx_mul_pr(eps_S0, sig6_S0);
    c6_S1       = gmx_mul_pr(eps_S1, sig6_S1);
#ifndef HALF_LJ
    c6_S2       = gmx_mul_pr(eps_S2, sig6_S2);
    c6_S3       = gmx_mul_pr(eps_S3, sig6_S3);
#endif
    c12_S0      = gmx_mul_pr(c6_S0, sig6_S0);
    c12_S1      = gmx_mul_pr(c6_S1, sig6_S1);
#ifndef HALF_LJ
    c12_S2      = gmx_mul_pr(c6_S2, sig6_S2);
    c12_S3      = gmx_mul_pr(c6_S3, sig6_S3);
#endif
#endif
#endif /* LJ_COMB_LB */

#endif /* CALC_LJ */

#ifdef CALC_ENERGIES
#ifdef ENERGY_GROUPS
    /* Extract the group pair index per j pair.
     * Energy groups are stored per i-cluster, so things get
     * complicated when the i- and j-cluster size don't match.
     */
    {
        int egps_j;
#if UNROLLJ == 2
        egps_j    = nbat->energrp[cj>>1];
        egp_jj[0] = ((egps_j >> ((cj & 1)*egps_jshift)) & egps_jmask)*egps_jstride;
#else
        /* We assume UNROLLI <= UNROLLJ */
        int jdi;
        for (jdi = 0; jdi < UNROLLJ/UNROLLI; jdi++)
        {
            int jj;
            egps_j = nbat->energrp[cj*(UNROLLJ/UNROLLI)+jdi];
            for (jj = 0; jj < (UNROLLI/2); jj++)
            {
                egp_jj[jdi*(UNROLLI/2)+jj] = ((egps_j >> (jj*egps_jshift)) & egps_jmask)*egps_jstride;
            }
        }
#endif
    }
#endif

#ifdef CALC_COULOMB
#ifndef ENERGY_GROUPS
    vctot_S      = gmx_add_pr(vctot_S, gmx_sum4_pr(vcoul_S0, vcoul_S1, vcoul_S2, vcoul_S3));
#else
    add_ener_grp(vcoul_S0, vctp[0], egp_jj);
    add_ener_grp(vcoul_S1, vctp[1], egp_jj);
    add_ener_grp(vcoul_S2, vctp[2], egp_jj);
    add_ener_grp(vcoul_S3, vctp[3], egp_jj);
#endif
#endif

#ifdef CALC_LJ
    /* Calculate the LJ energies */
    VLJ6_S0     = gmx_mul_pr(sixth_S, gmx_sub_pr(FrLJ6_S0, gmx_mul_pr(c6_S0, sh_invrc6_S)));
    VLJ6_S1     = gmx_mul_pr(sixth_S, gmx_sub_pr(FrLJ6_S1, gmx_mul_pr(c6_S1, sh_invrc6_S)));
#ifndef HALF_LJ
    VLJ6_S2     = gmx_mul_pr(sixth_S, gmx_sub_pr(FrLJ6_S2, gmx_mul_pr(c6_S2, sh_invrc6_S)));
    VLJ6_S3     = gmx_mul_pr(sixth_S, gmx_sub_pr(FrLJ6_S3, gmx_mul_pr(c6_S3, sh_invrc6_S)));
#endif
    VLJ12_S0    = gmx_mul_pr(twelveth_S, gmx_sub_pr(FrLJ12_S0, gmx_mul_pr(c12_S0, sh_invrc12_S)));
    VLJ12_S1    = gmx_mul_pr(twelveth_S, gmx_sub_pr(FrLJ12_S1, gmx_mul_pr(c12_S1, sh_invrc12_S)));
#ifndef HALF_LJ
    VLJ12_S2    = gmx_mul_pr(twelveth_S, gmx_sub_pr(FrLJ12_S2, gmx_mul_pr(c12_S2, sh_invrc12_S)));
    VLJ12_S3    = gmx_mul_pr(twelveth_S, gmx_sub_pr(FrLJ12_S3, gmx_mul_pr(c12_S3, sh_invrc12_S)));
#endif

    VLJ_S0      = gmx_sub_pr(VLJ12_S0, VLJ6_S0);
    VLJ_S1      = gmx_sub_pr(VLJ12_S1, VLJ6_S1);
#ifndef HALF_LJ
    VLJ_S2      = gmx_sub_pr(VLJ12_S2, VLJ6_S2);
    VLJ_S3      = gmx_sub_pr(VLJ12_S3, VLJ6_S3);
#endif
    /* The potential shift should be removed for pairs beyond cut-off */
    VLJ_S0      = gmx_blendzero_pr(VLJ_S0, wco_vdw_S0);
    VLJ_S1      = gmx_blendzero_pr(VLJ_S1, wco_vdw_S1);
#ifndef HALF_LJ
    VLJ_S2      = gmx_blendzero_pr(VLJ_S2, wco_vdw_S2);
    VLJ_S3      = gmx_blendzero_pr(VLJ_S3, wco_vdw_S3);
#endif
#ifdef CHECK_EXCLS
    /* The potential shift should be removed for excluded pairs */
    VLJ_S0      = gmx_blendzero_pr(VLJ_S0, interact_S0);
    VLJ_S1      = gmx_blendzero_pr(VLJ_S1, interact_S1);
#ifndef HALF_LJ
    VLJ_S2      = gmx_blendzero_pr(VLJ_S2, interact_S2);
    VLJ_S3      = gmx_blendzero_pr(VLJ_S3, interact_S3);
#endif
#endif
#ifndef ENERGY_GROUPS
    Vvdwtot_S   = gmx_add_pr(Vvdwtot_S,
#ifndef HALF_LJ
                             gmx_sum4_pr(VLJ_S0, VLJ_S1, VLJ_S2, VLJ_S3)
#else
                             gmx_add_pr(VLJ_S0, VLJ_S1)
#endif
                             );
#else
    add_ener_grp(VLJ_S0, vvdwtp[0], egp_jj);
    add_ener_grp(VLJ_S1, vvdwtp[1], egp_jj);
#ifndef HALF_LJ
    add_ener_grp(VLJ_S2, vvdwtp[2], egp_jj);
    add_ener_grp(VLJ_S3, vvdwtp[3], egp_jj);
#endif
#endif
#endif /* CALC_LJ */
#endif /* CALC_ENERGIES */

#ifdef CALC_LJ
    fscal_S0    = gmx_mul_pr(rinvsq_S0,
#ifdef CALC_COULOMB
                             gmx_add_pr(frcoul_S0,
#else
                             (
#endif
                                        gmx_sub_pr(FrLJ12_S0, FrLJ6_S0)));
    fscal_S1    = gmx_mul_pr(rinvsq_S1,
#ifdef CALC_COULOMB
                             gmx_add_pr(frcoul_S1,
#else
                             (
#endif
                                        gmx_sub_pr(FrLJ12_S1, FrLJ6_S1)));
#else
    fscal_S0    = gmx_mul_pr(rinvsq_S0, frcoul_S0);
    fscal_S1    = gmx_mul_pr(rinvsq_S1, frcoul_S1);
#endif /* CALC_LJ */
#if defined CALC_LJ && !defined HALF_LJ
    fscal_S2    = gmx_mul_pr(rinvsq_S2,
#ifdef CALC_COULOMB
                             gmx_add_pr(frcoul_S2,
#else
                             (
#endif
                                        gmx_sub_pr(FrLJ12_S2, FrLJ6_S2)));
    fscal_S3    = gmx_mul_pr(rinvsq_S3,
#ifdef CALC_COULOMB
                             gmx_add_pr(frcoul_S3,
#else
                             (
#endif
                                        gmx_sub_pr(FrLJ12_S3, FrLJ6_S3)));
#else
    /* Atom 2 and 3 don't have LJ, so only add Coulomb forces */
    fscal_S2    = gmx_mul_pr(rinvsq_S2, frcoul_S2);
    fscal_S3    = gmx_mul_pr(rinvsq_S3, frcoul_S3);
#endif

    /* Calculate temporary vectorial force */
    tx_S0       = gmx_mul_pr(fscal_S0, dx_S0);
    tx_S1       = gmx_mul_pr(fscal_S1, dx_S1);
    tx_S2       = gmx_mul_pr(fscal_S2, dx_S2);
    tx_S3       = gmx_mul_pr(fscal_S3, dx_S3);
    ty_S0       = gmx_mul_pr(fscal_S0, dy_S0);
    ty_S1       = gmx_mul_pr(fscal_S1, dy_S1);
    ty_S2       = gmx_mul_pr(fscal_S2, dy_S2);
    ty_S3       = gmx_mul_pr(fscal_S3, dy_S3);
    tz_S0       = gmx_mul_pr(fscal_S0, dz_S0);
    tz_S1       = gmx_mul_pr(fscal_S1, dz_S1);
    tz_S2       = gmx_mul_pr(fscal_S2, dz_S2);
    tz_S3       = gmx_mul_pr(fscal_S3, dz_S3);

    /* Increment i atom force */
    fix_S0      = gmx_add_pr(fix_S0, tx_S0);
    fix_S1      = gmx_add_pr(fix_S1, tx_S1);
    fix_S2      = gmx_add_pr(fix_S2, tx_S2);
    fix_S3      = gmx_add_pr(fix_S3, tx_S3);
    fiy_S0      = gmx_add_pr(fiy_S0, ty_S0);
    fiy_S1      = gmx_add_pr(fiy_S1, ty_S1);
    fiy_S2      = gmx_add_pr(fiy_S2, ty_S2);
    fiy_S3      = gmx_add_pr(fiy_S3, ty_S3);
    fiz_S0      = gmx_add_pr(fiz_S0, tz_S0);
    fiz_S1      = gmx_add_pr(fiz_S1, tz_S1);
    fiz_S2      = gmx_add_pr(fiz_S2, tz_S2);
    fiz_S3      = gmx_add_pr(fiz_S3, tz_S3);

    /* Decrement j atom force */
    gmx_store_pr(f+ajx,
                 gmx_sub_pr( gmx_load_pr(f+ajx), gmx_sum4_pr(tx_S0, tx_S1, tx_S2, tx_S3) ));
    gmx_store_pr(f+ajy,
                 gmx_sub_pr( gmx_load_pr(f+ajy), gmx_sum4_pr(ty_S0, ty_S1, ty_S2, ty_S3) ));
    gmx_store_pr(f+ajz,
                 gmx_sub_pr( gmx_load_pr(f+ajz), gmx_sum4_pr(tz_S0, tz_S1, tz_S2, tz_S3) ));
}

#undef  rinv_ex_S0
#undef  rinv_ex_S1
#undef  rinv_ex_S2
#undef  rinv_ex_S3

#undef  wco_vdw_S0
#undef  wco_vdw_S1
#undef  wco_vdw_S2
#undef  wco_vdw_S3

#undef  NBNXN_CUTOFF_USE_BLENDV

#undef  EXCL_FORCES
