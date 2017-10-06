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

/* Doxygen gets confused (buggy) about the block in this file in combination with
 * the  namespace prefix, and thinks store is documented here.
 * This will solve itself with the second-generation nbnxn kernels, so for now
 * we just tell Doxygen to stay out.
 */
#ifndef DOXYGEN

/* This is the innermost loop contents for the 4 x N atom simd kernel.
 * This flavor of the kernel calculates interactions of 4 i-atoms
 * with N j-atoms stored in N wide simd registers.
 */

/* When calculating RF or Ewald interactions we calculate the electrostatic/LJ
 * forces on excluded atom pairs here in the non-bonded loops.
 * But when energies and/or virial is required we calculate them
 * separately to as then it is easier to separate the energy and virial
 * contributions.
 */
#if defined CHECK_EXCLS && (defined CALC_COULOMB || defined LJ_EWALD_GEOM)
#define EXCL_FORCES
#endif

{
    int            cj, ajx, ajy, ajz;
    int gmx_unused aj;

#ifdef ENERGY_GROUPS
    /* Energy group indices for two atoms packed into one int */
    int        egp_jj[UNROLLJ/2];
#endif

#ifdef CHECK_EXCLS
    /* Interaction (non-exclusion) mask of all 1's or 0's */
    SimdBool  interact_S0;
    SimdBool  interact_S1;
    SimdBool  interact_S2;
    SimdBool  interact_S3;
#endif

    SimdReal  jx_S, jy_S, jz_S;
    SimdReal  dx_S0, dy_S0, dz_S0;
    SimdReal  dx_S1, dy_S1, dz_S1;
    SimdReal  dx_S2, dy_S2, dz_S2;
    SimdReal  dx_S3, dy_S3, dz_S3;
    SimdReal  tx_S0, ty_S0, tz_S0;
    SimdReal  tx_S1, ty_S1, tz_S1;
    SimdReal  tx_S2, ty_S2, tz_S2;
    SimdReal  tx_S3, ty_S3, tz_S3;
    SimdReal  rsq_S0, rinv_S0, rinvsq_S0;
    SimdReal  rsq_S1, rinv_S1, rinvsq_S1;
    SimdReal  rsq_S2, rinv_S2, rinvsq_S2;
    SimdReal  rsq_S3, rinv_S3, rinvsq_S3;

    /* wco: within cut-off, mask of all 1's or 0's */
    SimdBool  wco_S0;
    SimdBool  wco_S1;
    SimdBool  wco_S2;
    SimdBool  wco_S3;

#ifdef VDW_CUTOFF_CHECK
    SimdBool  wco_vdw_S0;
    SimdBool  wco_vdw_S1;
#ifndef HALF_LJ
    SimdBool  wco_vdw_S2;
    SimdBool  wco_vdw_S3;
#endif
#endif

#if (defined CALC_COULOMB && defined CALC_COUL_TAB) || defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
    SimdReal r_S0;
    SimdReal r_S1;
#if (defined CALC_COULOMB && defined CALC_COUL_TAB) || !defined HALF_LJ
    SimdReal r_S2;
    SimdReal r_S3;
#endif
#endif

#if defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
    SimdReal  rsw_S0, rsw2_S0;
    SimdReal  rsw_S1, rsw2_S1;
#ifndef HALF_LJ
    SimdReal  rsw_S2, rsw2_S2;
    SimdReal  rsw_S3, rsw2_S3;
#endif
#endif

#ifdef CALC_COULOMB
#ifdef CHECK_EXCLS
    /* 1/r masked with the interaction mask */
    SimdReal  rinv_ex_S0;
    SimdReal  rinv_ex_S1;
    SimdReal  rinv_ex_S2;
    SimdReal  rinv_ex_S3;
#endif
    SimdReal  jq_S;
    SimdReal  qq_S0;
    SimdReal  qq_S1;
    SimdReal  qq_S2;
    SimdReal  qq_S3;
#ifdef CALC_COUL_TAB
    /* The force (PME mesh force) we need to subtract from 1/r^2 */
    SimdReal  fsub_S0;
    SimdReal  fsub_S1;
    SimdReal  fsub_S2;
    SimdReal  fsub_S3;
#endif
#ifdef CALC_COUL_EWALD
    SimdReal  brsq_S0, brsq_S1, brsq_S2, brsq_S3;
    SimdReal  ewcorr_S0, ewcorr_S1, ewcorr_S2, ewcorr_S3;
#endif

    /* frcoul = (1/r - fsub)*r */
    SimdReal  frcoul_S0;
    SimdReal  frcoul_S1;
    SimdReal  frcoul_S2;
    SimdReal  frcoul_S3;
#ifdef CALC_COUL_TAB
    /* For tables: r, rs=r/sp, rf=floor(rs), frac=rs-rf */
    SimdReal         rs_S0, rf_S0, frac_S0;
    SimdReal         rs_S1, rf_S1, frac_S1;
    SimdReal         rs_S2, rf_S2, frac_S2;
    SimdReal         rs_S3, rf_S3, frac_S3;
    /* Table index: rs truncated to an int */
    SimdInt32        ti_S0, ti_S1, ti_S2, ti_S3;
    /* Linear force table values */
    SimdReal         ctab0_S0, ctab1_S0;
    SimdReal         ctab0_S1, ctab1_S1;
    SimdReal         ctab0_S2, ctab1_S2;
    SimdReal         ctab0_S3, ctab1_S3;
#ifdef CALC_ENERGIES
    /* Quadratic energy table value */
    SimdReal  ctabv_S0, dum_S0;
    SimdReal  ctabv_S1, dum_S1;
    SimdReal  ctabv_S2, dum_S2;
    SimdReal  ctabv_S3, dum_S3;
#endif
#endif
#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
    /* The potential (PME mesh) we need to subtract from 1/r */
    SimdReal  vc_sub_S0;
    SimdReal  vc_sub_S1;
    SimdReal  vc_sub_S2;
    SimdReal  vc_sub_S3;
#endif
#ifdef CALC_ENERGIES
    /* Electrostatic potential */
    SimdReal  vcoul_S0;
    SimdReal  vcoul_S1;
    SimdReal  vcoul_S2;
    SimdReal  vcoul_S3;
#endif
#endif
    /* The force times 1/r */
    SimdReal  fscal_S0;
    SimdReal  fscal_S1;
    SimdReal  fscal_S2;
    SimdReal  fscal_S3;

#ifdef CALC_LJ
#ifdef LJ_COMB_LB
    /* LJ sigma_j/2 and sqrt(epsilon_j) */
    SimdReal  hsig_j_S, seps_j_S;
    /* LJ sigma_ij and epsilon_ij */
    SimdReal  sig_S0, eps_S0;
    SimdReal  sig_S1, eps_S1;
#ifndef HALF_LJ
    SimdReal  sig_S2, eps_S2;
    SimdReal  sig_S3, eps_S3;
#endif
#ifdef CALC_ENERGIES
    SimdReal  sig2_S0, sig6_S0;
    SimdReal  sig2_S1, sig6_S1;
#ifndef HALF_LJ
    SimdReal  sig2_S2, sig6_S2;
    SimdReal  sig2_S3, sig6_S3;
#endif
#endif /* LJ_COMB_LB */
#endif /* CALC_LJ */

#ifdef LJ_COMB_GEOM
    SimdReal  c6s_j_S, c12s_j_S;
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB || defined LJ_EWALD_GEOM
    /* Index for loading LJ parameters, complicated when interleaving */
    int         aj2;
#endif

    /* Intermediate variables for LJ calculation */
#ifndef LJ_COMB_LB
    SimdReal  rinvsix_S0;
    SimdReal  rinvsix_S1;
#ifndef HALF_LJ
    SimdReal  rinvsix_S2;
    SimdReal  rinvsix_S3;
#endif
#endif
#ifdef LJ_COMB_LB
    SimdReal  sir_S0, sir2_S0, sir6_S0;
    SimdReal  sir_S1, sir2_S1, sir6_S1;
#ifndef HALF_LJ
    SimdReal  sir_S2, sir2_S2, sir6_S2;
    SimdReal  sir_S3, sir2_S3, sir6_S3;
#endif
#endif

    SimdReal  FrLJ6_S0, FrLJ12_S0, frLJ_S0;
    SimdReal  FrLJ6_S1, FrLJ12_S1, frLJ_S1;
#ifndef HALF_LJ
    SimdReal  FrLJ6_S2, FrLJ12_S2, frLJ_S2;
    SimdReal  FrLJ6_S3, FrLJ12_S3, frLJ_S3;
#endif
#endif /* CALC_LJ */

    /* j-cluster index */
    cj            = l_cj[cjind].cj;

    /* Atom indices (of the first atom in the cluster) */
    aj            = cj*UNROLLJ;
#if defined CALC_LJ && (defined LJ_COMB_GEOM || defined LJ_COMB_LB || defined LJ_EWALD_GEOM)
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
    gmx_load_simd_4xn_interactions(l_cj[cjind].excl,
                                   filter_S0, filter_S1,
                                   filter_S2, filter_S3,
                                   nbat->simd_interaction_array,
                                   &interact_S0, &interact_S1,
                                   &interact_S2, &interact_S3);
#endif /* CHECK_EXCLS */

    /* load j atom coordinates */
    jx_S        = load<SimdReal>(x+ajx);
    jy_S        = load<SimdReal>(x+ajy);
    jz_S        = load<SimdReal>(x+ajz);

    /* Calculate distance */
    dx_S0       = ix_S0 - jx_S;
    dy_S0       = iy_S0 - jy_S;
    dz_S0       = iz_S0 - jz_S;
    dx_S1       = ix_S1 - jx_S;
    dy_S1       = iy_S1 - jy_S;
    dz_S1       = iz_S1 - jz_S;
    dx_S2       = ix_S2 - jx_S;
    dy_S2       = iy_S2 - jy_S;
    dz_S2       = iz_S2 - jz_S;
    dx_S3       = ix_S3 - jx_S;
    dy_S3       = iy_S3 - jy_S;
    dz_S3       = iz_S3 - jz_S;

    /* rsq = dx*dx+dy*dy+dz*dz */
    rsq_S0      = norm2(dx_S0, dy_S0, dz_S0);
    rsq_S1      = norm2(dx_S1, dy_S1, dz_S1);
    rsq_S2      = norm2(dx_S2, dy_S2, dz_S2);
    rsq_S3      = norm2(dx_S3, dy_S3, dz_S3);

    /* Do the cut-off check */
    wco_S0      = (rsq_S0 < rc2_S);
    wco_S1      = (rsq_S1 < rc2_S);
    wco_S2      = (rsq_S2 < rc2_S);
    wco_S3      = (rsq_S3 < rc2_S);

#ifdef CHECK_EXCLS
#ifdef EXCL_FORCES
    /* Only remove the (sub-)diagonal to avoid double counting */
#if UNROLLJ == UNROLLI
    if (cj == ci_sh)
    {
        wco_S0  = wco_S0 && diagonal_mask_S0;
        wco_S1  = wco_S1 && diagonal_mask_S1;
        wco_S2  = wco_S2 && diagonal_mask_S2;
        wco_S3  = wco_S3 && diagonal_mask_S3;
    }
#else
#if UNROLLJ < UNROLLI
    if (cj == ci_sh*2)
    {
        wco_S0  = wco_S0 && diagonal_mask0_S0;
        wco_S1  = wco_S1 && diagonal_mask0_S1;
        wco_S2  = wco_S2 && diagonal_mask0_S2;
        wco_S3  = wco_S3 && diagonal_mask0_S3;
    }
    if (cj == ci_sh*2 + 1)
    {
        wco_S0  = wco_S0 && diagonal_mask1_S0;
        wco_S1  = wco_S1 && diagonal_mask1_S1;
        wco_S2  = wco_S2 && diagonal_mask1_S2;
        wco_S3  = wco_S3 && diagonal_mask1_S3;
    }
#else
    if (cj*2 == ci_sh)
    {
        wco_S0  = wco_S0 && diagonal_mask0_S0;
        wco_S1  = wco_S1 && diagonal_mask0_S1;
        wco_S2  = wco_S2 && diagonal_mask0_S2;
        wco_S3  = wco_S3 && diagonal_mask0_S3;
    }
    else if (cj*2 + 1 == ci_sh)
    {
        wco_S0  = wco_S0 && diagonal_mask1_S0;
        wco_S1  = wco_S1 && diagonal_mask1_S1;
        wco_S2  = wco_S2 && diagonal_mask1_S2;
        wco_S3  = wco_S3 && diagonal_mask1_S3;
    }
#endif
#endif
#else /* EXCL_FORCES */
      /* No exclusion forces: remove all excluded atom pairs from the list */
    wco_S0      = wco_S0 && interact_S0;
    wco_S1      = wco_S1 && interact_S1;
    wco_S2      = wco_S2 && interact_S2;
    wco_S3      = wco_S3 && interact_S3;
#endif
#endif

#ifdef COUNT_PAIRS
    {
        int  i, j;
        GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH)  tmp[2*GMX_SIMD_REAL_WIDTH];

        for (i = 0; i < UNROLLI; i++)
        {
            store(tmp, rc2_S - (i == 0 ? rsq_S0 : (i == 1 ? rsq_S1 : (i == 2 ? rsq_S2 : rsq_S3))));
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

    // Ensure the distances do not fall below the limit where r^-12 overflows.
    // This should never happen for normal interactions.
    rsq_S0      = max(rsq_S0, minRsq_S);
    rsq_S1      = max(rsq_S1, minRsq_S);
    rsq_S2      = max(rsq_S2, minRsq_S);
    rsq_S3      = max(rsq_S3, minRsq_S);

    /* Calculate 1/r */
#if !GMX_DOUBLE
    rinv_S0     = invsqrt(rsq_S0);
    rinv_S1     = invsqrt(rsq_S1);
    rinv_S2     = invsqrt(rsq_S2);
    rinv_S3     = invsqrt(rsq_S3);
#else
    invsqrtPair(rsq_S0, rsq_S1, &rinv_S0, &rinv_S1);
    invsqrtPair(rsq_S2, rsq_S3, &rinv_S2, &rinv_S3);
#endif

#ifdef CALC_COULOMB
    /* Load parameters for j atom */
    jq_S        = load<SimdReal>(q+aj);
    qq_S0       = iq_S0 * jq_S;
    qq_S1       = iq_S1 * jq_S;
    qq_S2       = iq_S2 * jq_S;
    qq_S3       = iq_S3 * jq_S;
#endif

#ifdef CALC_LJ
#if !defined LJ_COMB_GEOM && !defined LJ_COMB_LB && !defined FIX_LJ_C
    SimdReal c6_S0, c6_S1, c12_S0, c12_S1;
    gatherLoadTranspose<c_simdBestPairAlignment>(nbfp0, type+aj, &c6_S0, &c12_S0);
    gatherLoadTranspose<c_simdBestPairAlignment>(nbfp1, type+aj, &c6_S1, &c12_S1);
#ifndef HALF_LJ
    SimdReal c6_S2, c6_S3, c12_S2, c12_S3;
    gatherLoadTranspose<c_simdBestPairAlignment>(nbfp2, type+aj, &c6_S2, &c12_S2);
    gatherLoadTranspose<c_simdBestPairAlignment>(nbfp3, type+aj, &c6_S3, &c12_S3);
#endif
#endif /* not defined any LJ rule */

#ifdef LJ_COMB_GEOM
    c6s_j_S     = load<SimdReal>(ljc+aj2+0);
    c12s_j_S    = load<SimdReal>(ljc+aj2+STRIDE);
    SimdReal c6_S0  = c6s_S0 * c6s_j_S;
    SimdReal c6_S1  = c6s_S1 * c6s_j_S;
#ifndef HALF_LJ
    SimdReal c6_S2  = c6s_S2 * c6s_j_S;
    SimdReal c6_S3  = c6s_S3 * c6s_j_S;
#endif
    SimdReal c12_S0 = c12s_S0 * c12s_j_S;
    SimdReal c12_S1 = c12s_S1 * c12s_j_S;
#ifndef HALF_LJ
    SimdReal c12_S2 = c12s_S2 * c12s_j_S;
    SimdReal c12_S3 = c12s_S3 * c12s_j_S;
#endif
#endif /* LJ_COMB_GEOM */

#ifdef LJ_COMB_LB
    hsig_j_S    = load<SimdReal>(ljc+aj2+0);
    seps_j_S    = load<SimdReal>(ljc+aj2+STRIDE);

    sig_S0      = hsig_i_S0 + hsig_j_S;
    sig_S1      = hsig_i_S1 + hsig_j_S;
    eps_S0      = seps_i_S0 * seps_j_S;
    eps_S1      = seps_i_S1 * seps_j_S;
#ifndef HALF_LJ
    sig_S2      = hsig_i_S2 + hsig_j_S;
    sig_S3      = hsig_i_S3 + hsig_j_S;
    eps_S2      = seps_i_S2 * seps_j_S;
    eps_S3      = seps_i_S3 * seps_j_S;
#endif
#endif /* LJ_COMB_LB */

#endif /* CALC_LJ */

    /* Set rinv to zero for r beyond the cut-off */
    rinv_S0     = selectByMask(rinv_S0, wco_S0);
    rinv_S1     = selectByMask(rinv_S1, wco_S1);
    rinv_S2     = selectByMask(rinv_S2, wco_S2);
    rinv_S3     = selectByMask(rinv_S3, wco_S3);

    rinvsq_S0   = rinv_S0 * rinv_S0;
    rinvsq_S1   = rinv_S1 * rinv_S1;
    rinvsq_S2   = rinv_S2 * rinv_S2;
    rinvsq_S3   = rinv_S3 * rinv_S3;

#ifdef CALC_COULOMB
    /* Note that here we calculate force*r, not the usual force/r.
     * This allows avoiding masking the reaction-field contribution,
     * as frcoul is later multiplied by rinvsq which has been
     * masked with the cut-off check.
     */

#ifdef EXCL_FORCES
    /* Only add 1/r for non-excluded atom pairs */
    rinv_ex_S0  = selectByMask(rinv_S0, interact_S0);
    rinv_ex_S1  = selectByMask(rinv_S1, interact_S1);
    rinv_ex_S2  = selectByMask(rinv_S2, interact_S2);
    rinv_ex_S3  = selectByMask(rinv_S3, interact_S3);
#else
    /* No exclusion forces, we always need 1/r */
#define     rinv_ex_S0    rinv_S0
#define     rinv_ex_S1    rinv_S1
#define     rinv_ex_S2    rinv_S2
#define     rinv_ex_S3    rinv_S3
#endif

#ifdef CALC_COUL_RF
    /* Electrostatic interactions */
    frcoul_S0   = qq_S0 * fma(rsq_S0, mrc_3_S, rinv_ex_S0);
    frcoul_S1   = qq_S1 * fma(rsq_S1, mrc_3_S, rinv_ex_S1);
    frcoul_S2   = qq_S2 * fma(rsq_S2, mrc_3_S, rinv_ex_S2);
    frcoul_S3   = qq_S3 * fma(rsq_S3, mrc_3_S, rinv_ex_S3);

#ifdef CALC_ENERGIES
    vcoul_S0    = qq_S0 * (rinv_ex_S0 + fma(rsq_S0, hrc_3_S, moh_rc_S));
    vcoul_S1    = qq_S1 * (rinv_ex_S1 + fma(rsq_S1, hrc_3_S, moh_rc_S));
    vcoul_S2    = qq_S2 * (rinv_ex_S2 + fma(rsq_S2, hrc_3_S, moh_rc_S));
    vcoul_S3    = qq_S3 * (rinv_ex_S3 + fma(rsq_S3, hrc_3_S, moh_rc_S));
#endif
#endif

#ifdef CALC_COUL_EWALD
    /* We need to mask (or limit) rsq for the cut-off,
     * as large distances can cause an overflow in gmx_pmecorrF/V.
     */
    brsq_S0     = beta2_S * selectByMask(rsq_S0, wco_S0);
    brsq_S1     = beta2_S * selectByMask(rsq_S1, wco_S1);
    brsq_S2     = beta2_S * selectByMask(rsq_S2, wco_S2);
    brsq_S3     = beta2_S * selectByMask(rsq_S3, wco_S3);
    ewcorr_S0   = beta_S * pmeForceCorrection(brsq_S0);
    ewcorr_S1   = beta_S * pmeForceCorrection(brsq_S1);
    ewcorr_S2   = beta_S * pmeForceCorrection(brsq_S2);
    ewcorr_S3   = beta_S * pmeForceCorrection(brsq_S3);
    frcoul_S0   = qq_S0 * fma(ewcorr_S0, brsq_S0, rinv_ex_S0);
    frcoul_S1   = qq_S1 * fma(ewcorr_S1, brsq_S1, rinv_ex_S1);
    frcoul_S2   = qq_S2 * fma(ewcorr_S2, brsq_S2, rinv_ex_S2);
    frcoul_S3   = qq_S3 * fma(ewcorr_S3, brsq_S3, rinv_ex_S3);

#ifdef CALC_ENERGIES
    vc_sub_S0   = beta_S * pmePotentialCorrection(brsq_S0);
    vc_sub_S1   = beta_S * pmePotentialCorrection(brsq_S1);
    vc_sub_S2   = beta_S * pmePotentialCorrection(brsq_S2);
    vc_sub_S3   = beta_S * pmePotentialCorrection(brsq_S3);
#endif

#endif /* CALC_COUL_EWALD */

#ifdef CALC_COUL_TAB
    /* Electrostatic interactions */
    r_S0        = rsq_S0 * rinv_S0;
    r_S1        = rsq_S1 * rinv_S1;
    r_S2        = rsq_S2 * rinv_S2;
    r_S3        = rsq_S3 * rinv_S3;
    /* Convert r to scaled table units */
    rs_S0       = r_S0 * invtsp_S;
    rs_S1       = r_S1 * invtsp_S;
    rs_S2       = r_S2 * invtsp_S;
    rs_S3       = r_S3 * invtsp_S;
    /* Truncate scaled r to an int */
    ti_S0       = cvttR2I(rs_S0);
    ti_S1       = cvttR2I(rs_S1);
    ti_S2       = cvttR2I(rs_S2);
    ti_S3       = cvttR2I(rs_S3);

    rf_S0       = trunc(rs_S0);
    rf_S1       = trunc(rs_S1);
    rf_S2       = trunc(rs_S2);
    rf_S3       = trunc(rs_S3);

    frac_S0     = rs_S0 - rf_S0;
    frac_S1     = rs_S1 - rf_S1;
    frac_S2     = rs_S2 - rf_S2;
    frac_S3     = rs_S3 - rf_S3;

    /* Load and interpolate table forces and possibly energies.
     * Force and energy can be combined in one table, stride 4: FDV0
     * or in two separate tables with stride 1: F and V
     * Currently single precision uses FDV0, double F and V.
     */
#ifndef CALC_ENERGIES
#ifdef TAB_FDV0
    gatherLoadBySimdIntTranspose<4>(tab_coul_F, ti_S0, &ctab0_S0, &ctab1_S0);
    gatherLoadBySimdIntTranspose<4>(tab_coul_F, ti_S1, &ctab0_S1, &ctab1_S1);
    gatherLoadBySimdIntTranspose<4>(tab_coul_F, ti_S2, &ctab0_S2, &ctab1_S2);
    gatherLoadBySimdIntTranspose<4>(tab_coul_F, ti_S3, &ctab0_S3, &ctab1_S3);
#else
    gatherLoadUBySimdIntTranspose<1>(tab_coul_F, ti_S0, &ctab0_S0, &ctab1_S0);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_F, ti_S1, &ctab0_S1, &ctab1_S1);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_F, ti_S2, &ctab0_S2, &ctab1_S2);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_F, ti_S3, &ctab0_S3, &ctab1_S3);
    ctab1_S0 = ctab1_S0 - ctab0_S0;
    ctab1_S1 = ctab1_S1 - ctab0_S1;
    ctab1_S2 = ctab1_S2 - ctab0_S2;
    ctab1_S3 = ctab1_S3 - ctab0_S3;
#endif
#else
#ifdef TAB_FDV0
    gatherLoadBySimdIntTranspose<4>(tab_coul_F, ti_S0, &ctab0_S0, &ctab1_S0, &ctabv_S0, &dum_S0);
    gatherLoadBySimdIntTranspose<4>(tab_coul_F, ti_S1, &ctab0_S1, &ctab1_S1, &ctabv_S1, &dum_S1);
    gatherLoadBySimdIntTranspose<4>(tab_coul_F, ti_S2, &ctab0_S2, &ctab1_S2, &ctabv_S2, &dum_S2);
    gatherLoadBySimdIntTranspose<4>(tab_coul_F, ti_S3, &ctab0_S3, &ctab1_S3, &ctabv_S3, &dum_S3);
#else
    gatherLoadUBySimdIntTranspose<1>(tab_coul_F, ti_S0, &ctab0_S0, &ctab1_S0);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_V, ti_S0, &ctabv_S0, &dum_S0);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_F, ti_S1, &ctab0_S1, &ctab1_S1);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_V, ti_S1, &ctabv_S1, &dum_S1);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_F, ti_S2, &ctab0_S2, &ctab1_S2);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_V, ti_S2, &ctabv_S2, &dum_S2);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_F, ti_S3, &ctab0_S3, &ctab1_S3);
    gatherLoadUBySimdIntTranspose<1>(tab_coul_V, ti_S3, &ctabv_S3, &dum_S3);
    ctab1_S0 = ctab1_S0 - ctab0_S0;
    ctab1_S1 = ctab1_S1 - ctab0_S1;
    ctab1_S2 = ctab1_S2 - ctab0_S2;
    ctab1_S3 = ctab1_S3 - ctab0_S3;
#endif
#endif
    fsub_S0     = fma(frac_S0, ctab1_S0, ctab0_S0);
    fsub_S1     = fma(frac_S1, ctab1_S1, ctab0_S1);
    fsub_S2     = fma(frac_S2, ctab1_S2, ctab0_S2);
    fsub_S3     = fma(frac_S3, ctab1_S3, ctab0_S3);
    frcoul_S0   = qq_S0 * fnma(fsub_S0, r_S0, rinv_ex_S0);
    frcoul_S1   = qq_S1 * fnma(fsub_S1, r_S1, rinv_ex_S1);
    frcoul_S2   = qq_S2 * fnma(fsub_S2, r_S2, rinv_ex_S2);
    frcoul_S3   = qq_S3 * fnma(fsub_S3, r_S3, rinv_ex_S3);

#ifdef CALC_ENERGIES
    vc_sub_S0   = fma((mhalfsp_S * frac_S0), (ctab0_S0 + fsub_S0), ctabv_S0);
    vc_sub_S1   = fma((mhalfsp_S * frac_S1), (ctab0_S1 + fsub_S1), ctabv_S1);
    vc_sub_S2   = fma((mhalfsp_S * frac_S2), (ctab0_S2 + fsub_S2), ctabv_S2);
    vc_sub_S3   = fma((mhalfsp_S * frac_S3), (ctab0_S3 + fsub_S3), ctabv_S3);
#endif
#endif /* CALC_COUL_TAB */

#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
#ifndef NO_SHIFT_EWALD
    /* Add Ewald potential shift to vc_sub for convenience */
#ifdef CHECK_EXCLS
    vc_sub_S0   = vc_sub_S0 + selectByMask(sh_ewald_S, interact_S0);
    vc_sub_S1   = vc_sub_S1 + selectByMask(sh_ewald_S, interact_S1);
    vc_sub_S2   = vc_sub_S2 + selectByMask(sh_ewald_S, interact_S2);
    vc_sub_S3   = vc_sub_S3 + selectByMask(sh_ewald_S, interact_S3);
#else
    vc_sub_S0   = vc_sub_S0 + sh_ewald_S;
    vc_sub_S1   = vc_sub_S1 + sh_ewald_S;
    vc_sub_S2   = vc_sub_S2 + sh_ewald_S;
    vc_sub_S3   = vc_sub_S3 + sh_ewald_S;
#endif
#endif

    vcoul_S0    = qq_S0 * (rinv_ex_S0 - vc_sub_S0);
    vcoul_S1    = qq_S1 * (rinv_ex_S1 - vc_sub_S1);
    vcoul_S2    = qq_S2 * (rinv_ex_S2 - vc_sub_S2);
    vcoul_S3    = qq_S3 * (rinv_ex_S3 - vc_sub_S3);

#endif

#ifdef CALC_ENERGIES
    /* Mask energy for cut-off and diagonal */
    vcoul_S0    = selectByMask(vcoul_S0, wco_S0);
    vcoul_S1    = selectByMask(vcoul_S1, wco_S1);
    vcoul_S2    = selectByMask(vcoul_S2, wco_S2);
    vcoul_S3    = selectByMask(vcoul_S3, wco_S3);
#endif

#endif /* CALC_COULOMB */

#ifdef CALC_LJ
    /* Lennard-Jones interaction */

#ifdef VDW_CUTOFF_CHECK
    wco_vdw_S0  = (rsq_S0 < rcvdw2_S);
    wco_vdw_S1  = (rsq_S1 < rcvdw2_S);
#ifndef HALF_LJ
    wco_vdw_S2  = (rsq_S2 < rcvdw2_S);
    wco_vdw_S3  = (rsq_S3 < rcvdw2_S);
#endif
#else
    /* Same cut-off for Coulomb and VdW, reuse the registers */
#define     wco_vdw_S0    wco_S0
#define     wco_vdw_S1    wco_S1
#define     wco_vdw_S2    wco_S2
#define     wco_vdw_S3    wco_S3
#endif

#ifndef LJ_COMB_LB
    rinvsix_S0  = rinvsq_S0 * rinvsq_S0 * rinvsq_S0;
    rinvsix_S1  = rinvsq_S1 * rinvsq_S1 * rinvsq_S1;
#ifdef EXCL_FORCES
    rinvsix_S0  = selectByMask(rinvsix_S0, interact_S0);
    rinvsix_S1  = selectByMask(rinvsix_S1, interact_S1);
#endif
#ifndef HALF_LJ
    rinvsix_S2  = rinvsq_S2 * rinvsq_S2 * rinvsq_S2;
    rinvsix_S3  = rinvsq_S3 * rinvsq_S3 * rinvsq_S3;
#ifdef EXCL_FORCES
    rinvsix_S2  = selectByMask(rinvsix_S2, interact_S2);
    rinvsix_S3  = selectByMask(rinvsix_S3, interact_S3);
#endif
#endif

#if defined LJ_CUT || defined LJ_POT_SWITCH
    /* We have plain LJ or LJ-PME with simple C6/6 C12/12 coefficients */
    FrLJ6_S0    = c6_S0 * rinvsix_S0;
    FrLJ6_S1    = c6_S1 * rinvsix_S1;
#ifndef HALF_LJ
    FrLJ6_S2    = c6_S2 * rinvsix_S2;
    FrLJ6_S3    = c6_S3 * rinvsix_S3;
#endif
    FrLJ12_S0   = c12_S0 * rinvsix_S0 * rinvsix_S0;
    FrLJ12_S1   = c12_S1 * rinvsix_S1 * rinvsix_S1;
#ifndef HALF_LJ
    FrLJ12_S2   = c12_S2 * rinvsix_S2 * rinvsix_S2;
    FrLJ12_S3   = c12_S3 * rinvsix_S3 * rinvsix_S3;
#endif
#endif

#if defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
    /* We switch the LJ force */
    r_S0        = rsq_S0 * rinv_S0;
    rsw_S0      = max(r_S0 - rswitch_S, zero_S);
    rsw2_S0     = rsw_S0 * rsw_S0;
    r_S1        = rsq_S1 * rinv_S1;
    rsw_S1      = max(r_S1 - rswitch_S, zero_S);
    rsw2_S1     = rsw_S1 * rsw_S1;
#ifndef HALF_LJ
    r_S2        = rsq_S2 * rinv_S2;
    rsw_S2      = max(r_S2 - rswitch_S, zero_S);
    rsw2_S2     = rsw_S2 * rsw_S2;
    r_S3        = rsq_S3 * rinv_S3;
    rsw_S3      = max(r_S3 - rswitch_S, zero_S);
    rsw2_S3     = rsw_S3 * rsw_S3;
#endif
#endif

#ifdef LJ_FORCE_SWITCH

#define gmx_add_fr_switch(fr, rsw, rsw2_r, c2, c3) fma(fma(c3, rsw, c2), rsw2_r, (fr))
    SimdReal rsw2_r_S0 = rsw2_S0 * r_S0;
    SimdReal rsw2_r_S1 = rsw2_S1 * r_S1;
    FrLJ6_S0    = c6_S0 * gmx_add_fr_switch(rinvsix_S0, rsw_S0, rsw2_r_S0, p6_fc2_S, p6_fc3_S);
    FrLJ6_S1    = c6_S1 * gmx_add_fr_switch(rinvsix_S1, rsw_S1, rsw2_r_S1, p6_fc2_S, p6_fc3_S);
#ifndef HALF_LJ
    SimdReal rsw2_r_S2 = rsw2_S2 * r_S2;
    SimdReal rsw2_r_S3 = rsw2_S3 * r_S3;
    FrLJ6_S2    = c6_S2 * gmx_add_fr_switch(rinvsix_S2, rsw_S2, rsw2_r_S2, p6_fc2_S, p6_fc3_S);
    FrLJ6_S3    = c6_S3 * gmx_add_fr_switch(rinvsix_S3, rsw_S3, rsw2_r_S3, p6_fc2_S, p6_fc3_S);
#endif
    FrLJ12_S0   = c12_S0 * gmx_add_fr_switch((rinvsix_S0 * rinvsix_S0), rsw_S0, rsw2_r_S0, p12_fc2_S, p12_fc3_S);
    FrLJ12_S1   = c12_S1 * gmx_add_fr_switch((rinvsix_S1 * rinvsix_S1), rsw_S1, rsw2_r_S1, p12_fc2_S, p12_fc3_S);
#ifndef HALF_LJ
    FrLJ12_S2   = c12_S2 * gmx_add_fr_switch((rinvsix_S2 * rinvsix_S2), rsw_S2, rsw2_r_S2, p12_fc2_S, p12_fc3_S);
    FrLJ12_S3   = c12_S3 * gmx_add_fr_switch((rinvsix_S3 * rinvsix_S3), rsw_S3, rsw2_r_S3, p12_fc2_S, p12_fc3_S);
#endif
#undef gmx_add_fr_switch
#endif /* LJ_FORCE_SWITCH */

#endif /* not LJ_COMB_LB */

#ifdef LJ_COMB_LB
    sir_S0      = sig_S0 * rinv_S0;
    sir_S1      = sig_S1 * rinv_S1;
#ifndef HALF_LJ
    sir_S2      = sig_S2 * rinv_S2;
    sir_S3      = sig_S3 * rinv_S3;
#endif
    sir2_S0     = sir_S0 * sir_S0;
    sir2_S1     = sir_S1 * sir_S1;
#ifndef HALF_LJ
    sir2_S2     = sir_S2 * sir_S2;
    sir2_S3     = sir_S3 * sir_S3;
#endif
    sir6_S0     = sir2_S0 * sir2_S0 * sir2_S0;
    sir6_S1     = sir2_S1 * sir2_S1 * sir2_S1;
#ifdef EXCL_FORCES
    sir6_S0     = selectByMask(sir6_S0, interact_S0);
    sir6_S1     = selectByMask(sir6_S1, interact_S1);
#endif
#ifndef HALF_LJ
    sir6_S2     = sir2_S2 * sir2_S2 * sir2_S2;
    sir6_S3     = sir2_S3 * sir2_S3 * sir2_S3;
#ifdef EXCL_FORCES
    sir6_S2     = selectByMask(sir6_S2, interact_S2);
    sir6_S3     = selectByMask(sir6_S3, interact_S3);
#endif
#endif
#ifdef VDW_CUTOFF_CHECK
    sir6_S0     = selectByMask(sir6_S0, wco_vdw_S0);
    sir6_S1     = selectByMask(sir6_S1, wco_vdw_S1);
#ifndef HALF_LJ
    sir6_S2     = selectByMask(sir6_S2, wco_vdw_S2);
    sir6_S3     = selectByMask(sir6_S3, wco_vdw_S3);
#endif
#endif
    FrLJ6_S0    = eps_S0 * sir6_S0;
    FrLJ6_S1    = eps_S1 * sir6_S1;
#ifndef HALF_LJ
    FrLJ6_S2    = eps_S2 * sir6_S2;
    FrLJ6_S3    = eps_S3 * sir6_S3;
#endif
    FrLJ12_S0   = FrLJ6_S0 * sir6_S0;
    FrLJ12_S1   = FrLJ6_S1 * sir6_S1;
#ifndef HALF_LJ
    FrLJ12_S2   = FrLJ6_S2 * sir6_S2;
    FrLJ12_S3   = FrLJ6_S3 * sir6_S3;
#endif
#if defined CALC_ENERGIES
    /* We need C6 and C12 to calculate the LJ potential shift */
    sig2_S0     = sig_S0 * sig_S0;
    sig2_S1     = sig_S1 * sig_S1;
#ifndef HALF_LJ
    sig2_S2     = sig_S2 * sig_S2;
    sig2_S3     = sig_S3 * sig_S3;
#endif
    sig6_S0     = sig2_S0 * sig2_S0 * sig2_S0;
    sig6_S1     = sig2_S1 * sig2_S1 * sig2_S1;
#ifndef HALF_LJ
    sig6_S2     = sig2_S2 * sig2_S2 * sig2_S2;
    sig6_S3     = sig2_S3 * sig2_S3 * sig2_S3;
#endif
    SimdReal c6_S0  = eps_S0 * sig6_S0;
    SimdReal c6_S1  = eps_S1 * sig6_S1;
#ifndef HALF_LJ
    SimdReal c6_S2  = eps_S2 * sig6_S2;
    SimdReal c6_S3  = eps_S3 * sig6_S3;
#endif
    SimdReal c12_S0 = c6_S0 * sig6_S0;
    SimdReal c12_S1 = c6_S1 * sig6_S1;
#ifndef HALF_LJ
    SimdReal c12_S2 = c6_S2 * sig6_S2;
    SimdReal c12_S3 = c6_S3 * sig6_S3;
#endif
#endif
#endif /* LJ_COMB_LB */

    /* Determine the total scalar LJ force*r */
    frLJ_S0     = FrLJ12_S0 - FrLJ6_S0;
    frLJ_S1     = FrLJ12_S1 - FrLJ6_S1;
#ifndef HALF_LJ
    frLJ_S2     = FrLJ12_S2 - FrLJ6_S2;
    frLJ_S3     = FrLJ12_S3 - FrLJ6_S3;
#endif

#if (defined LJ_CUT || defined LJ_FORCE_SWITCH) && defined CALC_ENERGIES

#ifdef LJ_CUT
    /* Calculate the LJ energies, with constant potential shift */
    SimdReal VLJ6_S0  = sixth_S * fma(c6_S0, p6_cpot_S, FrLJ6_S0);
    SimdReal VLJ6_S1  = sixth_S * fma(c6_S1, p6_cpot_S, FrLJ6_S1);
#ifndef HALF_LJ
    SimdReal VLJ6_S2  = sixth_S * fma(c6_S2, p6_cpot_S, FrLJ6_S2);
    SimdReal VLJ6_S3  = sixth_S * fma(c6_S3, p6_cpot_S, FrLJ6_S3);
#endif
    SimdReal VLJ12_S0 = twelveth_S * fma(c12_S0, p12_cpot_S, FrLJ12_S0);
    SimdReal VLJ12_S1 = twelveth_S * fma(c12_S1, p12_cpot_S, FrLJ12_S1);
#ifndef HALF_LJ
    SimdReal VLJ12_S2 = twelveth_S * fma(c12_S2, p12_cpot_S, FrLJ12_S2);
    SimdReal VLJ12_S3 = twelveth_S * fma(c12_S3, p12_cpot_S, FrLJ12_S3);
#endif
#endif /* LJ_CUT */

#ifdef LJ_FORCE_SWITCH
#define v_fswitch_r(rsw, rsw2, c0, c3, c4) fma(fma(c4, rsw, c3), (rsw2 * rsw), c0)

    SimdReal VLJ6_S0  = c6_S0 * fma(sixth_S, rinvsix_S0, v_fswitch_r(rsw_S0, rsw2_S0, p6_6cpot_S, p6_vc3_S, p6_vc4_S));
    SimdReal VLJ6_S1  = c6_S1 * fma(sixth_S, rinvsix_S1, v_fswitch_r(rsw_S1, rsw2_S1, p6_6cpot_S, p6_vc3_S, p6_vc4_S));
#ifndef HALF_LJ
    SimdReal VLJ6_S2  = c6_S2 * fma(sixth_S, rinvsix_S2, v_fswitch_r(rsw_S2, rsw2_S2, p6_6cpot_S, p6_vc3_S, p6_vc4_S));
    SimdReal VLJ6_S3  = c6_S3 * fma(sixth_S, rinvsix_S3, v_fswitch_r(rsw_S3, rsw2_S3, p6_6cpot_S, p6_vc3_S, p6_vc4_S));
#endif
    SimdReal VLJ12_S0 = c12_S0 * fma(twelveth_S, rinvsix_S0 * rinvsix_S0, v_fswitch_r(rsw_S0, rsw2_S0, p12_12cpot_S, p12_vc3_S, p12_vc4_S));
    SimdReal VLJ12_S1 = c12_S1 * fma(twelveth_S, rinvsix_S1 * rinvsix_S1, v_fswitch_r(rsw_S1, rsw2_S1, p12_12cpot_S, p12_vc3_S, p12_vc4_S));
#ifndef HALF_LJ
    SimdReal VLJ12_S2 = c12_S2 * fma(twelveth_S, rinvsix_S2 * rinvsix_S2, v_fswitch_r(rsw_S2, rsw2_S2, p12_12cpot_S, p12_vc3_S, p12_vc4_S));
    SimdReal VLJ12_S3 = c12_S3 * fma(twelveth_S, rinvsix_S3 * rinvsix_S3, v_fswitch_r(rsw_S3, rsw2_S3, p12_12cpot_S, p12_vc3_S, p12_vc4_S));
#endif
#undef v_fswitch_r
#endif /* LJ_FORCE_SWITCH */

    /* Add up the repulsion and dispersion */
    SimdReal VLJ_S0 = VLJ12_S0 - VLJ6_S0;
    SimdReal VLJ_S1 = VLJ12_S1 - VLJ6_S1;
#ifndef HALF_LJ
    SimdReal VLJ_S2 = VLJ12_S2 - VLJ6_S2;
    SimdReal VLJ_S3 = VLJ12_S3 - VLJ6_S3;
#endif

#endif /* (LJ_CUT || LJ_FORCE_SWITCH) && CALC_ENERGIES */

#ifdef LJ_POT_SWITCH
    /* We always need the potential, since it is needed for the force */
    SimdReal VLJ_S0 = fnma(sixth_S, FrLJ6_S0, twelveth_S * FrLJ12_S0);
    SimdReal VLJ_S1 = fnma(sixth_S, FrLJ6_S1, twelveth_S * FrLJ12_S1);
#ifndef HALF_LJ
    SimdReal VLJ_S2 = fnma(sixth_S, FrLJ6_S2, twelveth_S * FrLJ12_S2);
    SimdReal VLJ_S3 = fnma(sixth_S, FrLJ6_S3, twelveth_S * FrLJ12_S3);
#endif

    {
        SimdReal sw_S0, dsw_S0;
        SimdReal sw_S1, dsw_S1;
#ifndef HALF_LJ
        SimdReal sw_S2, dsw_S2;
        SimdReal sw_S3, dsw_S3;
#endif

#define switch_r(rsw, rsw2, c3, c4, c5) fma(fma(fma(c5, rsw, c4), rsw, c3), (rsw2 * rsw), one_S)
#define dswitch_r(rsw, rsw2, c2, c3, c4) (fma(fma(c4, rsw, c3), rsw, c2) * rsw2)

        sw_S0  = switch_r(rsw_S0, rsw2_S0, swV3_S, swV4_S, swV5_S);
        dsw_S0 = dswitch_r(rsw_S0, rsw2_S0, swF2_S, swF3_S, swF4_S);
        sw_S1  = switch_r(rsw_S1, rsw2_S1, swV3_S, swV4_S, swV5_S);
        dsw_S1 = dswitch_r(rsw_S1, rsw2_S1, swF2_S, swF3_S, swF4_S);
#ifndef HALF_LJ
        sw_S2  = switch_r(rsw_S2, rsw2_S2, swV3_S, swV4_S, swV5_S);
        dsw_S2 = dswitch_r(rsw_S2, rsw2_S2, swF2_S, swF3_S, swF4_S);
        sw_S3  = switch_r(rsw_S3, rsw2_S3, swV3_S, swV4_S, swV5_S);
        dsw_S3 = dswitch_r(rsw_S3, rsw2_S3, swF2_S, swF3_S, swF4_S);
#endif
        frLJ_S0 = fnma(dsw_S0 * VLJ_S0, r_S0, sw_S0 * frLJ_S0);
        frLJ_S1 = fnma(dsw_S1 * VLJ_S1, r_S1, sw_S1 * frLJ_S1);
#ifndef HALF_LJ
        frLJ_S2 = fnma(dsw_S2 * VLJ_S2, r_S2, sw_S2 * frLJ_S2);
        frLJ_S3 = fnma(dsw_S3 * VLJ_S3, r_S3, sw_S3 * frLJ_S3);
#endif
#ifdef CALC_ENERGIES
        VLJ_S0  = sw_S0 * VLJ_S0;
        VLJ_S1  = sw_S1 * VLJ_S1;
#ifndef HALF_LJ
        VLJ_S2  = sw_S2 * VLJ_S2;
        VLJ_S3  = sw_S3 * VLJ_S3;
#endif
#endif

#undef switch_r
#undef dswitch_r
    }
#endif /* LJ_POT_SWITCH */

#if defined CALC_ENERGIES && defined CHECK_EXCLS
    /* The potential shift should be removed for excluded pairs */
    VLJ_S0      = selectByMask(VLJ_S0, interact_S0);
    VLJ_S1      = selectByMask(VLJ_S1, interact_S1);
#ifndef HALF_LJ
    VLJ_S2      = selectByMask(VLJ_S2, interact_S2);
    VLJ_S3      = selectByMask(VLJ_S3, interact_S3);
#endif
#endif

#ifdef LJ_EWALD_GEOM
    {
        SimdReal c6s_j_S;
        SimdReal c6grid_S0, rinvsix_nm_S0, cr2_S0, expmcr2_S0, poly_S0;
        SimdReal c6grid_S1, rinvsix_nm_S1, cr2_S1, expmcr2_S1, poly_S1;
#ifndef HALF_LJ
        SimdReal c6grid_S2, rinvsix_nm_S2, cr2_S2, expmcr2_S2, poly_S2;
        SimdReal c6grid_S3, rinvsix_nm_S3, cr2_S3, expmcr2_S3, poly_S3;
#endif
#ifdef CALC_ENERGIES
        SimdReal sh_mask_S0;
        SimdReal sh_mask_S1;
#ifndef HALF_LJ
        SimdReal sh_mask_S2;
        SimdReal sh_mask_S3;
#endif
#endif

        /* Determine C6 for the grid using the geometric combination rule */
        c6s_j_S         = load<SimdReal>(ljc+aj2+0);
        c6grid_S0       = c6s_S0 * c6s_j_S;
        c6grid_S1       = c6s_S1 * c6s_j_S;
#ifndef HALF_LJ
        c6grid_S2       = c6s_S2 * c6s_j_S;
        c6grid_S3       = c6s_S3 * c6s_j_S;
#endif

#ifdef CHECK_EXCLS
        /* Recalculate rinvsix without exclusion mask (compiler might optimize) */
        rinvsix_nm_S0 = rinvsq_S0 * rinvsq_S0 * rinvsq_S0;
        rinvsix_nm_S1 = rinvsq_S1 * rinvsq_S1 * rinvsq_S1;
#ifndef HALF_LJ
        rinvsix_nm_S2 = rinvsq_S2 * rinvsq_S2 * rinvsq_S2;
        rinvsix_nm_S3 = rinvsq_S3 * rinvsq_S3 * rinvsq_S3;
#endif
#else
        /* We didn't use a mask, so we can copy */
        rinvsix_nm_S0 = rinvsix_S0;
        rinvsix_nm_S1 = rinvsix_S1;
#ifndef HALF_LJ
        rinvsix_nm_S2 = rinvsix_S2;
        rinvsix_nm_S3 = rinvsix_S3;
#endif
#endif

        /* Mask for the cut-off to avoid overflow of cr2^2 */
        cr2_S0        = lje_c2_S * selectByMask(rsq_S0, wco_vdw_S0);
        cr2_S1        = lje_c2_S * selectByMask(rsq_S1, wco_vdw_S1);
#ifndef HALF_LJ
        cr2_S2        = lje_c2_S * selectByMask(rsq_S2, wco_vdw_S2);
        cr2_S3        = lje_c2_S * selectByMask(rsq_S3, wco_vdw_S3);
#endif
        // Unsafe version of our exp() should be fine, since these arguments should never
        // be smaller than -127 for any reasonable choice of cutoff or ewald coefficients.
        expmcr2_S0    = exp<MathOptimization::Unsafe>( -cr2_S0 );
        expmcr2_S1    = exp<MathOptimization::Unsafe>( -cr2_S1 );
#ifndef HALF_LJ
        expmcr2_S2    = exp<MathOptimization::Unsafe>( -cr2_S2 );
        expmcr2_S3    = exp<MathOptimization::Unsafe>( -cr2_S3 );
#endif

        /* 1 + cr2 + 1/2*cr2^2 */
        poly_S0       = fma(fma(half_S, cr2_S0, one_S), cr2_S0, one_S);
        poly_S1       = fma(fma(half_S, cr2_S1, one_S), cr2_S1, one_S);
#ifndef HALF_LJ
        poly_S2       = fma(fma(half_S, cr2_S2, one_S), cr2_S2, one_S);
        poly_S3       = fma(fma(half_S, cr2_S3, one_S), cr2_S3, one_S);
#endif

        /* We calculate LJ F*r = (6*C6)*(r^-6 - F_mesh/6), we use:
         * r^-6*cexp*(1 + cr2 + cr2^2/2 + cr2^3/6) = cexp*(r^-6*poly + c^6/6)
         */
        frLJ_S0       = fma(c6grid_S0, fnma(expmcr2_S0, fma(rinvsix_nm_S0, poly_S0, lje_c6_6_S), rinvsix_nm_S0), frLJ_S0);
        frLJ_S1       = fma(c6grid_S1, fnma(expmcr2_S1, fma(rinvsix_nm_S1, poly_S1, lje_c6_6_S), rinvsix_nm_S1), frLJ_S1);
#ifndef HALF_LJ
        frLJ_S2       = fma(c6grid_S2, fnma(expmcr2_S2, fma(rinvsix_nm_S2, poly_S2, lje_c6_6_S), rinvsix_nm_S2), frLJ_S2);
        frLJ_S3       = fma(c6grid_S3, fnma(expmcr2_S3, fma(rinvsix_nm_S3, poly_S3, lje_c6_6_S), rinvsix_nm_S3), frLJ_S3);
#endif

#ifdef CALC_ENERGIES
#ifdef CHECK_EXCLS
        sh_mask_S0    = selectByMask(lje_vc_S, interact_S0);
        sh_mask_S1    = selectByMask(lje_vc_S, interact_S1);
#ifndef HALF_LJ
        sh_mask_S2    = selectByMask(lje_vc_S, interact_S2);
        sh_mask_S3    = selectByMask(lje_vc_S, interact_S3);
#endif
#else
        sh_mask_S0    = lje_vc_S;
        sh_mask_S1    = lje_vc_S;
#ifndef HALF_LJ
        sh_mask_S2    = lje_vc_S;
        sh_mask_S3    = lje_vc_S;
#endif
#endif

        VLJ_S0        = fma(sixth_S * c6grid_S0, fma(rinvsix_nm_S0, fnma(expmcr2_S0, poly_S0, one_S), sh_mask_S0), VLJ_S0);
        VLJ_S1        = fma(sixth_S * c6grid_S1, fma(rinvsix_nm_S1, fnma(expmcr2_S1, poly_S1, one_S), sh_mask_S1), VLJ_S1);
#ifndef HALF_LJ
        VLJ_S2        = fma(sixth_S * c6grid_S2, fma(rinvsix_nm_S2, fnma(expmcr2_S2, poly_S2, one_S), sh_mask_S2), VLJ_S2);
        VLJ_S3        = fma(sixth_S * c6grid_S3, fma(rinvsix_nm_S3, fnma(expmcr2_S3, poly_S3, one_S), sh_mask_S3), VLJ_S3);
#endif
#endif /* CALC_ENERGIES */
    }
#endif /* LJ_EWALD_GEOM */

#if defined VDW_CUTOFF_CHECK
    /* frLJ is multiplied later by rinvsq, which is masked for the Coulomb
     * cut-off, but if the VdW cut-off is shorter, we need to mask with that.
     */
    frLJ_S0     = selectByMask(frLJ_S0, wco_vdw_S0);
    frLJ_S1     = selectByMask(frLJ_S1, wco_vdw_S1);
#ifndef HALF_LJ
    frLJ_S2     = selectByMask(frLJ_S2, wco_vdw_S2);
    frLJ_S3     = selectByMask(frLJ_S3, wco_vdw_S3);
#endif
#endif

#ifdef CALC_ENERGIES
    /* The potential shift should be removed for pairs beyond cut-off */
    VLJ_S0      = selectByMask(VLJ_S0, wco_vdw_S0);
    VLJ_S1      = selectByMask(VLJ_S1, wco_vdw_S1);
#ifndef HALF_LJ
    VLJ_S2      = selectByMask(VLJ_S2, wco_vdw_S2);
    VLJ_S3      = selectByMask(VLJ_S3, wco_vdw_S3);
#endif
#endif

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
    vctot_S      = vctot_S + vcoul_S0 + vcoul_S1 + vcoul_S2 + vcoul_S3;
#else
    add_ener_grp(vcoul_S0, vctp[0], egp_jj);
    add_ener_grp(vcoul_S1, vctp[1], egp_jj);
    add_ener_grp(vcoul_S2, vctp[2], egp_jj);
    add_ener_grp(vcoul_S3, vctp[3], egp_jj);
#endif
#endif

#ifdef CALC_LJ

#ifndef ENERGY_GROUPS
#ifndef HALF_LJ
    Vvdwtot_S   = Vvdwtot_S + VLJ_S0 + VLJ_S1 + VLJ_S2 + VLJ_S3;
#else
    Vvdwtot_S   = Vvdwtot_S + VLJ_S0 + VLJ_S1;
#endif
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
#ifdef CALC_COULOMB
    fscal_S0    = rinvsq_S0 * (frcoul_S0 + frLJ_S0);
#else
    fscal_S0    = rinvsq_S0 * frLJ_S0;
#endif
#ifdef CALC_COULOMB
    fscal_S1    = rinvsq_S1 * (frcoul_S1 + frLJ_S1);
#else
    fscal_S1    = rinvsq_S1 * frLJ_S1;
#endif
#else
    fscal_S0    = rinvsq_S0 * frcoul_S0;
    fscal_S1    = rinvsq_S1 * frcoul_S1;
#endif /* CALC_LJ */
#if defined CALC_LJ && !defined HALF_LJ
#ifdef CALC_COULOMB
    fscal_S2    = rinvsq_S2 * (frcoul_S2 + frLJ_S2);
    fscal_S3    = rinvsq_S3 * (frcoul_S3 + frLJ_S3);
#else
    fscal_S2    = rinvsq_S2 * frLJ_S2;
    fscal_S3    = rinvsq_S3 * frLJ_S3;
#endif
#else
    /* Atom 2 and 3 don't have LJ, so only add Coulomb forces */
    fscal_S2    = rinvsq_S2 * frcoul_S2;
    fscal_S3    = rinvsq_S3 * frcoul_S3;
#endif

    /* Calculate temporary vectorial force */
    tx_S0       = fscal_S0 * dx_S0;
    tx_S1       = fscal_S1 * dx_S1;
    tx_S2       = fscal_S2 * dx_S2;
    tx_S3       = fscal_S3 * dx_S3;
    ty_S0       = fscal_S0 * dy_S0;
    ty_S1       = fscal_S1 * dy_S1;
    ty_S2       = fscal_S2 * dy_S2;
    ty_S3       = fscal_S3 * dy_S3;
    tz_S0       = fscal_S0 * dz_S0;
    tz_S1       = fscal_S1 * dz_S1;
    tz_S2       = fscal_S2 * dz_S2;
    tz_S3       = fscal_S3 * dz_S3;

    /* Increment i atom force */
    fix_S0      = fix_S0 + tx_S0;
    fix_S1      = fix_S1 + tx_S1;
    fix_S2      = fix_S2 + tx_S2;
    fix_S3      = fix_S3 + tx_S3;
    fiy_S0      = fiy_S0 + ty_S0;
    fiy_S1      = fiy_S1 + ty_S1;
    fiy_S2      = fiy_S2 + ty_S2;
    fiy_S3      = fiy_S3 + ty_S3;
    fiz_S0      = fiz_S0 + tz_S0;
    fiz_S1      = fiz_S1 + tz_S1;
    fiz_S2      = fiz_S2 + tz_S2;
    fiz_S3      = fiz_S3 + tz_S3;

    /* Decrement j atom force */
    store(f+ajx, load<SimdReal>(f+ajx) - (tx_S0 + tx_S1 + tx_S2 + tx_S3));
    store(f+ajy, load<SimdReal>(f+ajy) - (ty_S0 + ty_S1 + ty_S2 + ty_S3));
    store(f+ajz, load<SimdReal>(f+ajz) - (tz_S0 + tz_S1 + tz_S2 + tz_S3));
}

#undef  rinv_ex_S0
#undef  rinv_ex_S1
#undef  rinv_ex_S2
#undef  rinv_ex_S3

#undef  wco_vdw_S0
#undef  wco_vdw_S1
#undef  wco_vdw_S2
#undef  wco_vdw_S3

#undef  EXCL_FORCES

#endif // !DOXYGEN
