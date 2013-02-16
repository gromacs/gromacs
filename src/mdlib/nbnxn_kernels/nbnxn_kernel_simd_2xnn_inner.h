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
 * modify it under the terms of the GNU LeSIMDr General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * LeSIMDr General Public License for more details.
 *
 * You should have received a copy of the GNU LeSIMDr General Public
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
 * This flavor of the kernel duplicates the data for N j-particles in
 * 2xN wide SIMD registers to do operate on 2 i-particles at once.
 * This leads to 4/2=2 sets of most instructions. Therefore we call
 * this kernel 2x(N+N) = 2xnn
 *
 * This 2xnn kernel is basically the 4xn equivalent with half the registers
 * and instructions removed.
 *
 * An alternative would be to load to different cluster of N j-particles
 * into SIMD registers, giving a 4x(N+N) kernel. This doubles the amount
 * of instructions, which could lead to better scheduling. But we actually
 * observed worse scheduling for the AVX-256 4x8 normal analytical PME
 * kernel, which has a lower pair throughput than 2x(4+4) with gcc 4.7.
 * It could be worth trying this option, but it takes some more effort.
 * This 2xnn kernel is basically the 4xn equivalent with
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
 * this can be faster with when we have defined gmx_blendv_pr, i.e. an instruction
 * that selects from two SIMD registers based on the contents of a third.
 */
#if !(defined CHECK_EXCLS || defined CALC_ENERGIES) && defined gmx_blendv_pr && !defined COUNT_PAIRS
/* The next two statements are hacks select optimization based on compiler and architecture indirectly based
 * on the instruction set supported. However, they only affect performance, not correctness,
 * and do not have to be altered for new architectures. They should be fixed eventually...
 */

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
    gmx_mm_pr  int_SIMD0;
    gmx_mm_pr  int_SIMD2;
#endif

    gmx_mm_pr  jxSIMD, jySIMD, jzSIMD;
    gmx_mm_pr  dx_SIMD0, dy_SIMD0, dz_SIMD0;
    gmx_mm_pr  dx_SIMD2, dy_SIMD2, dz_SIMD2;
    gmx_mm_pr  tx_SIMD0, ty_SIMD0, tz_SIMD0;
    gmx_mm_pr  tx_SIMD2, ty_SIMD2, tz_SIMD2;
    gmx_mm_pr  rsq_SIMD0, rinv_SIMD0, rinvsq_SIMD0;
    gmx_mm_pr  rsq_SIMD2, rinv_SIMD2, rinvsq_SIMD2;
#ifndef NBNXN_CUTOFF_USE_BLENDV
    /* wco: within cut-off, mask of all 1's or 0's */
    gmx_mm_pr  wco_SIMD0;
    gmx_mm_pr  wco_SIMD2;
#endif
#ifdef VDW_CUTOFF_CHECK
    gmx_mm_pr  wco_vdw_SIMD0;
#ifndef HALF_LJ
    gmx_mm_pr  wco_vdw_SIMD2;
#endif
#endif
#ifdef CALC_COULOMB
#ifdef CHECK_EXCLS
    /* 1/r masked with the interaction mask */
    gmx_mm_pr  rinv_ex_SIMD0;
    gmx_mm_pr  rinv_ex_SIMD2;
#endif
    gmx_mm_pr  jq_SIMD;
    gmx_mm_pr  qq_SIMD0;
    gmx_mm_pr  qq_SIMD2;
#ifdef CALC_COUL_TAB
    /* The force (PME mesh force) we need to subtract from 1/r^2 */
    gmx_mm_pr  fsub_SIMD0;
    gmx_mm_pr  fsub_SIMD2;
#endif
#ifdef CALC_COUL_EWALD
    gmx_mm_pr  brsq_SIMD0, brsq_SIMD2;
    gmx_mm_pr  ewcorr_SIMD0, ewcorr_SIMD2;
#endif

    /* frcoul = (1/r - fsub)*r */
    gmx_mm_pr  frcoul_SIMD0;
    gmx_mm_pr  frcoul_SIMD2;
#ifdef CALC_COUL_TAB
    /* For tables: r, rs=r/sp, rf=floor(rs), frac=rs-rf */
    gmx_mm_pr  r_SIMD0, rs_SIMD0, rf_SIMD0, frac_SIMD0;
    gmx_mm_pr  r_SIMD2, rs_SIMD2, rf_SIMD2, frac_SIMD2;
    /* Table index: rs truncated to an int */
    gmx_epi32  ti_SIMD0, ti_SIMD2;

    /* Linear force table values */
    gmx_mm_pr  ctab0_SIMD0, ctab1_SIMD0;
    gmx_mm_pr  ctab0_SIMD2, ctab1_SIMD2;
#ifdef CALC_ENERGIES
    /* Quadratic energy table value */
    gmx_mm_pr  ctabv_SIMD0;
    gmx_mm_pr  ctabv_SIMD2;
#endif
#endif
#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
    /* The potential (PME mesh) we need to subtract from 1/r */
    gmx_mm_pr  vc_sub_SIMD0;
    gmx_mm_pr  vc_sub_SIMD2;
#endif
#ifdef CALC_ENERGIES
    /* Electrostatic potential */
    gmx_mm_pr  vcoul_SIMD0;
    gmx_mm_pr  vcoul_SIMD2;
#endif
#endif
    /* The force times 1/r */
    gmx_mm_pr  fscal_SIMD0;
    gmx_mm_pr  fscal_SIMD2;

#ifdef CALC_LJ
#ifdef LJ_COMB_LB
    /* LJ sigma_j/2 and sqrt(epsilon_j) */
    gmx_mm_pr  hsig_j_SIMD, seps_j_SIMD;
    /* LJ sigma_ij and epsilon_ij */
    gmx_mm_pr  sig_SIMD0, eps_SIMD0;
#ifndef HALF_LJ
    gmx_mm_pr  sig_SIMD2, eps_SIMD2;
#endif
#ifdef CALC_ENERGIES
    gmx_mm_pr  sig2_SIMD0, sig6_SIMD0;
#ifndef HALF_LJ
    gmx_mm_pr  sig2_SIMD2, sig6_SIMD2;
#endif
#endif /* LJ_COMB_LB */
#endif /* CALC_LJ */

#ifdef LJ_COMB_GEOM
    gmx_mm_pr  c6s_j_SIMD, c12s_j_SIMD;
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB
    /* Index for loading LJ parameters, complicated when interleaving */
    int         aj2;
#endif

#ifndef FIX_LJ_C
    /* LJ C6 and C12 parameters, used with geometric comb. rule */
    gmx_mm_pr  c6_SIMD0, c12_SIMD0;
#ifndef HALF_LJ
    gmx_mm_pr  c6_SIMD2, c12_SIMD2;
#endif
#endif

    /* Intermediate variables for LJ calculation */
#ifndef LJ_COMB_LB
    gmx_mm_pr  rinvsix_SIMD0;
#ifndef HALF_LJ
    gmx_mm_pr  rinvsix_SIMD2;
#endif
#endif
#ifdef LJ_COMB_LB
    gmx_mm_pr  sir_SIMD0, sir2_SIMD0, sir6_SIMD0;
#ifndef HALF_LJ
    gmx_mm_pr  sir_SIMD2, sir2_SIMD2, sir6_SIMD2;
#endif
#endif

    gmx_mm_pr  FrLJ6_SIMD0, FrLJ12_SIMD0;
#ifndef HALF_LJ
    gmx_mm_pr  FrLJ6_SIMD2, FrLJ12_SIMD2;
#endif
#ifdef CALC_ENERGIES
    gmx_mm_pr  VLJ6_SIMD0, VLJ12_SIMD0, VLJ_SIMD0;
#ifndef HALF_LJ
    gmx_mm_pr  VLJ6_SIMD2, VLJ12_SIMD2, VLJ_SIMD2;
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
    {
        /* Load integer interaction mask */
        /* With AVX there are no integer operations, so cast to real */
        gmx_mm_pr mask_pr = gmx_mm_castsi256_pr(_mm256_set1_epi32(l_cj[cjind].excl));
        /* Intel Compiler version 12.1.3 20120130 is buggy: use cast.
         * With gcc we don't need the cast, but it's faster.
         */
        /* FIXME: Remove AVX-specific code in define to a separate routine in utils */
#define cast_cvt(x)  _mm256_cvtepi32_ps(_mm256_castps_si256(x))
        int_SIMD0  = gmx_cmpneq_pr(cast_cvt(gmx_and_pr(mask_pr, mask0)), zero_SIMD);
        int_SIMD2  = gmx_cmpneq_pr(cast_cvt(gmx_and_pr(mask_pr, mask2)), zero_SIMD);
#undef cast_cvt
    }
#endif
    /* load j atom coordinates */
    jxSIMD         = gmx_loaddh_pr(x+ajx);
    jySIMD         = gmx_loaddh_pr(x+ajy);
    jzSIMD         = gmx_loaddh_pr(x+ajz);

    /* Calculate distance */
    dx_SIMD0       = gmx_sub_pr(ix_SIMD0, jxSIMD);
    dy_SIMD0       = gmx_sub_pr(iy_SIMD0, jySIMD);
    dz_SIMD0       = gmx_sub_pr(iz_SIMD0, jzSIMD);
    dx_SIMD2       = gmx_sub_pr(ix_SIMD2, jxSIMD);
    dy_SIMD2       = gmx_sub_pr(iy_SIMD2, jySIMD);
    dz_SIMD2       = gmx_sub_pr(iz_SIMD2, jzSIMD);

    /* rsq = dx*dx+dy*dy+dz*dz */
    rsq_SIMD0      = gmx_calc_rsq_pr(dx_SIMD0, dy_SIMD0, dz_SIMD0);
    rsq_SIMD2      = gmx_calc_rsq_pr(dx_SIMD2, dy_SIMD2, dz_SIMD2);

#ifndef NBNXN_CUTOFF_USE_BLENDV
    wco_SIMD0      = gmx_cmplt_pr(rsq_SIMD0, rc2_SIMD);
    wco_SIMD2      = gmx_cmplt_pr(rsq_SIMD2, rc2_SIMD);
#endif

#ifdef CHECK_EXCLS
#ifdef EXCL_FORCES
    /* Only remove the (sub-)diagonal to avoid double counting */
#if UNROLLJ == UNROLLI
    if (cj == ci_sh)
    {
        wco_SIMD0  = gmx_and_pr(wco_SIMD0, diag_SIMD0);
        wco_SIMD2  = gmx_and_pr(wco_SIMD2, diag_SIMD2);
    }
#else
#error "only UNROLLJ == UNROLLI currently supported in the joined kernels"
#endif
#else /* EXCL_FORCES */
      /* Remove all excluded atom pairs from the list */
    wco_SIMD0      = gmx_and_pr(wco_SIMD0, int_SIMD0);
    wco_SIMD2      = gmx_and_pr(wco_SIMD2, int_SIMD2);
#endif
#endif

#ifdef COUNT_PAIRS
    {
        int  i, j;
        real tmp[UNROLLJ];
        for (i = 0; i < UNROLLI; i++)
        {
            gmx_storeu_pr(tmp, i == 0 ? wco_SIMD0 : (i == 1 ? wco_SIMD1 : (i == 2 ? wco_SIMD2 : wco_SIMD3)));
            for (j = 0; j < UNROLLJ; j++)
            {
                if (!(tmp[j] == 0))
                {
                    npair++;
                }
            }
        }
    }
#endif

#ifdef CHECK_EXCLS
    /* For excluded pairs add a small number to avoid r^-6 = NaN */
    rsq_SIMD0      = gmx_add_pr(rsq_SIMD0, gmx_andnot_pr(int_SIMD0, avoid_sing_SIMD));
    rsq_SIMD2      = gmx_add_pr(rsq_SIMD2, gmx_andnot_pr(int_SIMD2, avoid_sing_SIMD));
#endif

    /* Calculate 1/r */
    rinv_SIMD0     = gmx_invsqrt_pr(rsq_SIMD0);
    rinv_SIMD2     = gmx_invsqrt_pr(rsq_SIMD2);

#ifdef CALC_COULOMB
    /* Load parameters for j atom */
    jq_SIMD        = gmx_loaddh_pr(q+aj);
    qq_SIMD0       = gmx_mul_pr(iq_SIMD0, jq_SIMD);
    qq_SIMD2       = gmx_mul_pr(iq_SIMD2, jq_SIMD);
#endif

#ifdef CALC_LJ

#if !defined LJ_COMB_GEOM && !defined LJ_COMB_LB && !defined FIX_LJ_C
    load_lj_pair_params2(nbfp0, nbfp1, type, aj, c6_SIMD0, c12_SIMD0);
#ifndef HALF_LJ
    load_lj_pair_params2(nbfp2, nbfp3, type, aj, c6_SIMD2, c12_SIMD2);
#endif
#endif /* not defined any LJ rule */

#ifdef LJ_COMB_GEOM
    c6s_j_SIMD     = gmx_loaddh_pr(ljc+aj2+0);
    c12s_j_SIMD    = gmx_loaddh_pr(ljc+aj2+STRIDE);
    c6_SIMD0       = gmx_mul_pr(c6s_SIMD0, c6s_j_SIMD );
#ifndef HALF_LJ
    c6_SIMD2       = gmx_mul_pr(c6s_SIMD2, c6s_j_SIMD );
#endif
    c12_SIMD0      = gmx_mul_pr(c12s_SIMD0, c12s_j_SIMD);
#ifndef HALF_LJ
    c12_SIMD2      = gmx_mul_pr(c12s_SIMD2, c12s_j_SIMD);
#endif
#endif /* LJ_COMB_GEOM */

#ifdef LJ_COMB_LB
    hsig_j_SIMD    = gmx_loaddh_pr(ljc+aj2+0);
    seps_j_SIMD    = gmx_loaddh_pr(ljc+aj2+STRIDE);

    sig_SIMD0      = gmx_add_pr(hsig_i_SIMD0, hsig_j_SIMD);
    eps_SIMD0      = gmx_mul_pr(seps_i_SIMD0, seps_j_SIMD);
#ifndef HALF_LJ
    sig_SIMD2      = gmx_add_pr(hsig_i_SIMD2, hsig_j_SIMD);
    eps_SIMD2      = gmx_mul_pr(seps_i_SIMD2, seps_j_SIMD);
#endif
#endif /* LJ_COMB_LB */

#endif /* CALC_LJ */

#ifdef NBNXN_CUTOFF_USE_BLENDV
    /* We only need to mask for the cut-off: blendv is faster */
    rinv_SIMD0     = gmx_blendv_pr(rinv_SIMD0, zero_SIMD, gmx_sub_pr(rc2_SIMD, rsq_SIMD0));
    rinv_SIMD2     = gmx_blendv_pr(rinv_SIMD2, zero_SIMD, gmx_sub_pr(rc2_SIMD, rsq_SIMD2));
#else
    rinv_SIMD0     = gmx_and_pr(rinv_SIMD0, wco_SIMD0);
    rinv_SIMD2     = gmx_and_pr(rinv_SIMD2, wco_SIMD2);
#endif

    rinvsq_SIMD0   = gmx_mul_pr(rinv_SIMD0, rinv_SIMD0);
    rinvsq_SIMD2   = gmx_mul_pr(rinv_SIMD2, rinv_SIMD2);

#ifdef CALC_COULOMB
    /* Note that here we calculate force*r, not the usual force/r.
     * This allows avoiding masking the reaction-field contribution,
     * as frcoul is later multiplied by rinvsq which has been
     * masked with the cut-off check.
     */

#ifdef EXCL_FORCES
    /* Only add 1/r for non-excluded atom pairs */
    rinv_ex_SIMD0  = gmx_and_pr(rinv_SIMD0, int_SIMD0);
    rinv_ex_SIMD2  = gmx_and_pr(rinv_SIMD2, int_SIMD2);
#else
    /* No exclusion forces, we always need 1/r */
#define     rinv_ex_SIMD0    rinv_SIMD0
#define     rinv_ex_SIMD2    rinv_SIMD2
#endif

#ifdef CALC_COUL_RF
    /* Electrostatic interactions */
    frcoul_SIMD0   = gmx_mul_pr(qq_SIMD0, gmx_add_pr(rinv_ex_SIMD0, gmx_mul_pr(rsq_SIMD0, mrc_3_SIMD)));
    frcoul_SIMD2   = gmx_mul_pr(qq_SIMD2, gmx_add_pr(rinv_ex_SIMD2, gmx_mul_pr(rsq_SIMD2, mrc_3_SIMD)));

#ifdef CALC_ENERGIES
    vcoul_SIMD0    = gmx_mul_pr(qq_SIMD0, gmx_add_pr(rinv_ex_SIMD0, gmx_add_pr(gmx_mul_pr(rsq_SIMD0, hrc_3_SIMD), moh_rc_SIMD)));
    vcoul_SIMD2    = gmx_mul_pr(qq_SIMD2, gmx_add_pr(rinv_ex_SIMD2, gmx_add_pr(gmx_mul_pr(rsq_SIMD2, hrc_3_SIMD), moh_rc_SIMD)));
#endif
#endif

#ifdef CALC_COUL_EWALD
    /* We need to mask (or limit) rsq for the cut-off,
     * as large distances can cause an overflow in gmx_pmecorrF/V.
     */
#ifdef NBNXN_CUTOFF_USE_BLENDV
    /* Strangely, putting mul on a separate line is slower (icc 13) */
    brsq_SIMD0     = gmx_mul_pr(beta2_SIMD, gmx_blendv_pr(rsq_SIMD0, zero_SIMD, gmx_sub_pr(rc2_SIMD, rsq_SIMD0)));
    brsq_SIMD2     = gmx_mul_pr(beta2_SIMD, gmx_blendv_pr(rsq_SIMD2, zero_SIMD, gmx_sub_pr(rc2_SIMD, rsq_SIMD2)));
#else
    brsq_SIMD0     = gmx_mul_pr(beta2_SIMD, gmx_and_pr(rsq_SIMD0, wco_SIMD0));
    brsq_SIMD2     = gmx_mul_pr(beta2_SIMD, gmx_and_pr(rsq_SIMD2, wco_SIMD2));
#endif
    ewcorr_SIMD0   = gmx_mul_pr(gmx_pmecorrF_pr(brsq_SIMD0), beta_SIMD);
    ewcorr_SIMD2   = gmx_mul_pr(gmx_pmecorrF_pr(brsq_SIMD2), beta_SIMD);
    frcoul_SIMD0   = gmx_mul_pr(qq_SIMD0, gmx_add_pr(rinv_ex_SIMD0, gmx_mul_pr(ewcorr_SIMD0, brsq_SIMD0)));
    frcoul_SIMD2   = gmx_mul_pr(qq_SIMD2, gmx_add_pr(rinv_ex_SIMD2, gmx_mul_pr(ewcorr_SIMD2, brsq_SIMD2)));

#ifdef CALC_ENERGIES
    vc_sub_SIMD0   = gmx_mul_pr(gmx_pmecorrV_pr(brsq_SIMD0), beta_SIMD);
    vc_sub_SIMD2   = gmx_mul_pr(gmx_pmecorrV_pr(brsq_SIMD2), beta_SIMD);
#endif

#endif /* CALC_COUL_EWALD */

#ifdef CALC_COUL_TAB
    /* Electrostatic interactions */
    r_SIMD0        = gmx_mul_pr(rsq_SIMD0, rinv_SIMD0);
    r_SIMD2        = gmx_mul_pr(rsq_SIMD2, rinv_SIMD2);
    /* Convert r to scaled table units */
    rs_SIMD0       = gmx_mul_pr(r_SIMD0, invtsp_SIMD);
    rs_SIMD2       = gmx_mul_pr(r_SIMD2, invtsp_SIMD);
    /* Truncate scaled r to an int */
    ti_SIMD0       = gmx_cvttpr_epi32(rs_SIMD0);
    ti_SIMD2       = gmx_cvttpr_epi32(rs_SIMD2);
#ifdef gmx_floor_pr
    /* SIMD4.1 floor is faster than gmx_cvtepi32_ps int->float cast */
    rf_SIMD0       = gmx_floor_pr(rs_SIMD0);
    rf_SIMD2       = gmx_floor_pr(rs_SIMD2);
#else
    rf_SIMD0       = gmx_cvtepi32_pr(ti_SIMD0);
    rf_SIMD2       = gmx_cvtepi32_pr(ti_SIMD2);
#endif
    frac_SIMD0     = gmx_sub_pr(rs_SIMD0, rf_SIMD0);
    frac_SIMD2     = gmx_sub_pr(rs_SIMD2, rf_SIMD2);

    /* Load and interpolate table forces and possibly energies.
     * Force and energy can be combined in one table, stride 4: FDV0
     * or in two separate tables with stride 1: F and V
     * Currently single precision uses FDV0, double F and V.
     */
#ifdef CALC_ENERGIES
#ifdef TAB_FDV0
    load_table_f_v(tab_coul_F, ti_SIMD0, ti0, ctab0_SIMD0, ctab1_SIMD0, ctabv_SIMD0);
    load_table_f_v(tab_coul_F, ti_SIMD2, ti2, ctab0_SIMD2, ctab1_SIMD2, ctabv_SIMD2);
#else
    load_table_f_v(tab_coul_F, tab_coul_V, ti_SIMD0, ti0, ctab0_SIMD0, ctab1_SIMD0, ctabv_SIMD0);
    load_table_f_v(tab_coul_F, tab_coul_V, ti_SIMD2, ti2, ctab0_SIMD2, ctab1_SIMD2, ctabv_SIMD2);
#endif
#else
    load_table_f(tab_coul_F, ti_SIMD0, ti0, ctab0_SIMD0, ctab1_SIMD0);
    load_table_f(tab_coul_F, ti_SIMD2, ti2, ctab0_SIMD2, ctab1_SIMD2);
#endif
    fsub_SIMD0     = gmx_add_pr(ctab0_SIMD0, gmx_mul_pr(frac_SIMD0, ctab1_SIMD0));
    fsub_SIMD2     = gmx_add_pr(ctab0_SIMD2, gmx_mul_pr(frac_SIMD2, ctab1_SIMD2));
    frcoul_SIMD0   = gmx_mul_pr(qq_SIMD0, gmx_sub_pr(rinv_ex_SIMD0, gmx_mul_pr(fsub_SIMD0, r_SIMD0)));
    frcoul_SIMD2   = gmx_mul_pr(qq_SIMD2, gmx_sub_pr(rinv_ex_SIMD2, gmx_mul_pr(fsub_SIMD2, r_SIMD2)));

#ifdef CALC_ENERGIES
    vc_sub_SIMD0   = gmx_add_pr(ctabv_SIMD0, gmx_mul_pr(gmx_mul_pr(mhalfsp_SIMD, frac_SIMD0), gmx_add_pr(ctab0_SIMD0, fsub_SIMD0)));
    vc_sub_SIMD2   = gmx_add_pr(ctabv_SIMD2, gmx_mul_pr(gmx_mul_pr(mhalfsp_SIMD, frac_SIMD2), gmx_add_pr(ctab0_SIMD2, fsub_SIMD2)));
#endif
#endif /* CALC_COUL_TAB */

#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
#ifndef NO_SHIFT_EWALD
    /* Add Ewald potential shift to vc_sub for convenience */
#ifdef CHECK_EXCLS
    vc_sub_SIMD0   = gmx_add_pr(vc_sub_SIMD0, gmx_and_pr(sh_ewald_SIMD, int_SIMD0));
    vc_sub_SIMD2   = gmx_add_pr(vc_sub_SIMD2, gmx_and_pr(sh_ewald_SIMD, int_SIMD2));
#else
    vc_sub_SIMD0   = gmx_add_pr(vc_sub_SIMD0, sh_ewald_SIMD);
    vc_sub_SIMD2   = gmx_add_pr(vc_sub_SIMD2, sh_ewald_SIMD);
#endif
#endif

    vcoul_SIMD0    = gmx_mul_pr(qq_SIMD0, gmx_sub_pr(rinv_ex_SIMD0, vc_sub_SIMD0));
    vcoul_SIMD2    = gmx_mul_pr(qq_SIMD2, gmx_sub_pr(rinv_ex_SIMD2, vc_sub_SIMD2));
#endif

#ifdef CALC_ENERGIES
    /* Mask energy for cut-off and diagonal */
    vcoul_SIMD0    = gmx_and_pr(vcoul_SIMD0, wco_SIMD0);
    vcoul_SIMD2    = gmx_and_pr(vcoul_SIMD2, wco_SIMD2);
#endif

#endif /* CALC_COULOMB */

#ifdef CALC_LJ
    /* Lennard-Jones interaction */

#ifdef VDW_CUTOFF_CHECK
    wco_vdw_SIMD0  = gmx_cmplt_pr(rsq_SIMD0, rcvdw2_SIMD);
#ifndef HALF_LJ
    wco_vdw_SIMD2  = gmx_cmplt_pr(rsq_SIMD2, rcvdw2_SIMD);
#endif
#else
    /* Same cut-off for Coulomb and VdW, reuse the registers */
#define     wco_vdw_SIMD0    wco_SIMD0
#define     wco_vdw_SIMD2    wco_SIMD2
#endif

#ifndef LJ_COMB_LB
    rinvsix_SIMD0  = gmx_mul_pr(rinvsq_SIMD0, gmx_mul_pr(rinvsq_SIMD0, rinvsq_SIMD0));
#ifdef EXCL_FORCES
    rinvsix_SIMD0  = gmx_and_pr(rinvsix_SIMD0, int_SIMD0);
#endif
#ifndef HALF_LJ
    rinvsix_SIMD2  = gmx_mul_pr(rinvsq_SIMD2, gmx_mul_pr(rinvsq_SIMD2, rinvsq_SIMD2));
#ifdef EXCL_FORCES
    rinvsix_SIMD2  = gmx_and_pr(rinvsix_SIMD2, int_SIMD2);
#endif
#endif
#ifdef VDW_CUTOFF_CHECK
    rinvsix_SIMD0  = gmx_and_pr(rinvsix_SIMD0, wco_vdw_SIMD0);
#ifndef HALF_LJ
    rinvsix_SIMD2  = gmx_and_pr(rinvsix_SIMD2, wco_vdw_SIMD2);
#endif
#endif
    FrLJ6_SIMD0    = gmx_mul_pr(c6_SIMD0, rinvsix_SIMD0);
#ifndef HALF_LJ
    FrLJ6_SIMD2    = gmx_mul_pr(c6_SIMD2, rinvsix_SIMD2);
#endif
    FrLJ12_SIMD0   = gmx_mul_pr(c12_SIMD0, gmx_mul_pr(rinvsix_SIMD0, rinvsix_SIMD0));
#ifndef HALF_LJ
    FrLJ12_SIMD2   = gmx_mul_pr(c12_SIMD2, gmx_mul_pr(rinvsix_SIMD2, rinvsix_SIMD2));
#endif
#endif /* not LJ_COMB_LB */

#ifdef LJ_COMB_LB
    sir_SIMD0      = gmx_mul_pr(sig_SIMD0, rinv_SIMD0);
#ifndef HALF_LJ
    sir_SIMD2      = gmx_mul_pr(sig_SIMD2, rinv_SIMD2);
#endif
    sir2_SIMD0     = gmx_mul_pr(sir_SIMD0, sir_SIMD0);
#ifndef HALF_LJ
    sir2_SIMD2     = gmx_mul_pr(sir_SIMD2, sir_SIMD2);
#endif
    sir6_SIMD0     = gmx_mul_pr(sir2_SIMD0, gmx_mul_pr(sir2_SIMD0, sir2_SIMD0));
#ifdef EXCL_FORCES
    sir6_SIMD0     = gmx_and_pr(sir6_SIMD0, int_SIMD0);
#endif
#ifndef HALF_LJ
    sir6_SIMD2     = gmx_mul_pr(sir2_SIMD2, gmx_mul_pr(sir2_SIMD2, sir2_SIMD2));
#ifdef EXCL_FORCES
    sir6_SIMD2     = gmx_and_pr(sir6_SIMD2, int_SIMD2);
#endif
#endif
#ifdef VDW_CUTOFF_CHECK
    sir6_SIMD0     = gmx_and_pr(sir6_SIMD0, wco_vdw_SIMD0);
#ifndef HALF_LJ
    sir6_SIMD2     = gmx_and_pr(sir6_SIMD2, wco_vdw_SIMD2);
#endif
#endif
    FrLJ6_SIMD0    = gmx_mul_pr(eps_SIMD0, sir6_SIMD0);
#ifndef HALF_LJ
    FrLJ6_SIMD2    = gmx_mul_pr(eps_SIMD2, sir6_SIMD2);
#endif
    FrLJ12_SIMD0   = gmx_mul_pr(FrLJ6_SIMD0, sir6_SIMD0);
#ifndef HALF_LJ
    FrLJ12_SIMD2   = gmx_mul_pr(FrLJ6_SIMD2, sir6_SIMD2);
#endif
#if defined CALC_ENERGIES
    /* We need C6 and C12 to calculate the LJ potential shift */
    sig2_SIMD0     = gmx_mul_pr(sig_SIMD0, sig_SIMD0);
#ifndef HALF_LJ
    sig2_SIMD2     = gmx_mul_pr(sig_SIMD2, sig_SIMD2);
#endif
    sig6_SIMD0     = gmx_mul_pr(sig2_SIMD0, gmx_mul_pr(sig2_SIMD0, sig2_SIMD0));
#ifndef HALF_LJ
    sig6_SIMD2     = gmx_mul_pr(sig2_SIMD2, gmx_mul_pr(sig2_SIMD2, sig2_SIMD2));
#endif
    c6_SIMD0       = gmx_mul_pr(eps_SIMD0, sig6_SIMD0);
#ifndef HALF_LJ
    c6_SIMD2       = gmx_mul_pr(eps_SIMD2, sig6_SIMD2);
#endif
    c12_SIMD0      = gmx_mul_pr(c6_SIMD0, sig6_SIMD0);
#ifndef HALF_LJ
    c12_SIMD2      = gmx_mul_pr(c6_SIMD2, sig6_SIMD2);
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
#ifdef ENERGY_GROUPS
    add_ener_grp_halves(vcoul_SIMD0, vctp[0], vctp[1], egp_jj);
    add_ener_grp_halves(vcoul_SIMD2, vctp[2], vctp[3], egp_jj);
#else
    vctotSIMD      = gmx_add_pr(vctotSIMD, gmx_add_pr(vcoul_SIMD0, vcoul_SIMD2));
#endif
#endif

#ifdef CALC_LJ
    /* Calculate the LJ energies */
    VLJ6_SIMD0     = gmx_mul_pr(sixthSIMD, gmx_sub_pr(FrLJ6_SIMD0, gmx_mul_pr(c6_SIMD0, sh_invrc6_SIMD)));
#ifndef HALF_LJ
    VLJ6_SIMD2     = gmx_mul_pr(sixthSIMD, gmx_sub_pr(FrLJ6_SIMD2, gmx_mul_pr(c6_SIMD2, sh_invrc6_SIMD)));
#endif
    VLJ12_SIMD0    = gmx_mul_pr(twelvethSIMD, gmx_sub_pr(FrLJ12_SIMD0, gmx_mul_pr(c12_SIMD0, sh_invrc12_SIMD)));
#ifndef HALF_LJ
    VLJ12_SIMD2    = gmx_mul_pr(twelvethSIMD, gmx_sub_pr(FrLJ12_SIMD2, gmx_mul_pr(c12_SIMD2, sh_invrc12_SIMD)));
#endif

    VLJ_SIMD0      = gmx_sub_pr(VLJ12_SIMD0, VLJ6_SIMD0);
#ifndef HALF_LJ
    VLJ_SIMD2      = gmx_sub_pr(VLJ12_SIMD2, VLJ6_SIMD2);
#endif
    /* The potential shift should be removed for pairs beyond cut-off */
    VLJ_SIMD0      = gmx_and_pr(VLJ_SIMD0, wco_vdw_SIMD0);
#ifndef HALF_LJ
    VLJ_SIMD2      = gmx_and_pr(VLJ_SIMD2, wco_vdw_SIMD2);
#endif
#ifdef CHECK_EXCLS
    /* The potential shift should be removed for excluded pairs */
    VLJ_SIMD0      = gmx_and_pr(VLJ_SIMD0, int_SIMD0);
#ifndef HALF_LJ
    VLJ_SIMD2      = gmx_and_pr(VLJ_SIMD2, int_SIMD2);
#endif
#endif
#ifdef ENERGY_GROUPS
    add_ener_grp_halves(VLJ_SIMD0, vvdwtp[0], vvdwtp[1], egp_jj);
#ifndef HALF_LJ
    add_ener_grp_halves(VLJ_SIMD2, vvdwtp[2], vvdwtp[3], egp_jj);
#endif
#else
    VvdwtotSIMD    = gmx_add_pr(VvdwtotSIMD,
#ifdef HALF_LJ
                                VLJ_SIMD0
#else
                                gmx_add_pr(VLJ_SIMD0, VLJ_SIMD2)
#endif
                                );
#endif
#endif /* CALC_LJ */
#endif /* CALC_ENERGIES */

#ifdef CALC_LJ
    fscal_SIMD0    = gmx_mul_pr(rinvsq_SIMD0,
#ifdef CALC_COULOMB
                               gmx_add_pr(frcoul_SIMD0,
#else
                               (
#endif
                                          gmx_sub_pr(FrLJ12_SIMD0, FrLJ6_SIMD0)));
#else
    fscal_SIMD0    = gmx_mul_pr(rinvsq_SIMD0, frcoul_SIMD0);
#endif /* CALC_LJ */
#if defined CALC_LJ && !defined HALF_LJ
    fscal_SIMD2    = gmx_mul_pr(rinvsq_SIMD2,
#ifdef CALC_COULOMB
                               gmx_add_pr(frcoul_SIMD2,
#else
                               (
#endif
                                          gmx_sub_pr(FrLJ12_SIMD2, FrLJ6_SIMD2)));
#else
    /* Atom 2 and 3 don't have LJ, so only add Coulomb forces */
    fscal_SIMD2    = gmx_mul_pr(rinvsq_SIMD2, frcoul_SIMD2);
#endif

    /* Calculate temporary vectorial force */
    tx_SIMD0       = gmx_mul_pr(fscal_SIMD0, dx_SIMD0);
    tx_SIMD2       = gmx_mul_pr(fscal_SIMD2, dx_SIMD2);
    ty_SIMD0       = gmx_mul_pr(fscal_SIMD0, dy_SIMD0);
    ty_SIMD2       = gmx_mul_pr(fscal_SIMD2, dy_SIMD2);
    tz_SIMD0       = gmx_mul_pr(fscal_SIMD0, dz_SIMD0);
    tz_SIMD2       = gmx_mul_pr(fscal_SIMD2, dz_SIMD2);

    /* Increment i atom force */
    fix_SIMD0      = gmx_add_pr(fix_SIMD0, tx_SIMD0);
    fix_SIMD2      = gmx_add_pr(fix_SIMD2, tx_SIMD2);
    fiy_SIMD0      = gmx_add_pr(fiy_SIMD0, ty_SIMD0);
    fiy_SIMD2      = gmx_add_pr(fiy_SIMD2, ty_SIMD2);
    fiz_SIMD0      = gmx_add_pr(fiz_SIMD0, tz_SIMD0);
    fiz_SIMD2      = gmx_add_pr(fiz_SIMD2, tz_SIMD2);

    /* Decrement j atom force */
    gmx_store_hpr(f+ajx,
                  gmx_sub_hpr( gmx_load_hpr(f+ajx), gmx_sum4_hpr(tx_SIMD0, tx_SIMD2) ));
    gmx_store_hpr(f+ajy,
                  gmx_sub_hpr( gmx_load_hpr(f+ajy), gmx_sum4_hpr(ty_SIMD0, ty_SIMD2) ));
    gmx_store_hpr(f+ajz,
                  gmx_sub_hpr( gmx_load_hpr(f+ajz), gmx_sum4_hpr(tz_SIMD0, tz_SIMD2) ));
}

#undef  rinv_ex_SIMD0
#undef  rinv_ex_SIMD2

#undef  wco_vdw_SIMD0
#undef  wco_vdw_SIMD2

#undef  NBNXN_CUTOFF_USE_BLENDV

#undef  EXCL_FORCES
