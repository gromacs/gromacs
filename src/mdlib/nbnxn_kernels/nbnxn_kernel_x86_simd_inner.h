/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */

/* This is the innermost loop contents for the n vs n atom
 * SSE2 single precision kernels.
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

#if !(defined CHECK_EXCLS || defined CALC_ENERGIES) && defined GMX_X86_SSE4_1 && !defined COUNT_PAIRS && !defined __GNUC__
/* Without exclusions and energies we only need to mask the cut-off,
 * this is faster with blendv (only available with SSE4.1 and later).
 * With gcc blendv is much slower; tested with gcc 4.6.2 and 4.6.3.
 */
#define CUTOFF_BLENDV
#endif

        {
            int        cj,aj,ajx,ajy,ajz;

#ifdef ENERGY_GROUPS
            int        egps_j;
            int        egp_jj[UNROLLJ>>1];
            int        jj;
#endif

#ifdef CHECK_EXCLS
            /* Interaction (non-exclusion) mask of all 1's or 0's */
            gmx_mm_pr  int_SSE0;
            gmx_mm_pr  int_SSE1;
            gmx_mm_pr  int_SSE2;
            gmx_mm_pr  int_SSE3;
#endif

            gmx_mm_pr  jxSSE,jySSE,jzSSE;
            gmx_mm_pr  dx_SSE0,dy_SSE0,dz_SSE0;
            gmx_mm_pr  dx_SSE1,dy_SSE1,dz_SSE1;
            gmx_mm_pr  dx_SSE2,dy_SSE2,dz_SSE2;
            gmx_mm_pr  dx_SSE3,dy_SSE3,dz_SSE3;
            gmx_mm_pr  tx_SSE0,ty_SSE0,tz_SSE0;
            gmx_mm_pr  tx_SSE1,ty_SSE1,tz_SSE1;
            gmx_mm_pr  tx_SSE2,ty_SSE2,tz_SSE2;
            gmx_mm_pr  tx_SSE3,ty_SSE3,tz_SSE3;
            gmx_mm_pr  rsq_SSE0,rinv_SSE0,rinvsq_SSE0;
            gmx_mm_pr  rsq_SSE1,rinv_SSE1,rinvsq_SSE1;
            gmx_mm_pr  rsq_SSE2,rinv_SSE2,rinvsq_SSE2;
            gmx_mm_pr  rsq_SSE3,rinv_SSE3,rinvsq_SSE3;
#ifndef CUTOFF_BLENDV
            /* wco: within cut-off, mask of all 1's or 0's */
            gmx_mm_pr  wco_SSE0;
            gmx_mm_pr  wco_SSE1;
            gmx_mm_pr  wco_SSE2;
            gmx_mm_pr  wco_SSE3;
#endif
#ifdef CALC_COULOMB
#ifdef CHECK_EXCLS
            /* 1/r masked with the interaction mask */
            gmx_mm_pr  rinv_ex_SSE0;
            gmx_mm_pr  rinv_ex_SSE1;
            gmx_mm_pr  rinv_ex_SSE2;
            gmx_mm_pr  rinv_ex_SSE3;
#endif
            gmx_mm_pr  jq_SSE;
            gmx_mm_pr  qq_SSE0;
            gmx_mm_pr  qq_SSE1;
            gmx_mm_pr  qq_SSE2;
            gmx_mm_pr  qq_SSE3;
#ifdef CALC_COUL_TAB
            /* The force (PME mesh force) we need to subtract from 1/r^2 */
            gmx_mm_pr  fsub_SSE0;
            gmx_mm_pr  fsub_SSE1;
            gmx_mm_pr  fsub_SSE2;
            gmx_mm_pr  fsub_SSE3;
#endif
            /* frcoul = (1/r - fsub)*r */
            gmx_mm_pr  frcoul_SSE0;
            gmx_mm_pr  frcoul_SSE1;
            gmx_mm_pr  frcoul_SSE2;
            gmx_mm_pr  frcoul_SSE3;
#ifdef CALC_COUL_TAB
            /* For tables: r, rs=r/sp, rf=floor(rs), frac=rs-rf */
            gmx_mm_pr  r_SSE0,rs_SSE0,rf_SSE0,frac_SSE0;
            gmx_mm_pr  r_SSE1,rs_SSE1,rf_SSE1,frac_SSE1;
            gmx_mm_pr  r_SSE2,rs_SSE2,rf_SSE2,frac_SSE2;
            gmx_mm_pr  r_SSE3,rs_SSE3,rf_SSE3,frac_SSE3;
            /* Table index: rs converted to an int */ 
#if !(defined GMX_MM256_HERE && defined GMX_DOUBLE)
            gmx_epi32  ti_SSE0,ti_SSE1,ti_SSE2,ti_SSE3;
#else
            __m128i    ti_SSE0,ti_SSE1,ti_SSE2,ti_SSE3;
#endif
            /* Linear force table values */
            gmx_mm_pr  ctab0_SSE0,ctab1_SSE0;
            gmx_mm_pr  ctab0_SSE1,ctab1_SSE1;
            gmx_mm_pr  ctab0_SSE2,ctab1_SSE2;
            gmx_mm_pr  ctab0_SSE3,ctab1_SSE3;
#ifdef CALC_ENERGIES
            /* Quadratic energy table value */
            gmx_mm_pr  ctabv_SSE0;
            gmx_mm_pr  ctabv_SSE1;
            gmx_mm_pr  ctabv_SSE2;
            gmx_mm_pr  ctabv_SSE3;
            /* The potential (PME mesh) we need to subtract from 1/r */
            gmx_mm_pr  vc_sub_SSE0;
            gmx_mm_pr  vc_sub_SSE1;
            gmx_mm_pr  vc_sub_SSE2;
            gmx_mm_pr  vc_sub_SSE3;
#endif
#endif
#ifdef CALC_ENERGIES
            /* Electrostatic potential */
            gmx_mm_pr  vcoul_SSE0;
            gmx_mm_pr  vcoul_SSE1;
            gmx_mm_pr  vcoul_SSE2;
            gmx_mm_pr  vcoul_SSE3;
#endif
#endif
            /* The force times 1/r */
            gmx_mm_pr  fscal_SSE0;
            gmx_mm_pr  fscal_SSE1;
            gmx_mm_pr  fscal_SSE2;
            gmx_mm_pr  fscal_SSE3;

#ifdef CALC_LJ
#ifdef LJ_COMB_LB
            /* LJ sigma_j/2 and sqrt(epsilon_j) */
            gmx_mm_pr  hsig_j_SSE,seps_j_SSE;
            /* LJ sigma_ij and epsilon_ij */
            gmx_mm_pr  sig_SSE0,eps_SSE0;
            gmx_mm_pr  sig_SSE1,eps_SSE1;
#ifndef HALF_LJ
            gmx_mm_pr  sig_SSE2,eps_SSE2;
            gmx_mm_pr  sig_SSE3,eps_SSE3;
#endif
#ifdef CALC_ENERGIES
            gmx_mm_pr  sig6_SSE0;
            gmx_mm_pr  sig6_SSE1;
#ifndef HALF_LJ
            gmx_mm_pr  sig6_SSE2;
            gmx_mm_pr  sig6_SSE3;
#endif
#endif
#endif

#ifdef LJ_COMB_GEOM
            gmx_mm_pr  c6s_j_SSE,c12s_j_SSE;
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB
            /* Index for loading LJ parameters, complicated when interleaving */
            int         aj2;
#endif

#ifndef FIX_LJ_C
            /* LJ C6 and C12 parameters, used with geometric comb. rule */
            gmx_mm_pr  c6_SSE0,c12_SSE0;
            gmx_mm_pr  c6_SSE1,c12_SSE1;
#ifndef HALF_LJ
            gmx_mm_pr  c6_SSE2,c12_SSE2;
            gmx_mm_pr  c6_SSE3,c12_SSE3;
#endif
#endif

            /* Intermediate variables for LJ calculation */
#ifdef LJ_COMB_LB
            gmx_mm_pr  sir_SSE0,sir2_SSE0,sir6_SSE0;
            gmx_mm_pr  sir_SSE1,sir2_SSE1,sir6_SSE1;
#ifndef HALF_LJ
            gmx_mm_pr  sir_SSE2,sir2_SSE2,sir6_SSE2;
            gmx_mm_pr  sir_SSE3,sir2_SSE3,sir6_SSE3;
#endif
#else
            gmx_mm_pr  rinvsix_SSE0;
            gmx_mm_pr  rinvsix_SSE1;
#ifndef HALF_LJ
            gmx_mm_pr  rinvsix_SSE2;
            gmx_mm_pr  rinvsix_SSE3;
#endif
#endif

            gmx_mm_pr  FrLJ6_SSE0,FrLJ12_SSE0;
            gmx_mm_pr  FrLJ6_SSE1,FrLJ12_SSE1;
#ifndef HALF_LJ
            gmx_mm_pr  FrLJ6_SSE2,FrLJ12_SSE2;
            gmx_mm_pr  FrLJ6_SSE3,FrLJ12_SSE3;
#endif
#ifdef CALC_ENERGIES
            gmx_mm_pr  VLJ6_SSE0,VLJ12_SSE0,VLJ_SSE0;
            gmx_mm_pr  VLJ6_SSE1,VLJ12_SSE1,VLJ_SSE1;
#ifndef HALF_LJ
            gmx_mm_pr  VLJ6_SSE2,VLJ12_SSE2,VLJ_SSE2;
            gmx_mm_pr  VLJ6_SSE3,VLJ12_SSE3,VLJ_SSE3;
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
#ifndef GMX_MM256_HERE
            {
                /* Load integer interaction mask */
                __m128i mask_int = _mm_set1_epi32(l_cj[cjind].excl);

                /* The is no unequal sse instruction, so we need a not here */
                int_SSE0  = gmx_mm_castsi128_pr(_mm_cmpeq_epi32(_mm_andnot_si128(mask_int,mask0),zeroi_SSE));
                int_SSE1  = gmx_mm_castsi128_pr(_mm_cmpeq_epi32(_mm_andnot_si128(mask_int,mask1),zeroi_SSE));
                int_SSE2  = gmx_mm_castsi128_pr(_mm_cmpeq_epi32(_mm_andnot_si128(mask_int,mask2),zeroi_SSE));
                int_SSE3  = gmx_mm_castsi128_pr(_mm_cmpeq_epi32(_mm_andnot_si128(mask_int,mask3),zeroi_SSE));
            }
#else
            {
#ifndef GMX_DOUBLE
                /* Load integer interaction mask */
                /* With AVX there are no integer operations, so cast to real */
                gmx_mm_pr mask_pr = gmx_mm_castsi256_pr(_mm256_set1_epi32(l_cj[cjind].excl));
                /* We can't compare all 4*8=32 float bits: shift the mask */
                gmx_mm_pr masksh_pr = gmx_mm_castsi256_pr(_mm256_set1_epi32(l_cj[cjind].excl>>(2*UNROLLJ)));
                /* Intel Compiler version 12.1.3 20120130 is buggy: use cast.
                 * With gcc we don't need the cast, but it's faster.
                 */
#define cast_cvt(x)  _mm256_cvtepi32_ps(_mm256_castps_si256(x))
                int_SSE0  = gmx_cmpneq_pr(cast_cvt(gmx_and_pr(mask_pr,mask0)),zero_SSE);
                int_SSE1  = gmx_cmpneq_pr(cast_cvt(gmx_and_pr(mask_pr,mask1)),zero_SSE);
                int_SSE2  = gmx_cmpneq_pr(cast_cvt(gmx_and_pr(masksh_pr,mask0)),zero_SSE);
                int_SSE3  = gmx_cmpneq_pr(cast_cvt(gmx_and_pr(masksh_pr,mask1)),zero_SSE);
#undef cast_cvt
#else
                /* Load integer interaction mask */
                /* With AVX there are no integer operations,
                 * and there is no int to double conversion, so cast to float
                 */
                __m256 mask_ps = _mm256_castsi256_ps(_mm256_set1_epi32(l_cj[cjind].excl));
#define cast_cvt(x)  _mm256_castps_pd(_mm256_cvtepi32_ps(_mm256_castps_si256(x)))
                int_SSE0  = gmx_cmpneq_pr(cast_cvt(_mm256_and_ps(mask_ps,mask0)),zero_SSE);
                int_SSE1  = gmx_cmpneq_pr(cast_cvt(_mm256_and_ps(mask_ps,mask1)),zero_SSE);
                int_SSE2  = gmx_cmpneq_pr(cast_cvt(_mm256_and_ps(mask_ps,mask2)),zero_SSE);
                int_SSE3  = gmx_cmpneq_pr(cast_cvt(_mm256_and_ps(mask_ps,mask3)),zero_SSE);
#undef cast_cvt
#endif
            }
#endif
#endif
            /* load j atom coordinates */
            jxSSE         = gmx_load_pr(x+ajx);
            jySSE         = gmx_load_pr(x+ajy);
            jzSSE         = gmx_load_pr(x+ajz);
            
            /* Calculate distance */
            dx_SSE0       = gmx_sub_pr(ix_SSE0,jxSSE);
            dy_SSE0       = gmx_sub_pr(iy_SSE0,jySSE);
            dz_SSE0       = gmx_sub_pr(iz_SSE0,jzSSE);
            dx_SSE1       = gmx_sub_pr(ix_SSE1,jxSSE);
            dy_SSE1       = gmx_sub_pr(iy_SSE1,jySSE);
            dz_SSE1       = gmx_sub_pr(iz_SSE1,jzSSE);
            dx_SSE2       = gmx_sub_pr(ix_SSE2,jxSSE);
            dy_SSE2       = gmx_sub_pr(iy_SSE2,jySSE);
            dz_SSE2       = gmx_sub_pr(iz_SSE2,jzSSE);
            dx_SSE3       = gmx_sub_pr(ix_SSE3,jxSSE);
            dy_SSE3       = gmx_sub_pr(iy_SSE3,jySSE);
            dz_SSE3       = gmx_sub_pr(iz_SSE3,jzSSE);
            
            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0      = gmx_calc_rsq_pr(dx_SSE0,dy_SSE0,dz_SSE0);
            rsq_SSE1      = gmx_calc_rsq_pr(dx_SSE1,dy_SSE1,dz_SSE1);
            rsq_SSE2      = gmx_calc_rsq_pr(dx_SSE2,dy_SSE2,dz_SSE2);
            rsq_SSE3      = gmx_calc_rsq_pr(dx_SSE3,dy_SSE3,dz_SSE3);

#ifndef CUTOFF_BLENDV
            wco_SSE0      = gmx_cmplt_pr(rsq_SSE0,rc2_SSE);
            wco_SSE1      = gmx_cmplt_pr(rsq_SSE1,rc2_SSE);
            wco_SSE2      = gmx_cmplt_pr(rsq_SSE2,rc2_SSE);
            wco_SSE3      = gmx_cmplt_pr(rsq_SSE3,rc2_SSE);
#endif

#ifdef CHECK_EXCLS
#ifdef EXCL_FORCES
            /* Only remove the (sub-)diagonal to avoid double counting */
#if UNROLLJ == UNROLLI
            if (cj == ci_sh)
            {
                wco_SSE0  = gmx_and_pr(wco_SSE0,diag_SSE0);
                wco_SSE1  = gmx_and_pr(wco_SSE1,diag_SSE1);
                wco_SSE2  = gmx_and_pr(wco_SSE2,diag_SSE2);
                wco_SSE3  = gmx_and_pr(wco_SSE3,diag_SSE3);
            }
#else
#if UNROLLJ < UNROLLI
            if (cj == ci_sh*2)
            {
                wco_SSE0  = gmx_and_pr(wco_SSE0,diag0_SSE0);
                wco_SSE1  = gmx_and_pr(wco_SSE1,diag0_SSE1);
                wco_SSE2  = gmx_and_pr(wco_SSE2,diag0_SSE2);
                wco_SSE3  = gmx_and_pr(wco_SSE3,diag0_SSE3);
            }
            if (cj == ci_sh*2 + 1)
            { 
                wco_SSE0  = gmx_and_pr(wco_SSE0,diag1_SSE0);
                wco_SSE1  = gmx_and_pr(wco_SSE1,diag1_SSE1);
                wco_SSE2  = gmx_and_pr(wco_SSE2,diag1_SSE2);
                wco_SSE3  = gmx_and_pr(wco_SSE3,diag1_SSE3);
            }
#else
            if (cj*2 == ci_sh)
            {
                wco_SSE0  = gmx_and_pr(wco_SSE0,diag0_SSE0);
                wco_SSE1  = gmx_and_pr(wco_SSE1,diag0_SSE1);
                wco_SSE2  = gmx_and_pr(wco_SSE2,diag0_SSE2);
                wco_SSE3  = gmx_and_pr(wco_SSE3,diag0_SSE3);
            }
            else if (cj*2 + 1 == ci_sh)
            {
                wco_SSE0  = gmx_and_pr(wco_SSE0,diag1_SSE0);
                wco_SSE1  = gmx_and_pr(wco_SSE1,diag1_SSE1);
                wco_SSE2  = gmx_and_pr(wco_SSE2,diag1_SSE2);
                wco_SSE3  = gmx_and_pr(wco_SSE3,diag1_SSE3);
            }
#endif
#endif
#else
            /* Remove all excluded atom pairs from the list */
            wco_SSE0      = gmx_and_pr(wco_SSE0,int_SSE0);
            wco_SSE1      = gmx_and_pr(wco_SSE1,int_SSE1);
            wco_SSE2      = gmx_and_pr(wco_SSE2,int_SSE2);
            wco_SSE3      = gmx_and_pr(wco_SSE3,int_SSE3);
#endif
#endif

#ifdef COUNT_PAIRS
            {
                int i,j;
                real tmp[UNROLLJ];
                for(i=0; i<UNROLLI; i++)
                {
                    gmx_storeu_pr(tmp,i==0 ? wco_SSE0 : (i==1 ? wco_SSE1 : (i==2 ? wco_SSE2 : wco_SSE3)));
                    for(j=0; j<UNROLLJ; j++)
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
            rsq_SSE0      = gmx_add_pr(rsq_SSE0,gmx_andnot_pr(int_SSE0,avoid_sing_SSE));
            rsq_SSE1      = gmx_add_pr(rsq_SSE1,gmx_andnot_pr(int_SSE1,avoid_sing_SSE));
            rsq_SSE2      = gmx_add_pr(rsq_SSE2,gmx_andnot_pr(int_SSE2,avoid_sing_SSE));
            rsq_SSE3      = gmx_add_pr(rsq_SSE3,gmx_andnot_pr(int_SSE3,avoid_sing_SSE));
#endif

            /* Calculate 1/r */
#ifndef GMX_DOUBLE
            rinv_SSE0     = gmx_invsqrt_pr(rsq_SSE0);
            rinv_SSE1     = gmx_invsqrt_pr(rsq_SSE1);
            rinv_SSE2     = gmx_invsqrt_pr(rsq_SSE2);
            rinv_SSE3     = gmx_invsqrt_pr(rsq_SSE3);
#else
            GMX_MM_INVSQRT2_PD(rsq_SSE0,rsq_SSE1,rinv_SSE0,rinv_SSE1);
            GMX_MM_INVSQRT2_PD(rsq_SSE2,rsq_SSE3,rinv_SSE2,rinv_SSE3);
#endif

#ifdef CALC_COULOMB
            /* Load parameters for j atom */
            jq_SSE        = gmx_load_pr(q+aj);
            qq_SSE0       = gmx_mul_pr(iq_SSE0,jq_SSE);
            qq_SSE1       = gmx_mul_pr(iq_SSE1,jq_SSE);
            qq_SSE2       = gmx_mul_pr(iq_SSE2,jq_SSE);
            qq_SSE3       = gmx_mul_pr(iq_SSE3,jq_SSE);
#endif

#ifdef CALC_LJ
#ifdef LJ_COMB_LB
            hsig_j_SSE    = gmx_load_pr(ljc+aj2+0);
            seps_j_SSE    = gmx_load_pr(ljc+aj2+STRIDE);

            sig_SSE0      = gmx_add_pr(hsig_i_SSE0,hsig_j_SSE);
            sig_SSE1      = gmx_add_pr(hsig_i_SSE1,hsig_j_SSE);
            eps_SSE0      = gmx_mul_pr(seps_i_SSE0,seps_j_SSE);
            eps_SSE1      = gmx_mul_pr(seps_i_SSE1,seps_j_SSE);
#ifndef HALF_LJ
            sig_SSE2      = gmx_add_pr(hsig_i_SSE2,hsig_j_SSE);
            sig_SSE3      = gmx_add_pr(hsig_i_SSE3,hsig_j_SSE);
            eps_SSE2      = gmx_mul_pr(seps_i_SSE2,seps_j_SSE);
            eps_SSE3      = gmx_mul_pr(seps_i_SSE3,seps_j_SSE);
#endif
#else
#ifdef LJ_COMB_GEOM
            c6s_j_SSE     = gmx_load_pr(ljc+aj2+0);
            c12s_j_SSE    = gmx_load_pr(ljc+aj2+STRIDE);
            c6_SSE0       = gmx_mul_pr(c6s_SSE0 ,c6s_j_SSE );
            c6_SSE1       = gmx_mul_pr(c6s_SSE1 ,c6s_j_SSE );
#ifndef HALF_LJ
            c6_SSE2       = gmx_mul_pr(c6s_SSE2 ,c6s_j_SSE );
            c6_SSE3       = gmx_mul_pr(c6s_SSE3 ,c6s_j_SSE );
#endif
            c12_SSE0      = gmx_mul_pr(c12s_SSE0,c12s_j_SSE);
            c12_SSE1      = gmx_mul_pr(c12s_SSE1,c12s_j_SSE);
#ifndef HALF_LJ
            c12_SSE2      = gmx_mul_pr(c12s_SSE2,c12s_j_SSE);
            c12_SSE3      = gmx_mul_pr(c12s_SSE3,c12s_j_SSE);
#endif
#else
#ifndef FIX_LJ_C
            load_lj_pair_params(nbfp0,type,aj,c6_SSE0,c12_SSE0);
            load_lj_pair_params(nbfp1,type,aj,c6_SSE1,c12_SSE1);
#ifndef HALF_LJ
            load_lj_pair_params(nbfp2,type,aj,c6_SSE2,c12_SSE2);
            load_lj_pair_params(nbfp3,type,aj,c6_SSE3,c12_SSE3);
#endif
#endif /* FIX_LJ_C */
#endif /* LJ_COMB_GEOM */
#endif /* LJ_COMB_LB */
#endif /* CALC_LJ */

#ifndef CUTOFF_BLENDV
            rinv_SSE0     = gmx_and_pr(rinv_SSE0,wco_SSE0);
            rinv_SSE1     = gmx_and_pr(rinv_SSE1,wco_SSE1);
            rinv_SSE2     = gmx_and_pr(rinv_SSE2,wco_SSE2);
            rinv_SSE3     = gmx_and_pr(rinv_SSE3,wco_SSE3);
#else
            /* We only need to mask for the cut-off: blendv is faster */
            rinv_SSE0     = gmx_blendv_pr(rinv_SSE0,zero_SSE,gmx_sub_pr(rc2_SSE,rsq_SSE0));
            rinv_SSE1     = gmx_blendv_pr(rinv_SSE1,zero_SSE,gmx_sub_pr(rc2_SSE,rsq_SSE1));
            rinv_SSE2     = gmx_blendv_pr(rinv_SSE2,zero_SSE,gmx_sub_pr(rc2_SSE,rsq_SSE2));
            rinv_SSE3     = gmx_blendv_pr(rinv_SSE3,zero_SSE,gmx_sub_pr(rc2_SSE,rsq_SSE3));
#endif

            rinvsq_SSE0   = gmx_mul_pr(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1   = gmx_mul_pr(rinv_SSE1,rinv_SSE1);
            rinvsq_SSE2   = gmx_mul_pr(rinv_SSE2,rinv_SSE2);
            rinvsq_SSE3   = gmx_mul_pr(rinv_SSE3,rinv_SSE3);

#ifdef CALC_COULOMB

#ifdef EXCL_FORCES
            /* Only add 1/r for non-excluded atom pairs */
            rinv_ex_SSE0  = gmx_and_pr(rinv_SSE0,int_SSE0);
            rinv_ex_SSE1  = gmx_and_pr(rinv_SSE1,int_SSE1);
            rinv_ex_SSE2  = gmx_and_pr(rinv_SSE2,int_SSE2);
            rinv_ex_SSE3  = gmx_and_pr(rinv_SSE3,int_SSE3);
#else
            /* No exclusion forces, we always need 1/r */
#define     rinv_ex_SSE0    rinv_SSE0
#define     rinv_ex_SSE1    rinv_SSE1
#define     rinv_ex_SSE2    rinv_SSE2
#define     rinv_ex_SSE3    rinv_SSE3
#endif

#ifdef CALC_COUL_RF
            /* Electrostatic interactions */
            frcoul_SSE0   = gmx_mul_pr(qq_SSE0,gmx_add_pr(rinv_ex_SSE0,gmx_mul_pr(rsq_SSE0,mrc_3_SSE)));
            frcoul_SSE1   = gmx_mul_pr(qq_SSE1,gmx_add_pr(rinv_ex_SSE1,gmx_mul_pr(rsq_SSE1,mrc_3_SSE)));
            frcoul_SSE2   = gmx_mul_pr(qq_SSE2,gmx_add_pr(rinv_ex_SSE2,gmx_mul_pr(rsq_SSE2,mrc_3_SSE)));
            frcoul_SSE3   = gmx_mul_pr(qq_SSE3,gmx_add_pr(rinv_ex_SSE3,gmx_mul_pr(rsq_SSE3,mrc_3_SSE)));

#ifdef CALC_ENERGIES
            vcoul_SSE0    = gmx_mul_pr(qq_SSE0,gmx_add_pr(rinv_ex_SSE0,gmx_add_pr(gmx_mul_pr(rsq_SSE0,hrc_3_SSE),moh_rc_SSE)));
            vcoul_SSE1    = gmx_mul_pr(qq_SSE1,gmx_add_pr(rinv_ex_SSE1,gmx_add_pr(gmx_mul_pr(rsq_SSE1,hrc_3_SSE),moh_rc_SSE)));
            vcoul_SSE2    = gmx_mul_pr(qq_SSE2,gmx_add_pr(rinv_ex_SSE2,gmx_add_pr(gmx_mul_pr(rsq_SSE2,hrc_3_SSE),moh_rc_SSE)));
            vcoul_SSE3    = gmx_mul_pr(qq_SSE3,gmx_add_pr(rinv_ex_SSE3,gmx_add_pr(gmx_mul_pr(rsq_SSE3,hrc_3_SSE),moh_rc_SSE)));
#endif
#endif

#ifdef CALC_COUL_TAB
            /* Electrostatic interactions */
            r_SSE0        = gmx_mul_pr(rsq_SSE0,rinv_SSE0);
            r_SSE1        = gmx_mul_pr(rsq_SSE1,rinv_SSE1);
            r_SSE2        = gmx_mul_pr(rsq_SSE2,rinv_SSE2);
            r_SSE3        = gmx_mul_pr(rsq_SSE3,rinv_SSE3);
            /* Convert r to scaled table units */
            rs_SSE0       = gmx_mul_pr(r_SSE0,invtsp_SSE);
            rs_SSE1       = gmx_mul_pr(r_SSE1,invtsp_SSE);
            rs_SSE2       = gmx_mul_pr(r_SSE2,invtsp_SSE);
            rs_SSE3       = gmx_mul_pr(r_SSE3,invtsp_SSE);
            /* Truncate scaled r to an int */
            ti_SSE0       = gmx_cvttpr_epi32(rs_SSE0);
            ti_SSE1       = gmx_cvttpr_epi32(rs_SSE1);
            ti_SSE2       = gmx_cvttpr_epi32(rs_SSE2);
            ti_SSE3       = gmx_cvttpr_epi32(rs_SSE3);
#ifdef GMX_X86_SSE4_1
            /* SSE4.1 floor is faster than gmx_cvtepi32_ps int->float cast */
            rf_SSE0       = gmx_floor_pr(rs_SSE0);
            rf_SSE1       = gmx_floor_pr(rs_SSE1);
            rf_SSE2       = gmx_floor_pr(rs_SSE2);
            rf_SSE3       = gmx_floor_pr(rs_SSE3);
#else
            rf_SSE0       = gmx_cvtepi32_pr(ti_SSE0);
            rf_SSE1       = gmx_cvtepi32_pr(ti_SSE1);
            rf_SSE2       = gmx_cvtepi32_pr(ti_SSE2);
            rf_SSE3       = gmx_cvtepi32_pr(ti_SSE3);
#endif
            frac_SSE0     = gmx_sub_pr(rs_SSE0,rf_SSE0);
            frac_SSE1     = gmx_sub_pr(rs_SSE1,rf_SSE1);
            frac_SSE2     = gmx_sub_pr(rs_SSE2,rf_SSE2);
            frac_SSE3     = gmx_sub_pr(rs_SSE3,rf_SSE3);

            /* Load and interpolate table forces and possibly energies.
             * Force and energy can be combined in one table, stride 4: FDV0
             * or in two separate tables with stride 1: F and V
             * Currently single precision uses FDV0, double F and V.
             */
#ifndef CALC_ENERGIES
            load_table_f(tab_coul_F,ti_SSE0,ti0,ctab0_SSE0,ctab1_SSE0);
            load_table_f(tab_coul_F,ti_SSE1,ti1,ctab0_SSE1,ctab1_SSE1);
            load_table_f(tab_coul_F,ti_SSE2,ti2,ctab0_SSE2,ctab1_SSE2);
            load_table_f(tab_coul_F,ti_SSE3,ti3,ctab0_SSE3,ctab1_SSE3);
#else
#ifdef TAB_FDV0
            load_table_f_v(tab_coul_F,ti_SSE0,ti0,ctab0_SSE0,ctab1_SSE0,ctabv_SSE0);
            load_table_f_v(tab_coul_F,ti_SSE1,ti1,ctab0_SSE1,ctab1_SSE1,ctabv_SSE1);
            load_table_f_v(tab_coul_F,ti_SSE2,ti2,ctab0_SSE2,ctab1_SSE2,ctabv_SSE2);
            load_table_f_v(tab_coul_F,ti_SSE3,ti3,ctab0_SSE3,ctab1_SSE3,ctabv_SSE3);
#else
            load_table_f_v(tab_coul_F,tab_coul_V,ti_SSE0,ti0,ctab0_SSE0,ctab1_SSE0,ctabv_SSE0);
            load_table_f_v(tab_coul_F,tab_coul_V,ti_SSE1,ti1,ctab0_SSE1,ctab1_SSE1,ctabv_SSE1);
            load_table_f_v(tab_coul_F,tab_coul_V,ti_SSE2,ti2,ctab0_SSE2,ctab1_SSE2,ctabv_SSE2);
            load_table_f_v(tab_coul_F,tab_coul_V,ti_SSE3,ti3,ctab0_SSE3,ctab1_SSE3,ctabv_SSE3);
#endif
#endif
            fsub_SSE0     = gmx_add_pr(ctab0_SSE0,gmx_mul_pr(frac_SSE0,ctab1_SSE0));
            fsub_SSE1     = gmx_add_pr(ctab0_SSE1,gmx_mul_pr(frac_SSE1,ctab1_SSE1));
            fsub_SSE2     = gmx_add_pr(ctab0_SSE2,gmx_mul_pr(frac_SSE2,ctab1_SSE2));
            fsub_SSE3     = gmx_add_pr(ctab0_SSE3,gmx_mul_pr(frac_SSE3,ctab1_SSE3));
            frcoul_SSE0   = gmx_mul_pr(qq_SSE0,gmx_sub_pr(rinv_ex_SSE0,gmx_mul_pr(fsub_SSE0,r_SSE0)));
            frcoul_SSE1   = gmx_mul_pr(qq_SSE1,gmx_sub_pr(rinv_ex_SSE1,gmx_mul_pr(fsub_SSE1,r_SSE1)));
            frcoul_SSE2   = gmx_mul_pr(qq_SSE2,gmx_sub_pr(rinv_ex_SSE2,gmx_mul_pr(fsub_SSE2,r_SSE2)));
            frcoul_SSE3   = gmx_mul_pr(qq_SSE3,gmx_sub_pr(rinv_ex_SSE3,gmx_mul_pr(fsub_SSE3,r_SSE3)));

#ifdef CALC_ENERGIES
            vc_sub_SSE0   = gmx_add_pr(ctabv_SSE0,gmx_mul_pr(gmx_mul_pr(mhalfsp_SSE,frac_SSE0),gmx_add_pr(ctab0_SSE0,fsub_SSE0)));
            vc_sub_SSE1   = gmx_add_pr(ctabv_SSE1,gmx_mul_pr(gmx_mul_pr(mhalfsp_SSE,frac_SSE1),gmx_add_pr(ctab0_SSE1,fsub_SSE1)));
            vc_sub_SSE2   = gmx_add_pr(ctabv_SSE2,gmx_mul_pr(gmx_mul_pr(mhalfsp_SSE,frac_SSE2),gmx_add_pr(ctab0_SSE2,fsub_SSE2)));
            vc_sub_SSE3   = gmx_add_pr(ctabv_SSE3,gmx_mul_pr(gmx_mul_pr(mhalfsp_SSE,frac_SSE3),gmx_add_pr(ctab0_SSE3,fsub_SSE3)));

#ifndef NO_SHIFT_EWALD
            /* Add Ewald potential shift to vc_sub for convenience */
            vc_sub_SSE0   = gmx_add_pr(vc_sub_SSE0,gmx_and_pr(sh_ewald_SSE,wco_SSE0));
            vc_sub_SSE1   = gmx_add_pr(vc_sub_SSE1,gmx_and_pr(sh_ewald_SSE,wco_SSE1));
            vc_sub_SSE2   = gmx_add_pr(vc_sub_SSE2,gmx_and_pr(sh_ewald_SSE,wco_SSE2));
            vc_sub_SSE3   = gmx_add_pr(vc_sub_SSE3,gmx_and_pr(sh_ewald_SSE,wco_SSE3));
#endif
            
            vcoul_SSE0    = gmx_mul_pr(qq_SSE0,gmx_sub_pr(rinv_ex_SSE0,vc_sub_SSE0));
            vcoul_SSE1    = gmx_mul_pr(qq_SSE1,gmx_sub_pr(rinv_ex_SSE1,vc_sub_SSE1));
            vcoul_SSE2    = gmx_mul_pr(qq_SSE2,gmx_sub_pr(rinv_ex_SSE2,vc_sub_SSE2));
            vcoul_SSE3    = gmx_mul_pr(qq_SSE3,gmx_sub_pr(rinv_ex_SSE3,vc_sub_SSE3));

#endif
#endif

#ifdef CALC_ENERGIES
            /* Mask energy for cut-off and diagonal */
            vcoul_SSE0    = gmx_and_pr(vcoul_SSE0,wco_SSE0);
            vcoul_SSE1    = gmx_and_pr(vcoul_SSE1,wco_SSE1);
            vcoul_SSE2    = gmx_and_pr(vcoul_SSE2,wco_SSE2);
            vcoul_SSE3    = gmx_and_pr(vcoul_SSE3,wco_SSE3);
#endif

#endif /* CALC_COULOMB */

#ifdef CALC_LJ
            /* Lennard-Jones interaction */
#ifdef LJ_COMB_LB
            sir_SSE0      = gmx_mul_pr(sig_SSE0,rinv_SSE0);
            sir_SSE1      = gmx_mul_pr(sig_SSE1,rinv_SSE1);
#ifndef HALF_LJ
            sir_SSE2      = gmx_mul_pr(sig_SSE2,rinv_SSE2);
            sir_SSE3      = gmx_mul_pr(sig_SSE3,rinv_SSE3);
#endif
            sir2_SSE0     = gmx_mul_pr(sir_SSE0,sir_SSE0);
            sir2_SSE1     = gmx_mul_pr(sir_SSE1,sir_SSE1);
#ifndef HALF_LJ
            sir2_SSE2     = gmx_mul_pr(sir_SSE2,sir_SSE2);
            sir2_SSE3     = gmx_mul_pr(sir_SSE3,sir_SSE3);
#endif
            sir6_SSE0     = gmx_mul_pr(sir2_SSE0,gmx_mul_pr(sir2_SSE0,sir2_SSE0));
            sir6_SSE1     = gmx_mul_pr(sir2_SSE1,gmx_mul_pr(sir2_SSE1,sir2_SSE1));
#ifdef EXCL_FORCES
            sir6_SSE0     = gmx_and_pr(sir6_SSE0,int_SSE0);
            sir6_SSE1     = gmx_and_pr(sir6_SSE1,int_SSE1);
#endif
#ifndef HALF_LJ
            sir6_SSE2     = gmx_mul_pr(sir2_SSE2,gmx_mul_pr(sir2_SSE2,sir2_SSE2));
            sir6_SSE3     = gmx_mul_pr(sir2_SSE3,gmx_mul_pr(sir2_SSE3,sir2_SSE3));
#ifdef EXCL_FORCES
            sir6_SSE2     = gmx_and_pr(sir6_SSE2,int_SSE2);
            sir6_SSE3     = gmx_and_pr(sir6_SSE3,int_SSE3);
#endif
#endif
            FrLJ6_SSE0    = gmx_mul_pr(eps_SSE0,sir6_SSE0);
            FrLJ6_SSE1    = gmx_mul_pr(eps_SSE1,sir6_SSE1);
#ifndef HALF_LJ
            FrLJ6_SSE2    = gmx_mul_pr(eps_SSE2,sir6_SSE2);
            FrLJ6_SSE3    = gmx_mul_pr(eps_SSE3,sir6_SSE3);
#endif
            FrLJ12_SSE0   = gmx_mul_pr(FrLJ6_SSE0,sir6_SSE0);
            FrLJ12_SSE1   = gmx_mul_pr(FrLJ6_SSE1,sir6_SSE1);
#ifndef HALF_LJ
            FrLJ12_SSE2   = gmx_mul_pr(FrLJ6_SSE2,sir6_SSE2);
            FrLJ12_SSE3   = gmx_mul_pr(FrLJ6_SSE3,sir6_SSE3);
#endif
#if defined CALC_ENERGIES && !defined NO_LJ_SHIFT
            /* We need C6 and C12 to calculate the LJ potential shift */
            sig6_SSE0     = gmx_mul_pr(sig_SSE0,sig_SSE0);
            sig6_SSE1     = gmx_mul_pr(sig_SSE1,sig_SSE1);
#ifndef HALF_LJ
            sig6_SSE2     = gmx_mul_pr(sig_SSE2,sig_SSE2);
            sig6_SSE3     = gmx_mul_pr(sig_SSE3,sig_SSE3);
#endif
            sig6_SSE0     = gmx_mul_pr(sig6_SSE0,gmx_mul_pr(sig6_SSE0,sig6_SSE0));
            sig6_SSE1     = gmx_mul_pr(sig6_SSE1,gmx_mul_pr(sig6_SSE1,sig6_SSE1));
#ifndef HALF_LJ
            sig6_SSE2     = gmx_mul_pr(sig6_SSE2,gmx_mul_pr(sig6_SSE2,sig6_SSE2));
            sig6_SSE3     = gmx_mul_pr(sig6_SSE3,gmx_mul_pr(sig6_SSE3,sig6_SSE3));
#endif
            c6_SSE0       = gmx_mul_pr(eps_SSE0,sig6_SSE0);
            c6_SSE1       = gmx_mul_pr(eps_SSE1,sig6_SSE1);
#ifndef HALF_LJ
            c6_SSE2       = gmx_mul_pr(eps_SSE2,sig6_SSE2);
            c6_SSE3       = gmx_mul_pr(eps_SSE3,sig6_SSE3);
#endif
            c12_SSE0      = gmx_mul_pr(c6_SSE0,sig6_SSE0);
            c12_SSE1      = gmx_mul_pr(c6_SSE1,sig6_SSE1);
#ifndef HALF_LJ
            c12_SSE2      = gmx_mul_pr(c6_SSE2,sig6_SSE2);
            c12_SSE3      = gmx_mul_pr(c6_SSE3,sig6_SSE3);
#endif
#endif
#else /* ifdef LJ_COMB_LB */
            rinvsix_SSE0  = gmx_mul_pr(rinvsq_SSE0,gmx_mul_pr(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1  = gmx_mul_pr(rinvsq_SSE1,gmx_mul_pr(rinvsq_SSE1,rinvsq_SSE1));
#ifdef EXCL_FORCES
            rinvsix_SSE0  = gmx_and_pr(rinvsix_SSE0,int_SSE0);
            rinvsix_SSE1  = gmx_and_pr(rinvsix_SSE1,int_SSE1);
#endif
#ifndef HALF_LJ
            rinvsix_SSE2  = gmx_mul_pr(rinvsq_SSE2,gmx_mul_pr(rinvsq_SSE2,rinvsq_SSE2));
            rinvsix_SSE3  = gmx_mul_pr(rinvsq_SSE3,gmx_mul_pr(rinvsq_SSE3,rinvsq_SSE3));
#ifdef EXCL_FORCES
            rinvsix_SSE2  = gmx_and_pr(rinvsix_SSE2,int_SSE2);
            rinvsix_SSE3  = gmx_and_pr(rinvsix_SSE3,int_SSE3);
#endif
#endif
            FrLJ6_SSE0    = gmx_mul_pr(c6_SSE0,rinvsix_SSE0);
            FrLJ6_SSE1    = gmx_mul_pr(c6_SSE1,rinvsix_SSE1);
#ifndef HALF_LJ
            FrLJ6_SSE2    = gmx_mul_pr(c6_SSE2,rinvsix_SSE2);
            FrLJ6_SSE3    = gmx_mul_pr(c6_SSE3,rinvsix_SSE3);
#endif
            FrLJ12_SSE0   = gmx_mul_pr(c12_SSE0,gmx_mul_pr(rinvsix_SSE0,rinvsix_SSE0));
            FrLJ12_SSE1   = gmx_mul_pr(c12_SSE1,gmx_mul_pr(rinvsix_SSE1,rinvsix_SSE1));
#ifndef HALF_LJ
            FrLJ12_SSE2   = gmx_mul_pr(c12_SSE2,gmx_mul_pr(rinvsix_SSE2,rinvsix_SSE2));
            FrLJ12_SSE3   = gmx_mul_pr(c12_SSE3,gmx_mul_pr(rinvsix_SSE3,rinvsix_SSE3));
#endif
#endif /* LJ_COMB_LB */
#endif /* CALC_LJ */
            
#ifdef CALC_ENERGIES
#ifdef ENERGY_GROUPS
            /* Extract the group pair index per j pair */
#if UNROLLJ == 2
            egps_j        = nbat->energrp[cj>>1];
            egp_jj[0]     = ((egps_j >> ((cj & 1)*egps_jshift)) & egps_jmask)*egps_jstride;
#else
            egps_j        = nbat->energrp[cj];
            for(jj=0; jj<(UNROLLJ>>1); jj++)
            {
                egp_jj[jj]  = ((egps_j >> (jj*egps_jshift)) & egps_jmask)*egps_jstride;
            }
#endif
#endif

#ifdef CALC_COULOMB
#ifndef ENERGY_GROUPS
            vctotSSE      = gmx_add_pr(vctotSSE, gmx_sum4_pr(vcoul_SSE0,vcoul_SSE1,vcoul_SSE2,vcoul_SSE3));
#else
            add_ener_grp(vcoul_SSE0,vctp[0],egp_jj);
            add_ener_grp(vcoul_SSE1,vctp[1],egp_jj);
            add_ener_grp(vcoul_SSE2,vctp[2],egp_jj);
            add_ener_grp(vcoul_SSE3,vctp[3],egp_jj);
#endif
#endif

#ifdef CALC_LJ
            /* Calculate the LJ energies */
#ifdef NO_LJ_SHIFT
            VLJ6_SSE0     = gmx_mul_pr(sixthSSE,FrLJ6_SSE0);
            VLJ6_SSE1     = gmx_mul_pr(sixthSSE,FrLJ6_SSE1);
#ifndef HALF_LJ
            VLJ6_SSE2     = gmx_mul_pr(sixthSSE,FrLJ6_SSE2);
            VLJ6_SSE3     = gmx_mul_pr(sixthSSE,FrLJ6_SSE3);
#endif
            VLJ12_SSE0    = gmx_mul_pr(twelvethSSE,FrLJ12_SSE0);
            VLJ12_SSE1    = gmx_mul_pr(twelvethSSE,FrLJ12_SSE1);
#ifndef HALF_LJ
            VLJ12_SSE2    = gmx_mul_pr(twelvethSSE,FrLJ12_SSE2);
            VLJ12_SSE3    = gmx_mul_pr(twelvethSSE,FrLJ12_SSE3);
#endif
#else /* ifdef NO_LJ_SHIFT */
            VLJ6_SSE0     = gmx_mul_pr(sixthSSE,gmx_sub_pr(FrLJ6_SSE0,gmx_mul_pr(c6_SSE0,sh_invrc6_SSE)));
            VLJ6_SSE1     = gmx_mul_pr(sixthSSE,gmx_sub_pr(FrLJ6_SSE1,gmx_mul_pr(c6_SSE1,sh_invrc6_SSE)));
#ifndef HALF_LJ
            VLJ6_SSE2     = gmx_mul_pr(sixthSSE,gmx_sub_pr(FrLJ6_SSE2,gmx_mul_pr(c6_SSE2,sh_invrc6_SSE)));
            VLJ6_SSE3     = gmx_mul_pr(sixthSSE,gmx_sub_pr(FrLJ6_SSE3,gmx_mul_pr(c6_SSE3,sh_invrc6_SSE)));
#endif
            VLJ12_SSE0    = gmx_mul_pr(twelvethSSE,gmx_sub_pr(FrLJ12_SSE0,gmx_mul_pr(c12_SSE0,sh_invrc12_SSE)));
            VLJ12_SSE1    = gmx_mul_pr(twelvethSSE,gmx_sub_pr(FrLJ12_SSE1,gmx_mul_pr(c12_SSE1,sh_invrc12_SSE)));
#ifndef HALF_LJ
            VLJ12_SSE2    = gmx_mul_pr(twelvethSSE,gmx_sub_pr(FrLJ12_SSE2,gmx_mul_pr(c12_SSE2,sh_invrc12_SSE)));
            VLJ12_SSE3    = gmx_mul_pr(twelvethSSE,gmx_sub_pr(FrLJ12_SSE3,gmx_mul_pr(c12_SSE3,sh_invrc12_SSE)));
#endif
#endif /* ifdef NO_LJ_SHIFT */

            VLJ_SSE0      = gmx_sub_pr(VLJ12_SSE0,VLJ6_SSE0);
            VLJ_SSE1      = gmx_sub_pr(VLJ12_SSE1,VLJ6_SSE1);
#ifndef HALF_LJ
            VLJ_SSE2      = gmx_sub_pr(VLJ12_SSE2,VLJ6_SSE2);
            VLJ_SSE3      = gmx_sub_pr(VLJ12_SSE3,VLJ6_SSE3);
#endif
#ifndef NO_LJ_SHIFT
            /* The potential shift should be removed non-interacting pairs */
            VLJ_SSE0      = gmx_and_pr(VLJ_SSE0,wco_SSE0);
            VLJ_SSE1      = gmx_and_pr(VLJ_SSE1,wco_SSE1);
#ifndef HALF_LJ
            VLJ_SSE2      = gmx_and_pr(VLJ_SSE2,wco_SSE2);
            VLJ_SSE3      = gmx_and_pr(VLJ_SSE3,wco_SSE3);
#endif
#endif
#ifndef ENERGY_GROUPS
            VvdwtotSSE    = gmx_add_pr(VvdwtotSSE,
#ifndef HALF_LJ
                                       gmx_sum4_pr(VLJ_SSE0,VLJ_SSE1,VLJ_SSE2,VLJ_SSE3)
#else
                                       gmx_add_pr(VLJ_SSE0,VLJ_SSE1)
#endif
                                      );
#else
            add_ener_grp(VLJ_SSE0,vvdwtp[0],egp_jj);
            add_ener_grp(VLJ_SSE1,vvdwtp[1],egp_jj);
#ifndef HALF_LJ
            add_ener_grp(VLJ_SSE2,vvdwtp[2],egp_jj);
            add_ener_grp(VLJ_SSE3,vvdwtp[3],egp_jj);
#endif
#endif
#endif /* CALC_LJ */
#endif /* CALC_ENERGIES */

#ifdef CALC_LJ
            fscal_SSE0    = gmx_mul_pr(rinvsq_SSE0,
#ifdef CALC_COULOMB
                                                   gmx_add_pr(frcoul_SSE0,
#else
                                                   (
#endif
                                                    gmx_sub_pr(FrLJ12_SSE0,FrLJ6_SSE0)));
            fscal_SSE1    = gmx_mul_pr(rinvsq_SSE1,
#ifdef CALC_COULOMB
                                                   gmx_add_pr(frcoul_SSE1,
#else
                                                   (
#endif
                                                    gmx_sub_pr(FrLJ12_SSE1,FrLJ6_SSE1)));
#else
            fscal_SSE0    = gmx_mul_pr(rinvsq_SSE0,frcoul_SSE0);
            fscal_SSE1    = gmx_mul_pr(rinvsq_SSE1,frcoul_SSE1);
#endif /* CALC_LJ */
#if defined CALC_LJ && !defined HALF_LJ
            fscal_SSE2    = gmx_mul_pr(rinvsq_SSE2,
#ifdef CALC_COULOMB
                                                   gmx_add_pr(frcoul_SSE2,
#else
                                                   (
#endif
                                                    gmx_sub_pr(FrLJ12_SSE2,FrLJ6_SSE2)));
            fscal_SSE3    = gmx_mul_pr(rinvsq_SSE3,
#ifdef CALC_COULOMB
                                                   gmx_add_pr(frcoul_SSE3,
#else
                                                   (
#endif
                                                    gmx_sub_pr(FrLJ12_SSE3,FrLJ6_SSE3)));
#else
            /* Atom 2 and 3 don't have LJ, so only add Coulomb forces */
            fscal_SSE2    = gmx_mul_pr(rinvsq_SSE2,frcoul_SSE2);
            fscal_SSE3    = gmx_mul_pr(rinvsq_SSE3,frcoul_SSE3);
#endif
            
            /* Calculate temporary vectorial force */
            tx_SSE0       = gmx_mul_pr(fscal_SSE0,dx_SSE0);
            tx_SSE1       = gmx_mul_pr(fscal_SSE1,dx_SSE1);
            tx_SSE2       = gmx_mul_pr(fscal_SSE2,dx_SSE2);
            tx_SSE3       = gmx_mul_pr(fscal_SSE3,dx_SSE3);
            ty_SSE0       = gmx_mul_pr(fscal_SSE0,dy_SSE0);
            ty_SSE1       = gmx_mul_pr(fscal_SSE1,dy_SSE1);
            ty_SSE2       = gmx_mul_pr(fscal_SSE2,dy_SSE2);
            ty_SSE3       = gmx_mul_pr(fscal_SSE3,dy_SSE3);
            tz_SSE0       = gmx_mul_pr(fscal_SSE0,dz_SSE0);
            tz_SSE1       = gmx_mul_pr(fscal_SSE1,dz_SSE1);
            tz_SSE2       = gmx_mul_pr(fscal_SSE2,dz_SSE2);
            tz_SSE3       = gmx_mul_pr(fscal_SSE3,dz_SSE3);
            
            /* Increment i atom force */
            fix_SSE0      = gmx_add_pr(fix_SSE0,tx_SSE0);
            fix_SSE1      = gmx_add_pr(fix_SSE1,tx_SSE1);
            fix_SSE2      = gmx_add_pr(fix_SSE2,tx_SSE2);
            fix_SSE3      = gmx_add_pr(fix_SSE3,tx_SSE3);
            fiy_SSE0      = gmx_add_pr(fiy_SSE0,ty_SSE0);
            fiy_SSE1      = gmx_add_pr(fiy_SSE1,ty_SSE1);
            fiy_SSE2      = gmx_add_pr(fiy_SSE2,ty_SSE2);
            fiy_SSE3      = gmx_add_pr(fiy_SSE3,ty_SSE3);
            fiz_SSE0      = gmx_add_pr(fiz_SSE0,tz_SSE0);
            fiz_SSE1      = gmx_add_pr(fiz_SSE1,tz_SSE1);
            fiz_SSE2      = gmx_add_pr(fiz_SSE2,tz_SSE2);
            fiz_SSE3      = gmx_add_pr(fiz_SSE3,tz_SSE3);
            
            /* Decrement j atom force */
            gmx_store_pr(f+ajx,
                         gmx_sub_pr( gmx_load_pr(f+ajx), gmx_sum4_pr(tx_SSE0,tx_SSE1,tx_SSE2,tx_SSE3) ));
            gmx_store_pr(f+ajy,
                         gmx_sub_pr( gmx_load_pr(f+ajy), gmx_sum4_pr(ty_SSE0,ty_SSE1,ty_SSE2,ty_SSE3) ));
            gmx_store_pr(f+ajz,
                         gmx_sub_pr( gmx_load_pr(f+ajz), gmx_sum4_pr(tz_SSE0,tz_SSE1,tz_SSE2,tz_SSE3) ));
        }

#undef  rinv_ex_SSE0
#undef  rinv_ex_SSE1
#undef  rinv_ex_SSE2
#undef  rinv_ex_SSE3

#undef  CUTOFF_BLENDV

#undef  EXCL_FORCES
