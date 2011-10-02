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

/* This is the innermost loop contents for the ns cell vs cell
 * SSE2 single precision kernels.
 */

        {
            int        sj,ssj,ssjx,ssjy,ssjz;

            __m128     mask_int;
            __m128     int_SSE0;
            __m128     int_SSE1;
            __m128     int_SSE2;
            __m128     int_SSE3;

            __m128     jxSSE,jySSE,jzSSE;
            __m128     dx_SSE0,dy_SSE0,dz_SSE0;
            __m128     dx_SSE1,dy_SSE1,dz_SSE1;
            __m128     dx_SSE2,dy_SSE2,dz_SSE2;
            __m128     dx_SSE3,dy_SSE3,dz_SSE3;
            __m128     tx_SSE0,ty_SSE0,tz_SSE0;
            __m128     tx_SSE1,ty_SSE1,tz_SSE1;
            __m128     tx_SSE2,ty_SSE2,tz_SSE2;
            __m128     tx_SSE3,ty_SSE3,tz_SSE3;
            __m128     rsq_SSE0,rinv_SSE0,rinvsq_SSE0;
            __m128     rsq_SSE1,rinv_SSE1,rinvsq_SSE1;
            __m128     rsq_SSE2,rinv_SSE2,rinvsq_SSE2;
            __m128     rsq_SSE3,rinv_SSE3,rinvsq_SSE3;
            __m128     wco_SSE0;
            __m128     wco_SSE1;
            __m128     wco_SSE2;
            __m128     wco_SSE3;
            __m128     wco_any_SSE01,wco_any_SSE23,wco_any_SSE;
#ifdef CALC_COULOMB
            __m128     jq_SSE;
            __m128     qq_SSE0;
            __m128     qq_SSE1;
            __m128     qq_SSE2;
            __m128     qq_SSE3;
#ifndef CALC_COUL_RF
            __m128     fcoul_SSE0;
            __m128     fcoul_SSE1;
            __m128     fcoul_SSE2;
            __m128     fcoul_SSE3;
#endif
            __m128     frcoul_SSE0;
            __m128     frcoul_SSE1;
            __m128     frcoul_SSE2;
            __m128     frcoul_SSE3;
#ifndef CALC_COUL_RF
            __m128     r_SSE0,rs_SSE0,rf_SSE0,frac_SSE0;
            __m128     r_SSE1,rs_SSE1,rf_SSE1,frac_SSE1;
            __m128     r_SSE2,rs_SSE2,rf_SSE2,frac_SSE2;
            __m128     r_SSE3,rs_SSE3,rf_SSE3,frac_SSE3;
            __m128i    ti_SSE0;
            __m128i    ti_SSE1;
            __m128i    ti_SSE2;
            __m128i    ti_SSE3;
            int        e;
            __m128     ctab_SSE0[4],ctab0_SSE0,ctab1_SSE0;
            __m128     ctab_SSE1[4],ctab0_SSE1,ctab1_SSE1;
            __m128     ctab_SSE2[4],ctab0_SSE2,ctab1_SSE2;
            __m128     ctab_SSE3[4],ctab0_SSE3,ctab1_SSE3;
#ifdef CALC_ENERGIES
            __m128     ctabv_SSE0;
            __m128     ctabv_SSE1;
            __m128     ctabv_SSE2;
            __m128     ctabv_SSE3;
#endif
#endif
#ifdef CALC_ENERGIES
            __m128     vcoul_SSE0;
            __m128     vcoul_SSE1;
            __m128     vcoul_SSE2;
            __m128     vcoul_SSE3;
#ifdef ENERGY_GROUPS
            __m128     vct_SSE0;
            __m128     vct_SSE1;
            __m128     vct_SSE2;
            __m128     vct_SSE3;
#endif
#endif
#endif
            __m128     fscal_SSE0;
            __m128     fscal_SSE1;
            __m128     fscal_SSE2;
            __m128     fscal_SSE3;

#ifdef LJ_COMB_LB
            __m128     hsig_j_SSE,seps_j_SSE;
#endif

#ifdef LJ_COMB_GEOM
            __m128     c6s_j_SSE,c12s_j_SSE;
#endif

#if !defined LJ_COMB_GEOM && !defined LJ_COMB_LB
            int         p;
            __m128      clj_SSE0[4];
            __m128      clj_SSE1[4];
#ifndef HALF_LJ
            __m128      clj_SSE2[4];
            __m128      clj_SSE3[4];
#endif
#endif

#ifndef FIX_LJ_C
            __m128     c6_SSE0,c12_SSE0;
            __m128     c6_SSE1,c12_SSE1;
#ifndef HALF_LJ
            __m128     c6_SSE2,c12_SSE2;
            __m128     c6_SSE3,c12_SSE3;
#endif
#endif

#ifdef LJ_COMB_LB
            __m128     sir_SSE0,sir2_SSE0,sir6_SSE0;
            __m128     sir_SSE1,sir2_SSE1,sir6_SSE1;
#ifndef HALF_LJ
            __m128     sir_SSE2,sir2_SSE2,sir6_SSE2;
            __m128     sir_SSE3,sir2_SSE3,sir6_SSE3;
#endif
#else
            __m128     rinvsix_SSE0;
            __m128     rinvsix_SSE1;
#ifndef HALF_LJ
            __m128     rinvsix_SSE2;
            __m128     rinvsix_SSE3;
#endif
#endif

            __m128     Vvdw6_SSE0,Vvdw12_SSE0;
            __m128     Vvdw6_SSE1,Vvdw12_SSE1;
#ifndef HALF_LJ
            __m128     Vvdw6_SSE2,Vvdw12_SSE2;
            __m128     Vvdw6_SSE3,Vvdw12_SSE3;
#endif
#ifdef CALC_ENERGIES
            __m128     Vvdw_SSE0,Vvdw_SSE1;
#ifndef HALF_LJ
            __m128     Vvdw_SSE2,Vvdw_SSE3;
#endif
#ifdef ENERGY_GROUPS
            __m128     vvdwt_SSE0;
            __m128     vvdwt_SSE1;
#ifndef HALF_LJ
            __m128     vvdwt_SSE2;
            __m128     vvdwt_SSE3;
#endif
#endif
#endif

            sj               = cj[sjind].c;

            ssj              = sj*SIMD_WIDTH;
            ssjx             = ssj*DIM;
            ssjy             = ssj*DIM + SIMD_WIDTH;
            ssjz             = ssj*DIM + 2*SIMD_WIDTH;

#ifdef CHECK_EXCLS
            /* Load integer interaction mask as float to avoid sse casts */
            mask_int         = _mm_load1_ps((float *)&cj[sjind].excl);
            int_SSE0         = _mm_cmpneq_ps(_mm_and_ps(mask_int,mask0),zero_SSE);
            int_SSE1         = _mm_cmpneq_ps(_mm_and_ps(mask_int,mask1),zero_SSE);
            int_SSE2         = _mm_cmpneq_ps(_mm_and_ps(mask_int,mask2),zero_SSE);
            int_SSE3         = _mm_cmpneq_ps(_mm_and_ps(mask_int,mask3),zero_SSE);
#endif
            /* load j atom coordinates */
            jxSSE            = _mm_load_ps(x+ssjx);
            jySSE            = _mm_load_ps(x+ssjy);
            jzSSE            = _mm_load_ps(x+ssjz);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0,jxSSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0,jySSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0,jzSSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1,jxSSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1,jySSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1,jzSSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2,jxSSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2,jySSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2,jzSSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3,jxSSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3,jySSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3,jzSSE);
            
            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0,dy_SSE0,dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1,dy_SSE1,dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2,dy_SSE2,dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3,dy_SSE3,dz_SSE3);

            wco_SSE0           = _mm_cmplt_ps(rsq_SSE0,rc2_SSE);
            wco_SSE1           = _mm_cmplt_ps(rsq_SSE1,rc2_SSE);
            wco_SSE2           = _mm_cmplt_ps(rsq_SSE2,rc2_SSE);
            wco_SSE3           = _mm_cmplt_ps(rsq_SSE3,rc2_SSE);

#ifdef CHECK_EXCLS
            wco_SSE0           = _mm_and_ps(wco_SSE0,int_SSE0);
            wco_SSE1           = _mm_and_ps(wco_SSE1,int_SSE1);
            wco_SSE2           = _mm_and_ps(wco_SSE2,int_SSE2);
            wco_SSE3           = _mm_and_ps(wco_SSE3,int_SSE3);
#endif

#ifdef COUNT_PAIRS
            {
                int i,j;
                float tmp[4];
                for(i=0; i<SIMD_WIDTH; i++)
                {
                    _mm_storeu_ps(tmp,i==0 ? wco_SSE0 : (i==1 ? wco_SSE1 : (i==2 ? wco_SSE2 : wco_SSE3)));
                    for(j=0; j<SIMD_WIDTH; j++)
                    {
                        if (!(tmp[j] == 0))
                        {
                            npair++;
                        }
                    }
                }
            }
#endif

#ifdef SSE_CONDITIONAL
            wco_any_SSE01      = _mm_or_ps(wco_SSE0,wco_SSE1);
            wco_any_SSE23      = _mm_or_ps(wco_SSE2,wco_SSE3);
            wco_any_SSE        = _mm_or_ps(wco_any_SSE01,wco_any_SSE23);

            _mm_store_ps(wco_any_align,wco_any_SSE);
            
            if (wco_any_align[0] == 0 &&
                wco_any_align[1] == 0 && 
                wco_any_align[2] == 0 &&
                wco_any_align[3] == 0)
            {
                continue;
            }
#endif

            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);

#ifdef CALC_COULOMB
            /* Load parameters for j atom */
            jq_SSE             = _mm_load_ps(q+ssj);
            qq_SSE0            = _mm_mul_ps(iq_SSE0,jq_SSE);
            qq_SSE1            = _mm_mul_ps(iq_SSE1,jq_SSE);
            qq_SSE2            = _mm_mul_ps(iq_SSE2,jq_SSE);
            qq_SSE3            = _mm_mul_ps(iq_SSE3,jq_SSE);
#endif

#ifdef LJ_COMB_LB
            hsig_j_SSE         = _mm_load_ps(ljc+ssj*2+0);
            seps_j_SSE         = _mm_load_ps(ljc+ssj*2+4);
#else
#ifdef LJ_COMB_GEOM
            c6s_j_SSE          = _mm_load_ps(ljc+ssj*2+0);
            c12s_j_SSE         = _mm_load_ps(ljc+ssj*2+4);
            c6_SSE0            = _mm_mul_ps(c6s_SSE0 ,c6s_j_SSE );
            c6_SSE1            = _mm_mul_ps(c6s_SSE1 ,c6s_j_SSE );
#ifndef HALF_LJ
            c6_SSE2            = _mm_mul_ps(c6s_SSE2 ,c6s_j_SSE );
            c6_SSE3            = _mm_mul_ps(c6s_SSE3 ,c6s_j_SSE );
#endif
            c12_SSE0           = _mm_mul_ps(c12s_SSE0,c12s_j_SSE);
            c12_SSE1           = _mm_mul_ps(c12s_SSE1,c12s_j_SSE);
#ifndef HALF_LJ
            c12_SSE2           = _mm_mul_ps(c12s_SSE2,c12s_j_SSE);
            c12_SSE3           = _mm_mul_ps(c12s_SSE3,c12s_j_SSE);
#endif
#else
#ifndef FIX_LJ_C
            for(p=0; p<4; p++)
            {
                clj_SSE0[p]   = _mm_load_ps(nbfp0+type[ssj+p]*4);
            }
            for(p=0; p<4; p++)
            {
                clj_SSE1[p]   = _mm_load_ps(nbfp1+type[ssj+p]*4);
            }
#ifndef HALF_LJ
            for(p=0; p<4; p++)
            {
                clj_SSE2[p]   = _mm_load_ps(nbfp2+type[ssj+p]*4);
            }
            for(p=0; p<4; p++)
            {
                clj_SSE3[p]   = _mm_load_ps(nbfp3+type[ssj+p]*4);
            }
#endif
            GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SSE0[0],clj_SSE0[1],clj_SSE0[2],clj_SSE0[3],c6_SSE0,c12_SSE0);
            GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SSE1[0],clj_SSE1[1],clj_SSE1[2],clj_SSE1[3],c6_SSE1,c12_SSE1);
#ifndef HALF_LJ
            GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SSE2[0],clj_SSE2[1],clj_SSE2[2],clj_SSE2[3],c6_SSE2,c12_SSE2);
            GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SSE3[0],clj_SSE3[1],clj_SSE3[2],clj_SSE3[3],c6_SSE3,c12_SSE3);
#endif
#endif /* FIX_LJ_C */
#endif /* LJ_COMB_GEOM */
#endif /* LJ_COMB_LB */

            
            rinv_SSE0          = _mm_and_ps(rinv_SSE0,wco_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1,wco_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2,wco_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3,wco_SSE3);

            rinvsq_SSE0        = _mm_mul_ps(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1        = _mm_mul_ps(rinv_SSE1,rinv_SSE1);
#if !defined HALF_LJ || defined CALC_COUL_RF
            rinvsq_SSE2        = _mm_mul_ps(rinv_SSE2,rinv_SSE2);
            rinvsq_SSE3        = _mm_mul_ps(rinv_SSE3,rinv_SSE3);
#endif

#ifdef CALC_COULOMB
#ifdef CALC_COUL_RF
            /* Coulomb interaction */
            frcoul_SSE0        = _mm_mul_ps(qq_SSE0,_mm_add_ps(rinv_SSE0,_mm_mul_ps(rsq_SSE0,mrc_3_SSE)));
            frcoul_SSE1        = _mm_mul_ps(qq_SSE1,_mm_add_ps(rinv_SSE1,_mm_mul_ps(rsq_SSE1,mrc_3_SSE)));
            frcoul_SSE2        = _mm_mul_ps(qq_SSE2,_mm_add_ps(rinv_SSE2,_mm_mul_ps(rsq_SSE2,mrc_3_SSE)));
            frcoul_SSE3        = _mm_mul_ps(qq_SSE3,_mm_add_ps(rinv_SSE3,_mm_mul_ps(rsq_SSE3,mrc_3_SSE)));

#ifdef CALC_ENERGIES
            vcoul_SSE0         = _mm_mul_ps(qq_SSE0,_mm_add_ps(rinv_SSE0,_mm_add_ps(_mm_mul_ps(rsq_SSE0,hrc_3_SSE),moh_rc_SSE)));
            vcoul_SSE1         = _mm_mul_ps(qq_SSE1,_mm_add_ps(rinv_SSE1,_mm_add_ps(_mm_mul_ps(rsq_SSE1,hrc_3_SSE),moh_rc_SSE)));
            vcoul_SSE2         = _mm_mul_ps(qq_SSE2,_mm_add_ps(rinv_SSE2,_mm_add_ps(_mm_mul_ps(rsq_SSE2,hrc_3_SSE),moh_rc_SSE)));
            vcoul_SSE3         = _mm_mul_ps(qq_SSE3,_mm_add_ps(rinv_SSE3,_mm_add_ps(_mm_mul_ps(rsq_SSE3,hrc_3_SSE),moh_rc_SSE)));
            vcoul_SSE0         = _mm_and_ps(vcoul_SSE0,wco_SSE0);
            vcoul_SSE1         = _mm_and_ps(vcoul_SSE1,wco_SSE1);
            vcoul_SSE2         = _mm_and_ps(vcoul_SSE2,wco_SSE2);
            vcoul_SSE3         = _mm_and_ps(vcoul_SSE3,wco_SSE3);
#endif
#else
            r_SSE0             = _mm_mul_ps(rsq_SSE0,rinv_SSE0);
            r_SSE1             = _mm_mul_ps(rsq_SSE1,rinv_SSE1);
            r_SSE2             = _mm_mul_ps(rsq_SSE2,rinv_SSE2);
            r_SSE3             = _mm_mul_ps(rsq_SSE3,rinv_SSE3);
            /* Convert r to scaled table units */
            rs_SSE0            = _mm_mul_ps(r_SSE0,invtsp_SSE);
            rs_SSE1            = _mm_mul_ps(r_SSE1,invtsp_SSE);
            rs_SSE2            = _mm_mul_ps(r_SSE2,invtsp_SSE);
            rs_SSE3            = _mm_mul_ps(r_SSE3,invtsp_SSE);
            /* Truncate scaled r to an int */
            ti_SSE0            = _mm_cvttps_epi32(rs_SSE0);
            ti_SSE1            = _mm_cvttps_epi32(rs_SSE1);
            ti_SSE2            = _mm_cvttps_epi32(rs_SSE2);
            ti_SSE3            = _mm_cvttps_epi32(rs_SSE3);
#ifdef GMX_SSE4_1
            /* SSE4.1 floor is faster than _mm_cvtepi32_ps int->float cast */
            rf_SSE0            = _mm_floor_ps(rs_SSE0);
            rf_SSE1            = _mm_floor_ps(rs_SSE1);
            rf_SSE2            = _mm_floor_ps(rs_SSE2);
            rf_SSE3            = _mm_floor_ps(rs_SSE3);
#else
            rf_SSE0            = _mm_cvtepi32_ps(ti_SSE0);
            rf_SSE1            = _mm_cvtepi32_ps(ti_SSE1);
            rf_SSE2            = _mm_cvtepi32_ps(ti_SSE2);
            rf_SSE3            = _mm_cvtepi32_ps(ti_SSE3);
#endif
            frac_SSE0          = _mm_sub_ps(rs_SSE0,rf_SSE0);
            frac_SSE1          = _mm_sub_ps(rs_SSE1,rf_SSE1);
            frac_SSE2          = _mm_sub_ps(rs_SSE2,rf_SSE2);
            frac_SSE3          = _mm_sub_ps(rs_SSE3,rf_SSE3);

            /* Load 4 table floats, 2 are used with force only, 3 with energy */
            _mm_store_si128((__m128i *)ti0,ti_SSE0);
            for(e=0; e<4; e++)
            {
                ctab_SSE0[e]   = _mm_load_ps(tab_coul_FDV0+ti0[e]*4);

            }
            _mm_store_si128((__m128i *)ti1,ti_SSE1);
            for(e=0; e<4; e++)
            {
                ctab_SSE1[e]   = _mm_load_ps(tab_coul_FDV0+ti1[e]*4);

            }
            _mm_store_si128((__m128i *)ti2,ti_SSE2);
            for(e=0; e<4; e++)
            {
                ctab_SSE2[e]   = _mm_load_ps(tab_coul_FDV0+ti2[e]*4);

            }
            _mm_store_si128((__m128i *)ti3,ti_SSE3);
            for(e=0; e<4; e++)
            {
                ctab_SSE3[e]   = _mm_load_ps(tab_coul_FDV0+ti3[e]*4);

            }
            /* Shuffle the force table entries to a convenient order */
            GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE0[0],ctab_SSE0[1],ctab_SSE0[2],ctab_SSE0[3],ctab0_SSE0,ctab1_SSE0);
            GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE1[0],ctab_SSE1[1],ctab_SSE1[2],ctab_SSE1[3],ctab0_SSE1,ctab1_SSE1);
            GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE2[0],ctab_SSE2[1],ctab_SSE2[2],ctab_SSE2[3],ctab0_SSE2,ctab1_SSE2);
            GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE3[0],ctab_SSE3[1],ctab_SSE3[2],ctab_SSE3[3],ctab0_SSE3,ctab1_SSE3);
            fcoul_SSE0         = _mm_add_ps(ctab0_SSE0,_mm_mul_ps(frac_SSE0,ctab1_SSE0));
            fcoul_SSE1         = _mm_add_ps(ctab0_SSE1,_mm_mul_ps(frac_SSE1,ctab1_SSE1));
            fcoul_SSE2         = _mm_add_ps(ctab0_SSE2,_mm_mul_ps(frac_SSE2,ctab1_SSE2));
            fcoul_SSE3         = _mm_add_ps(ctab0_SSE3,_mm_mul_ps(frac_SSE3,ctab1_SSE3));
            frcoul_SSE0        = _mm_mul_ps(qq_SSE0,_mm_mul_ps(fcoul_SSE0,r_SSE0));
            frcoul_SSE1        = _mm_mul_ps(qq_SSE1,_mm_mul_ps(fcoul_SSE1,r_SSE1));
#ifndef HALF_LJ
            frcoul_SSE2        = _mm_mul_ps(qq_SSE2,_mm_mul_ps(fcoul_SSE2,r_SSE2));
            frcoul_SSE3        = _mm_mul_ps(qq_SSE3,_mm_mul_ps(fcoul_SSE3,r_SSE3));
#endif
#ifdef CALC_ENERGIES
            GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SSE0[0],ctab_SSE0[1],ctab_SSE0[2],ctab_SSE0[3],ctabv_SSE0);
            GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SSE1[0],ctab_SSE1[1],ctab_SSE1[2],ctab_SSE1[3],ctabv_SSE1);
            GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SSE2[0],ctab_SSE2[1],ctab_SSE2[2],ctab_SSE2[3],ctabv_SSE2);
            GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SSE3[0],ctab_SSE3[1],ctab_SSE3[2],ctab_SSE3[3],ctabv_SSE3);
            vcoul_SSE0 = _mm_mul_ps(qq_SSE0,_mm_add_ps(ctabv_SSE0,_mm_mul_ps(_mm_mul_ps(mhalfsp_SSE,frac_SSE0),_mm_add_ps(ctab0_SSE0,fcoul_SSE0))));
            vcoul_SSE1 = _mm_mul_ps(qq_SSE1,_mm_add_ps(ctabv_SSE1,_mm_mul_ps(_mm_mul_ps(mhalfsp_SSE,frac_SSE1),_mm_add_ps(ctab0_SSE1,fcoul_SSE1))));
            vcoul_SSE2 = _mm_mul_ps(qq_SSE2,_mm_add_ps(ctabv_SSE2,_mm_mul_ps(_mm_mul_ps(mhalfsp_SSE,frac_SSE2),_mm_add_ps(ctab0_SSE2,fcoul_SSE2))));
            vcoul_SSE3 = _mm_mul_ps(qq_SSE3,_mm_add_ps(ctabv_SSE3,_mm_mul_ps(_mm_mul_ps(mhalfsp_SSE,frac_SSE3),_mm_add_ps(ctab0_SSE3,fcoul_SSE3))));
            /*
            {
                float r[4];
                float rs[4];
                float rf[4];
                float fr[4];
                float tf[4];
                float tv[4];
                float f[4];
                float v[4];
                _mm_storeu_ps(r,r_SSE0);
                _mm_storeu_ps(rs,rs_SSE0);
                _mm_storeu_ps(rf,rf_SSE0);
                _mm_storeu_ps(fr,frac_SSE0);
                _mm_storeu_ps(tf,ctab_SSE0[0]);
                _mm_storeu_ps(tv,ctabv_SSE0);
                _mm_storeu_ps(f,fcoul_SSE0);
                _mm_storeu_ps(v,vcoul_SSE0);
                printf("r %5.3f rs %6.2f rf %5.1f fr %5.3f ri %3d tf %5.3f %5.3f tv %5.2f f %5.3f v %5.2f\n",
                       r[0],rs[0],rf[0],fr[0],ti0[0],tf[0],tf[1],tv[0],f[0],v[0]);
            }
            */
#endif
#endif
#endif

            /* Lennard-Jones interaction */
#ifdef LJ_COMB_LB
            sir_SSE0           = _mm_mul_ps(_mm_add_ps(hsig_i_SSE0,hsig_j_SSE),rinv_SSE0);
            sir_SSE1           = _mm_mul_ps(_mm_add_ps(hsig_i_SSE1,hsig_j_SSE),rinv_SSE1);
#ifndef HALF_LJ
            sir_SSE2           = _mm_mul_ps(_mm_add_ps(hsig_i_SSE2,hsig_j_SSE),rinv_SSE2);
            sir_SSE3           = _mm_mul_ps(_mm_add_ps(hsig_i_SSE3,hsig_j_SSE),rinv_SSE3);
#endif
            sir2_SSE0          = _mm_mul_ps(sir_SSE0,sir_SSE0);
            sir2_SSE1          = _mm_mul_ps(sir_SSE1,sir_SSE1);
#ifndef HALF_LJ
            sir2_SSE2          = _mm_mul_ps(sir_SSE2,sir_SSE2);
            sir2_SSE3          = _mm_mul_ps(sir_SSE3,sir_SSE3);
#endif
            sir6_SSE0          = _mm_mul_ps(sir2_SSE0,_mm_mul_ps(sir2_SSE0,sir2_SSE0));
            sir6_SSE1          = _mm_mul_ps(sir2_SSE1,_mm_mul_ps(sir2_SSE1,sir2_SSE1));
#ifndef HALF_LJ
            sir6_SSE2          = _mm_mul_ps(sir2_SSE2,_mm_mul_ps(sir2_SSE2,sir2_SSE2));
            sir6_SSE3          = _mm_mul_ps(sir2_SSE3,_mm_mul_ps(sir2_SSE3,sir2_SSE3));
#endif
            Vvdw6_SSE0         = _mm_mul_ps(_mm_mul_ps(seps_i_SSE0,seps_j_SSE),sir6_SSE0);
            Vvdw6_SSE1         = _mm_mul_ps(_mm_mul_ps(seps_i_SSE1,seps_j_SSE),sir6_SSE1);
#ifndef HALF_LJ
            Vvdw6_SSE2         = _mm_mul_ps(_mm_mul_ps(seps_i_SSE2,seps_j_SSE),sir6_SSE2);
            Vvdw6_SSE3         = _mm_mul_ps(_mm_mul_ps(seps_i_SSE3,seps_j_SSE),sir6_SSE3);
#endif
            Vvdw12_SSE0        = _mm_mul_ps(Vvdw6_SSE0,sir6_SSE0);
            Vvdw12_SSE1        = _mm_mul_ps(Vvdw6_SSE1,sir6_SSE1);
#ifndef HALF_LJ
            Vvdw12_SSE2        = _mm_mul_ps(Vvdw6_SSE2,sir6_SSE2);
            Vvdw12_SSE3        = _mm_mul_ps(Vvdw6_SSE3,sir6_SSE3);
#endif
#else /* LJ_COMB_LB */
            rinvsix_SSE0       = _mm_mul_ps(rinvsq_SSE0,_mm_mul_ps(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1       = _mm_mul_ps(rinvsq_SSE1,_mm_mul_ps(rinvsq_SSE1,rinvsq_SSE1));
#ifndef HALF_LJ
            rinvsix_SSE2       = _mm_mul_ps(rinvsq_SSE2,_mm_mul_ps(rinvsq_SSE2,rinvsq_SSE2));
            rinvsix_SSE3       = _mm_mul_ps(rinvsq_SSE3,_mm_mul_ps(rinvsq_SSE3,rinvsq_SSE3));
#endif
            Vvdw6_SSE0         = _mm_mul_ps(c6_SSE0,rinvsix_SSE0);
            Vvdw6_SSE1         = _mm_mul_ps(c6_SSE1,rinvsix_SSE1);
#ifndef HALF_LJ
            Vvdw6_SSE2         = _mm_mul_ps(c6_SSE2,rinvsix_SSE2);
            Vvdw6_SSE3         = _mm_mul_ps(c6_SSE3,rinvsix_SSE3);
#endif
            Vvdw12_SSE0        = _mm_mul_ps(c12_SSE0,_mm_mul_ps(rinvsix_SSE0,rinvsix_SSE0));
            Vvdw12_SSE1        = _mm_mul_ps(c12_SSE1,_mm_mul_ps(rinvsix_SSE1,rinvsix_SSE1));
#ifndef HALF_LJ
            Vvdw12_SSE2        = _mm_mul_ps(c12_SSE2,_mm_mul_ps(rinvsix_SSE2,rinvsix_SSE2));
            Vvdw12_SSE3        = _mm_mul_ps(c12_SSE3,_mm_mul_ps(rinvsix_SSE3,rinvsix_SSE3));
#endif
#endif /* LJ_COMB_LB */
            
#ifdef CALC_ENERGIES
#ifdef ENERGY_GROUPS
            egps_j4            = nbat->energrp[sj]*4;
#endif

#ifdef CALC_COULOMB
#ifndef ENERGY_GROUPS
            vctotSSE           = _mm_add_ps(vctotSSE, gmx_mm_sum4_ps(vcoul_SSE0,vcoul_SSE1,vcoul_SSE2,vcoul_SSE3));
#else
            vct_SSE0           = _mm_load_ps(vctp[0]+egps_j4);
            _mm_store_ps(vctp[0]+egps_j4,_mm_add_ps(vct_SSE0,vcoul_SSE0));
            vct_SSE1           = _mm_load_ps(vctp[1]+egps_j4);
            _mm_store_ps(vctp[1]+egps_j4,_mm_add_ps(vct_SSE1,vcoul_SSE1));
            vct_SSE2           = _mm_load_ps(vctp[2]+egps_j4);
            _mm_store_ps(vctp[2]+egps_j4,_mm_add_ps(vct_SSE2,vcoul_SSE2));
            vct_SSE3           = _mm_load_ps(vctp[3]+egps_j4);
            _mm_store_ps(vctp[3]+egps_j4,_mm_add_ps(vct_SSE3,vcoul_SSE3));
#endif
#endif
            Vvdw_SSE0          = _mm_sub_ps(_mm_mul_ps(twelvethSSE,Vvdw12_SSE0),_mm_mul_ps(sixthSSE,Vvdw6_SSE0));
            Vvdw_SSE1          = _mm_sub_ps(_mm_mul_ps(twelvethSSE,Vvdw12_SSE1),_mm_mul_ps(sixthSSE,Vvdw6_SSE1));
#ifndef HALF_LJ
            Vvdw_SSE2          = _mm_sub_ps(_mm_mul_ps(twelvethSSE,Vvdw12_SSE2),_mm_mul_ps(sixthSSE,Vvdw6_SSE2));
            Vvdw_SSE3          = _mm_sub_ps(_mm_mul_ps(twelvethSSE,Vvdw12_SSE3),_mm_mul_ps(sixthSSE,Vvdw6_SSE3));
#endif
#ifndef ENERGY_GROUPS
            VvdwtotSSE         = _mm_add_ps(VvdwtotSSE,
#ifndef HALF_LJ
                                            gmx_mm_sum4_ps(
#else
                                                _mm_add_ps(
#endif
                                                           Vvdw_SSE0,Vvdw_SSE1
#ifndef HALF_LJ
                                                           ,
                                                           Vvdw_SSE2,Vvdw_SSE3
#endif
                                                           ));
#else
            vvdwt_SSE0         = _mm_load_ps(vvdwtp[0]+egps_j4);
            _mm_store_ps(vvdwtp[0]+egps_j4,_mm_add_ps(vvdwt_SSE0,Vvdw_SSE0));
            vvdwt_SSE1         = _mm_load_ps(vvdwtp[1]+egps_j4);
            _mm_store_ps(vvdwtp[1]+egps_j4,_mm_add_ps(vvdwt_SSE1,Vvdw_SSE1));
#ifndef HALF_LJ
            vvdwt_SSE2         = _mm_load_ps(vvdwtp[2]+egps_j4);
            _mm_store_ps(vvdwtp[2]+egps_j4,_mm_add_ps(vvdwt_SSE2,Vvdw_SSE2));
            vvdwt_SSE3         = _mm_load_ps(vvdwtp[3]+egps_j4);
            _mm_store_ps(vvdwtp[3]+egps_j4,_mm_add_ps(vvdwt_SSE3,Vvdw_SSE3));
#endif
#endif
#endif
                                                                            
            fscal_SSE0         = _mm_mul_ps(rinvsq_SSE0,
#ifdef CALC_COULOMB
                                           _mm_add_ps(frcoul_SSE0,
#else
                                                     (
#endif
                                                      _mm_sub_ps(Vvdw12_SSE0,Vvdw6_SSE0)));
            fscal_SSE1         = _mm_mul_ps(rinvsq_SSE1,
#ifdef CALC_COULOMB
                                           _mm_add_ps(frcoul_SSE1,
#else
                                                     (
#endif
                                                      _mm_sub_ps(Vvdw12_SSE1,Vvdw6_SSE1)));
#ifndef HALF_LJ
            fscal_SSE2         = _mm_mul_ps(rinvsq_SSE2,
#ifdef CALC_COULOMB
                                           _mm_add_ps(frcoul_SSE2,
#else
                                                     (
#endif
                                                      _mm_sub_ps(Vvdw12_SSE2,Vvdw6_SSE2)));
            fscal_SSE3         = _mm_mul_ps(rinvsq_SSE3,
#ifdef CALC_COULOMB
                                           _mm_add_ps(frcoul_SSE3,
#else
                                                     (
#endif
                                                      _mm_sub_ps(Vvdw12_SSE3,Vvdw6_SSE3)));
#else
            /* Atom 2 and 3 don't have LJ, so only add Coulomb forces */
#ifdef CALC_COUL_RF
            fscal_SSE2         = _mm_mul_ps(rinvsq_SSE2,frcoul_SSE2);
            fscal_SSE3         = _mm_mul_ps(rinvsq_SSE3,frcoul_SSE3);
#else
            fscal_SSE2         = _mm_mul_ps(qq_SSE2,_mm_mul_ps(rinv_SSE2,fcoul_SSE2));
            fscal_SSE3         = _mm_mul_ps(qq_SSE3,_mm_mul_ps(rinv_SSE3,fcoul_SSE3));
#endif
#endif
            
            /* Calculate temporary vectorial force */
            tx_SSE0            = _mm_mul_ps(fscal_SSE0,dx_SSE0);
            tx_SSE1            = _mm_mul_ps(fscal_SSE1,dx_SSE1);
            tx_SSE2            = _mm_mul_ps(fscal_SSE2,dx_SSE2);
            tx_SSE3            = _mm_mul_ps(fscal_SSE3,dx_SSE3);
            ty_SSE0            = _mm_mul_ps(fscal_SSE0,dy_SSE0);
            ty_SSE1            = _mm_mul_ps(fscal_SSE1,dy_SSE1);
            ty_SSE2            = _mm_mul_ps(fscal_SSE2,dy_SSE2);
            ty_SSE3            = _mm_mul_ps(fscal_SSE3,dy_SSE3);
            tz_SSE0            = _mm_mul_ps(fscal_SSE0,dz_SSE0);
            tz_SSE1            = _mm_mul_ps(fscal_SSE1,dz_SSE1);
            tz_SSE2            = _mm_mul_ps(fscal_SSE2,dz_SSE2);
            tz_SSE3            = _mm_mul_ps(fscal_SSE3,dz_SSE3);
            
            /* Increment i atom force */
            fix_SSE0          = _mm_add_ps(fix_SSE0,tx_SSE0);
            fix_SSE1          = _mm_add_ps(fix_SSE1,tx_SSE1);
            fix_SSE2          = _mm_add_ps(fix_SSE2,tx_SSE2);
            fix_SSE3          = _mm_add_ps(fix_SSE3,tx_SSE3);
            fiy_SSE0          = _mm_add_ps(fiy_SSE0,ty_SSE0);
            fiy_SSE1          = _mm_add_ps(fiy_SSE1,ty_SSE1);
            fiy_SSE2          = _mm_add_ps(fiy_SSE2,ty_SSE2);
            fiy_SSE3          = _mm_add_ps(fiy_SSE3,ty_SSE3);
            fiz_SSE0          = _mm_add_ps(fiz_SSE0,tz_SSE0);
            fiz_SSE1          = _mm_add_ps(fiz_SSE1,tz_SSE1);
            fiz_SSE2          = _mm_add_ps(fiz_SSE2,tz_SSE2);
            fiz_SSE3          = _mm_add_ps(fiz_SSE3,tz_SSE3);
            
            /* Decrement j atom force */
            _mm_store_ps(f+ssjx,
                         _mm_sub_ps( _mm_load_ps(f+ssjx), gmx_mm_sum4_ps(tx_SSE0,tx_SSE1,tx_SSE2,tx_SSE3) ));
            _mm_store_ps(f+ssjy,
                         _mm_sub_ps( _mm_load_ps(f+ssjy), gmx_mm_sum4_ps(ty_SSE0,ty_SSE1,ty_SSE2,ty_SSE3) ));
            _mm_store_ps(f+ssjz,
                         _mm_sub_ps( _mm_load_ps(f+ssjz), gmx_mm_sum4_ps(tz_SSE0,tz_SSE1,tz_SSE2,tz_SSE3) ));
        }
