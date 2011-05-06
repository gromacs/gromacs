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

            sj               = nbl->cj[sjind].c;

            ssj              = sj*SIMD_WIDTH;
            ssjx             = ssj*DIM;
            ssjy             = ssj*DIM + SIMD_WIDTH;
            ssjz             = ssj*DIM + 2*SIMD_WIDTH;

#ifdef CHECK_EXCLS
            /* Load integer interaction mask as float to avoid sse casts */
            mask_int         = _mm_load1_ps((float *)&nbl->cj[sjind].excl);
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
            jqSSE              = _mm_load_ps(q+ssj);
            qq_SSE0            = _mm_mul_ps(iq_SSE0,jqSSE);
            qq_SSE1            = _mm_mul_ps(iq_SSE1,jqSSE);
            qq_SSE2            = _mm_mul_ps(iq_SSE2,jqSSE);
            qq_SSE3            = _mm_mul_ps(iq_SSE3,jqSSE);
#endif
            
#ifndef FIX_LJ_C
            for(jp=0; jp<UNROLLJ; jp++)
            {
                pvdw_c6 [0*UNROLLJ+jp] = nbfp0[type[ssj+0]*2];
                pvdw_c6 [1*UNROLLJ+jp] = nbfp1[type[ssj+1]*2];
#ifndef HALF_LJ
                pvdw_c6 [2*UNROLLJ+jp] = nbfp2[type[ssj+2]*2];
                pvdw_c6 [3*UNROLLJ+jp] = nbfp3[type[ssj+3]*2];
#endif
                pvdw_c12[0*UNROLLJ+jp] = nbfp0[type[ssj+0]*2+1];
                pvdw_c12[1*UNROLLJ+jp] = nbfp1[type[ssj+1]*2+1];
#ifndef HALF_LJ
                pvdw_c12[2*UNROLLJ+jp] = nbfp2[type[ssj+2]*2+1];
                pvdw_c12[3*UNROLLJ+jp] = nbfp3[type[ssj+3]*2+1];
#endif
            }

            c6_SSE0            = _mm_load_ps(pvdw_c6 +0*UNROLLJ);
            c6_SSE1            = _mm_load_ps(pvdw_c6 +1*UNROLLJ);
#ifndef HALF_LJ
            c6_SSE2            = _mm_load_ps(pvdw_c6 +2*UNROLLJ);
            c6_SSE3            = _mm_load_ps(pvdw_c6 +3*UNROLLJ);
#endif
            c12_SSE0           = _mm_load_ps(pvdw_c12+0*UNROLLJ);
            c12_SSE1           = _mm_load_ps(pvdw_c12+1*UNROLLJ);
#ifndef HALF_LJ
            c12_SSE2           = _mm_load_ps(pvdw_c12+2*UNROLLJ);
            c12_SSE3           = _mm_load_ps(pvdw_c12+3*UNROLLJ);
#endif
#endif /* FIX_LJ_C */
            
            rinv_SSE0          = _mm_and_ps(rinv_SSE0,wco_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1,wco_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2,wco_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3,wco_SSE3);

            rinvsq_SSE0        = _mm_mul_ps(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1        = _mm_mul_ps(rinv_SSE1,rinv_SSE1);
            rinvsq_SSE2        = _mm_mul_ps(rinv_SSE2,rinv_SSE2);
            rinvsq_SSE3        = _mm_mul_ps(rinv_SSE3,rinv_SSE3);

#ifdef CALC_COULOMB
            /* Coulomb interaction */
            vcoul_SSE0         = _mm_mul_ps(qq_SSE0,rinv_SSE0);
            vcoul_SSE1         = _mm_mul_ps(qq_SSE1,rinv_SSE1);
            vcoul_SSE2         = _mm_mul_ps(qq_SSE2,rinv_SSE2);
            vcoul_SSE3         = _mm_mul_ps(qq_SSE3,rinv_SSE3);
#endif            

            /* Lennard-Jones interaction */
            rinvsix_SSE0       = _mm_mul_ps(rinvsq_SSE0,_mm_mul_ps(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1       = _mm_mul_ps(rinvsq_SSE1,_mm_mul_ps(rinvsq_SSE1,rinvsq_SSE1));
            rinvsix_SSE2       = _mm_mul_ps(rinvsq_SSE2,_mm_mul_ps(rinvsq_SSE2,rinvsq_SSE2));
            rinvsix_SSE3       = _mm_mul_ps(rinvsq_SSE3,_mm_mul_ps(rinvsq_SSE3,rinvsq_SSE3));
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
            
#ifdef CALC_ENERGIES
#ifdef CALC_COULOMB
            vctotSSE           = _mm_add_ps(vctotSSE, gmx_mm_sum4_ps(vcoul_SSE0,vcoul_SSE1,vcoul_SSE2,vcoul_SSE3));
#endif
            VvdwtotSSE         = _mm_add_ps(VvdwtotSSE, gmx_mm_sum4_ps(_mm_sub_ps(Vvdw12_SSE0,Vvdw6_SSE0),
                                                                    _mm_sub_ps(Vvdw12_SSE1,Vvdw6_SSE1),
                                                                    _mm_sub_ps(Vvdw12_SSE2,Vvdw6_SSE2),
                                                                    _mm_sub_ps(Vvdw12_SSE3,Vvdw6_SSE3)));
#endif
                                                                            
            fscal_SSE0         = _mm_mul_ps(rinvsq_SSE0,
#ifdef CALC_COULOMB
                                           _mm_add_ps(vcoul_SSE0,
#else
                                                     (
#endif
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw12_SSE0),
                                                                 _mm_mul_ps(sixSSE,Vvdw6_SSE0))));
            fscal_SSE1         = _mm_mul_ps(rinvsq_SSE1,
#ifdef CALC_COULOMB
                                           _mm_add_ps(vcoul_SSE1,
#else
                                                     (
#endif
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw12_SSE1),
                                                                 _mm_mul_ps(sixSSE,Vvdw6_SSE1))));
#ifndef HALF_LJ
            fscal_SSE2         = _mm_mul_ps(rinvsq_SSE2,
#ifdef CALC_COULOMB
                                           _mm_add_ps(vcoul_SSE2,
#else
                                                     (
#endif
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw12_SSE2),
                                                                 _mm_mul_ps(sixSSE,Vvdw6_SSE2))));
            fscal_SSE3         = _mm_mul_ps(rinvsq_SSE3,
#ifdef CALC_COULOMB
                                           _mm_add_ps(vcoul_SSE3,
#else
                                                     (
#endif
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw12_SSE3),
                                                                 _mm_mul_ps(sixSSE,Vvdw6_SSE3))));
#else
            fscal_SSE2         = _mm_mul_ps(rinvsq_SSE2,vcoul_SSE2);
            fscal_SSE3         = _mm_mul_ps(rinvsq_SSE3,vcoul_SSE3);
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
            
            /* Inner loop uses 38 flops/iteration */
