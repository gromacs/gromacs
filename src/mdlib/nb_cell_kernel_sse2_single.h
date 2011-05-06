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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "typedefs.h"

#include "gmx_sse2_single.h"


#include <xmmintrin.h>
#include <emmintrin.h>

#define SIMD_WIDTH 4
#define UNROLLI    4
#define UNROLLJ    4


/* All functionality defines are set here, except for:
 * CALC_ENERGIES, which is set before calling these functions,
 * CHECK_EXCLS, which is set just before including the inner loop contents.
 */

/* Calculate Coulomb interactions */
//#define CALC_COULOMB

/* Assumes only the first half of the particles in each cell have LJ */
//#define HALF_LJ

/* Assumes are LJ parameters are indentical */
#define FIX_LJ_C

static void
#ifndef CALC_ENERGIES
nb_cell_kernel_sse2_single_noener
#else
nb_cell_kernel_sse2_single_ener
#endif
                            (const gmx_nblist_t         *nbl,
                             const gmx_nb_atomdata_t    *nbat,
                             const t_forcerec *         fr,
                             real                       tabscale,  
                             const real *               VFtab,
                             real *                     f
#ifdef CALC_SHIFTFORCES
                             ,
                             real *                     fshift
#endif
#ifdef CALC_ENERGIES
                             ,
                             real *                     Vc,
                             real *                     Vvdw
#endif
                            )
{
    const gmx_nbl_ci_t *nbln;
    const int          *type;
    const real         *q;
    const real         *shiftvec;
    const real         *x;
    const real         *nbfp0,*nbfp1,*nbfp2,*nbfp3;
    real       facel;
    real       *nbfp_i;
    int        n,si;
    int        ish3;
    //real       shX,shY,shZ;
	int        ssi,ssix,ssiy,ssiz;
    int        sjind0,sjind1,sjind,sj,ssj,ssjx,ssjy,ssjz;
    int        ip,jp;
    real       pvdw_array[2*UNROLLI*UNROLLJ+3];
    real       *pvdw_c6,*pvdw_c12;
    
    __m128     shX_SSE;
    __m128     shY_SSE;
    __m128     shZ_SSE;
    __m128     ix_SSE0,iy_SSE0,iz_SSE0;
    __m128     ix_SSE1,iy_SSE1,iz_SSE1;
    __m128     ix_SSE2,iy_SSE2,iz_SSE2;
    __m128     ix_SSE3,iy_SSE3,iz_SSE3;
	__m128     fix_SSE0,fiy_SSE0,fiz_SSE0;
	__m128     fix_SSE1,fiy_SSE1,fiz_SSE1;
	__m128     fix_SSE2,fiy_SSE2,fiz_SSE2;
	__m128     fix_SSE3,fiy_SSE3,fiz_SSE3;

    __m128     mask0 = gmx_mm_castsi128_ps( _mm_set_epi32(0x0008, 0x0004, 0x0002, 0x0001) );
    __m128     mask1 = gmx_mm_castsi128_ps( _mm_set_epi32(0x0080, 0x0040, 0x0020, 0x0010) );
    __m128     mask2 = gmx_mm_castsi128_ps( _mm_set_epi32(0x0800, 0x0400, 0x0200, 0x0100) );
    __m128     mask3 = gmx_mm_castsi128_ps( _mm_set_epi32(0x8000, 0x4000, 0x2000, 0x1000) );
    __m128     zero_SSE = gmx_mm_castsi128_ps( _mm_set_epi32(0x0, 0x0, 0x0, 0x0) );
    __m128     mask_int;
    __m128     int_SSE0;
    __m128     int_SSE1;
    __m128     int_SSE2;
    __m128     int_SSE3;

	__m128     jxSSE,jySSE,jzSSE,jqSSE;
	__m128     dx_SSE0,dy_SSE0,dz_SSE0;
	__m128     dx_SSE1,dy_SSE1,dz_SSE1;
	__m128     dx_SSE2,dy_SSE2,dz_SSE2;
	__m128     dx_SSE3,dy_SSE3,dz_SSE3;
	__m128     tx_SSE0,ty_SSE0,tz_SSE0;
	__m128     tx_SSE1,ty_SSE1,tz_SSE1;
	__m128     tx_SSE2,ty_SSE2,tz_SSE2;
	__m128     tx_SSE3,ty_SSE3,tz_SSE3;
	__m128     rsq_SSE0,rinv_SSE0,rinvsq_SSE0,rinvsix_SSE0;
	__m128     rsq_SSE1,rinv_SSE1,rinvsq_SSE1,rinvsix_SSE1;
	__m128     rsq_SSE2,rinv_SSE2,rinvsq_SSE2,rinvsix_SSE2;
	__m128     rsq_SSE3,rinv_SSE3,rinvsq_SSE3,rinvsix_SSE3;
    __m128     wco_SSE0;
    __m128     wco_SSE1;
    __m128     wco_SSE2;
    __m128     wco_SSE3;
    __m128     wco_any_SSE01,wco_any_SSE23,wco_any_SSE;
	__m128     qq_SSE0,iq_SSE0;
	__m128     qq_SSE1,iq_SSE1;
	__m128     qq_SSE2,iq_SSE2;
	__m128     qq_SSE3,iq_SSE3;
	__m128     vcoul_SSE0,Vvdw6_SSE0,Vvdw12_SSE0,fscal_SSE0;
	__m128     vcoul_SSE1,Vvdw6_SSE1,Vvdw12_SSE1,fscal_SSE1;
	__m128     vcoul_SSE2,Vvdw6_SSE2,Vvdw12_SSE2,fscal_SSE2;
	__m128     vcoul_SSE3,Vvdw6_SSE3,Vvdw12_SSE3,fscal_SSE3;
    __m128     c6_SSE0,c12_SSE0;
    __m128     c6_SSE1,c12_SSE1;
    __m128     c6_SSE2,c12_SSE2;
    __m128     c6_SSE3,c12_SSE3;
    
    __m128     vctotSSE,VvdwtotSSE;
    __m128     sixSSE,twelveSSE;
    __m128i    ikSSE,imSSE,ifourSSE;
    __m128i    ioneSSE;

    __m128     rc2_SSE;

    float      wco_any_array[7],*wco_any_align;
    float      mask_array[7],*mask_align;
#ifdef CALC_ENERGIES
    float      tmpsum_array[7],*tmpsum;
#endif
#ifdef CALC_SHIFTFORCES
    float      shf_array[7],*shf;
#endif

    int ninner;

#ifdef COUNT_PAIRS
    int npair=0;
#endif

    pvdw_c6  = (float *)(((size_t)(pvdw_array+3)) & (~((size_t)15)));
    pvdw_c12 = pvdw_c6 + UNROLLI*UNROLLJ;

    q                   = nbat->q;
    type                = nbat->type;
    facel               = fr->epsfac;
    shiftvec            = fr->shift_vec[0];
    x                   = nbat->x;
    
    sixSSE    = _mm_set1_ps(6.0);
    twelveSSE = _mm_set1_ps(12.0);
    ifourSSE  = _mm_set1_epi32(4);
    ioneSSE   = _mm_set1_epi32(1);

    rc2_SSE   = _mm_set1_ps(nbl->rcut*nbl->rcut);

    wco_any_align = (float *)(((size_t)(wco_any_array+3)) & (~((size_t)15)));

    mask_align = (float *)(((size_t)(mask_array+3)) & (~((size_t)15)));

#ifdef CALC_ENERGIES
    tmpsum = (float *)(((size_t)(tmpsum_array+3)) & (~((size_t)15)));
#endif
#ifdef CALC_ENERGIES
    shf = (float *)(((size_t)(shf_array+3)) & (~((size_t)15)));
#endif

#ifdef FIX_LJ_C
    for(jp=0; jp<UNROLLJ; jp++)
    {
        pvdw_c6 [0*UNROLLJ+jp] = nbat->nbfp[0*2];
        pvdw_c6 [1*UNROLLJ+jp] = nbat->nbfp[0*2];
#ifndef HALF_LJ
        pvdw_c6 [2*UNROLLJ+jp] = nbat->nbfp[0*2];
        pvdw_c6 [3*UNROLLJ+jp] = nbat->nbfp[0*2];
#else
        pvdw_c6 [2*UNROLLJ+jp] = 0;
        pvdw_c6 [3*UNROLLJ+jp] = 0;
#endif
        pvdw_c12[0*UNROLLJ+jp] = nbat->nbfp[0*2+1];
        pvdw_c12[1*UNROLLJ+jp] = nbat->nbfp[0*2+1];
#ifndef HALF_LJ
        pvdw_c12[2*UNROLLJ+jp] = nbat->nbfp[0*2+1];
        pvdw_c12[3*UNROLLJ+jp] = nbat->nbfp[0*2+1];
#else
        pvdw_c12[2*UNROLLJ+jp] = 0;
        pvdw_c12[3*UNROLLJ+jp] = 0;
#endif
    }
    c6_SSE0            = _mm_load_ps(pvdw_c6 +0*UNROLLJ);
    c6_SSE1            = _mm_load_ps(pvdw_c6 +1*UNROLLJ);
    c6_SSE2            = _mm_load_ps(pvdw_c6 +2*UNROLLJ);
    c6_SSE3            = _mm_load_ps(pvdw_c6 +3*UNROLLJ);
    c12_SSE0           = _mm_load_ps(pvdw_c12+0*UNROLLJ);
    c12_SSE1           = _mm_load_ps(pvdw_c12+1*UNROLLJ);
    c12_SSE2           = _mm_load_ps(pvdw_c12+2*UNROLLJ);
    c12_SSE3           = _mm_load_ps(pvdw_c12+3*UNROLLJ);
#endif

    ninner = 0;
    for(n=0; n<nbl->nci; n++)
    {
        nbln = &nbl->ci[n];

        ish3             = 3*nbln->shift;     
        sjind0           = nbln->cj_ind_start;      
        sjind1           = nbln->cj_ind_end;    
        /* Currently only works super-cells equal to sub-cells */
        si               = nbln->ci;

        shX_SSE = _mm_load1_ps(shiftvec+ish3);
        shY_SSE = _mm_load1_ps(shiftvec+ish3+1);
        shZ_SSE = _mm_load1_ps(shiftvec+ish3+2);
       
		/* Load i atom data */
        ssi              = si*SIMD_WIDTH;
        ssix             = ssi*DIM;
        ssiy             = ssix + SIMD_WIDTH;
        ssiz             = ssiy + SIMD_WIDTH;
		ix_SSE0          = _mm_add_ps(_mm_load1_ps(x+ssix)  ,shX_SSE);
		ix_SSE1          = _mm_add_ps(_mm_load1_ps(x+ssix+1),shX_SSE);
		ix_SSE2          = _mm_add_ps(_mm_load1_ps(x+ssix+2),shX_SSE);
		ix_SSE3          = _mm_add_ps(_mm_load1_ps(x+ssix+3),shX_SSE);
		iy_SSE0          = _mm_add_ps(_mm_load1_ps(x+ssiy)  ,shY_SSE);
		iy_SSE1          = _mm_add_ps(_mm_load1_ps(x+ssiy+1),shY_SSE);
		iy_SSE2          = _mm_add_ps(_mm_load1_ps(x+ssiy+2),shY_SSE);
		iy_SSE3          = _mm_add_ps(_mm_load1_ps(x+ssiy+3),shY_SSE);
        iz_SSE0          = _mm_add_ps(_mm_load1_ps(x+ssiz)  ,shZ_SSE);
		iz_SSE1          = _mm_add_ps(_mm_load1_ps(x+ssiz+1),shZ_SSE);
		iz_SSE2          = _mm_add_ps(_mm_load1_ps(x+ssiz+2),shZ_SSE);
		iz_SSE3          = _mm_add_ps(_mm_load1_ps(x+ssiz+3),shZ_SSE);
#ifdef CALC_COULOMB
		iq_SSE0          = _mm_set1_ps(facel*q[ssi]);
		iq_SSE1          = _mm_set1_ps(facel*q[ssi+1]);
		iq_SSE2          = _mm_set1_ps(facel*q[ssi+2]);
		iq_SSE3          = _mm_set1_ps(facel*q[ssi+3]);
#endif

        nbfp0 = nbat->nbfp + type[ssi  ]*nbat->ntype*2;
        nbfp1 = nbat->nbfp + type[ssi+1]*nbat->ntype*2;
        nbfp2 = nbat->nbfp + type[ssi+2]*nbat->ntype*2;
        nbfp3 = nbat->nbfp + type[ssi+3]*nbat->ntype*2;

		/* Zero the potential energy for this list */
		VvdwtotSSE       = _mm_setzero_ps();
		vctotSSE         = _mm_setzero_ps();

		/* Clear i atom forces */
		fix_SSE0           = _mm_setzero_ps();
		fix_SSE1           = _mm_setzero_ps();
		fix_SSE2           = _mm_setzero_ps();
		fix_SSE3           = _mm_setzero_ps();
		fiy_SSE0           = _mm_setzero_ps();
		fiy_SSE1           = _mm_setzero_ps();
		fiy_SSE2           = _mm_setzero_ps();
		fiy_SSE3           = _mm_setzero_ps();
		fiz_SSE0           = _mm_setzero_ps();
		fiz_SSE1           = _mm_setzero_ps();
		fiz_SSE2           = _mm_setzero_ps();
		fiz_SSE3           = _mm_setzero_ps();

        sjind = sjind0;
        while (sjind < sjind1 && nbl->cj[sjind].excl != 0xffff)
        {
#define CHECK_EXCLS
#include "nb_cell_kernel_sse2_single_inner.h"
#undef CHECK_EXCLS
            sjind++;
        }

        for(; (sjind<sjind1); sjind++)
        {
#include "nb_cell_kernel_sse2_single_inner.h"
        }
        ninner += sjind1 - sjind0;

	/* Add i forces to mem and shiftedn force list */
        _MM_TRANSPOSE4_PS(fix_SSE0,fix_SSE1,fix_SSE2,fix_SSE3);
        fix_SSE0 = _mm_add_ps(fix_SSE0,fix_SSE1);
        fix_SSE2 = _mm_add_ps(fix_SSE2,fix_SSE3);
        fix_SSE0 = _mm_add_ps(fix_SSE0,fix_SSE2);
        _mm_store_ps(f+ssix, _mm_add_ps(fix_SSE0, _mm_load_ps(f+ssix)));
        
        _MM_TRANSPOSE4_PS(fiy_SSE0,fiy_SSE1,fiy_SSE2,fiy_SSE3);
        fiy_SSE0 = _mm_add_ps(fiy_SSE0,fiy_SSE1);
        fiy_SSE2 = _mm_add_ps(fiy_SSE2,fiy_SSE3);
        fiy_SSE0 = _mm_add_ps(fiy_SSE0,fiy_SSE2);
        _mm_store_ps(f+ssiy, _mm_add_ps(fiy_SSE0, _mm_load_ps(f+ssiy)));
        
        _MM_TRANSPOSE4_PS(fiz_SSE0,fiz_SSE1,fiz_SSE2,fiz_SSE3);
        fiz_SSE0 = _mm_add_ps(fiz_SSE0,fiz_SSE1);
        fiz_SSE2 = _mm_add_ps(fiz_SSE2,fiz_SSE3);
        fiz_SSE0 = _mm_add_ps(fiz_SSE0,fiz_SSE2);
        _mm_store_ps(f+ssiz, _mm_add_ps(fiz_SSE0, _mm_load_ps(f+ssiz)));

#ifdef CALC_SHIFTFORCES
        _mm_store_ps(shf,fix_SSE0);
        fshift[ish3+0] += shf[0]+shf[1]+shf[2]+shf[3]; 
        _mm_store_ps(shf,fiy_SSE0);
        fshift[ish3+1] += shf[0]+shf[1]+shf[2]+shf[3]; 
        _mm_store_ps(shf,fiz_SSE0);
        fshift[ish3+2] += shf[0]+shf[1]+shf[2]+shf[3]; 
#endif
		
#ifdef CALC_ENERGIES
#ifdef CALC_COULOMB
        _mm_store_ps(tmpsum,vctotSSE);
        *Vc += tmpsum[0]+tmpsum[1]+tmpsum[2]+tmpsum[3];
#endif
		
        _mm_store_ps(tmpsum,VvdwtotSSE);
        *Vvdw += tmpsum[0]+tmpsum[1]+tmpsum[2]+tmpsum[3];
#endif
		
		/* Outer loop uses 6 flops/iteration */
	}

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n",npair);
#endif
}

#undef SIMD_WIDTH
#undef UNROLLI   
#undef UNROLLJ   
