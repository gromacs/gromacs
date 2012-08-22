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

/* GMX_MM128_HERE or GMX_MM256_HERE should be set before including this file */
#include "gmx_x86_simd_macros.h"

#define SUM_SIMD4(x) (x[0]+x[1]+x[2]+x[3])

#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    GMX_X86_SIMD_WIDTH_HERE

#if defined GMX_MM128_HERE || defined GMX_DOUBLE
#define STRIDE     4
#endif
#if defined GMX_MM256_HERE && !defined GMX_DOUBLE
#define STRIDE     8
#endif 

#ifdef GMX_MM128_HERE
#ifndef GMX_DOUBLE
/* SSE single precision 4x4 kernel */
#define SUM_SIMD(x) SUM_SIMD4(x)
#define TAB_FDV0
#else
/* SSE double precision 4x2 kernel */
#define SUM_SIMD(x) (x[0]+x[1])
#endif
#endif

#ifdef GMX_MM256_HERE
#ifndef GMX_DOUBLE
/* AVX single precision 4x8 kernel */
#define SUM_SIMD(x) (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7])
#define TAB_FDV0
#else
/* AVX double precision 4x4 kernel */
#define SUM_SIMD(x) SUM_SIMD4(x)
#endif
#endif

#define SIMD_MASK_ALL   0xffffffff

#include "nbnxn_kernel_x86_simd_utils.h"

/* All functionality defines are set here, except for:
 * CALC_ENERGIES, ENERGY_GROUPS which are defined before.
 * CHECK_EXCLS, which is set just before including the inner loop contents.
 * The combination rule defines, LJ_COMB_GEOM or LJ_COMB_LB are currently
 * set before calling the kernel function. We might want to move that
 * to inside the n-loop and have a different combination rule for different
 * ci's, as no combination rule gives a 50% performance hit for LJ.
 */

/* We always calculate shift forces, because it's cheap anyhow */
#define CALC_SHIFTFORCES

/* Assumes all LJ parameters are identical */
/* #define FIX_LJ_C */

#define NBK_FUNC_NAME_C_LJC(b,s,c,ljc,e) b##_##s##_##c##_comb_##ljc##_##e

#if defined LJ_COMB_GEOM
#define NBK_FUNC_NAME_C(b,s,c,e) NBK_FUNC_NAME_C_LJC(b,s,c,geom,e)
#else
#if defined LJ_COMB_LB
#define NBK_FUNC_NAME_C(b,s,c,e) NBK_FUNC_NAME_C_LJC(b,s,c,lb,e)
#else
#define NBK_FUNC_NAME_C(b,s,c,e) NBK_FUNC_NAME_C_LJC(b,s,c,none,e)
#endif
#endif

#ifdef CALC_COUL_RF
#define NBK_FUNC_NAME(b,s,e) NBK_FUNC_NAME_C(b,s,rf,e)
#endif
#ifdef CALC_COUL_TAB
#define NBK_FUNC_NAME(b,s,e) NBK_FUNC_NAME_C(b,s,tab,e)
#endif

#ifdef GMX_MM128_HERE
#define NBK_FUNC_NAME_S128_OR_S256(b,e) NBK_FUNC_NAME(b,x86_simd128,e)
#endif
#ifdef GMX_MM256_HERE
#define NBK_FUNC_NAME_S128_OR_S256(b,e) NBK_FUNC_NAME(b,x86_simd256,e)
#endif

static void
#ifndef CALC_ENERGIES
NBK_FUNC_NAME_S128_OR_S256(nbnxn_kernel,noener)
#else
#ifndef ENERGY_GROUPS
NBK_FUNC_NAME_S128_OR_S256(nbnxn_kernel,ener)
#else
NBK_FUNC_NAME_S128_OR_S256(nbnxn_kernel,energrp)
#endif
#endif
#undef NBK_FUNC_NAME
#undef NBK_FUNC_NAME_C
#undef NBK_FUNC_NAME_C_LJC
                            (const nbnxn_pairlist_t     *nbl,
                             const nbnxn_atomdata_t     *nbat,
                             const interaction_const_t  *ic,
                             rvec                       *shift_vec, 
                             real                       *f
#ifdef CALC_SHIFTFORCES
                             ,
                             real                       *fshift
#endif
#ifdef CALC_ENERGIES
                             ,
                             real                       *Vvdw,
                             real                       *Vc
#endif
                            )
{
    const nbnxn_ci_t   *nbln;
    const nbnxn_cj_t   *l_cj;
    const int          *type;
    const real         *q;
    const real         *shiftvec;
    const real         *x;
    const real         *nbfp0,*nbfp1,*nbfp2=NULL,*nbfp3=NULL;
    real       facel;
    real       *nbfp_ptr;
    int        nbfp_stride;
    int        n,ci,ci_sh;
    int        ish,ish3;
    gmx_bool   half_LJ,do_coul;
    int        sci,scix,sciy,sciz,sci2;
    int        cjind0,cjind1,cjind;
    int        ip,jp;

#ifdef ENERGY_GROUPS
    int        Vstride_i;
    int        egps_ishift,egps_imask;
    int        egps_jshift,egps_jmask,egps_jstride;
    int        egps_i;
    real       *vvdwtp[UNROLLI];
    real       *vctp[UNROLLI];
#endif
    
    gmx_mm_pr  shX_SSE;
    gmx_mm_pr  shY_SSE;
    gmx_mm_pr  shZ_SSE;
    gmx_mm_pr  ix_SSE0,iy_SSE0,iz_SSE0;
    gmx_mm_pr  ix_SSE1,iy_SSE1,iz_SSE1;
    gmx_mm_pr  ix_SSE2,iy_SSE2,iz_SSE2;
    gmx_mm_pr  ix_SSE3,iy_SSE3,iz_SSE3;
    gmx_mm_pr  fix_SSE0,fiy_SSE0,fiz_SSE0;
    gmx_mm_pr  fix_SSE1,fiy_SSE1,fiz_SSE1;
    gmx_mm_pr  fix_SSE2,fiy_SSE2,fiz_SSE2;
    gmx_mm_pr  fix_SSE3,fiy_SSE3,fiz_SSE3;
#if UNROLLJ >= 4
#ifndef GMX_DOUBLE
    __m128     fix_SSE,fiy_SSE,fiz_SSE;
#else
    __m256d    fix_SSE,fiy_SSE,fiz_SSE;
#endif
#else
    __m128d    fix0_SSE,fiy0_SSE,fiz0_SSE;
    __m128d    fix2_SSE,fiy2_SSE,fiz2_SSE;
#endif

#ifndef GMX_MM256_HERE
#ifndef GMX_DOUBLE
    __m128i    mask0 = _mm_set_epi32( 0x0008, 0x0004, 0x0002, 0x0001 );
    __m128i    mask1 = _mm_set_epi32( 0x0080, 0x0040, 0x0020, 0x0010 );
    __m128i    mask2 = _mm_set_epi32( 0x0800, 0x0400, 0x0200, 0x0100 );
    __m128i    mask3 = _mm_set_epi32( 0x8000, 0x4000, 0x2000, 0x1000 );
#else
    /* For double precision we need to set two 32bit ints for one double */
    __m128i    mask0 = _mm_set_epi32( 0x0002, 0x0002, 0x0001, 0x0001 );
    __m128i    mask1 = _mm_set_epi32( 0x0008, 0x0008, 0x0004, 0x0004 );
    __m128i    mask2 = _mm_set_epi32( 0x0020, 0x0020, 0x0010, 0x0010 );
    __m128i    mask3 = _mm_set_epi32( 0x0080, 0x0080, 0x0040, 0x0040 );
#endif
#else
    /* AVX: use floating point masks, as there are no integer instructions */
#ifndef GMX_DOUBLE
    gmx_mm_pr  mask0 = _mm256_castsi256_ps(_mm256_set_epi32( 0x0080, 0x0040, 0x0020, 0x0010, 0x0008, 0x0004, 0x0002, 0x0001 ));
    gmx_mm_pr  mask1 = _mm256_castsi256_ps(_mm256_set_epi32( 0x8000, 0x4000, 0x2000, 0x1000, 0x0800, 0x0400, 0x0200, 0x0100 ));
#else
    /* There is no 256-bit int to double conversion, so we use float here */
    __m256     mask0 = _mm256_castsi256_ps(_mm256_set_epi32( 0x0008, 0x0008, 0x0004, 0x0004, 0x0002, 0x0002, 0x0001, 0x0001 ));
    __m256     mask1 = _mm256_castsi256_ps(_mm256_set_epi32( 0x0080, 0x0080, 0x0040, 0x0040, 0x0020, 0x0020, 0x0010, 0x0010 ));
    __m256     mask2 = _mm256_castsi256_ps(_mm256_set_epi32( 0x0800, 0x0800, 0x0400, 0x0400, 0x0200, 0x0200, 0x0100, 0x0100 ));
    __m256     mask3 = _mm256_castsi256_ps(_mm256_set_epi32( 0x8000, 0x8000, 0x4000, 0x4000, 0x2000, 0x2000, 0x1000, 0x1000 ));
#endif
#endif

#ifndef GMX_MM256_HERE
#ifndef GMX_DOUBLE
    __m128     diag_SSE0 = gmx_mm_castsi128_pr( _mm_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000 ));
    __m128     diag_SSE1 = gmx_mm_castsi128_pr( _mm_set_epi32( 0xffffffff, 0xffffffff, 0x00000000, 0x00000000 ));
    __m128     diag_SSE2 = gmx_mm_castsi128_pr( _mm_set_epi32( 0xffffffff, 0x00000000, 0x00000000, 0x00000000 ));
    __m128     diag_SSE3 = gmx_mm_castsi128_pr( _mm_set_epi32( 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
#else
    __m128d    diag0_SSE0 = gmx_mm_castsi128_pd( _mm_set_epi32( 0xffffffff, 0xffffffff, 0x00000000, 0x00000000 ));
    __m128d    diag0_SSE1 = gmx_mm_castsi128_pd( _mm_set_epi32( 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    __m128d    diag0_SSE2 = gmx_mm_castsi128_pd( _mm_set_epi32( 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    __m128d    diag0_SSE3 = gmx_mm_castsi128_pd( _mm_set_epi32( 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    __m128d    diag1_SSE0 = gmx_mm_castsi128_pd( _mm_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff ));
    __m128d    diag1_SSE1 = gmx_mm_castsi128_pd( _mm_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff ));
    __m128d    diag1_SSE2 = gmx_mm_castsi128_pd( _mm_set_epi32( 0xffffffff, 0xffffffff, 0x00000000, 0x00000000 ));
    __m128d    diag1_SSE3 = gmx_mm_castsi128_pd( _mm_set_epi32( 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
#endif
#else /* GMX_MM256_HERE */
#ifndef GMX_DOUBLE
    gmx_mm_pr  diag0_SSE0 = _mm256_castsi256_ps( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000 ));
    gmx_mm_pr  diag0_SSE1 = _mm256_castsi256_ps( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag0_SSE2 = _mm256_castsi256_ps( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag0_SSE3 = _mm256_castsi256_ps( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag1_SSE0 = _mm256_castsi256_ps( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag1_SSE1 = _mm256_castsi256_ps( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag1_SSE2 = _mm256_castsi256_ps( _mm256_set_epi32( 0xffffffff, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag1_SSE3 = _mm256_castsi256_ps( _mm256_set_epi32( 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
#else
    gmx_mm_pr  diag_SSE0 = _mm256_castsi256_pd( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag_SSE1 = _mm256_castsi256_pd( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag_SSE2 = _mm256_castsi256_pd( _mm256_set_epi32( 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
    gmx_mm_pr  diag_SSE3 = _mm256_castsi256_pd( _mm256_set_epi32( 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 ));
#endif
#endif

#ifndef GMX_MM256_HERE
    __m128i    zeroi_SSE = _mm_setzero_si128();
#endif
#ifdef GMX_X86_SSE4_1
    gmx_mm_pr  zero_SSE = gmx_set1_pr(0);
#endif

    gmx_mm_pr  one_SSE=gmx_set1_pr(1.0);
    gmx_mm_pr  iq_SSE0=gmx_setzero_pr();
    gmx_mm_pr  iq_SSE1=gmx_setzero_pr();
    gmx_mm_pr  iq_SSE2=gmx_setzero_pr();
    gmx_mm_pr  iq_SSE3=gmx_setzero_pr();
    gmx_mm_pr  mrc_3_SSE;
#ifdef CALC_ENERGIES
    gmx_mm_pr  hrc_3_SSE,moh_rc_SSE;
#endif
#ifdef CALC_COUL_TAB
    /* Coulomb table variables */
    gmx_mm_pr  invtsp_SSE;
    const real *tab_coul_F;
#ifndef TAB_FDV0
    const real *tab_coul_V;
#endif
#ifdef GMX_MM256_HERE
    int        ti0_array[2*UNROLLJ-1],*ti0;
    int        ti1_array[2*UNROLLJ-1],*ti1;
    int        ti2_array[2*UNROLLJ-1],*ti2;
    int        ti3_array[2*UNROLLJ-1],*ti3;
#endif
#ifdef CALC_ENERGIES
    gmx_mm_pr  mhalfsp_SSE;
    gmx_mm_pr  sh_ewald_SSE;
#endif
#endif

#ifdef LJ_COMB_LB
    const real *ljc;

    gmx_mm_pr  hsig_i_SSE0,seps_i_SSE0;
    gmx_mm_pr  hsig_i_SSE1,seps_i_SSE1;
    gmx_mm_pr  hsig_i_SSE2,seps_i_SSE2;
    gmx_mm_pr  hsig_i_SSE3,seps_i_SSE3;
#else
#ifdef FIX_LJ_C
    real       pvdw_array[2*UNROLLI*UNROLLJ+3];
    real       *pvdw_c6,*pvdw_c12;
    gmx_mm_pr  c6_SSE0,c12_SSE0;
    gmx_mm_pr  c6_SSE1,c12_SSE1;
    gmx_mm_pr  c6_SSE2,c12_SSE2;
    gmx_mm_pr  c6_SSE3,c12_SSE3;
#endif

#ifdef LJ_COMB_GEOM
    const real *ljc;

    gmx_mm_pr  c6s_SSE0,c12s_SSE0;
    gmx_mm_pr  c6s_SSE1,c12s_SSE1;
    gmx_mm_pr  c6s_SSE2=gmx_setzero_pr(),c12s_SSE2=gmx_setzero_pr();
    gmx_mm_pr  c6s_SSE3=gmx_setzero_pr(),c12s_SSE3=gmx_setzero_pr();
#endif
#endif /* LJ_COMB_LB */

    gmx_mm_pr  vctotSSE,VvdwtotSSE;
    gmx_mm_pr  sixthSSE,twelvethSSE;

    gmx_mm_pr  avoid_sing_SSE;
    gmx_mm_pr  rc2_SSE;

#ifdef CALC_ENERGIES
    gmx_mm_pr  sh_invrc6_SSE,sh_invrc12_SSE;

    /* cppcheck-suppress unassignedVariable */
    real       tmpsum_array[15],*tmpsum;
#endif
#ifdef CALC_SHIFTFORCES
    /* cppcheck-suppress unassignedVariable */
    real       shf_array[15],*shf;
#endif

    int ninner;

#ifdef COUNT_PAIRS
    int npair=0;
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB
    ljc = nbat->lj_comb;
#else
    /* No combination rule used */
#ifndef GMX_DOUBLE
    nbfp_ptr    = nbat->nbfp_s4;
    nbfp_stride = 4;
#else
    nbfp_ptr    = nbat->nbfp;
    nbfp_stride = 2;
#endif
#endif

#ifdef CALC_COUL_TAB
#ifdef GMX_MM256_HERE
    /* Generate aligned table pointers */
    ti0 = (int *)(((size_t)(ti0_array+UNROLLJ-1)) & (~((size_t)(UNROLLJ*sizeof(real)-1))));
    ti1 = (int *)(((size_t)(ti1_array+UNROLLJ-1)) & (~((size_t)(UNROLLJ*sizeof(real)-1))));
    ti2 = (int *)(((size_t)(ti2_array+UNROLLJ-1)) & (~((size_t)(UNROLLJ*sizeof(real)-1))));
    ti3 = (int *)(((size_t)(ti3_array+UNROLLJ-1)) & (~((size_t)(UNROLLJ*sizeof(real)-1))));
#endif

    invtsp_SSE  = gmx_set1_pr(ic->tabq_scale);
#ifdef CALC_ENERGIES
    mhalfsp_SSE = gmx_set1_pr(-0.5/ic->tabq_scale);

    sh_ewald_SSE = gmx_set1_pr(ic->sh_ewald);
#endif

#ifdef TAB_FDV0
    tab_coul_F = ic->tabq_coul_FDV0;
#else
    tab_coul_F = ic->tabq_coul_F;
    tab_coul_V = ic->tabq_coul_V;
#endif
#endif

    q                   = nbat->q;
    type                = nbat->type;
    facel               = ic->epsfac;
    shiftvec            = shift_vec[0];
    x                   = nbat->x;

    avoid_sing_SSE = gmx_set1_pr(NBNXN_AVOID_SING_R2_INC);

    /* These kernels only support rvdw = rcoulomb */
    rc2_SSE   = gmx_set1_pr(ic->rvdw*ic->rvdw);

#ifdef CALC_ENERGIES
    sixthSSE    = gmx_set1_pr(1.0/6.0);
    twelvethSSE = gmx_set1_pr(1.0/12.0);

    sh_invrc6_SSE  = gmx_set1_pr(ic->sh_invrc6);
    sh_invrc12_SSE = gmx_set1_pr(ic->sh_invrc6*ic->sh_invrc6);
#endif

    mrc_3_SSE = gmx_set1_pr(-2*ic->k_rf);

#ifdef CALC_ENERGIES
    hrc_3_SSE = gmx_set1_pr(ic->k_rf);
    
    moh_rc_SSE = gmx_set1_pr(-ic->c_rf); 
#endif

#ifdef CALC_ENERGIES
    tmpsum = (real *)(((size_t)(tmpsum_array+7)) & (~((size_t)31)));
#endif
#ifdef CALC_SHIFTFORCES
    shf = (real *)(((size_t)(shf_array+7)) & (~((size_t)31)));
#endif

#ifdef FIX_LJ_C
    pvdw_c6  = (real *)(((size_t)(pvdw_array+3)) & (~((size_t)15)));
    pvdw_c12 = pvdw_c6 + UNROLLI*UNROLLJ;

    for(jp=0; jp<UNROLLJ; jp++)
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
    c6_SSE0            = gmx_load_pr(pvdw_c6 +0*UNROLLJ);
    c6_SSE1            = gmx_load_pr(pvdw_c6 +1*UNROLLJ);
    c6_SSE2            = gmx_load_pr(pvdw_c6 +2*UNROLLJ);
    c6_SSE3            = gmx_load_pr(pvdw_c6 +3*UNROLLJ);

    c12_SSE0           = gmx_load_pr(pvdw_c12+0*UNROLLJ);
    c12_SSE1           = gmx_load_pr(pvdw_c12+1*UNROLLJ);
    c12_SSE2           = gmx_load_pr(pvdw_c12+2*UNROLLJ);
    c12_SSE3           = gmx_load_pr(pvdw_c12+3*UNROLLJ);
#endif /* FIX_LJ_C */

#ifdef ENERGY_GROUPS
    egps_ishift  = nbat->neg_2log;
    egps_imask   = (1<<egps_ishift) - 1;
    egps_jshift  = 2*nbat->neg_2log;
    egps_jmask   = (1<<egps_jshift) - 1;
    egps_jstride = (UNROLLJ>>1)*UNROLLJ;
    /* Major division is over i-particles: divide nVS by 4 for i-stride */
    Vstride_i    = nbat->nenergrp*(1<<nbat->neg_2log)*egps_jstride;
#endif

    l_cj = nbl->cj;

    ninner = 0;
    for(n=0; n<nbl->nci; n++)
    {
        nbln = &nbl->ci[n];

        ish              = (nbln->shift & NBNXN_CI_SHIFT);
        ish3             = ish*3;
        cjind0           = nbln->cj_ind_start;      
        cjind1           = nbln->cj_ind_end;    
        /* Currently only works super-cells equal to sub-cells */
        ci               = nbln->ci;
        ci_sh            = (ish == CENTRAL ? ci : -1);

        shX_SSE = gmx_load1_pr(shiftvec+ish3);
        shY_SSE = gmx_load1_pr(shiftvec+ish3+1);
        shZ_SSE = gmx_load1_pr(shiftvec+ish3+2);

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

        half_LJ = (nbln->shift & NBNXN_CI_HALF_LJ(0));
        do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));

#ifdef ENERGY_GROUPS
        egps_i = nbat->energrp[ci];
        {
            int ia,egp_ia;

            for(ia=0; ia<4; ia++)
            {
                egp_ia = (egps_i >> (ia*egps_ishift)) & egps_imask;
                vvdwtp[ia] = Vvdw + egp_ia*Vstride_i;
                vctp[ia]   = Vc   + egp_ia*Vstride_i;
            }
        }
#endif
#if defined CALC_ENERGIES
#if UNROLLJ == 4
        if (do_coul && l_cj[nbln->cj_ind_start].cj == ci_sh)
#endif
#if UNROLLJ == 2
        if (do_coul && l_cj[nbln->cj_ind_start].cj == (ci_sh<<1))
#endif
#if UNROLLJ == 8
        if (do_coul && l_cj[nbln->cj_ind_start].cj == (ci_sh>>1))
#endif
        {
            int  ia;
            real Vc_sub_self;

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

            for(ia=0; ia<UNROLLI; ia++)
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
#endif

		/* Load i atom data */
        sciy             = scix + STRIDE;
        sciz             = sciy + STRIDE;
        ix_SSE0          = gmx_add_pr(gmx_load1_pr(x+scix)  ,shX_SSE);
        ix_SSE1          = gmx_add_pr(gmx_load1_pr(x+scix+1),shX_SSE);
        ix_SSE2          = gmx_add_pr(gmx_load1_pr(x+scix+2),shX_SSE);
        ix_SSE3          = gmx_add_pr(gmx_load1_pr(x+scix+3),shX_SSE);
        iy_SSE0          = gmx_add_pr(gmx_load1_pr(x+sciy)  ,shY_SSE);
        iy_SSE1          = gmx_add_pr(gmx_load1_pr(x+sciy+1),shY_SSE);
        iy_SSE2          = gmx_add_pr(gmx_load1_pr(x+sciy+2),shY_SSE);
        iy_SSE3          = gmx_add_pr(gmx_load1_pr(x+sciy+3),shY_SSE);
        iz_SSE0          = gmx_add_pr(gmx_load1_pr(x+sciz)  ,shZ_SSE);
        iz_SSE1          = gmx_add_pr(gmx_load1_pr(x+sciz+1),shZ_SSE);
        iz_SSE2          = gmx_add_pr(gmx_load1_pr(x+sciz+2),shZ_SSE);
        iz_SSE3          = gmx_add_pr(gmx_load1_pr(x+sciz+3),shZ_SSE);

        /* With half_LJ we currently always calculate Coulomb interactions */
        if (do_coul || half_LJ)
        {
            iq_SSE0      = gmx_set1_pr(facel*q[sci]);
            iq_SSE1      = gmx_set1_pr(facel*q[sci+1]);
            iq_SSE2      = gmx_set1_pr(facel*q[sci+2]);
            iq_SSE3      = gmx_set1_pr(facel*q[sci+3]);
        }

#ifdef LJ_COMB_LB
        hsig_i_SSE0      = gmx_load1_pr(ljc+sci2+0);
        hsig_i_SSE1      = gmx_load1_pr(ljc+sci2+1);
        hsig_i_SSE2      = gmx_load1_pr(ljc+sci2+2);
        hsig_i_SSE3      = gmx_load1_pr(ljc+sci2+3);
        seps_i_SSE0      = gmx_load1_pr(ljc+sci2+STRIDE+0);
        seps_i_SSE1      = gmx_load1_pr(ljc+sci2+STRIDE+1);
        seps_i_SSE2      = gmx_load1_pr(ljc+sci2+STRIDE+2);
        seps_i_SSE3      = gmx_load1_pr(ljc+sci2+STRIDE+3);
#else
#ifdef LJ_COMB_GEOM
        c6s_SSE0         = gmx_load1_pr(ljc+sci2+0);
        c6s_SSE1         = gmx_load1_pr(ljc+sci2+1);
        if (!half_LJ)
        {
            c6s_SSE2     = gmx_load1_pr(ljc+sci2+2);
            c6s_SSE3     = gmx_load1_pr(ljc+sci2+3);
        }
        c12s_SSE0        = gmx_load1_pr(ljc+sci2+STRIDE+0);
        c12s_SSE1        = gmx_load1_pr(ljc+sci2+STRIDE+1);
        if (!half_LJ)
        {
            c12s_SSE2    = gmx_load1_pr(ljc+sci2+STRIDE+2);
            c12s_SSE3    = gmx_load1_pr(ljc+sci2+STRIDE+3);
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

        /* Zero the potential energy for this list */
        VvdwtotSSE       = gmx_setzero_pr();
        vctotSSE         = gmx_setzero_pr();

        /* Clear i atom forces */
        fix_SSE0           = gmx_setzero_pr();
        fix_SSE1           = gmx_setzero_pr();
        fix_SSE2           = gmx_setzero_pr();
        fix_SSE3           = gmx_setzero_pr();
        fiy_SSE0           = gmx_setzero_pr();
        fiy_SSE1           = gmx_setzero_pr();
        fiy_SSE2           = gmx_setzero_pr();
        fiy_SSE3           = gmx_setzero_pr();
        fiz_SSE0           = gmx_setzero_pr();
        fiz_SSE1           = gmx_setzero_pr();
        fiz_SSE2           = gmx_setzero_pr();
        fiz_SSE3           = gmx_setzero_pr();

        cjind = cjind0;

        /* Currently all kernels use (at least half) LJ */
#define CALC_LJ
        if (half_LJ)
        {
#define CALC_COULOMB
#define HALF_LJ
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != SIMD_MASK_ALL)
            {
#include "nbnxn_kernel_x86_simd_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for(; (cjind<cjind1); cjind++)
            {
#include "nbnxn_kernel_x86_simd_inner.h"
            }
#undef HALF_LJ
#undef CALC_COULOMB
        }
        else if (do_coul)
        {
#define CALC_COULOMB
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != SIMD_MASK_ALL)
            {
#include "nbnxn_kernel_x86_simd_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for(; (cjind<cjind1); cjind++)
            {
#include "nbnxn_kernel_x86_simd_inner.h"
            }
#undef CALC_COULOMB
        }
        else
        {
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != SIMD_MASK_ALL)
            {
#include "nbnxn_kernel_x86_simd_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for(; (cjind<cjind1); cjind++)
            {
#include "nbnxn_kernel_x86_simd_inner.h"
            }
        }
#undef CALC_LJ
        ninner += cjind1 - cjind0;

        /* Add accumulated i-forces to the force array */
#if UNROLLJ >= 4
#ifndef GMX_DOUBLE
#define gmx_load_ps4  _mm_load_ps
#define gmx_store_ps4 _mm_store_ps
#define gmx_add_ps4   _mm_add_ps
#else
#define gmx_load_ps4  _mm256_load_pd
#define gmx_store_ps4 _mm256_store_pd
#define gmx_add_ps4   _mm256_add_pd
#endif
        GMX_MM_TRANSPOSE_SUM4_PR(fix_SSE0,fix_SSE1,fix_SSE2,fix_SSE3,fix_SSE);
        gmx_store_ps4(f+scix, gmx_add_ps4(fix_SSE, gmx_load_ps4(f+scix)));

        GMX_MM_TRANSPOSE_SUM4_PR(fiy_SSE0,fiy_SSE1,fiy_SSE2,fiy_SSE3,fiy_SSE);
        gmx_store_ps4(f+sciy, gmx_add_ps4(fiy_SSE, gmx_load_ps4(f+sciy)));

        GMX_MM_TRANSPOSE_SUM4_PR(fiz_SSE0,fiz_SSE1,fiz_SSE2,fiz_SSE3,fiz_SSE);
        gmx_store_ps4(f+sciz, gmx_add_ps4(fiz_SSE, gmx_load_ps4(f+sciz)));

#ifdef CALC_SHIFTFORCES
        gmx_store_ps4(shf,fix_SSE);
        fshift[ish3+0] += SUM_SIMD4(shf);
        gmx_store_ps4(shf,fiy_SSE);
        fshift[ish3+1] += SUM_SIMD4(shf);
        gmx_store_ps4(shf,fiz_SSE);
        fshift[ish3+2] += SUM_SIMD4(shf);
#endif
#else
        GMX_MM_TRANSPOSE_SUM2_PD(fix_SSE0,fix_SSE1,fix0_SSE);
        _mm_store_pd(f+scix, _mm_add_pd(fix0_SSE, _mm_load_pd(f+scix)));
        GMX_MM_TRANSPOSE_SUM2_PD(fix_SSE2,fix_SSE3,fix2_SSE);
        _mm_store_pd(f+scix+2, _mm_add_pd(fix2_SSE, _mm_load_pd(f+scix+2)));

        GMX_MM_TRANSPOSE_SUM2_PD(fiy_SSE0,fiy_SSE1,fiy0_SSE);
        _mm_store_pd(f+sciy, _mm_add_pd(fiy0_SSE, _mm_load_pd(f+sciy)));
        GMX_MM_TRANSPOSE_SUM2_PD(fiy_SSE2,fiy_SSE3,fiy2_SSE);
        _mm_store_pd(f+sciy+2, _mm_add_pd(fiy2_SSE, _mm_load_pd(f+sciy+2)));

        GMX_MM_TRANSPOSE_SUM2_PD(fiz_SSE0,fiz_SSE1,fiz0_SSE);
        _mm_store_pd(f+sciz, _mm_add_pd(fiz0_SSE, _mm_load_pd(f+sciz)));
        GMX_MM_TRANSPOSE_SUM2_PD(fiz_SSE2,fiz_SSE3,fiz2_SSE);
        _mm_store_pd(f+sciz+2, _mm_add_pd(fiz2_SSE, _mm_load_pd(f+sciz+2)));

#ifdef CALC_SHIFTFORCES
        _mm_store_pd(shf,_mm_add_pd(fix0_SSE,fix2_SSE));
        fshift[ish3+0] += shf[0] + shf[1];
        _mm_store_pd(shf,_mm_add_pd(fiy0_SSE,fiy2_SSE));
        fshift[ish3+1] += shf[0] + shf[1];
        _mm_store_pd(shf,_mm_add_pd(fiz0_SSE,fiz2_SSE));
        fshift[ish3+2] += shf[0] + shf[1];
#endif
#endif
		
#ifdef CALC_ENERGIES
        if (do_coul)
        {
            gmx_store_pr(tmpsum,vctotSSE);
            *Vc += SUM_SIMD(tmpsum);
        }
		
        gmx_store_pr(tmpsum,VvdwtotSSE);
        *Vvdw += SUM_SIMD(tmpsum);
#endif
		
		/* Outer loop uses 6 flops/iteration */
	}

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n",npair);
#endif
}

#undef gmx_load_ps4
#undef gmx_store_ps4
#undef gmx_store_ps4

#undef CALC_SHIFTFORCES

#undef UNROLLI   
#undef UNROLLJ   
#undef STRIDE
#undef TAB_FDV0
