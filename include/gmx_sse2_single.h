/*
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
#ifndef _gmx_sse2_single_h_
#define _gmx_sse2_single_h_

/* We require SSE2 now! */

#include <math.h> 


#include <xmmintrin.h> /* SSE */
#include <emmintrin.h> /* SSE2 */

#ifdef GMX_SSE3
#  include <pmmintrin.h> /* SSE3 */
#endif
#ifdef GMX_SSE4
#  include <smmintrin.h> /* SSE4.1 */
#endif

#include <stdio.h>

#include "types/simple.h"


/***************************************************
 *                                                 *
 * COMPILER RANT WARNING:                          *
 *                                                 *
 * Ideally, this header would be filled with       *
 * simple static inline functions. Unfortunately,  *
 * many vendors provide really braindead compilers *
 * that either cannot handle more than 1-2 SSE     *
 * function parameters, and some cannot handle     *
 * pointers to SSE __m128 datatypes as parameters  *
 * at all. Thus, for portability we have had to    *
 * implement all but the simplest routines as      *
 * macros instead...                               *
 *                                                 *
 ***************************************************/


/***************************************************
 *                                                 *
 *   Wrappers/replacements for some instructions   *
 *   not available in all SSE versions.            *
 *                                                 *
 ***************************************************/

#ifdef GMX_SSE4
#  define gmx_mm_extract_epi32(x, imm) _mm_extract_epi32(x,imm)
#else
#  define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))
#endif

/*
 * Some compilers require a cast to change the interpretation
 * of a register from FP to Int and vice versa, and not all of
 * the provide instructions to do this. Roll our own wrappers...
 */

#if (defined (_MSC_VER) || defined(__INTEL_COMPILER))
#  define gmx_mm_castsi128_ps(a) _mm_castsi128_ps(a)
#  define gmx_mm_castps_si128(a) _mm_castps_si128(a)
#  define gmx_mm_castps_ps128(a) (a)
#elif defined(__GNUC__)
#  define gmx_mm_castsi128_ps(a) ((__m128)(a))
#  define gmx_mm_castps_si128(a) ((__m128i)(a))
#  define gmx_mm_castps_ps128(a) ((__m128)(a))
#else
static __m128  gmx_mm_castsi128_ps(__m128i a) { return *(__m128 *) &a;  } 
static __m128i gmx_mm_castps_si128(__m128 a)  { return *(__m128i *) &a; } 
static __m128  gmx_mm_castps_ps128(__m128 a) { return *(__m128 *) &a;  } 
#endif



/* IO functions, just for debugging */

static void
printxmm(const char *s,__m128 xmm)
{
	float f[4];
	
	_mm_storeu_ps(f,xmm);
	printf("%s: %8.5g %8.5g %8.5g %8.5g\n",s,f[0],f[1],f[2],f[3]);	
}


static void
printxmmsum(const char *s,__m128 xmm)
{
	float f[4];
	
	_mm_storeu_ps(f,xmm);
	printf("%s (sum): %15.10g\n",s,f[0]+f[1]+f[2]+f[3]);	
}


static void
printxmmi(const char *s,__m128i xmmi)
{
    int i[4];
    
    _mm_storeu_si128((__m128i *)i,xmmi);
    printf("%10s: %2d %2d %2d %2d\n",s,i[0],i[1],i[2],i[3]);      
}


static int gmx_mm_check_and_reset_overflow(void)
{
    int MXCSR;
    int sse_overflow;

    MXCSR = _mm_getcsr();
    /* The overflow flag is bit 3 in the register */
    if (MXCSR & 0x0008)
    {
        sse_overflow = 1;
        /* Set the overflow flag to zero */
        MXCSR = MXCSR & 0xFFF7;
        _mm_setcsr(MXCSR);
    }
    else
    {
       sse_overflow = 0;
    }

    return sse_overflow;
}


/************************
 *                      *
 * Simple math routines *
 *                      *
 ************************/

static inline __m128
gmx_mm_invsqrt_ps(__m128 x)
{
    const __m128 half  = _mm_set_ps(0.5,0.5,0.5,0.5);
    const __m128 three = _mm_set_ps(3.0,3.0,3.0,3.0);
    
    __m128 lu = _mm_rsqrt_ps(x);
    
    return _mm_mul_ps(half,_mm_mul_ps(_mm_sub_ps(three,_mm_mul_ps(_mm_mul_ps(lu,lu),x)),lu));
}

static inline __m128
gmx_mm_sqrt_ps(__m128 x)
{
    __m128 mask;
    __m128 res;
    
    mask = _mm_cmpeq_ps(x,_mm_setzero_ps());
    res  = _mm_andnot_ps(mask,gmx_mm_invsqrt_ps(x));
    
    res  = _mm_mul_ps(x,res);
    
    return res;
}

static inline __m128
gmx_mm_inv_ps(__m128 x)
{
	const __m128 two = _mm_set_ps(2.0f,2.0f,2.0f,2.0f);
    
    __m128 lu = _mm_rcp_ps(x);
    
	return _mm_mul_ps(lu,_mm_sub_ps(two,_mm_mul_ps(lu,x)));
}


static inline __m128
gmx_mm_calc_rsq_ps(__m128 dx, __m128 dy, __m128 dz)
{
    return _mm_add_ps( _mm_add_ps( _mm_mul_ps(dx,dx), _mm_mul_ps(dy,dy) ), _mm_mul_ps(dz,dz) );
}

/* Normal sum of four xmm registers */
#define gmx_mm_sum4_ps(t0,t1,t2,t3)  _mm_add_ps(_mm_add_ps(t0,t1),_mm_add_ps(t2,t3))

static __m128
gmx_mm_log_ps(__m128 x)
{
    /* Same algorithm as cephes library */
	const __m128  expmask    = gmx_mm_castsi128_ps( _mm_set_epi32(0x7F800000, 0x7F800000, 0x7F800000, 0x7F800000) );
    const __m128i expbase_m1 = _mm_set1_epi32(127-1); /* We want non-IEEE format */
    const __m128  half       = _mm_set1_ps(0.5f);
    const __m128  one        = _mm_set1_ps(1.0f);
    const __m128  invsq2     = _mm_set1_ps(1.0f/sqrt(2.0f));
    const __m128  corr1      = _mm_set1_ps(-2.12194440e-4f);
    const __m128  corr2      = _mm_set1_ps(0.693359375f);
    
    const __m128 CA_1        = _mm_set1_ps(0.070376836292f);
    const __m128 CB_0        = _mm_set1_ps(1.6714950086782716f);
    const __m128 CB_1        = _mm_set1_ps(-2.452088066061482f);
    const __m128 CC_0        = _mm_set1_ps(1.5220770854701728f);
    const __m128 CC_1        = _mm_set1_ps(-1.3422238433233642f);
    const __m128 CD_0        = _mm_set1_ps(1.386218787509749f);
    const __m128 CD_1        = _mm_set1_ps(0.35075468953796346f);
    const __m128 CE_0        = _mm_set1_ps(1.3429983063133937f);
    const __m128 CE_1        = _mm_set1_ps(1.807420826584643f);
    
    __m128  fexp,fexp1;
    __m128i iexp;
    __m128  mask;
    __m128  x1,x2;
    __m128  y;
    __m128  pA,pB,pC,pD,pE,tB,tC,tD,tE;
    
    /* Separate x into exponent and mantissa, with a mantissa in the range [0.5..1[ (not IEEE754 standard!) */
    fexp  = _mm_and_ps(x,expmask);
    iexp  = gmx_mm_castps_si128(fexp);
    iexp  = _mm_srli_epi32(iexp,23);
    iexp  = _mm_sub_epi32(iexp,expbase_m1);
    
    x     = _mm_andnot_ps(expmask,x);
    x     = _mm_or_ps(x,one);
    x     = _mm_mul_ps(x,half);
    
    mask  = _mm_cmplt_ps(x,invsq2);
    
    x     = _mm_add_ps(x,_mm_and_ps(mask,x));
    x     = _mm_sub_ps(x,one);
    iexp  = _mm_add_epi32(iexp,gmx_mm_castps_si128(mask)); /* 0xFFFFFFFF = -1 as int */
    
    x2    = _mm_mul_ps(x,x);
    
    pA    = _mm_mul_ps(CA_1,x);
    pB    = _mm_mul_ps(CB_1,x);
    pC    = _mm_mul_ps(CC_1,x);
    pD    = _mm_mul_ps(CD_1,x);
    pE    = _mm_mul_ps(CE_1,x);
    tB    = _mm_add_ps(CB_0,x2);
    tC    = _mm_add_ps(CC_0,x2);
    tD    = _mm_add_ps(CD_0,x2);
    tE    = _mm_add_ps(CE_0,x2);
    pB    = _mm_add_ps(pB,tB);
    pC    = _mm_add_ps(pC,tC);
    pD    = _mm_add_ps(pD,tD);
    pE    = _mm_add_ps(pE,tE);
    
    pA    = _mm_mul_ps(pA,pB);
    pC    = _mm_mul_ps(pC,pD);
    pE    = _mm_mul_ps(pE,x2);
    pA    = _mm_mul_ps(pA,pC);
    y     = _mm_mul_ps(pA,pE);
    
    fexp  = _mm_cvtepi32_ps(iexp);
    y     = _mm_add_ps(y,_mm_mul_ps(fexp,corr1));
    
    y     = _mm_sub_ps(y, _mm_mul_ps(half,x2));
    x2    = _mm_add_ps(x,y);
    
    x2    = _mm_add_ps(x2,_mm_mul_ps(fexp,corr2));
    
    return x2;
}


/* 
 * Exponential function.
 * 
 * Exp(x) is calculate from the relation Exp(x)=2^(y), where y=log2(e)*x
 * Thus, the contents of this routine is mostly about calculating 2^y.
 *
 * This is done by separating y=z+w, where z=[y] is an integer. For technical reasons it is easiest
 * for us to round to the _nearest_ integer and have w in [-0.5,0.5] rather than always rounding down.
 * (It is not until SSE4 there was an efficient operation to do rounding towards -infinity).
 *
 * With this we get 2^y=2^z*2^w
 *
 * Since we have IEEE fp representation, we can easily calculate 2^z by adding the FP exponent bias
 * (127 in single), and shifting the integer to the exponent field of the FP number (23 bits up).
 *
 * The 2^w term is calculated from a (5,0)-th order (no denominator) Minimax polynomia on the interval
 * [-0.5,0.5]. The coefficiencts of this was derived in Mathematica using the command:
 *
 * MiniMaxApproximation[(2^x), {x, {-0.5, 0.5}, 5, 0}, WorkingPrecision -> 15]
 *
 * The lowest exponent we can represent in IEEE single-precision binary format is 2^-126; below that 
 * it will wrap around and lead to very large positive numbers. This corresponds to a lower bound
 * on the argument for exp(x) of roughly -87.33. For smaller arguments the return value will be 0.0.
 *
 * There appears to be a slight loss of precision for large arguments (~50), where the largest relative
 * error reaches ~3e-6. However, since the actual value for that argument is around 10^21, it might
 * not matter for typical single precision workloads. This is likely caused by the polynomial evaluation,
 * and the only way around would then be a table-based version, which I haven't managed to get the
 * same performance from.
 * 
 * The _average_ accuracy is 22.7 bits in the range [-10,10], and the worst roughly 1 bit worse.
 */
static __m128 
gmx_mm_exp_ps(__m128 x)
{
    const __m128  argscale = _mm_set1_ps(1.442695040888963f);
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m128  arglimit = _mm_set1_ps(-126.0f/1.442695040888963f);
    
    const __m128i expbase  = _mm_set1_epi32(127);
    const __m128  CA0      = _mm_set1_ps(0.00132764719920600f);
    const __m128  CB0      = _mm_set1_ps(3.17196359322f);
    const __m128  CC0      = _mm_set1_ps(20.36135752425f);
    const __m128  CC1      = _mm_set1_ps(-0.681627790451f);
    const __m128  CD0      = _mm_set1_ps(11.66225206128f);
    const __m128  CD1      = _mm_set1_ps(4.79739947827f);
    
    
    __m128  valuemask;
    __m128i iexppart;
    __m128  fexppart;
    __m128  intpart;
    __m128  z,z2;
    __m128  factB,factC,factD;
    
	z         = _mm_mul_ps(x,argscale);     
    iexppart  = _mm_cvtps_epi32(z);
#if GMX_SSE4
    /* This reduces latency and speeds up the code by roughly 5% when supported */
    intpart   = _mm_round_ps(z,0);
#else
    intpart   = _mm_cvtepi32_ps(iexppart);
#endif
    iexppart  = _mm_slli_epi32(_mm_add_epi32(iexppart,expbase),23);    
    valuemask = _mm_cmpgt_ps(x,arglimit);
    
	z         = _mm_sub_ps(z,intpart);
    z2        = _mm_mul_ps(z,z);  
    
    fexppart  = _mm_and_ps(valuemask,gmx_mm_castsi128_ps(iexppart));
    
    /* Since SSE floating-point has relatively high latency it is faster to do 
     * factorized polynomial summation with independent terms than using alternating add/multiply, i.e.
     * p(z) = A0 * (B0 + z) * (C0 + C1*z + z^2) * (D0 + D1*z + z^2) 
     */
    factB     = _mm_add_ps(CB0,z);
    factC     = _mm_add_ps(CC0,_mm_mul_ps(CC1,z) );
    factC     = _mm_add_ps(factC,z2);
    factD     = _mm_add_ps(CD0,_mm_mul_ps(CD1,z) );
    factD     = _mm_add_ps(factD,z2);
    
    z         = _mm_mul_ps(CA0,fexppart);
    factB     = _mm_mul_ps(factB,factC);
    z         = _mm_mul_ps(z,factD);
    z         = _mm_mul_ps(z,factB);
    
	/* Currently uses 22 actual (real, not including casts) SSE instructions */
	return  z;
}



static int 
gmx_mm_sincos_ps(__m128 x,
                 __m128 *sinval,
                 __m128 *cosval)
{
    const __m128 _sincosf_two_over_pi = _mm_set_ps(2.0/M_PI,2.0/M_PI,2.0/M_PI,2.0/M_PI); 
    const __m128 _sincosf_half        = _mm_set_ps(0.5,0.5,0.5,0.5);
    const __m128 _sincosf_one         = _mm_set_ps(1.0,1.0,1.0,1.0);
    
    const __m128i _sincosf_izero      = _mm_set1_epi32(0);                                                   
    const __m128i _sincosf_ione       = _mm_set1_epi32(1);                                                   
    const __m128i _sincosf_itwo       = _mm_set1_epi32(2);                                                   
    const __m128i _sincosf_ithree     = _mm_set1_epi32(3);                                                   
    
    const __m128 _sincosf_kc1 = _mm_set_ps(1.57079625129,1.57079625129,1.57079625129,1.57079625129);                   
    const __m128 _sincosf_kc2 = _mm_set_ps(7.54978995489e-8,7.54978995489e-8,7.54978995489e-8,7.54978995489e-8);       
    const __m128 _sincosf_cc0 = _mm_set_ps(-0.0013602249,-0.0013602249,-0.0013602249,-0.0013602249);                   
    const __m128 _sincosf_cc1 = _mm_set_ps(0.0416566950,0.0416566950,0.0416566950,0.0416566950);                       
    const __m128 _sincosf_cc2 = _mm_set_ps(-0.4999990225,-0.4999990225,-0.4999990225,-0.4999990225);                   
    const __m128 _sincosf_sc0 = _mm_set_ps(-0.0001950727,-0.0001950727,-0.0001950727,-0.0001950727);                   
    const __m128 _sincosf_sc1 = _mm_set_ps(0.0083320758,0.0083320758,0.0083320758,0.0083320758);                       
    const __m128 _sincosf_sc2 = _mm_set_ps(-0.1666665247,-0.1666665247,-0.1666665247,-0.1666665247);                   
    
    __m128 _sincosf_signbit           = gmx_mm_castsi128_ps( _mm_set1_epi32(0x80000000) );                   
    __m128 _sincosf_tiny              = gmx_mm_castsi128_ps( _mm_set1_epi32(0x3e400000) );                   
    
    __m128 _sincosf_xl;                                                                                      
    __m128 _sincosf_xl2;                                                                                     
    __m128 _sincosf_xl3;                                                                                     
    __m128 _sincosf_qf;                                                                                      
    __m128 _sincosf_absxl;                                                                                   
    __m128 _sincosf_p1;                                                                                      
    __m128 _sincosf_cx;                                                                                      
    __m128 _sincosf_sx;                                                                                      
    __m128 _sincosf_ts;                                                                                      
    __m128 _sincosf_tc;                                                                                      
    __m128 _sincosf_tsn;                                                                                     
    __m128 _sincosf_tcn;                                                                                     
    __m128i _sincosf_q;                                                                                      
    __m128i _sincosf_offsetSin;                                                                              
    __m128i _sincosf_offsetCos;                                                                              
    __m128 _sincosf_sinMask;                                                                                 
    __m128 _sincosf_cosMask;                                                                                 
    __m128 _sincosf_isTiny;                                                                                  
    __m128 _sincosf_ct0;                                                                                     
    __m128 _sincosf_ct1;                                                                                     
    __m128 _sincosf_ct2;                                                                                     
    __m128 _sincosf_st1;                                                                                     
    __m128 _sincosf_st2;                                                                                     
    
    _sincosf_xl        = _mm_mul_ps(x,_sincosf_two_over_pi);                                                 
    
    _sincosf_xl        = _mm_add_ps(_sincosf_xl,_mm_or_ps(_mm_and_ps(_sincosf_xl,_sincosf_signbit),_sincosf_half)); 
    
    _sincosf_q         = _mm_cvttps_epi32(_sincosf_xl);                                                      
    _sincosf_qf        = _mm_cvtepi32_ps(_sincosf_q);                                                        
    
    _sincosf_offsetSin   = _mm_and_si128(_sincosf_q,_sincosf_ithree);                                        
    _sincosf_offsetCos   = _mm_add_epi32(_sincosf_offsetSin,_sincosf_ione);                                  
    
    _sincosf_p1 = _mm_mul_ps(_sincosf_qf,_sincosf_kc1);                                                      
    _sincosf_xl = _mm_mul_ps(_sincosf_qf,_sincosf_kc2);                                                      
    _sincosf_p1 = _mm_sub_ps(x,_sincosf_p1);                                                                 
    _sincosf_xl = _mm_sub_ps(_sincosf_p1,_sincosf_xl);                                                       
    
    _sincosf_absxl  = _mm_andnot_ps(_sincosf_signbit,_sincosf_xl);                                           
    _sincosf_isTiny = _mm_cmpgt_ps(_sincosf_tiny,_sincosf_absxl);                                            
    
    _sincosf_xl2    = _mm_mul_ps(_sincosf_xl,_sincosf_xl);                                                   
    _sincosf_xl3    = _mm_mul_ps(_sincosf_xl2,_sincosf_xl);                                                  
    
    _sincosf_ct1    = _mm_mul_ps(_sincosf_cc0,_sincosf_xl2);                                                 
    _sincosf_ct1    = _mm_add_ps(_sincosf_ct1,_sincosf_cc1);                                                 
    _sincosf_st1    = _mm_mul_ps(_sincosf_sc0,_sincosf_xl2);                                                 
    _sincosf_st1    = _mm_add_ps(_sincosf_st1,_sincosf_sc1);                                                 
    _sincosf_ct2    = _mm_mul_ps(_sincosf_ct1,_sincosf_xl2);                                                 
    _sincosf_ct2    = _mm_add_ps(_sincosf_ct2,_sincosf_cc2);                                                 
    _sincosf_st2    = _mm_mul_ps(_sincosf_st1,_sincosf_xl2);                                                 
    _sincosf_st2    = _mm_add_ps(_sincosf_st2,_sincosf_sc2);                                                 
    
    _sincosf_cx     = _mm_mul_ps(_sincosf_ct2,_sincosf_xl2);                                                 
    _sincosf_cx     = _mm_add_ps(_sincosf_cx,_sincosf_one);                                                  
    
    _sincosf_sx     = _mm_mul_ps(_sincosf_st2,_sincosf_xl3);                                                 
    _sincosf_sx     = _mm_add_ps(_sincosf_sx,_sincosf_xl);                                                   
    
    _sincosf_sinMask = gmx_mm_castsi128_ps( _mm_cmpeq_epi32( _mm_and_si128(_sincosf_offsetSin,_sincosf_ione), _sincosf_izero) ); 
    _sincosf_cosMask = gmx_mm_castsi128_ps( _mm_cmpeq_epi32( _mm_and_si128(_sincosf_offsetCos,_sincosf_ione), _sincosf_izero) ); 
    
    _sincosf_ts     = _mm_or_ps( _mm_and_ps(_sincosf_sinMask,_sincosf_sx) , _mm_andnot_ps(_sincosf_sinMask,_sincosf_cx) ); 
    _sincosf_tc     = _mm_or_ps( _mm_and_ps(_sincosf_cosMask,_sincosf_sx) , _mm_andnot_ps(_sincosf_cosMask,_sincosf_cx) ); 
    
    _sincosf_sinMask = gmx_mm_castsi128_ps(  _mm_cmpeq_epi32( _mm_and_si128(_sincosf_offsetSin,_sincosf_itwo), _sincosf_izero) );
    _sincosf_tsn    = _mm_xor_ps(_sincosf_signbit,_sincosf_ts);                                              
    _sincosf_ts     = _mm_or_ps( _mm_and_ps(_sincosf_sinMask,_sincosf_ts) , _mm_andnot_ps(_sincosf_sinMask,_sincosf_tsn) ); 
    
    _sincosf_cosMask = gmx_mm_castsi128_ps(  _mm_cmpeq_epi32( _mm_and_si128(_sincosf_offsetCos,_sincosf_itwo), _sincosf_izero) ); 
    _sincosf_tcn    = _mm_xor_ps(_sincosf_signbit,_sincosf_tc);                                              
    _sincosf_tc     = _mm_or_ps( _mm_and_ps(_sincosf_cosMask,_sincosf_tc) , _mm_andnot_ps(_sincosf_cosMask,_sincosf_tcn) ); 
    
    *sinval = _sincosf_ts;                                                                                    
    *cosval = _sincosf_tc;      
    
    return 0;
}

static __m128
gmx_mm_tan_ps(__m128 x)
{
    __m128 sinval,cosval;
    __m128 tanval;
    
    gmx_mm_sincos_ps(x,&sinval,&cosval);
    
    tanval = _mm_mul_ps(sinval,gmx_mm_inv_ps(cosval));
    
    return tanval;
}


static __m128
gmx_mm_asin_ps(__m128 x)
{
    /* Same algorithm as cephes library */
    const __m128 signmask  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x7FFFFFFF) );
    const __m128 limitlow  = _mm_set1_ps(1e-4f);
    const __m128 half      = _mm_set1_ps(0.5f);
    const __m128 one       = _mm_set1_ps(1.0f);
    const __m128 halfpi    = _mm_set1_ps(M_PI/2.0f);
    
    const __m128 CC5        = _mm_set1_ps(4.2163199048E-2f);
    const __m128 CC4        = _mm_set1_ps(2.4181311049E-2f);
    const __m128 CC3        = _mm_set1_ps(4.5470025998E-2f);
    const __m128 CC2        = _mm_set1_ps(7.4953002686E-2f);
    const __m128 CC1        = _mm_set1_ps(1.6666752422E-1f);
    
    __m128 sign;
    __m128 mask;
    __m128 xabs;
    __m128 z,z1,z2,q,q1,q2;
    __m128 pA,pB;
    
    sign  = _mm_andnot_ps(signmask,x);
    xabs  = _mm_and_ps(x,signmask);
    
    mask  = _mm_cmpgt_ps(xabs,half);
    
    z1    = _mm_mul_ps(half, _mm_sub_ps(one,xabs));
    q1    = _mm_mul_ps(z1,gmx_mm_invsqrt_ps(z1));
    q1    = _mm_andnot_ps(_mm_cmpeq_ps(xabs,one),q1);
    
    q2    = xabs;
    z2    = _mm_mul_ps(q2,q2);
    
    z     = _mm_or_ps( _mm_and_ps(mask,z1) , _mm_andnot_ps(mask,z2) );
    q     = _mm_or_ps( _mm_and_ps(mask,q1) , _mm_andnot_ps(mask,q2) );
    
    z2    = _mm_mul_ps(z,z);
    
    pA    = _mm_mul_ps(CC5,z2);
    pB    = _mm_mul_ps(CC4,z2);
    
    pA    = _mm_add_ps(pA,CC3);
    pB    = _mm_add_ps(pB,CC2);
    
    pA    = _mm_mul_ps(pA,z2);
    pB    = _mm_mul_ps(pB,z2);
    
    pA    = _mm_add_ps(pA,CC1);
    pA    = _mm_mul_ps(pA,z);
    
    z     = _mm_add_ps(pA,pB);
    z     = _mm_mul_ps(z,q);
    z     = _mm_add_ps(z,q);
    
    q2    = _mm_sub_ps(halfpi,z);
    q2    = _mm_sub_ps(q2,z);
    
    z     = _mm_or_ps( _mm_and_ps(mask,q2) , _mm_andnot_ps(mask,z) );
    
    mask  = _mm_cmpgt_ps(xabs,limitlow);
    z     = _mm_or_ps( _mm_and_ps(mask,z) , _mm_andnot_ps(mask,xabs) );
    
    z = _mm_xor_ps(z,sign);
    
    return z;
}


static __m128
gmx_mm_acos_ps(__m128 x)
{
    const __m128 signmask  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x7FFFFFFF) );
    const __m128 one_ps    = _mm_set1_ps(1.0f);
    const __m128 half_ps   = _mm_set1_ps(0.5f);
    const __m128 pi_ps     = _mm_set1_ps(M_PI);
    const __m128 halfpi_ps = _mm_set1_ps(M_PI/2.0f);
    
    __m128 mask1;
    __m128 mask2;
    __m128 xabs;
    __m128 z,z1,z2,z3;
    
    xabs  = _mm_and_ps(x,signmask);    
    mask1 = _mm_cmpgt_ps(xabs,half_ps);
    mask2 = _mm_cmpgt_ps(x,_mm_setzero_ps());
    
    z     = _mm_mul_ps(half_ps,_mm_sub_ps(one_ps,xabs));
    z     = _mm_mul_ps(z,gmx_mm_invsqrt_ps(z));
    z     = _mm_andnot_ps(_mm_cmpeq_ps(xabs,one_ps),z);
    
    z     = _mm_or_ps( _mm_and_ps(mask1,z) , _mm_andnot_ps(mask1,x) );
    z     = gmx_mm_asin_ps(z);
    
    z2    = _mm_add_ps(z,z);
    z1    = _mm_sub_ps(pi_ps,z2);
    z3    = _mm_sub_ps(halfpi_ps,z);    
    
    z     = _mm_or_ps( _mm_and_ps(mask2,z2) , _mm_andnot_ps(mask2,z1) );
    z     = _mm_or_ps( _mm_and_ps(mask1,z) , _mm_andnot_ps(mask1,z3) );
    
    return z;
}


static __m128
gmx_mm_atan_ps(__m128 x)
{
    /* Same algorithm as cephes library */
    const __m128 signmask  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x7FFFFFFF) );
    const __m128 limit1    = _mm_set1_ps(0.414213562373095f);
    const __m128 limit2    = _mm_set1_ps(2.414213562373095f);
    const __m128 quarterpi = _mm_set1_ps(0.785398163397448f);
    const __m128 halfpi    = _mm_set1_ps(1.570796326794896f);
    const __m128 mone      = _mm_set1_ps(-1.0f);
    const __m128 CC3       = _mm_set1_ps(-3.33329491539E-1f);
    const __m128 CC5       = _mm_set1_ps(1.99777106478E-1f);
    const __m128 CC7       = _mm_set1_ps(-1.38776856032E-1);
    const __m128 CC9       = _mm_set1_ps(8.05374449538e-2f);
    
    __m128 sign;
    __m128 mask1,mask2;
    __m128 y,z1,z2;
    __m128 x2,x4;
    __m128 sum1,sum2;
    
    sign  = _mm_andnot_ps(signmask,x);
    x     = _mm_and_ps(x,signmask);
    
    mask1 = _mm_cmpgt_ps(x,limit1);
    mask2 = _mm_cmpgt_ps(x,limit2);
    
    z1    = _mm_mul_ps(_mm_add_ps(x,mone),gmx_mm_inv_ps(_mm_sub_ps(x,mone)));
    z2    = _mm_mul_ps(mone,gmx_mm_inv_ps(x));
    
    y     = _mm_and_ps(mask1,quarterpi);
    y     = _mm_or_ps( _mm_and_ps(mask2,halfpi) , _mm_andnot_ps(mask2,y) );
    
    x     = _mm_or_ps( _mm_and_ps(mask1,z1) , _mm_andnot_ps(mask1,x) );
    x     = _mm_or_ps( _mm_and_ps(mask2,z2) , _mm_andnot_ps(mask2,x) );
    
    x2    = _mm_mul_ps(x,x);
    x4    = _mm_mul_ps(x2,x2);
    
    sum1  = _mm_mul_ps(CC9,x4);    
    sum2  = _mm_mul_ps(CC7,x4);    
    sum1  = _mm_add_ps(sum1,CC5);
    sum2  = _mm_add_ps(sum2,CC3); 
    sum1  = _mm_mul_ps(sum1,x4);
    sum2  = _mm_mul_ps(sum2,x2);
    
    sum1  = _mm_add_ps(sum1,sum2);
    sum1  = _mm_sub_ps(sum1,mone);
    sum1  = _mm_mul_ps(sum1,x);
    y     = _mm_add_ps(y,sum1);
    
    y     = _mm_xor_ps(y,sign);
    
    return y;
}


static __m128
gmx_mm_atan2_ps(__m128 y, __m128 x)
{
    const __m128 pi          = _mm_set1_ps(M_PI);
    const __m128 minuspi     = _mm_set1_ps(-M_PI);
    const __m128 halfpi      = _mm_set1_ps(M_PI/2.0);
    const __m128 minushalfpi = _mm_set1_ps(-M_PI/2.0);
    
    __m128 z,z1,z3,z4;
    __m128 w;
    __m128 maskx_lt,maskx_eq;
    __m128 masky_lt,masky_eq;
    __m128 mask1,mask2,mask3,mask4,maskall;
    
    maskx_lt  = _mm_cmplt_ps(x,_mm_setzero_ps());
    masky_lt  = _mm_cmplt_ps(y,_mm_setzero_ps());
    maskx_eq  = _mm_cmpeq_ps(x,_mm_setzero_ps());
    masky_eq  = _mm_cmpeq_ps(y,_mm_setzero_ps());
    
    z         = _mm_mul_ps(y,gmx_mm_inv_ps(x));
    z         = gmx_mm_atan_ps(z);
    
    mask1     = _mm_and_ps(maskx_eq,masky_lt);
    mask2     = _mm_andnot_ps(maskx_lt,masky_eq);
    mask3     = _mm_andnot_ps( _mm_or_ps(masky_lt,masky_eq) , maskx_eq);
    mask4     = _mm_and_ps(masky_eq,maskx_lt);
    
    maskall   = _mm_or_ps( _mm_or_ps(mask1,mask2), _mm_or_ps(mask3,mask4) );
    
    z         = _mm_andnot_ps(maskall,z);
    z1        = _mm_and_ps(mask1,minushalfpi);
    z3        = _mm_and_ps(mask3,halfpi);
    z4        = _mm_and_ps(mask4,pi);
    
    z         = _mm_or_ps( _mm_or_ps(z,z1), _mm_or_ps(z3,z4) );
    
    mask1     = _mm_andnot_ps(masky_lt,maskx_lt);
    mask2     = _mm_and_ps(maskx_lt,masky_lt);
    
    w         = _mm_or_ps( _mm_and_ps(mask1,pi), _mm_and_ps(mask2,minuspi) );    
    w         = _mm_andnot_ps(maskall,w);
    
    z         = _mm_add_ps(z,w);
    
    return z;
}

/* Load a single value from 1-4 places, merge into xmm register */

#define GMX_MM_LOAD_4VALUES_PS(ptr1,ptr2,ptr3,ptr4,xmm1) \
{                                                         \
    __m128 _txmm2,_txmm3,_txmm4;                          \
    xmm1           = _mm_load_ss(ptr1);                   \
    _txmm2         = _mm_load_ss(ptr2);                   \
    _txmm3         = _mm_load_ss(ptr3);                   \
    _txmm4         = _mm_load_ss(ptr4);                   \
    xmm1           = _mm_unpacklo_ps(xmm1,_txmm3);        \
    _txmm2         = _mm_unpacklo_ps(_txmm2,_txmm4);      \
    xmm1           = _mm_unpacklo_ps(xmm1,_txmm2);        \
}


#define GMX_MM_LOAD_3VALUES_PS(ptr1,ptr2,ptr3,xmm1) \
{                                                    \
    __m128 _txmm2,_txmm3;                            \
    xmm1           = _mm_load_ss(ptr1);              \
    _txmm2         = _mm_load_ss(ptr2);              \
    _txmm3         = _mm_load_ss(ptr3);              \
    xmm1           = _mm_unpacklo_ps(xmm1,_txmm3);   \
    xmm1           = _mm_unpacklo_ps(xmm1,_txmm2);   \
}


#define GMX_MM_LOAD_2VALUES_PS(ptr1,ptr2,xmm1)    \
{                                                  \
    __m128 _txmm2;                                 \
    xmm1           = _mm_load_ss(ptr1);            \
    _txmm2         = _mm_load_ss(ptr2);            \
    xmm1           = _mm_unpacklo_ps(xmm1,_txmm2); \
}


#define GMX_MM_LOAD_1VALUE_PS(ptr1,xmm1) \
{                                         \
      xmm1           = _mm_load_ss(ptr1); \
}

/* Store data in an xmm register into 1-4 different places */
#define GMX_MM_STORE_4VALUES_PS(ptr1,ptr2,ptr3,ptr4,xmm1)             \
{                                                                      \
    __m128 _txmm2,_txmm3,_txmm4;                                       \
    _txmm3       = _mm_movehl_ps(_mm_setzero_ps(),xmm1);               \
    _txmm2       = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(1,1,1,1));     \
    _txmm4       = _mm_shuffle_ps(_txmm3,_txmm3,_MM_SHUFFLE(1,1,1,1)); \
    _mm_store_ss(ptr1,xmm1);                                           \
    _mm_store_ss(ptr2,_txmm2);                                         \
    _mm_store_ss(ptr3,_txmm3);                                         \
    _mm_store_ss(ptr4,_txmm4);                                         \
}


#define GMX_MM_STORE_3VALUES_PS(ptr1,ptr2,ptr3,xmm1)              \
{                                                                  \
    __m128 _txmm2,_txmm3;                                          \
    _txmm3       = _mm_movehl_ps(_mm_setzero_ps(),xmm1);           \
    _txmm2       = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(1,1,1,1)); \
    _mm_store_ss(ptr1,xmm1);                                       \
    _mm_store_ss(ptr2,_txmm2);                                     \
    _mm_store_ss(ptr3,_txmm3);                                     \
}


#define GMX_MM_STORE_2VALUES_PS(ptr1,ptr2,xmm1)                   \
{                                                                  \
    __m128 _txmm2;                                                 \
    _txmm2       = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(1,1,1,1)); \
    _mm_store_ss(ptr1,xmm1);                                       \
    _mm_store_ss(ptr2,_txmm2);                                     \
}


#define GMX_MM_STORE_1VALUE_PS(ptr1,xmm1) \
{                                          \
    _mm_store_ss(ptr1,xmm1);               \
}


/* Similar to store, but increments value in memory */
#define GMX_MM_INCREMENT_8VALUES_PS(ptr1,ptr2,ptr3,ptr4,ptr5,ptr6,ptr7,ptr8,xmm1,xmm2)    \
{                                                                  \
    __m128 _tincr1,_tincr2;                                        \
    GMX_MM_LOAD_4VALUES_PS(ptr1,ptr2,ptr3,ptr4,_tincr1);          \
    GMX_MM_LOAD_4VALUES_PS(ptr5,ptr6,ptr7,ptr8,_tincr2);          \
    _tincr1 = _mm_add_ps(_tincr1,xmm1);                            \
    _tincr2 = _mm_add_ps(_tincr2,xmm2);                            \
    GMX_MM_STORE_4VALUES_PS(ptr1,ptr2,ptr3,ptr4,_tincr1);         \
    GMX_MM_STORE_4VALUES_PS(ptr5,ptr6,ptr7,ptr8,_tincr2);         \
}

#define GMX_MM_INCREMENT_4VALUES_PS(ptr1,ptr2,ptr3,ptr4,xmm1)    \
{                                                                 \
    __m128 _tincr;                                                \
    GMX_MM_LOAD_4VALUES_PS(ptr1,ptr2,ptr3,ptr4,_tincr);          \
    _tincr = _mm_add_ps(_tincr,xmm1);                             \
    GMX_MM_STORE_4VALUES_PS(ptr1,ptr2,ptr3,ptr4,_tincr);         \
}

#define GMX_MM_INCREMENT_3VALUES_PS(ptr1,ptr2,ptr3,xmm1)         \
{                                                                 \
    __m128 _tincr;                                                \
    GMX_MM_LOAD_3VALUES_PS(ptr1,ptr2,ptr3,_tincr);               \
    _tincr = _mm_add_ps(_tincr,xmm1);                             \
    GMX_MM_STORE_3VALUES_PS(ptr1,ptr2,ptr3,_tincr);              \
}

#define GMX_MM_INCREMENT_2VALUES_PS(ptr1,ptr2,xmm1)         \
{                                                            \
    __m128 _tincr;                                           \
    GMX_MM_LOAD_2VALUES_PS(ptr1,ptr2,_tincr);               \
    _tincr = _mm_add_ps(_tincr,xmm1);                        \
    GMX_MM_STORE_2VALUES_PS(ptr1,ptr2,_tincr);              \
}

#define GMX_MM_INCREMENT_1VALUE_PS(ptr1,xmm1)         \
{                                                      \
    __m128 _tincr;                                     \
    GMX_MM_LOAD_1VALUE_PS(ptr1,_tincr);               \
    _tincr = _mm_add_ss(_tincr,xmm1);                  \
    GMX_MM_STORE_1VALUE_PS(ptr1,_tincr);              \
}



/* Routines to load pairs from 1-4 places, put in two separate xmm registers. Useful to load LJ parameters! */
#define GMX_MM_LOAD_4PAIRS_PS(ptr1,ptr2,ptr3,ptr4,c6,c12)     \
{                                                             \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4;                           \
    _tmp1  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1));  \
    _tmp2  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2));  \
    _tmp3  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3));  \
    _tmp4  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr4));  \
    _tmp1  = _mm_unpacklo_ps(_tmp1,_tmp3);                    \
    _tmp2  = _mm_unpacklo_ps(_tmp2,_tmp4);                    \
    c6     = _mm_unpacklo_ps(_tmp1,_tmp2);                    \
    c12    = _mm_unpackhi_ps(_tmp1,_tmp2);                    \
}

#define GMX_MM_LOAD_3PAIRS_PS(ptr1,ptr2,ptr3,c6,c12)        \
{                                                            \
    __m128 _tmp1,_tmp2,_tmp3;                                \
    _tmp1  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1)); \
    _tmp2  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2)); \
    _tmp3  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3)); \
    _tmp1  = _mm_unpacklo_ps(_tmp1,_tmp3);                   \
    _tmp2  = _mm_unpacklo_ps(_tmp2,_mm_setzero_ps());        \
    c6     = _mm_unpacklo_ps(_tmp1,_tmp2);                   \
    c12    = _mm_unpackhi_ps(_tmp1,_tmp2);                   \
}


#define GMX_MM_LOAD_2PAIRS_PS(ptr1,ptr2,c6,c12)             \
{                                                            \
    __m128 _tmp1,_tmp2;                                      \
    _tmp1  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1)); \
    _tmp2  = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2)); \
    c6     = _mm_unpacklo_ps(_tmp1,_tmp2);                   \
    c12    = _mm_movehl_ps(c12,c6);                          \
}

#define GMX_MM_LOAD_1PAIR_PS(ptr1,c6,c12)                   \
{                                                            \
    c6     = _mm_load_ss(ptr1);                              \
    c12    = _mm_load_ss(ptr1+1);                            \
}


/* Routines to load 1-4 rvecs from 1-4 places. 
 * We mainly use these to load coordinates. The extra routines
 * are very efficient for the water-water loops, since we e.g.
 * know that a TIP4p water has 4 atoms, so we should load 12 floats+shuffle.
 */
#define GMX_MM_LOAD_1RVEC_1POINTER_PS(ptr1,jx1,jy1,jz1) {             \
	 jx1            = _mm_load_ss(ptr1);                                \
     jy1            = _mm_load_ss((ptr1)+1);                            \
     jz1            = _mm_load_ss((ptr1)+2);                            \
}

#define GMX_MM_LOAD_2RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) {      \
	 jx1            = _mm_load_ss(ptr1);                                      \
     jy1            = _mm_load_ss((ptr1)+1);                                  \
     jz1            = _mm_load_ss((ptr1)+2);                                  \
	 jx2            = _mm_load_ss((ptr1)+3);                                  \
     jy2            = _mm_load_ss((ptr1)+4);                                  \
     jz2            = _mm_load_ss((ptr1)+5);                                  \
}


#define GMX_MM_LOAD_3RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
	 jx1            = _mm_load_ss(ptr1);                                    \
     jy1            = _mm_load_ss((ptr1)+1);                                \
     jz1            = _mm_load_ss((ptr1)+2);                                \
	 jx2            = _mm_load_ss((ptr1)+3);                                \
     jy2            = _mm_load_ss((ptr1)+4);                                \
     jz2            = _mm_load_ss((ptr1)+5);                                \
	 jx3            = _mm_load_ss((ptr1)+6);                                \
     jy3            = _mm_load_ss((ptr1)+7);                                \
     jz3            = _mm_load_ss((ptr1)+8);                                \
}


#define GMX_MM_LOAD_4RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
	 jx1            = _mm_load_ss(ptr1);                                    \
     jy1            = _mm_load_ss((ptr1)+1);                                \
     jz1            = _mm_load_ss((ptr1)+2);                                \
	 jx2            = _mm_load_ss((ptr1)+3);                                \
     jy2            = _mm_load_ss((ptr1)+4);                                \
     jz2            = _mm_load_ss((ptr1)+5);                                \
	 jx3            = _mm_load_ss((ptr1)+6);                                \
     jy3            = _mm_load_ss((ptr1)+7);                                \
     jz3            = _mm_load_ss((ptr1)+8);                                \
	 jx4            = _mm_load_ss((ptr1)+9);                                \
     jy4            = _mm_load_ss((ptr1)+10);                               \
     jz4            = _mm_load_ss((ptr1)+11);                               \
}


#define GMX_MM_LOAD_1RVEC_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1) {  \
      __m128 _tmp1,_tmp2;                                           \
      _tmp1           = _mm_load_ss(ptr1);                          \
	  _tmp2           = _mm_load_ss(ptr2);                          \
      _tmp1           = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));      \
      _tmp2           = _mm_loadh_pi(_tmp2,(__m64 *)(ptr2+1));      \
      jx1             = _mm_unpacklo_ps(_tmp1,_tmp2);               \
      jy1             = _mm_unpackhi_ps(_tmp1,_tmp2);               \
	  jx1             = _mm_unpacklo_ps(_tmp1,_tmp2);               \
      jz1             = _mm_movehl_ps(jz1,jy1);                     \
}

#define GMX_MM_LOAD_2RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1, _tmp2;                                                      \
	 _tmp1          = _mm_loadu_ps(ptr1);                                      \
     jy1            = _mm_loadu_ps(ptr2);                                      \
     jy2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));        \
     _tmp2          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2+4));        \
     jx1            = _mm_unpacklo_ps(_tmp1,jy1);                              \
     jz1            = _mm_unpackhi_ps(_tmp1,jy1);                              \
     jy2            = _mm_unpacklo_ps(jy2,_tmp2);                              \
     jy1            = _mm_movehl_ps(jx1,jx1);                                  \
     jx2            = _mm_movehl_ps(jz1,jz1);                                  \
     jz2            = _mm_movehl_ps(jy2,jy2);                                  \
}


#define GMX_MM_LOAD_3RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1, _tmp2, _tmp3;                                                           \
	 _tmp1          = _mm_loadu_ps(ptr1);                                                  \
     jy1            = _mm_loadu_ps(ptr2);                                                  \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                \
     jz2            = _mm_loadu_ps(ptr2+4);                                                \
     jz3            = _mm_load_ss(ptr1+8);                                                 \
     _tmp3          = _mm_load_ss(ptr2+8);                                                 \
     jx1            = _mm_unpacklo_ps(_tmp1,jy1);                                          \
     jz1            = _mm_unpackhi_ps(_tmp1,jy1);                                          \
     jy2            = _mm_unpacklo_ps(_tmp2,jz2);                                          \
     jx3            = _mm_unpackhi_ps(_tmp2,jz2);                                          \
     jy1            = _mm_movehl_ps(jx1,jx1);                                              \
     jx2            = _mm_movehl_ps(jz1,jz1);                                              \
     jz2            = _mm_movehl_ps(jy2,jy2);                                              \
     jy3            = _mm_movehl_ps(jx3,jx3);                                              \
     jz3            = _mm_unpacklo_ps(jz3,_tmp3);                                          \
}


#define GMX_MM_LOAD_4RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1, _tmp2, _tmp3,_tmp4;                                                                 \
	 _tmp1          = _mm_loadu_ps(ptr1);                                                              \
     jy1            = _mm_loadu_ps(ptr2);                                                              \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                            \
     jz2            = _mm_loadu_ps(ptr2+4);                                                            \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                                            \
     _tmp4          = _mm_loadu_ps(ptr2+8);                                                            \
     jx1            = _mm_unpacklo_ps(_tmp1,jy1);                                                      \
     jz1            = _mm_unpackhi_ps(_tmp1,jy1);                                                      \
     jy2            = _mm_unpacklo_ps(_tmp2,jz2);                                                      \
     jx3            = _mm_unpackhi_ps(_tmp2,jz2);                                                      \
     jz3            = _mm_unpacklo_ps(_tmp3,_tmp4);                                                    \
     jy4            = _mm_unpackhi_ps(_tmp3,_tmp4);                                                    \
     jy1            = _mm_movehl_ps(jx1,jx1);                                                          \
     jx2            = _mm_movehl_ps(jz1,jz1);                                                          \
     jz2            = _mm_movehl_ps(jy2,jy2);                                                          \
     jy3            = _mm_movehl_ps(jx3,jx3);                                                          \
     jx4            = _mm_movehl_ps(jz3,jz3);                                                          \
     jz4            = _mm_movehl_ps(jy4,jy4);                                                          \
}


#define GMX_MM_LOAD_1RVEC_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1) { \
     __m128 _tmp1,_tmp3,_tmp4;                                         \
	 jx1            = _mm_load_ss(ptr1);                               \
     jy1            = _mm_load_ss(ptr2);                               \
     jz1            = _mm_load_ss(ptr3);                               \
	 jx1            = _mm_loadh_pi(jx1,(__m64 *)(ptr1+1));             \
     jy1            = _mm_loadh_pi(jy1,(__m64 *)(ptr2+1));             \
     jz1            = _mm_loadh_pi(jz1,(__m64 *)(ptr3+1));             \
     _tmp1          = _mm_unpacklo_ps(jx1,jy1);                        \
     _tmp3          = _mm_unpackhi_ps(jx1,jy1);                        \
     _tmp4          = _mm_unpackhi_ps(jz1,jz1);                        \
     jx1            = _mm_movelh_ps(_tmp1,jz1);                        \
     jy1            = _mm_movelh_ps(_tmp3,_tmp4);                      \
     jz1            = _mm_movehl_ps(_tmp4,_tmp3);                      \
}


#define GMX_MM_LOAD_2RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1, _tmp2;                                                           \
	 jx1            = _mm_loadu_ps(ptr1);                                           \
     jy1            = _mm_loadu_ps(ptr2);                                           \
     jz1            = _mm_loadu_ps(ptr3);                                           \
     jx2            = _mm_setzero_ps();                                             \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                            \
     _tmp1          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));             \
     jz2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2+4));             \
     _tmp2          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3+4));             \
     _tmp1          = _mm_unpacklo_ps(_tmp1,_tmp2);                                 \
     jz2            = _mm_unpacklo_ps(jz2,_mm_setzero_ps());                        \
     jy2            = _mm_unpacklo_ps(_tmp1,jz2);                                   \
     jz2            = _mm_unpackhi_ps(_tmp1,jz2);                                   \
}


#define GMX_MM_LOAD_3RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1, _tmp2;                                                                       \
	 jx1            = _mm_loadu_ps(ptr1);                                                       \
     jy1            = _mm_loadu_ps(ptr2);                                                       \
     jz1            = _mm_loadu_ps(ptr3);                                                       \
     jx2            = _mm_setzero_ps();                                                         \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                        \
     jy2            = _mm_loadu_ps(ptr1+4);                                                     \
     jz2            = _mm_loadu_ps(ptr2+4);                                                     \
     jx3            = _mm_loadu_ps(ptr3+4);                                                     \
     jy3            = _mm_setzero_ps();                                                         \
     _MM_TRANSPOSE4_PS(jy2,jz2,jx3,jy3);                                                        \
     jz3            = _mm_load_ss(ptr1+8);                                                      \
     _tmp1          = _mm_load_ss(ptr2+8);                                                      \
     _tmp2          = _mm_load_ss(ptr3+8);                                                      \
     jz3            = _mm_unpacklo_ps(jz3,_tmp2);                                               \
     _tmp1          = _mm_unpacklo_ps(_tmp1,_mm_setzero_ps());                                  \
     jz3            = _mm_unpacklo_ps(jz3,_tmp1);                                               \
}


#define GMX_MM_LOAD_4RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
	 jx1            = _mm_loadu_ps(ptr1);                                                                   \
     jy1            = _mm_loadu_ps(ptr2);                                                                   \
     jz1            = _mm_loadu_ps(ptr3);                                                                   \
     jx2            = _mm_setzero_ps();                                                                     \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                                    \
     jy2            = _mm_loadu_ps(ptr1+4);                                                                 \
     jz2            = _mm_loadu_ps(ptr2+4);                                                                 \
     jx3            = _mm_loadu_ps(ptr3+4);                                                                 \
     jy3            = _mm_setzero_ps();                                                                     \
     _MM_TRANSPOSE4_PS(jy2,jz2,jx3,jy3);                                                                    \
     jz3            = _mm_loadu_ps(ptr1+8);                                                                 \
     jx4            = _mm_loadu_ps(ptr2+8);                                                                 \
     jy4            = _mm_loadu_ps(ptr3+8);                                                                 \
     jz4            = _mm_setzero_ps();                                                                     \
     _MM_TRANSPOSE4_PS(jz3,jx4,jy4,jz4);                                                                    \
}



#define GMX_MM_LOAD_1RVEC_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1) {  \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5;                                   \
	 jx1            = _mm_load_ss(ptr1);                                     \
     _tmp1          = _mm_load_ss(ptr2);                                     \
     jy1            = _mm_load_ss(ptr3);                                     \
     jz1            = _mm_load_ss(ptr4);                                     \
	 jx1            = _mm_loadh_pi(jx1,(__m64 *)(ptr1+1));                   \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr2+1));                 \
     jy1            = _mm_loadh_pi(jy1,(__m64 *)(ptr3+1));                   \
     jz1            = _mm_loadh_pi(jz1,(__m64 *)(ptr4+1));                   \
     _tmp2          = _mm_unpacklo_ps(jx1,_tmp1);                            \
     _tmp3          = _mm_unpacklo_ps(jy1,jz1);                              \
     _tmp4          = _mm_unpackhi_ps(jx1,_tmp1);                            \
     _tmp5          = _mm_unpackhi_ps(jy1,jz1);                              \
     jx1            = _mm_movelh_ps(_tmp2,_tmp3);                            \
     jy1            = _mm_movelh_ps(_tmp4,_tmp5);                            \
     jz1            = _mm_movehl_ps(_tmp5,_tmp4);                            \
}


#define GMX_MM_LOAD_2RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1, _tmp2;                                                                \
	 jx1            = _mm_loadu_ps(ptr1);                                                \
     jy1            = _mm_loadu_ps(ptr2);                                                \
     jz1            = _mm_loadu_ps(ptr3);                                                \
     jx2            = _mm_loadu_ps(ptr4);                                                \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                 \
     jy2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));                  \
     jz2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2+4));                  \
     _tmp1          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3+4));                  \
     _tmp2          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr4+4));                  \
     _tmp1          = _mm_unpacklo_ps(jy2,_tmp1);                                        \
     _tmp2          = _mm_unpacklo_ps(jz2,_tmp2);                                        \
     jy2            = _mm_unpacklo_ps(_tmp1,_tmp2);                                      \
     jz2            = _mm_unpackhi_ps(_tmp1,_tmp2);                                      \
}


#define GMX_MM_LOAD_3RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1, _tmp2, _tmp3;                                                                     \
	 jx1            = _mm_loadu_ps(ptr1);                                                            \
     jy1            = _mm_loadu_ps(ptr2);                                                            \
     jz1            = _mm_loadu_ps(ptr3);                                                            \
     jx2            = _mm_loadu_ps(ptr4);                                                            \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                             \
     jy2            = _mm_loadu_ps(ptr1+4);                                                          \
     jz2            = _mm_loadu_ps(ptr2+4);                                                          \
     jx3            = _mm_loadu_ps(ptr3+4);                                                          \
     jy3            = _mm_loadu_ps(ptr4+4);                                                          \
     _MM_TRANSPOSE4_PS(jy2,jz2,jx3,jy3);                                                             \
     jz3            = _mm_load_ss(ptr1+8);                                                           \
     _tmp1          = _mm_load_ss(ptr2+8);                                                           \
     _tmp2          = _mm_load_ss(ptr3+8);                                                           \
     _tmp3          = _mm_load_ss(ptr4+8);                                                           \
     jz3            = _mm_unpacklo_ps(jz3,_tmp2);                                                    \
     _tmp1          = _mm_unpacklo_ps(_tmp1,_tmp3);                                                  \
     jz3            = _mm_unpacklo_ps(jz3,_tmp1);                                                    \
}


#define GMX_MM_LOAD_4RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
	 jx1            = _mm_loadu_ps(ptr1);                                                                        \
     jy1            = _mm_loadu_ps(ptr2);                                                                        \
     jz1            = _mm_loadu_ps(ptr3);                                                                        \
     jx2            = _mm_loadu_ps(ptr4);                                                                        \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                                         \
     jy2            = _mm_loadu_ps(ptr1+4);                                                                      \
     jz2            = _mm_loadu_ps(ptr2+4);                                                                      \
     jx3            = _mm_loadu_ps(ptr3+4);                                                                      \
     jy3            = _mm_loadu_ps(ptr4+4);                                                                      \
     _MM_TRANSPOSE4_PS(jy2,jz2,jx3,jy3);                                                                         \
     jz3            = _mm_loadu_ps(ptr1+8);                                                                      \
     jx4            = _mm_loadu_ps(ptr2+8);                                                                      \
     jy4            = _mm_loadu_ps(ptr3+8);                                                                      \
     jz4            = _mm_loadu_ps(ptr4+8);                                                                      \
     _MM_TRANSPOSE4_PS(jz3,jx4,jy4,jz4);                                                                         \
}


/* Routines to increment rvecs in memory, typically use for j particle force updates */
#define GMX_MM_INCREMENT_1RVEC_1POINTER_PS(ptr1,jx1,jy1,jz1) {      \
     __m128 _tmp1;                                                    \
     jy1            = _mm_unpacklo_ps(jy1,jz1);                       \
     jx1            = _mm_movelh_ps(jx1,jy1);                         \
     _tmp1          = _mm_load_ss(ptr1);                              \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));          \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                          \
     _mm_store_ss(ptr1,_tmp1);                                        \
     _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                          \
}


#define GMX_MM_INCREMENT_2RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1, _tmp2;                                                     \
     _tmp1          = _mm_loadu_ps(ptr1);                                     \
     _tmp2          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));       \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                               \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                               \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                               \
     jx1            = _mm_movelh_ps(jx1,jz1);                                 \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                                  \
     _tmp2          = _mm_add_ps(_tmp2,jy2);                                  \
     _mm_storeu_ps(ptr1,_tmp1);                                               \
     _mm_storel_pi((__m64 *)(ptr1+4),_tmp2);                                  \
}


#define GMX_MM_INCREMENT_3RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1, _tmp2, _tmp3;                                                          \
     _tmp1          = _mm_loadu_ps(ptr1);                                                 \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                               \
     _tmp3          = _mm_load_ss(ptr1+8);                                                \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                           \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                           \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                           \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                           \
     jx1            = _mm_movelh_ps(jx1,jz1);                                             \
     jy2            = _mm_movelh_ps(jy2,jx3);                                             \
     _tmp1           = _mm_add_ps(_tmp1,jx1);                                             \
     _tmp2           = _mm_add_ps(_tmp2,jy2);                                             \
     _tmp3           = _mm_add_ss(_tmp3,jz3);                                             \
     _mm_storeu_ps(ptr1,_tmp1);                                                           \
     _mm_storeu_ps(ptr1+4,_tmp2);                                                         \
     _mm_store_ss(ptr1+8,_tmp3);                                                          \
}


#define GMX_MM_INCREMENT_4RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1, _tmp2, _tmp3;                                                                      \
     _tmp1          = _mm_loadu_ps(ptr1);                                                             \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                           \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                                           \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                                       \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                                       \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                                       \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                                       \
     jz3            = _mm_unpacklo_ps(jz3,jx4);                                                       \
     jy4            = _mm_unpacklo_ps(jy4,jz4);                                                       \
     jx1            = _mm_movelh_ps(jx1,jz1);                                                         \
     jy2            = _mm_movelh_ps(jy2,jx3);                                                         \
     jz3            = _mm_movelh_ps(jz3,jy4);                                                         \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                                                          \
     _tmp2          = _mm_add_ps(_tmp2,jy2);                                                          \
     _tmp3          = _mm_add_ps(_tmp3,jz3);                                                          \
     _mm_storeu_ps(ptr1,_tmp1);                                                                       \
     _mm_storeu_ps(ptr1+4,_tmp2);                                                                     \
     _mm_storeu_ps(ptr1+8,_tmp3);                                                                     \
}


#define GMX_MM_INCREMENT_1RVEC_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1) {        \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4;                                          \
     _tmp1          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1));         \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr2));                    \
     _tmp2          = _mm_load_ss(ptr1+2);                                    \
     _tmp3          = _mm_load_ss(ptr2+2);                                    \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                               \
     _tmp4          = _mm_shuffle_ps(jz1,jz1,_MM_SHUFFLE(0,0,0,1));           \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                                  \
     _mm_storel_pi((__m64 *)(ptr1),_tmp1);                                    \
     _mm_storeh_pi((__m64 *)(ptr2),_tmp1);                                    \
     _mm_store_ss(ptr1+2,_mm_add_ss(_tmp2,jz1));                              \
	 _mm_store_ss(ptr2+2,_mm_add_ss(_tmp3,_tmp4));                            \
}


#define GMX_MM_INCREMENT_2RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) {  \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5;                                           \
     _tmp1          = _mm_loadu_ps(ptr1);                                            \
     _tmp2          = _mm_loadu_ps(ptr2);                                            \
     _tmp3          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));              \
     _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr2+4));                         \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                      \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                      \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                      \
     _tmp4          = _mm_movelh_ps(jx1,jz1);                                        \
     _tmp5          = _mm_movehl_ps(jz1,jx1);                                        \
     _tmp1          = _mm_add_ps(_tmp1,_tmp4);                                       \
     _tmp2          = _mm_add_ps(_tmp2,_tmp5);                                       \
     _tmp3          = _mm_add_ps(_tmp3,jy2);                                         \
     _mm_storeu_ps(ptr1,_tmp1);                                                      \
     _mm_storeu_ps(ptr2,_tmp2);                                                      \
     _mm_storel_pi((__m64 *)(ptr1+4),_tmp3);                                         \
	 _mm_storeh_pi((__m64 *)(ptr2+4),_tmp3);                                         \
}


#define GMX_MM_INCREMENT_3RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;                \
     _tmp1          = _mm_loadu_ps(ptr1);                                                       \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                     \
     _tmp3          = _mm_load_ss(ptr1+8);                                                      \
     _tmp4          = _mm_loadu_ps(ptr2);                                                       \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                                     \
     _tmp6          = _mm_load_ss(ptr2+8);                                                      \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                                 \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                                 \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                                 \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                                 \
     _tmp7          = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));                             \
     _tmp8          = _mm_movelh_ps(jx1,jz1);                                                   \
     _tmp9          = _mm_movehl_ps(jz1,jx1);                                                   \
     _tmp10         = _mm_movelh_ps(jy2,jx3);                                                   \
     _tmp11         = _mm_movehl_ps(jx3,jy2);                                                   \
     _tmp1          = _mm_add_ps(_tmp1,_tmp8);                                                  \
     _tmp2          = _mm_add_ps(_tmp2,_tmp10);                                                 \
     _tmp3          = _mm_add_ss(_tmp3,jz3);                                                    \
     _tmp4          = _mm_add_ps(_tmp4,_tmp9);                                                  \
     _tmp5          = _mm_add_ps(_tmp5,_tmp11);                                                 \
     _tmp6          = _mm_add_ss(_tmp6,_tmp7);                                                  \
     _mm_storeu_ps(ptr1,_tmp1);                                                                 \
     _mm_storeu_ps(ptr1+4,_tmp2);                                                               \
     _mm_store_ss(ptr1+8,_tmp3);                                                                \
     _mm_storeu_ps(ptr2,_tmp4);                                                                 \
     _mm_storeu_ps(ptr2+4,_tmp5);                                                               \
     _mm_store_ss(ptr2+8,_tmp6);                                                                \
}


#define GMX_MM_INCREMENT_4RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11,_tmp12,_tmp13;              \
     _tmp1          = _mm_loadu_ps(ptr1);                                                                   \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                                 \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                                                 \
     _tmp4          = _mm_loadu_ps(ptr2);                                                                   \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                                                 \
     _tmp6          = _mm_loadu_ps(ptr2+8);                                                                 \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                                             \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                                             \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                                             \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                                             \
     jz3            = _mm_unpacklo_ps(jz3,jx4);                                                             \
     jy4            = _mm_unpacklo_ps(jy4,jz4);                                                             \
     _tmp8          = _mm_movelh_ps(jx1,jz1);                                                               \
     _tmp9          = _mm_movehl_ps(jz1,jx1);                                                               \
     _tmp10         = _mm_movelh_ps(jy2,jx3);                                                               \
     _tmp11         = _mm_movehl_ps(jx3,jy2);                                                               \
     _tmp12         = _mm_movelh_ps(jz3,jy4);                                                               \
     _tmp13         = _mm_movehl_ps(jy4,jz3);                                                               \
     _tmp1          = _mm_add_ps(_tmp1,_tmp8);                                                              \
     _tmp2          = _mm_add_ps(_tmp2,_tmp10);                                                             \
     _tmp3          = _mm_add_ps(_tmp3,_tmp12);                                                             \
     _tmp4          = _mm_add_ps(_tmp4,_tmp9);                                                              \
     _tmp5          = _mm_add_ps(_tmp5,_tmp11);                                                             \
     _tmp6          = _mm_add_ps(_tmp6,_tmp13);                                                             \
     _mm_storeu_ps(ptr1,_tmp1);                                                                             \
     _mm_storeu_ps(ptr1+4,_tmp2);                                                                           \
     _mm_storeu_ps(ptr1+8,_tmp3);                                                                           \
     _mm_storeu_ps(ptr2,_tmp4);                                                                             \
     _mm_storeu_ps(ptr2+4,_tmp5);                                                                           \
     _mm_storeu_ps(ptr2+8,_tmp6);                                                                           \
}


#define GMX_MM_INCREMENT_1RVEC_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1) {   \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7;                        \
     _tmp1          = _mm_load_ss(ptr1);                                      \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));                  \
     _tmp2          = _mm_load_ss(ptr2);                                      \
     _tmp2          = _mm_loadh_pi(_tmp2,(__m64 *)(ptr2+1));                  \
     _tmp3          = _mm_load_ss(ptr3);                                      \
     _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr3+1));                  \
     _tmp4          = _mm_unpacklo_ps(jy1,jz1);                               \
     _tmp5          = _mm_unpackhi_ps(jy1,jz1);                               \
     _tmp6          = _mm_shuffle_ps(jx1,_tmp4,_MM_SHUFFLE(3,2,0,1));         \
     _tmp7          = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,2));           \
     jx1            = _mm_movelh_ps(jx1,_tmp4);                               \
     _tmp7          = _mm_movelh_ps(_tmp7,_tmp5);                             \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                                  \
     _tmp2          = _mm_add_ps(_tmp2,_tmp6);                                \
     _tmp3          = _mm_add_ps(_tmp3,_tmp7);                                \
     _mm_store_ss(ptr1,_tmp1);                                                \
     _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                                  \
     _mm_store_ss(ptr2,_tmp2);                                                \
     _mm_storeh_pi((__m64 *)(ptr2+1),_tmp2);                                  \
     _mm_store_ss(ptr3,_tmp3);                                                \
     _mm_storeh_pi((__m64 *)(ptr3+1),_tmp3);                                  \
}


#define GMX_MM_INCREMENT_2RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;                \
     _tmp1          = _mm_loadu_ps(ptr1);                                                \
     _tmp2          = _mm_loadu_ps(ptr2);                                                \
     _tmp3          = _mm_loadu_ps(ptr3);                                                \
     _tmp4          = _mm_loadl_pi(_tmp4,(__m64 *)(ptr1+4));                             \
     _tmp4          = _mm_loadh_pi(_tmp4,(__m64 *)(ptr2+4));                             \
     _tmp5          = _mm_loadl_pi(_tmp5,(__m64 *)(ptr3+4));                             \
     _tmp6          = _mm_unpackhi_ps(jx1,jy1);                                          \
	 jx1            = _mm_unpacklo_ps(jx1,jy1);                                          \
     _tmp7          = _mm_unpackhi_ps(jz1,jx2);                                          \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                          \
     _tmp8          = _mm_unpackhi_ps(jy2,jz2);                                          \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                          \
     _tmp9          = _mm_movelh_ps(jx1,jz1);                                            \
     _tmp10         = _mm_movehl_ps(jz1,jx1);                                            \
     _tmp6          = _mm_movelh_ps(_tmp6,_tmp7);                                        \
     _tmp1          = _mm_add_ps(_tmp1,_tmp9);                                           \
     _tmp2          = _mm_add_ps(_tmp2,_tmp10);                                          \
     _tmp3          = _mm_add_ps(_tmp3,_tmp6);                                           \
     _tmp4          = _mm_add_ps(_tmp4,jy2);                                             \
     _tmp5          = _mm_add_ps(_tmp5,_tmp8);                                           \
     _mm_storeu_ps(ptr1,_tmp1);                                                          \
     _mm_storeu_ps(ptr2,_tmp2);                                                          \
     _mm_storeu_ps(ptr3,_tmp3);                                                          \
     _mm_storel_pi((__m64 *)(ptr1+4),_tmp4);                                             \
     _mm_storeh_pi((__m64 *)(ptr2+4),_tmp4);                                             \
	 _mm_storel_pi((__m64 *)(ptr3+4),_tmp5);                                             \
}


#define GMX_MM_INCREMENT_3RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;                            \
     __m128 _tmp11,_tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19;                          \
     _tmp1          = _mm_loadu_ps(ptr1);                                                            \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                          \
     _tmp3          = _mm_load_ss(ptr1+8);                                                           \
     _tmp4          = _mm_loadu_ps(ptr2);                                                            \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                                          \
     _tmp6          = _mm_load_ss(ptr2+8);                                                           \
     _tmp7          = _mm_loadu_ps(ptr3);                                                            \
     _tmp8          = _mm_loadu_ps(ptr3+4);                                                          \
     _tmp9          = _mm_load_ss(ptr3+8);                                                           \
     _tmp10         = _mm_unpackhi_ps(jx1,jy1);                                                      \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                                      \
     _tmp11         = _mm_unpackhi_ps(jz1,jx2);                                                      \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                                      \
     _tmp12         = _mm_unpackhi_ps(jy2,jz2);                                                      \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                                      \
     _tmp13         = _mm_unpackhi_ps(jx3,jy3);                                                      \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                                      \
     _tmp14         = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));                                  \
     _tmp15         = _mm_movehl_ps(jz3,jz3);                                                        \
     _tmp16         = _mm_movelh_ps(jx1,jz1);                                                        \
     _tmp17         = _mm_movehl_ps(jz1,jx1);                                                        \
     _tmp10         = _mm_movelh_ps(_tmp10,_tmp11);                                                  \
     _tmp18         = _mm_movelh_ps(jy2,jx3);                                                        \
     _tmp19         = _mm_movehl_ps(jx3,jy2);                                                        \
     _tmp12         = _mm_movelh_ps(_tmp12,_tmp13);                                                  \
     _tmp1          = _mm_add_ps(_tmp1,_tmp16);                                                      \
     _tmp2          = _mm_add_ps(_tmp2,_tmp18);                                                      \
     _tmp3          = _mm_add_ss(_tmp3,jz3);                                                         \
     _tmp4          = _mm_add_ps(_tmp4,_tmp17);                                                      \
     _tmp5          = _mm_add_ps(_tmp5,_tmp19);                                                      \
     _tmp6          = _mm_add_ss(_tmp6,_tmp14);                                                      \
     _tmp7          = _mm_add_ps(_tmp7,_tmp10);                                                      \
     _tmp8          = _mm_add_ps(_tmp8,_tmp12);                                                      \
     _tmp9          = _mm_add_ss(_tmp9,_tmp15);                                                      \
     _mm_storeu_ps(ptr1,_tmp1);                                                                      \
     _mm_storeu_ps(ptr1+4,_tmp2);                                                                    \
     _mm_store_ss(ptr1+8,_tmp3);                                                                     \
     _mm_storeu_ps(ptr2,_tmp4);                                                                      \
     _mm_storeu_ps(ptr2+4,_tmp5);                                                                    \
     _mm_store_ss(ptr2+8,_tmp6);                                                                     \
     _mm_storeu_ps(ptr3,_tmp7);                                                                      \
     _mm_storeu_ps(ptr3+4,_tmp8);                                                                    \
     _mm_store_ss(ptr3+8,_tmp9);                                                                     \
}


#define GMX_MM_INCREMENT_4RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;                                 \
     __m128 _tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19,_tmp20,_tmp21;                               \
     _tmp1          = _mm_loadu_ps(ptr1);                                                                        \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                                      \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                                                      \
     _tmp4          = _mm_loadu_ps(ptr2);                                                                        \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                                                      \
     _tmp6          = _mm_loadu_ps(ptr2+8);                                                                      \
     _tmp7          = _mm_loadu_ps(ptr3);                                                                        \
     _tmp8          = _mm_loadu_ps(ptr3+4);                                                                      \
     _tmp9          = _mm_loadu_ps(ptr3+8);                                                                      \
     _tmp10         = _mm_unpackhi_ps(jx1,jy1);                                                                  \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                                                  \
     _tmp11         = _mm_unpackhi_ps(jz1,jx2);                                                                  \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                                                  \
     _tmp12         = _mm_unpackhi_ps(jy2,jz2);                                                                  \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                                                  \
     _tmp13         = _mm_unpackhi_ps(jx3,jy3);                                                                  \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                                                  \
     _tmp14         = _mm_unpackhi_ps(jz3,jx4);                                                                  \
     jz3            = _mm_unpacklo_ps(jz3,jx4);                                                                  \
     _tmp15         = _mm_unpackhi_ps(jy4,jz4);                                                                  \
     jy4            = _mm_unpacklo_ps(jy4,jz4);                                                                  \
     _tmp16         = _mm_movelh_ps(jx1,jz1);                                                                    \
     _tmp17         = _mm_movehl_ps(jz1,jx1);                                                                    \
     _tmp10         = _mm_movelh_ps(_tmp10,_tmp11);                                                              \
     _tmp18         = _mm_movelh_ps(jy2,jx3);                                                                    \
     _tmp19         = _mm_movehl_ps(jx3,jy2);                                                                    \
     _tmp12         = _mm_movelh_ps(_tmp12,_tmp13);                                                              \
     _tmp20         = _mm_movelh_ps(jz3,jy4);                                                                    \
     _tmp21         = _mm_movehl_ps(jy4,jz3);                                                                    \
     _tmp14         = _mm_movelh_ps(_tmp14,_tmp15);                                                              \
     _tmp1          = _mm_add_ps(_tmp1,_tmp16);                                                                  \
     _tmp2          = _mm_add_ps(_tmp2,_tmp18);                                                                  \
     _tmp3          = _mm_add_ps(_tmp3,_tmp20);                                                                  \
     _tmp4          = _mm_add_ps(_tmp4,_tmp17);                                                                  \
     _tmp5          = _mm_add_ps(_tmp5,_tmp19);                                                                  \
     _tmp6          = _mm_add_ps(_tmp6,_tmp21);                                                                  \
     _tmp7          = _mm_add_ps(_tmp7,_tmp10);                                                                  \
     _tmp8          = _mm_add_ps(_tmp8,_tmp12);                                                                  \
     _tmp9          = _mm_add_ps(_tmp9,_tmp14);                                                                  \
     _mm_storeu_ps(ptr1,_tmp1);                                                                                  \
     _mm_storeu_ps(ptr1+4,_tmp2);                                                                                \
     _mm_storeu_ps(ptr1+8,_tmp3);                                                                                \
     _mm_storeu_ps(ptr2,_tmp4);                                                                                  \
     _mm_storeu_ps(ptr2+4,_tmp5);                                                                                \
     _mm_storeu_ps(ptr2+8,_tmp6);                                                                                \
     _mm_storeu_ps(ptr3,_tmp7);                                                                                  \
     _mm_storeu_ps(ptr3+4,_tmp8);                                                                                \
     _mm_storeu_ps(ptr3+8,_tmp9);                                                                                \
}



#define GMX_MM_INCREMENT_1RVEC_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;        \
     _tmp1          = _mm_load_ss(ptr1);                                         \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));                     \
     _tmp2          = _mm_load_ss(ptr2);                                         \
     _tmp2          = _mm_loadh_pi(_tmp2,(__m64 *)(ptr2+1));                     \
     _tmp3          = _mm_load_ss(ptr3);                                         \
     _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr3+1));                     \
     _tmp4          = _mm_load_ss(ptr4);                                         \
     _tmp4          = _mm_loadh_pi(_tmp4,(__m64 *)(ptr4+1));                     \
     _tmp5          = _mm_unpacklo_ps(jy1,jz1);                                  \
     _tmp6          = _mm_unpackhi_ps(jy1,jz1);                                  \
     _tmp7          = _mm_shuffle_ps(jx1,_tmp5,_MM_SHUFFLE(1,0,0,0));            \
     _tmp8          = _mm_shuffle_ps(jx1,_tmp5,_MM_SHUFFLE(3,2,0,1));            \
     _tmp9          = _mm_shuffle_ps(jx1,_tmp6,_MM_SHUFFLE(1,0,0,2));            \
     _tmp10         = _mm_shuffle_ps(jx1,_tmp6,_MM_SHUFFLE(3,2,0,3));            \
     _tmp1          = _mm_add_ps(_tmp1,_tmp7);                                   \
     _tmp2          = _mm_add_ps(_tmp2,_tmp8);                                   \
     _tmp3          = _mm_add_ps(_tmp3,_tmp9);                                   \
     _tmp4          = _mm_add_ps(_tmp4,_tmp10);                                  \
     _mm_store_ss(ptr1,_tmp1);                                                   \
     _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                                     \
     _mm_store_ss(ptr2,_tmp2);                                                   \
     _mm_storeh_pi((__m64 *)(ptr2+1),_tmp2);                                     \
     _mm_store_ss(ptr3,_tmp3);                                                   \
     _mm_storeh_pi((__m64 *)(ptr3+1),_tmp3);                                     \
     _mm_store_ss(ptr4,_tmp4);                                                   \
     _mm_storeh_pi((__m64 *)(ptr4+1),_tmp4);                                     \
}


#define GMX_MM_INCREMENT_2RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2) {  \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11,_tmp12,_tmp13; \
     _tmp1          = _mm_loadu_ps(ptr1);                                                      \
     _tmp2          = _mm_loadu_ps(ptr2);                                                      \
     _tmp3          = _mm_loadu_ps(ptr3);                                                      \
     _tmp4          = _mm_loadu_ps(ptr4);                                                      \
     _tmp5          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));                        \
     _tmp5          = _mm_loadh_pi(_tmp5,(__m64 *)(ptr2+4));                                   \
     _tmp6          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3+4));                        \
     _tmp6          = _mm_loadh_pi(_tmp6,(__m64 *)(ptr4+4));                                   \
     _tmp7          = _mm_unpackhi_ps(jx1,jy1);                                                \
	 jx1            = _mm_unpacklo_ps(jx1,jy1);                                                \
     _tmp8          = _mm_unpackhi_ps(jz1,jx2);                                                \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                                \
     _tmp9          = _mm_unpackhi_ps(jy2,jz2);                                                \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                                \
     _tmp10         = _mm_movelh_ps(jx1,jz1);                                                  \
     _tmp11         = _mm_movehl_ps(jz1,jx1);                                                  \
     _tmp12         = _mm_movelh_ps(_tmp7,_tmp8);                                              \
     _tmp13         = _mm_movehl_ps(_tmp8,_tmp7);                                              \
     _tmp1          = _mm_add_ps(_tmp1,_tmp10);                                                \
     _tmp2          = _mm_add_ps(_tmp2,_tmp11);                                                \
     _tmp3          = _mm_add_ps(_tmp3,_tmp12);                                                \
     _tmp4          = _mm_add_ps(_tmp4,_tmp13);                                                \
     _tmp5          = _mm_add_ps(_tmp5,jy2);                                                   \
     _tmp6          = _mm_add_ps(_tmp6,_tmp9);                                                 \
     _mm_storeu_ps(ptr1,_tmp1);                                                                \
     _mm_storeu_ps(ptr2,_tmp2);                                                                \
     _mm_storeu_ps(ptr3,_tmp3);                                                                \
     _mm_storeu_ps(ptr4,_tmp4);                                                                \
     _mm_storel_pi((__m64 *)(ptr1+4),_tmp5);                                                   \
     _mm_storeh_pi((__m64 *)(ptr2+4),_tmp5);                                                   \
	 _mm_storel_pi((__m64 *)(ptr3+4),_tmp6);                                                   \
	 _mm_storeh_pi((__m64 *)(ptr4+4),_tmp6);                                                   \
}


#define GMX_MM_INCREMENT_3RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;                                 \
     __m128 _tmp11,_tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19;                               \
     __m128 _tmp20,_tmp21,_tmp22,_tmp23,_tmp24,_tmp25;                                                    \
     _tmp1          = _mm_loadu_ps(ptr1);                                                                 \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                               \
     _tmp3          = _mm_load_ss(ptr1+8);                                                                \
     _tmp4          = _mm_loadu_ps(ptr2);                                                                 \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                                               \
     _tmp6          = _mm_load_ss(ptr2+8);                                                                \
     _tmp7          = _mm_loadu_ps(ptr3);                                                                 \
     _tmp8          = _mm_loadu_ps(ptr3+4);                                                               \
     _tmp9          = _mm_load_ss(ptr3+8);                                                                \
     _tmp10         = _mm_loadu_ps(ptr4);                                                                 \
     _tmp11         = _mm_loadu_ps(ptr4+4);                                                               \
     _tmp12         = _mm_load_ss(ptr4+8);                                                                \
     _tmp13         = _mm_unpackhi_ps(jx1,jy1);                                                           \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                                           \
     _tmp14         = _mm_unpackhi_ps(jz1,jx2);                                                           \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                                           \
     _tmp15         = _mm_unpackhi_ps(jy2,jz2);                                                           \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                                           \
     _tmp16         = _mm_unpackhi_ps(jx3,jy3);                                                           \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                                           \
     _tmp17         = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));                                       \
     _tmp18         = _mm_movehl_ps(jz3,jz3);                                                             \
     _tmp19         = _mm_shuffle_ps(_tmp18,_tmp18,_MM_SHUFFLE(0,0,0,1));                                 \
     _tmp20         = _mm_movelh_ps(jx1,jz1);                                                             \
     _tmp21         = _mm_movehl_ps(jz1,jx1);                                                             \
     _tmp22         = _mm_movelh_ps(_tmp13,_tmp14);                                                       \
     _tmp14         = _mm_movehl_ps(_tmp14,_tmp13);                                                       \
     _tmp23         = _mm_movelh_ps(jy2,jx3);                                                             \
     _tmp24         = _mm_movehl_ps(jx3,jy2);                                                             \
     _tmp25         = _mm_movelh_ps(_tmp15,_tmp16);                                                       \
     _tmp16         = _mm_movehl_ps(_tmp16,_tmp15);                                                       \
     _tmp1          = _mm_add_ps(_tmp1,_tmp20);                                                           \
     _tmp2          = _mm_add_ps(_tmp2,_tmp23);                                                           \
     _tmp3          = _mm_add_ss(_tmp3,jz3);                                                              \
     _tmp4          = _mm_add_ps(_tmp4,_tmp21);                                                           \
     _tmp5          = _mm_add_ps(_tmp5,_tmp24);                                                           \
     _tmp6          = _mm_add_ss(_tmp6,_tmp17);                                                           \
     _tmp7          = _mm_add_ps(_tmp7,_tmp22);                                                           \
     _tmp8          = _mm_add_ps(_tmp8,_tmp25);                                                           \
     _tmp9          = _mm_add_ss(_tmp9,_tmp18);                                                           \
     _tmp10         = _mm_add_ps(_tmp10,_tmp14);                                                          \
     _tmp11         = _mm_add_ps(_tmp11,_tmp16);                                                          \
     _tmp12         = _mm_add_ss(_tmp12,_tmp19);                                                          \
     _mm_storeu_ps(ptr1,_tmp1);                                                                           \
     _mm_storeu_ps(ptr1+4,_tmp2);                                                                         \
     _mm_store_ss(ptr1+8,_tmp3);                                                                          \
     _mm_storeu_ps(ptr2,_tmp4);                                                                           \
     _mm_storeu_ps(ptr2+4,_tmp5);                                                                         \
     _mm_store_ss(ptr2+8,_tmp6);                                                                          \
     _mm_storeu_ps(ptr3,_tmp7);                                                                           \
     _mm_storeu_ps(ptr3+4,_tmp8);                                                                         \
     _mm_store_ss(ptr3+8,_tmp9);                                                                          \
     _mm_storeu_ps(ptr4,_tmp10);                                                                          \
     _mm_storeu_ps(ptr4+4,_tmp11);                                                                        \
     _mm_store_ss(ptr4+8,_tmp12);                                                                         \
}


#define GMX_MM_INCREMENT_4RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;                                      \
     __m128 _tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19,_tmp20,_tmp21,_tmp22;                             \
     __m128 _tmp23,_tmp24;                                                                                            \
     _tmp1          = _mm_loadu_ps(ptr1);                                                                             \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                                                           \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                                                           \
     _tmp4          = _mm_loadu_ps(ptr2);                                                                             \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                                                           \
     _tmp6          = _mm_loadu_ps(ptr2+8);                                                                           \
     _tmp7          = _mm_loadu_ps(ptr3);                                                                             \
     _tmp8          = _mm_loadu_ps(ptr3+4);                                                                           \
     _tmp9          = _mm_loadu_ps(ptr3+8);                                                                           \
     _tmp10         = _mm_loadu_ps(ptr4);                                                                             \
     _tmp11         = _mm_loadu_ps(ptr4+4);                                                                           \
     _tmp12         = _mm_loadu_ps(ptr4+8);                                                                           \
     _tmp13         = _mm_unpackhi_ps(jx1,jy1);                                                                       \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                                                       \
     _tmp14         = _mm_unpackhi_ps(jz1,jx2);                                                                       \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                                                       \
     _tmp15         = _mm_unpackhi_ps(jy2,jz2);                                                                       \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                                                       \
     _tmp16         = _mm_unpackhi_ps(jx3,jy3);                                                                       \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                                                       \
     _tmp17         = _mm_unpackhi_ps(jz3,jx4);                                                                       \
     jz3            = _mm_unpacklo_ps(jz3,jx4);                                                                       \
     _tmp18         = _mm_unpackhi_ps(jy4,jz4);                                                                       \
     jy4            = _mm_unpacklo_ps(jy4,jz4);                                                                       \
     _tmp19         = _mm_movelh_ps(jx1,jz1);                                                                         \
     jz1            = _mm_movehl_ps(jz1,jx1);                                                                         \
     _tmp20         = _mm_movelh_ps(_tmp13,_tmp14);                                                                   \
     _tmp14         = _mm_movehl_ps(_tmp14,_tmp13);                                                                   \
     _tmp21         = _mm_movelh_ps(jy2,jx3);                                                                         \
     jx3            = _mm_movehl_ps(jx3,jy2);                                                                         \
     _tmp22         = _mm_movelh_ps(_tmp15,_tmp16);                                                                   \
     _tmp16         = _mm_movehl_ps(_tmp16,_tmp15);                                                                   \
     _tmp23         = _mm_movelh_ps(jz3,jy4);                                                                         \
     jy4            = _mm_movehl_ps(jy4,jz3);                                                                         \
     _tmp24         = _mm_movelh_ps(_tmp17,_tmp18);                                                                   \
     _tmp18         = _mm_movehl_ps(_tmp18,_tmp17);                                                                   \
     _tmp1          = _mm_add_ps(_tmp1,_tmp19);                                                                       \
     _tmp2          = _mm_add_ps(_tmp2,_tmp21);                                                                       \
     _tmp3          = _mm_add_ps(_tmp3,_tmp23);                                                                       \
     _tmp4          = _mm_add_ps(_tmp4,jz1);                                                                          \
     _tmp5          = _mm_add_ps(_tmp5,jx3);                                                                          \
     _tmp6          = _mm_add_ps(_tmp6,jy4);                                                                          \
     _tmp7          = _mm_add_ps(_tmp7,_tmp20);                                                                       \
     _tmp8          = _mm_add_ps(_tmp8,_tmp22);                                                                       \
     _tmp9          = _mm_add_ps(_tmp9,_tmp24);                                                                       \
     _tmp10         = _mm_add_ps(_tmp10,_tmp14);                                                                      \
     _tmp11         = _mm_add_ps(_tmp11,_tmp16);                                                                      \
     _tmp12         = _mm_add_ps(_tmp12,_tmp18);                                                                      \
     _mm_storeu_ps(ptr1,_tmp1);                                                                                       \
     _mm_storeu_ps(ptr1+4,_tmp2);                                                                                     \
     _mm_storeu_ps(ptr1+8,_tmp3);                                                                                     \
     _mm_storeu_ps(ptr2,_tmp4);                                                                                       \
     _mm_storeu_ps(ptr2+4,_tmp5);                                                                                     \
     _mm_storeu_ps(ptr2+8,_tmp6);                                                                                     \
     _mm_storeu_ps(ptr3,_tmp7);                                                                                       \
     _mm_storeu_ps(ptr3+4,_tmp8);                                                                                     \
     _mm_storeu_ps(ptr3+8,_tmp9);                                                                                     \
     _mm_storeu_ps(ptr4,_tmp10);                                                                                      \
     _mm_storeu_ps(ptr4+4,_tmp11);                                                                                    \
     _mm_storeu_ps(ptr4+8,_tmp12);                                                                                    \
}



#define GMX_MM_DECREMENT_1RVEC_1POINTER_PS(ptr1,jx1,jy1,jz1) {     \
    __m128 _tmp1;                                                    \
    jy1            = _mm_unpacklo_ps(jy1,jz1);                       \
    jx1            = _mm_movelh_ps(jx1,jy1);                         \
    _tmp1          = _mm_load_ss(ptr1);                              \
    _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));          \
    _tmp1          = _mm_sub_ps(_tmp1,jx1);                          \
    _mm_store_ss(ptr1,_tmp1);                                        \
    _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                          \
}


#define GMX_MM_DECREMENT_2RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) { \
    __m128 _tmp1, _tmp2;                                                      \
    _tmp1          = _mm_loadu_ps(ptr1);                                      \
    _tmp2          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));        \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
    jx1            = _mm_movelh_ps(jx1,jz1);                                  \
    _tmp1          = _mm_sub_ps(_tmp1,jx1);                                   \
    _tmp2          = _mm_sub_ps(_tmp2,jy2);                                   \
    _mm_storeu_ps(ptr1,_tmp1);                                                \
    _mm_storel_pi((__m64 *)(ptr1+4),_tmp2);                                   \
}


#define GMX_MM_DECREMENT_3RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
    __m128 _tmp1, _tmp2, _tmp3;                                                           \
    _tmp1          = _mm_loadu_ps(ptr1);                                                  \
    _tmp2          = _mm_loadu_ps(ptr1+4);                                                \
    _tmp3          = _mm_load_ss(ptr1+8);                                                 \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                            \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                            \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                            \
    jx3            = _mm_unpacklo_ps(jx3,jy3);                                            \
    jx1            = _mm_movelh_ps(jx1,jz1);                                              \
    jy2            = _mm_movelh_ps(jy2,jx3);                                              \
    _tmp1          = _mm_sub_ps(_tmp1,jx1);                                               \
    _tmp2          = _mm_sub_ps(_tmp2,jy2);                                               \
    _tmp3          = _mm_sub_ss(_tmp3,jz3);                                               \
    _mm_storeu_ps(ptr1,_tmp1);                                                            \
    _mm_storeu_ps(ptr1+4,_tmp2);                                                          \
    _mm_store_ss(ptr1+8,_tmp3);                                                           \
}


#define GMX_MM_DECREMENT_4RVECS_1POINTER_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
    __m128 _tmp1, _tmp2, _tmp3;                                                                       \
    _tmp1          = _mm_loadu_ps(ptr1);                                                              \
    _tmp2          = _mm_loadu_ps(ptr1+4);                                                            \
    _tmp3          = _mm_loadu_ps(ptr1+8);                                                            \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                                        \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                                        \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                                        \
    jx3            = _mm_unpacklo_ps(jx3,jy3);                                                        \
    jz3            = _mm_unpacklo_ps(jz3,jx4);                                                        \
    jy4            = _mm_unpacklo_ps(jy4,jz4);                                                        \
    jx1            = _mm_movelh_ps(jx1,jz1);                                                          \
    jy2            = _mm_movelh_ps(jy2,jx3);                                                          \
    jz3            = _mm_movelh_ps(jz3,jy4);                                                          \
    _tmp1          = _mm_sub_ps(_tmp1,jx1);                                                           \
    _tmp2          = _mm_sub_ps(_tmp2,jy2);                                                           \
    _tmp3          = _mm_sub_ps(_tmp3,jz3);                                                           \
    _mm_storeu_ps(ptr1,_tmp1);                                                                        \
    _mm_storeu_ps(ptr1+4,_tmp2);                                                                      \
    _mm_storeu_ps(ptr1+8,_tmp3);                                                                      \
}


#define GMX_MM_DECREMENT_1RVEC_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1) {        \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4;                                           \
    _tmp1          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1));          \
    _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr2));                     \
    _tmp2          = _mm_load_ss(ptr1+2);                                     \
    _tmp3          = _mm_load_ss(ptr2+2);                                     \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
    _tmp4          = _mm_shuffle_ps(jz1,jz1,_MM_SHUFFLE(0,0,0,1));            \
    _tmp1          = _mm_sub_ps(_tmp1,jx1);                                   \
    _mm_storel_pi((__m64 *)(ptr1),_tmp1);                                     \
    _mm_storeh_pi((__m64 *)(ptr2),_tmp1);                                     \
    _mm_store_ss(ptr1+2,_mm_sub_ss(_tmp2,jz1));                               \
    _mm_store_ss(ptr2+2,_mm_sub_ss(_tmp3,_tmp4));                             \
}


#define GMX_MM_DECREMENT_2RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5;                                           \
    _tmp1          = _mm_loadu_ps(ptr1);                                            \
    _tmp2          = _mm_loadu_ps(ptr2);                                            \
    _tmp3          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));              \
    _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr2+4));                         \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                      \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                      \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                      \
    _tmp4          = _mm_movelh_ps(jx1,jz1);                                        \
    _tmp5          = _mm_movehl_ps(jz1,jx1);                                        \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp4);                                       \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp5);                                       \
    _tmp3          = _mm_sub_ps(_tmp3,jy2);                                         \
    _mm_storeu_ps(ptr1,_tmp1);                                                      \
    _mm_storeu_ps(ptr2,_tmp2);                                                      \
    _mm_storel_pi((__m64 *)(ptr1+4),_tmp3);                                         \
    _mm_storeh_pi((__m64 *)(ptr2+4),_tmp3);                                         \
}


#define GMX_MM_DECREMENT_3RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) {\
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;                \
    _tmp1          = _mm_loadu_ps(ptr1);                                                       \
    _tmp2          = _mm_loadu_ps(ptr1+4);                                                     \
    _tmp3          = _mm_load_ss(ptr1+8);                                                      \
    _tmp4          = _mm_loadu_ps(ptr2);                                                       \
    _tmp5          = _mm_loadu_ps(ptr2+4);                                                     \
    _tmp6          = _mm_load_ss(ptr2+8);                                                      \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                                 \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                                 \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                                 \
    jx3            = _mm_unpacklo_ps(jx3,jy3);                                                 \
    _tmp7          = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));                             \
    _tmp8          = _mm_movelh_ps(jx1,jz1);                                                   \
    _tmp9          = _mm_movehl_ps(jz1,jx1);                                                   \
    _tmp10         = _mm_movelh_ps(jy2,jx3);                                                   \
    _tmp11         = _mm_movehl_ps(jx3,jy2);                                                   \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp8);                                                  \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp10);                                                 \
    _tmp3          = _mm_sub_ss(_tmp3,jz3);                                                    \
    _tmp4          = _mm_sub_ps(_tmp4,_tmp9);                                                  \
    _tmp5          = _mm_sub_ps(_tmp5,_tmp11);                                                 \
    _tmp6          = _mm_sub_ss(_tmp6,_tmp7);                                                  \
    _mm_storeu_ps(ptr1,_tmp1);                                                                 \
    _mm_storeu_ps(ptr1+4,_tmp2);                                                               \
    _mm_store_ss(ptr1+8,_tmp3);                                                                \
    _mm_storeu_ps(ptr2,_tmp4);                                                                 \
    _mm_storeu_ps(ptr2+4,_tmp5);                                                               \
    _mm_store_ss(ptr2+8,_tmp6);                                                                \
}


#define GMX_MM_DECREMENT_4RVECS_2POINTERS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) {\
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11,_tmp12,_tmp13;              \
    _tmp1          = _mm_loadu_ps(ptr1);                                                                   \
    _tmp2          = _mm_loadu_ps(ptr1+4);                                                                 \
    _tmp3          = _mm_loadu_ps(ptr1+8);                                                                 \
    _tmp4          = _mm_loadu_ps(ptr2);                                                                   \
    _tmp5          = _mm_loadu_ps(ptr2+4);                                                                 \
    _tmp6          = _mm_loadu_ps(ptr2+8);                                                                 \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                                             \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                                             \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                                             \
    jx3            = _mm_unpacklo_ps(jx3,jy3);                                                             \
    jz3            = _mm_unpacklo_ps(jz3,jx4);                                                             \
    jy4            = _mm_unpacklo_ps(jy4,jz4);                                                             \
    _tmp8          = _mm_movelh_ps(jx1,jz1);                                                               \
    _tmp9          = _mm_movehl_ps(jz1,jx1);                                                               \
    _tmp10         = _mm_movelh_ps(jy2,jx3);                                                               \
    _tmp11         = _mm_movehl_ps(jx3,jy2);                                                               \
    _tmp12         = _mm_movelh_ps(jz3,jy4);                                                               \
    _tmp13         = _mm_movehl_ps(jy4,jz3);                                                               \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp8);                                                              \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp10);                                                             \
    _tmp3          = _mm_sub_ps(_tmp3,_tmp12);                                                             \
    _tmp4          = _mm_sub_ps(_tmp4,_tmp9);                                                              \
    _tmp5          = _mm_sub_ps(_tmp5,_tmp11);                                                             \
    _tmp6          = _mm_sub_ps(_tmp6,_tmp13);                                                             \
    _mm_storeu_ps(ptr1,_tmp1);                                                                             \
    _mm_storeu_ps(ptr1+4,_tmp2);                                                                           \
    _mm_storeu_ps(ptr1+8,_tmp3);                                                                           \
    _mm_storeu_ps(ptr2,_tmp4);                                                                             \
    _mm_storeu_ps(ptr2+4,_tmp5);                                                                           \
    _mm_storeu_ps(ptr2+8,_tmp6);                                                                           \
}


#define GMX_MM_DECREMENT_1RVEC_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7;                       \
    _tmp1          = _mm_load_ss(ptr1);                                     \
    _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));                 \
    _tmp2          = _mm_load_ss(ptr2);                                     \
    _tmp2          = _mm_loadh_pi(_tmp2,(__m64 *)(ptr2+1));                 \
    _tmp3          = _mm_load_ss(ptr3);                                     \
    _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr3+1));                 \
    _tmp4          = _mm_unpacklo_ps(jy1,jz1);                              \
    _tmp5          = _mm_unpackhi_ps(jy1,jz1);                              \
    _tmp6          = _mm_shuffle_ps(jx1,_tmp4,_MM_SHUFFLE(3,2,0,1));        \
    _tmp7          = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,2));          \
    jx1            = _mm_movelh_ps(jx1,_tmp4);                              \
    _tmp7          = _mm_movelh_ps(_tmp7,_tmp5);                            \
    _tmp1          = _mm_sub_ps(_tmp1,jx1);                                 \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp6);                               \
    _tmp3          = _mm_sub_ps(_tmp3,_tmp7);                               \
    _mm_store_ss(ptr1,_tmp1);                                               \
    _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                                 \
    _mm_store_ss(ptr2,_tmp2);                                               \
    _mm_storeh_pi((__m64 *)(ptr2+1),_tmp2);                                 \
    _mm_store_ss(ptr3,_tmp3);                                               \
    _mm_storeh_pi((__m64 *)(ptr3+1),_tmp3);                                 \
}


#define GMX_MM_DECREMENT_2RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;                 \
    _tmp1          = _mm_loadu_ps(ptr1);                                                 \
    _tmp2          = _mm_loadu_ps(ptr2);                                                 \
    _tmp3          = _mm_loadu_ps(ptr3);                                                 \
    _tmp4          = _mm_loadl_pi(_tmp4,(__m64 *)(ptr1+4));                              \
    _tmp4          = _mm_loadh_pi(_tmp4,(__m64 *)(ptr2+4));                              \
    _tmp5          = _mm_loadl_pi(_tmp5,(__m64 *)(ptr3+4));                              \
    _tmp6          = _mm_unpackhi_ps(jx1,jy1);                                           \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                           \
    _tmp7          = _mm_unpackhi_ps(jz1,jx2);                                           \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                           \
    _tmp8          = _mm_unpackhi_ps(jy2,jz2);                                           \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                           \
    _tmp9          = _mm_movelh_ps(jx1,jz1);                                             \
    _tmp10         = _mm_movehl_ps(jz1,jx1);                                             \
    _tmp6          = _mm_movelh_ps(_tmp6,_tmp7);                                         \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp9);                                            \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp10);                                           \
    _tmp3          = _mm_sub_ps(_tmp3,_tmp6);                                            \
    _tmp4          = _mm_sub_ps(_tmp4,jy2);                                              \
    _tmp5          = _mm_sub_ps(_tmp5,_tmp8);                                            \
    _mm_storeu_ps(ptr1,_tmp1);                                                           \
    _mm_storeu_ps(ptr2,_tmp2);                                                           \
    _mm_storeu_ps(ptr3,_tmp3);                                                           \
    _mm_storel_pi((__m64 *)(ptr1+4),_tmp4);                                              \
    _mm_storeh_pi((__m64 *)(ptr2+4),_tmp4);                                              \
    _mm_storel_pi((__m64 *)(ptr3+4),_tmp5);                                              \
}


#define GMX_MM_DECREMENT_3RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;      \
    __m128 _tmp11,_tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19;    \
    _tmp1          = _mm_loadu_ps(ptr1);                                      \
    _tmp2          = _mm_loadu_ps(ptr1+4);                                    \
    _tmp3          = _mm_load_ss(ptr1+8);                                     \
    _tmp4          = _mm_loadu_ps(ptr2);                                      \
    _tmp5          = _mm_loadu_ps(ptr2+4);                                    \
    _tmp6          = _mm_load_ss(ptr2+8);                                     \
    _tmp7          = _mm_loadu_ps(ptr3);                                      \
    _tmp8          = _mm_loadu_ps(ptr3+4);                                    \
    _tmp9          = _mm_load_ss(ptr3+8);                                     \
    _tmp10         = _mm_unpackhi_ps(jx1,jy1);                                \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
    _tmp11         = _mm_unpackhi_ps(jz1,jx2);                                \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
    _tmp12         = _mm_unpackhi_ps(jy2,jz2);                                \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
    _tmp13         = _mm_unpackhi_ps(jx3,jy3);                                \
    jx3            = _mm_unpacklo_ps(jx3,jy3);                                \
    _tmp14         = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));            \
    _tmp15         = _mm_movehl_ps(jz3,jz3);                                  \
    _tmp16         = _mm_movelh_ps(jx1,jz1);                                  \
    _tmp17         = _mm_movehl_ps(jz1,jx1);                                  \
    _tmp10         = _mm_movelh_ps(_tmp10,_tmp11);                            \
    _tmp18         = _mm_movelh_ps(jy2,jx3);                                  \
    _tmp19         = _mm_movehl_ps(jx3,jy2);                                  \
    _tmp12         = _mm_movelh_ps(_tmp12,_tmp13);                            \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp16);                                \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp18);                                \
    _tmp3          = _mm_sub_ss(_tmp3,jz3);                                   \
    _tmp4          = _mm_sub_ps(_tmp4,_tmp17);                                \
    _tmp5          = _mm_sub_ps(_tmp5,_tmp19);                                \
    _tmp6          = _mm_sub_ss(_tmp6,_tmp14);                                \
    _tmp7          = _mm_sub_ps(_tmp7,_tmp10);                                \
    _tmp8          = _mm_sub_ps(_tmp8,_tmp12);                                \
    _tmp9          = _mm_sub_ss(_tmp9,_tmp15);                                \
    _mm_storeu_ps(ptr1,_tmp1);                                                \
    _mm_storeu_ps(ptr1+4,_tmp2);                                              \
    _mm_store_ss(ptr1+8,_tmp3);                                               \
    _mm_storeu_ps(ptr2,_tmp4);                                                \
    _mm_storeu_ps(ptr2+4,_tmp5);                                              \
    _mm_store_ss(ptr2+8,_tmp6);                                               \
    _mm_storeu_ps(ptr3,_tmp7);                                                \
    _mm_storeu_ps(ptr3+4,_tmp8);                                              \
    _mm_store_ss(ptr3+8,_tmp9);                                               \
}


#define GMX_MM_DECREMENT_4RVECS_3POINTERS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;                                  \
    __m128 _tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19,_tmp20,_tmp21;                                \
    _tmp1          = _mm_loadu_ps(ptr1);                                      \
    _tmp2          = _mm_loadu_ps(ptr1+4);                                    \
    _tmp3          = _mm_loadu_ps(ptr1+8);                                    \
    _tmp4          = _mm_loadu_ps(ptr2);                                      \
    _tmp5          = _mm_loadu_ps(ptr2+4);                                    \
    _tmp6          = _mm_loadu_ps(ptr2+8);                                    \
    _tmp7          = _mm_loadu_ps(ptr3);                                      \
    _tmp8          = _mm_loadu_ps(ptr3+4);                                    \
    _tmp9          = _mm_loadu_ps(ptr3+8);                                    \
    _tmp10         = _mm_unpackhi_ps(jx1,jy1);                                \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
    _tmp11         = _mm_unpackhi_ps(jz1,jx2);                                \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
    _tmp12         = _mm_unpackhi_ps(jy2,jz2);                                \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
    _tmp13         = _mm_unpackhi_ps(jx3,jy3);                                \
    jx3            = _mm_unpacklo_ps(jx3,jy3);                                \
    _tmp14         = _mm_unpackhi_ps(jz3,jx4);                                \
    jz3            = _mm_unpacklo_ps(jz3,jx4);                                \
    _tmp15         = _mm_unpackhi_ps(jy4,jz4);                                \
    jy4            = _mm_unpacklo_ps(jy4,jz4);                                \
    _tmp16         = _mm_movelh_ps(jx1,jz1);                                  \
    _tmp17         = _mm_movehl_ps(jz1,jx1);                                  \
    _tmp10         = _mm_movelh_ps(_tmp10,_tmp11);                            \
    _tmp18         = _mm_movelh_ps(jy2,jx3);                                  \
    _tmp19         = _mm_movehl_ps(jx3,jy2);                                  \
    _tmp12         = _mm_movelh_ps(_tmp12,_tmp13);                            \
    _tmp20         = _mm_movelh_ps(jz3,jy4);                                  \
    _tmp21         = _mm_movehl_ps(jy4,jz3);                                  \
    _tmp14         = _mm_movelh_ps(_tmp14,_tmp15);                            \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp16);                                \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp18);                                \
    _tmp3          = _mm_sub_ps(_tmp3,_tmp20);                                \
    _tmp4          = _mm_sub_ps(_tmp4,_tmp17);                                \
    _tmp5          = _mm_sub_ps(_tmp5,_tmp19);                                \
    _tmp6          = _mm_sub_ps(_tmp6,_tmp21);                                \
    _tmp7          = _mm_sub_ps(_tmp7,_tmp10);                                \
    _tmp8          = _mm_sub_ps(_tmp8,_tmp12);                                \
    _tmp9          = _mm_sub_ps(_tmp9,_tmp14);                                \
    _mm_storeu_ps(ptr1,_tmp1);                                                \
    _mm_storeu_ps(ptr1+4,_tmp2);                                              \
    _mm_storeu_ps(ptr1+8,_tmp3);                                              \
    _mm_storeu_ps(ptr2,_tmp4);                                                \
    _mm_storeu_ps(ptr2+4,_tmp5);                                              \
    _mm_storeu_ps(ptr2+8,_tmp6);                                              \
    _mm_storeu_ps(ptr3,_tmp7);                                                \
    _mm_storeu_ps(ptr3+4,_tmp8);                                              \
    _mm_storeu_ps(ptr3+8,_tmp9);                                              \
}




#define GMX_MM_DECREMENT_1RVEC_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;         \
    _tmp1          = _mm_load_ss(ptr1);                              \
    _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));          \
    _tmp2          = _mm_load_ss(ptr2);                              \
    _tmp2          = _mm_loadh_pi(_tmp2,(__m64 *)(ptr2+1));          \
    _tmp3          = _mm_load_ss(ptr3);                              \
    _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr3+1));          \
    _tmp4          = _mm_load_ss(ptr4);                              \
    _tmp4          = _mm_loadh_pi(_tmp4,(__m64 *)(ptr4+1));          \
    _tmp5          = _mm_unpacklo_ps(jy1,jz1);                       \
    _tmp6          = _mm_unpackhi_ps(jy1,jz1);                       \
    _tmp7          = _mm_shuffle_ps(jx1,_tmp5,_MM_SHUFFLE(1,0,0,0)); \
    _tmp8          = _mm_shuffle_ps(jx1,_tmp5,_MM_SHUFFLE(3,2,0,1)); \
    _tmp9          = _mm_shuffle_ps(jx1,_tmp6,_MM_SHUFFLE(1,0,0,2)); \
    _tmp10         = _mm_shuffle_ps(jx1,_tmp6,_MM_SHUFFLE(3,2,0,3)); \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp7);                        \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp8);                        \
    _tmp3          = _mm_sub_ps(_tmp3,_tmp9);                        \
    _tmp4          = _mm_sub_ps(_tmp4,_tmp10);                       \
    _mm_store_ss(ptr1,_tmp1);                                        \
    _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                          \
    _mm_store_ss(ptr2,_tmp2);                                        \
    _mm_storeh_pi((__m64 *)(ptr2+1),_tmp2);                          \
    _mm_store_ss(ptr3,_tmp3);                                        \
    _mm_storeh_pi((__m64 *)(ptr3+1),_tmp3);                          \
    _mm_store_ss(ptr4,_tmp4);                                        \
    _mm_storeh_pi((__m64 *)(ptr4+1),_tmp4);                          \
}



#define GMX_MM_DECREMENT_2RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11,_tmp12,_tmp13; \
    _tmp1          = _mm_loadu_ps(ptr1);                                       \
    _tmp2          = _mm_loadu_ps(ptr2);                                       \
    _tmp3          = _mm_loadu_ps(ptr3);                                       \
    _tmp4          = _mm_loadu_ps(ptr4);                                       \
    _tmp5          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));         \
    _tmp5          = _mm_loadh_pi(_tmp5,(__m64 *)(ptr2+4));                    \
    _tmp6          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3+4));         \
    _tmp6          = _mm_loadh_pi(_tmp6,(__m64 *)(ptr4+4));                    \
    _tmp7          = _mm_unpackhi_ps(jx1,jy1);                                 \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                 \
    _tmp8          = _mm_unpackhi_ps(jz1,jx2);                                 \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                 \
    _tmp9          = _mm_unpackhi_ps(jy2,jz2);                                 \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                 \
    _tmp10         = _mm_movelh_ps(jx1,jz1);                                   \
    _tmp11         = _mm_movehl_ps(jz1,jx1);                                   \
    _tmp12         = _mm_movelh_ps(_tmp7,_tmp8);                               \
    _tmp13         = _mm_movehl_ps(_tmp8,_tmp7);                               \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp10);                                 \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp11);                                 \
    _tmp3          = _mm_sub_ps(_tmp3,_tmp12);                                 \
    _tmp4          = _mm_sub_ps(_tmp4,_tmp13);                                 \
    _tmp5          = _mm_sub_ps(_tmp5,jy2);                                    \
    _tmp6          = _mm_sub_ps(_tmp6,_tmp9);                                  \
    _mm_storeu_ps(ptr1,_tmp1);                                                 \
    _mm_storeu_ps(ptr2,_tmp2);                                                 \
    _mm_storeu_ps(ptr3,_tmp3);                                                 \
    _mm_storeu_ps(ptr4,_tmp4);                                                 \
    _mm_storel_pi((__m64 *)(ptr1+4),_tmp5);                                    \
    _mm_storeh_pi((__m64 *)(ptr2+4),_tmp5);                                    \
    _mm_storel_pi((__m64 *)(ptr3+4),_tmp6);                                    \
    _mm_storeh_pi((__m64 *)(ptr4+4),_tmp6);                                    \
}


#define GMX_MM_DECREMENT_3RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;       \
    __m128 _tmp11,_tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19;     \
    __m128 _tmp20,_tmp21,_tmp22,_tmp23,_tmp24,_tmp25;                          \
    _tmp1          = _mm_loadu_ps(ptr1);                                       \
    _tmp2          = _mm_loadu_ps(ptr1+4);                                     \
    _tmp3          = _mm_load_ss(ptr1+8);                                      \
    _tmp4          = _mm_loadu_ps(ptr2);                                       \
    _tmp5          = _mm_loadu_ps(ptr2+4);                                     \
    _tmp6          = _mm_load_ss(ptr2+8);                                      \
    _tmp7          = _mm_loadu_ps(ptr3);                                       \
    _tmp8          = _mm_loadu_ps(ptr3+4);                                     \
    _tmp9          = _mm_load_ss(ptr3+8);                                      \
    _tmp10         = _mm_loadu_ps(ptr4);                                       \
    _tmp11         = _mm_loadu_ps(ptr4+4);                                     \
    _tmp12         = _mm_load_ss(ptr4+8);                                      \
    _tmp13         = _mm_unpackhi_ps(jx1,jy1);                                 \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                 \
    _tmp14         = _mm_unpackhi_ps(jz1,jx2);                                 \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                 \
    _tmp15         = _mm_unpackhi_ps(jy2,jz2);                                 \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                 \
    _tmp16         = _mm_unpackhi_ps(jx3,jy3);                                 \
    jx3            = _mm_unpacklo_ps(jx3,jy3);                                 \
    _tmp17         = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));             \
    _tmp18         = _mm_movehl_ps(jz3,jz3);                                   \
    _tmp19         = _mm_shuffle_ps(_tmp18,_tmp18,_MM_SHUFFLE(0,0,0,1));       \
    _tmp20         = _mm_movelh_ps(jx1,jz1);                                   \
    _tmp21         = _mm_movehl_ps(jz1,jx1);                                   \
    _tmp22         = _mm_movelh_ps(_tmp13,_tmp14);                             \
    _tmp14         = _mm_movehl_ps(_tmp14,_tmp13);                             \
    _tmp23         = _mm_movelh_ps(jy2,jx3);                                   \
    _tmp24         = _mm_movehl_ps(jx3,jy2);                                   \
    _tmp25         = _mm_movelh_ps(_tmp15,_tmp16);                             \
    _tmp16         = _mm_movehl_ps(_tmp16,_tmp15);                             \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp20);                                 \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp23);                                 \
    _tmp3          = _mm_sub_ss(_tmp3,jz3);                                    \
    _tmp4          = _mm_sub_ps(_tmp4,_tmp21);                                 \
    _tmp5          = _mm_sub_ps(_tmp5,_tmp24);                                 \
    _tmp6          = _mm_sub_ss(_tmp6,_tmp17);                                 \
    _tmp7          = _mm_sub_ps(_tmp7,_tmp22);                                 \
    _tmp8          = _mm_sub_ps(_tmp8,_tmp25);                                 \
    _tmp9          = _mm_sub_ss(_tmp9,_tmp18);                                 \
    _tmp10         = _mm_sub_ps(_tmp10,_tmp14);                                \
    _tmp11         = _mm_sub_ps(_tmp11,_tmp16);                                \
    _tmp12         = _mm_sub_ss(_tmp12,_tmp19);                                \
    _mm_storeu_ps(ptr1,_tmp1);                                                 \
    _mm_storeu_ps(ptr1+4,_tmp2);                                               \
    _mm_store_ss(ptr1+8,_tmp3);                                                \
    _mm_storeu_ps(ptr2,_tmp4);                                                 \
    _mm_storeu_ps(ptr2+4,_tmp5);                                               \
    _mm_store_ss(ptr2+8,_tmp6);                                                \
    _mm_storeu_ps(ptr3,_tmp7);                                                 \
    _mm_storeu_ps(ptr3+4,_tmp8);                                               \
    _mm_store_ss(ptr3+8,_tmp9);                                                \
    _mm_storeu_ps(ptr4,_tmp10);                                                \
    _mm_storeu_ps(ptr4+4,_tmp11);                                              \
    _mm_store_ss(ptr4+8,_tmp12);                                               \
}


#define GMX_MM_DECREMENT_4RVECS_4POINTERS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
    __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;         \
    __m128 _tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19,_tmp20,_tmp21,_tmp22;\
    __m128 _tmp23,_tmp24;                                                     \
    _tmp1          = _mm_loadu_ps(ptr1);                                      \
    _tmp2          = _mm_loadu_ps(ptr1+4);                                    \
    _tmp3          = _mm_loadu_ps(ptr1+8);                                    \
    _tmp4          = _mm_loadu_ps(ptr2);                                      \
    _tmp5          = _mm_loadu_ps(ptr2+4);                                    \
    _tmp6          = _mm_loadu_ps(ptr2+8);                                    \
    _tmp7          = _mm_loadu_ps(ptr3);                                      \
    _tmp8          = _mm_loadu_ps(ptr3+4);                                    \
    _tmp9          = _mm_loadu_ps(ptr3+8);                                    \
    _tmp10         = _mm_loadu_ps(ptr4);                                      \
    _tmp11         = _mm_loadu_ps(ptr4+4);                                    \
    _tmp12         = _mm_loadu_ps(ptr4+8);                                    \
    _tmp13         = _mm_unpackhi_ps(jx1,jy1);                                \
    jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
    _tmp14         = _mm_unpackhi_ps(jz1,jx2);                                \
    jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
    _tmp15         = _mm_unpackhi_ps(jy2,jz2);                                \
    jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
    _tmp16         = _mm_unpackhi_ps(jx3,jy3);                                \
    jx3            = _mm_unpacklo_ps(jx3,jy3);                                \
    _tmp17         = _mm_unpackhi_ps(jz3,jx4);                                \
    jz3            = _mm_unpacklo_ps(jz3,jx4);                                \
    _tmp18         = _mm_unpackhi_ps(jy4,jz4);                                \
    jy4            = _mm_unpacklo_ps(jy4,jz4);                                \
    _tmp19         = _mm_movelh_ps(jx1,jz1);                                  \
    jz1            = _mm_movehl_ps(jz1,jx1);                                  \
    _tmp20         = _mm_movelh_ps(_tmp13,_tmp14);                            \
    _tmp14         = _mm_movehl_ps(_tmp14,_tmp13);                            \
    _tmp21         = _mm_movelh_ps(jy2,jx3);                                  \
    jx3            = _mm_movehl_ps(jx3,jy2);                                  \
    _tmp22         = _mm_movelh_ps(_tmp15,_tmp16);                            \
    _tmp16         = _mm_movehl_ps(_tmp16,_tmp15);                            \
    _tmp23         = _mm_movelh_ps(jz3,jy4);                                  \
    jy4            = _mm_movehl_ps(jy4,jz3);                                  \
    _tmp24         = _mm_movelh_ps(_tmp17,_tmp18);                            \
    _tmp18         = _mm_movehl_ps(_tmp18,_tmp17);                            \
    _tmp1          = _mm_sub_ps(_tmp1,_tmp19);                                \
    _tmp2          = _mm_sub_ps(_tmp2,_tmp21);                                \
    _tmp3          = _mm_sub_ps(_tmp3,_tmp23);                                \
    _tmp4          = _mm_sub_ps(_tmp4,jz1);                                   \
    _tmp5          = _mm_sub_ps(_tmp5,jx3);                                   \
    _tmp6          = _mm_sub_ps(_tmp6,jy4);                                   \
    _tmp7          = _mm_sub_ps(_tmp7,_tmp20);                                \
    _tmp8          = _mm_sub_ps(_tmp8,_tmp22);                                \
    _tmp9          = _mm_sub_ps(_tmp9,_tmp24);                                \
    _tmp10         = _mm_sub_ps(_tmp10,_tmp14);                               \
    _tmp11         = _mm_sub_ps(_tmp11,_tmp16);                               \
    _tmp12         = _mm_sub_ps(_tmp12,_tmp18);                               \
    _mm_storeu_ps(ptr1,_tmp1);                                                \
    _mm_storeu_ps(ptr1+4,_tmp2);                                              \
    _mm_storeu_ps(ptr1+8,_tmp3);                                              \
    _mm_storeu_ps(ptr2,_tmp4);                                                \
    _mm_storeu_ps(ptr2+4,_tmp5);                                              \
    _mm_storeu_ps(ptr2+8,_tmp6);                                              \
    _mm_storeu_ps(ptr3,_tmp7);                                                \
    _mm_storeu_ps(ptr3+4,_tmp8);                                              \
    _mm_storeu_ps(ptr3+8,_tmp9);                                              \
    _mm_storeu_ps(ptr4,_tmp10);                                               \
    _mm_storeu_ps(ptr4+4,_tmp11);                                             \
    _mm_storeu_ps(ptr4+8,_tmp12);                                             \
}






/* Routine to be called with rswitch/rcut at the beginning of a kernel
 * to set up the 7 constants used for analytic 5th order switch calculations.
 */
#define GMX_MM_SETUP_SWITCH5_PS(rswitch,rcut,switch_C3,switch_C4,switch_C5,switch_D2,switch_D3,switch_D4) {  \
	const __m128  _swsetup_cm6  = _mm_set_ps( -6.0, -6.0, -6.0, -6.0);                                                 \
	const __m128 _swsetup_cm10  = _mm_set_ps(-10.0,-10.0,-10.0,-10.0);                                                 \
	const __m128  _swsetup_c15  = _mm_set_ps( 15.0, 15.0, 15.0, 15.0);                                                 \
	const __m128 _swsetup_cm30  = _mm_set_ps(-30.0,-30.0,-30.0,-30.0);                                                 \
	const __m128  _swsetup_c60  = _mm_set_ps( 60.0, 60.0, 60.0, 60.0);                                                 \
                                                                                                             \
	__m128 d,dinv,dinv2,dinv3,dinv4,dinv5;                                                                   \
	                                                                                                         \
	d       = _mm_sub_ps(rcut,rswitch);                                                                      \
	dinv    = gmx_mm_inv_ps(d);                                                                              \
	dinv2   = _mm_mul_ps(dinv,dinv);                                                                         \
	dinv3   = _mm_mul_ps(dinv2,dinv);                                                                        \
	dinv4   = _mm_mul_ps(dinv2,dinv2);                                                                       \
	dinv5   = _mm_mul_ps(dinv3,dinv2);                                                                       \
	                                                                                                         \
	switch_C3 = _mm_mul_ps(_swsetup_cm10,dinv3);                                                             \
	switch_C4 = _mm_mul_ps(_swsetup_c15,dinv4);                                                              \
	switch_C5 = _mm_mul_ps(_swsetup_cm6,dinv5);                                                              \
	switch_D2 = _mm_mul_ps(_swsetup_cm30,dinv3);                                                             \
	switch_D3 = _mm_mul_ps(_swsetup_c60,dinv4);                                                              \
	switch_D4 = _mm_mul_ps(_swsetup_cm30,dinv5);                                                             \
}


#define GMX_MM_EVALUATE_SWITCH5_PS(r,rswitch,rcut,sw,dsw,sw_C3,sw_C4,sw_C5,sw_D2,sw_D3,sw_D4) { \
    const __m128  _sw_one  = _mm_set_ps(  1.0,  1.0,  1.0,  1.0);                                         \
    __m128 d,d2;                                                                                \
    d     = _mm_max_ps(r,rswitch);                                                              \
    d     = _mm_min_ps(d,rcut);                                                                 \
    d     = _mm_sub_ps(d,rswitch);                                                              \
    d2    = _mm_mul_ps(d,d);                                                                    \
    sw    = _mm_mul_ps(d,sw_C5);                                                                \
    dsw   = _mm_mul_ps(d,sw_D4);                                                                \
    sw    = _mm_add_ps(sw,sw_C4);                                                               \
    dsw   = _mm_add_ps(dsw,sw_D3);                                                              \
    sw    = _mm_mul_ps(sw,d);                                                                   \
    dsw   = _mm_mul_ps(dsw,d);                                                                  \
    sw    = _mm_add_ps(sw,sw_C3);                                                               \
    dsw   = _mm_add_ps(dsw,sw_D2);                                                              \
    sw    = _mm_mul_ps(sw,_mm_mul_ps(d,d2));                                                    \
    dsw   = _mm_mul_ps(dsw,d2);                                                                 \
    sw    = _mm_add_ps(sw,_sw_one);                                                             \
}



static inline void
gmx_mm_update_iforce_1atom_ps(__m128 *fix1, __m128 *fiy1, __m128 *fiz1,
                              float *fptr,
                              float *fshiftptr)
{
	__m128 t1,t2,t3;
	
#ifdef GMX_SSE3
	*fix1 = _mm_hadd_ps(*fix1,*fix1);   
	*fiy1 = _mm_hadd_ps(*fiy1,*fiz1);
	
	*fix1 = _mm_hadd_ps(*fix1,*fiy1); /* fiz1 fiy1 fix1 fix1 */ 
#else
	/* SSE2 */
	/* transpose data */
	t1 = *fix1;
	_MM_TRANSPOSE4_PS(*fix1,t1,*fiy1,*fiz1);  
	*fix1 = _mm_add_ps(_mm_add_ps(*fix1,t1), _mm_add_ps(*fiy1,*fiz1));
#endif
	t2 = _mm_load_ss(fptr);
	t2 = _mm_loadh_pi(t2,(__m64 *)(fptr+1));
	t3 = _mm_load_ss(fshiftptr);
	t3 = _mm_loadh_pi(t3,(__m64 *)(fshiftptr+1));
	
	t2 = _mm_add_ps(t2,*fix1);
	t3 = _mm_add_ps(t3,*fix1);
	
	_mm_store_ss(fptr,t2);
	_mm_storeh_pi((__m64 *)(fptr+1),t2);
	_mm_store_ss(fshiftptr,t3);
	_mm_storeh_pi((__m64 *)(fshiftptr+1),t3);
}

static inline void
gmx_mm_update_iforce_2atoms_ps(__m128 *fix1, __m128 *fiy1, __m128 *fiz1,
                               __m128 *fix2, __m128 *fiy2, __m128 *fiz2,
                               float *fptr,
                               float *fshiftptr)
{
	__m128 t1,t2,t4;
	
#ifdef GMX_SSE3
	*fix1 = _mm_hadd_ps(*fix1,*fiy1);   
	*fiz1 = _mm_hadd_ps(*fiz1,*fix2);
	*fiy2 = _mm_hadd_ps(*fiy2,*fiz2);
	
	*fix1 = _mm_hadd_ps(*fix1,*fiz1); /* fix2 fiz1 fiy1 fix1 */ 
	*fiy2 = _mm_hadd_ps(*fiy2,*fiy2); /*  -    -   fiz2 fiy2 */
#else
	/* SSE2 */
	/* transpose data */
	_MM_TRANSPOSE4_PS(*fix1,*fiy1,*fiz1,*fix2);  
	t1 = _mm_unpacklo_ps(*fiy2,*fiz2);
	t2 = _mm_unpackhi_ps(*fiy2,*fiz2);
		
	*fix1 = _mm_add_ps(_mm_add_ps(*fix1,*fiy1), _mm_add_ps(*fiz1,*fix2));
	t1   = _mm_add_ps(t1,t2);
	t2   = _mm_movehl_ps(t2,t1);
	*fiy2 = _mm_add_ps(t1,t2);
#endif
	_mm_storeu_ps(fptr,   _mm_add_ps(*fix1,_mm_loadu_ps(fptr)  ));
	t1 = _mm_loadl_pi(t1,(__m64 *)(fptr+4));
	_mm_storel_pi((__m64 *)(fptr+4), _mm_add_ps(*fiy2,t1));
	
	t4 = _mm_load_ss(fshiftptr+2);
	t4 = _mm_loadh_pi(t4,(__m64 *)(fshiftptr));
	
	t1 = _mm_shuffle_ps(*fix1,*fiy2,_MM_SHUFFLE(0,0,3,2));   /* fiy2  -   fix2 fiz1 */
	t1 = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(3,1,0,0));       /* fiy2 fix2  -   fiz1 */
	t2 = _mm_shuffle_ps(*fiy2,*fix1,_MM_SHUFFLE(1,0,0,1));   /* fiy1 fix1  -   fiz2 */

	t1 = _mm_add_ps(t1,t2);
	t1 = _mm_add_ps(t1,t4); /* y x - z */
	
	_mm_store_ss(fshiftptr+2,t1);
	_mm_storeh_pi((__m64 *)(fshiftptr),t1);
}



static inline void
gmx_mm_update_iforce_3atoms_ps(__m128 *fix1, __m128 *fiy1, __m128 *fiz1,
                               __m128 *fix2, __m128 *fiy2, __m128 *fiz2,
                               __m128 *fix3, __m128 *fiy3, __m128 *fiz3,
                               float *fptr,
                               float *fshiftptr)
{
	__m128 t1,t2,t3,t4;
	
#ifdef GMX_SSE3
	*fix1 = _mm_hadd_ps(*fix1,*fiy1);   
	*fiz1 = _mm_hadd_ps(*fiz1,*fix2);
	*fiy2 = _mm_hadd_ps(*fiy2,*fiz2);
	*fix3 = _mm_hadd_ps(*fix3,*fiy3);
	*fiz3 = _mm_hadd_ps(*fiz3,*fiz3);
	
	*fix1 = _mm_hadd_ps(*fix1,*fiz1); /* fix2 fiz1 fiy1 fix1 */ 
	*fiy2 = _mm_hadd_ps(*fiy2,*fix3); /* fiy3 fix3 fiz2 fiy2 */
	*fiz3 = _mm_hadd_ps(*fiz3,*fiz3); /*  -    -    -   fiz3 */
#else
	/* SSE2 */
	/* transpose data */
	_MM_TRANSPOSE4_PS(*fix1,*fiy1,*fiz1,*fix2);  
	_MM_TRANSPOSE4_PS(*fiy2,*fiz2,*fix3,*fiy3);
	t2   = _mm_movehl_ps(_mm_setzero_ps(),*fiz3);
	t1   = _mm_shuffle_ps(*fiz3,*fiz3,_MM_SHUFFLE(0,0,0,1));
	t3   = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(0,0,0,1));
	
	*fix1 = _mm_add_ps(_mm_add_ps(*fix1,*fiy1), _mm_add_ps(*fiz1,*fix2));
	*fiy2 = _mm_add_ps(_mm_add_ps(*fiy2,*fiz2), _mm_add_ps(*fix3,*fiy3));
	*fiz3 = _mm_add_ss(_mm_add_ps(*fiz3,t1)  , _mm_add_ps(t2,t3));
#endif
	_mm_storeu_ps(fptr,  _mm_add_ps(*fix1,_mm_loadu_ps(fptr)  ));
	_mm_storeu_ps(fptr+4,_mm_add_ps(*fiy2,_mm_loadu_ps(fptr+4)));
	_mm_store_ss (fptr+8,_mm_add_ss(*fiz3,_mm_load_ss(fptr+8) ));
	
	t4 = _mm_load_ss(fshiftptr+2);
	t4 = _mm_loadh_pi(t4,(__m64 *)(fshiftptr));
	
	t1 = _mm_shuffle_ps(*fiz3,*fix1,_MM_SHUFFLE(1,0,0,0));   /* fiy1 fix1  -   fiz3 */
	t2 = _mm_shuffle_ps(*fix1,*fiy2,_MM_SHUFFLE(3,2,2,2));   /* fiy3 fix3  -   fiz1 */
	t3 = _mm_shuffle_ps(*fiy2,*fix1,_MM_SHUFFLE(3,3,0,1));   /* fix2 fix2 fiy2 fiz2 */
	t3 = _mm_shuffle_ps(t3  ,t3  ,_MM_SHUFFLE(1,2,0,0));   /* fiy2 fix2  -   fiz2 */

	t1 = _mm_add_ps(t1,t2);
	t3 = _mm_add_ps(t3,t4);
	t1 = _mm_add_ps(t1,t3); /* y x - z */
	
	_mm_store_ss(fshiftptr+2,t1);
	_mm_storeh_pi((__m64 *)(fshiftptr),t1);
}


static inline void
gmx_mm_update_iforce_4atoms_ps(__m128 *fix1, __m128 *fiy1, __m128 *fiz1,
                               __m128 *fix2, __m128 *fiy2, __m128 *fiz2,
                               __m128 *fix3, __m128 *fiy3, __m128 *fiz3,
                               __m128 *fix4, __m128 *fiy4, __m128 *fiz4,
                               float *fptr,
                               float *fshiftptr)
{
	__m128 t1,t2,t3,t4,t5;
	
#ifdef GMX_SSE3
	*fix1 = _mm_hadd_ps(*fix1,*fiy1);   
	*fiz1 = _mm_hadd_ps(*fiz1,*fix2);
	*fiy2 = _mm_hadd_ps(*fiy2,*fiz2);
	*fix3 = _mm_hadd_ps(*fix3,*fiy3);
	*fiz3 = _mm_hadd_ps(*fiz3,*fix4);
	*fiy4 = _mm_hadd_ps(*fiy4,*fiz4);
	
	*fix1 = _mm_hadd_ps(*fix1,*fiz1); /* fix2 fiz1 fiy1 fix1 */ 
	*fiy2 = _mm_hadd_ps(*fiy2,*fix3); /* fiy3 fix3 fiz2 fiy2 */
	*fiz3 = _mm_hadd_ps(*fiz3,*fiy4); /* fiz4 fiy4 fix4 fiz3 */
#else
	/* SSE2 */
	/* transpose data */
	_MM_TRANSPOSE4_PS(*fix1,*fiy1,*fiz1,*fix2);  
	_MM_TRANSPOSE4_PS(*fiy2,*fiz2,*fix3,*fiy3);
	_MM_TRANSPOSE4_PS(*fiz3,*fix4,*fiy4,*fiz4);
	
	*fix1 = _mm_add_ps(_mm_add_ps(*fix1,*fiy1), _mm_add_ps(*fiz1,*fix2));
	*fiy2 = _mm_add_ps(_mm_add_ps(*fiy2,*fiz2), _mm_add_ps(*fix3,*fiy3));
	*fiz3 = _mm_add_ps(_mm_add_ps(*fiz3,*fix4), _mm_add_ps(*fiy4,*fiz4));
#endif
	_mm_storeu_ps(fptr,  _mm_add_ps(*fix1,_mm_loadu_ps(fptr)  ));
	_mm_storeu_ps(fptr+4,_mm_add_ps(*fiy2,_mm_loadu_ps(fptr+4)));
	_mm_storeu_ps(fptr+8,_mm_add_ps(*fiz3,_mm_loadu_ps(fptr+8)));
	
	t5 = _mm_load_ss(fshiftptr+2);
	t5 = _mm_loadh_pi(t5,(__m64 *)(fshiftptr));
	
	t1 = _mm_shuffle_ps(*fix1,*fix1,_MM_SHUFFLE(1,0,2,2));   /* fiy1 fix1  -   fiz1 */
	t2 = _mm_shuffle_ps(*fiy2,*fiy2,_MM_SHUFFLE(3,2,1,1));   /* fiy3 fix3  -   fiz2 */
	t3 = _mm_shuffle_ps(*fiz3,*fiz3,_MM_SHUFFLE(2,1,0,0));   /* fiy4 fix4  -   fiz3 */
	t4 = _mm_shuffle_ps(*fix1,*fiy2,_MM_SHUFFLE(0,0,3,3));   /* fiy2 fiy2 fix2 fix2 */
	t4 = _mm_shuffle_ps(*fiz3,t4  ,_MM_SHUFFLE(2,0,3,3));   /* fiy2 fix2  -   fiz4 */
	
	t1 = _mm_add_ps(t1,t2);
	t3 = _mm_add_ps(t3,t4);
	t1 = _mm_add_ps(t1,t3); /* y x - z */
	t5 = _mm_add_ps(t5,t1);
	
	_mm_store_ss(fshiftptr+2,t5);
	_mm_storeh_pi((__m64 *)(fshiftptr),t5);
}


#ifdef GMX_SSE3

#define GMX_MM_UPDATE_1POT_PS(pot1,ptr1)                                    \
{                                                                           \
    pot1 = _mm_hadd_ps(pot1,pot1);                                          \
    pot1 = _mm_hadd_ps(pot1,pot1);                                          \
    _mm_store_ss(ptr1,_mm_add_ss(pot1,_mm_load_ss(ptr1)));                  \
}

#define GMX_MM_UPDATE_2POT_PS(pot1,ptr1,pot2,ptr2)                          \
{                                                                           \
    pot1 = _mm_hadd_ps(pot1,pot2);                                          \
    pot1 = _mm_hadd_ps(pot1,pot1);                                          \
    pot2 = _mm_shuffle_ps(pot1,pot1,_MM_SHUFFLE(0,0,0,1));                  \
	_mm_store_ss(ptr1,_mm_add_ss(pot1,_mm_load_ss(ptr1)));                  \
	_mm_store_ss(ptr2,_mm_add_ss(pot2,_mm_load_ss(ptr2)));                  \
}

#else

#define GMX_MM_UPDATE_1POT_PS(pot1,ptr1)                                    \
{                                                                           \
    pot1 = _mm_add_ps(pot1,_mm_movehl_ps(pot1,pot1));                       \
    pot1 = _mm_add_ps(pot1,_mm_shuffle_ps(pot1,pot1,_MM_SHUFFLE(0,0,0,1))); \
    _mm_store_ss(ptr1,_mm_add_ss(pot1,_mm_load_ss(ptr1)));                  \
}

#define GMX_MM_UPDATE_2POT_PS(pot1,ptr1,pot2,ptr2)                          \
{                                                                           \
	__m128 _updt1_,_updt2;                                                  \
	_updt1 = _mm_movehl_ps(pot2,pot1); /* 2d 2c 1d 1c */                    \
	_updt2 = _mm_movelh_ps(pot1,pot2); /* 2b 2a 1b 1a */                    \
	_updt1 = _mm_add_ps(_updt1,_updt2);       /* 2  2  1  1  */             \
	_updt2 = _mm_shuffle_ps(_updt1,_updt1,_MM_SHUFFLE(3,3,1,1));            \
	pot1   = _mm_add_ps(_updt1,_updt2);       /* -  2  -  1  */             \
	pot2   = _mm_movehl_ps(_updt2,pot1);    /* -  -  -  2  */               \
	_mm_store_ss(ptr1,_mm_add_ss(pot1,_mm_load_ss(ptr1)));                  \
	_mm_store_ss(ptr2,_mm_add_ss(pot2,_mm_load_ss(ptr2)));                  \
}

#endif


#define GMX_MM_UPDATE_4POT_PS(pot1,ptr1,pot2,ptr2,pot3,ptr3,pot4,ptr4)      \
{                                                                           \
    _MM_TRANSPOSE4_PS(pot1,pot2,pot3,pot4);                                 \
    pot1 = _mm_add_ps(_mm_add_ps(pot1,pot2),_mm_add_ps(pot3,pot4));         \
    pot2 = _mm_shuffle_ps(pot1,pot1,_MM_SHUFFLE(1,1,1,1));                  \
    pot3 = _mm_shuffle_ps(pot1,pot1,_MM_SHUFFLE(2,2,2,2));                  \
    pot4 = _mm_shuffle_ps(pot1,pot1,_MM_SHUFFLE(3,3,3,3));                  \
	_mm_store_ss(ptr1,_mm_add_ss(pot1,_mm_load_ss(ptr1)));                  \
	_mm_store_ss(ptr2,_mm_add_ss(pot2,_mm_load_ss(ptr2)));                  \
	_mm_store_ss(ptr3,_mm_add_ss(pot3,_mm_load_ss(ptr3)));                  \
	_mm_store_ss(ptr4,_mm_add_ss(pot4,_mm_load_ss(ptr4)));                  \
}


#endif /* _gmx_sse2_single_h_ */
