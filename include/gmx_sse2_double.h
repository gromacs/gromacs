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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#  define gmx_mm_castsi128_pd(a) _mm_castsi128_pd(a)
#  define gmx_mm_castpd_si128(a) _mm_castpd_si128(a)
#  define gmx_mm_castpd_pd128(a) (a)
#elif defined(__GNUC__)
#  define gmx_mm_castsi128_pd(a) ((__m128d)(a))
#  define gmx_mm_castpd_si128(a) ((__m128i)(a))
#  define gmx_mm_castpd_pd128(a) ((__m128d)(a))
#else
static __m128d gmx_mm_castsi128_pd(__m128i a) { return *(__m128d *) &a; } 
static __m128i gmx_mm_castpd_si128(__m128d a) { return *(__m128i *) &a; } 
static __m128d gmx_mm_castpd_pd128(__m128d a) { return *(__m128d *) &a; } 
#endif




static inline void
gmx_printxmm_pd(const char *s,__m128d xmm)
{
	double f[2];
	
	_mm_storeu_pd(f,xmm);
	printf("%s: %15.10g %15.10g\n",s,f[0],f[1]);	
}


static inline __m128d
gmx_mm_inv_pd(__m128d x)
{
	const __m128d two  = (const __m128d) {2.0,2.0};
	
	/* Lookup instruction only exists in single precision, convert back and forth... */
	__m128d lu = _mm_cvtps_pd(_mm_rcp_ps( _mm_cvtpd_ps(x)));
	
	/* Perform two N-R steps for double precision */
	lu         = _mm_mul_pd(lu,_mm_sub_pd(two,_mm_mul_pd(x,lu)));
	return (__m128d) _mm_mul_pd(lu,_mm_sub_pd(two,_mm_mul_pd(x,lu)));
}



static inline __m128d
gmx_mm_invsqrt_pd(__m128d x)
{
	const __m128d half  = (const __m128d) {0.5,0.5};
	const __m128d three = (const __m128d) {3.0,3.0};
    
	/* Lookup instruction only exists in single precision, convert back and forth... */
	__m128d lu = _mm_cvtps_pd(_mm_rsqrt_ps( _mm_cvtpd_ps(x)));
	
	lu = _mm_mul_pd(half,_mm_mul_pd(_mm_sub_pd(three,_mm_mul_pd(_mm_mul_pd(lu,lu),x)),lu));
	return (__m128d) _mm_mul_pd(half,_mm_mul_pd(_mm_sub_pd(three,_mm_mul_pd(_mm_mul_pd(lu,lu),x)),lu));
}

static inline __m128d
gmx_mm_sqrt_pd(__m128d x)
{
    __m128d mask;
    __m128d res;
    
    mask = _mm_cmpeq_pd(x,_mm_setzero_pd());
    res  = _mm_andnot_pd(mask,gmx_mm_invsqrt_pd(x));
    
    res  = _mm_mul_pd(x,res);
    
    return res;
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
 * (1023 in double), and shifting the integer to the exponent field of the FP number (52 bits up).
 *
 * The 2^w term is calculated from a (10,0)-th order (no denominator) Minimax polynomia on the interval
 * [-0.5,0.5]. The coefficiencts of this was derived in Mathematica using the command:
 *
 * MiniMaxApproximation[(2^x), {x, {-0.5, 0.5}, 10, 0}, WorkingPrecision -> 20]
 *
 * The lowest exponent we can represent in IEEE double-precision binary format is 2^-1022; below that 
 * it will wrap around and lead to very large positive numbers. This corresponds to a lower bound
 * on the argument for exp(x) of roughly -708.39. For smaller arguments the return value will be 0.0.
 *
 * There appears to be a slight loss of precision for large arguments (~250), where the largest relative
 * error reaches ~2e-14. However, since the actual value for that argument is around 1E100, it might
 * not matter for typical workloads. This is likely caused by the polynomial evaluation,
 * and the only way around would then be a table-based version, which I haven't managed to get the
 * same performance from.
 * 
 * The _average_ accuracy is about 51 bits in the range [-20,20], and the worst roughly 1 bit worse.
 */
static __m128d 
gmx_mm_exp_pd(__m128d x)
{
    const __m128d argscale = _mm_set1_pd(1.442695040888963387);
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m128d arglimit = _mm_set1_pd(-1022.0/1.442695040888963387);
    const __m128i expbase  = _mm_set1_epi32(1023);
    
    const __m128d CA0       = _mm_set1_pd(7.0372789822689374920e-9);
    const __m128d CB0       = _mm_set1_pd(89.491964762085371);
    const __m128d CB1       = _mm_set1_pd(-9.7373870675164587);
    const __m128d CC0       = _mm_set1_pd(51.247261867992408);
    const __m128d CC1       = _mm_set1_pd(-0.184020268133945);
    const __m128d CD0       = _mm_set1_pd(36.82070153762337);
    const __m128d CD1       = _mm_set1_pd(5.416849282638991);
    const __m128d CE0       = _mm_set1_pd(30.34003452248759);
    const __m128d CE1       = _mm_set1_pd(8.726173289493301);
    const __m128d CF0       = _mm_set1_pd(27.73526969472330);
    const __m128d CF1       = _mm_set1_pd(10.284755658866532);    
    
    __m128d valuemask;
    __m128i iexppart;
    __m128d fexppart;
    __m128d intpart;
    __m128d z,z2;
    __m128d factB,factC,factD,factE,factF;
    
	z         = _mm_mul_pd(x,argscale);  
    iexppart  = _mm_cvtpd_epi32(z);    
    
#if GMX_SSE4
    /* This reduces latency and speeds up the code by roughly 5% when supported */
    intpart   = _mm_round_pd(z,0);
#else
    intpart   = _mm_cvtepi32_pd(iexppart);
#endif
    /* The two lowest elements of iexppart now contains 32-bit numbers with a correctly biased exponent.
     * To be able to shift it into the exponent for a double precision number we first need to
     * shuffle so that the lower half contains the first element, and the upper half the second.
     * This should really be done as a zero-extension, but since the next instructions will shift
     * the registers left by 52 bits it doesn't matter what we put there - it will be shifted out.
     * (thus we just use element 2 from iexppart).
     */
    iexppart  = _mm_shuffle_epi32(iexppart,_MM_SHUFFLE(2,1,2,0));  
    
    /* Do the shift operation on the 64-bit registers */
    iexppart  = _mm_add_epi32(iexppart,expbase);
    iexppart  = _mm_slli_epi64(iexppart,52);
    valuemask = _mm_cmpgt_pd(x,arglimit);
    
	z         = _mm_sub_pd(z,intpart);
    z2        = _mm_mul_pd(z,z);  
    
    fexppart  = _mm_and_pd(valuemask,gmx_mm_castsi128_pd(iexppart));
    
    /* Since SSE doubleing-point has relatively high latency it is faster to do 
     * factorized polynomial summation with independent terms than using alternating add/multiply, i.e.
     * p(z) = A0 * (B0 + z) * (C0 + C1*z + z^2) * (D0 + D1*z + z^2) * (E0 + E1*z + z^2) * (F0 + F1*z + z^2) 
     */
    
    factB     = _mm_add_pd(CB0,_mm_mul_pd(CB1,z) );
    factB     = _mm_add_pd(factB,z2);
    factC     = _mm_add_pd(CC0,_mm_mul_pd(CC1,z) );
    factC     = _mm_add_pd(factC,z2);
    factD     = _mm_add_pd(CD0,_mm_mul_pd(CD1,z) );
    factD     = _mm_add_pd(factD,z2);
    factE     = _mm_add_pd(CE0,_mm_mul_pd(CE1,z) );
    factE     = _mm_add_pd(factE,z2);
    factF     = _mm_add_pd(CF0,_mm_mul_pd(CF1,z) );
    factF     = _mm_add_pd(factF,z2);
    
    z         = _mm_mul_pd(CA0,fexppart);
    factB     = _mm_mul_pd(factB,factC);
    factD     = _mm_mul_pd(factD,factE);
    z         = _mm_mul_pd(z,factF);
    factB     = _mm_mul_pd(factB,factD);
    z         = _mm_mul_pd(z,factB);
    
	/* Currently uses 32 actual (real, not including casts) SSE instructions */
	return  z;
}

static __m128d
gmx_mm_log_pd(__m128d x)
{
    /* Same algorithm as cephes library */
	const __m128d expmask    = gmx_mm_castsi128_pd( _mm_set_epi32(0x7FF00000, 0x00000000, 0x7FF00000, 0x00000000) );
    
    const __m128i expbase_m1 = _mm_set1_epi32(1023-1); /* We want non-IEEE format */
    const __m128d half       = _mm_set1_pd(0.5);
    const __m128d one        = _mm_set1_pd(1.0);
    const __m128i itwo       = _mm_set1_epi32(2);
    const __m128d invsq2     = _mm_set1_pd(1.0/sqrt(2.0));
    
    const __m128d corr1      = _mm_set1_pd(-2.121944400546905827679e-4);
    const __m128d corr2      = _mm_set1_pd(0.693359375);
    
    const __m128d P5         = _mm_set1_pd(1.01875663804580931796E-4);
    const __m128d P4         = _mm_set1_pd(4.97494994976747001425E-1);
    const __m128d P3         = _mm_set1_pd(4.70579119878881725854E0);
    const __m128d P2         = _mm_set1_pd(1.44989225341610930846E1);
    const __m128d P1         = _mm_set1_pd(1.79368678507819816313E1);
    const __m128d P0         = _mm_set1_pd(7.70838733755885391666E0);
    
    const __m128d Q4         = _mm_set1_pd(1.12873587189167450590E1);
    const __m128d Q3         = _mm_set1_pd(4.52279145837532221105E1);
    const __m128d Q2         = _mm_set1_pd(8.29875266912776603211E1);
    const __m128d Q1         = _mm_set1_pd(7.11544750618563894466E1);
    const __m128d Q0         = _mm_set1_pd(2.31251620126765340583E1);
    
    const __m128d R2         = _mm_set1_pd(-7.89580278884799154124E-1);
    const __m128d R1         = _mm_set1_pd(1.63866645699558079767E1);
    const __m128d R0         = _mm_set1_pd(-6.41409952958715622951E1);
    
    const __m128d S2         = _mm_set1_pd(-3.56722798256324312549E1);
    const __m128d S1         = _mm_set1_pd(3.12093766372244180303E2);
    const __m128d S0         = _mm_set1_pd(-7.69691943550460008604E2);
    
    __m128d fexp;
    __m128i iexp,iexp1,signbit,iexpabs;
    __m128i imask1;
    __m128d mask1,mask2;
    __m128d corr,t1,t2,q;
    __m128d zA,yA,xA,zB,yB,xB,z;
    __m128d polyR,polyS;
    __m128d polyP1,polyP2,polyQ1,polyQ2;
    
    /* Separate x into exponent and mantissa, with a mantissa in the range [0.5..1[ (not IEEE754 standard!) */
    fexp   = _mm_and_pd(x,expmask);
    iexp   = gmx_mm_castpd_si128(fexp);
    iexp   = _mm_srli_epi64(iexp,52);
    iexp   = _mm_sub_epi32(iexp,expbase_m1);
    iexp   = _mm_shuffle_epi32(iexp, _MM_SHUFFLE(1,1,2,0) );
    
    x      = _mm_andnot_pd(expmask,x);
    x      = _mm_or_pd(x,one);
    x      = _mm_mul_pd(x,half);
    
    signbit = _mm_srai_epi32(iexp,31);
    iexpabs = _mm_sub_epi32(_mm_xor_si128(iexp,signbit),signbit);
    
    imask1 = _mm_cmpgt_epi32( iexpabs, itwo );
    mask2  = _mm_cmplt_pd(x,invsq2);
    
    fexp   = _mm_cvtepi32_pd(iexp);    
    corr   = _mm_and_pd(mask2,one);
    fexp   = _mm_sub_pd(fexp,corr);
    
    /* If mask1 is set ('A') */
    zA     = _mm_sub_pd(x,half);
    t1     = _mm_or_pd( _mm_andnot_pd(mask2,zA), _mm_and_pd(mask2,x) );
    zA     = _mm_sub_pd(t1,half);
    t2     = _mm_or_pd( _mm_andnot_pd(mask2,x), _mm_and_pd(mask2,zA) );
    yA     = _mm_mul_pd(half,_mm_add_pd(t2,one));
    
    xA     = _mm_mul_pd(zA,gmx_mm_inv_pd(yA));
    zA     = _mm_mul_pd(xA,xA);
    
    /* EVALUATE POLY */
    polyR  = _mm_mul_pd(R2,zA);
    polyR  = _mm_add_pd(polyR,R1);
    polyR  = _mm_mul_pd(polyR,zA);
    polyR  = _mm_add_pd(polyR,R0);
    
    polyS  = _mm_add_pd(zA,S2);
    polyS  = _mm_mul_pd(polyS,zA);
    polyS  = _mm_add_pd(polyS,S1);
    polyS  = _mm_mul_pd(polyS,zA);
    polyS  = _mm_add_pd(polyS,S0);
    
    q      = _mm_mul_pd(polyR,gmx_mm_inv_pd(polyS));
    zA     = _mm_mul_pd(_mm_mul_pd(xA,zA),q);
    
    zA     = _mm_add_pd(zA,_mm_mul_pd(corr1,fexp));
    zA     = _mm_add_pd(zA,xA);
    zA     = _mm_add_pd(zA,_mm_mul_pd(corr2,fexp));
    
    
    /* If mask1 is not set ('B') */
    corr   = _mm_and_pd(mask2,x);
    xB     = _mm_add_pd(x,corr);
    xB     = _mm_sub_pd(xB,one);
    zB     = _mm_mul_pd(xB,xB);
    
    polyP1 = _mm_mul_pd(P5,zB);
    polyP2 = _mm_mul_pd(P4,zB);
    polyP1 = _mm_add_pd(polyP1,P3);
    polyP2 = _mm_add_pd(polyP2,P2);
    polyP1 = _mm_mul_pd(polyP1,zB);
    polyP2 = _mm_mul_pd(polyP2,zB);
    polyP1 = _mm_add_pd(polyP1,P1);
    polyP2 = _mm_add_pd(polyP2,P0);
    polyP1 = _mm_mul_pd(polyP1,xB);
    polyP1 = _mm_add_pd(polyP1,polyP2);
    
    polyQ2 = _mm_mul_pd(Q4,zB);
    polyQ1 = _mm_add_pd(zB,Q3);
    polyQ2 = _mm_add_pd(polyQ2,Q2);
    polyQ1 = _mm_mul_pd(polyQ1,zB);
    polyQ2 = _mm_mul_pd(polyQ2,zB);
    polyQ1 = _mm_add_pd(polyQ1,Q1);
    polyQ2 = _mm_add_pd(polyQ2,Q0);
    polyQ1 = _mm_mul_pd(polyQ1,xB);
    polyQ1 = _mm_add_pd(polyQ1,polyQ2);
    
    q      = _mm_mul_pd(polyP1,gmx_mm_inv_pd(polyQ1));
    yB     = _mm_mul_pd(_mm_mul_pd(xB,zB),q);
    
    yB     = _mm_add_pd(yB,_mm_mul_pd(corr1,fexp));
    yB     = _mm_sub_pd(yB,_mm_mul_pd(half,zB));
    zB     = _mm_add_pd(xB,yB);
    zB     = _mm_add_pd(zB,_mm_mul_pd(corr2,fexp));
    
    
    mask1  = gmx_mm_castsi128_pd( _mm_shuffle_epi32(imask1, _MM_SHUFFLE(1,1,0,0)) );
    z      = _mm_or_pd( _mm_and_pd(mask1,zA), _mm_andnot_pd(mask1,zB) );
    
    return z;
}


static int
gmx_mm_sincos_pd(__m128d x,
                 __m128d *sinval,
                 __m128d *cosval)
{
    const __m128d sintable[17] =
    {
        _mm_set_pd( sin(  0.0 * (M_PI/2.0) / 16.0) , cos(  0.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  1.0 * (M_PI/2.0) / 16.0) , cos(  1.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  2.0 * (M_PI/2.0) / 16.0) , cos(  2.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  3.0 * (M_PI/2.0) / 16.0) , cos(  3.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  4.0 * (M_PI/2.0) / 16.0) , cos(  4.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  5.0 * (M_PI/2.0) / 16.0) , cos(  5.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  6.0 * (M_PI/2.0) / 16.0) , cos(  6.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  7.0 * (M_PI/2.0) / 16.0) , cos(  7.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  8.0 * (M_PI/2.0) / 16.0) , cos(  8.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin(  9.0 * (M_PI/2.0) / 16.0) , cos(  9.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin( 10.0 * (M_PI/2.0) / 16.0) , cos( 10.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin( 11.0 * (M_PI/2.0) / 16.0) , cos( 11.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin( 12.0 * (M_PI/2.0) / 16.0) , cos( 12.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin( 13.0 * (M_PI/2.0) / 16.0) , cos( 13.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin( 14.0 * (M_PI/2.0) / 16.0) , cos( 14.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin( 15.0 * (M_PI/2.0) / 16.0) , cos( 15.0 * (M_PI/2.0) / 16.0) ),
        _mm_set_pd( sin( 16.0 * (M_PI/2.0) / 16.0) , cos( 16.0 * (M_PI/2.0) / 16.0) )
    };
    
    const __m128d signmask    = _mm_castsi128_pd( _mm_set_epi32(0x7FFFFFFF,0xFFFFFFFF,0x7FFFFFFF,0xFFFFFFFF) );
    const __m128d tabscale    = _mm_set1_pd(32.0/M_PI);
    const __m128d invtabscale = _mm_set1_pd(M_PI/32.0);
    const __m128d one         = _mm_set1_pd(1.0);
    const __m128d half        = _mm_set1_pd(0.5);
    const __m128i i32         = _mm_set1_epi32(32);
    const __m128i i16         = _mm_set1_epi32(16);
    const __m128i tabmask     = _mm_set1_epi32(0x3F);
    const __m128d sinP7       = _mm_set1_pd(-1.9835958374985742404167983310657359E-4);
    const __m128d sinP5       = _mm_set1_pd(8.3333330133863188710910904896670286347068944E-3);
    const __m128d sinP3       = _mm_set1_pd(-1.66666666666049927649121240304808461825E-1);
    const __m128d sinP1       = _mm_set1_pd(9.99999999999999814240922423058227089729828E-1);
    
    const __m128d cosP6       = _mm_set1_pd(-1.3884108697213500852322439746988228E-3);
    const __m128d cosP4       = _mm_set1_pd(4.16666637872444585215198198619627956E-2);
    const __m128d cosP2       = _mm_set1_pd(-4.999999999944495616220671189143471E-1);
    const __m128d cosP0       = _mm_set1_pd(9.999999999999983282334075742852867E-1);
    
    __m128d scalex;
    __m128i tabidx,corridx;
    __m128d xabs,z,z2,polySin,polyCos;
    __m128d xpoint;
    __m128d ypoint0,ypoint1;
    __m128d sinpoint,cospoint;
    __m128d xsign,ssign,csign;
    __m128i imask,sswapsign,cswapsign;
    __m128d minusone;
    
    xsign    = _mm_andnot_pd(signmask,x);
    xabs     = _mm_and_pd(x,signmask);
    
    scalex   = _mm_mul_pd(tabscale,xabs);
    tabidx   = _mm_cvtpd_epi32(scalex);
    xpoint   = _mm_cvtepi32_pd(tabidx);
    xpoint   = _mm_mul_pd(xpoint,invtabscale);
    z        = _mm_sub_pd(xabs, xpoint);
    
    /* Range reduction to 0..2*Pi */
    tabidx   = _mm_and_si128(tabidx,tabmask);
    
    /* tabidx is now in range [0,..,64] */
    imask     = _mm_cmpgt_epi32(tabidx,i32);
    sswapsign = imask;
    cswapsign = imask;
    corridx   = _mm_and_si128(imask,i32);
    tabidx    = _mm_sub_epi32(tabidx,corridx);
    
    /* tabidx is now in range [0..32] */
    imask     = _mm_cmpgt_epi32(tabidx,i16);
    cswapsign = _mm_xor_si128(cswapsign,imask);
    corridx   = _mm_sub_epi32(i32,tabidx);
    tabidx    = _mm_or_si128( _mm_and_si128(imask,corridx), _mm_andnot_si128(imask,tabidx) );
    
    /* tabidx is now in range [0..16] */
    sswapsign = _mm_shuffle_epi32(sswapsign,_MM_SHUFFLE(2,2,0,0));
    cswapsign = _mm_shuffle_epi32(cswapsign,_MM_SHUFFLE(2,2,0,0));
    minusone  = _mm_sub_pd(_mm_setzero_pd(),one);
    
    ssign     = _mm_or_pd(_mm_and_pd( _mm_castsi128_pd(sswapsign),minusone ),
                          _mm_andnot_pd( _mm_castsi128_pd(sswapsign),one ));
    csign     = _mm_or_pd(_mm_and_pd( _mm_castsi128_pd(cswapsign),minusone ),
                          _mm_andnot_pd( _mm_castsi128_pd(cswapsign),one ));
    
    /* First lookup into table */
    ypoint0  = sintable[gmx_mm_extract_epi32(tabidx,0)];
    ypoint1  = sintable[gmx_mm_extract_epi32(tabidx,1)];
    sinpoint = _mm_unpackhi_pd(ypoint0,ypoint1);
    cospoint = _mm_unpacklo_pd(ypoint0,ypoint1);
    
    sinpoint = _mm_mul_pd(sinpoint,ssign);
    cospoint = _mm_mul_pd(cospoint,csign);
    
    z2       = _mm_mul_pd(z,z);
    
    polySin  = _mm_mul_pd(sinP7,z2);
    polySin  = _mm_add_pd(polySin,sinP5);
    polySin  = _mm_mul_pd(polySin,z2);
    polySin  = _mm_add_pd(polySin,sinP3);
    polySin  = _mm_mul_pd(polySin,z2);
    polySin  = _mm_add_pd(polySin,sinP1);
    polySin  = _mm_mul_pd(polySin,z);
    
    polyCos  = _mm_mul_pd(cosP6,z2);
    polyCos  = _mm_add_pd(polyCos,cosP4);
    polyCos  = _mm_mul_pd(polyCos,z2);
    polyCos  = _mm_add_pd(polyCos,cosP2);
    polyCos  = _mm_mul_pd(polyCos,z2);
    polyCos  = _mm_add_pd(polyCos,cosP0);
    
    *sinval  = _mm_xor_pd(_mm_add_pd( _mm_mul_pd(sinpoint,polyCos) , _mm_mul_pd(cospoint,polySin) ),xsign);
    *cosval  = _mm_sub_pd( _mm_mul_pd(cospoint,polyCos) , _mm_mul_pd(sinpoint,polySin) );
    
    return 0;
}



static __m128d
gmx_mm_tan_pd(__m128d x)
{
    __m128d sinval,cosval;
    __m128d tanval;
    
    gmx_mm_sincos_pd(x,&sinval,&cosval);
    
    tanval = _mm_mul_pd(sinval,gmx_mm_inv_pd(cosval));
    
    return tanval;
}



static __m128d
gmx_mm_asin_pd(__m128d x)
{
    /* Same algorithm as cephes library */
    const __m128d signmask  = gmx_mm_castsi128_pd( _mm_set_epi32(0x7FFFFFFF,0xFFFFFFFF,0x7FFFFFFF,0xFFFFFFFF) );
    const __m128d limit1    = _mm_set1_pd(0.625);
    const __m128d limit2    = _mm_set1_pd(1e-8);
    const __m128d one       = _mm_set1_pd(1.0);
    const __m128d halfpi    = _mm_set1_pd(M_PI/2.0);
    const __m128d quarterpi = _mm_set1_pd(M_PI/4.0);
    const __m128d morebits  = _mm_set1_pd(6.123233995736765886130E-17);
    
    const __m128d P5        = _mm_set1_pd(4.253011369004428248960E-3);
    const __m128d P4        = _mm_set1_pd(-6.019598008014123785661E-1);
    const __m128d P3        = _mm_set1_pd(5.444622390564711410273E0);
    const __m128d P2        = _mm_set1_pd(-1.626247967210700244449E1);
    const __m128d P1        = _mm_set1_pd(1.956261983317594739197E1);
    const __m128d P0        = _mm_set1_pd(-8.198089802484824371615E0);
    
    const __m128d Q4        = _mm_set1_pd(-1.474091372988853791896E1);
    const __m128d Q3        = _mm_set1_pd(7.049610280856842141659E1);
    const __m128d Q2        = _mm_set1_pd(-1.471791292232726029859E2);
    const __m128d Q1        = _mm_set1_pd(1.395105614657485689735E2);
    const __m128d Q0        = _mm_set1_pd(-4.918853881490881290097E1);
    
    const __m128d R4        = _mm_set1_pd(2.967721961301243206100E-3);
    const __m128d R3        = _mm_set1_pd(-5.634242780008963776856E-1);
    const __m128d R2        = _mm_set1_pd(6.968710824104713396794E0);
    const __m128d R1        = _mm_set1_pd(-2.556901049652824852289E1);
    const __m128d R0        = _mm_set1_pd(2.853665548261061424989E1);
    
    const __m128d S3        = _mm_set1_pd(-2.194779531642920639778E1);
    const __m128d S2        = _mm_set1_pd(1.470656354026814941758E2);
    const __m128d S1        = _mm_set1_pd(-3.838770957603691357202E2);
    const __m128d S0        = _mm_set1_pd(3.424398657913078477438E2);
    
    __m128d sign;
    __m128d mask;
    __m128d xabs;
    __m128d zz,ww,z,q,w,y,zz2,ww2;
    __m128d PA,PB;
    __m128d QA,QB;
    __m128d RA,RB;
    __m128d SA,SB;
    __m128d nom,denom;
    
    sign  = _mm_andnot_pd(signmask,x);
    xabs  = _mm_and_pd(x,signmask);
    
    mask  = _mm_cmpgt_pd(xabs,limit1);
    
    zz    = _mm_sub_pd(one,xabs);
    ww    = _mm_mul_pd(xabs,xabs);
    zz2   = _mm_mul_pd(zz,zz);
    ww2   = _mm_mul_pd(ww,ww);
    
    /* R */
    RA    = _mm_mul_pd(R4,zz2);  
    RB    = _mm_mul_pd(R3,zz2);
    RA    = _mm_add_pd(RA,R2); 
    RB    = _mm_add_pd(RB,R1);
    RA    = _mm_mul_pd(RA,zz2); 
    RB    = _mm_mul_pd(RB,zz);
    RA    = _mm_add_pd(RA,R0);
    RA    = _mm_add_pd(RA,RB);
    
    /* S, SA = zz2 */
    SB    = _mm_mul_pd(S3,zz2); 
    SA    = _mm_add_pd(zz2,S2);  
    SB    = _mm_add_pd(SB,S1);
    SA    = _mm_mul_pd(SA,zz2);
    SB    = _mm_mul_pd(SB,zz);
    SA    = _mm_add_pd(SA,S0);
    SA    = _mm_add_pd(SA,SB);
    
    /* P */
    PA    = _mm_mul_pd(P5,ww2);  
    PB    = _mm_mul_pd(P4,ww2);
    PA    = _mm_add_pd(PA,P3); 
    PB    = _mm_add_pd(PB,P2);
    PA    = _mm_mul_pd(PA,ww2); 
    PB    = _mm_mul_pd(PB,ww2);
    PA    = _mm_add_pd(PA,P1); 
    PB    = _mm_add_pd(PB,P0);
    PA    = _mm_mul_pd(PA,ww);
    PA    = _mm_add_pd(PA,PB);
    
    /* Q, QA = ww2 */
    QB    = _mm_mul_pd(Q4,ww2);  
    QA    = _mm_add_pd(ww2,Q3);  
    QB    = _mm_add_pd(QB,Q2);
    QA    = _mm_mul_pd(QA,ww2); 
    QB    = _mm_mul_pd(QB,ww2);
    QA    = _mm_add_pd(QA,Q1); 
    QB    = _mm_add_pd(QB,Q0); 
    QA    = _mm_mul_pd(QA,ww);  
    QA    = _mm_add_pd(QA,QB);
    
    RA    = _mm_mul_pd(RA,zz);
    PA    = _mm_mul_pd(PA,ww);
    
    nom   = _mm_or_pd( _mm_and_pd(mask,RA) , _mm_andnot_pd(mask,PA) );
    denom = _mm_or_pd( _mm_and_pd(mask,SA) , _mm_andnot_pd(mask,QA) );
    
    q     = _mm_mul_pd( nom, gmx_mm_inv_pd(denom) );
    
    zz    = _mm_add_pd(zz,zz);
    zz    = gmx_mm_sqrt_pd(zz);
    z     = _mm_sub_pd(quarterpi,zz);
    zz    = _mm_mul_pd(zz,q);
    zz    = _mm_sub_pd(zz,morebits);
    z     = _mm_sub_pd(z,zz);
    z     = _mm_add_pd(z,quarterpi);
    
    w     = _mm_mul_pd(xabs,q);
    w     = _mm_add_pd(w,xabs);
    
    z     = _mm_or_pd( _mm_and_pd(mask,z) , _mm_andnot_pd(mask,w) );
    
    mask  = _mm_cmpgt_pd(xabs,limit2);
    z     = _mm_or_pd( _mm_and_pd(mask,z) , _mm_andnot_pd(mask,xabs) );
    
    z = _mm_xor_pd(z,sign);
    
    return z;
}


static __m128d
gmx_mm_acos_pd(__m128d x)
{
    const __m128d signmask  = gmx_mm_castsi128_pd( _mm_set_epi32(0x7FFFFFFF,0xFFFFFFFF,0x7FFFFFFF,0xFFFFFFFF) );
    const __m128d one_pd    = _mm_set1_pd(1.0);
    const __m128d half_pd   = _mm_set1_pd(0.5);
    const __m128d pi_pd     = _mm_set1_pd(M_PI);
    const __m128d halfpi_pd = _mm_set1_pd(M_PI/2.0);
    
    __m128d mask1;
    __m128d mask2;
    __m128d xabs;
    __m128d z,z1,z2,z3;
    
    xabs  = _mm_and_pd(x,signmask);    
    mask1 = _mm_cmpgt_pd(xabs,half_pd);
    mask2 = _mm_cmpgt_pd(x,_mm_setzero_pd());
    
    z     = _mm_mul_pd(half_pd,_mm_sub_pd(one_pd,xabs));
    z     = _mm_mul_pd(z,gmx_mm_invsqrt_pd(z));
    z     = _mm_andnot_pd(_mm_cmpeq_pd(xabs,one_pd),z);
    
    z     = _mm_or_pd( _mm_and_pd(mask1,z) , _mm_andnot_pd(mask1,x) );
    z     = gmx_mm_asin_pd(z);
    
    z2    = _mm_add_pd(z,z);
    z1    = _mm_sub_pd(pi_pd,z2);
    z3    = _mm_sub_pd(halfpi_pd,z);    
    
    z     = _mm_or_pd( _mm_and_pd(mask2,z2) , _mm_andnot_pd(mask2,z1) );
    z     = _mm_or_pd( _mm_and_pd(mask1,z) , _mm_andnot_pd(mask1,z3) );
    
    return z;
}

static __m128d
gmx_mm_atan_pd(__m128d x)
{
    /* Same algorithm as cephes library */
    const __m128d signmask  = gmx_mm_castsi128_pd( _mm_set_epi32(0x7FFFFFFF,0xFFFFFFFF,0x7FFFFFFF,0xFFFFFFFF) );
    const __m128d limit1    = _mm_set1_pd(0.66);
    const __m128d limit2    = _mm_set1_pd(2.41421356237309504880);
    const __m128d quarterpi = _mm_set1_pd(M_PI/4.0);
    const __m128d halfpi    = _mm_set1_pd(M_PI/2.0);
    const __m128d mone      = _mm_set1_pd(-1.0);
    const __m128d morebits1 = _mm_set1_pd(0.5*6.123233995736765886130E-17);
    const __m128d morebits2 = _mm_set1_pd(6.123233995736765886130E-17);
    
    const __m128d P4        = _mm_set1_pd(-8.750608600031904122785E-1);
    const __m128d P3        = _mm_set1_pd(-1.615753718733365076637E1);
    const __m128d P2        = _mm_set1_pd(-7.500855792314704667340E1);
    const __m128d P1        = _mm_set1_pd(-1.228866684490136173410E2);
    const __m128d P0        = _mm_set1_pd(-6.485021904942025371773E1);
    
    const __m128d Q4        = _mm_set1_pd(2.485846490142306297962E1);
    const __m128d Q3        = _mm_set1_pd(1.650270098316988542046E2);
    const __m128d Q2        = _mm_set1_pd(4.328810604912902668951E2);
    const __m128d Q1        = _mm_set1_pd(4.853903996359136964868E2);
    const __m128d Q0        = _mm_set1_pd(1.945506571482613964425E2);
    
    __m128d sign;
    __m128d mask1,mask2;
    __m128d y,t1,t2;
    __m128d z,z2;
    __m128d P_A,P_B,Q_A,Q_B;
    
    sign   = _mm_andnot_pd(signmask,x);
    x      = _mm_and_pd(x,signmask);
    
    mask1  = _mm_cmpgt_pd(x,limit1);
    mask2  = _mm_cmpgt_pd(x,limit2);
    
    t1     = _mm_mul_pd(_mm_add_pd(x,mone),gmx_mm_inv_pd(_mm_sub_pd(x,mone)));
    t2     = _mm_mul_pd(mone,gmx_mm_inv_pd(x));
    
    y      = _mm_and_pd(mask1,quarterpi);
    y      = _mm_or_pd( _mm_and_pd(mask2,halfpi) , _mm_andnot_pd(mask2,y) );
    
    x      = _mm_or_pd( _mm_and_pd(mask1,t1) , _mm_andnot_pd(mask1,x) );
    x      = _mm_or_pd( _mm_and_pd(mask2,t2) , _mm_andnot_pd(mask2,x) );
    
    z      = _mm_mul_pd(x,x);
    z2     = _mm_mul_pd(z,z);
    
    P_A    = _mm_mul_pd(P4,z2);  
    P_B    = _mm_mul_pd(P3,z2);  
    P_A    = _mm_add_pd(P_A,P2); 
    P_B    = _mm_add_pd(P_B,P1); 
    P_A    = _mm_mul_pd(P_A,z2); 
    P_B    = _mm_mul_pd(P_B,z);  
    P_A    = _mm_add_pd(P_A,P0); 
    P_A    = _mm_add_pd(P_A,P_B);
    
    /* Q_A = z2 */
    Q_B    = _mm_mul_pd(Q4,z2);  
    Q_A    = _mm_add_pd(z2,Q3);  
    Q_B    = _mm_add_pd(Q_B,Q2); 
    Q_A    = _mm_mul_pd(Q_A,z2); 
    Q_B    = _mm_mul_pd(Q_B,z2); 
    Q_A    = _mm_add_pd(Q_A,Q1); 
    Q_B    = _mm_add_pd(Q_B,Q0); 
    Q_A    = _mm_mul_pd(Q_A,z);  
    Q_A    = _mm_add_pd(Q_A,Q_B);
    
    z      = _mm_mul_pd(z,P_A);
    z      = _mm_mul_pd(z,gmx_mm_inv_pd(Q_A));
    z      = _mm_mul_pd(z,x);
    z      = _mm_add_pd(z,x);
    
    t1     = _mm_and_pd(mask1,morebits1);
    t1     = _mm_or_pd( _mm_and_pd(mask2,morebits2) , _mm_andnot_pd(mask2,t1) );
    
    z      = _mm_add_pd(z,t1);
    y      = _mm_add_pd(y,z);
    
    y      = _mm_xor_pd(y,sign);
    
    return y;
}

static __m128d
gmx_mm_atan2_pd(__m128d y, __m128d x)
{
    const __m128d pi          = _mm_set1_pd(M_PI);
    const __m128d minuspi     = _mm_set1_pd(-M_PI);
    const __m128d halfpi      = _mm_set1_pd(M_PI/2.0);
    const __m128d minushalfpi = _mm_set1_pd(-M_PI/2.0);
    
    __m128d z,z1,z3,z4;
    __m128d w;
    __m128d maskx_lt,maskx_eq;
    __m128d masky_lt,masky_eq;
    __m128d mask1,mask2,mask3,mask4,maskall;
    
    maskx_lt  = _mm_cmplt_pd(x,_mm_setzero_pd());
    masky_lt  = _mm_cmplt_pd(y,_mm_setzero_pd());
    maskx_eq  = _mm_cmpeq_pd(x,_mm_setzero_pd());
    masky_eq  = _mm_cmpeq_pd(y,_mm_setzero_pd());
    
    z         = _mm_mul_pd(y,gmx_mm_inv_pd(x));
    z         = gmx_mm_atan_pd(z);
    
    mask1     = _mm_and_pd(maskx_eq,masky_lt);
    mask2     = _mm_andnot_pd(maskx_lt,masky_eq);
    mask3     = _mm_andnot_pd( _mm_or_pd(masky_lt,masky_eq) , maskx_eq);
    mask4     = _mm_and_pd(masky_eq,maskx_lt);
    
    maskall   = _mm_or_pd( _mm_or_pd(mask1,mask2), _mm_or_pd(mask3,mask4) );
    
    z         = _mm_andnot_pd(maskall,z);
    z1        = _mm_and_pd(mask1,minushalfpi);
    z3        = _mm_and_pd(mask3,halfpi);
    z4        = _mm_and_pd(mask4,pi);
    
    z         = _mm_or_pd( _mm_or_pd(z,z1), _mm_or_pd(z3,z4) );
    
    mask1     = _mm_andnot_pd(masky_lt,maskx_lt);
    mask2     = _mm_and_pd(maskx_lt,masky_lt);
    
    w         = _mm_or_pd( _mm_and_pd(mask1,pi), _mm_and_pd(mask2,minuspi) );    
    w         = _mm_andnot_pd(maskall,w);
    
    z         = _mm_add_pd(z,w);
    
    return z;
}


#define GMX_MM_TRANSPOSE2_PD(row0, row1) {            \
     __m128d __gmx_t1 = row0;                         \
     row0           = _mm_unpacklo_pd(row0,row1);     \
     row1           = _mm_unpackhi_pd(__gmx_t1,row1); \
}


#define GMX_MM_LOAD_1JCOORD_1ATOM_PD(ptr1,jx1,jy1,jz1) {               \
	 jx1            = _mm_loadu_pd(ptr1);                              \
     jy1            = _mm_unpackhi_pd(jx1,jx1);                        \
     jz1            = _mm_load_sd(ptr1+2);                             \
}


#define GMX_MM_LOAD_1JCOORD_2ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) { \
	 jx1            = _mm_loadu_pd(ptr1);                             \
     jy1            = _mm_unpackhi_pd(jx1,jx1);                       \
     jz1            = _mm_loadu_pd(ptr1+2);                           \
     jx2            = _mm_unpackhi_pd(jz1,jz1);                       \
     jy2            = _mm_loadu_pd(ptr1+4);                           \
     jz2            = _mm_unpackhi_pd(jy2,jy2);                       \
}


#define GMX_MM_LOAD_1JCOORD_3ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
	 jx1            = _mm_loadu_pd(ptr1);                             \
     jy1            = _mm_unpackhi_pd(jx1,jx1);                       \
     jz1            = _mm_loadu_pd(ptr1+2);                           \
     jx2            = _mm_unpackhi_pd(jz1,jz1);                       \
     jy2            = _mm_loadu_pd(ptr1+4);                           \
     jz2            = _mm_unpackhi_pd(jy2,jy2);                       \
	 jx3            = _mm_loadu_pd(ptr1+6);                           \
     jy3            = _mm_unpackhi_pd(jx3,jx3);                       \
     jz3            = _mm_load_sd(ptr1+8);                            \
}


#define GMX_MM_LOAD_1JCOORD_4ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
	 jx1            = _mm_loadu_pd(ptr1);                             \
     jy1            = _mm_unpackhi_pd(jx1,jx1);                       \
     jz1            = _mm_loadu_pd(ptr1+2);                           \
     jx2            = _mm_unpackhi_pd(jz1,jz1);                       \
     jy2            = _mm_loadu_pd(ptr1+4);                           \
     jz2            = _mm_unpackhi_pd(jy2,jy2);                       \
	 jx3            = _mm_loadu_pd(ptr1+6);                           \
     jy3            = _mm_unpackhi_pd(jx3,jx3);                       \
     jz3            = _mm_loadu_pd(ptr1+8);                           \
     jx4            = _mm_unpackhi_pd(jz3,jz3);                       \
     jy4            = _mm_loadu_pd(ptr1+10);                          \
     jz4            = _mm_unpackhi_pd(jy4,jy4);                       \
}

#define GMX_MM_LOAD_2JCOORD_1ATOM_PD(ptr1,ptr2,jx1,jy1,jz1) {   \
     __m128d _tmp1;                                             \
	 _tmp1           = _mm_loadu_pd(ptr1);                      \
     jy1             = _mm_loadu_pd(ptr2);                      \
     jz1             = _mm_load_sd(ptr1+2);                     \
     jx1             = _mm_unpacklo_pd(_tmp1,jy1);              \
     jy1             = _mm_unpackhi_pd(_tmp1,jy1);              \
     jz1             = _mm_loadh_pd(jz1,ptr2+2);                \
}


#define GMX_MM_LOAD_2JCOORD_2ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) {  \
     __m128d _tmp1, _tmp2,_tmp3;                                            \
	 _tmp1          = _mm_loadu_pd(ptr1);                                   \
	 jy1            = _mm_loadu_pd(ptr2);                                   \
	 _tmp2          = _mm_loadu_pd(ptr1+2);                                 \
	 jx2            = _mm_loadu_pd(ptr2+2);                                 \
	 _tmp3          = _mm_loadu_pd(ptr1+4);                                 \
	 jz2            = _mm_loadu_pd(ptr2+4);                                 \
     jx1            = _mm_unpacklo_pd(_tmp1,jy1);                           \
     jy1            = _mm_unpackhi_pd(_tmp1,jy1);                           \
     jz1            = _mm_unpacklo_pd(_tmp2,jx2);                           \
     jx2            = _mm_unpackhi_pd(_tmp2,jx2);                           \
     jy2            = _mm_unpacklo_pd(_tmp3,jz2);                           \
     jz2            = _mm_unpackhi_pd(_tmp3,jz2);                           \
}


#define GMX_MM_LOAD_2JCOORD_3ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128d _tmp1, _tmp2, _tmp3, _tmp4, _tmp5;                             \
	 _tmp1          = _mm_loadu_pd(ptr1);                                   \
	 jy1            = _mm_loadu_pd(ptr2);                                   \
	 _tmp2          = _mm_loadu_pd(ptr1+2);                                 \
	 jx2            = _mm_loadu_pd(ptr2+2);                                 \
	 _tmp3          = _mm_loadu_pd(ptr1+4);                                 \
	 jz2            = _mm_loadu_pd(ptr2+4);                                 \
	 _tmp4          = _mm_loadu_pd(ptr1+6);                                 \
	 jy3            = _mm_loadu_pd(ptr2+6);                                 \
	 jz3            = _mm_load_sd(ptr1+8);                                  \
	 _tmp5          = _mm_load_sd(ptr2+8);                                  \
     jx1            = _mm_unpacklo_pd(_tmp1,jy1);                           \
     jy1            = _mm_unpackhi_pd(_tmp1,jy1);                           \
     jz1            = _mm_unpacklo_pd(_tmp2,jx2);                           \
     jx2            = _mm_unpackhi_pd(_tmp2,jx2);                           \
     jy2            = _mm_unpacklo_pd(_tmp3,jz2);                           \
     jz2            = _mm_unpackhi_pd(_tmp3,jz2);                           \
     jx3            = _mm_unpacklo_pd(_tmp4,jy3);                           \
     jy3            = _mm_unpackhi_pd(_tmp4,jy3);                           \
     jz3            = _mm_unpacklo_pd(jz2,_tmp5);                           \
}


#define GMX_MM_LOAD_2JCOORD_4ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128d _tmp1, _tmp2,_tmp3, _tmp4, _tmp5, _tmp6;                       \
	 _tmp1          = _mm_loadu_pd(ptr1);                                   \
	 jy1            = _mm_loadu_pd(ptr2);                                   \
	 _tmp2          = _mm_loadu_pd(ptr1+2);                                 \
	 jx2            = _mm_loadu_pd(ptr2+2);                                 \
	 _tmp3          = _mm_loadu_pd(ptr1+4);                                 \
	 jz2            = _mm_loadu_pd(ptr2+4);                                 \
	 _tmp4          = _mm_loadu_pd(ptr1+6);                                 \
	 jy3            = _mm_loadu_pd(ptr2+6);                                 \
	 _tmp5          = _mm_loadu_pd(ptr1+8);                                 \
	 jx4            = _mm_loadu_pd(ptr2+8);                                 \
	 _tmp6          = _mm_loadu_pd(ptr1+10);                                \
	 jz4            = _mm_loadu_pd(ptr2+10);                                \
     jx1            = _mm_unpacklo_pd(_tmp1,jy1);                           \
     jy1            = _mm_unpackhi_pd(_tmp1,jy1);                           \
     jz1            = _mm_unpacklo_pd(_tmp2,jx2);                           \
     jx2            = _mm_unpackhi_pd(_tmp2,jx2);                           \
     jy2            = _mm_unpacklo_pd(_tmp3,jz2);                           \
     jz2            = _mm_unpackhi_pd(_tmp3,jz2);                           \
     jx3            = _mm_unpacklo_pd(_tmp4,jy3);                           \
     jy3            = _mm_unpackhi_pd(_tmp4,jy3);                           \
     jz3            = _mm_unpacklo_pd(_tmp5,jx4);                           \
     jx4            = _mm_unpackhi_pd(_tmp5,jx4);                           \
     jy4            = _mm_unpacklo_pd(_tmp6,jz4);                           \
     jz4            = _mm_unpackhi_pd(_tmp6,jz4);                           \
}


#define GMX_MM_UPDATE_1JCOORD_1ATOM_PD(ptr1,jx1,jy1,jz1) {            \
     jx1            = _mm_unpacklo_pd(jx1,jy1);                       \
     _mm_storeu_pd(ptr1, _mm_add_pd( _mm_loadu_pd(ptr1), jx1 ));      \
     _mm_store_sd(ptr1+2, _mm_add_sd( _mm_load_sd(ptr1+2), jz1));     \
}


#define GMX_MM_UPDATE_1JCOORD_2ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) { \
     jx1            = _mm_unpacklo_pd(jx1,jy1);                         \
     jz1            = _mm_unpacklo_pd(jz1,jx2);                         \
     jy2            = _mm_unpacklo_pd(jy2,jz2);                         \
     _mm_storeu_pd(ptr1, _mm_add_pd( _mm_loadu_pd(ptr1), jx1 ));        \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), jz1 ));    \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), jy2 ));    \
}


#define GMX_MM_UPDATE_1JCOORD_3ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     jx1            = _mm_unpacklo_pd(jx1,jy1);                         \
     jz1            = _mm_unpacklo_pd(jz1,jx2);                         \
     jy2            = _mm_unpacklo_pd(jy2,jz2);                         \
     jx3            = _mm_unpacklo_pd(jx3,jy3);                         \
     _mm_storeu_pd(ptr1, _mm_add_pd( _mm_loadu_pd(ptr1), jx1 ));        \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), jz1 ));    \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), jy2 ));    \
     _mm_storeu_pd(ptr1+6, _mm_add_pd( _mm_loadu_pd(ptr1+6), jx3 ));    \
     _mm_store_sd(ptr1+8, _mm_add_pd( _mm_load_sd(ptr1+8), jz3 ));      \
}


#define GMX_MM_UPDATE_1JCOORD_4ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     jx1            = _mm_unpacklo_pd(jx1,jy1);                         \
     jz1            = _mm_unpacklo_pd(jz1,jx2);                         \
     jy2            = _mm_unpacklo_pd(jy2,jz2);                         \
     jx3            = _mm_unpacklo_pd(jx3,jy3);                         \
     jz3            = _mm_unpacklo_pd(jz3,jx4);                         \
     jy4            = _mm_unpacklo_pd(jy4,jz4);                         \
     _mm_storeu_pd(ptr1,    _mm_add_pd( _mm_loadu_pd(ptr1), jx1 ));     \
     _mm_storeu_pd(ptr1+2,  _mm_add_pd( _mm_loadu_pd(ptr1+2), jz1 ));   \
     _mm_storeu_pd(ptr1+4,  _mm_add_pd( _mm_loadu_pd(ptr1+4), jy2 ));   \
     _mm_storeu_pd(ptr1+6,  _mm_add_pd( _mm_loadu_pd(ptr1+6), jx3 ));   \
     _mm_storeu_pd(ptr1+8,  _mm_add_pd( _mm_loadu_pd(ptr1+8), jz3 ));   \
     _mm_storeu_pd(ptr1+10, _mm_add_pd( _mm_loadu_pd(ptr1+10), jy4 ));  \
}


#define GMX_MM_UPDATE_2JCOORD_1ATOM_PD(ptr1,ptr2,jx1,jy1,jz1) {         \
     __m128d _tmp1,_tmp2,_tmp3;                                         \
     _tmp1          = _mm_unpacklo_pd(jx1,jy1);                         \
     _tmp2          = _mm_unpackhi_pd(jx1,jy1);                         \
     _tmp3          = _mm_unpackhi_pd(jz1,jz1);                         \
     _mm_storeu_pd(ptr1,  _mm_add_pd( _mm_loadu_pd(ptr1), _tmp1 ));     \
     _mm_storeu_pd(ptr2,  _mm_add_pd( _mm_loadu_pd(ptr2), _tmp2 ));     \
     _mm_store_sd(ptr1+2, _mm_add_pd( _mm_load_sd(ptr1+2), jz1 ));      \
     _mm_store_sd(ptr2+2, _mm_add_sd( _mm_load_sd(ptr2+2), _tmp3 ));    \
}


#define GMX_MM_UPDATE_2JCOORD_2ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) {  \
     __m128d _tmp1,_tmp2,_tmp3;                                         \
     _tmp1          = _mm_unpacklo_pd(jx1,jy1);                         \
     jy1            = _mm_unpackhi_pd(jx1,jy1);                         \
     _tmp2          = _mm_unpacklo_pd(jz1,jx2);                         \
     jx2            = _mm_unpackhi_pd(jz1,jx2);                         \
     _tmp3          = _mm_unpacklo_pd(jy2,jz2);                         \
     jz2            = _mm_unpackhi_pd(jy2,jz2);                         \
     _mm_storeu_pd(ptr1,   _mm_add_pd( _mm_loadu_pd(ptr1),   _tmp1 ));  \
     _mm_storeu_pd(ptr2,   _mm_add_pd( _mm_loadu_pd(ptr2),     jy1 ));  \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), _tmp2 ));  \
     _mm_storeu_pd(ptr2+2, _mm_add_pd( _mm_loadu_pd(ptr2+2),   jx2 ));  \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), _tmp3 ));  \
     _mm_storeu_pd(ptr2+4, _mm_add_pd( _mm_loadu_pd(ptr2+4),   jz2 ));  \
}


#define GMX_MM_UPDATE_2JCOORD_3ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128d _tmp1,_tmp2,_tmp3,_tmp4,_tmp5;                             \
     _tmp1          = _mm_unpacklo_pd(jx1,jy1);                         \
     jy1            = _mm_unpackhi_pd(jx1,jy1);                         \
     _tmp2          = _mm_unpacklo_pd(jz1,jx2);                         \
     jx2            = _mm_unpackhi_pd(jz1,jx2);                         \
     _tmp3          = _mm_unpacklo_pd(jy2,jz2);                         \
     jz2            = _mm_unpackhi_pd(jy2,jz2);                         \
     _tmp4          = _mm_unpacklo_pd(jx3,jy3);                         \
     jy3            = _mm_unpackhi_pd(jx3,jy3);                         \
     _tmp5          = _mm_unpackhi_pd(jz3,jz3);                         \
     _mm_storeu_pd(ptr1,   _mm_add_pd( _mm_loadu_pd(ptr1),   _tmp1 ));  \
     _mm_storeu_pd(ptr2,   _mm_add_pd( _mm_loadu_pd(ptr2),     jy1 ));  \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), _tmp2 ));  \
     _mm_storeu_pd(ptr2+2, _mm_add_pd( _mm_loadu_pd(ptr2+2),   jx2 ));  \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), _tmp3 ));  \
     _mm_storeu_pd(ptr2+4, _mm_add_pd( _mm_loadu_pd(ptr2+4),   jz2 ));  \
     _mm_storeu_pd(ptr1+6, _mm_add_pd( _mm_loadu_pd(ptr1+6), _tmp4 ));  \
     _mm_storeu_pd(ptr2+6, _mm_add_pd( _mm_loadu_pd(ptr2+6),   jy3 ));  \
     _mm_store_sd(ptr1+8, _mm_add_sd( _mm_load_sd(ptr1+8),   jz3 ));    \
     _mm_store_sd(ptr2+8, _mm_add_sd( _mm_load_sd(ptr2+8), _tmp5 ));    \
}


#define GMX_MM_UPDATE_2JCOORD_4ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128d _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6;                       \
     _tmp1          = _mm_unpacklo_pd(jx1,jy1);                         \
     jy1            = _mm_unpackhi_pd(jx1,jy1);                         \
     _tmp2          = _mm_unpacklo_pd(jz1,jx2);                         \
     jx2            = _mm_unpackhi_pd(jz1,jx2);                         \
     _tmp3          = _mm_unpacklo_pd(jy2,jz2);                         \
     jz2            = _mm_unpackhi_pd(jy2,jz2);                         \
     _tmp4          = _mm_unpacklo_pd(jx3,jy3);                         \
     jy3            = _mm_unpackhi_pd(jx3,jy3);                         \
     _tmp5          = _mm_unpacklo_pd(jz3,jx4);                         \
     jx4            = _mm_unpackhi_pd(jz3,jx4);                         \
     _tmp6          = _mm_unpacklo_pd(jy4,jz4);                         \
     jz4            = _mm_unpackhi_pd(jy4,jz4);                         \
     _mm_storeu_pd(ptr1,   _mm_add_pd( _mm_loadu_pd(ptr1),   _tmp1 ));  \
     _mm_storeu_pd(ptr2,   _mm_add_pd( _mm_loadu_pd(ptr2),     jy1 ));  \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), _tmp2 ));  \
     _mm_storeu_pd(ptr2+2, _mm_add_pd( _mm_loadu_pd(ptr2+2),   jx2 ));  \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), _tmp3 ));  \
     _mm_storeu_pd(ptr2+4, _mm_add_pd( _mm_loadu_pd(ptr2+4),   jz2 ));  \
     _mm_storeu_pd(ptr1+6, _mm_add_pd( _mm_loadu_pd(ptr1+6), _tmp4 ));  \
     _mm_storeu_pd(ptr2+6, _mm_add_pd( _mm_loadu_pd(ptr2+6),   jy3 ));  \
     _mm_storeu_pd(ptr1+8, _mm_add_pd( _mm_loadu_pd(ptr1+8), _tmp5 ));  \
     _mm_storeu_pd(ptr2+8, _mm_add_pd( _mm_loadu_pd(ptr2+8),   jx4 ));  \
     _mm_storeu_pd(ptr1+10,_mm_add_pd( _mm_loadu_pd(ptr1+10),_tmp6 ));  \
     _mm_storeu_pd(ptr2+10,_mm_add_pd( _mm_loadu_pd(ptr2+10),  jz4 ));  \
}


static inline __m128d
gmx_mm_scalarprod_pd(__m128d x, __m128d y, __m128d z)
{
	return (__m128d) _mm_add_pd(_mm_add_pd(_mm_mul_pd(x,x),_mm_mul_pd(y,y)),_mm_mul_pd(z,z));
}

static inline __m128d
gmx_mm_calc_rsq_pd(__m128d dx, __m128d dy, __m128d dz)
{
    return _mm_add_pd( _mm_add_pd( _mm_mul_pd(dx,dx), _mm_mul_pd(dy,dy) ), _mm_mul_pd(dz,dz) );
}



/* Routine to be called with rswitch/rcut at the beginning of a kernel
 * to set up the 7 constants used for analytic 5th order switch calculations.
 */
static inline void
gmx_mm_setup_switch5_constants_pd(__m128d rswitch, __m128d rcut,
								  __m128d *switch_C3, __m128d *switch_C4, __m128d *switch_C5,
								  __m128d *switch_D2, __m128d *switch_D3, __m128d *switch_D4)
{
	const __m128d  cm6  = (const __m128d) { -6.0, -6.0};
	const __m128d cm10  = (const __m128d) {-10.0,-10.0};
	const __m128d  c15  = (const __m128d) { 15.0, 15.0};
	const __m128d cm30  = (const __m128d) {-30.0,-30.0};
	const __m128d  c60  = (const __m128d) { 60.0, 60.0};

	__m128d d,dinv,dinv2,dinv3,dinv4,dinv5;
	
	d       = _mm_sub_pd(rcut,rswitch);
	dinv    = gmx_mm_inv_pd(d);
	dinv2   = _mm_mul_pd(dinv,dinv);
	dinv3   = _mm_mul_pd(dinv2,dinv);
	dinv4   = _mm_mul_pd(dinv2,dinv2);
	dinv5   = _mm_mul_pd(dinv3,dinv2);
	
	*switch_C3 = _mm_mul_pd(cm10,dinv3);
	*switch_C4 = _mm_mul_pd(c15,dinv4);
	*switch_C5 = _mm_mul_pd(cm6,dinv5);
	*switch_D2 = _mm_mul_pd(cm30,dinv3);
	*switch_D3 = _mm_mul_pd(c60,dinv4);
	*switch_D4 = _mm_mul_pd(cm30,dinv5);
}


#define GMX_MM_SET_SWITCH5_PD(r,rswitch,rcut,sw,dsw,sw_C3,sw_C4,sw_C5,sw_D2,sw_D3,sw_D4) {   \
    const __m128d  _sw_one  = (const __m128d) {1.0,1.0};                                     \
    __m128d d,d2;                                                                            \
    d     = _mm_max_pd(r,rswitch);                                                           \
    d     = _mm_min_pd(d,rcut);                                                              \
    d     = _mm_sub_pd(d,rswitch);                                                           \
    d2    = _mm_mul_pd(d,d);                                                                 \
    sw    = _mm_mul_pd(d,sw_C5);                                                             \
    dsw   = _mm_mul_pd(d,sw_D4);                                                             \
    sw    = _mm_add_pd(sw,sw_C4);                                                            \
    dsw   = _mm_add_pd(dsw,sw_D3);                                                           \
    sw    = _mm_mul_pd(sw,d);                                                                \
    dsw   = _mm_mul_pd(dsw,d);                                                               \
    sw    = _mm_add_pd(sw,sw_C3);                                                            \
    dsw   = _mm_add_pd(dsw,sw_D2);                                                           \
    sw    = _mm_mul_pd(sw,_mm_mul_pd(d,d2));                                                 \
    dsw   = _mm_mul_pd(dsw,d2);                                                              \
    sw    = _mm_add_pd(sw,_sw_one);                                                          \
}


/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128d
gmx_mm_int_coulomb_pd(__m128d rinv, __m128d qq,__m128d *vctot)
{
	__m128d vcoul = _mm_mul_pd(qq,rinv);
	*vctot   = _mm_add_pd(*vctot,vcoul);
	return (__m128d) vcoul;
}


static inline void
gmx_mm_int_coulomb_noforce_pd(__m128d rinv, __m128d qq,__m128d *vctot)
{
	__m128d vcoul = _mm_mul_pd(qq,rinv);
	*vctot   = _mm_add_pd(*vctot,vcoul);
	return;
}

/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128d
gmx_mm_int_coulombrf_pd(__m128d rinv, __m128d rsq, __m128d krf, __m128d crf, __m128d qq,__m128d *vctot)
{
	const __m128d two  = (const __m128d) {2.0,2.0};
	__m128d vcoul,krsq;
	
	krsq   = _mm_mul_pd(krf,rsq);
	vcoul  = _mm_mul_pd(qq, _mm_sub_pd(_mm_add_pd(rinv,krsq),crf));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	return (__m128d) _mm_mul_pd(qq, _mm_sub_pd(rinv, _mm_mul_pd(two,krsq)));
}


static inline void
gmx_mm_int_coulombrf_noforce_pd(__m128d rinv, __m128d rsq, __m128d krf, __m128d crf, __m128d qq,__m128d *vctot)
{
	__m128d vcoul,krsq;
	
	krsq   = _mm_mul_pd(krf,rsq);
	vcoul  = _mm_mul_pd(qq, _mm_sub_pd(_mm_add_pd(rinv,krsq),crf));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	return;
}



/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128d
gmx_mm_int_lj_pd(__m128d rinvsq, __m128d c6, __m128d c12, __m128d *vvdwtot)
{
	const __m128d six    = (const __m128d) {6.0,6.0};
	const __m128d twelve = (const __m128d) {12.0,12.0};
	
	__m128d rinvsix,vvdw6,vvdw12;
		
	rinvsix  = _mm_mul_pd(_mm_mul_pd(rinvsq,rinvsq),rinvsq);
	vvdw6    = _mm_mul_pd(c6,rinvsix);  
	vvdw12   = _mm_mul_pd(c12, _mm_mul_pd(rinvsix,rinvsix));
	*vvdwtot = _mm_add_pd(*vvdwtot , _mm_sub_pd(vvdw12,vvdw6));
	
	return (__m128d) _mm_sub_pd( _mm_mul_pd(twelve,vvdw12),_mm_mul_pd(six,vvdw6));
}
		   

static inline void
gmx_mm_int_lj_potonly_pd(__m128d rinvsq, __m128d c6, __m128d c12, __m128d *vvdwtot)
{
	__m128d rinvsix,vvdw6,vvdw12;
	
	rinvsix  = _mm_mul_pd(_mm_mul_pd(rinvsq,rinvsq),rinvsq);
	vvdw6    = _mm_mul_pd(c6,rinvsix);  
	vvdw12   = _mm_mul_pd(c12, _mm_mul_pd(rinvsix,rinvsix));
	*vvdwtot = _mm_add_pd(*vvdwtot , _mm_sub_pd(vvdw12,vvdw6));
	
	return;
}

#ifdef GMX_SSE4
#  define gmx_mm_extract_epi32(x, imm) _mm_extract_epi32(x,imm)
#else
#  define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))
#endif



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_2_table_coulomb_pd(__m128d r, __m128d tabscale, double * VFtab, __m128d qq, __m128d *vctot)
{
    __m128d  rt,eps,eps2,Y,F,G,H,vcoul;
	__m128i  n0;
	int     n_a,n_b;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	Y        = _mm_load_pd(VFtab + 4* n_a);
	F        = _mm_load_pd(VFtab + 4* n_b);
	G        = _mm_load_pd(VFtab + 4* n_a + 2);
	H        = _mm_load_pd(VFtab + 4* n_b + 2);
	GMX_MM_TRANSPOSE2_PD(Y,F);
	GMX_MM_TRANSPOSE2_PD(G,H);
	H        = _mm_mul_pd(H,eps2);              /* Heps2 */
	G        = _mm_mul_pd(G,eps);               /* Geps  */
	F        = _mm_add_pd(F, _mm_add_pd(G,H));  /* Fp    */
	vcoul    = _mm_mul_pd(qq, _mm_add_pd(Y, _mm_mul_pd(eps,F)));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	F        = _mm_mul_pd(qq, _mm_add_pd(F, _mm_add_pd(G, _mm_add_pd(H,H))));
	
	return (__m128d) _mm_mul_pd(F,tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_2_table_lj_pd(__m128d r, __m128d tabscale, double * VFtab, int offset, __m128d c6, __m128d c12, __m128d *vvdwtot)
{
    __m128d  rt,eps,eps2,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i  n0;
	int     n_a,n_b;
	
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	
	/* For a few cases, like TIP4p waters, there are particles with LJ-only interactions in a loop where
	 * the table data might contain both coulomb and LJ. To handle this case, we use an offset value of 0
	 * if the data is an LJ-only table, and 1 if it is actually a mixed coul+lj table.
	 */
	Yd       = _mm_load_pd(VFtab + 4*(offset+2)* n_a + 4*offset);
	Fd       = _mm_load_pd(VFtab + 4*(offset+2)* n_b + 4*offset);
	Gd       = _mm_load_pd(VFtab + 4*(offset+2)* n_a + 4*offset + 2);
	Hd       = _mm_load_pd(VFtab + 4*(offset+2)* n_b + 4*offset + 2);
	Yr       = _mm_load_pd(VFtab + 4*(offset+2)* n_a + 4*offset + 4);
	Fr       = _mm_load_pd(VFtab + 4*(offset+2)* n_b + 4*offset + 4);
	Gr       = _mm_load_pd(VFtab + 4*(offset+2)* n_a + 4*offset + 6);
	Hr       = _mm_load_pd(VFtab + 4*(offset+2)* n_b + 4*offset + 6);
	GMX_MM_TRANSPOSE2_PD(Yd,Fd);
	GMX_MM_TRANSPOSE2_PD(Gd,Hd);
	GMX_MM_TRANSPOSE2_PD(Yr,Fr);
	GMX_MM_TRANSPOSE2_PD(Gr,Hr);
	Hd       = _mm_mul_pd(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_pd(Gd,eps);               /* Geps  */
	Fd       = _mm_add_pd(Fd, _mm_add_pd(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_pd(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_pd(Gr,eps);               /* Geps  */
	Fr       = _mm_add_pd(Fr, _mm_add_pd(Gr,Hr));  /* Fp    */
	vvdw6    = _mm_mul_pd(c6,  _mm_add_pd(Yd, _mm_mul_pd(eps,Fd)));
	vvdw12   = _mm_mul_pd(c12, _mm_add_pd(Yr, _mm_mul_pd(eps,Fr)));
	*vvdwtot = _mm_add_pd(*vvdwtot, _mm_add_pd(vvdw6,vvdw12));
	
	Fd        = _mm_mul_pd(c6,  _mm_add_pd(Fd, _mm_add_pd(Gd, _mm_add_pd(Hd,Hd))));
	Fr        = _mm_mul_pd(c12, _mm_add_pd(Fr, _mm_add_pd(Gr, _mm_add_pd(Hr,Hr))));
	
	return (__m128d) _mm_mul_pd( _mm_add_pd(Fd,Fr),tabscale);
}


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_2_table_coulomb_and_lj_pd(__m128d r, __m128d tabscale, double * VFtab, __m128d qq, __m128d c6, __m128d c12, 
								  __m128d *vctot, __m128d *vvdwtot)
{
    __m128d  rt,eps,eps2,vcoul,Yc,Fc,Gc,Hc,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i  n0;
	int     n_a,n_b;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	
	Yc       = _mm_load_pd(VFtab + 12* n_a);
	Fc       = _mm_load_pd(VFtab + 12* n_b);
	Gc       = _mm_load_pd(VFtab + 12* n_a + 2);
	Hc       = _mm_load_pd(VFtab + 12* n_b + 2);
	Yd       = _mm_load_pd(VFtab + 12* n_a + 4);
	Fd       = _mm_load_pd(VFtab + 12* n_b + 4);
	Gd       = _mm_load_pd(VFtab + 12* n_a + 6);
	Hd       = _mm_load_pd(VFtab + 12* n_b + 6);
	Yr       = _mm_load_pd(VFtab + 12* n_a + 8);
	Fr       = _mm_load_pd(VFtab + 12* n_b + 8);
	Gr       = _mm_load_pd(VFtab + 12* n_a + 10);
	Hr       = _mm_load_pd(VFtab + 12* n_b + 10);
	GMX_MM_TRANSPOSE2_PD(Yc,Fc);
	GMX_MM_TRANSPOSE2_PD(Gc,Hc);
	GMX_MM_TRANSPOSE2_PD(Yd,Fd);
	GMX_MM_TRANSPOSE2_PD(Gd,Hd);
	GMX_MM_TRANSPOSE2_PD(Yr,Fr);
	GMX_MM_TRANSPOSE2_PD(Gr,Hr);
	Hc       = _mm_mul_pd(Hc,eps2);              /* Heps2 */
	Gc       = _mm_mul_pd(Gc,eps);               /* Geps  */
	Fc       = _mm_add_pd(Fc, _mm_add_pd(Gc,Hc));  /* Fp    */
	Hd       = _mm_mul_pd(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_pd(Gd,eps);               /* Geps  */
	Fd       = _mm_add_pd(Fd, _mm_add_pd(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_pd(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_pd(Gr,eps);               /* Geps  */
	Fr       = _mm_add_pd(Fr, _mm_add_pd(Gr,Hr));  /* Fp    */
	
	vcoul    = _mm_mul_pd(qq, _mm_add_pd(Yc, _mm_mul_pd(eps,Fc)));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	vvdw6    = _mm_mul_pd(c6,  _mm_add_pd(Yd, _mm_mul_pd(eps,Fd)));
	vvdw12   = _mm_mul_pd(c12, _mm_add_pd(Yr, _mm_mul_pd(eps,Fr)));
	*vvdwtot = _mm_add_pd(*vvdwtot, _mm_add_pd(vvdw6,vvdw12));
	
	Fc       = _mm_mul_pd(qq, _mm_add_pd(Fc, _mm_add_pd(Gc, _mm_add_pd(Hc,Hc))));
	Fd       = _mm_mul_pd(c6,  _mm_add_pd(Fd, _mm_add_pd(Gd, _mm_add_pd(Hd,Hd))));
	Fr       = _mm_mul_pd(c12, _mm_add_pd(Fr, _mm_add_pd(Gr, _mm_add_pd(Hr,Hr))));
	
	return (__m128d) _mm_mul_pd( _mm_add_pd(Fc,_mm_add_pd(Fd,Fr)),tabscale);
}




/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_1_table_coulomb_pd(__m128d r, __m128d tabscale, double * VFtab, __m128d qq, __m128d *vctot)
{
    __m128d  rt,eps,eps2,Y,F,G,H,vcoul;
	__m128i  n0;
	int     n_a;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	Y        = _mm_load_pd(VFtab + 4* n_a);
	F        = _mm_setzero_pd();
	G        = _mm_load_pd(VFtab + 4* n_a + 2);
	H        = _mm_setzero_pd();
	GMX_MM_TRANSPOSE2_PD(Y,F);
	GMX_MM_TRANSPOSE2_PD(G,H);

	H        = _mm_mul_pd(H,eps2);              /* Heps2 */
	G        = _mm_mul_pd(G,eps);               /* Geps  */
	F        = _mm_add_pd(F, _mm_add_pd(G,H));  /* Fp    */
	vcoul    = _mm_mul_pd(qq, _mm_add_pd(Y, _mm_mul_pd(eps,F)));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	F        = _mm_mul_pd(qq, _mm_add_pd(F, _mm_add_pd(G, _mm_add_pd(H,H))));
	
	return (__m128d) _mm_mul_pd(F,tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_1_table_lj_pd(__m128d r, __m128d tabscale, double * VFtab, int offset, __m128d c6, __m128d c12, __m128d *vvdwtot)
{
    __m128d  rt,eps,eps2,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i  n0;
	int     n_a;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	/* For a few cases, like TIP4p waters, there are particles with LJ-only interactions in a loop where
	 * the table data might contain both coulomb and LJ. To handle this case, we use an offset value of 0
	 * if the data is an LJ-only table, and 1 if it is actually a mixed coul+lj table.
	 */
	Yd       = _mm_load_pd(VFtab + 4*(offset+2)*n_a + 4*offset);
	Fd       = _mm_setzero_pd();
	Gd       = _mm_load_pd(VFtab + 4*(offset+2)*n_a + 4*offset + 2);
	Hd       = _mm_setzero_pd();
	Yr       = _mm_load_pd(VFtab + 4*(offset+2)*n_a + 4*offset + 4);
	Fr       = _mm_setzero_pd();
	Gr       = _mm_load_pd(VFtab + 4*(offset+2)*n_a + 4*offset + 6);
	Hr       = _mm_setzero_pd();
	GMX_MM_TRANSPOSE2_PD(Yd,Fd);
	GMX_MM_TRANSPOSE2_PD(Gd,Hd);
	GMX_MM_TRANSPOSE2_PD(Yr,Fr);
	GMX_MM_TRANSPOSE2_PD(Gr,Hr);
	
	Hd       = _mm_mul_pd(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_pd(Gd,eps);               /* Geps  */
	Fd       = _mm_add_pd(Fd, _mm_add_pd(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_pd(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_pd(Gr,eps);               /* Geps  */
	Fr       = _mm_add_pd(Fr, _mm_add_pd(Gr,Hr));  /* Fp    */
	vvdw6    = _mm_mul_pd(c6,  _mm_add_pd(Yd, _mm_mul_pd(eps,Fd)));
	vvdw12   = _mm_mul_pd(c12, _mm_add_pd(Yr, _mm_mul_pd(eps,Fr)));
	*vvdwtot = _mm_add_pd(*vvdwtot, _mm_add_pd(vvdw6,vvdw12));
	
	Fd        = _mm_mul_pd(c6,  _mm_add_pd(Fd, _mm_add_pd(Gd, _mm_add_pd(Hd,Hd))));
	Fr        = _mm_mul_pd(c12, _mm_add_pd(Fr, _mm_add_pd(Gr, _mm_add_pd(Hr,Hr))));
	
	return (__m128d) _mm_mul_pd( _mm_add_pd(Fd,Fr),tabscale);
}


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_1_table_coulomb_and_lj_pd(__m128d r, __m128d tabscale, double * VFtab, __m128d qq, __m128d c6, __m128d c12, 
									 __m128d *vctot, __m128d *vvdwtot)
{
    __m128d  rt,eps,eps2,vcoul,Yc,Fc,Gc,Hc,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i  n0;
	int     n_a;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	Yc       = _mm_load_pd(VFtab + 12* n_a);
	Fc       = _mm_setzero_pd();
	Gc       = _mm_load_pd(VFtab + 12* n_a + 2);
	Hc       = _mm_setzero_pd();
	Yd       = _mm_load_pd(VFtab + 12* n_a + 4);
	Fd       = _mm_setzero_pd();
	Gd       = _mm_load_pd(VFtab + 12* n_a + 6);
	Hd       = _mm_setzero_pd();
	Yr       = _mm_load_pd(VFtab + 12* n_a + 8);
	Fr       = _mm_setzero_pd();
	Gr       = _mm_load_pd(VFtab + 12* n_a + 10);
	Hr       = _mm_setzero_pd();
	GMX_MM_TRANSPOSE2_PD(Yc,Fc);
	GMX_MM_TRANSPOSE2_PD(Gc,Hc);
	GMX_MM_TRANSPOSE2_PD(Yd,Fd);
	GMX_MM_TRANSPOSE2_PD(Gd,Hd);
	GMX_MM_TRANSPOSE2_PD(Yr,Fr);
	GMX_MM_TRANSPOSE2_PD(Gr,Hr);	
	Hc       = _mm_mul_pd(Hc,eps2);              /* Heps2 */
	Gc       = _mm_mul_pd(Gc,eps);               /* Geps  */
	Fc       = _mm_add_pd(Fc, _mm_add_pd(Gc,Hc));  /* Fp    */
	Hd       = _mm_mul_pd(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_pd(Gd,eps);               /* Geps  */
	Fd       = _mm_add_pd(Fd, _mm_add_pd(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_pd(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_pd(Gr,eps);               /* Geps  */
	Fr       = _mm_add_pd(Fr, _mm_add_pd(Gr,Hr));  /* Fp    */
	
	vcoul    = _mm_mul_pd(qq, _mm_add_pd(Yc, _mm_mul_pd(eps,Fc)));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	vvdw6    = _mm_mul_pd(c6,  _mm_add_pd(Yd, _mm_mul_pd(eps,Fd)));
	vvdw12   = _mm_mul_pd(c12, _mm_add_pd(Yr, _mm_mul_pd(eps,Fr)));
	*vvdwtot = _mm_add_pd(*vvdwtot, _mm_add_pd(vvdw6,vvdw12));
	
	Fc       = _mm_mul_pd(qq, _mm_add_pd(Fc, _mm_add_pd(Gc, _mm_add_pd(Hc,Hc))));
	Fd       = _mm_mul_pd(c6,  _mm_add_pd(Fd, _mm_add_pd(Gd, _mm_add_pd(Hd,Hd))));
	Fr       = _mm_mul_pd(c12, _mm_add_pd(Fr, _mm_add_pd(Gr, _mm_add_pd(Hr,Hr))));
	
	return (__m128d) _mm_mul_pd( _mm_add_pd(Fc,_mm_add_pd(Fd,Fr)),tabscale);
}



/* Return force should be multiplied by +rinv to get fscal */
static inline __m128d
gmx_mm_int_2_genborn_pd(__m128d r, __m128d isai, 
						double * isaj1, double *isaj2, 
						__m128d gbtabscale, double * GBtab, __m128d qq, __m128d *dvdasum, 
						double *dvdaj1, double *dvdaj2,
						__m128d *vgbtot)
{
	const __m128d half  = (const __m128d) {0.5,0.5};
	
    __m128d  rt,eps,eps2,Y,F,G,H,VV,FF,ftmp,isaprod,t2,t3,t4,isaj,vgb,dvdatmp;
	__m128i  n0;
	int     n_a,n_b,n_c,n_d;
	
	/* Assemble isaj */
	isaj     = _mm_load_sd(isaj1);
	t2       = _mm_load_sd(isaj2);
	isaj     = _mm_unpacklo_pd(isaj,t2);  /* t2 t1 */
	
	isaprod     = _mm_mul_pd(isai,isaj);
	qq          = _mm_mul_pd(qq,isaprod);
	gbtabscale  = _mm_mul_pd( isaprod, gbtabscale );
	
	rt       = _mm_mul_pd(r,gbtabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	Y        = _mm_load_pd(GBtab + 4* n_a);
	F        = _mm_load_pd(GBtab + 4* n_b);
	G        = _mm_load_pd(GBtab + 4* n_a + 2);
	H        = _mm_load_pd(GBtab + 4* n_b + 2);
	GMX_MM_TRANSPOSE2_PD(Y,F);
	GMX_MM_TRANSPOSE2_PD(G,H);
	G        = _mm_mul_pd(G,eps);               /* Geps  */
	H        = _mm_mul_pd(H,eps2);              /* Heps2 */
	F        = _mm_add_pd(_mm_add_pd(F,G),H);  /* Fp    */
	
	VV       = _mm_add_pd(Y, _mm_mul_pd(eps,F));
	FF       = _mm_add_pd(_mm_add_pd(F,G), _mm_add_pd(H,H));
	
	vgb      = _mm_mul_pd(qq, VV);
	*vgbtot  = _mm_sub_pd(*vgbtot,vgb); /* Yes, the sign is correct */
	
	ftmp     = _mm_mul_pd(_mm_mul_pd(qq, FF), gbtabscale);
	
	dvdatmp  = _mm_mul_pd(half,_mm_add_pd(vgb,_mm_mul_pd(ftmp,r)));
	
	*dvdasum = _mm_add_pd(*dvdasum,dvdatmp);
	
	dvdatmp  = _mm_mul_pd(_mm_mul_pd(dvdatmp,isaj), isaj);
	
	/* Update 2 dada[j] values */
	Y        = _mm_load_sd(dvdaj1);
	F        = _mm_load_sd(dvdaj2);
	t2       = _mm_unpackhi_pd(dvdatmp,dvdatmp);
	
	_mm_store_sd( dvdaj1 , _mm_add_sd( Y, dvdatmp ) );
	_mm_store_sd( dvdaj2 , _mm_add_sd( F, t2 ) );
	
	return (__m128d) ftmp;
}

/* Return force should be multiplied by +rinv to get fscal */
static inline __m128d
gmx_mm_int_1_genborn_pd(__m128d r, __m128d isai, 
						double * isaj1, 
						__m128d gbtabscale, double * GBtab, __m128d qq, __m128d *dvdasum, 
						double *dvdaj1, 
						__m128d *vgbtot)
{
	const __m128d half  = (const __m128d) {0.5,0.5};
	
    __m128d  rt,eps,eps2,Y,F,G,H,VV,FF,ftmp,isaprod,t2,t3,t4,isaj,vgb,dvdatmp;
	__m128i  n0;
	int     n_a,n_b,n_c,n_d;
	
	/* Assemble isaj */
	isaj     = _mm_load_sd(isaj1);
	
	isaprod     = _mm_mul_pd(isai,isaj);
	qq          = _mm_mul_pd(qq,isaprod);
	gbtabscale  = _mm_mul_pd( isaprod, gbtabscale );
	
	rt       = _mm_mul_pd(r,gbtabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	Y        = _mm_load_pd(GBtab + 4* n_a);
	F        = _mm_setzero_pd();
	G        = _mm_load_pd(GBtab + 4* n_a + 2);
	H        = _mm_setzero_pd();
	GMX_MM_TRANSPOSE2_PD(Y,F);
	GMX_MM_TRANSPOSE2_PD(G,H);

	G        = _mm_mul_pd(G,eps);               /* Geps  */
	H        = _mm_mul_pd(H,eps2);              /* Heps2 */
	F        = _mm_add_pd(_mm_add_pd(F,G),H);  /* Fp    */
	
	VV       = _mm_add_pd(Y, _mm_mul_pd(eps,F));
	FF       = _mm_add_pd(_mm_add_pd(F,G), _mm_add_pd(H,H));
	
	vgb      = _mm_mul_pd(qq, VV);
	*vgbtot  = _mm_sub_pd(*vgbtot,vgb); /* Yes, the sign is correct */
	
	ftmp     = _mm_mul_pd(_mm_mul_pd(qq, FF), gbtabscale);
	
	dvdatmp  = _mm_mul_pd(half,_mm_add_pd(vgb,_mm_mul_pd(ftmp,r)));
	
	*dvdasum = _mm_add_pd(*dvdasum,dvdatmp);
	
	dvdatmp  = _mm_mul_pd(_mm_mul_pd(dvdatmp,isaj), isaj);
	
	/* Update 1 dada[j] values */
	Y        = _mm_load_sd(dvdaj1);
	
	_mm_store_sd( dvdaj1 , _mm_add_sd( Y, dvdatmp ) );
	
	return (__m128d) ftmp;
}




static inline void
gmx_mm_update_iforce_1atom_pd(__m128d fix1, __m128d fiy1, __m128d fiz1,
							  double *fptr,
							  double *fshiftptr)
{
	__m128d t1,t2,t3;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_pd(fix1,fiy1);   
	/* fiz1 is fine as it is */
#else
	/* SSE2 */
	/* transpose data */
	t1 = fix1;
	fix1 = _mm_unpacklo_pd(fix1,fiy1); /* y0 x0 */
	fiy1 = _mm_unpackhi_pd(t1,fiy1);   /* y1 x1 */
	
	fix1 = _mm_add_pd(fix1,fiy1);
	fiz1 = _mm_add_sd( fiz1, _mm_unpackhi_pd(fiz1,fiz1 ));
#endif
	_mm_storeu_pd( fptr, _mm_add_pd( _mm_loadu_pd(fptr), fix1 ));
	_mm_store_sd( fptr+2, _mm_add_sd( _mm_load_sd(fptr+2), fiz1 ));
	
	_mm_storeu_pd( fshiftptr, _mm_add_pd( _mm_loadu_pd(fshiftptr), fix1 ));
	_mm_store_sd( fshiftptr+2, _mm_add_sd( _mm_load_sd(fshiftptr+2), fiz1 ));
}

static inline void
gmx_mm_update_iforce_2atoms_pd(__m128d fix1, __m128d fiy1, __m128d fiz1,
							   __m128d fix2, __m128d fiy2, __m128d fiz2,
							   double *fptr,
							   double *fshiftptr)
{
	__m128d t1;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_pd(fix1,fiy1);   
	fiz1 = _mm_hadd_pd(fiz1,fix2);   
	fiy2 = _mm_hadd_pd(fiy2,fiz2);   
#else
	/* SSE2 */
	/* transpose data */
	GMX_MM_TRANSPOSE2_PD(fix1,fiy1);
	GMX_MM_TRANSPOSE2_PD(fiz1,fix2);
	GMX_MM_TRANSPOSE2_PD(fiy2,fiz2);
	
	fix1 = _mm_add_pd(fix1,fiy1);
	fiz1 = _mm_add_pd(fiz1,fix2);
	fiy2 = _mm_add_pd(fiy2,fiz2);
#endif
	_mm_storeu_pd( fptr, _mm_add_pd( _mm_loadu_pd(fptr), fix1 ));
	_mm_storeu_pd( fptr+2, _mm_add_pd( _mm_loadu_pd(fptr+2), fiz1 ));
	_mm_storeu_pd( fptr+4, _mm_add_pd( _mm_loadu_pd(fptr+4), fiy2 ));
	
	t1 = _mm_shuffle_pd(fiz1,fiy2,_MM_SHUFFLE2(0,1)); 
	fix1 = _mm_add_pd(fix1,t1); /* x and y sums */	
	fiz1 = _mm_add_sd(fiz1, _mm_unpackhi_pd(fiy2,fiy2)); /* z sum */

	_mm_storeu_pd( fshiftptr, _mm_add_pd( _mm_loadu_pd(fshiftptr), fix1 ));
	_mm_store_sd( fshiftptr+2, _mm_add_sd( _mm_load_sd(fshiftptr+2), fiz1 ));
}



static inline void
gmx_mm_update_iforce_3atoms_pd(__m128d fix1, __m128d fiy1, __m128d fiz1,
							   __m128d fix2, __m128d fiy2, __m128d fiz2,
							   __m128d fix3, __m128d fiy3, __m128d fiz3,
							   double *fptr,
							   double *fshiftptr)
{
	__m128d t1,t2;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_pd(fix1,fiy1);   
	fiz1 = _mm_hadd_pd(fiz1,fix2);   
	fiy2 = _mm_hadd_pd(fiy2,fiz2);   
	fix3 = _mm_hadd_pd(fix3,fiy3);   
	/* fiz3 is fine as it is */
#else
	/* SSE2 */
	/* transpose data */
	GMX_MM_TRANSPOSE2_PD(fix1,fiy1);
	GMX_MM_TRANSPOSE2_PD(fiz1,fix2);
	GMX_MM_TRANSPOSE2_PD(fiy2,fiz2);
	t1 = fix3;
	fix3 = _mm_unpacklo_pd(fix3,fiy3); /* y0 x0 */
	fiy3 = _mm_unpackhi_pd(t1,fiy3);   /* y1 x1 */
	
	fix1 = _mm_add_pd(fix1,fiy1);
	fiz1 = _mm_add_pd(fiz1,fix2);
	fiy2 = _mm_add_pd(fiy2,fiz2);
	
	fix3 = _mm_add_pd(fix3,fiy3);
	fiz3 = _mm_add_sd( fiz3, _mm_unpackhi_pd(fiz3,fiz3));
#endif
	_mm_storeu_pd( fptr, _mm_add_pd( _mm_loadu_pd(fptr), fix1 ));
	_mm_storeu_pd( fptr+2, _mm_add_pd( _mm_loadu_pd(fptr+2), fiz1 ));
	_mm_storeu_pd( fptr+4, _mm_add_pd( _mm_loadu_pd(fptr+4), fiy2 ));
	_mm_storeu_pd( fptr+6, _mm_add_pd( _mm_loadu_pd(fptr+6), fix3 ));
	_mm_store_sd( fptr+8, _mm_add_sd( _mm_load_sd(fptr+8), fiz3 ));
	
	fix1 = _mm_add_pd(fix1,fix3);
	t1   = _mm_shuffle_pd(fiz1,fiy2,_MM_SHUFFLE2(0,1)); 
	fix1 = _mm_add_pd(fix1,t1); /* x and y sums */	

	t2   = _mm_shuffle_pd(fiy2,fiy2,_MM_SHUFFLE2(1,1));
	fiz1 = _mm_add_sd(fiz1,fiz3);
	fiz1 = _mm_add_sd(fiz1,t2); /* z sum */
	
	_mm_storeu_pd( fshiftptr, _mm_add_pd( _mm_loadu_pd(fshiftptr), fix1 ));
	_mm_store_sd( fshiftptr+2, _mm_add_sd( _mm_load_sd(fshiftptr+2), fiz1 ));
}



static inline void
gmx_mm_update_iforce_4atoms_pd(__m128d fix1, __m128d fiy1, __m128d fiz1,
							   __m128d fix2, __m128d fiy2, __m128d fiz2,
							   __m128d fix3, __m128d fiy3, __m128d fiz3,
							   __m128d fix4, __m128d fiy4, __m128d fiz4,
							   double *fptr,
							   double *fshiftptr)
{
	__m128d t1,t2;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_pd(fix1,fiy1);   
	fiz1 = _mm_hadd_pd(fiz1,fix2);   
	fiy2 = _mm_hadd_pd(fiy2,fiz2);   
	fix3 = _mm_hadd_pd(fix3,fiy3);   
	fiz3 = _mm_hadd_pd(fiz3,fix4);   
	fiy4 = _mm_hadd_pd(fiy4,fiz4);   
#else
	/* SSE2 */
	/* transpose data */
	GMX_MM_TRANSPOSE2_PD(fix1,fiy1);
	GMX_MM_TRANSPOSE2_PD(fiz1,fix2);
	GMX_MM_TRANSPOSE2_PD(fiy2,fiz2);
	GMX_MM_TRANSPOSE2_PD(fix3,fiy3);
	GMX_MM_TRANSPOSE2_PD(fiz3,fix4);
	GMX_MM_TRANSPOSE2_PD(fiy4,fiz4);
	
	fix1 = _mm_add_pd(fix1,fiy1);
	fiz1 = _mm_add_pd(fiz1,fix2);
	fiy2 = _mm_add_pd(fiy2,fiz2);
	fix3 = _mm_add_pd(fix3,fiy3);
	fiz3 = _mm_add_pd(fiz3,fix4);
	fiy4 = _mm_add_pd(fiy4,fiz4);
#endif
	_mm_storeu_pd( fptr, _mm_add_pd( _mm_loadu_pd(fptr), fix1 ));
	_mm_storeu_pd( fptr+2, _mm_add_pd( _mm_loadu_pd(fptr+2), fiz1 ));
	_mm_storeu_pd( fptr+4, _mm_add_pd( _mm_loadu_pd(fptr+4), fiy2 ));
	_mm_storeu_pd( fptr+6, _mm_add_pd( _mm_loadu_pd(fptr+6), fix3 ));
	_mm_storeu_pd( fptr+8, _mm_add_pd( _mm_loadu_pd(fptr+8), fiz3 ));
	_mm_storeu_pd( fptr+10, _mm_add_pd( _mm_loadu_pd(fptr+10), fiy4 ));
	
	t1 = _mm_shuffle_pd(fiz1,fiy2,_MM_SHUFFLE2(0,1)); 
	fix1 = _mm_add_pd(fix1,t1); 
	t2 = _mm_shuffle_pd(fiz3,fiy4,_MM_SHUFFLE2(0,1)); 
	fix3 = _mm_add_pd(fix3,t2); 
	fix1 = _mm_add_pd(fix1,fix3); /* x and y sums */
	
	
	fiz1 = _mm_add_sd(fiz1, _mm_unpackhi_pd(fiy2,fiy2)); 
	fiz3 = _mm_add_sd(fiz3, _mm_unpackhi_pd(fiy4,fiy4)); 
	fiz1 = _mm_add_sd(fiz1,fiz3); /* z sum */
	
	
	_mm_storeu_pd( fshiftptr, _mm_add_pd( _mm_loadu_pd(fshiftptr), fix1 ));
	_mm_store_sd( fshiftptr+2, _mm_add_sd( _mm_load_sd(fshiftptr+2), fiz1 ));
}


static inline void
gmx_mm_update_1pot_pd(__m128d pot1, double *ptr1)
{
#ifdef GMX_SSE3
	pot1 = _mm_hadd_pd(pot1,pot1);
#else
	/* SSE2 */
	pot1 = _mm_add_pd(pot1, _mm_unpackhi_pd(pot1,pot1));
#endif
	_mm_store_sd(ptr1,_mm_add_sd(pot1,_mm_load_sd(ptr1)));
}
				   
static inline void
gmx_mm_update_2pot_pd(__m128d pot1, double *ptr1, __m128d pot2, double *ptr2)
{
#ifdef GMX_SSE3
	pot1 = _mm_hadd_pd(pot1,pot2); 
#else
	/* SSE2 */
	GMX_MM_TRANSPOSE2_PD(pot1,pot2);
	pot1 = _mm_add_pd(pot1,pot2);
#endif
	pot2 = _mm_unpackhi_pd(pot1,pot1);
	
	_mm_store_sd(ptr1,_mm_add_sd(pot1,_mm_load_sd(ptr1)));
	_mm_store_sd(ptr2,_mm_add_sd(pot2,_mm_load_sd(ptr2)));
}


