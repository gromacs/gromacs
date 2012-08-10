/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of GROMACS.
 * Copyright (c) 2012-  
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _gmx_math_x86_avx_128_fma_double_h_
#define _gmx_math_x86_avx_128_fma_double_h_

#include <math.h>

#include "gmx_x86_avx_128_fma.h"


#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288
#endif


/************************
 *                      *
 * Simple math routines *
 *                      *
 ************************/

/* 1.0/sqrt(x) */
static gmx_inline __m128d
gmx_mm_invsqrt_pd(__m128d x)
{
    const __m128d half  = _mm_set1_pd(0.5);
    const __m128d three = _mm_set1_pd(3.0);

    /* Lookup instruction only exists in single precision, convert back and forth... */
    __m128d lu = _mm_cvtps_pd(_mm_rsqrt_ps( _mm_cvtpd_ps(x)));

    lu = _mm_mul_pd(_mm_mul_pd(half,lu),_mm_nmacc_pd(_mm_mul_pd(lu,lu),x,three));
    return _mm_mul_pd(_mm_mul_pd(half,lu),_mm_nmacc_pd(_mm_mul_pd(lu,lu),x,three));
}

/* 1.0/sqrt(x), done for a pair of arguments to improve throughput */
static void
gmx_mm_invsqrt_pair_pd(__m128d x1, __m128d x2, __m128d *invsqrt1, __m128d *invsqrt2)
{
    const __m128d half   = _mm_set1_pd(0.5);
    const __m128d three  = _mm_set1_pd(3.0);
    const __m128  halff  = _mm_set1_ps(0.5f);
    const __m128  threef = _mm_set1_ps(3.0f);
    
    __m128 xf,luf;
    __m128d lu1,lu2;
    
    /* Do first N-R step in float for 2x throughput */
    xf  = _mm_shuffle_ps(_mm_cvtpd_ps(x1),_mm_cvtpd_ps(x2),_MM_SHUFFLE(1,0,1,0));
    luf = _mm_rsqrt_ps(xf);
    
    luf = _mm_mul_ps(_mm_mul_ps(halff,luf),_mm_nmacc_ps(_mm_mul_ps(luf,luf),xf,threef));

    
    lu2 = _mm_cvtps_pd(_mm_shuffle_ps(luf,luf,_MM_SHUFFLE(3,2,3,2)));
    lu1 = _mm_cvtps_pd(luf);
    
    *invsqrt1 = _mm_mul_pd(_mm_mul_pd(half,lu1),_mm_nmacc_pd(_mm_mul_pd(lu1,lu1),x1,three));
    *invsqrt2 = _mm_mul_pd(_mm_mul_pd(half,lu2),_mm_nmacc_pd(_mm_mul_pd(lu2,lu2),x2,three));
}

/* sqrt(x) - Do NOT use this (but rather invsqrt) if you actually need 1.0/sqrt(x) */
static gmx_inline __m128d
gmx_mm_sqrt_pd(__m128d x)
{
    __m128d mask;
    __m128d res;

    mask = _mm_cmpeq_pd(x,_mm_setzero_pd());
    res  = _mm_andnot_pd(mask,gmx_mm_invsqrt_pd(x));

    res  = _mm_mul_pd(x,res);

    return res;
}

/* 1.0/x */
static gmx_inline __m128d
gmx_mm_inv_pd(__m128d x)
{
    const __m128d two  = _mm_set1_pd(2.0);

    /* Lookup instruction only exists in single precision, convert back and forth... */
    __m128d lu = _mm_cvtps_pd(_mm_rcp_ps( _mm_cvtpd_ps(x)));

    /* Perform two N-R steps for double precision */
    lu         = _mm_mul_pd(lu,_mm_nmacc_pd(lu,x,two));
    return _mm_mul_pd(lu,_mm_nmacc_pd(lu,x,two));
}

static gmx_inline __m128d
gmx_mm_abs_pd(__m128d x)
{
    const __m128d signmask  = gmx_mm_castsi128_pd( _mm_set_epi32(0x7FFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF) );

    return _mm_and_pd(x,signmask);
}


/*
 * 2^x function.
 *
 * The 2^w term is calculated from a (6,0)-th order (no denominator) Minimax polynomia on the interval
 * [-0.5,0.5].
 *
 * The approximation on [-0.5,0.5] is a rational Padé approximation, 1+2*P(x^2)/(Q(x^2)-P(x^2)),
 * according to the same algorithm as used in the Cephes/netlib math routines.
 */
static __m128d
gmx_mm_exp2_pd(__m128d x)
{
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m128d arglimit = _mm_set1_pd(1022.0);
    const __m128i expbase  = _mm_set1_epi32(1023);

    const __m128d P2       = _mm_set1_pd(2.30933477057345225087e-2);
    const __m128d P1       = _mm_set1_pd(2.02020656693165307700e1);
    const __m128d P0       = _mm_set1_pd(1.51390680115615096133e3);
    /* Q2 == 1.0 */
    const __m128d Q1       = _mm_set1_pd(2.33184211722314911771e2);
    const __m128d Q0       = _mm_set1_pd(4.36821166879210612817e3);
    const __m128d one      = _mm_set1_pd(1.0);
    const __m128d two      = _mm_set1_pd(2.0);

    __m128d valuemask;
    __m128i iexppart;
    __m128d fexppart;
    __m128d intpart;
    __m128d z,z2;
    __m128d PolyP,PolyQ;

    iexppart  = _mm_cvtpd_epi32(x);
    intpart   = _mm_round_pd(x,_MM_FROUND_TO_NEAREST_INT);

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

    valuemask = _mm_cmpge_pd(arglimit,gmx_mm_abs_pd(x));
    fexppart  = _mm_and_pd(valuemask,gmx_mm_castsi128_pd(iexppart));

    z         = _mm_sub_pd(x,intpart);
    z2        = _mm_mul_pd(z,z);

    PolyP     = _mm_macc_pd(P2,z2,P1);
    PolyQ     = _mm_add_pd(z2,Q1);
    PolyP     = _mm_macc_pd(PolyP,z2,P0);
    PolyQ     = _mm_macc_pd(PolyQ,z2,Q0);
    PolyP     = _mm_mul_pd(PolyP,z);

    z         = _mm_mul_pd(PolyP,gmx_mm_inv_pd(_mm_sub_pd(PolyQ,PolyP)));
    z         = _mm_macc_pd(two,z,one);

    z         = _mm_mul_pd(z,fexppart);

    return  z;
}

/* Exponential function. This could be calculated from 2^x as Exp(x)=2^(y), where y=log2(e)*x,
 * but there will then be a small rounding error since we lose some precision due to the
 * multiplication. This will then be magnified a lot by the exponential.
 *
 * Instead, we calculate the fractional part directly as a Padé approximation of
 * Exp(z) on [-0.5,0.5]. We use extended precision arithmetics to calculate the fraction
 * remaining after 2^y, which avoids the precision-loss.
 */
static __m128d
gmx_mm_exp_pd(__m128d exparg)
{
    const __m128d argscale = _mm_set1_pd(1.4426950408889634073599);
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m128d arglimit = _mm_set1_pd(1022.0);
    const __m128i expbase  = _mm_set1_epi32(1023);

    const __m128d invargscale0  = _mm_set1_pd(6.93145751953125e-1);
    const __m128d invargscale1  = _mm_set1_pd(1.42860682030941723212e-6);

    const __m128d P2       = _mm_set1_pd(1.26177193074810590878e-4);
    const __m128d P1       = _mm_set1_pd(3.02994407707441961300e-2);
    /* P0 == 1.0 */
    const __m128d Q3       = _mm_set1_pd(3.00198505138664455042E-6);
    const __m128d Q2       = _mm_set1_pd(2.52448340349684104192E-3);
    const __m128d Q1       = _mm_set1_pd(2.27265548208155028766E-1);
    /* Q0 == 2.0 */
    const __m128d one      = _mm_set1_pd(1.0);
    const __m128d two      = _mm_set1_pd(2.0);

    __m128d valuemask;
    __m128i iexppart;
    __m128d fexppart;
    __m128d intpart;
    __m128d x,z,z2;
    __m128d PolyP,PolyQ;

    x             = _mm_mul_pd(exparg,argscale);

    iexppart  = _mm_cvtpd_epi32(x);
    intpart   = _mm_round_pd(x,_MM_FROUND_TO_NEAREST_INT);

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

    valuemask = _mm_cmpge_pd(arglimit,gmx_mm_abs_pd(x));
    fexppart  = _mm_and_pd(valuemask,gmx_mm_castsi128_pd(iexppart));

    z         = _mm_sub_pd(exparg,_mm_mul_pd(invargscale0,intpart));
    z         = _mm_sub_pd(z,_mm_mul_pd(invargscale1,intpart));

    z2        = _mm_mul_pd(z,z);

    PolyQ     = _mm_macc_pd(Q3,z2,Q2);
    PolyP     = _mm_macc_pd(P2,z2,P1);
    PolyQ     = _mm_macc_pd(PolyQ,z2,Q1);

    PolyP     = _mm_macc_pd(PolyP,z2,one);
    PolyQ     = _mm_macc_pd(PolyQ,z2,two);

    PolyP     = _mm_mul_pd(PolyP,z);

    z         = _mm_mul_pd(PolyP,gmx_mm_inv_pd(_mm_sub_pd(PolyQ,PolyP)));
    z         = _mm_macc_pd(two,z,one);

    z         = _mm_mul_pd(z,fexppart);

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
    const __m128d two        = _mm_set1_pd(2.0);
    const __m128d invsq2     = _mm_set1_pd(1.0/sqrt(2.0));

    const __m128d corr1      = _mm_set1_pd(-2.121944400546905827679e-4);
    const __m128d corr2      = _mm_set1_pd(0.693359375);

    const __m128d P5         = _mm_set1_pd(1.01875663804580931796e-4);
    const __m128d P4         = _mm_set1_pd(4.97494994976747001425e-1);
    const __m128d P3         = _mm_set1_pd(4.70579119878881725854e0);
    const __m128d P2         = _mm_set1_pd(1.44989225341610930846e1);
    const __m128d P1         = _mm_set1_pd(1.79368678507819816313e1);
    const __m128d P0         = _mm_set1_pd(7.70838733755885391666e0);

    const __m128d Q4         = _mm_set1_pd(1.12873587189167450590e1);
    const __m128d Q3         = _mm_set1_pd(4.52279145837532221105e1);
    const __m128d Q2         = _mm_set1_pd(8.29875266912776603211e1);
    const __m128d Q1         = _mm_set1_pd(7.11544750618563894466e1);
    const __m128d Q0         = _mm_set1_pd(2.31251620126765340583e1);

    const __m128d R2         = _mm_set1_pd(-7.89580278884799154124e-1);
    const __m128d R1         = _mm_set1_pd(1.63866645699558079767e1);
    const __m128d R0         = _mm_set1_pd(-6.41409952958715622951e1);

    const __m128d S2         = _mm_set1_pd(-3.56722798256324312549E1);
    const __m128d S1         = _mm_set1_pd(3.12093766372244180303E2);
    const __m128d S0         = _mm_set1_pd(-7.69691943550460008604E2);

    __m128d fexp;
    __m128i iexp;

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
    fexp   = _mm_cvtepi32_pd(iexp);

    x      = _mm_andnot_pd(expmask,x);
    x      = _mm_or_pd(x,one);
    x      = _mm_mul_pd(x,half);

    mask1     = _mm_cmpgt_pd(gmx_mm_abs_pd(fexp),two);
    mask2     = _mm_cmplt_pd(x,invsq2);

    fexp   = _mm_sub_pd(fexp,_mm_and_pd(mask2,one));

    /* If mask1 is set ('A') */
    zA     = _mm_sub_pd(x,half);
    t1     = _mm_blendv_pd( zA,x,mask2 );
    zA     = _mm_sub_pd(t1,half);
    t2     = _mm_blendv_pd( x,zA,mask2 );
    yA     = _mm_mul_pd(half,_mm_add_pd(t2,one));

    xA     = _mm_mul_pd(zA,gmx_mm_inv_pd(yA));
    zA     = _mm_mul_pd(xA,xA);

    /* EVALUATE POLY */
    polyR  = _mm_macc_pd(R2,zA,R1);
    polyR  = _mm_macc_pd(polyR,zA,R0);

    polyS  = _mm_add_pd(zA,S2);
    polyS  = _mm_macc_pd(polyS,zA,S1);
    polyS  = _mm_macc_pd(polyS,zA,S0);

    q      = _mm_mul_pd(polyR,gmx_mm_inv_pd(polyS));
    zA     = _mm_mul_pd(_mm_mul_pd(xA,zA),q);

    zA     = _mm_macc_pd(corr1,fexp,zA);
    zA     = _mm_add_pd(zA,xA);
    zA     = _mm_macc_pd(corr2,fexp,zA);

    /* If mask1 is not set ('B') */
    corr   = _mm_and_pd(mask2,x);
    xB     = _mm_add_pd(x,corr);
    xB     = _mm_sub_pd(xB,one);
    zB     = _mm_mul_pd(xB,xB);

    polyP1 = _mm_macc_pd(P5,zB,P3);
    polyP2 = _mm_macc_pd(P4,zB,P2);
    polyP1 = _mm_macc_pd(polyP1,zB,P1);
    polyP2 = _mm_macc_pd(polyP2,zB,P0);
    polyP1 = _mm_macc_pd(polyP1,xB,polyP2);

    polyQ2 = _mm_macc_pd(Q4,zB,Q2);
    polyQ1 = _mm_add_pd(zB,Q3);
    polyQ1 = _mm_macc_pd(polyQ1,zB,Q1);
    polyQ2 = _mm_macc_pd(polyQ2,zB,Q0);
    polyQ1 = _mm_macc_pd(polyQ1,xB,polyQ2);

    fexp   = _mm_and_pd(fexp,_mm_cmpneq_pd(fexp,_mm_setzero_pd()));

    q      = _mm_mul_pd(polyP1,gmx_mm_inv_pd(polyQ1));
    yB     = _mm_macc_pd(_mm_mul_pd(xB,zB),q,_mm_mul_pd(corr1,fexp));

    yB     = _mm_nmacc_pd(half,zB,yB);
    zB     = _mm_add_pd(xB,yB);
    zB     = _mm_macc_pd(corr2,fexp,zB);

    z      = _mm_blendv_pd( zB,zA,mask1 );

    return z;
}



static __m128d
gmx_mm_erf_pd(__m128d x)
{
    /* Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75] */
    const __m128d CAP4      = _mm_set1_pd(-0.431780540597889301512e-4);
    const __m128d CAP3      = _mm_set1_pd(-0.00578562306260059236059);
    const __m128d CAP2      = _mm_set1_pd(-0.028593586920219752446);
    const __m128d CAP1      = _mm_set1_pd(-0.315924962948621698209);
    const __m128d CAP0      = _mm_set1_pd(0.14952975608477029151);

    const __m128d CAQ5      = _mm_set1_pd(-0.374089300177174709737e-5);
    const __m128d CAQ4      = _mm_set1_pd(0.00015126584532155383535);
    const __m128d CAQ3      = _mm_set1_pd(0.00536692680669480725423);
    const __m128d CAQ2      = _mm_set1_pd(0.0668686825594046122636);
    const __m128d CAQ1      = _mm_set1_pd(0.402604990869284362773);
    /* CAQ0 == 1.0 */
    const __m128d CAoffset  = _mm_set1_pd(0.9788494110107421875);

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5] */
    const __m128d CBP6      = _mm_set1_pd(2.49650423685462752497647637088e-10);
    const __m128d CBP5      = _mm_set1_pd(0.00119770193298159629350136085658);
    const __m128d CBP4      = _mm_set1_pd(0.0164944422378370965881008942733);
    const __m128d CBP3      = _mm_set1_pd(0.0984581468691775932063932439252);
    const __m128d CBP2      = _mm_set1_pd(0.317364595806937763843589437418);
    const __m128d CBP1      = _mm_set1_pd(0.554167062641455850932670067075);
    const __m128d CBP0      = _mm_set1_pd(0.427583576155807163756925301060);
    const __m128d CBQ7      = _mm_set1_pd(0.00212288829699830145976198384930);
    const __m128d CBQ6      = _mm_set1_pd(0.0334810979522685300554606393425);
    const __m128d CBQ5      = _mm_set1_pd(0.2361713785181450957579508850717);
    const __m128d CBQ4      = _mm_set1_pd(0.955364736493055670530981883072);
    const __m128d CBQ3      = _mm_set1_pd(2.36815675631420037315349279199);
    const __m128d CBQ2      = _mm_set1_pd(3.55261649184083035537184223542);
    const __m128d CBQ1      = _mm_set1_pd(2.93501136050160872574376997993);
    /* CBQ0 == 1.0 */

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf] */
    const __m128d CCP6      = _mm_set1_pd(-2.8175401114513378771);
    const __m128d CCP5      = _mm_set1_pd(-3.22729451764143718517);
    const __m128d CCP4      = _mm_set1_pd(-2.5518551727311523996);
    const __m128d CCP3      = _mm_set1_pd(-0.687717681153649930619);
    const __m128d CCP2      = _mm_set1_pd(-0.212652252872804219852);
    const __m128d CCP1      = _mm_set1_pd(0.0175389834052493308818);
    const __m128d CCP0      = _mm_set1_pd(0.00628057170626964891937);

    const __m128d CCQ6      = _mm_set1_pd(5.48409182238641741584);
    const __m128d CCQ5      = _mm_set1_pd(13.5064170191802889145);
    const __m128d CCQ4      = _mm_set1_pd(22.9367376522880577224);
    const __m128d CCQ3      = _mm_set1_pd(15.930646027911794143);
    const __m128d CCQ2      = _mm_set1_pd(11.0567237927800161565);
    const __m128d CCQ1      = _mm_set1_pd(2.79257750980575282228);
    /* CCQ0 == 1.0 */
    const __m128d CCoffset  = _mm_set1_pd(0.5579090118408203125);

    const __m128d one       = _mm_set1_pd(1.0);
    const __m128d two       = _mm_set1_pd(2.0);

    const __m128d signbit   = gmx_mm_castsi128_pd( _mm_set_epi32(0x80000000,0x00000000,0x80000000,0x00000000) );

    __m128d xabs,x2,x4,t,t2,w,w2;
    __m128d PolyAP0,PolyAP1,PolyAQ0,PolyAQ1;
    __m128d PolyBP0,PolyBP1,PolyBQ0,PolyBQ1;
    __m128d PolyCP0,PolyCP1,PolyCQ0,PolyCQ1;
    __m128d res_erf,res_erfcB,res_erfcC,res_erfc,res;
    __m128d mask,expmx2;

    /* Calculate erf() */
    xabs     = gmx_mm_abs_pd(x);
    x2       = _mm_mul_pd(x,x);
    x4       = _mm_mul_pd(x2,x2);

    PolyAP0  = _mm_macc_pd(CAP4,x4,CAP2);
    PolyAP1  = _mm_macc_pd(CAP3,x4,CAP1);
    PolyAP0  = _mm_macc_pd(PolyAP0,x4,CAP0);
    PolyAP0  = _mm_macc_pd(PolyAP1,x2,PolyAP0);

    PolyAQ1  = _mm_macc_pd(CAQ5,x4,CAQ3);
    PolyAQ0  = _mm_macc_pd(CAQ4,x4,CAQ2);
    PolyAQ1  = _mm_macc_pd(PolyAQ1,x4,CAQ1);
    PolyAQ0  = _mm_macc_pd(PolyAQ0,x4,one);
    PolyAQ0  = _mm_macc_pd(PolyAQ1,x2,PolyAQ0);

    res_erf  = _mm_macc_pd(PolyAP0,gmx_mm_inv_pd(PolyAQ0),CAoffset);
    res_erf  = _mm_mul_pd(x,res_erf);

    /* Calculate erfc() in range [1,4.5] */
    t       = _mm_sub_pd(xabs,one);
    t2      = _mm_mul_pd(t,t);

    PolyBP0  = _mm_macc_pd(CBP6,t2,CBP4);
    PolyBP1  = _mm_macc_pd(CBP5,t2,CBP3);
    PolyBP0  = _mm_macc_pd(PolyBP0,t2,CBP2);
    PolyBP1  = _mm_macc_pd(PolyBP1,t2,CBP1);
    PolyBP0  = _mm_macc_pd(PolyBP0,t2,CBP0);
    PolyBP0  = _mm_macc_pd(PolyBP1,t,PolyBP0);

    PolyBQ1 = _mm_macc_pd(CBQ7,t2,CBQ5);
    PolyBQ0 = _mm_macc_pd(CBQ6,t2,CBQ4);
    PolyBQ1 = _mm_macc_pd(PolyBQ1,t2,CBQ3);
    PolyBQ0 = _mm_macc_pd(PolyBQ0,t2,CBQ2);
    PolyBQ1 = _mm_macc_pd(PolyBQ1,t2,CBQ1);
    PolyBQ0 = _mm_macc_pd(PolyBQ0,t2,one);
    PolyBQ0 = _mm_macc_pd(PolyBQ1,t,PolyBQ0);

    res_erfcB = _mm_mul_pd(PolyBP0,gmx_mm_inv_pd(PolyBQ0));

    res_erfcB = _mm_mul_pd(res_erfcB,xabs);

    /* Calculate erfc() in range [4.5,inf] */
    w       = gmx_mm_inv_pd(xabs);
    w2      = _mm_mul_pd(w,w);

    PolyCP0  = _mm_macc_pd(CCP6,w2,CCP4);
    PolyCP1  = _mm_macc_pd(CCP5,w2,CCP3);
    PolyCP0  = _mm_macc_pd(PolyCP0,w2,CCP2);
    PolyCP1  = _mm_macc_pd(PolyCP1,w2,CCP1);
    PolyCP0  = _mm_macc_pd(PolyCP0,w2,CCP0);
    PolyCP0  = _mm_macc_pd(PolyCP1,w,PolyCP0);

    PolyCQ0  = _mm_macc_pd(CCQ6,w2,CCQ4);
    PolyCQ1  = _mm_macc_pd(CCQ5,w2,CCQ3);
    PolyCQ0  = _mm_macc_pd(PolyCQ0,w2,CCQ2);
    PolyCQ1  = _mm_macc_pd(PolyCQ1,w2,CCQ1);
    PolyCQ0  = _mm_macc_pd(PolyCQ0,w2,one);
    PolyCQ0  = _mm_macc_pd(PolyCQ1,w,PolyCQ0);

    expmx2   = gmx_mm_exp_pd( _mm_or_pd(signbit, x2) );

    res_erfcC = _mm_macc_pd(PolyCP0,gmx_mm_inv_pd(PolyCQ0),CCoffset);
    res_erfcC = _mm_mul_pd(res_erfcC,w);

    mask = _mm_cmpgt_pd(xabs,_mm_set1_pd(4.5));
    res_erfc = _mm_blendv_pd(res_erfcB,res_erfcC,mask);

    res_erfc = _mm_mul_pd(res_erfc,expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask = _mm_cmplt_pd(x,_mm_setzero_pd());
    res_erfc = _mm_blendv_pd(res_erfc,_mm_sub_pd(two,res_erfc),mask);

    /* Select erf() or erfc() */
    mask = _mm_cmplt_pd(xabs,one);
    res  = _mm_blendv_pd(_mm_sub_pd(one,res_erfc),res_erf,mask);

    return res;
}


static __m128d
gmx_mm_erfc_pd(__m128d x)
{
    /* Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75] */
    const __m128d CAP4      = _mm_set1_pd(-0.431780540597889301512e-4);
    const __m128d CAP3      = _mm_set1_pd(-0.00578562306260059236059);
    const __m128d CAP2      = _mm_set1_pd(-0.028593586920219752446);
    const __m128d CAP1      = _mm_set1_pd(-0.315924962948621698209);
    const __m128d CAP0      = _mm_set1_pd(0.14952975608477029151);

    const __m128d CAQ5      = _mm_set1_pd(-0.374089300177174709737e-5);
    const __m128d CAQ4      = _mm_set1_pd(0.00015126584532155383535);
    const __m128d CAQ3      = _mm_set1_pd(0.00536692680669480725423);
    const __m128d CAQ2      = _mm_set1_pd(0.0668686825594046122636);
    const __m128d CAQ1      = _mm_set1_pd(0.402604990869284362773);
    /* CAQ0 == 1.0 */
    const __m128d CAoffset  = _mm_set1_pd(0.9788494110107421875);

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5] */
    const __m128d CBP6      = _mm_set1_pd(2.49650423685462752497647637088e-10);
    const __m128d CBP5      = _mm_set1_pd(0.00119770193298159629350136085658);
    const __m128d CBP4      = _mm_set1_pd(0.0164944422378370965881008942733);
    const __m128d CBP3      = _mm_set1_pd(0.0984581468691775932063932439252);
    const __m128d CBP2      = _mm_set1_pd(0.317364595806937763843589437418);
    const __m128d CBP1      = _mm_set1_pd(0.554167062641455850932670067075);
    const __m128d CBP0      = _mm_set1_pd(0.427583576155807163756925301060);
    const __m128d CBQ7      = _mm_set1_pd(0.00212288829699830145976198384930);
    const __m128d CBQ6      = _mm_set1_pd(0.0334810979522685300554606393425);
    const __m128d CBQ5      = _mm_set1_pd(0.2361713785181450957579508850717);
    const __m128d CBQ4      = _mm_set1_pd(0.955364736493055670530981883072);
    const __m128d CBQ3      = _mm_set1_pd(2.36815675631420037315349279199);
    const __m128d CBQ2      = _mm_set1_pd(3.55261649184083035537184223542);
    const __m128d CBQ1      = _mm_set1_pd(2.93501136050160872574376997993);
    /* CBQ0 == 1.0 */

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf] */
    const __m128d CCP6      = _mm_set1_pd(-2.8175401114513378771);
    const __m128d CCP5      = _mm_set1_pd(-3.22729451764143718517);
    const __m128d CCP4      = _mm_set1_pd(-2.5518551727311523996);
    const __m128d CCP3      = _mm_set1_pd(-0.687717681153649930619);
    const __m128d CCP2      = _mm_set1_pd(-0.212652252872804219852);
    const __m128d CCP1      = _mm_set1_pd(0.0175389834052493308818);
    const __m128d CCP0      = _mm_set1_pd(0.00628057170626964891937);

    const __m128d CCQ6      = _mm_set1_pd(5.48409182238641741584);
    const __m128d CCQ5      = _mm_set1_pd(13.5064170191802889145);
    const __m128d CCQ4      = _mm_set1_pd(22.9367376522880577224);
    const __m128d CCQ3      = _mm_set1_pd(15.930646027911794143);
    const __m128d CCQ2      = _mm_set1_pd(11.0567237927800161565);
    const __m128d CCQ1      = _mm_set1_pd(2.79257750980575282228);
    /* CCQ0 == 1.0 */
    const __m128d CCoffset  = _mm_set1_pd(0.5579090118408203125);

    const __m128d one       = _mm_set1_pd(1.0);
    const __m128d two       = _mm_set1_pd(2.0);

    const __m128d signbit   = gmx_mm_castsi128_pd( _mm_set_epi32(0x80000000,0x00000000,0x80000000,0x00000000) );

    __m128d xabs,x2,x4,t,t2,w,w2;
    __m128d PolyAP0,PolyAP1,PolyAQ0,PolyAQ1;
    __m128d PolyBP0,PolyBP1,PolyBQ0,PolyBQ1;
    __m128d PolyCP0,PolyCP1,PolyCQ0,PolyCQ1;
    __m128d res_erf,res_erfcB,res_erfcC,res_erfc,res;
    __m128d mask,expmx2;

    /* Calculate erf() */
    xabs     = gmx_mm_abs_pd(x);
    x2       = _mm_mul_pd(x,x);
    x4       = _mm_mul_pd(x2,x2);

    PolyAP0  = _mm_macc_pd(CAP4,x4,CAP2);
    PolyAP1  = _mm_macc_pd(CAP3,x4,CAP1);
    PolyAP0  = _mm_macc_pd(PolyAP0,x4,CAP0);
    PolyAP0  = _mm_macc_pd(PolyAP1,x2,PolyAP0);

    PolyAQ1  = _mm_macc_pd(CAQ5,x4,CAQ3);
    PolyAQ0  = _mm_macc_pd(CAQ4,x4,CAQ2);
    PolyAQ1  = _mm_macc_pd(PolyAQ1,x4,CAQ1);
    PolyAQ0  = _mm_macc_pd(PolyAQ0,x4,one);
    PolyAQ0  = _mm_macc_pd(PolyAQ1,x2,PolyAQ0);

    res_erf  = _mm_macc_pd(PolyAP0,gmx_mm_inv_pd(PolyAQ0),CAoffset);
    res_erf  = _mm_mul_pd(x,res_erf);

    /* Calculate erfc() in range [1,4.5] */
    t       = _mm_sub_pd(xabs,one);
    t2      = _mm_mul_pd(t,t);

    PolyBP0  = _mm_macc_pd(CBP6,t2,CBP4);
    PolyBP1  = _mm_macc_pd(CBP5,t2,CBP3);
    PolyBP0  = _mm_macc_pd(PolyBP0,t2,CBP2);
    PolyBP1  = _mm_macc_pd(PolyBP1,t2,CBP1);
    PolyBP0  = _mm_macc_pd(PolyBP0,t2,CBP0);
    PolyBP0  = _mm_macc_pd(PolyBP1,t,PolyBP0);

    PolyBQ1 = _mm_macc_pd(CBQ7,t2,CBQ5);
    PolyBQ0 = _mm_macc_pd(CBQ6,t2,CBQ4);
    PolyBQ1 = _mm_macc_pd(PolyBQ1,t2,CBQ3);
    PolyBQ0 = _mm_macc_pd(PolyBQ0,t2,CBQ2);
    PolyBQ1 = _mm_macc_pd(PolyBQ1,t2,CBQ1);
    PolyBQ0 = _mm_macc_pd(PolyBQ0,t2,one);
    PolyBQ0 = _mm_macc_pd(PolyBQ1,t,PolyBQ0);

    res_erfcB = _mm_mul_pd(PolyBP0,gmx_mm_inv_pd(PolyBQ0));

    res_erfcB = _mm_mul_pd(res_erfcB,xabs);

    /* Calculate erfc() in range [4.5,inf] */
    w       = gmx_mm_inv_pd(xabs);
    w2      = _mm_mul_pd(w,w);

    PolyCP0  = _mm_macc_pd(CCP6,w2,CCP4);
    PolyCP1  = _mm_macc_pd(CCP5,w2,CCP3);
    PolyCP0  = _mm_macc_pd(PolyCP0,w2,CCP2);
    PolyCP1  = _mm_macc_pd(PolyCP1,w2,CCP1);
    PolyCP0  = _mm_macc_pd(PolyCP0,w2,CCP0);
    PolyCP0  = _mm_macc_pd(PolyCP1,w,PolyCP0);

    PolyCQ0  = _mm_macc_pd(CCQ6,w2,CCQ4);
    PolyCQ1  = _mm_macc_pd(CCQ5,w2,CCQ3);
    PolyCQ0  = _mm_macc_pd(PolyCQ0,w2,CCQ2);
    PolyCQ1  = _mm_macc_pd(PolyCQ1,w2,CCQ1);
    PolyCQ0  = _mm_macc_pd(PolyCQ0,w2,one);
    PolyCQ0  = _mm_macc_pd(PolyCQ1,w,PolyCQ0);

    expmx2   = gmx_mm_exp_pd( _mm_or_pd(signbit, x2) );

    res_erfcC = _mm_macc_pd(PolyCP0,gmx_mm_inv_pd(PolyCQ0),CCoffset);
    res_erfcC = _mm_mul_pd(res_erfcC,w);

    mask = _mm_cmpgt_pd(xabs,_mm_set1_pd(4.5));
    res_erfc = _mm_blendv_pd(res_erfcB,res_erfcC,mask);

    res_erfc = _mm_mul_pd(res_erfc,expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask = _mm_cmplt_pd(x,_mm_setzero_pd());
    res_erfc = _mm_blendv_pd(res_erfc,_mm_sub_pd(two,res_erfc),mask);

    /* Select erf() or erfc() */
    mask = _mm_cmplt_pd(xabs,one);
    res  = _mm_blendv_pd(res_erfc,_mm_sub_pd(one,res_erf),mask);

    return res;
}




static int
gmx_mm_sincos_pd(__m128d x,
                 __m128d *sinval,
                 __m128d *cosval)
{
#ifdef _MSC_VER
    __declspec(align(16))
    const double sintable[34] =
    {
        1.00000000000000000e+00 , 0.00000000000000000e+00 ,
        9.95184726672196929e-01 , 9.80171403295606036e-02 ,
        9.80785280403230431e-01 , 1.95090322016128248e-01 ,
        9.56940335732208824e-01 , 2.90284677254462331e-01 ,
        9.23879532511286738e-01 , 3.82683432365089782e-01 ,
        8.81921264348355050e-01 , 4.71396736825997642e-01 ,
        8.31469612302545236e-01 , 5.55570233019602178e-01 ,
        7.73010453362736993e-01 , 6.34393284163645488e-01 ,
        7.07106781186547573e-01 , 7.07106781186547462e-01 ,
        6.34393284163645599e-01 , 7.73010453362736882e-01 ,
        5.55570233019602289e-01 , 8.31469612302545125e-01 ,
        4.71396736825997809e-01 , 8.81921264348354939e-01 ,
        3.82683432365089837e-01 , 9.23879532511286738e-01 ,
        2.90284677254462276e-01 , 9.56940335732208935e-01 ,
        1.95090322016128304e-01 , 9.80785280403230431e-01 ,
        9.80171403295607702e-02 , 9.95184726672196818e-01 ,
        0.0 , 1.00000000000000000e+00
    };
#else
    const __m128d sintable[17] =
    {
        _mm_set_pd( 0.0 , 1.0 ),
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
        _mm_set_pd(  1.0 , 0.0 )
    };
#endif

    const __m128d signmask    = gmx_mm_castsi128_pd( _mm_set_epi32(0x7FFFFFFF,0xFFFFFFFF,0x7FFFFFFF,0xFFFFFFFF) );
    const __m128i signbit_epi32  = _mm_set1_epi32(0x80000000);

    const __m128d tabscale      = _mm_set1_pd(32.0/M_PI);
    const __m128d invtabscale0  = _mm_set1_pd(9.81747508049011230469e-02);
    const __m128d invtabscale1  = _mm_set1_pd(1.96197799156550576057e-08);
    const __m128i ione          = _mm_set1_epi32(1);
    const __m128i i32           = _mm_set1_epi32(32);
    const __m128i i16           = _mm_set1_epi32(16);
    const __m128i tabmask       = _mm_set1_epi32(0x3F);
    const __m128d sinP7         = _mm_set1_pd(-1.0/5040.0);
    const __m128d sinP5         = _mm_set1_pd(1.0/120.0);
    const __m128d sinP3         = _mm_set1_pd(-1.0/6.0);
    const __m128d sinP1         = _mm_set1_pd(1.0);

    const __m128d cosP6         = _mm_set1_pd(-1.0/720.0);
    const __m128d cosP4         = _mm_set1_pd(1.0/24.0);
    const __m128d cosP2         = _mm_set1_pd(-1.0/2.0);
    const __m128d cosP0         = _mm_set1_pd(1.0);

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

    xpoint   = _mm_round_pd(scalex,_MM_FROUND_TO_NEAREST_INT);

    /* Extended precision arithmetics */
    z        = _mm_nmacc_pd(invtabscale0,xpoint,xabs);
    z        = _mm_nmacc_pd(invtabscale1,xpoint,z);

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
    tabidx    = _mm_blendv_epi8(tabidx,corridx,imask);
    /* tabidx is now in range [0..16] */
    ssign     = _mm_cvtepi32_pd( _mm_or_si128( sswapsign , ione ) );
    csign     = _mm_cvtepi32_pd( _mm_or_si128( cswapsign , ione ) );

#ifdef _MSC_VER
    ypoint0  = _mm_load_pd(sintable + 2*_mm_extract_epi32(tabidx,0));
    ypoint1  = _mm_load_pd(sintable + 2*_mm_extract_epi32(tabidx,1));
#else
    ypoint0  = sintable[_mm_extract_epi32(tabidx,0)];
    ypoint1  = sintable[_mm_extract_epi32(tabidx,1)];
#endif
    sinpoint = _mm_unpackhi_pd(ypoint0,ypoint1);
    cospoint = _mm_unpacklo_pd(ypoint0,ypoint1);

    sinpoint = _mm_mul_pd(sinpoint,ssign);
    cospoint = _mm_mul_pd(cospoint,csign);

    z2       = _mm_mul_pd(z,z);

    polySin  = _mm_macc_pd(sinP7,z2,sinP5);
    polySin  = _mm_macc_pd(polySin,z2,sinP3);
    polySin  = _mm_macc_pd(polySin,z2,sinP1);
    polySin  = _mm_mul_pd(polySin,z);

    polyCos  = _mm_macc_pd(cosP6,z2,cosP4);
    polyCos  = _mm_macc_pd(polyCos,z2,cosP2);
    polyCos  = _mm_macc_pd(polyCos,z2,cosP0);

    *sinval  = _mm_xor_pd(_mm_add_pd( _mm_mul_pd(sinpoint,polyCos) , _mm_mul_pd(cospoint,polySin) ),xsign);
    *cosval  = _mm_sub_pd( _mm_mul_pd(cospoint,polyCos) , _mm_mul_pd(sinpoint,polySin) );

    return 0;
}

/*
 * IMPORTANT: Do NOT call both sin & cos if you need both results, since each of them
 * will then call the sincos() routine and waste a factor 2 in performance!
 */
static __m128d
gmx_mm_sin_pd(__m128d x)
{
    __m128d s,c;
    gmx_mm_sincos_pd(x,&s,&c);
    return s;
}

/*
 * IMPORTANT: Do NOT call both sin & cos if you need both results, since each of them
 * will then call the sincos() routine and waste a factor 2 in performance!
 */
static __m128d
gmx_mm_cos_pd(__m128d x)
{
    __m128d s,c;
    gmx_mm_sincos_pd(x,&s,&c);
    return c;
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
    const __m128d morebits  = _mm_set1_pd(6.123233995736765886130e-17);

    const __m128d P5        = _mm_set1_pd(4.253011369004428248960e-3);
    const __m128d P4        = _mm_set1_pd(-6.019598008014123785661e-1);
    const __m128d P3        = _mm_set1_pd(5.444622390564711410273e0);
    const __m128d P2        = _mm_set1_pd(-1.626247967210700244449e1);
    const __m128d P1        = _mm_set1_pd(1.956261983317594739197e1);
    const __m128d P0        = _mm_set1_pd(-8.198089802484824371615e0);

    const __m128d Q4        = _mm_set1_pd(-1.474091372988853791896e1);
    const __m128d Q3        = _mm_set1_pd(7.049610280856842141659e1);
    const __m128d Q2        = _mm_set1_pd(-1.471791292232726029859e2);
    const __m128d Q1        = _mm_set1_pd(1.395105614657485689735e2);
    const __m128d Q0        = _mm_set1_pd(-4.918853881490881290097e1);

    const __m128d R4        = _mm_set1_pd(2.967721961301243206100e-3);
    const __m128d R3        = _mm_set1_pd(-5.634242780008963776856e-1);
    const __m128d R2        = _mm_set1_pd(6.968710824104713396794e0);
    const __m128d R1        = _mm_set1_pd(-2.556901049652824852289e1);
    const __m128d R0        = _mm_set1_pd(2.853665548261061424989e1);

    const __m128d S3        = _mm_set1_pd(-2.194779531642920639778e1);
    const __m128d S2        = _mm_set1_pd(1.470656354026814941758e2);
    const __m128d S1        = _mm_set1_pd(-3.838770957603691357202e2);
    const __m128d S0        = _mm_set1_pd(3.424398657913078477438e2);

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
    RA    = _mm_macc_pd(R4,zz2,R2);
    RB    = _mm_macc_pd(R3,zz2,R1);
    RA    = _mm_macc_pd(RA,zz2,R0);
    RA    = _mm_macc_pd(RB,zz,RA);

    /* S, SA = zz2 */
    SB    = _mm_macc_pd(S3,zz2,S1);
    SA    = _mm_add_pd(zz2,S2);
    SA    = _mm_macc_pd(SA,zz2,S0);
    SA    = _mm_macc_pd(SB,zz,SA);

    /* P */
    PA    = _mm_macc_pd(P5,ww2,P3);
    PB    = _mm_macc_pd(P4,ww2,P2);
    PA    = _mm_macc_pd(PA,ww2,P1);
    PB    = _mm_macc_pd(PB,ww2,P0);
    PA    = _mm_macc_pd(PA,ww,PB);

    /* Q, QA = ww2 */
    QB    = _mm_macc_pd(Q4,ww2,Q2);
    QA    = _mm_add_pd(ww2,Q3);
    QA    = _mm_macc_pd(QA,ww2,Q1);
    QB    = _mm_macc_pd(QB,ww2,Q0);
    QA    = _mm_macc_pd(QA,ww,QB);

    RA    = _mm_mul_pd(RA,zz);
    PA    = _mm_mul_pd(PA,ww);

    nom   = _mm_blendv_pd( PA,RA,mask );
    denom = _mm_blendv_pd( QA,SA,mask );

    q     = _mm_mul_pd( nom, gmx_mm_inv_pd(denom) );

    zz    = _mm_add_pd(zz,zz);
    zz    = gmx_mm_sqrt_pd(zz);
    z     = _mm_sub_pd(quarterpi,zz);
    zz    = _mm_mul_pd(zz,q);
    zz    = _mm_sub_pd(zz,morebits);
    z     = _mm_sub_pd(z,zz);
    z     = _mm_add_pd(z,quarterpi);

    w     = _mm_macc_pd(xabs,q,xabs);

    z     = _mm_blendv_pd( w,z,mask );

    mask  = _mm_cmpgt_pd(xabs,limit2);
    z     = _mm_blendv_pd( xabs,z,mask );

    z = _mm_xor_pd(z,sign);

    return z;
}


static __m128d
gmx_mm_acos_pd(__m128d x)
{
    const __m128d signmask  = gmx_mm_castsi128_pd( _mm_set_epi32(0x7FFFFFFF,0xFFFFFFFF,0x7FFFFFFF,0xFFFFFFFF) );
    const __m128d one        = _mm_set1_pd(1.0);
    const __m128d half       = _mm_set1_pd(0.5);
    const __m128d pi         = _mm_set1_pd(M_PI);
    const __m128d quarterpi0 = _mm_set1_pd(7.85398163397448309616e-1);
    const __m128d quarterpi1 = _mm_set1_pd(6.123233995736765886130e-17);


    __m128d mask1;

    __m128d z,z1,z2;

    mask1 = _mm_cmpgt_pd(x,half);
    z1    = _mm_mul_pd(half,_mm_sub_pd(one,x));
    z1    = gmx_mm_sqrt_pd(z1);
    z     = _mm_blendv_pd( x,z1,mask1 );

    z     = gmx_mm_asin_pd(z);

    z1    = _mm_add_pd(z,z);

    z2    = _mm_sub_pd(quarterpi0,z);
    z2    = _mm_add_pd(z2,quarterpi1);
    z2    = _mm_add_pd(z2,quarterpi0);

    z     = _mm_blendv_pd(z2,z1,mask1);

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

    P_A    = _mm_macc_pd(P4,z2,P2);
    P_B    = _mm_macc_pd(P3,z2,P1);
    P_A    = _mm_macc_pd(P_A,z2,P0);
    P_A    = _mm_macc_pd(P_B,z,P_A);

    /* Q_A = z2 */
    Q_B    = _mm_macc_pd(Q4,z2,Q2);
    Q_A    = _mm_add_pd(z2,Q3);
    Q_A    = _mm_macc_pd(Q_A,z2,Q1);
    Q_B    = _mm_macc_pd(Q_B,z2,Q0);
    Q_A    = _mm_macc_pd(Q_A,z,Q_B);

    z      = _mm_mul_pd(z,P_A);
    z      = _mm_mul_pd(z,gmx_mm_inv_pd(Q_A));
    z      = _mm_macc_pd(z,x,x);

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

    w         = _mm_blendv_pd(pi,minuspi,masky_lt);
    w         = _mm_and_pd(w,maskx_lt);

    w         = _mm_andnot_pd(maskall,w);

    z         = _mm_add_pd(z,w);

    return z;
}

#endif /*_gmx_math_x86_avx_128_fma_double_h_ */
