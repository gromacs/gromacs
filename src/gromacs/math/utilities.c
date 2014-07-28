/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "utilities.h"

#include "config.h"

#include <assert.h>
#include <limits.h>
#include <math.h>

#ifdef HAVE__FINITE
#include <float.h>
#endif
#ifndef GMX_NATIVE_WINDOWS
#include <fenv.h>
#endif

int gmx_nint(real a)
{
    const real half = .5;
    int        result;

    result = (a < 0.) ? ((int)(a - half)) : ((int)(a + half));
    return result;
}

real cuberoot(real x)
{
    if (x < 0)
    {
        return (-pow(-x, 1.0/3.0));
    }
    else
    {
        return (pow(x, 1.0/3.0));
    }
}

real sign(real x, real y)
{
    if (y < 0)
    {
        return -fabs(x);
    }
    else
    {
        return +fabs(x);
    }
}

/* Double and single precision erf() and erfc() from
 * the Sun Freely Distributable Math Library FDLIBM.
 * See http://www.netlib.org/fdlibm
 * Specific file used: s_erf.c, version 1.3 95/01/18
 */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

static const double
    tiny        = 1e-300,
    half        =  5.00000000000000000000e-01, /* 0x3FE00000, 0x00000000 */
    one         =  1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
    two         =  2.00000000000000000000e+00, /* 0x40000000, 0x00000000 */
/* c = (float)0.84506291151 */
    erx =  8.45062911510467529297e-01,         /* 0x3FEB0AC1, 0x60000000 */
/*
 * Coefficients for approximation to  erf on [0,0.84375]
 */
    efx  =  1.28379167095512586316e-01,        /* 0x3FC06EBA, 0x8214DB69 */
    efx8 =  1.02703333676410069053e+00,        /* 0x3FF06EBA, 0x8214DB69 */
    pp0  =  1.28379167095512558561e-01,        /* 0x3FC06EBA, 0x8214DB68 */
    pp1  = -3.25042107247001499370e-01,        /* 0xBFD4CD7D, 0x691CB913 */
    pp2  = -2.84817495755985104766e-02,        /* 0xBF9D2A51, 0xDBD7194F */
    pp3  = -5.77027029648944159157e-03,        /* 0xBF77A291, 0x236668E4 */
    pp4  = -2.37630166566501626084e-05,        /* 0xBEF8EAD6, 0x120016AC */
    qq1  =  3.97917223959155352819e-01,        /* 0x3FD97779, 0xCDDADC09 */
    qq2  =  6.50222499887672944485e-02,        /* 0x3FB0A54C, 0x5536CEBA */
    qq3  =  5.08130628187576562776e-03,        /* 0x3F74D022, 0xC4D36B0F */
    qq4  =  1.32494738004321644526e-04,        /* 0x3F215DC9, 0x221C1A10 */
    qq5  = -3.96022827877536812320e-06,        /* 0xBED09C43, 0x42A26120 */
/*
 * Coefficients for approximation to  erf  in [0.84375,1.25]
 */
    pa0  = -2.36211856075265944077e-03,        /* 0xBF6359B8, 0xBEF77538 */
    pa1  =  4.14856118683748331666e-01,        /* 0x3FDA8D00, 0xAD92B34D */
    pa2  = -3.72207876035701323847e-01,        /* 0xBFD7D240, 0xFBB8C3F1 */
    pa3  =  3.18346619901161753674e-01,        /* 0x3FD45FCA, 0x805120E4 */
    pa4  = -1.10894694282396677476e-01,        /* 0xBFBC6398, 0x3D3E28EC */
    pa5  =  3.54783043256182359371e-02,        /* 0x3FA22A36, 0x599795EB */
    pa6  = -2.16637559486879084300e-03,        /* 0xBF61BF38, 0x0A96073F */
    qa1  =  1.06420880400844228286e-01,        /* 0x3FBB3E66, 0x18EEE323 */
    qa2  =  5.40397917702171048937e-01,        /* 0x3FE14AF0, 0x92EB6F33 */
    qa3  =  7.18286544141962662868e-02,        /* 0x3FB2635C, 0xD99FE9A7 */
    qa4  =  1.26171219808761642112e-01,        /* 0x3FC02660, 0xE763351F */
    qa5  =  1.36370839120290507362e-02,        /* 0x3F8BEDC2, 0x6B51DD1C */
    qa6  =  1.19844998467991074170e-02,        /* 0x3F888B54, 0x5735151D */
/*
 * Coefficients for approximation to  erfc in [1.25,1/0.35]
 */
    ra0  = -9.86494403484714822705e-03,        /* 0xBF843412, 0x600D6435 */
    ra1  = -6.93858572707181764372e-01,        /* 0xBFE63416, 0xE4BA7360 */
    ra2  = -1.05586262253232909814e+01,        /* 0xC0251E04, 0x41B0E726 */
    ra3  = -6.23753324503260060396e+01,        /* 0xC04F300A, 0xE4CBA38D */
    ra4  = -1.62396669462573470355e+02,        /* 0xC0644CB1, 0x84282266 */
    ra5  = -1.84605092906711035994e+02,        /* 0xC067135C, 0xEBCCABB2 */
    ra6  = -8.12874355063065934246e+01,        /* 0xC0545265, 0x57E4D2F2 */
    ra7  = -9.81432934416914548592e+00,        /* 0xC023A0EF, 0xC69AC25C */
    sa1  =  1.96512716674392571292e+01,        /* 0x4033A6B9, 0xBD707687 */
    sa2  =  1.37657754143519042600e+02,        /* 0x4061350C, 0x526AE721 */
    sa3  =  4.34565877475229228821e+02,        /* 0x407B290D, 0xD58A1A71 */
    sa4  =  6.45387271733267880336e+02,        /* 0x40842B19, 0x21EC2868 */
    sa5  =  4.29008140027567833386e+02,        /* 0x407AD021, 0x57700314 */
    sa6  =  1.08635005541779435134e+02,        /* 0x405B28A3, 0xEE48AE2C */
    sa7  =  6.57024977031928170135e+00,        /* 0x401A47EF, 0x8E484A93 */
    sa8  = -6.04244152148580987438e-02,        /* 0xBFAEEFF2, 0xEE749A62 */
/*
 * Coefficients for approximation to  erfc in [1/.35,28]
 */
    rb0  = -9.86494292470009928597e-03,        /* 0xBF843412, 0x39E86F4A */
    rb1  = -7.99283237680523006574e-01,        /* 0xBFE993BA, 0x70C285DE */
    rb2  = -1.77579549177547519889e+01,        /* 0xC031C209, 0x555F995A */
    rb3  = -1.60636384855821916062e+02,        /* 0xC064145D, 0x43C5ED98 */
    rb4  = -6.37566443368389627722e+02,        /* 0xC083EC88, 0x1375F228 */
    rb5  = -1.02509513161107724954e+03,        /* 0xC0900461, 0x6A2E5992 */
    rb6  = -4.83519191608651397019e+02,        /* 0xC07E384E, 0x9BDC383F */
    sb1  =  3.03380607434824582924e+01,        /* 0x403E568B, 0x261D5190 */
    sb2  =  3.25792512996573918826e+02,        /* 0x40745CAE, 0x221B9F0A */
    sb3  =  1.53672958608443695994e+03,        /* 0x409802EB, 0x189D5118 */
    sb4  =  3.19985821950859553908e+03,        /* 0x40A8FFB7, 0x688C246A */
    sb5  =  2.55305040643316442583e+03,        /* 0x40A3F219, 0xCEDF3BE6 */
    sb6  =  4.74528541206955367215e+02,        /* 0x407DA874, 0xE79FE763 */
    sb7  = -2.24409524465858183362e+01;        /* 0xC03670E2, 0x42712D62 */

double gmx_erfd(double x)
{
#ifdef GMX_FLOAT_FORMAT_IEEE754
    gmx_int32_t hx, ix, i;
    double      R, S, P, Q, s, y, z, r;

    union
    {
        double d;
        int    i[2];
    }
    conv;

    conv.d = x;

#ifdef GMX_IEEE754_BIG_ENDIAN_WORD_ORDER
    hx = conv.i[0];
#else
    hx = conv.i[1];
#endif

    ix = hx&0x7fffffff;
    if (ix >= 0x7ff00000)
    {
        /* erf(nan)=nan */
        i = ((gmx_uint32_t)hx>>31)<<1;
        return (double)(1-i)+one/x; /* erf(+-inf)=+-1 */
    }

    if (ix < 0x3feb0000)
    {
        /* |x|<0.84375 */
        if (ix < 0x3e300000)
        {
            /* |x|<2**-28 */
            if (ix < 0x00800000)
            {
                return 0.125*(8.0*x+efx8*x);  /*avoid underflow */
            }
            return x + efx*x;
        }
        z = x*x;
        r = pp0+z*(pp1+z*(pp2+z*(pp3+z*pp4)));
        s = one+z*(qq1+z*(qq2+z*(qq3+z*(qq4+z*qq5))));
        y = r/s;
        return x + x*y;
    }
    if (ix < 0x3ff40000)
    {
        /* 0.84375 <= |x| < 1.25 */
        s = fabs(x)-one;
        P = pa0+s*(pa1+s*(pa2+s*(pa3+s*(pa4+s*(pa5+s*pa6)))));
        Q = one+s*(qa1+s*(qa2+s*(qa3+s*(qa4+s*(qa5+s*qa6)))));
        if (hx >= 0)
        {
            return erx + P/Q;
        }
        else
        {
            return -erx - P/Q;
        }
    }
    if (ix >= 0x40180000)
    {
        /* inf>|x|>=6 */
        if (hx >= 0)
        {
            return one-tiny;
        }
        else
        {
            return tiny-one;
        }
    }
    x = fabs(x);
    s = one/(x*x);
    if (ix < 0x4006DB6E)
    {
        /* |x| < 1/0.35 */
        R = ra0+s*(ra1+s*(ra2+s*(ra3+s*(ra4+s*(ra5+s*(ra6+s*ra7))))));
        S = one+s*(sa1+s*(sa2+s*(sa3+s*(sa4+s*(sa5+s*(sa6+s*(sa7+s*sa8)))))));
    }
    else
    {
        /* |x| >= 1/0.35 */
        R = rb0+s*(rb1+s*(rb2+s*(rb3+s*(rb4+s*(rb5+s*rb6)))));
        S = one+s*(sb1+s*(sb2+s*(sb3+s*(sb4+s*(sb5+s*(sb6+s*sb7))))));
    }

    conv.d = x;

#ifdef GMX_IEEE754_BIG_ENDIAN_WORD_ORDER
    conv.i[1] = 0;
#else
    conv.i[0] = 0;
#endif

    z = conv.d;

    r  =  exp(-z*z-0.5625)*exp((z-x)*(z+x)+R/S);
    if (hx >= 0)
    {
        return one-r/x;
    }
    else
    {
        return r/x-one;
    }
#else
    /* No IEEE754 information. We need to trust that the OS provides erf(). */
    return erf(x);
#endif
}


double gmx_erfcd(double x)
{
#ifdef GMX_FLOAT_FORMAT_IEEE754
    gmx_int32_t hx, ix;
    double      R, S, P, Q, s, y, z, r;

    union
    {
        double d;
        int    i[2];
    }
    conv;

    conv.d = x;

#ifdef GMX_IEEE754_BIG_ENDIAN_WORD_ORDER
    hx = conv.i[0];
#else
    hx = conv.i[1];
#endif

    ix = hx&0x7fffffff;
    if (ix >= 0x7ff00000)
    {
        /* erfc(nan)=nan */
        /* erfc(+-inf)=0,2 */
        return (double)(((gmx_uint32_t)hx>>31)<<1)+one/x;
    }

    if (ix < 0x3feb0000)
    {
        /* |x|<0.84375 */
        double r1, r2, s1, s2, s3, z2, z4;
        if (ix < 0x3c700000)     /* |x|<2**-56 */
        {
            return one-x;
        }
        z = x*x;
        r = pp0+z*(pp1+z*(pp2+z*(pp3+z*pp4)));
        s = one+z*(qq1+z*(qq2+z*(qq3+z*(qq4+z*qq5))));
        y = r/s;
        if (hx < 0x3fd00000)
        {
            /* x<1/4 */
            return one-(x+x*y);
        }
        else
        {
            r  = x*y;
            r += (x-half);
            return half - r;
        }
    }

    if (ix < 0x3ff40000)
    {
        /* 0.84375 <= |x| < 1.25 */
        s = fabs(x)-one;
        P = pa0+s*(pa1+s*(pa2+s*(pa3+s*(pa4+s*(pa5+s*pa6)))));
        Q = one+s*(qa1+s*(qa2+s*(qa3+s*(qa4+s*(qa5+s*qa6)))));
        if (hx >= 0)
        {
            z  = one-erx; return z - P/Q;
        }
        else
        {
            z = erx+P/Q; return one+z;
        }
    }
    if (ix < 0x403c0000)
    {
        /* |x|<28 */
        x = fabs(x);
        s = one/(x*x);
        if (ix < 0x4006DB6D)
        {
            /* |x| < 1/.35 ~ 2.857143*/
            R = ra0+s*(ra1+s*(ra2+s*(ra3+s*(ra4+s*(ra5+s*(ra6+s*ra7))))));
            S = one+s*(sa1+s*(sa2+s*(sa3+s*(sa4+s*(sa5+s*(sa6+s*(sa7+s*sa8)))))));
        }
        else
        {
            /* |x| >= 1/.35 ~ 2.857143 */
            if (hx < 0 && ix >= 0x40180000)
            {
                return two-tiny; /* x < -6 */
            }
            R = rb0+s*(rb1+s*(rb2+s*(rb3+s*(rb4+s*(rb5+s*rb6)))));
            S = one+s*(sb1+s*(sb2+s*(sb3+s*(sb4+s*(sb5+s*(sb6+s*sb7))))));
        }

        conv.d = x;

#ifdef GMX_IEEE754_BIG_ENDIAN_WORD_ORDER
        conv.i[1] = 0;
#else
        conv.i[0] = 0;
#endif

        z = conv.d;

        r  =  exp(-z*z-0.5625)*exp((z-x)*(z+x)+R/S);

        if (hx > 0)
        {
            return r/x;
        }
        else
        {
            return two-r/x;
        }
    }
    else
    {
        if (hx > 0)
        {
            return tiny*tiny;
        }
        else
        {
            return two-tiny;
        }
    }
#else
    /* No IEEE754 information. We need to trust that the OS provides erfc(). */
    return erfc(x);
#endif
}


static const float
    tinyf =  1e-30,
    halff =  5.0000000000e-01, /* 0x3F000000 */
    onef  =  1.0000000000e+00, /* 0x3F800000 */
    twof  =  2.0000000000e+00, /* 0x40000000 */
/* c = (subfloat)0.84506291151 */
    erxf =  8.4506291151e-01,  /* 0x3f58560b */
/*
 * Coefficients for approximation to  erf on [0,0.84375]
 */
    efxf  =  1.2837916613e-01, /* 0x3e0375d4 */
    efx8f =  1.0270333290e+00, /* 0x3f8375d4 */
    pp0f  =  1.2837916613e-01, /* 0x3e0375d4 */
    pp1f  = -3.2504209876e-01, /* 0xbea66beb */
    pp2f  = -2.8481749818e-02, /* 0xbce9528f */
    pp3f  = -5.7702702470e-03, /* 0xbbbd1489 */
    pp4f  = -2.3763017452e-05, /* 0xb7c756b1 */
    qq1f  =  3.9791721106e-01, /* 0x3ecbbbce */
    qq2f  =  6.5022252500e-02, /* 0x3d852a63 */
    qq3f  =  5.0813062117e-03, /* 0x3ba68116 */
    qq4f  =  1.3249473704e-04, /* 0x390aee49 */
    qq5f  = -3.9602282413e-06, /* 0xb684e21a */
/*
 * Coefficients for approximation to  erf  in [0.84375,1.25]
 */
    pa0f = -2.3621185683e-03,  /* 0xbb1acdc6 */
    pa1f =  4.1485610604e-01,  /* 0x3ed46805 */
    pa2f = -3.7220788002e-01,  /* 0xbebe9208 */
    pa3f =  3.1834661961e-01,  /* 0x3ea2fe54 */
    pa4f = -1.1089469492e-01,  /* 0xbde31cc2 */
    pa5f =  3.5478305072e-02,  /* 0x3d1151b3 */
    pa6f = -2.1663755178e-03,  /* 0xbb0df9c0 */
    qa1f =  1.0642088205e-01,  /* 0x3dd9f331 */
    qa2f =  5.4039794207e-01,  /* 0x3f0a5785 */
    qa3f =  7.1828655899e-02,  /* 0x3d931ae7 */
    qa4f =  1.2617121637e-01,  /* 0x3e013307 */
    qa5f =  1.3637083583e-02,  /* 0x3c5f6e13 */
    qa6f =  1.1984500103e-02,  /* 0x3c445aa3 */
/*
 * Coefficients for approximation to  erfc in [1.25,1/0.35]
 */
    ra0f = -9.8649440333e-03,  /* 0xbc21a093 */
    ra1f = -6.9385856390e-01,  /* 0xbf31a0b7 */
    ra2f = -1.0558626175e+01,  /* 0xc128f022 */
    ra3f = -6.2375331879e+01,  /* 0xc2798057 */
    ra4f = -1.6239666748e+02,  /* 0xc322658c */
    ra5f = -1.8460508728e+02,  /* 0xc3389ae7 */
    ra6f = -8.1287437439e+01,  /* 0xc2a2932b */
    ra7f = -9.8143291473e+00,  /* 0xc11d077e */
    sa1f =  1.9651271820e+01,  /* 0x419d35ce */
    sa2f =  1.3765776062e+02,  /* 0x4309a863 */
    sa3f =  4.3456588745e+02,  /* 0x43d9486f */
    sa4f =  6.4538726807e+02,  /* 0x442158c9 */
    sa5f =  4.2900814819e+02,  /* 0x43d6810b */
    sa6f =  1.0863500214e+02,  /* 0x42d9451f */
    sa7f =  6.5702495575e+00,  /* 0x40d23f7c */
    sa8f = -6.0424413532e-02,  /* 0xbd777f97 */
/*
 * Coefficients for approximation to  erfc in [1/.35,28]
 */
    rb0f = -9.8649431020e-03,  /* 0xbc21a092 */
    rb1f = -7.9928326607e-01,  /* 0xbf4c9dd4 */
    rb2f = -1.7757955551e+01,  /* 0xc18e104b */
    rb3f = -1.6063638306e+02,  /* 0xc320a2ea */
    rb4f = -6.3756646729e+02,  /* 0xc41f6441 */
    rb5f = -1.0250950928e+03,  /* 0xc480230b */
    rb6f = -4.8351919556e+02,  /* 0xc3f1c275 */
    sb1f =  3.0338060379e+01,  /* 0x41f2b459 */
    sb2f =  3.2579251099e+02,  /* 0x43a2e571 */
    sb3f =  1.5367296143e+03,  /* 0x44c01759 */
    sb4f =  3.1998581543e+03,  /* 0x4547fdbb */
    sb5f =  2.5530502930e+03,  /* 0x451f90ce */
    sb6f =  4.7452853394e+02,  /* 0x43ed43a7 */
    sb7f = -2.2440952301e+01;  /* 0xc1b38712 */


typedef union
{
    float         value;
    gmx_uint32_t  word;
} ieee_float_shape_type;

#define GET_FLOAT_WORD(i, d)                 \
    do {                                \
        ieee_float_shape_type gf_u;                   \
        gf_u.value = (d);                     \
        (i)        = gf_u.word;                      \
    } while (0)


#define SET_FLOAT_WORD(d, i)                 \
    do {                                \
        ieee_float_shape_type sf_u;                   \
        sf_u.word = (i);                      \
        (d)       = sf_u.value;                     \
    } while (0)


float gmx_erff(float x)
{
    gmx_int32_t hx, ix, i;
    float       R, S, P, Q, s, y, z, r;

    union
    {
        float  f;
        int    i;
    }
    conv;

    conv.f = x;
    hx     = conv.i;

    ix = hx&0x7fffffff;
    if (ix >= 0x7f800000)
    {
        /* erf(nan)=nan */
        i = ((gmx_uint32_t)hx>>31)<<1;
        return (float)(1-i)+onef/x; /* erf(+-inf)=+-1 */
    }

    if (ix < 0x3f580000)
    {
        /* |x|<0.84375 */
        if (ix < 0x31800000)
        {
            /* |x|<2**-28 */
            if (ix < 0x04000000)
            {
                return (float)0.125*((float)8.0*x+efx8f*x);             /*avoid underflow */
            }
            return x + efxf*x;
        }
        z = x*x;
        r = pp0f+z*(pp1f+z*(pp2f+z*(pp3f+z*pp4f)));
        s = onef+z*(qq1f+z*(qq2f+z*(qq3f+z*(qq4f+z*qq5f))));
        y = r/s;
        return x + x*y;
    }
    if (ix < 0x3fa00000)
    {
        /* 0.84375 <= |x| < 1.25 */
        s = fabs(x)-onef;
        P = pa0f+s*(pa1f+s*(pa2f+s*(pa3f+s*(pa4f+s*(pa5f+s*pa6f)))));
        Q = onef+s*(qa1f+s*(qa2f+s*(qa3f+s*(qa4f+s*(qa5f+s*qa6f)))));
        if (hx >= 0)
        {
            return erxf + P/Q;
        }
        else
        {
            return -erxf - P/Q;
        }
    }
    if (ix >= 0x40c00000)
    {
        /* inf>|x|>=6 */
        if (hx >= 0)
        {
            return onef-tinyf;
        }
        else
        {
            return tinyf-onef;
        }
    }
    x = fabs(x);
    s = onef/(x*x);
    if (ix < 0x4036DB6E)
    {
        /* |x| < 1/0.35 */
        R = ra0f+s*(ra1f+s*(ra2f+s*(ra3f+s*(ra4f+s*(ra5f+s*(ra6f+s*ra7f))))));
        S = onef+s*(sa1f+s*(sa2f+s*(sa3f+s*(sa4f+s*(sa5f+s*(sa6f+s*(sa7f+s*sa8f)))))));
    }
    else
    {
        /* |x| >= 1/0.35 */
        R = rb0f+s*(rb1f+s*(rb2f+s*(rb3f+s*(rb4f+s*(rb5f+s*rb6f)))));
        S = onef+s*(sb1f+s*(sb2f+s*(sb3f+s*(sb4f+s*(sb5f+s*(sb6f+s*sb7f))))));
    }

    conv.f = x;
    conv.i = conv.i & 0xfffff000;
    z      = conv.f;

    r  =  exp(-z*z-(float)0.5625)*exp((z-x)*(z+x)+R/S);
    if (hx >= 0)
    {
        return onef-r/x;
    }
    else
    {
        return r/x-onef;
    }
}

float gmx_erfcf(float x)
{
    gmx_int32_t hx, ix;
    float       R, S, P, Q, s, y, z, r;

    union
    {
        float  f;
        int    i;
    }
    conv;

    conv.f = x;
    hx     = conv.i;

    ix = hx&0x7fffffff;
    if (ix >= 0x7f800000)
    {
        /* erfc(nan)=nan */
        /* erfc(+-inf)=0,2 */
        return (float)(((gmx_uint32_t)hx>>31)<<1)+onef/x;
    }

    if (ix < 0x3f580000)
    {
        /* |x|<0.84375 */
        if (ix < 0x23800000)
        {
            return onef-x;  /* |x|<2**-56 */
        }
        z = x*x;
        r = pp0f+z*(pp1f+z*(pp2f+z*(pp3f+z*pp4f)));
        s = onef+z*(qq1f+z*(qq2f+z*(qq3f+z*(qq4f+z*qq5f))));
        y = r/s;
        if (hx < 0x3e800000)
        {
            /* x<1/4 */
            return onef-(x+x*y);
        }
        else
        {
            r  = x*y;
            r += (x-halff);
            return halff - r;
        }
    }
    if (ix < 0x3fa00000)
    {
        /* 0.84375 <= |x| < 1.25 */
        s = fabs(x)-onef;
        P = pa0f+s*(pa1f+s*(pa2f+s*(pa3f+s*(pa4f+s*(pa5f+s*pa6f)))));
        Q = onef+s*(qa1f+s*(qa2f+s*(qa3f+s*(qa4f+s*(qa5f+s*qa6f)))));
        if (hx >= 0)
        {
            z  = onef-erxf; return z - P/Q;
        }
        else
        {
            z = erxf+P/Q; return onef+z;
        }
    }
    if (ix < 0x41e00000)
    {
        /* |x|<28 */
        x = fabs(x);
        s = onef/(x*x);
        if (ix < 0x4036DB6D)
        {
            /* |x| < 1/.35 ~ 2.857143*/
            R = ra0f+s*(ra1f+s*(ra2f+s*(ra3f+s*(ra4f+s*(ra5f+s*(ra6f+s*ra7f))))));
            S = onef+s*(sa1f+s*(sa2f+s*(sa3f+s*(sa4f+s*(sa5f+s*(sa6f+s*(sa7f+s*sa8f)))))));
        }
        else
        {
            /* |x| >= 1/.35 ~ 2.857143 */
            if (hx < 0 && ix >= 0x40c00000)
            {
                return twof-tinyf;                     /* x < -6 */
            }
            R = rb0f+s*(rb1f+s*(rb2f+s*(rb3f+s*(rb4f+s*(rb5f+s*rb6f)))));
            S = onef+s*(sb1f+s*(sb2f+s*(sb3f+s*(sb4f+s*(sb5f+s*(sb6f+s*sb7f))))));
        }

        conv.f = x;
        conv.i = conv.i & 0xfffff000;
        z      = conv.f;

        r  =  exp(-z*z-(float)0.5625)*exp((z-x)*(z+x)+R/S);
        if (hx > 0)
        {
            return r/x;
        }
        else
        {
            return twof-r/x;
        }
    }
    else
    {
        if (hx > 0)
        {
            return tinyf*tinyf;
        }
        else
        {
            return twof-tinyf;
        }
    }
}


gmx_bool gmx_isfinite(real gmx_unused x)
{
    gmx_bool returnval;

#ifdef HAVE__FINITE
    returnval = _finite(x);
#elif defined HAVE_ISFINITE
    returnval = isfinite(x);
#elif defined HAVE__ISFINITE
    returnval = _isfinite(x);
#else
    /* If no suitable function was found, assume the value is
     * finite. */
    returnval = TRUE;
#endif
    return returnval;
}

gmx_bool gmx_isnan(real x)
{
    return x != x;
}

int
gmx_within_tol(double   f1,
               double   f2,
               double   tol)
{
    /* The or-equal is important - otherwise we return false if f1==f2==0 */
    if (fabs(f1-f2) <= tol*0.5*(fabs(f1)+fabs(f2)) )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int
gmx_numzero(double a)
{
    return gmx_within_tol(a, 0.0, GMX_REAL_MIN/GMX_REAL_EPS);
}

unsigned int
gmx_log2i(unsigned int n)
{
    assert(n != 0); /* behavior differs for 0 */
#if defined(__INTEL_COMPILER)
    return _bit_scan_reverse(n);
#elif defined(__GNUC__) && UINT_MAX == 4294967295U /*also for clang*/
    return __builtin_clz(n) ^ 31U;                 /* xor gets optimized out */
#elif defined(_MSC_VER) && _MSC_VER >= 1400
    {
        unsigned long i;
        _BitScanReverse(&i, n);
        return i;
    }
#elif defined(__xlC__)
    return 31 - __cntlz4(n);
#else
    /* http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup */
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    static const char     LogTable256[256] = {
        -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
        LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
        LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
    };
#undef LT

    unsigned int r;     /* r will be lg(n) */
    unsigned int t, tt; /* temporaries */

    if ((tt = n >> 16) != 0)
    {
        r = ((t = tt >> 8) != 0) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
    }
    else
    {
        r = ((t = n >> 8) != 0) ? 8 + LogTable256[t] : LogTable256[n];
    }
    return r;
#endif
}

gmx_bool
check_int_multiply_for_overflow(gmx_int64_t  a,
                                gmx_int64_t  b,
                                gmx_int64_t *result)
{
    gmx_int64_t sign = 1;
    if ((0 == a) || (0 == b))
    {
        *result = 0;
        return TRUE;
    }
    if (a < 0)
    {
        a    = -a;
        sign = -sign;
    }
    if (b < 0)
    {
        b    = -b;
        sign = -sign;
    }
    if (GMX_INT64_MAX / b < a)
    {
        *result = (sign > 0) ? GMX_INT64_MAX : GMX_INT64_MIN;
        return FALSE;
    }
    *result = sign * a * b;
    return TRUE;
}

int gmx_greatest_common_divisor(int p, int q)
{
    int tmp;
    while (q != 0)
    {
        tmp = q;
        q   = p % q;
        p   = tmp;
    }
    return p;
}

int gmx_feenableexcept()
{
#ifdef HAVE_FEENABLEEXCEPT
    return feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#elif (defined(__i386__) || defined(__x86_64__)) && defined(__APPLE__)
    /* Author:  David N. Williams
     * License:  Public Domain
     *
     * Might also work on non-Apple Unix. But should be tested
     * before enabling.
     */
    unsigned int  excepts = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW;
    static fenv_t fenv;
    unsigned int  new_excepts = excepts & FE_ALL_EXCEPT,
                  old_excepts; // previous masks

    if (fegetenv (&fenv) )
    {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);

    return ( fesetenv (&fenv) ? -1 : old_excepts );
#else
    return -1;
#endif
}
