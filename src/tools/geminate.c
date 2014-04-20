/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBGSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_version.h>
#endif

#include <math.h>
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "geminate.h"
#include "gmx_omp.h"


/* The first few sections of this file contain functions that were adopted,
 * and to some extent modified, by Erik Marklund (erikm[aT]xray.bmc.uu.se,
 * http://folding.bmc.uu.se) from code written by Omer Markovitch (email, url).
 * This is also the case with the function eq10v2().
 *
 * The parts menetioned in the previous paragraph were contributed under the BSD license.
 */


/* This first part is from complex.c which I recieved from Omer Markowitch.
 * - Erik Marklund
 *
 * ------------- from complex.c ------------- */

/* Complexation of a paired number (r,i)                                     */
static gem_complex gem_cmplx(double r, double i)
{
    gem_complex value;
    value.r = r;
    value.i = i;
    return value;
}

/* Complexation of a real number, x */
static gem_complex gem_c(double x)
{
    gem_complex value;
    value.r = x;
    value.i = 0;
    return value;
}

/* Magnitude of a complex number z                                           */
static double gem_cx_abs(gem_complex z)
{
    return (sqrt(z.r*z.r+z.i*z.i));
}

/* Addition of two complex numbers z1 and z2                                 */
static gem_complex gem_cxadd(gem_complex z1, gem_complex z2)
{
    gem_complex value;
    value.r = z1.r+z2.r;
    value.i = z1.i+z2.i;
    return value;
}

/* Addition of a complex number z1 and a real number r */
static gem_complex gem_cxradd(gem_complex z, double r)
{
    gem_complex value;
    value.r = z.r + r;
    value.i = z.i;
    return value;
}

/* Subtraction of two complex numbers z1 and z2                              */
static gem_complex gem_cxsub(gem_complex z1, gem_complex z2)
{
    gem_complex value;
    value.r = z1.r-z2.r;
    value.i = z1.i-z2.i;
    return value;
}

/* Multiplication of two complex numbers z1 and z2                           */
static gem_complex gem_cxmul(gem_complex z1, gem_complex z2)
{
    gem_complex value;
    value.r = z1.r*z2.r-z1.i*z2.i;
    value.i = z1.r*z2.i+z1.i*z2.r;
    return value;
}

/* Square of a complex number z                                              */
static gem_complex gem_cxsq(gem_complex z)
{
    gem_complex value;
    value.r = z.r*z.r-z.i*z.i;
    value.i = z.r*z.i*2.;
    return value;
}

/* multiplication of a complex number z and a real number r */
static gem_complex gem_cxrmul(gem_complex z, double r)
{
    gem_complex value;
    value.r = z.r*r;
    value.i = z.i*r;
    return value;
}

/* Division of two complex numbers z1 and z2                                 */
static gem_complex gem_cxdiv(gem_complex z1, gem_complex z2)
{
    gem_complex value;
    double      num;
    num = z2.r*z2.r+z2.i*z2.i;
    if (num == 0.)
    {
        fprintf(stderr, "ERROR in gem_cxdiv function\n");
    }
    value.r = (z1.r*z2.r+z1.i*z2.i)/num; value.i = (z1.i*z2.r-z1.r*z2.i)/num;
    return value;
}

/* division of a complex z number by a real number x */
static gem_complex gem_cxrdiv(gem_complex z, double r)
{
    gem_complex value;
    value.r = z.r/r;
    value.i = z.i/r;
    return value;
}

/* division of a real number r by a complex number x */
static gem_complex gem_rxcdiv(double r, gem_complex z)
{
    gem_complex value;
    double      f;
    f       = r/(z.r*z.r+z.i*z.i);
    value.r = f*z.r;
    value.i = -f*z.i;
    return value;
}

/* Exponential of a complex number-- exp (z)=|exp(z.r)|*{cos(z.i)+I*sin(z.i)}*/
static gem_complex gem_cxdexp(gem_complex z)
{
    gem_complex value;
    double      exp_z_r;
    exp_z_r = exp(z.r);
    value.r = exp_z_r*cos(z.i);
    value.i = exp_z_r*sin(z.i);
    return value;
}

/* Logarithm of a complex number -- log(z)=log|z|+I*Arg(z),                  */
/*                                  where -PI < Arg(z) < PI                  */
static gem_complex gem_cxlog(gem_complex z)
{
    gem_complex value;
    double      mag2;
    mag2 = z.r*z.r+z.i*z.i;
    if (mag2 < 0.)
    {
        fprintf(stderr, "ERROR in gem_cxlog func\n");
    }
    value.r = log(sqrt(mag2));
    if (z.r == 0.)
    {
        value.i = PI/2.;
        if (z.i < 0.)
        {
            value.i = -value.i;
        }
    }
    else
    {
        value.i = atan2(z.i, z.r);
    }
    return value;
}

/* Square root of a complex number z = |z| exp(I*the) -- z^(1/2)             */
/*                               z^(1/2)=|z|^(1/2)*[cos(the/2)+I*sin(the/2)] */
/*                               where 0 < the < 2*PI                        */
static gem_complex gem_cxdsqrt(gem_complex z)
{
    gem_complex value;
    double      sq;
    sq      = gem_cx_abs(z);
    value.r = sqrt(fabs((sq+z.r)*0.5)); /* z'.r={|z|*[1+cos(the)]/2}^(1/2) */
    value.i = sqrt(fabs((sq-z.r)*0.5)); /* z'.i={|z|*[1-cos(the)]/2}^(1/2) */
    if (z.i < 0.)
    {
        value.r = -value.r;
    }
    return value;
}

/* Complex power of a complex number z1^z2                                   */
static gem_complex gem_cxdpow(gem_complex z1, gem_complex z2)
{
    gem_complex value;
    value = gem_cxdexp(gem_cxmul(gem_cxlog(z1), z2));
    return value;
}

/* ------------ end of complex.c ------------ */

/* This next part was derived from cubic.c, also received from Omer Markovitch.
 * ------------- from cubic.c ------------- */

/* Solver for a cubic equation: x^3-a*x^2+b*x-c=0                            */
static void gem_solve(gem_complex *al, gem_complex *be, gem_complex *gam,
                      double a, double b, double c)
{
    double      t1, t2, two_3, temp;
    gem_complex ctemp, ct3;

    two_3 = pow(2., 1./3.); t1 = -a*a+3.*b; t2 = 2.*a*a*a-9.*a*b+27.*c;
    temp  = 4.*t1*t1*t1+t2*t2;

    ctemp = gem_cmplx(temp, 0.);   ctemp = gem_cxadd(gem_cmplx(t2, 0.), gem_cxdsqrt(ctemp));
    ct3   = gem_cxdpow(ctemp, gem_cmplx(1./3., 0.));

    ctemp = gem_rxcdiv(-two_3*t1/3., ct3);
    ctemp = gem_cxadd(ctemp, gem_cxrdiv(ct3, 3.*two_3));

    *gam = gem_cxadd(gem_cmplx(a/3., 0.), ctemp);

    ctemp = gem_cxmul(gem_cxsq(*gam), gem_cxsq(gem_cxsub(*gam, gem_cmplx(a, 0.))));
    ctemp = gem_cxadd(ctemp, gem_cxmul(gem_cmplx(-4.*c, 0.), *gam));
    ctemp = gem_cxdiv(gem_cxdsqrt(ctemp), *gam);
    *al   = gem_cxrmul(gem_cxsub(gem_cxsub(gem_cmplx(a, 0.), *gam), ctemp), 0.5);
    *be   = gem_cxrmul(gem_cxadd(gem_cxsub(gem_cmplx(a, 0.), *gam), ctemp), 0.5);
}

/* ------------ end of cubic.c ------------ */

/* This next part was derived from cerror.c and rerror.c, also received from Omer Markovitch.
 * ------------- from [cr]error.c ------------- */

/************************************************************/
/* Real valued error function and related functions         */
/* x, y     : real variables                                */
/* erf(x)   : error function                                */
/* erfc(x)  : complementary error function                  */
/* omega(x) : exp(x*x)*erfc(x)                              */
/* W(x,y)   : exp(-x*x)*omega(x+y)=exp(2*x*y+y^2)*erfc(x+y) */
/************************************************************/

/*---------------------------------------------------------------------------*/
/* Utilzed the series approximation of erf(x)                                */
/* Relative error=|err(x)|/erf(x)<EPS                                        */
/* Handbook of Mathematical functions, Abramowitz, p 297                     */
/* Note: When x>=6 series sum deteriorates -> Asymptotic series used instead */
/*---------------------------------------------------------------------------*/
/*  This gives erfc(z) correctly only upto >10-15 */

static double gem_erf(double x)
{
    double n, sum, temp, exp2, xm, x2, x4, x6, x8, x10, x12;
    temp = x;
    sum  = temp;
    xm   = 26.;
    x2   = x*x;
    x4   = x2*x2;
    x6   = x4*x2;
    x8   = x6*x2;
    x10  = x8*x2;
    x12  = x10*x2;
    exp2 = exp(-x2);
    if (x <= xm)
    {
        for (n = 1.; n <= 2000.; n += 1.)
        {
            temp *= 2.*x2/(2.*n+1.);
            sum  += temp;
            if (fabs(temp/sum) < 1.E-16)
            {
                break;
            }
        }

        if (n >= 2000.)
        {
            fprintf(stderr, "In Erf calc - iteration exceeds %lg\n", n);
        }
        sum *= 2./sPI*exp2;
    }
    else
    {
        /* from the asymptotic expansion of experfc(x) */
        sum = (1. - 0.5/x2 + 0.75/x4
               - 1.875/x6 + 6.5625/x8
               - 29.53125/x10 + 162.421875/x12)
            / sPI/x;
        sum *= exp2; /* now sum is erfc(x) */
        sum  = -sum+1.;
    }
    return x >= 0.0 ? sum : -sum;
}

/* Result --> Alex's code for erfc and experfc behaves better than this      */
/* complementray error function.                Returns 1.-erf(x)            */
static double gem_erfc(double x)
{
    double t, z, ans;
    z = fabs(x);
    t = 1.0/(1.0+0.5*z);

    ans = t * exp(-z*z - 1.26551223 +
                  t*(1.00002368 +
                     t*(0.37409196 +
                        t*(0.09678418 +
                           t*(-0.18628806 +
                              t*(0.27886807 +
                                 t*(-1.13520398 +
                                    t*(1.48851587 +
                                       t*(-0.82215223 +
                                          t*0.17087277)))))))));

    return x >= 0.0 ? ans : 2.0-ans;
}

/* omega(x)=exp(x^2)erfc(x)                                                  */
static double gem_omega(double x)
{
    double xm, ans, xx, x4, x6, x8, x10, x12;
    xm  = 26;
    xx  = x*x;
    x4  = xx*xx;
    x6  = x4*xx;
    x8  = x6*xx;
    x10 = x8*xx;
    x12 = x10*xx;

    if (x <= xm)
    {
        ans = exp(xx)*gem_erfc(x);
    }
    else
    {
        /* Asymptotic expansion */
        ans = (1. - 0.5/xx + 0.75/x4 - 1.875/x6 + 6.5625/x8 - 29.53125/x10 + 162.421875/x12) / sPI/x;
    }
    return ans;
}

/*---------------------------------------------------------------------------*/
/* Utilzed the series approximation of erf(z=x+iy)                           */
/* Relative error=|err(z)|/|erf(z)|<EPS                                      */
/* Handbook of Mathematical functions, Abramowitz, p 299                     */
/* comega(z=x+iy)=cexp(z^2)*cerfc(z)                                         */
/*---------------------------------------------------------------------------*/
static gem_complex gem_comega(gem_complex z)
{
    gem_complex value;
    double      x, y;
    double      sumr, sumi, n, n2, f, temp, temp1;
    double      x2, y2, cos_2xy, sin_2xy, cosh_2xy, sinh_2xy, cosh_ny, sinh_ny, exp_y2;

    x        = z.r;
    y        = z.i;
    x2       = x*x;
    y2       = y*y;
    sumr     = 0.;
    sumi     = 0.;
    cos_2xy  = cos(2.*x*y);
    sin_2xy  = sin(2.*x*y);
    cosh_2xy = cosh(2.*x*y);
    sinh_2xy = sinh(2.*x*y);
    exp_y2   = exp(-y2);

    for (n = 1.0, temp = 0.; n <= 2000.; n += 1.0)
    {
        n2      = n*n;
        cosh_ny = cosh(n*y);
        sinh_ny = sinh(n*y);
        f       = exp(-n2/4.)/(n2+4.*x2);
        /* if(f<1.E-200) break; */
        sumr += (2.*x - 2.*x*cosh_ny*cos_2xy + n*sinh_ny*sin_2xy)*f;
        sumi += (2.*x*cosh_ny*sin_2xy + n*sinh_ny*cos_2xy)*f;
        temp1 = sqrt(sumr*sumr+sumi*sumi);
        if (fabs((temp1-temp)/temp1) < 1.E-16)
        {
            break;
        }
        temp = temp1;
    }
    if (n == 2000.)
    {
        fprintf(stderr, "iteration exceeds %lg\n", n);
    }
    sumr *= 2./PI;
    sumi *= 2./PI;

    if (x != 0.)
    {
        f = 1./2./PI/x;
    }
    else
    {
        f = 0.;
    }
    value.r = gem_omega(x)-(f*(1.-cos_2xy)+sumr);
    value.i = -(f*sin_2xy+sumi);
    value   = gem_cxmul(value, gem_cmplx(exp_y2*cos_2xy, exp_y2*sin_2xy));
    return (value);
}

/* ------------ end of [cr]error.c ------------ */

/*_ REVERSIBLE GEMINATE RECOMBINATION
 *
 * Here are the functions for reversible geminate recombination. */

/* Changes the unit from square cm per s to square Ångström per ps,
 * since Omers code uses the latter units while g_mds outputs the former.
 * g_hbond expects a diffusion coefficent given in square cm per s. */
static double sqcm_per_s_to_sqA_per_ps (real D)
{
    fprintf(stdout, "Diffusion coefficient is %f A^2/ps\n", D*1e4);
    return (double)(D*1e4);
}


static double eq10v2(double theoryCt[], double time[], int manytimes,
                     double ka, double kd, t_gemParams *params)
{
    /* Finding the 3 roots */
    int    i;
    double kD, D, r, a, b, c, tsqrt, sumimaginary;
    gem_complex
           alpha, beta, gamma,
        c1, c2, c3, c4,
        oma, omb, omc,
        part1, part2, part3, part4;

    kD = params->kD;
    D  = params->D;
    r  = params->sigma;
    a  = (1.0 + ka/kD) * sqrt(D)/r;
    b  = kd;
    c  = kd * sqrt(D)/r;

    gem_solve(&alpha, &beta, &gamma, a, b, c);
    /* Finding the 3 roots */

    sumimaginary = 0;
    part1        = gem_cxmul(alpha, gem_cxmul(gem_cxadd(beta,  gamma), gem_cxsub(beta,  gamma)));                 /* 1(2+3)(2-3) */
    part2        = gem_cxmul(beta,  gem_cxmul(gem_cxadd(gamma, alpha), gem_cxsub(gamma, alpha)));                 /* 2(3+1)(3-1) */
    part3        = gem_cxmul(gamma, gem_cxmul(gem_cxadd(alpha, beta), gem_cxsub(alpha, beta)));                   /* 3(1+2)(1-2) */
    part4        = gem_cxmul(gem_cxsub(gamma, alpha), gem_cxmul(gem_cxsub(alpha, beta), gem_cxsub(beta, gamma))); /* (3-1)(1-2)(2-3) */

#pragma omp parallel for \
    private(i, tsqrt, oma, omb, omc, c1, c2, c3, c4) \
    reduction(+:sumimaginary) \
    default(shared) \
    schedule(guided)
    for (i = 0; i < manytimes; i++)
    {
        tsqrt         = sqrt(time[i]);
        oma           = gem_comega(gem_cxrmul(alpha, tsqrt));
        c1            = gem_cxmul(oma, gem_cxdiv(part1, part4));
        omb           = gem_comega(gem_cxrmul(beta, tsqrt));
        c2            = gem_cxmul(omb, gem_cxdiv(part2, part4));
        omc           = gem_comega(gem_cxrmul(gamma, tsqrt));
        c3            = gem_cxmul(omc, gem_cxdiv(part3, part4));
        c4.r          = c1.r+c2.r+c3.r;
        c4.i          = c1.i+c2.i+c3.i;
        theoryCt[i]   = c4.r;
        sumimaginary += c4.i * c4.i;
    }

    return sumimaginary;

} /* eq10v2 */

/* This returns the real-valued index(!) to an ACF, equidistant on a log scale. */
static double getLogIndex(const int i, const t_gemParams *params)
{
    return (exp(((double)(i)) * params->logQuota) -1);
}

extern t_gemParams *init_gemParams(const double sigma, const double D,
                                   const real *t, const int len, const int nFitPoints,
                                   const real begFit, const real endFit,
                                   const real ballistic, const int nBalExp, const gmx_bool bDt)
{
    double       tDelta;
    t_gemParams *p;

    snew(p, 1);

    /* A few hardcoded things here. For now anyway. */
/*   p->ka_min   = 0; */
/*   p->ka_max   = 100; */
/*   p->dka      = 10; */
/*   p->kd_min   = 0; */
/*   p->kd_max   = 2; */
/*   p->dkd      = 0.2; */
    p->ka       = 0;
    p->kd       = 0;
/*   p->lsq      = -1; */
/*   p->lifetime = 0; */
    p->sigma    = sigma*10; /* Omer uses Å, not nm */
/*   p->lsq_old  = 99999; */
    p->D        = sqcm_per_s_to_sqA_per_ps(D);
    p->kD       = 4*acos(-1.0)*sigma*p->D;


    /* Parameters used by calcsquare(). Better to calculate them
     * here than in calcsquare every time it's called. */
    p->len = len;
/*   p->logAfterTime = logAfterTime; */
    tDelta       = (t[len-1]-t[0]) / len;
    if (tDelta <= 0)
    {
        gmx_fatal(FARGS, "Time between frames is non-positive!");
    }
    else
    {
        p->tDelta = tDelta;
    }

    p->nExpFit      = nBalExp;
/*   p->nLin         = logAfterTime / tDelta; */
    p->nFitPoints   = nFitPoints;
    p->begFit       = begFit;
    p->endFit       = endFit;
    p->logQuota     = (double)(log(p->len))/(p->nFitPoints-1);
/*   if (p->nLin <= 0) { */
/*     fprintf(stderr, "Number of data points in the linear regime is non-positive!\n"); */
/*     sfree(p); */
/*     return NULL; */
/*   } */
/* We want the same number of data points in the log regime. Not crucial, but seems like a good idea. */
/* p->logDelta = log(((float)len)/p->nFitPoints) / p->nFitPoints;/\* log(((float)len)/p->nLin) / p->nLin; *\/ */
/*   p->logPF    = p->nFitPoints*p->nFitPoints/(float)len; /\* p->nLin*p->nLin/(float)len; *\/ */
/* logPF and logDelta are stitched together with the macro GETLOGINDEX defined in geminate.h */

/*   p->logMult      = pow((float)len, 1.0/nLin);/\* pow(t[len-1]-t[0], 1.0/p->nLin); *\/ */
    p->ballistic    =  ballistic;
    return p;
}

/* There was a misunderstanding regarding the fitting. From our
 * recent correspondence it appears that Omer's code require
 * the ACF data on a log-scale and does not operate on the raw data.
 * This needs to be redone in gemFunc_residual() as well as in the
 * t_gemParams structure. */
#ifdef HAVE_LIBGSL
static double gemFunc_residual2(const gsl_vector *p, void *data)
{
    gemFitData *GD;
    int         i, iLog, nLin, nFitPoints, nData;
    double      r, residual2, *ctTheory, *y;

    GD         = (gemFitData *)data;
    nLin       = GD->params->nLin;
    nFitPoints = GD->params->nFitPoints;
    nData      = GD->nData;
    residual2  = 0;
    ctTheory   = GD->ctTheory;
    y          = GD->y;

    /* Now, we'd save a lot of time by not calculating eq10v2 for all
     * time points, but only those that we sample to calculate the mse.
     */

    eq10v2(GD->ctTheory, GD->doubleLogTime /* GD->time */, nFitPoints /* GD->nData */,
           gsl_vector_get(p, 0), gsl_vector_get(p, 1),
           GD->params);

    fixGemACF(GD->ctTheory, nFitPoints);

    /* Removing a bunch of points from the log-part. */
#pragma omp parallel for schedule(dynamic) \
    firstprivate(nData, ctTheory, y, nFitPoints) \
    private (i, iLog, r) \
    reduction(+:residual2) \
    default(shared)
    for (i = 0; i < nFitPoints; i++)
    {
        iLog       = GD->logtime[i];
        r          = log(ctTheory[i /* iLog */]);
        residual2 += sqr(r-log(y[iLog]));
    }
    residual2 /= nFitPoints; /* Not really necessary for the fitting, but gives more meaning to the returned data. */
    /* printf("ka=%3.5f\tkd=%3.5f\trmse=%3.5f\n", gsl_vector_get(p,0), gsl_vector_get(p,1), residual2); */
    return residual2;

/*   for (i=0; i<nLin; i++) { */
/*     /\* Linear part ----------*\/ */
/*     r = ctTheory[i]; */
/*     residual2 += sqr(r-y[i]); */
/*     /\* Log part -------------*\/ */
/*     iLog = GETLOGINDEX(i, GD->params); */
/*     if (iLog >= nData) */
/*       gmx_fatal(FARGS, "log index out of bounds: %i", iLog); */
/*     r = ctTheory[iLog]; */
/*     residual2 += sqr(r-y[iLog]); */

/*   } */
/*   residual2 /= GD->n; */
/*   return residual2; */
}
#endif

static real* d2r(const double *d, const int nn);

extern real fitGemRecomb(double *ct, double *time, double **ctFit,
                         const int nData, t_gemParams *params)
{

    int         nThreads, i, iter, status, maxiter;
    real        size, d2, tol, *dumpdata;
    size_t      p, n;
    gemFitData *GD;
    char       *dumpstr, dumpname[128];

    /* nmsimplex2 had convergence problems prior to gsl v1.14,
     * but it's O(N) instead of O(N) operations, so let's use it if v >= 1.14 */
#ifdef HAVE_LIBGSL
    gsl_multimin_fminimizer *s;
    gsl_vector              *x, *dx; /* parameters and initial step size */
    gsl_multimin_function    fitFunc;
#ifdef GSL_MAJOR_VERSION
#ifdef GSL_MINOR_VERSION
#if ((GSL_MAJOR_VERSION == 1 && GSL_MINOR_VERSION >= 14) || \
    (GSL_MAJOR_VERSION > 1))
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
#else
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
#endif /* #if ... */
#endif /* GSL_MINOR_VERSION */
#else
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
#endif /* GSL_MAJOR_VERSION */
    fprintf(stdout, "Will fit ka and kd to the ACF according to the reversible geminate recombination model.\n");
#else  /* HAVE_LIBGSL */
    fprintf(stderr, "Sorry, can't do reversible geminate recombination without gsl. "
            "Recompile using --with-gsl.\n");
    return -1;
#endif /* HAVE_LIBGSL */

#ifdef HAVE_LIBGSL
#ifdef GMX_OPENMP
    nThreads = gmx_omp_get_max_threads();
    gmx_omp_set_num_threads(nThreads);
    fprintf(stdout, "We will be using %i threads.\n", nThreads);
#endif

    iter    = 0;
    status  = 0;
    maxiter = 100;
    tol     = 1e-10;

    p = 2;                                        /* Number of parameters to fit. ka and kd.  */
    n = params->nFitPoints; /* params->nLin*2 */; /* Number of points in the reduced dataset  */

    if (params->D <= 0)
    {
        fprintf(stderr, "Fitting of D is not implemented yet. It must be provided on the command line.\n");
        return -1;
    }

/*   if (nData<n) { */
/*     fprintf(stderr, "Reduced data set larger than the complete data set!\n"); */
/*     n=nData; */
/*   } */
    snew(dumpdata, nData);
    snew(GD, 1);

    GD->n        = n;
    GD->y        = ct;
    GD->ctTheory = NULL;
    snew(GD->ctTheory, nData);
    GD->LinLog = NULL;
    snew(GD->LinLog, n);
    GD->time   = time;
    GD->ka     = 0;
    GD->kd     = 0;
    GD->tDelta = time[1]-time[0];
    GD->nData  = nData;
    GD->params = params;
    snew(GD->logtime, params->nFitPoints);
    snew(GD->doubleLogTime, params->nFitPoints);

    for (i = 0; i < params->nFitPoints; i++)
    {
        GD->doubleLogTime[i]  = (double)(getLogIndex(i, params));
        GD->logtime[i]        = (int)(GD->doubleLogTime[i]);
        GD->doubleLogTime[i] *= GD->tDelta;

        if (GD->logtime[i] >= nData)
        {
            fprintf(stderr, "Ayay. It seems we're indexing out of bounds.\n");
            params->nFitPoints = i;
        }
    }

    fitFunc.f      = &gemFunc_residual2;
    fitFunc.n      = 2;
    fitFunc.params = (void*)GD;

    x  = gsl_vector_alloc (fitFunc.n);
    dx = gsl_vector_alloc (fitFunc.n);
    gsl_vector_set (x,  0, 25);
    gsl_vector_set (x,  1, 0.5);
    gsl_vector_set (dx, 0, 0.1);
    gsl_vector_set (dx, 1, 0.01);


    s = gsl_multimin_fminimizer_alloc (T, fitFunc.n);
    gsl_multimin_fminimizer_set (s, &fitFunc, x, dx);
    gsl_vector_free (x);
    gsl_vector_free (dx);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);

        if (status != 0)
        {
            gmx_fatal(FARGS, "Something went wrong in the iteration in minimizer %s:\n \"%s\"\n",
                      gsl_multimin_fminimizer_name(s), gsl_strerror(status));
        }

        d2         = gsl_multimin_fminimizer_minimum(s);
        size       = gsl_multimin_fminimizer_size(s);
        params->ka = gsl_vector_get (s->x, 0);
        params->kd = gsl_vector_get (s->x, 1);

        if (status)
        {
            fprintf(stderr, "%s\n", gsl_strerror(status));
            break;
        }

        status = gsl_multimin_test_size(size, tol);

        if (status == GSL_SUCCESS)
        {
            fprintf(stdout, "Converged to minimum at\n");
        }

        printf ("iter %5d: ka = %2.5f  kd = %2.5f  f() = %7.3f  size = %.3f  chi2 = %2.5f\n",
                iter,
                params->ka,
                params->kd,
                s->fval, size, d2);

        if (iter%1 == 0)
        {
            eq10v2(GD->ctTheory, time, nData, params->ka, params->kd, params);
            /* fixGemACF(GD->ctTheory, nFitPoints); */
            sprintf(dumpname, "Iter_%i.xvg", iter);
            for (i = 0; i < GD->nData; i++)
            {
                dumpdata[i] = (real)(GD->ctTheory[i]);
                if (!gmx_isfinite(dumpdata[i]))
                {
                    gmx_fatal(FARGS, "Non-finite value in acf.");
                }
            }
            dumpN(dumpdata, GD->nData, dumpname);
        }
    }
    while ((status == GSL_CONTINUE) && (iter < maxiter));

    /*   /\* Calculate the theoretical ACF from the parameters one last time. *\/ */
    eq10v2(GD->ctTheory, time, nData, params->ka, params->kd, params);
    *ctFit = GD->ctTheory;

    sfree(GD);
    gsl_multimin_fminimizer_free (s);


    return d2;

#endif /* HAVE_LIBGSL */
}

#ifdef HAVE_LIBGSL
static int balFunc_f(const gsl_vector *x, void *data, gsl_vector *f)
{
    /* C + sum{ A_i * exp(-B_i * t) }*/

    balData *BD;
    int      n, i, j, nexp;
    double  *y, *A, *B, C,     /* There are the parameters to be optimized. */
             t, ct;

    BD   = (balData *)data;
    n    = BD->n;
    nexp = BD->nexp;
    y    = BD->y,
    snew(A, nexp);
    snew(B, nexp);

    for (i = 0; i < nexp; i++)
    {
        A[i] = gsl_vector_get(x, i*2);
        B[i] = gsl_vector_get(x, i*2+1);
    }
    C = gsl_vector_get(x, nexp*2);

    for (i = 0; i < n; i++)
    {
        t  = i*BD->tDelta;
        ct = 0;
        for (j = 0; j < nexp; j++)
        {
            ct += A[j] * exp(B[j] * t);
        }
        ct += C;
        gsl_vector_set (f, i, ct - y[i]);
    }
    return GSL_SUCCESS;
}

/* The derivative stored in jacobian form (J)*/
static int balFunc_df(const gsl_vector *params, void *data, gsl_matrix *J)
{
    balData *BD;
    size_t   n, i, j;
    double  *y, *A, *B, C, /* There are the parameters. */
             t;
    int      nexp;

    BD   = (balData*)data;
    n    = BD->n;
    y    = BD->y;
    nexp = BD->nexp;

    snew(A, nexp);
    snew(B, nexp);

    for (i = 0; i < nexp; i++)
    {
        A[i] = gsl_vector_get(params, i*2);
        B[i] = gsl_vector_get(params, i*2+1);
    }
    C = gsl_vector_get(params, nexp*2);
    for (i = 0; i < n; i++)
    {
        t = i*BD->tDelta;
        for (j = 0; j < nexp; j++)
        {
            gsl_matrix_set (J, i, j*2,   exp(B[j]*t));        /* df(t)/dA_j */
            gsl_matrix_set (J, i, j*2+1, A[j]*t*exp(B[j]*t)); /* df(t)/dB_j */
        }
        gsl_matrix_set (J, i, nexp*2, 1);                     /* df(t)/dC */
    }
    return GSL_SUCCESS;
}

/* Calculation of the function and its derivative */
static int balFunc_fdf(const gsl_vector *params, void *data,
                       gsl_vector *f, gsl_matrix *J)
{
    balFunc_f(params, data, f);
    balFunc_df(params, data, J);
    return GSL_SUCCESS;
}
#endif /* HAVE_LIBGSL */

/* Removes the ballistic term from the beginning of the ACF,
 * just like in Omer's paper.
 */
extern void takeAwayBallistic(double *ct, double *t, int len, real tMax, int nexp, gmx_bool bDerivative)
{

    /* Use nonlinear regression with GSL instead.
     * Fit with 4 exponentials and one constant term,
     * subtract the fatest exponential. */

    int      nData, i, status, iter;
    balData *BD;
    double  *guess,           /* Initial guess. */
    *A,                       /* The fitted parameters. (A1, B1, A2, B2,... C) */
             a[2],
             ddt[2];
    gmx_bool sorted;
    size_t   n;
    size_t   p;

    nData = 0;
    do
    {
        nData++;
    }
    while (t[nData] < tMax+t[0] && nData < len);

    p = nexp*2+1;            /* Number of parameters. */

#ifdef HAVE_LIBGSL
    const gsl_multifit_fdfsolver_type *T
        = gsl_multifit_fdfsolver_lmsder;

    gsl_multifit_fdfsolver   *s;           /* The solver itself. */
    gsl_multifit_function_fdf fitFunction; /* The function to be fitted. */
    gsl_matrix               *covar;       /* Covariance matrix for the parameters.
                                            * We'll not use the result, though. */
    gsl_vector_view           theParams;

    nData = 0;
    do
    {
        nData++;
    }
    while (t[nData] < tMax+t[0] && nData < len);

    guess = NULL;
    n     = nData;

    snew(guess, p);
    snew(A, p);
    covar = gsl_matrix_alloc (p, p);

    /* Set up an initial gess for the parameters.
     * The solver is somewhat sensitive to the initial guess,
     * but this worked fine for a TIP3P box with -geminate dd
     * EDIT: In fact, this seems like a good starting pont for other watermodels too. */
    for (i = 0; i < nexp; i++)
    {
        guess[i*2]   = 0.1;
        guess[i*2+1] = -0.5 + (((double)i)/nexp - 0.5)*0.3;
    }
    guess[nexp * 2] = 0.01;

    theParams = gsl_vector_view_array(guess, p);

    snew(BD, 1);
    BD->n      = n;
    BD->y      = ct;
    BD->tDelta = t[1]-t[0];
    BD->nexp   = nexp;

    fitFunction.f      =  &balFunc_f;
    fitFunction.df     =  &balFunc_df;
    fitFunction.fdf    =  &balFunc_fdf;
    fitFunction.n      =  nData;
    fitFunction.p      =  p;
    fitFunction.params =  BD;

    s = gsl_multifit_fdfsolver_alloc (T, nData, p);
    if (s == NULL)
    {
        gmx_fatal(FARGS, "Could not set up the nonlinear solver.");
    }

    gsl_multifit_fdfsolver_set(s, &fitFunction, &theParams.vector);

    /* \=============================================/ */

    iter = 0;
    do
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);

        if (status)
        {
            break;
        }
        status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
    while (iter < 5000 && status == GSL_CONTINUE);

    if (iter == 5000)
    {
        fprintf(stderr, "The non-linear fitting did not converge in 5000 steps.\n"
                "Check the quality of the fit!\n");
    }
    else
    {
        fprintf(stderr, "Non-linear fitting of ballistic term converged in %d steps.\n\n", (int)iter);
    }
    for (i = 0; i < nexp; i++)
    {
        fprintf(stdout, "%c * exp(%c * t) + ", 'A'+(char)i*2, 'B'+(char)i*2);
    }

    fprintf(stdout, "%c\n", 'A'+(char)nexp*2);
    fprintf(stdout, "Here are the actual numbers for A-%c:\n", 'A'+nexp*2);

    for (i = 0; i < nexp; i++)
    {
        A[i*2]   = gsl_vector_get(s->x, i*2);
        A[i*2+1] = gsl_vector_get(s->x, i*2+1);
        fprintf(stdout, " %g*exp(%g * x) +", A[i*2], A[i*2+1]);
    }
    A[i*2] = gsl_vector_get(s->x, i*2);        /* The last and constant term */
    fprintf(stdout, " %g\n", A[i*2]);

    fflush(stdout);

    /* Implement some check for parameter quality */
    for (i = 0; i < nexp; i++)
    {
        if (A[i*2] < 0 || A[i*2] > 1)
        {
            fprintf(stderr, "WARNING: ----------------------------------\n"
                    " | A coefficient does not lie within [0,1].\n"
                    " | This may or may not be a problem.\n"
                    " | Double check the quality of the fit!\n");
        }
        if (A[i*2+1] > 0)
        {
            fprintf(stderr, "WARNING: ----------------------------------\n"
                    " | One factor in the exponent is positive.\n"
                    " | This could be a problem if the coefficient\n"
                    " | is large. Double check the quality of the fit!\n");
        }
    }
    if (A[i*2] < 0 || A[i*2] > 1)
    {
        fprintf(stderr, "WARNING: ----------------------------------\n"
                " | The constant term does not lie within [0,1].\n"
                " | This may or may not be a problem.\n"
                " | Double check the quality of the fit!\n");
    }

    /* Sort the terms */
    sorted = (nexp > 1) ?  FALSE : TRUE;
    while (!sorted)
    {
        sorted = TRUE;
        for (i = 0; i < nexp-1; i++)
        {
            ddt[0] = A[i*2] * A[i*2+1];
            ddt[1] = A[i*2+2] * A[i*2+3];

            if ((bDerivative && (ddt[0] < 0 && ddt[1] < 0 && ddt[0] > ddt[1])) || /* Compare derivative at t=0... */
                (!bDerivative && (A[i*2+1] > A[i*2+3])))                          /* Or just the coefficient in the exponent */
            {
                sorted = FALSE;
                a[0]   = A[i*2];   /* coefficient */
                a[1]   = A[i*2+1]; /* parameter in the exponent */

                A[i*2]   = A[i*2+2];
                A[i*2+1] = A[i*2+3];

                A[i*2+2] = a[0];
                A[i*2+3] = a[1];
            }
        }
    }

    /* Subtract the fastest component */
    fprintf(stdout, "Fastest component is %g * exp(%g * t)\n"
            "Subtracting fastest component from ACF.\n", A[0], A[1]);

    for (i = 0; i < len; i++)
    {
        ct[i] = (ct[i] - A[0] * exp(A[1] * i*BD->tDelta)) / (1-A[0]);
    }

    sfree(guess);
    sfree(A);

    gsl_multifit_fdfsolver_free(s);
    gsl_matrix_free(covar);
    fflush(stdout);

#else
    /* We have no gsl. */
    fprintf(stderr, "Sorry, can't take away ballistic component without gsl. "
            "Recompile using --with-gsl.\n");
    return;
#endif /* HAVE_LIBGSL */

}

extern void dumpN(const real *e, const int nn, char *fn)
{
    /* For debugging only */
    int   i;
    FILE *f;
    char  standardName[] = "Nt.xvg";
    if (fn == NULL)
    {
        fn = standardName;
    }

    f = fopen(fn, "w");
    fprintf(f,
            "@ type XY\n"
            "@ xaxis label \"Frame\"\n"
            "@ yaxis label \"N\"\n"
            "@ s0 line type 3\n");

    for (i = 0; i < nn; i++)
    {
        fprintf(f, "%-10i %-g\n", i, e[i]);
    }

    fclose(f);
}

static real* d2r(const double *d, const int nn)
{
    real *r;
    int   i;

    snew(r, nn);
    for (i = 0; i < nn; i++)
    {
        r[i] = (real)d[i];
    }

    return r;
}

static void _patchBad(double *ct, int n, double dy)
{
    /* Just do lin. interpolation for now. */
    int i;

    for (i = 1; i < n; i++)
    {
        ct[i] = ct[0]+i*dy;
    }
}

static void patchBadPart(double *ct, int n)
{
    _patchBad(ct, n, (ct[n] - ct[0])/n);
}

static void patchBadTail(double *ct, int n)
{
    _patchBad(ct+1, n-1, ct[1]-ct[0]);

}

extern void fixGemACF(double *ct, int len)
{
    int      i, j, b, e;
    gmx_bool bBad;

    /* Let's separate two things:
     * - identification of bad parts
     * - patching of bad parts.
     */

    b    = 0; /* Start of a bad stretch */
    e    = 0; /* End of a bad stretch */
    bBad = FALSE;

    /* An acf of binary data must be one at t=0. */
    if (fabs(ct[0]-1.0) > 1e-6)
    {
        ct[0] = 1.0;
        fprintf(stderr, "|ct[0]-1.0| = %1.6f. Setting ct[0] to 1.0.\n", fabs(ct[0]-1.0));
    }

    for (i = 0; i < len; i++)
    {

#ifdef HAS_ISFINITE
        if (isfinite(ct[i]))
#elif defined(HAS__ISFINITE)
        if (_isfinite(ct[i]))
#else
        if (1)
#endif
        {
            if (!bBad)
            {
                /* Still on a good stretch. Proceed.*/
                continue;
            }

            /* Patch up preceding bad stretch. */
            if (i == (len-1))
            {
                /* It's the tail */
                if (b <= 1)
                {
                    gmx_fatal(FARGS, "The ACF is mostly NaN or Inf. Aborting.");
                }
                patchBadTail(&(ct[b-2]), (len-b)+1);
            }

            e = i;
            patchBadPart(&(ct[b-1]), (e-b)+1);
            bBad = FALSE;
        }
        else
        {
            if (!bBad)
            {
                b = i;

                bBad = TRUE;
            }
        }
    }
}
