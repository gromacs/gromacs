/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _maths_h
#define _maths_h

#include <math.h>
#include "types/simple.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef M_PI
#define	M_PI		3.14159265358979323846
#endif

#ifndef M_PI_2
#define	M_PI_2		1.57079632679489661923
#endif

#ifndef M_2PI
#define	M_2PI		6.28318530717958647692
#endif
    
#ifndef M_SQRT2
#define M_SQRT2 sqrt(2.0)
#endif
    
/* Suzuki-Yoshida Constants, for n=3 and n=5, for symplectic integration  */
/* for n=1, w0 = 1 */
/* for n=3, w0 = w2 = 1/(2-2^-(1/3)), w1 = 1-2*w0 */
/* for n=5, w0 = w1 = w3 = w4 = 1/(4-4^-(1/3)), w1 = 1-4*w0 */
    
#define MAX_SUZUKI_YOSHIDA_NUM 5
#define SUZUKI_YOSHIDA_NUM  5

static const double sy_const_1[] = { 1. };
static const double sy_const_3[] = { 0.828981543588751,-0.657963087177502,0.828981543588751 };
static const double sy_const_5[] = { 0.2967324292201065,0.2967324292201065,-0.186929716880426,0.2967324292201065,0.2967324292201065 };

static const double* sy_const[] = {
    NULL,
    sy_const_1,
    NULL,
    sy_const_3,
    NULL,
    sy_const_5
};

/*
static const double sy_const[MAX_SUZUKI_YOSHIDA_NUM+1][MAX_SUZUKI_YOSHIDA_NUM+1] = {
    {},
    {1},
    {},
    {0.828981543588751,-0.657963087177502,0.828981543588751},
    {},
    {0.2967324292201065,0.2967324292201065,-0.186929716880426,0.2967324292201065,0.2967324292201065}
};*/

int		gmx_nint(real a);
real    sign(real x,real y);

int		gmx_nint(real a);
real    sign(real x,real y);
real    cuberoot (real a);
real    gmx_erf(real x);
real    gmx_erfc(real x);

/*! \brief Check if two numbers are within a tolerance
 *
 *  This routine checks if the relative difference between two numbers is
 *  approximately within the given tolerance, defined as
 *  fabs(f1-f2)<=tolerance*fabs(f1+f2).
 *
 *  To check if two floating-point numbers are almost identical, use this routine 
 *  with the tolerance GMX_REAL_EPS, or GMX_DOUBLE_EPS if the check should be
 *  done in double regardless of Gromacs precision.
 *  
 *  To check if two algorithms produce similar results you will normally need
 *  to relax the tolerance significantly since many operations (e.g. summation)
 *  accumulate floating point errors.
 *
 *  \param f1  First number to compare
 *  \param f2  Second number to compare
 *  \param tol Tolerance to use
 *
 *  \return 1 if the relative difference is within tolerance, 0 if not.
 */
static int
gmx_within_tol(double   f1,
               double   f2,
               double   tol)
{
    /* The or-equal is important - otherwise we return false if f1==f2==0 */
    if( fabs(f1-f2) <= tol*0.5*(fabs(f1)+fabs(f2)) )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}



/** 
 * Check if a number is smaller than some preset safe minimum
 * value, currently defined as GMX_REAL_MIN/GMX_REAL_EPS.
 *
 * If a number is smaller than this value we risk numerical overflow
 * if any number larger than 1.0/GMX_REAL_EPS is divided by it.
 *
 * \return 1  if 'almost' numerically zero, 0 otherwise.
 */
static int
gmx_numzero(double a)
{
  return gmx_within_tol(a,0.0,GMX_REAL_MIN/GMX_REAL_EPS);
}


static real
gmx_log2(real x)
{
  const real iclog2 = 1.0/log( 2.0 );

    return log( x ) * iclog2;
}


#ifdef __cplusplus
}
#endif

#endif	/* _maths_h */
