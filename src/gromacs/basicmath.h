/*
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
 */
#ifndef GMX_BASICMATH_H
#define GMX_BASICMATH_H

#include <cmath>

#include "types/simple.h"

namespace gmx
{

template <typename T> static inline
T sqr(const T &value)
{
    return value * value;
}

/*! \brief
 * Check if two numbers are within a tolerance
 *
 * This routine checks if the relative difference between two numbers is
 * approximately within the given tolerance, defined as
 * fabs(f1-f2)<=tolerance*fabs(f1+f2).
 *
 * To check if two floating-point numbers are almost identical, use this routine
 * with the tolerance GMX_REAL_EPS, or GMX_DOUBLE_EPS if the check should be
 * done in double regardless of Gromacs precision.
 *
 * To check if two algorithms produce similar results you will normally need
 * to relax the tolerance significantly since many operations (e.g. summation)
 * accumulate floating point errors.
 *
 * \param f1  First number to compare
 * \param f2  Second number to compare
 * \param tol Tolerance to use
 *
 * \return 1 if the relative difference is within tolerance, 0 if not.
 */
static int
gmx_within_tol(double   f1,
               double   f2,
               double   tol)
{
    /* The or-equal is important - otherwise we return false if f1==f2==0 */
    if (std::fabs(f1-f2) <= tol * 0.5 * (std::fabs(f1) + std::fabs(f2)))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}



/*! \brief
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

}

#endif
