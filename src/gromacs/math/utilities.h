/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#ifndef GMX_MATH_UTILITIES_H
#define GMX_MATH_UTILITIES_H

#include <limits.h>
#include <math.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923
#endif

#ifndef M_2PI
#define M_2PI       6.28318530717958647692
#endif

#ifndef M_SQRT2
#define M_SQRT2 sqrt(2.0)
#endif

#ifndef M_1_PI
#define M_1_PI      0.31830988618379067154
#endif

#ifndef M_FLOAT_1_SQRTPI /* used in GPU kernels */
/* 1.0 / sqrt(M_PI) */
#define M_FLOAT_1_SQRTPI 0.564189583547756f
#endif

#ifndef M_1_SQRTPI
/* 1.0 / sqrt(M_PI) */
#define M_1_SQRTPI 0.564189583547756
#endif

#ifndef M_2_SQRTPI
/* 2.0 / sqrt(M_PI) */
#define M_2_SQRTPI  1.128379167095513
#endif

int     gmx_nint(real a);
real    sign(real x, real y);

real    cuberoot (real a);
double  gmx_erfd(double x);
double  gmx_erfcd(double x);
float   gmx_erff(float x);
float   gmx_erfcf(float x);
#ifdef GMX_DOUBLE
#define gmx_erf(x)   gmx_erfd(x)
#define gmx_erfc(x)  gmx_erfcd(x)
#else
#define gmx_erf(x)   gmx_erff(x)
#define gmx_erfc(x)  gmx_erfcf(x)
#endif

#if defined(_MSC_VER) && _MSC_VER < 1800
#define gmx_expm1(x) (exp(x)-1)
#define gmx_log1p(x) log(1+x)
#else
#define gmx_expm1 expm1
#define gmx_log1p log1p
#endif

gmx_bool gmx_isfinite(real x);
gmx_bool gmx_isnan(real x);

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
int
gmx_within_tol(double   f1,
               double   f2,
               double   tol);

/*!
 * \brief Check if a number is smaller than some preset safe minimum
 * value, currently defined as GMX_REAL_MIN/GMX_REAL_EPS.
 *
 * If a number is smaller than this value we risk numerical overflow
 * if any number larger than 1.0/GMX_REAL_EPS is divided by it.
 *
 * \return 1  if 'almost' numerically zero, 0 otherwise.
 */
int
gmx_numzero(double a);

/*! \brief Compute floor of logarithm to base 2
 *
 * \return log2(x)
 */
unsigned int
gmx_log2i(unsigned int x);

/*! \brief Multiply two large ints
 *
 * \return False iff overflow occured
 */
gmx_bool
check_int_multiply_for_overflow(gmx_int64_t  a,
                                gmx_int64_t  b,
                                gmx_int64_t *result);

/*! \brief Find greatest common divisor of two numbers
 *
 * \return GCD of the two inputs
 */
int
gmx_greatest_common_divisor(int p, int q);


/*! \brief Enable floating-point exceptions if supported on OS
 *
 * Enables division-by-zero, invalid, and overflow.
 */
int gmx_feenableexcept();

#ifdef __cplusplus
}
#endif

#endif
