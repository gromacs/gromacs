/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_MATH_UTILITIES_H
#define GMX_MATH_UTILITIES_H

#include <limits.h>
#include <stdint.h>

#include <cmath>

/*! \brief Enum to select safe or highly unsafe (faster) math functions.
 *
 *  Normally all the Gromacs math functions should apply reasonable care with
 *  input arguments. While we do not necessarily adhere strictly to IEEE
 *  (in particular not for arguments that might result in NaN, inf, etc.), the
 *  functions should return reasonable values or e.g. clamp results to zero.
 *
 *  However, in a few cases where we are extremely performance-sensitive it
 *  makes sense to forego these checks too in cases where we know the exact
 *  properties if the input data, and we really need to save every cycle we can.
 *
 *  This class is typically used as a template parameter to such calls to enable
 *  the caller to select the level of aggressiveness. We should always use the
 *  safe alternative as the default value, and document carefully what might
 *  happen with the unsafe alternative.
 */
enum class MathOptimization
{
    Safe,  //!< Don't do unsafe optimizations. This should always be default.
    Unsafe //!< Allow optimizations that can be VERY dangerous for general code.
};

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
bool gmx_within_tol(double f1, double f2, double tol);

/*!
 * \brief Check if a number is smaller than some preset safe minimum
 * value, currently defined as GMX_REAL_MIN/GMX_REAL_EPS.
 *
 * If a number is smaller than this value we risk numerical overflow
 * if any number larger than 1.0/GMX_REAL_EPS is divided by it.
 *
 * \return True if 'almost' numerically zero, false otherwise.
 */
bool gmx_numzero(double a);

/*! \brief Multiply two large ints
 *
 * \return False iff overflow occurred
 */
bool check_int_multiply_for_overflow(int64_t a, int64_t b, int64_t* result);

/*! \brief Enable floating-point exceptions if supported on OS
 *
 * Enables division-by-zero, invalid value, and overflow.
 *
 * \returns 0 if successful in enabling exceptions, anything else in case of failure/unsupported OS.
 */
int gmx_feenableexcept();

/*! \brief Disable floating-point exceptions if supported on OS
 *
 * Disables division-by-zero, invalid value, and overflow.
 *
 * \returns 0 if successful in disabling exceptions, anything else in case of failure/unsupported OS.
 */
int gmx_fedisableexcept();

/*! \brief Return true if the current build should enable floating-point exceptions by default.
 *
 * Currently, it returns true unless any of the following conditions are met:
 * - release build,
 * - SYCL build (Intel IGC, at least 1.0.5964, raises FP exceptions in JIT compilation),
 * - - See https://github.com/intel/intel-graphics-compiler/issues/164
 * - compilers with known buggy FP exception support (clang with any optimization)
 *   or suspected buggy FP exception support (gcc 7.* with optimization).
 *
 * Note that this function does not check whether the build/OS supports FP exceptions.
 *
 * \returns true if we should enable FP exceptions by default.
 */
bool gmxShouldEnableFPExceptions();

#endif
