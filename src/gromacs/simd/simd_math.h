/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
#ifndef GMX_SIMD_SIMD_MATH_H_
#define GMX_SIMD_SIMD_MATH_H_

/*! \libinternal \file
 *
 * \brief Math functions for SIMD datatypes.
 *
 * \attention This file is generic for all SIMD architectures, so you cannot
 * assume that any of the optional SIMD features (as defined in simd.h) are
 * present. In particular, this means you cannot assume support for integers,
 * logical operations (neither on floating-point nor integer values), shifts,
 * and the architecture might only have SIMD for either float or double.
 * Second, to keep this file clean and general, any additions to this file
 * must work for all possible SIMD architectures in both single and double
 * precision (if they support it), and you cannot make any assumptions about
 * SIMD width.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

#include "config.h"

#include <math.h>

#include "gromacs/math/utilities.h"
#include "gromacs/simd/simd.h"

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name Implementation accuracy settings
 *  \{
 */

/*! \} */

#ifdef GMX_SIMD_HAVE_FLOAT

/*! \name Single precision SIMD math functions
 *
 *  \note In most cases you should use the real-precision functions instead.
 *  \{
 */

/****************************************
 * SINGLE PRECISION SIMD MATH FUNCTIONS *
 ****************************************/

/*! \brief SIMD float utility to sum a+b+c+d.
 *
 * You should normally call the real-precision routine \ref gmx_simd_sum4_r.
 *
 * \param a term 1 (multiple values)
 * \param b term 2 (multiple values)
 * \param c term 3 (multiple values)
 * \param d term 4 (multiple values)
 * \return sum of terms 1-4 (multiple values)
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_sum4_f(gmx_simd_float_t a, gmx_simd_float_t b,
                gmx_simd_float_t c, gmx_simd_float_t d)
{
    return gmx_simd_add_f(gmx_simd_add_f(a, b), gmx_simd_add_f(c, d));
}

/*! \brief Return -a if b is negative, SIMD float.
 *
 * You should normally call the real-precision routine \ref gmx_simd_xor_sign_r.
 *
 * \param a Values to set sign for
 * \param b Values used to set sign
 * \return if b is negative, the sign of a will be changed.
 *
 * This is equivalent to doing an xor operation on a with the sign bit of b,
 * with the exception that negative zero is not considered to be negative
 * on architectures where \ref GMX_SIMD_HAVE_LOGICAL is not set.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_xor_sign_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
#ifdef GMX_SIMD_HAVE_LOGICAL
    return gmx_simd_xor_f(a, gmx_simd_and_f(gmx_simd_set1_f(GMX_FLOAT_NEGZERO), b));
#else
    return gmx_simd_blendv_f(a, gmx_simd_fneg_f(a), gmx_simd_cmplt_f(b, gmx_simd_setzero_f()));
#endif
}

#ifndef gmx_simd_rsqrt_iter_f
/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD float.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the inverse square root.
 *
 *  \param lu Approximation of 1/sqrt(x), typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/sqrt(x).
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_rsqrt_iter_f(gmx_simd_float_t lu, gmx_simd_float_t x)
{
#    ifdef GMX_SIMD_HAVE_FMA
    return gmx_simd_fmadd_f(gmx_simd_fnmadd_f(x, gmx_simd_mul_f(lu, lu), gmx_simd_set1_f(1.0f)), gmx_simd_mul_f(lu, gmx_simd_set1_f(0.5f)), lu);
#    else
    return gmx_simd_mul_f(gmx_simd_set1_f(0.5f), gmx_simd_mul_f(gmx_simd_sub_f(gmx_simd_set1_f(3.0f), gmx_simd_mul_f(gmx_simd_mul_f(lu, lu), x)), lu));
#    endif
}
#endif

/*! \brief Calculate 1/sqrt(x) for SIMD float.
 *
 * You should normally call the real-precision routine \ref gmx_simd_invsqrt_r.
 *
 *  \param x Argument that must be >0. This routine does not check arguments.
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_invsqrt_f(gmx_simd_float_t x)
{
    gmx_simd_float_t lu = gmx_simd_rsqrt_f(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rsqrt_iter_f(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rsqrt_iter_f(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rsqrt_iter_f(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for masked entries of SIMD float.
 *
 * Identical to gmx_simd_invsqrt_f but avoids fp-exception for non-masked entries.
 * The result for the non-masked entries is undefined and the user has to use blend
 * with the same mask to obtain a defined result.
 *
 *  \param x Argument that must be >0 for masked entries
 *  \param m Masked entries
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or entry was not masked.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_invsqrt_maskfpe_f(gmx_simd_float_t x, gmx_simd_fbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_invsqrt_f(x);
#else
    return gmx_simd_invsqrt_f(gmx_simd_blendv_f(gmx_simd_set1_f(1.0f), x, m));
#endif
}

/*! \brief Calculate 1/sqrt(x) for non-masked entries of SIMD float.
 *
 * Identical to gmx_simd_invsqrt_f but avoids fp-exception for masked entries.
 * The result for the non-masked entries is undefined and the user has to use blend
 * with the same mask to obtain a defined result.
 *
 *  \param x Argument that must be >0 for non-masked entries
 *  \param m Masked entries
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or entry was masked.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_invsqrt_notmaskfpe_f(gmx_simd_float_t x, gmx_simd_fbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_invsqrt_f(x);
#else
    return gmx_simd_invsqrt_f(gmx_simd_blendv_f(x, gmx_simd_set1_f(1.0f), m));
#endif
}

/*! \brief Calculate 1/sqrt(x) for two SIMD floats.
 *
 * You should normally call the real-precision routine \ref gmx_simd_invsqrt_pair_r.
 *
 * \param x0  First set of arguments, x0 must be positive - no argument checking.
 * \param x1  Second set of arguments, x1 must be positive - no argument checking.
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 *  In particular for double precision we can sometimes calculate square root
 *  pairs slightly faster by using single precision until the very last step.
 */
static gmx_inline void gmx_simdcall
gmx_simd_invsqrt_pair_f(gmx_simd_float_t x0,    gmx_simd_float_t x1,
                        gmx_simd_float_t *out0, gmx_simd_float_t *out1)
{
    *out0 = gmx_simd_invsqrt_f(x0);
    *out1 = gmx_simd_invsqrt_f(x1);
}

#ifndef gmx_simd_rcp_iter_f
/*! \brief Perform one Newton-Raphson iteration to improve 1/x for SIMD float.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the reciprocal.
 *
 *  \param lu Approximation of 1/x, typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/x.
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_rcp_iter_f(gmx_simd_float_t lu, gmx_simd_float_t x)
{
    return gmx_simd_mul_f(lu, gmx_simd_fnmadd_f(lu, x, gmx_simd_set1_f(2.0f)));
}
#endif

/*! \brief Calculate 1/x for SIMD float.
 *
 * You should normally call the real-precision routine \ref gmx_simd_inv_r.
 *
 *  \param x Argument that must be nonzero. This routine does not check arguments.
 *  \return 1/x. Result is undefined if your argument was invalid.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_inv_f(gmx_simd_float_t x)
{
    gmx_simd_float_t lu = gmx_simd_rcp_f(x);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rcp_iter_f(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rcp_iter_f(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rcp_iter_f(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/x for masked entries of SIMD float.
 *
 * Identical to gmx_simd_inv_f but avoids fp-exception for non-masked entries.
 * The result for the non-masked entries is undefined and the user has to use blend
 * with the same mask to obtain a defined result.
 *
 *  \param x Argument that must be nonzero for masked entries
 *  \param m Masked entries
 *  \return 1/x. Result is undefined if your argument was invalid or entry was not masked.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_inv_maskfpe_f(gmx_simd_float_t x, gmx_simd_fbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_inv_f(x);
#else
    return gmx_simd_inv_f(gmx_simd_blendv_f(gmx_simd_set1_f(1.0f), x, m));
#endif
}


/*! \brief Calculate 1/x for non-masked entries of SIMD float.
 *
 * Identical to gmx_simd_inv_f but avoids fp-exception for masked entries.
 * The result for the non-masked entries is undefined and the user has to use blend
 * with the same mask to obtain a defined result.
 *
 *  \param x Argument that must be nonzero for non-masked entries
 *  \param m Masked entries
 *  \return 1/x. Result is undefined if your argument was invalid or entry was masked.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_inv_notmaskfpe_f(gmx_simd_float_t x, gmx_simd_fbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_inv_f(x);
#else
    return gmx_simd_inv_f(gmx_simd_blendv_f(x, gmx_simd_set1_f(1.0f), m));
#endif
}

/*! \brief Calculate sqrt(x) correctly for SIMD floats, including argument 0.0.
 *
 * You should normally call the real-precision routine \ref gmx_simd_sqrt_r.
 *
 *  \param x Argument that must be >=0.
 *  \return sqrt(x). If x=0, the result will correctly be set to 0.
 *          The result is undefined if the input value is negative.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_sqrt_f(gmx_simd_float_t x)
{
    gmx_simd_fbool_t  mask;
    gmx_simd_float_t  res;

    mask = gmx_simd_cmpeq_f(x, gmx_simd_setzero_f());
    res  = gmx_simd_blendnotzero_f(gmx_simd_invsqrt_notmaskfpe_f(x, mask), mask);
    return gmx_simd_mul_f(res, x);
}

/*! \brief SIMD float log(x). This is the natural logarithm.
 *
 * You should normally call the real-precision routine \ref gmx_simd_log_r.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 */
#ifndef gmx_simd_log_f
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_log_f(gmx_simd_float_t x)
{
    const gmx_simd_float_t  half       = gmx_simd_set1_f(0.5f);
    const gmx_simd_float_t  one        = gmx_simd_set1_f(1.0f);
    const gmx_simd_float_t  sqrt2      = gmx_simd_set1_f(sqrt(2.0f));
    const gmx_simd_float_t  corr       = gmx_simd_set1_f(0.693147180559945286226764f);
    const gmx_simd_float_t  CL9        = gmx_simd_set1_f(0.2371599674224853515625f);
    const gmx_simd_float_t  CL7        = gmx_simd_set1_f(0.285279005765914916992188f);
    const gmx_simd_float_t  CL5        = gmx_simd_set1_f(0.400005519390106201171875f);
    const gmx_simd_float_t  CL3        = gmx_simd_set1_f(0.666666567325592041015625f);
    const gmx_simd_float_t  CL1        = gmx_simd_set1_f(2.0f);
    gmx_simd_float_t        fexp, x2, p;
    gmx_simd_fbool_t        mask;

    fexp  = gmx_simd_get_exponent_f(x);
    x     = gmx_simd_get_mantissa_f(x);

    mask  = gmx_simd_cmplt_f(sqrt2, x);
    /* Adjust to non-IEEE format for x>sqrt(2): exponent += 1, mantissa *= 0.5 */
    fexp  = gmx_simd_add_f(fexp, gmx_simd_blendzero_f(one, mask));
    x     = gmx_simd_mul_f(x, gmx_simd_blendv_f(one, half, mask));

    x     = gmx_simd_mul_f( gmx_simd_sub_f(x, one), gmx_simd_inv_f( gmx_simd_add_f(x, one) ) );
    x2    = gmx_simd_mul_f(x, x);

    p     = gmx_simd_fmadd_f(CL9, x2, CL7);
    p     = gmx_simd_fmadd_f(p, x2, CL5);
    p     = gmx_simd_fmadd_f(p, x2, CL3);
    p     = gmx_simd_fmadd_f(p, x2, CL1);
    p     = gmx_simd_fmadd_f(p, x, gmx_simd_mul_f(corr, fexp));

    return p;
}
#endif

#ifndef gmx_simd_exp2_f
/*! \brief SIMD float 2^x.
 *
 * You should normally call the real-precision routine \ref gmx_simd_exp2_r.
 *
 * \param x Argument.
 * \result 2^x. Undefined if input argument caused overflow.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_exp2_f(gmx_simd_float_t x)
{
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const gmx_simd_float_t  arglimit = gmx_simd_set1_f(126.0f);
    const gmx_simd_float_t  CC6      = gmx_simd_set1_f(0.0001534581200287996416911311);
    const gmx_simd_float_t  CC5      = gmx_simd_set1_f(0.001339993121934088894618990);
    const gmx_simd_float_t  CC4      = gmx_simd_set1_f(0.009618488957115180159497841);
    const gmx_simd_float_t  CC3      = gmx_simd_set1_f(0.05550328776964726865751735);
    const gmx_simd_float_t  CC2      = gmx_simd_set1_f(0.2402264689063408646490722);
    const gmx_simd_float_t  CC1      = gmx_simd_set1_f(0.6931472057372680777553816);
    const gmx_simd_float_t  one      = gmx_simd_set1_f(1.0f);

    gmx_simd_float_t        fexppart;
    gmx_simd_float_t        intpart;
    gmx_simd_float_t        p;
    gmx_simd_fbool_t        valuemask;

    fexppart  = gmx_simd_set_exponent_f(x);
    intpart   = gmx_simd_round_f(x);
    valuemask = gmx_simd_cmple_f(gmx_simd_fabs_f(x), arglimit);
    fexppart  = gmx_simd_blendzero_f(fexppart, valuemask);
    x         = gmx_simd_sub_f(x, intpart);

    p         = gmx_simd_fmadd_f(CC6, x, CC5);
    p         = gmx_simd_fmadd_f(p, x, CC4);
    p         = gmx_simd_fmadd_f(p, x, CC3);
    p         = gmx_simd_fmadd_f(p, x, CC2);
    p         = gmx_simd_fmadd_f(p, x, CC1);
    p         = gmx_simd_fmadd_f(p, x, one);
    x         = gmx_simd_mul_f(p, fexppart);
    return x;
}
#endif

#ifndef gmx_simd_exp_f
/*! \brief SIMD float exp(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_exp_r.
 *
 * In addition to scaling the argument for 2^x this routine correctly does
 * extended precision arithmetics to improve accuracy.
 *
 * \param x Argument.
 * \result exp(x). Undefined if input argument caused overflow,
 * which can happen if abs(x) \> 7e13.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_exp_f(gmx_simd_float_t x)
{
    const gmx_simd_float_t  argscale     = gmx_simd_set1_f(1.44269504088896341f);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const gmx_simd_float_t  arglimit     = gmx_simd_set1_f(126.0f);
    const gmx_simd_float_t  invargscale0 = gmx_simd_set1_f(-0.693145751953125f);
    const gmx_simd_float_t  invargscale1 = gmx_simd_set1_f(-1.428606765330187045e-06f);
    const gmx_simd_float_t  CC4          = gmx_simd_set1_f(0.00136324646882712841033936f);
    const gmx_simd_float_t  CC3          = gmx_simd_set1_f(0.00836596917361021041870117f);
    const gmx_simd_float_t  CC2          = gmx_simd_set1_f(0.0416710823774337768554688f);
    const gmx_simd_float_t  CC1          = gmx_simd_set1_f(0.166665524244308471679688f);
    const gmx_simd_float_t  CC0          = gmx_simd_set1_f(0.499999850988388061523438f);
    const gmx_simd_float_t  one          = gmx_simd_set1_f(1.0f);
    gmx_simd_float_t        fexppart;
    gmx_simd_float_t        intpart;
    gmx_simd_float_t        y, p;
    gmx_simd_fbool_t        valuemask;

    y         = gmx_simd_mul_f(x, argscale);
    fexppart  = gmx_simd_set_exponent_f(y);  /* rounds to nearest int internally */
    intpart   = gmx_simd_round_f(y);         /* use same rounding algorithm here */
    valuemask = gmx_simd_cmple_f(gmx_simd_fabs_f(y), arglimit);
    fexppart  = gmx_simd_blendzero_f(fexppart, valuemask);

    /* Extended precision arithmetics */
    x         = gmx_simd_fmadd_f(invargscale0, intpart, x);
    x         = gmx_simd_fmadd_f(invargscale1, intpart, x);

    p         = gmx_simd_fmadd_f(CC4, x, CC3);
    p         = gmx_simd_fmadd_f(p, x, CC2);
    p         = gmx_simd_fmadd_f(p, x, CC1);
    p         = gmx_simd_fmadd_f(p, x, CC0);
    p         = gmx_simd_fmadd_f(gmx_simd_mul_f(x, x), p, x);
    p         = gmx_simd_add_f(p, one);
    x         = gmx_simd_mul_f(p, fexppart);
    return x;
}
#endif

/*! \brief SIMD float erf(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_erf_r.
 *
 * \param x The value to calculate erf(x) for.
 * \result erf(x)
 *
 * This routine achieves very close to full precision, but we do not care about
 * the last bit or the subnormal result range.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_erf_f(gmx_simd_float_t x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const gmx_simd_float_t  CA6      = gmx_simd_set1_f(7.853861353153693e-5f);
    const gmx_simd_float_t  CA5      = gmx_simd_set1_f(-8.010193625184903e-4f);
    const gmx_simd_float_t  CA4      = gmx_simd_set1_f(5.188327685732524e-3f);
    const gmx_simd_float_t  CA3      = gmx_simd_set1_f(-2.685381193529856e-2f);
    const gmx_simd_float_t  CA2      = gmx_simd_set1_f(1.128358514861418e-1f);
    const gmx_simd_float_t  CA1      = gmx_simd_set1_f(-3.761262582423300e-1f);
    const gmx_simd_float_t  CA0      = gmx_simd_set1_f(1.128379165726710f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const gmx_simd_float_t  CB9      = gmx_simd_set1_f(-0.0018629930017603923f);
    const gmx_simd_float_t  CB8      = gmx_simd_set1_f(0.003909821287598495f);
    const gmx_simd_float_t  CB7      = gmx_simd_set1_f(-0.0052094582210355615f);
    const gmx_simd_float_t  CB6      = gmx_simd_set1_f(0.005685614362160572f);
    const gmx_simd_float_t  CB5      = gmx_simd_set1_f(-0.0025367682853477272f);
    const gmx_simd_float_t  CB4      = gmx_simd_set1_f(-0.010199799682318782f);
    const gmx_simd_float_t  CB3      = gmx_simd_set1_f(0.04369575504816542f);
    const gmx_simd_float_t  CB2      = gmx_simd_set1_f(-0.11884063474674492f);
    const gmx_simd_float_t  CB1      = gmx_simd_set1_f(0.2732120154030589f);
    const gmx_simd_float_t  CB0      = gmx_simd_set1_f(0.42758357702025784f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const gmx_simd_float_t  CC10     = gmx_simd_set1_f(-0.0445555913112064f);
    const gmx_simd_float_t  CC9      = gmx_simd_set1_f(0.21376355144663348f);
    const gmx_simd_float_t  CC8      = gmx_simd_set1_f(-0.3473187200259257f);
    const gmx_simd_float_t  CC7      = gmx_simd_set1_f(0.016690861551248114f);
    const gmx_simd_float_t  CC6      = gmx_simd_set1_f(0.7560973182491192f);
    const gmx_simd_float_t  CC5      = gmx_simd_set1_f(-1.2137903600145787f);
    const gmx_simd_float_t  CC4      = gmx_simd_set1_f(0.8411872321232948f);
    const gmx_simd_float_t  CC3      = gmx_simd_set1_f(-0.08670413896296343f);
    const gmx_simd_float_t  CC2      = gmx_simd_set1_f(-0.27124782687240334f);
    const gmx_simd_float_t  CC1      = gmx_simd_set1_f(-0.0007502488047806069f);
    const gmx_simd_float_t  CC0      = gmx_simd_set1_f(0.5642114853803148f);
    const gmx_simd_float_t  one      = gmx_simd_set1_f(1.0f);
    const gmx_simd_float_t  two      = gmx_simd_set1_f(2.0f);

    gmx_simd_float_t        x2, x4, y;
    gmx_simd_float_t        t, t2, w, w2;
    gmx_simd_float_t        pA0, pA1, pB0, pB1, pC0, pC1;
    gmx_simd_float_t        expmx2;
    gmx_simd_float_t        res_erf, res_erfc, res;
    gmx_simd_fbool_t        mask, msk_erf;

    /* Calculate erf() */
    x2   = gmx_simd_mul_f(x, x);
    x4   = gmx_simd_mul_f(x2, x2);

    pA0  = gmx_simd_fmadd_f(CA6, x4, CA4);
    pA1  = gmx_simd_fmadd_f(CA5, x4, CA3);
    pA0  = gmx_simd_fmadd_f(pA0, x4, CA2);
    pA1  = gmx_simd_fmadd_f(pA1, x4, CA1);
    pA0  = gmx_simd_mul_f(pA0, x4);
    pA0  = gmx_simd_fmadd_f(pA1, x2, pA0);
    /* Constant term must come last for precision reasons */
    pA0  = gmx_simd_add_f(pA0, CA0);

    res_erf = gmx_simd_mul_f(x, pA0);

    /* Calculate erfc */
    y       = gmx_simd_fabs_f(x);
    msk_erf = gmx_simd_cmplt_f(y, gmx_simd_set1_f(0.75f));
    t       = gmx_simd_inv_notmaskfpe_f(y, msk_erf);
    w       = gmx_simd_sub_f(t, one);
    t2      = gmx_simd_mul_f(t, t);
    w2      = gmx_simd_mul_f(w, w);

    /* No need for a floating-point sieve here (as in erfc), since erf()
     * will never return values that are extremely small for large args.
     */
    expmx2  = gmx_simd_exp_f( gmx_simd_fneg_f( gmx_simd_mul_f(y, y)));

    pB1  = gmx_simd_fmadd_f(CB9, w2, CB7);
    pB0  = gmx_simd_fmadd_f(CB8, w2, CB6);
    pB1  = gmx_simd_fmadd_f(pB1, w2, CB5);
    pB0  = gmx_simd_fmadd_f(pB0, w2, CB4);
    pB1  = gmx_simd_fmadd_f(pB1, w2, CB3);
    pB0  = gmx_simd_fmadd_f(pB0, w2, CB2);
    pB1  = gmx_simd_fmadd_f(pB1, w2, CB1);
    pB0  = gmx_simd_fmadd_f(pB0, w2, CB0);
    pB0  = gmx_simd_fmadd_f(pB1, w, pB0);

    pC0  = gmx_simd_fmadd_f(CC10, t2, CC8);
    pC1  = gmx_simd_fmadd_f(CC9, t2, CC7);
    pC0  = gmx_simd_fmadd_f(pC0, t2, CC6);
    pC1  = gmx_simd_fmadd_f(pC1, t2, CC5);
    pC0  = gmx_simd_fmadd_f(pC0, t2, CC4);
    pC1  = gmx_simd_fmadd_f(pC1, t2, CC3);
    pC0  = gmx_simd_fmadd_f(pC0, t2, CC2);
    pC1  = gmx_simd_fmadd_f(pC1, t2, CC1);

    pC0  = gmx_simd_fmadd_f(pC0, t2, CC0);
    pC0  = gmx_simd_fmadd_f(pC1, t, pC0);
    pC0  = gmx_simd_mul_f(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = gmx_simd_cmplt_f(two, y);
    res_erfc = gmx_simd_blendv_f(pB0, pC0, mask);
    res_erfc = gmx_simd_mul_f(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = gmx_simd_cmplt_f(x, gmx_simd_setzero_f());
    res_erfc = gmx_simd_blendv_f(res_erfc, gmx_simd_sub_f(two, res_erfc), mask);

    /* Select erf() or erfc() */
    res  = gmx_simd_blendv_f(gmx_simd_sub_f(one, res_erfc), res_erf, msk_erf);

    return res;
}

/*! \brief SIMD float erfc(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_erfc_r.
 *
 * \param x The value to calculate erfc(x) for.
 * \result erfc(x)
 *
 * This routine achieves full precision (bar the last bit) over most of the
 * input range, but for large arguments where the result is getting close
 * to the minimum representable numbers we accept slightly larger errors
 * (think results that are in the ballpark of 10^-30 for single precision,
 * or 10^-200 for double) since that is not relevant for MD.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_erfc_f(gmx_simd_float_t x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const gmx_simd_float_t  CA6      = gmx_simd_set1_f(7.853861353153693e-5f);
    const gmx_simd_float_t  CA5      = gmx_simd_set1_f(-8.010193625184903e-4f);
    const gmx_simd_float_t  CA4      = gmx_simd_set1_f(5.188327685732524e-3f);
    const gmx_simd_float_t  CA3      = gmx_simd_set1_f(-2.685381193529856e-2f);
    const gmx_simd_float_t  CA2      = gmx_simd_set1_f(1.128358514861418e-1f);
    const gmx_simd_float_t  CA1      = gmx_simd_set1_f(-3.761262582423300e-1f);
    const gmx_simd_float_t  CA0      = gmx_simd_set1_f(1.128379165726710f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const gmx_simd_float_t  CB9      = gmx_simd_set1_f(-0.0018629930017603923f);
    const gmx_simd_float_t  CB8      = gmx_simd_set1_f(0.003909821287598495f);
    const gmx_simd_float_t  CB7      = gmx_simd_set1_f(-0.0052094582210355615f);
    const gmx_simd_float_t  CB6      = gmx_simd_set1_f(0.005685614362160572f);
    const gmx_simd_float_t  CB5      = gmx_simd_set1_f(-0.0025367682853477272f);
    const gmx_simd_float_t  CB4      = gmx_simd_set1_f(-0.010199799682318782f);
    const gmx_simd_float_t  CB3      = gmx_simd_set1_f(0.04369575504816542f);
    const gmx_simd_float_t  CB2      = gmx_simd_set1_f(-0.11884063474674492f);
    const gmx_simd_float_t  CB1      = gmx_simd_set1_f(0.2732120154030589f);
    const gmx_simd_float_t  CB0      = gmx_simd_set1_f(0.42758357702025784f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const gmx_simd_float_t  CC10     = gmx_simd_set1_f(-0.0445555913112064f);
    const gmx_simd_float_t  CC9      = gmx_simd_set1_f(0.21376355144663348f);
    const gmx_simd_float_t  CC8      = gmx_simd_set1_f(-0.3473187200259257f);
    const gmx_simd_float_t  CC7      = gmx_simd_set1_f(0.016690861551248114f);
    const gmx_simd_float_t  CC6      = gmx_simd_set1_f(0.7560973182491192f);
    const gmx_simd_float_t  CC5      = gmx_simd_set1_f(-1.2137903600145787f);
    const gmx_simd_float_t  CC4      = gmx_simd_set1_f(0.8411872321232948f);
    const gmx_simd_float_t  CC3      = gmx_simd_set1_f(-0.08670413896296343f);
    const gmx_simd_float_t  CC2      = gmx_simd_set1_f(-0.27124782687240334f);
    const gmx_simd_float_t  CC1      = gmx_simd_set1_f(-0.0007502488047806069f);
    const gmx_simd_float_t  CC0      = gmx_simd_set1_f(0.5642114853803148f);
    /* Coefficients for expansion of exp(x) in [0,0.1] */
    /* CD0 and CD1 are both 1.0, so no need to declare them separately */
    const gmx_simd_float_t  CD2      = gmx_simd_set1_f(0.5000066608081202f);
    const gmx_simd_float_t  CD3      = gmx_simd_set1_f(0.1664795422874624f);
    const gmx_simd_float_t  CD4      = gmx_simd_set1_f(0.04379839977652482f);
    const gmx_simd_float_t  one      = gmx_simd_set1_f(1.0f);
    const gmx_simd_float_t  two      = gmx_simd_set1_f(2.0f);

    /* We need to use a small trick here, since we cannot assume all SIMD
     * architectures support integers, and the flag we want (0xfffff000) would
     * evaluate to NaN (i.e., it cannot be expressed as a floating-point num).
     * Instead, we represent the flags 0xf0f0f000 and 0x0f0f0000 as valid
     * fp numbers, and perform a logical or. Since the expression is constant,
     * we can at least hope it is evaluated at compile-time.
     */
#ifdef GMX_SIMD_HAVE_LOGICAL
    const gmx_simd_float_t  sieve    = gmx_simd_or_f(gmx_simd_set1_f(-5.965323564e+29f), gmx_simd_set1_f(7.05044434e-30f));
#else
    const int               isieve   = 0xFFFFF000;
    float                   mem[GMX_SIMD_REAL_WIDTH*2];
    float *                 pmem = gmx_simd_align_f(mem);
    union {
        float f; int i;
    } conv;
    int                     i;
#endif

    gmx_simd_float_t        x2, x4, y;
    gmx_simd_float_t        q, z, t, t2, w, w2;
    gmx_simd_float_t        pA0, pA1, pB0, pB1, pC0, pC1;
    gmx_simd_float_t        expmx2, corr;
    gmx_simd_float_t        res_erf, res_erfc, res;
    gmx_simd_fbool_t        mask, msk_erf;

    /* Calculate erf() */
    x2     = gmx_simd_mul_f(x, x);
    x4     = gmx_simd_mul_f(x2, x2);

    pA0  = gmx_simd_fmadd_f(CA6, x4, CA4);
    pA1  = gmx_simd_fmadd_f(CA5, x4, CA3);
    pA0  = gmx_simd_fmadd_f(pA0, x4, CA2);
    pA1  = gmx_simd_fmadd_f(pA1, x4, CA1);
    pA1  = gmx_simd_mul_f(pA1, x2);
    pA0  = gmx_simd_fmadd_f(pA0, x4, pA1);
    /* Constant term must come last for precision reasons */
    pA0  = gmx_simd_add_f(pA0, CA0);

    res_erf = gmx_simd_mul_f(x, pA0);

    /* Calculate erfc */
    y       = gmx_simd_fabs_f(x);
    msk_erf = gmx_simd_cmplt_f(y, gmx_simd_set1_f(0.75f));
    t       = gmx_simd_inv_notmaskfpe_f(y, msk_erf);
    w       = gmx_simd_sub_f(t, one);
    t2      = gmx_simd_mul_f(t, t);
    w2      = gmx_simd_mul_f(w, w);
    /*
     * We cannot simply calculate exp(-y2) directly in single precision, since
     * that will lose a couple of bits of precision due to the multiplication.
     * Instead, we introduce y=z+w, where the last 12 bits of precision are in w.
     * Then we get exp(-y2) = exp(-z2)*exp((z-y)*(z+y)).
     *
     * The only drawback with this is that it requires TWO separate exponential
     * evaluations, which would be horrible performance-wise. However, the argument
     * for the second exp() call is always small, so there we simply use a
     * low-order minimax expansion on [0,0.1].
     *
     * However, this neat idea requires support for logical ops (and) on
     * FP numbers, which some vendors decided isn't necessary in their SIMD
     * instruction sets (Hi, IBM VSX!). In principle we could use some tricks
     * in double, but we still need memory as a backup when that is not available,
     * and this case is rare enough that we go directly there...
     */
#ifdef GMX_SIMD_HAVE_LOGICAL
    z       = gmx_simd_and_f(y, sieve);
#else
    gmx_simd_store_f(pmem, y);
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f  = pmem[i];
        conv.i  = conv.i & isieve;
        pmem[i] = conv.f;
    }
    z = gmx_simd_load_f(pmem);
#endif
    q       = gmx_simd_mul_f( gmx_simd_sub_f(z, y), gmx_simd_add_f(z, y) );
    corr    = gmx_simd_fmadd_f(CD4, q, CD3);
    corr    = gmx_simd_fmadd_f(corr, q, CD2);
    corr    = gmx_simd_fmadd_f(corr, q, one);
    corr    = gmx_simd_fmadd_f(corr, q, one);

    expmx2  = gmx_simd_exp_f( gmx_simd_fneg_f( gmx_simd_mul_f(z, z) ) );
    expmx2  = gmx_simd_mul_f(expmx2, corr);

    pB1  = gmx_simd_fmadd_f(CB9, w2, CB7);
    pB0  = gmx_simd_fmadd_f(CB8, w2, CB6);
    pB1  = gmx_simd_fmadd_f(pB1, w2, CB5);
    pB0  = gmx_simd_fmadd_f(pB0, w2, CB4);
    pB1  = gmx_simd_fmadd_f(pB1, w2, CB3);
    pB0  = gmx_simd_fmadd_f(pB0, w2, CB2);
    pB1  = gmx_simd_fmadd_f(pB1, w2, CB1);
    pB0  = gmx_simd_fmadd_f(pB0, w2, CB0);
    pB0  = gmx_simd_fmadd_f(pB1, w, pB0);

    pC0  = gmx_simd_fmadd_f(CC10, t2, CC8);
    pC1  = gmx_simd_fmadd_f(CC9, t2, CC7);
    pC0  = gmx_simd_fmadd_f(pC0, t2, CC6);
    pC1  = gmx_simd_fmadd_f(pC1, t2, CC5);
    pC0  = gmx_simd_fmadd_f(pC0, t2, CC4);
    pC1  = gmx_simd_fmadd_f(pC1, t2, CC3);
    pC0  = gmx_simd_fmadd_f(pC0, t2, CC2);
    pC1  = gmx_simd_fmadd_f(pC1, t2, CC1);

    pC0  = gmx_simd_fmadd_f(pC0, t2, CC0);
    pC0  = gmx_simd_fmadd_f(pC1, t, pC0);
    pC0  = gmx_simd_mul_f(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = gmx_simd_cmplt_f(two, y);
    res_erfc = gmx_simd_blendv_f(pB0, pC0, mask);
    res_erfc = gmx_simd_mul_f(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = gmx_simd_cmplt_f(x, gmx_simd_setzero_f());
    res_erfc = gmx_simd_blendv_f(res_erfc, gmx_simd_sub_f(two, res_erfc), mask);

    /* Select erf() or erfc() */
    res  = gmx_simd_blendv_f(res_erfc, gmx_simd_sub_f(one, res_erf), msk_erf);

    return res;
}

/*! \brief SIMD float sin \& cos.
 *
 * You should normally call the real-precision routine \ref gmx_simd_sincos_r.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 * This version achieves close to machine precision, but for very large
 * magnitudes of the argument we inherently begin to lose accuracy due to the
 * argument reduction, despite using extended precision arithmetics internally.
 */
static gmx_inline void gmx_simdcall
gmx_simd_sincos_f(gmx_simd_float_t x, gmx_simd_float_t *sinval, gmx_simd_float_t *cosval)
{
    /* Constants to subtract Pi/4*x from y while minimizing precision loss */
    const gmx_simd_float_t  argred0         = gmx_simd_set1_f(-1.5703125);
    const gmx_simd_float_t  argred1         = gmx_simd_set1_f(-4.83751296997070312500e-04f);
    const gmx_simd_float_t  argred2         = gmx_simd_set1_f(-7.54953362047672271729e-08f);
    const gmx_simd_float_t  argred3         = gmx_simd_set1_f(-2.56334406825708960298e-12f);
    const gmx_simd_float_t  two_over_pi     = gmx_simd_set1_f(2.0f/M_PI);
    const gmx_simd_float_t  const_sin2      = gmx_simd_set1_f(-1.9515295891e-4f);
    const gmx_simd_float_t  const_sin1      = gmx_simd_set1_f( 8.3321608736e-3f);
    const gmx_simd_float_t  const_sin0      = gmx_simd_set1_f(-1.6666654611e-1f);
    const gmx_simd_float_t  const_cos2      = gmx_simd_set1_f( 2.443315711809948e-5f);
    const gmx_simd_float_t  const_cos1      = gmx_simd_set1_f(-1.388731625493765e-3f);
    const gmx_simd_float_t  const_cos0      = gmx_simd_set1_f( 4.166664568298827e-2f);
    const gmx_simd_float_t  half            = gmx_simd_set1_f(0.5f);
    const gmx_simd_float_t  one             = gmx_simd_set1_f(1.0f);
    gmx_simd_float_t        ssign, csign;
    gmx_simd_float_t        x2, y, z, psin, pcos, sss, ccc;
    gmx_simd_fbool_t        mask;
#if (defined GMX_SIMD_HAVE_FINT32) && (defined GMX_SIMD_HAVE_FINT32_ARITHMETICS) && (defined GMX_SIMD_HAVE_LOGICAL)
    const gmx_simd_fint32_t ione            = gmx_simd_set1_fi(1);
    const gmx_simd_fint32_t itwo            = gmx_simd_set1_fi(2);
    gmx_simd_fint32_t       iy;

    z       = gmx_simd_mul_f(x, two_over_pi);
    iy      = gmx_simd_cvt_f2i(z);
    y       = gmx_simd_round_f(z);

    mask    = gmx_simd_cvt_fib2fb(gmx_simd_cmpeq_fi(gmx_simd_and_fi(iy, ione), gmx_simd_setzero_fi()));
    ssign   = gmx_simd_blendzero_f(gmx_simd_set1_f(GMX_FLOAT_NEGZERO), gmx_simd_cvt_fib2fb(gmx_simd_cmpeq_fi(gmx_simd_and_fi(iy, itwo), itwo)));
    csign   = gmx_simd_blendzero_f(gmx_simd_set1_f(GMX_FLOAT_NEGZERO), gmx_simd_cvt_fib2fb(gmx_simd_cmpeq_fi(gmx_simd_and_fi(gmx_simd_add_fi(iy, ione), itwo), itwo)));
#else
    const gmx_simd_float_t  quarter         = gmx_simd_set1_f(0.25f);
    const gmx_simd_float_t  minusquarter    = gmx_simd_set1_f(-0.25f);
    gmx_simd_float_t        q;
    gmx_simd_fbool_t        m1, m2, m3;

    /* The most obvious way to find the arguments quadrant in the unit circle
     * to calculate the sign is to use integer arithmetic, but that is not
     * present in all SIMD implementations. As an alternative, we have devised a
     * pure floating-point algorithm that uses truncation for argument reduction
     * so that we get a new value 0<=q<1 over the unit circle, and then
     * do floating-point comparisons with fractions. This is likely to be
     * slightly slower (~10%) due to the longer latencies of floating-point, so
     * we only use it when integer SIMD arithmetic is not present.
     */
    ssign   = x;
    x       = gmx_simd_fabs_f(x);
    /* It is critical that half-way cases are rounded down */
    z       = gmx_simd_fmadd_f(x, two_over_pi, half);
    y       = gmx_simd_trunc_f(z);
    q       = gmx_simd_mul_f(z, quarter);
    q       = gmx_simd_sub_f(q, gmx_simd_trunc_f(q));
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q       = gmx_simd_sub_f(q, half);
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#    ifdef GMX_SIMD_HAVE_LOGICAL
    ssign   = gmx_simd_and_f(ssign, gmx_simd_set1_f(GMX_FLOAT_NEGZERO));
    csign   = gmx_simd_andnot_f(q, gmx_simd_set1_f(GMX_FLOAT_NEGZERO));
    ssign   = gmx_simd_xor_f(ssign, csign);
#    else
    csign   = gmx_simd_xor_sign_f(gmx_simd_set1_f(-1.0f), q);
    // ALT: csign = gmx_simd_fneg_f(gmx_simd_copysign(gmx_simd_set1_f(1.0),q));

    ssign   = gmx_simd_xor_sign_f(ssign, csign);    /* swap ssign if csign was set. */
#    endif
    /* Check if y had value 1 or 3 (remember we subtracted 0.5 from q) */
    m1      = gmx_simd_cmplt_f(q, minusquarter);
    m2      = gmx_simd_cmple_f(gmx_simd_setzero_f(), q);
    m3      = gmx_simd_cmplt_f(q, quarter);
    m2      = gmx_simd_and_fb(m2, m3);
    mask    = gmx_simd_or_fb(m1, m2);
    /* where mask is FALSE, set sign. */
    csign   = gmx_simd_xor_sign_f(csign, gmx_simd_blendv_f(gmx_simd_set1_f(-1.0f), one, mask));
#endif
    x       = gmx_simd_fmadd_f(y, argred0, x);
    x       = gmx_simd_fmadd_f(y, argred1, x);
    x       = gmx_simd_fmadd_f(y, argred2, x);
    x       = gmx_simd_fmadd_f(y, argred3, x);
    x2      = gmx_simd_mul_f(x, x);

    psin    = gmx_simd_fmadd_f(const_sin2, x2, const_sin1);
    psin    = gmx_simd_fmadd_f(psin, x2, const_sin0);
    psin    = gmx_simd_fmadd_f(psin, gmx_simd_mul_f(x, x2), x);
    pcos    = gmx_simd_fmadd_f(const_cos2, x2, const_cos1);
    pcos    = gmx_simd_fmadd_f(pcos, x2, const_cos0);
    pcos    = gmx_simd_fmsub_f(pcos, x2, half);
    pcos    = gmx_simd_fmadd_f(pcos, x2, one);

    sss     = gmx_simd_blendv_f(pcos, psin, mask);
    ccc     = gmx_simd_blendv_f(psin, pcos, mask);
    /* See comment for GMX_SIMD_HAVE_LOGICAL section above. */
#ifdef GMX_SIMD_HAVE_LOGICAL
    *sinval = gmx_simd_xor_f(sss, ssign);
    *cosval = gmx_simd_xor_f(ccc, csign);
#else
    *sinval = gmx_simd_xor_sign_f(sss, ssign);
    *cosval = gmx_simd_xor_sign_f(ccc, csign);
#endif
}

/*! \brief SIMD float sin(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_sin_r.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref gmx_simd_sincos_r and waste a factor 2 in performance.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_sin_f(gmx_simd_float_t x)
{
    gmx_simd_float_t s, c;
    gmx_simd_sincos_f(x, &s, &c);
    return s;
}

/*! \brief SIMD float cos(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_cos_r.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref gmx_simd_sincos_r and waste a factor 2 in performance.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_cos_f(gmx_simd_float_t x)
{
    gmx_simd_float_t s, c;
    gmx_simd_sincos_f(x, &s, &c);
    return c;
}

/*! \brief SIMD float tan(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_tan_r.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_tan_f(gmx_simd_float_t x)
{
    const gmx_simd_float_t  argred0         = gmx_simd_set1_f(-1.5703125);
    const gmx_simd_float_t  argred1         = gmx_simd_set1_f(-4.83751296997070312500e-04f);
    const gmx_simd_float_t  argred2         = gmx_simd_set1_f(-7.54953362047672271729e-08f);
    const gmx_simd_float_t  argred3         = gmx_simd_set1_f(-2.56334406825708960298e-12f);
    const gmx_simd_float_t  two_over_pi     = gmx_simd_set1_f(2.0f/M_PI);
    const gmx_simd_float_t  CT6             = gmx_simd_set1_f(0.009498288995810566122993911);
    const gmx_simd_float_t  CT5             = gmx_simd_set1_f(0.002895755790837379295226923);
    const gmx_simd_float_t  CT4             = gmx_simd_set1_f(0.02460087336161924491836265);
    const gmx_simd_float_t  CT3             = gmx_simd_set1_f(0.05334912882656359828045988);
    const gmx_simd_float_t  CT2             = gmx_simd_set1_f(0.1333989091464957704418495);
    const gmx_simd_float_t  CT1             = gmx_simd_set1_f(0.3333307599244198227797507);

    gmx_simd_float_t        x2, p, y, z;
    gmx_simd_fbool_t        mask;

#if (defined GMX_SIMD_HAVE_FINT32) && (defined GMX_SIMD_HAVE_FINT32_ARITHMETICS) && (defined GMX_SIMD_HAVE_LOGICAL)
    gmx_simd_fint32_t  iy;
    gmx_simd_fint32_t  ione = gmx_simd_set1_fi(1);

    z       = gmx_simd_mul_f(x, two_over_pi);
    iy      = gmx_simd_cvt_f2i(z);
    y       = gmx_simd_round_f(z);
    mask    = gmx_simd_cvt_fib2fb(gmx_simd_cmpeq_fi(gmx_simd_and_fi(iy, ione), ione));

    x       = gmx_simd_fmadd_f(y, argred0, x);
    x       = gmx_simd_fmadd_f(y, argred1, x);
    x       = gmx_simd_fmadd_f(y, argred2, x);
    x       = gmx_simd_fmadd_f(y, argred3, x);
    x       = gmx_simd_xor_f(gmx_simd_blendzero_f(gmx_simd_set1_f(GMX_FLOAT_NEGZERO), mask), x);
#else
    const gmx_simd_float_t  quarter         = gmx_simd_set1_f(0.25f);
    const gmx_simd_float_t  half            = gmx_simd_set1_f(0.5f);
    const gmx_simd_float_t  threequarter    = gmx_simd_set1_f(0.75f);
    gmx_simd_float_t        w, q;
    gmx_simd_fbool_t        m1, m2, m3;

    w       = gmx_simd_fabs_f(x);
    z       = gmx_simd_fmadd_f(w, two_over_pi, half);
    y       = gmx_simd_trunc_f(z);
    q       = gmx_simd_mul_f(z, quarter);
    q       = gmx_simd_sub_f(q, gmx_simd_trunc_f(q));
    m1      = gmx_simd_cmple_f(quarter, q);
    m2      = gmx_simd_cmplt_f(q, half);
    m3      = gmx_simd_cmple_f(threequarter, q);
    m1      = gmx_simd_and_fb(m1, m2);
    mask    = gmx_simd_or_fb(m1, m3);
    w       = gmx_simd_fmadd_f(y, argred0, w);
    w       = gmx_simd_fmadd_f(y, argred1, w);
    w       = gmx_simd_fmadd_f(y, argred2, w);
    w       = gmx_simd_fmadd_f(y, argred3, w);

    w       = gmx_simd_blendv_f(w, gmx_simd_fneg_f(w), mask);
    x       = gmx_simd_xor_sign_f(w, x);
#endif
    x2      = gmx_simd_mul_f(x, x);
    p       = gmx_simd_fmadd_f(CT6, x2, CT5);
    p       = gmx_simd_fmadd_f(p, x2, CT4);
    p       = gmx_simd_fmadd_f(p, x2, CT3);
    p       = gmx_simd_fmadd_f(p, x2, CT2);
    p       = gmx_simd_fmadd_f(p, x2, CT1);
    p       = gmx_simd_fmadd_f(x2, gmx_simd_mul_f(p, x), x);

    p       = gmx_simd_blendv_f( p, gmx_simd_inv_maskfpe_f(p, mask), mask);
    return p;
}

/*! \brief SIMD float asin(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_asin_r.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_asin_f(gmx_simd_float_t x)
{
    const gmx_simd_float_t limitlow   = gmx_simd_set1_f(1e-4f);
    const gmx_simd_float_t half       = gmx_simd_set1_f(0.5f);
    const gmx_simd_float_t one        = gmx_simd_set1_f(1.0f);
    const gmx_simd_float_t halfpi     = gmx_simd_set1_f((float)M_PI/2.0f);
    const gmx_simd_float_t CC5        = gmx_simd_set1_f(4.2163199048E-2f);
    const gmx_simd_float_t CC4        = gmx_simd_set1_f(2.4181311049E-2f);
    const gmx_simd_float_t CC3        = gmx_simd_set1_f(4.5470025998E-2f);
    const gmx_simd_float_t CC2        = gmx_simd_set1_f(7.4953002686E-2f);
    const gmx_simd_float_t CC1        = gmx_simd_set1_f(1.6666752422E-1f);
    gmx_simd_float_t       xabs;
    gmx_simd_float_t       z, z1, z2, q, q1, q2;
    gmx_simd_float_t       pA, pB;
    gmx_simd_fbool_t       mask, mask2;

    xabs  = gmx_simd_fabs_f(x);
    mask  = gmx_simd_cmplt_f(half, xabs);
    z1    = gmx_simd_mul_f(half, gmx_simd_sub_f(one, xabs));
    mask2 = gmx_simd_cmpeq_f(xabs, one);
    q1    = gmx_simd_mul_f(z1, gmx_simd_invsqrt_notmaskfpe_f(z1, mask2));
    q1    = gmx_simd_blendnotzero_f(q1, mask2);
    q2    = xabs;
    z2    = gmx_simd_mul_f(q2, q2);
    z     = gmx_simd_blendv_f(z2, z1, mask);
    q     = gmx_simd_blendv_f(q2, q1, mask);

    z2    = gmx_simd_mul_f(z, z);
    pA    = gmx_simd_fmadd_f(CC5, z2, CC3);
    pB    = gmx_simd_fmadd_f(CC4, z2, CC2);
    pA    = gmx_simd_fmadd_f(pA, z2, CC1);
    pA    = gmx_simd_mul_f(pA, z);
    z     = gmx_simd_fmadd_f(pB, z2, pA);
    z     = gmx_simd_fmadd_f(z, q, q);
    q2    = gmx_simd_sub_f(halfpi, z);
    q2    = gmx_simd_sub_f(q2, z);
    z     = gmx_simd_blendv_f(z, q2, mask);

    mask  = gmx_simd_cmplt_f(limitlow, xabs);
    z     = gmx_simd_blendv_f( xabs, z, mask );
    z     = gmx_simd_xor_sign_f(z, x);

    return z;
}

/*! \brief SIMD float acos(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_acos_r.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_acos_f(gmx_simd_float_t x)
{
    const gmx_simd_float_t one       = gmx_simd_set1_f(1.0f);
    const gmx_simd_float_t half      = gmx_simd_set1_f(0.5f);
    const gmx_simd_float_t pi        = gmx_simd_set1_f((float)M_PI);
    const gmx_simd_float_t halfpi    = gmx_simd_set1_f((float)M_PI/2.0f);
    gmx_simd_float_t       xabs;
    gmx_simd_float_t       z, z1, z2, z3;
    gmx_simd_fbool_t       mask1, mask2, mask3;

    xabs  = gmx_simd_fabs_f(x);
    mask1 = gmx_simd_cmplt_f(half, xabs);
    mask2 = gmx_simd_cmplt_f(gmx_simd_setzero_f(), x);

    z     = gmx_simd_mul_f(half, gmx_simd_sub_f(one, xabs));
    mask3 = gmx_simd_cmpeq_f(xabs, one);
    z     = gmx_simd_mul_f(z, gmx_simd_invsqrt_notmaskfpe_f(z, mask3));
    z     = gmx_simd_blendnotzero_f(z, mask3);
    z     = gmx_simd_blendv_f(x, z, mask1);
    z     = gmx_simd_asin_f(z);

    z2    = gmx_simd_add_f(z, z);
    z1    = gmx_simd_sub_f(pi, z2);
    z3    = gmx_simd_sub_f(halfpi, z);
    z     = gmx_simd_blendv_f(z1, z2, mask2);
    z     = gmx_simd_blendv_f(z3, z, mask1);

    return z;
}

/*! \brief SIMD float asin(x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_atan_r.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x), same argument/value range as standard math library.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_atan_f(gmx_simd_float_t x)
{
    const gmx_simd_float_t halfpi    = gmx_simd_set1_f(M_PI/2);
    const gmx_simd_float_t CA17      = gmx_simd_set1_f(0.002823638962581753730774f);
    const gmx_simd_float_t CA15      = gmx_simd_set1_f(-0.01595690287649631500244f);
    const gmx_simd_float_t CA13      = gmx_simd_set1_f(0.04250498861074447631836f);
    const gmx_simd_float_t CA11      = gmx_simd_set1_f(-0.07489009201526641845703f);
    const gmx_simd_float_t CA9       = gmx_simd_set1_f(0.1063479334115982055664f);
    const gmx_simd_float_t CA7       = gmx_simd_set1_f(-0.1420273631811141967773f);
    const gmx_simd_float_t CA5       = gmx_simd_set1_f(0.1999269574880599975585f);
    const gmx_simd_float_t CA3       = gmx_simd_set1_f(-0.3333310186862945556640f);
    const gmx_simd_float_t one       = gmx_simd_set1_f(1.0f);
    gmx_simd_float_t       x2, x3, x4, pA, pB;
    gmx_simd_fbool_t       mask, mask2;

    mask  = gmx_simd_cmplt_f(x, gmx_simd_setzero_f());
    x     = gmx_simd_fabs_f(x);
    mask2 = gmx_simd_cmplt_f(one, x);
    x     = gmx_simd_blendv_f(x, gmx_simd_inv_maskfpe_f(x, mask2), mask2);

    x2    = gmx_simd_mul_f(x, x);
    x3    = gmx_simd_mul_f(x2, x);
    x4    = gmx_simd_mul_f(x2, x2);
    pA    = gmx_simd_fmadd_f(CA17, x4, CA13);
    pB    = gmx_simd_fmadd_f(CA15, x4, CA11);
    pA    = gmx_simd_fmadd_f(pA, x4, CA9);
    pB    = gmx_simd_fmadd_f(pB, x4, CA7);
    pA    = gmx_simd_fmadd_f(pA, x4, CA5);
    pB    = gmx_simd_fmadd_f(pB, x4, CA3);
    pA    = gmx_simd_fmadd_f(pA, x2, pB);
    pA    = gmx_simd_fmadd_f(pA, x3, x);

    pA    = gmx_simd_blendv_f(pA, gmx_simd_sub_f(halfpi, pA), mask2);
    pA    = gmx_simd_blendv_f(pA, gmx_simd_fneg_f(pA), mask);

    return pA;
}

/*! \brief SIMD float atan2(y,x).
 *
 * You should normally call the real-precision routine \ref gmx_simd_atan2_r.
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x), same argument/value range as standard math library.
 *
 * \note This routine should provide correct results for all finite
 * non-zero or positive-zero arguments. However, negative zero arguments will
 * be treated as positive zero, which means the return value will deviate from
 * the standard math library atan2(y,x) for those cases. That should not be
 * of any concern in Gromacs, and in particular it will not affect calculations
 * of angles from vectors.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_atan2_f(gmx_simd_float_t y, gmx_simd_float_t x)
{
    const gmx_simd_float_t pi          = gmx_simd_set1_f(M_PI);
    const gmx_simd_float_t halfpi      = gmx_simd_set1_f(M_PI/2.0);
    gmx_simd_float_t       xinv, p, aoffset;
    gmx_simd_fbool_t       mask_x0, mask_y0, mask_xlt0, mask_ylt0;

    mask_x0   = gmx_simd_cmpeq_f(x, gmx_simd_setzero_f());
    mask_y0   = gmx_simd_cmpeq_f(y, gmx_simd_setzero_f());
    mask_xlt0 = gmx_simd_cmplt_f(x, gmx_simd_setzero_f());
    mask_ylt0 = gmx_simd_cmplt_f(y, gmx_simd_setzero_f());

    aoffset   = gmx_simd_blendzero_f(halfpi, mask_x0);
    aoffset   = gmx_simd_blendnotzero_f(aoffset, mask_y0);

    aoffset   = gmx_simd_blendv_f(aoffset, pi, mask_xlt0);
    aoffset   = gmx_simd_blendv_f(aoffset, gmx_simd_fneg_f(aoffset), mask_ylt0);

    xinv      = gmx_simd_blendnotzero_f(gmx_simd_inv_notmaskfpe_f(x, mask_x0), mask_x0);
    p         = gmx_simd_mul_f(y, xinv);
    p         = gmx_simd_atan_f(p);
    p         = gmx_simd_add_f(p, aoffset);

    return p;
}

/*! \brief Calculate the force correction due to PME analytically in SIMD float.
 *
 * You should normally call the real-precision routine \ref gmx_simd_pmecorrF_r.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb force - see below for details.
 *
 * This routine is meant to enable analytical evaluation of the
 * direct-space PME electrostatic force to avoid tables.
 *
 * The direct-space potential should be \f$ \mbox{erfc}(\beta r)/r\f$, but there
 * are some problems evaluating that:
 *
 * First, the error function is difficult (read: expensive) to
 * approxmiate accurately for intermediate to large arguments, and
 * this happens already in ranges of \f$(\beta r)\f$ that occur in simulations.
 * Second, we now try to avoid calculating potentials in Gromacs but
 * use forces directly.
 *
 * We can simply things slight by noting that the PME part is really
 * a correction to the normal Coulomb force since \f$\mbox{erfc}(z)=1-\mbox{erf}(z)\f$, i.e.
 * \f[
 * V = \frac{1}{r} - \frac{\mbox{erf}(\beta r)}{r}
 * \f]
 * The first term we already have from the inverse square root, so
 * that we can leave out of this routine.
 *
 * For pme tolerances of 1e-3 to 1e-8 and cutoffs of 0.5nm to 1.8nm,
 * the argument \f$beta r\f$ will be in the range 0.15 to ~4, which is
 * the range used for the minimax fit. Use your favorite plotting program
 * to realize how well-behaved \f$\frac{\mbox{erf}(z)}{z}\f$ is in this range!
 *
 * We approximate \f$f(z)=\mbox{erf}(z)/z\f$ with a rational minimax polynomial.
 * However, it turns out it is more efficient to approximate \f$f(z)/z\f$ and
 * then only use even powers. This is another minor optimization, since
 * we actually \a want \f$f(z)/z\f$, because it is going to be multiplied by
 * the vector between the two atoms to get the vectorial force. The
 * fastest flops are the ones we can avoid calculating!
 *
 * So, here's how it should be used:
 *
 * 1. Calculate \f$r^2\f$.
 * 2. Multiply by \f$\beta^2\f$, so you get \f$z^2=(\beta r)^2\f$.
 * 3. Evaluate this routine with \f$z^2\f$ as the argument.
 * 4. The return value is the expression:
 *
 * \f[
 *    \frac{2 \exp{-z^2}}{\sqrt{\pi} z^2}-\frac{\mbox{erf}(z)}{z^3}
 * \f]
 *
 * 5. Multiply the entire expression by \f$\beta^3\f$. This will get you
 *
 *  \f[
 *    \frac{2 \beta^3 \exp(-z^2)}{\sqrt{\pi} z^2} - \frac{\beta^3 \mbox{erf}(z)}{z^3}
 *  \f]
 *
 *    or, switching back to \f$r\f$ (since \f$z=r \beta\f$):
 *
 *  \f[
 *    \frac{2 \beta \exp(-r^2 \beta^2)}{\sqrt{\pi} r^2} - \frac{\mbox{erf}(r \beta)}{r^3}
 *  \f]
 *
 *    With a bit of math exercise you should be able to confirm that
 *    this is exactly
 *
 *  \f[
 *   \frac{\frac{d}{dr}\left( \frac{\mbox{erf}(\beta r)}{r} \right)}{r}
 *  \f]
 *
 * 6. Add the result to \f$r^{-3}\f$, multiply by the product of the charges,
 *    and you have your force (divided by \f$r\f$). A final multiplication
 *    with the vector connecting the two particles and you have your
 *    vectorial force to add to the particles.
 *
 * This approximation achieves an error slightly lower than 1e-6
 * in single precision and 1e-11 in double precision
 * for arguments smaller than 16 (\f$\beta r \leq 4 \f$);
 * when added to \f$1/r\f$ the error will be insignificant.
 * For \f$\beta r \geq 7206\f$ the return value can be inf or NaN.
 *
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_pmecorrF_f(gmx_simd_float_t z2)
{
    const gmx_simd_float_t  FN6      = gmx_simd_set1_f(-1.7357322914161492954e-8f);
    const gmx_simd_float_t  FN5      = gmx_simd_set1_f(1.4703624142580877519e-6f);
    const gmx_simd_float_t  FN4      = gmx_simd_set1_f(-0.000053401640219807709149f);
    const gmx_simd_float_t  FN3      = gmx_simd_set1_f(0.0010054721316683106153f);
    const gmx_simd_float_t  FN2      = gmx_simd_set1_f(-0.019278317264888380590f);
    const gmx_simd_float_t  FN1      = gmx_simd_set1_f(0.069670166153766424023f);
    const gmx_simd_float_t  FN0      = gmx_simd_set1_f(-0.75225204789749321333f);

    const gmx_simd_float_t  FD4      = gmx_simd_set1_f(0.0011193462567257629232f);
    const gmx_simd_float_t  FD3      = gmx_simd_set1_f(0.014866955030185295499f);
    const gmx_simd_float_t  FD2      = gmx_simd_set1_f(0.11583842382862377919f);
    const gmx_simd_float_t  FD1      = gmx_simd_set1_f(0.50736591960530292870f);
    const gmx_simd_float_t  FD0      = gmx_simd_set1_f(1.0f);

    gmx_simd_float_t        z4;
    gmx_simd_float_t        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = gmx_simd_mul_f(z2, z2);

    polyFD0        = gmx_simd_fmadd_f(FD4, z4, FD2);
    polyFD1        = gmx_simd_fmadd_f(FD3, z4, FD1);
    polyFD0        = gmx_simd_fmadd_f(polyFD0, z4, FD0);
    polyFD0        = gmx_simd_fmadd_f(polyFD1, z2, polyFD0);

    polyFD0        = gmx_simd_inv_f(polyFD0);

    polyFN0        = gmx_simd_fmadd_f(FN6, z4, FN4);
    polyFN1        = gmx_simd_fmadd_f(FN5, z4, FN3);
    polyFN0        = gmx_simd_fmadd_f(polyFN0, z4, FN2);
    polyFN1        = gmx_simd_fmadd_f(polyFN1, z4, FN1);
    polyFN0        = gmx_simd_fmadd_f(polyFN0, z4, FN0);
    polyFN0        = gmx_simd_fmadd_f(polyFN1, z2, polyFN0);

    return gmx_simd_mul_f(polyFN0, polyFD0);
}



/*! \brief Calculate the potential correction due to PME analytically in SIMD float.
 *
 * You should normally call the real-precision routine \ref gmx_simd_pmecorrV_r.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb potential - see below for details.
 *
 * See \ref gmx_simd_pmecorrF_f for details about the approximation.
 *
 * This routine calculates \f$\mbox{erf}(z)/z\f$, although you should provide \f$z^2\f$
 * as the input argument.
 *
 * Here's how it should be used:
 *
 * 1. Calculate \f$r^2\f$.
 * 2. Multiply by \f$\beta^2\f$, so you get \f$z^2=\beta^2*r^2\f$.
 * 3. Evaluate this routine with z^2 as the argument.
 * 4. The return value is the expression:
 *
 *  \f[
 *   \frac{\mbox{erf}(z)}{z}
 *  \f]
 *
 * 5. Multiply the entire expression by beta and switching back to \f$r\f$ (since \f$z=r \beta\f$):
 *
 *  \f[
 *    \frac{\mbox{erf}(r \beta)}{r}
 *  \f]
 *
 * 6. Subtract the result from \f$1/r\f$, multiply by the product of the charges,
 *    and you have your potential.
 *
 * This approximation achieves an error slightly lower than 1e-6
 * in single precision and 4e-11 in double precision
 * for arguments smaller than 16 (\f$ 0.15 \leq \beta r \leq 4 \f$);
 * for \f$ \beta r \leq 0.15\f$ the error can be twice as high;
 * when added to \f$1/r\f$ the error will be insignificant.
 * For \f$\beta r \geq 7142\f$ the return value can be inf or NaN.
 */
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_pmecorrV_f(gmx_simd_float_t z2)
{
    const gmx_simd_float_t  VN6      = gmx_simd_set1_f(1.9296833005951166339e-8f);
    const gmx_simd_float_t  VN5      = gmx_simd_set1_f(-1.4213390571557850962e-6f);
    const gmx_simd_float_t  VN4      = gmx_simd_set1_f(0.000041603292906656984871f);
    const gmx_simd_float_t  VN3      = gmx_simd_set1_f(-0.00013134036773265025626f);
    const gmx_simd_float_t  VN2      = gmx_simd_set1_f(0.038657983986041781264f);
    const gmx_simd_float_t  VN1      = gmx_simd_set1_f(0.11285044772717598220f);
    const gmx_simd_float_t  VN0      = gmx_simd_set1_f(1.1283802385263030286f);

    const gmx_simd_float_t  VD3      = gmx_simd_set1_f(0.0066752224023576045451f);
    const gmx_simd_float_t  VD2      = gmx_simd_set1_f(0.078647795836373922256f);
    const gmx_simd_float_t  VD1      = gmx_simd_set1_f(0.43336185284710920150f);
    const gmx_simd_float_t  VD0      = gmx_simd_set1_f(1.0f);

    gmx_simd_float_t        z4;
    gmx_simd_float_t        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = gmx_simd_mul_f(z2, z2);

    polyVD1        = gmx_simd_fmadd_f(VD3, z4, VD1);
    polyVD0        = gmx_simd_fmadd_f(VD2, z4, VD0);
    polyVD0        = gmx_simd_fmadd_f(polyVD1, z2, polyVD0);

    polyVD0        = gmx_simd_inv_f(polyVD0);

    polyVN0        = gmx_simd_fmadd_f(VN6, z4, VN4);
    polyVN1        = gmx_simd_fmadd_f(VN5, z4, VN3);
    polyVN0        = gmx_simd_fmadd_f(polyVN0, z4, VN2);
    polyVN1        = gmx_simd_fmadd_f(polyVN1, z4, VN1);
    polyVN0        = gmx_simd_fmadd_f(polyVN0, z4, VN0);
    polyVN0        = gmx_simd_fmadd_f(polyVN1, z2, polyVN0);

    return gmx_simd_mul_f(polyVN0, polyVD0);
}
#endif

/*! \} */

#ifdef GMX_SIMD_HAVE_DOUBLE

/*! \name Double precision SIMD math functions
 *
 *  \note In most cases you should use the real-precision functions instead.
 *  \{
 */

/****************************************
 * DOUBLE PRECISION SIMD MATH FUNCTIONS *
 ****************************************/

/*! \brief SIMD utility function to sum a+b+c+d for SIMD doubles.
 *
 * \copydetails gmx_simd_sum4_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_sum4_d(gmx_simd_double_t a, gmx_simd_double_t b,
                gmx_simd_double_t c, gmx_simd_double_t d)
{
    return gmx_simd_add_d(gmx_simd_add_d(a, b), gmx_simd_add_d(c, d));
}

/*! \brief Return -a if b is negative, SIMD double.
 *
 * You should normally call the real-precision routine \ref gmx_simd_xor_sign_r.
 *
 * \param a Values to set sign for
 * \param b Values used to set sign
 * \return if b is negative, the sign of a will be changed.
 *
 * This is equivalent to doing an xor operation on a with the sign bit of b,
 * with the exception that negative zero is not considered to be negative
 * on architectures where \ref GMX_SIMD_HAVE_LOGICAL is not set.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_xor_sign_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
#ifdef GMX_SIMD_HAVE_LOGICAL
    return gmx_simd_xor_d(a, gmx_simd_and_d(gmx_simd_set1_d(GMX_DOUBLE_NEGZERO), b));
#else
    return gmx_simd_blendv_d(a, gmx_simd_fneg_d(a), gmx_simd_cmplt_d(b, gmx_simd_setzero_d()));
#endif
}

#ifndef gmx_simd_rsqrt_iter_d
/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD double.
 *
 * \copydetails gmx_simd_rsqrt_iter_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_rsqrt_iter_d(gmx_simd_double_t lu, gmx_simd_double_t x)
{
#ifdef GMX_SIMD_HAVE_FMA
    return gmx_simd_fmadd_d(gmx_simd_fnmadd_d(x, gmx_simd_mul_d(lu, lu), gmx_simd_set1_d(1.0)), gmx_simd_mul_d(lu, gmx_simd_set1_d(0.5)), lu);
#else
    return gmx_simd_mul_d(gmx_simd_set1_d(0.5), gmx_simd_mul_d(gmx_simd_sub_d(gmx_simd_set1_d(3.0), gmx_simd_mul_d(gmx_simd_mul_d(lu, lu), x)), lu));
#endif
}
#endif

/*! \brief Calculate 1/sqrt(x) for SIMD double
 *
 * \copydetails gmx_simd_invsqrt_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_invsqrt_d(gmx_simd_double_t x)
{
    gmx_simd_double_t lu = gmx_simd_rsqrt_d(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd_rsqrt_iter_d(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for masked entries of SIMD double.
 *
 * \copydetails gmx_simd_invsqrt_maskfpe_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_invsqrt_maskfpe_d(gmx_simd_double_t x, gmx_simd_dbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_invsqrt_d(x);
#else
    return gmx_simd_invsqrt_d(gmx_simd_blendv_d(gmx_simd_set1_d(1.0), x, m));
#endif
}

/*! \brief Calculate 1/sqrt(x) for non-masked entries of SIMD double.
 *
 * \copydetails gmx_simd_invsqrt_notmaskfpe_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_invsqrt_notmaskfpe_d(gmx_simd_double_t x, gmx_simd_dbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_invsqrt_d(x);
#else
    return gmx_simd_invsqrt_d(gmx_simd_blendv_d(x, gmx_simd_set1_d(1.0), m));
#endif
}

/*! \brief Calculate 1/sqrt(x) for two SIMD doubles.
 *
 * \copydetails gmx_simd_invsqrt_pair_f
 */
static gmx_inline void gmx_simdcall
gmx_simd_invsqrt_pair_d(gmx_simd_double_t x0,    gmx_simd_double_t x1,
                        gmx_simd_double_t *out0, gmx_simd_double_t *out1)
{
#if (defined GMX_SIMD_HAVE_FLOAT) && (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH) && (GMX_SIMD_RSQRT_BITS < 22)
    gmx_simd_float_t  xf  = gmx_simd_cvt_dd2f(x0, x1);
    gmx_simd_float_t  luf = gmx_simd_rsqrt_f(xf);
    gmx_simd_double_t lu0, lu1;
    /* Intermediate target is single - mantissa+1 bits */
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = gmx_simd_rsqrt_iter_f(luf, xf);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = gmx_simd_rsqrt_iter_f(luf, xf);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = gmx_simd_rsqrt_iter_f(luf, xf);
#endif
    gmx_simd_cvt_f2dd(luf, &lu0, &lu1);
    /* Last iteration(s) performed in double - if we had 22 bits, this gets us to 44 (~1e-15) */
#if (GMX_SIMD_ACCURACY_BITS_SINGLE < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu0 = gmx_simd_rsqrt_iter_d(lu0, x0);
    lu1 = gmx_simd_rsqrt_iter_d(lu1, x1);
#endif
#if (GMX_SIMD_ACCURACY_BITS_SINGLE*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu0 = gmx_simd_rsqrt_iter_d(lu0, x0);
    lu1 = gmx_simd_rsqrt_iter_d(lu1, x1);
#endif
    *out0 = lu0;
    *out1 = lu1;
#else
    *out0 = gmx_simd_invsqrt_d(x0);
    *out1 = gmx_simd_invsqrt_d(x1);
#endif
}

#ifndef gmx_simd_rcp_iter_d
/*! \brief Perform one Newton-Raphson iteration to improve 1/x for SIMD double.
 *
 * \copydetails gmx_simd_rcp_iter_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_rcp_iter_d(gmx_simd_double_t lu, gmx_simd_double_t x)
{
    return gmx_simd_mul_d(lu, gmx_simd_fnmadd_d(lu, x, gmx_simd_set1_d(2.0)));
}
#endif

/*! \brief Calculate 1/x for SIMD double.
 *
 * \copydetails gmx_simd_inv_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_inv_d(gmx_simd_double_t x)
{
    gmx_simd_double_t lu = gmx_simd_rcp_d(x);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd_rcp_iter_d(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd_rcp_iter_d(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd_rcp_iter_d(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd_rcp_iter_d(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/x for masked entries of SIMD double.
 *
 * \copydetails gmx_simd_inv_maskfpe_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_inv_maskfpe_d(gmx_simd_double_t x, gmx_simd_dbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_inv_d(x);
#else
    return gmx_simd_inv_d(gmx_simd_blendv_d(gmx_simd_set1_d(1.0), x, m));
#endif
}

/*! \brief Calculate 1/x for non-masked entries of SIMD double.
 *
 * \copydetails gmx_simd_inv_notmaskfpe_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_inv_notmaskfpe_d(gmx_simd_double_t x, gmx_simd_dbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_inv_d(x);
#else
    return gmx_simd_inv_d(gmx_simd_blendv_d(x, gmx_simd_set1_d(1.0), m));
#endif
}

/*! \brief Calculate sqrt(x) correctly for SIMD doubles, including argument 0.0.
 *
 * \copydetails gmx_simd_sqrt_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_sqrt_d(gmx_simd_double_t x)
{
    gmx_simd_dbool_t   mask;
    gmx_simd_double_t  res;

    mask = gmx_simd_cmpeq_d(x, gmx_simd_setzero_d());
    res  = gmx_simd_blendnotzero_d(gmx_simd_invsqrt_notmaskfpe_d(x, mask), mask);
    return gmx_simd_mul_d(res, x);
}

/*! \brief SIMD double log(x). This is the natural logarithm.
 *
 * \copydetails gmx_simd_log_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_log_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t  half       = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t  one        = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t  sqrt2      = gmx_simd_set1_d(sqrt(2.0));
    const gmx_simd_double_t  corr       = gmx_simd_set1_d(0.693147180559945286226764);
    const gmx_simd_double_t  CL15       = gmx_simd_set1_d(0.148197055177935105296783);
    const gmx_simd_double_t  CL13       = gmx_simd_set1_d(0.153108178020442575739679);
    const gmx_simd_double_t  CL11       = gmx_simd_set1_d(0.181837339521549679055568);
    const gmx_simd_double_t  CL9        = gmx_simd_set1_d(0.22222194152736701733275);
    const gmx_simd_double_t  CL7        = gmx_simd_set1_d(0.285714288030134544449368);
    const gmx_simd_double_t  CL5        = gmx_simd_set1_d(0.399999999989941956712869);
    const gmx_simd_double_t  CL3        = gmx_simd_set1_d(0.666666666666685503450651);
    const gmx_simd_double_t  CL1        = gmx_simd_set1_d(2.0);
    gmx_simd_double_t        fexp, x2, p;
    gmx_simd_dbool_t         mask;

    fexp  = gmx_simd_get_exponent_d(x);
    x     = gmx_simd_get_mantissa_d(x);

    mask  = gmx_simd_cmplt_d(sqrt2, x);
    /* Adjust to non-IEEE format for x>sqrt(2): exponent += 1, mantissa *= 0.5 */
    fexp  = gmx_simd_add_d(fexp, gmx_simd_blendzero_d(one, mask));
    x     = gmx_simd_mul_d(x, gmx_simd_blendv_d(one, half, mask));

    x     = gmx_simd_mul_d( gmx_simd_sub_d(x, one), gmx_simd_inv_d( gmx_simd_add_d(x, one) ) );
    x2    = gmx_simd_mul_d(x, x);

    p     = gmx_simd_fmadd_d(CL15, x2, CL13);
    p     = gmx_simd_fmadd_d(p, x2, CL11);
    p     = gmx_simd_fmadd_d(p, x2, CL9);
    p     = gmx_simd_fmadd_d(p, x2, CL7);
    p     = gmx_simd_fmadd_d(p, x2, CL5);
    p     = gmx_simd_fmadd_d(p, x2, CL3);
    p     = gmx_simd_fmadd_d(p, x2, CL1);
    p     = gmx_simd_fmadd_d(p, x, gmx_simd_mul_d(corr, fexp));

    return p;
}

/*! \brief SIMD double 2^x.
 *
 * \copydetails gmx_simd_exp2_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_exp2_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t  arglimit      = gmx_simd_set1_d(1022.0);
    const gmx_simd_double_t  CE11          = gmx_simd_set1_d(4.435280790452730022081181e-10);
    const gmx_simd_double_t  CE10          = gmx_simd_set1_d(7.074105630863314448024247e-09);
    const gmx_simd_double_t  CE9           = gmx_simd_set1_d(1.017819803432096698472621e-07);
    const gmx_simd_double_t  CE8           = gmx_simd_set1_d(1.321543308956718799557863e-06);
    const gmx_simd_double_t  CE7           = gmx_simd_set1_d(0.00001525273348995851746990884);
    const gmx_simd_double_t  CE6           = gmx_simd_set1_d(0.0001540353046251466849082632);
    const gmx_simd_double_t  CE5           = gmx_simd_set1_d(0.001333355814678995257307880);
    const gmx_simd_double_t  CE4           = gmx_simd_set1_d(0.009618129107588335039176502);
    const gmx_simd_double_t  CE3           = gmx_simd_set1_d(0.05550410866481992147457793);
    const gmx_simd_double_t  CE2           = gmx_simd_set1_d(0.2402265069591015620470894);
    const gmx_simd_double_t  CE1           = gmx_simd_set1_d(0.6931471805599453304615075);
    const gmx_simd_double_t  one           = gmx_simd_set1_d(1.0);
    gmx_simd_double_t        fexppart;
    gmx_simd_double_t        intpart;
    gmx_simd_double_t        p;
    gmx_simd_dbool_t         valuemask;

    fexppart  = gmx_simd_set_exponent_d(x);  /* rounds to nearest int internally */
    intpart   = gmx_simd_round_d(x);         /* use same rounding mode here */
    valuemask = gmx_simd_cmple_d(gmx_simd_fabs_d(x), arglimit);
    fexppart  = gmx_simd_blendzero_d(fexppart, valuemask);
    x         = gmx_simd_sub_d(x, intpart);

    p         = gmx_simd_fmadd_d(CE11, x, CE10);
    p         = gmx_simd_fmadd_d(p, x, CE9);
    p         = gmx_simd_fmadd_d(p, x, CE8);
    p         = gmx_simd_fmadd_d(p, x, CE7);
    p         = gmx_simd_fmadd_d(p, x, CE6);
    p         = gmx_simd_fmadd_d(p, x, CE5);
    p         = gmx_simd_fmadd_d(p, x, CE4);
    p         = gmx_simd_fmadd_d(p, x, CE3);
    p         = gmx_simd_fmadd_d(p, x, CE2);
    p         = gmx_simd_fmadd_d(p, x, CE1);
    p         = gmx_simd_fmadd_d(p, x, one);
    x         = gmx_simd_mul_d(p, fexppart);
    return x;
}

/*! \brief SIMD double exp(x).
 *
 * \copydetails gmx_simd_exp_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_exp_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t  argscale      = gmx_simd_set1_d(1.44269504088896340735992468100);
    const gmx_simd_double_t  arglimit      = gmx_simd_set1_d(1022.0);
    const gmx_simd_double_t  invargscale0  = gmx_simd_set1_d(-0.69314718055966295651160180568695068359375);
    const gmx_simd_double_t  invargscale1  = gmx_simd_set1_d(-2.8235290563031577122588448175013436025525412068e-13);
    const gmx_simd_double_t  CE12          = gmx_simd_set1_d(2.078375306791423699350304e-09);
    const gmx_simd_double_t  CE11          = gmx_simd_set1_d(2.518173854179933105218635e-08);
    const gmx_simd_double_t  CE10          = gmx_simd_set1_d(2.755842049600488770111608e-07);
    const gmx_simd_double_t  CE9           = gmx_simd_set1_d(2.755691815216689746619849e-06);
    const gmx_simd_double_t  CE8           = gmx_simd_set1_d(2.480158383706245033920920e-05);
    const gmx_simd_double_t  CE7           = gmx_simd_set1_d(0.0001984127043518048611841321);
    const gmx_simd_double_t  CE6           = gmx_simd_set1_d(0.001388888889360258341755930);
    const gmx_simd_double_t  CE5           = gmx_simd_set1_d(0.008333333332907368102819109);
    const gmx_simd_double_t  CE4           = gmx_simd_set1_d(0.04166666666663836745814631);
    const gmx_simd_double_t  CE3           = gmx_simd_set1_d(0.1666666666666796929434570);
    const gmx_simd_double_t  CE2           = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t  one           = gmx_simd_set1_d(1.0);
    gmx_simd_double_t        fexppart;
    gmx_simd_double_t        intpart;
    gmx_simd_double_t        y, p;
    gmx_simd_dbool_t         valuemask;

    y         = gmx_simd_mul_d(x, argscale);
    fexppart  = gmx_simd_set_exponent_d(y);  /* rounds to nearest int internally */
    intpart   = gmx_simd_round_d(y);         /* use same rounding mode here */
    valuemask = gmx_simd_cmple_d(gmx_simd_fabs_d(y), arglimit);
    fexppart  = gmx_simd_blendzero_d(fexppart, valuemask);

    /* Extended precision arithmetics */
    x         = gmx_simd_fmadd_d(invargscale0, intpart, x);
    x         = gmx_simd_fmadd_d(invargscale1, intpart, x);

    p         = gmx_simd_fmadd_d(CE12, x, CE11);
    p         = gmx_simd_fmadd_d(p, x, CE10);
    p         = gmx_simd_fmadd_d(p, x, CE9);
    p         = gmx_simd_fmadd_d(p, x, CE8);
    p         = gmx_simd_fmadd_d(p, x, CE7);
    p         = gmx_simd_fmadd_d(p, x, CE6);
    p         = gmx_simd_fmadd_d(p, x, CE5);
    p         = gmx_simd_fmadd_d(p, x, CE4);
    p         = gmx_simd_fmadd_d(p, x, CE3);
    p         = gmx_simd_fmadd_d(p, x, CE2);
    p         = gmx_simd_fmadd_d(p, gmx_simd_mul_d(x, x), gmx_simd_add_d(x, one));
    x         = gmx_simd_mul_d(p, fexppart);
    return x;
}

/*! \brief SIMD double erf(x).
 *
 * \copydetails gmx_simd_erf_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_erf_d(gmx_simd_double_t x)
{
    /* Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75] */
    const gmx_simd_double_t CAP4      = gmx_simd_set1_d(-0.431780540597889301512e-4);
    const gmx_simd_double_t CAP3      = gmx_simd_set1_d(-0.00578562306260059236059);
    const gmx_simd_double_t CAP2      = gmx_simd_set1_d(-0.028593586920219752446);
    const gmx_simd_double_t CAP1      = gmx_simd_set1_d(-0.315924962948621698209);
    const gmx_simd_double_t CAP0      = gmx_simd_set1_d(0.14952975608477029151);

    const gmx_simd_double_t CAQ5      = gmx_simd_set1_d(-0.374089300177174709737e-5);
    const gmx_simd_double_t CAQ4      = gmx_simd_set1_d(0.00015126584532155383535);
    const gmx_simd_double_t CAQ3      = gmx_simd_set1_d(0.00536692680669480725423);
    const gmx_simd_double_t CAQ2      = gmx_simd_set1_d(0.0668686825594046122636);
    const gmx_simd_double_t CAQ1      = gmx_simd_set1_d(0.402604990869284362773);
    /* CAQ0 == 1.0 */
    const gmx_simd_double_t CAoffset  = gmx_simd_set1_d(0.9788494110107421875);

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5] */
    const gmx_simd_double_t CBP6      = gmx_simd_set1_d(2.49650423685462752497647637088e-10);
    const gmx_simd_double_t CBP5      = gmx_simd_set1_d(0.00119770193298159629350136085658);
    const gmx_simd_double_t CBP4      = gmx_simd_set1_d(0.0164944422378370965881008942733);
    const gmx_simd_double_t CBP3      = gmx_simd_set1_d(0.0984581468691775932063932439252);
    const gmx_simd_double_t CBP2      = gmx_simd_set1_d(0.317364595806937763843589437418);
    const gmx_simd_double_t CBP1      = gmx_simd_set1_d(0.554167062641455850932670067075);
    const gmx_simd_double_t CBP0      = gmx_simd_set1_d(0.427583576155807163756925301060);
    const gmx_simd_double_t CBQ7      = gmx_simd_set1_d(0.00212288829699830145976198384930);
    const gmx_simd_double_t CBQ6      = gmx_simd_set1_d(0.0334810979522685300554606393425);
    const gmx_simd_double_t CBQ5      = gmx_simd_set1_d(0.2361713785181450957579508850717);
    const gmx_simd_double_t CBQ4      = gmx_simd_set1_d(0.955364736493055670530981883072);
    const gmx_simd_double_t CBQ3      = gmx_simd_set1_d(2.36815675631420037315349279199);
    const gmx_simd_double_t CBQ2      = gmx_simd_set1_d(3.55261649184083035537184223542);
    const gmx_simd_double_t CBQ1      = gmx_simd_set1_d(2.93501136050160872574376997993);
    /* CBQ0 == 1.0 */

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf] */
    const gmx_simd_double_t CCP6      = gmx_simd_set1_d(-2.8175401114513378771);
    const gmx_simd_double_t CCP5      = gmx_simd_set1_d(-3.22729451764143718517);
    const gmx_simd_double_t CCP4      = gmx_simd_set1_d(-2.5518551727311523996);
    const gmx_simd_double_t CCP3      = gmx_simd_set1_d(-0.687717681153649930619);
    const gmx_simd_double_t CCP2      = gmx_simd_set1_d(-0.212652252872804219852);
    const gmx_simd_double_t CCP1      = gmx_simd_set1_d(0.0175389834052493308818);
    const gmx_simd_double_t CCP0      = gmx_simd_set1_d(0.00628057170626964891937);

    const gmx_simd_double_t CCQ6      = gmx_simd_set1_d(5.48409182238641741584);
    const gmx_simd_double_t CCQ5      = gmx_simd_set1_d(13.5064170191802889145);
    const gmx_simd_double_t CCQ4      = gmx_simd_set1_d(22.9367376522880577224);
    const gmx_simd_double_t CCQ3      = gmx_simd_set1_d(15.930646027911794143);
    const gmx_simd_double_t CCQ2      = gmx_simd_set1_d(11.0567237927800161565);
    const gmx_simd_double_t CCQ1      = gmx_simd_set1_d(2.79257750980575282228);
    /* CCQ0 == 1.0 */
    const gmx_simd_double_t CCoffset  = gmx_simd_set1_d(0.5579090118408203125);

    const gmx_simd_double_t one       = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t two       = gmx_simd_set1_d(2.0);

    gmx_simd_double_t       xabs, x2, x4, t, t2, w, w2;
    gmx_simd_double_t       PolyAP0, PolyAP1, PolyAQ0, PolyAQ1;
    gmx_simd_double_t       PolyBP0, PolyBP1, PolyBQ0, PolyBQ1;
    gmx_simd_double_t       PolyCP0, PolyCP1, PolyCQ0, PolyCQ1;
    gmx_simd_double_t       res_erf, res_erfcB, res_erfcC, res_erfc, res;
    gmx_simd_double_t       expmx2;
    gmx_simd_dbool_t        mask, mask_erf;

    /* Calculate erf() */
    xabs     = gmx_simd_fabs_d(x);
    mask_erf = gmx_simd_cmplt_d(xabs, one);
    x2       = gmx_simd_mul_d(x, x);
    x4       = gmx_simd_mul_d(x2, x2);

    PolyAP0  = gmx_simd_mul_d(CAP4, x4);
    PolyAP1  = gmx_simd_mul_d(CAP3, x4);
    PolyAP0  = gmx_simd_add_d(PolyAP0, CAP2);
    PolyAP1  = gmx_simd_add_d(PolyAP1, CAP1);
    PolyAP0  = gmx_simd_mul_d(PolyAP0, x4);
    PolyAP1  = gmx_simd_mul_d(PolyAP1, x2);
    PolyAP0  = gmx_simd_add_d(PolyAP0, CAP0);
    PolyAP0  = gmx_simd_add_d(PolyAP0, PolyAP1);

    PolyAQ1  = gmx_simd_mul_d(CAQ5, x4);
    PolyAQ0  = gmx_simd_mul_d(CAQ4, x4);
    PolyAQ1  = gmx_simd_add_d(PolyAQ1, CAQ3);
    PolyAQ0  = gmx_simd_add_d(PolyAQ0, CAQ2);
    PolyAQ1  = gmx_simd_mul_d(PolyAQ1, x4);
    PolyAQ0  = gmx_simd_mul_d(PolyAQ0, x4);
    PolyAQ1  = gmx_simd_add_d(PolyAQ1, CAQ1);
    PolyAQ0  = gmx_simd_add_d(PolyAQ0, one);
    PolyAQ1  = gmx_simd_mul_d(PolyAQ1, x2);
    PolyAQ0  = gmx_simd_add_d(PolyAQ0, PolyAQ1);

    res_erf  = gmx_simd_mul_d(PolyAP0, gmx_simd_inv_maskfpe_d(PolyAQ0, mask_erf));
    res_erf  = gmx_simd_add_d(CAoffset, res_erf);
    res_erf  = gmx_simd_mul_d(x, res_erf);

    /* Calculate erfc() in range [1,4.5] */
    t       = gmx_simd_sub_d(xabs, one);
    t2      = gmx_simd_mul_d(t, t);

    PolyBP0  = gmx_simd_mul_d(CBP6, t2);
    PolyBP1  = gmx_simd_mul_d(CBP5, t2);
    PolyBP0  = gmx_simd_add_d(PolyBP0, CBP4);
    PolyBP1  = gmx_simd_add_d(PolyBP1, CBP3);
    PolyBP0  = gmx_simd_mul_d(PolyBP0, t2);
    PolyBP1  = gmx_simd_mul_d(PolyBP1, t2);
    PolyBP0  = gmx_simd_add_d(PolyBP0, CBP2);
    PolyBP1  = gmx_simd_add_d(PolyBP1, CBP1);
    PolyBP0  = gmx_simd_mul_d(PolyBP0, t2);
    PolyBP1  = gmx_simd_mul_d(PolyBP1, t);
    PolyBP0  = gmx_simd_add_d(PolyBP0, CBP0);
    PolyBP0  = gmx_simd_add_d(PolyBP0, PolyBP1);

    PolyBQ1 = gmx_simd_mul_d(CBQ7, t2);
    PolyBQ0 = gmx_simd_mul_d(CBQ6, t2);
    PolyBQ1 = gmx_simd_add_d(PolyBQ1, CBQ5);
    PolyBQ0 = gmx_simd_add_d(PolyBQ0, CBQ4);
    PolyBQ1 = gmx_simd_mul_d(PolyBQ1, t2);
    PolyBQ0 = gmx_simd_mul_d(PolyBQ0, t2);
    PolyBQ1 = gmx_simd_add_d(PolyBQ1, CBQ3);
    PolyBQ0 = gmx_simd_add_d(PolyBQ0, CBQ2);
    PolyBQ1 = gmx_simd_mul_d(PolyBQ1, t2);
    PolyBQ0 = gmx_simd_mul_d(PolyBQ0, t2);
    PolyBQ1 = gmx_simd_add_d(PolyBQ1, CBQ1);
    PolyBQ0 = gmx_simd_add_d(PolyBQ0, one);
    PolyBQ1 = gmx_simd_mul_d(PolyBQ1, t);
    PolyBQ0 = gmx_simd_add_d(PolyBQ0, PolyBQ1);

    res_erfcB = gmx_simd_mul_d(PolyBP0, gmx_simd_inv_notmaskfpe_d(PolyBQ0, mask_erf));

    res_erfcB = gmx_simd_mul_d(res_erfcB, xabs);

    /* Calculate erfc() in range [4.5,inf] */
    w       = gmx_simd_inv_notmaskfpe_d(xabs, mask_erf);
    w2      = gmx_simd_mul_d(w, w);

    PolyCP0  = gmx_simd_mul_d(CCP6, w2);
    PolyCP1  = gmx_simd_mul_d(CCP5, w2);
    PolyCP0  = gmx_simd_add_d(PolyCP0, CCP4);
    PolyCP1  = gmx_simd_add_d(PolyCP1, CCP3);
    PolyCP0  = gmx_simd_mul_d(PolyCP0, w2);
    PolyCP1  = gmx_simd_mul_d(PolyCP1, w2);
    PolyCP0  = gmx_simd_add_d(PolyCP0, CCP2);
    PolyCP1  = gmx_simd_add_d(PolyCP1, CCP1);
    PolyCP0  = gmx_simd_mul_d(PolyCP0, w2);
    PolyCP1  = gmx_simd_mul_d(PolyCP1, w);
    PolyCP0  = gmx_simd_add_d(PolyCP0, CCP0);
    PolyCP0  = gmx_simd_add_d(PolyCP0, PolyCP1);

    PolyCQ0  = gmx_simd_mul_d(CCQ6, w2);
    PolyCQ1  = gmx_simd_mul_d(CCQ5, w2);
    PolyCQ0  = gmx_simd_add_d(PolyCQ0, CCQ4);
    PolyCQ1  = gmx_simd_add_d(PolyCQ1, CCQ3);
    PolyCQ0  = gmx_simd_mul_d(PolyCQ0, w2);
    PolyCQ1  = gmx_simd_mul_d(PolyCQ1, w2);
    PolyCQ0  = gmx_simd_add_d(PolyCQ0, CCQ2);
    PolyCQ1  = gmx_simd_add_d(PolyCQ1, CCQ1);
    PolyCQ0  = gmx_simd_mul_d(PolyCQ0, w2);
    PolyCQ1  = gmx_simd_mul_d(PolyCQ1, w);
    PolyCQ0  = gmx_simd_add_d(PolyCQ0, one);
    PolyCQ0  = gmx_simd_add_d(PolyCQ0, PolyCQ1);

    expmx2   = gmx_simd_exp_d( gmx_simd_fneg_d(x2) );

    res_erfcC = gmx_simd_mul_d(PolyCP0, gmx_simd_inv_notmaskfpe_d(PolyCQ0, mask_erf));
    res_erfcC = gmx_simd_add_d(res_erfcC, CCoffset);
    res_erfcC = gmx_simd_mul_d(res_erfcC, w);

    mask     = gmx_simd_cmplt_d(gmx_simd_set1_d(4.5), xabs);
    res_erfc = gmx_simd_blendv_d(res_erfcB, res_erfcC, mask);

    res_erfc = gmx_simd_mul_d(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = gmx_simd_cmplt_d(x, gmx_simd_setzero_d());
    res_erfc = gmx_simd_blendv_d(res_erfc, gmx_simd_sub_d(two, res_erfc), mask);

    /* Select erf() or erfc() */
    res  = gmx_simd_blendv_d(gmx_simd_sub_d(one, res_erfc), res_erf, mask_erf);

    return res;
}

/*! \brief SIMD double erfc(x).
 *
 * \copydetails gmx_simd_erfc_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_erfc_d(gmx_simd_double_t x)
{
    /* Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75] */
    const gmx_simd_double_t CAP4      = gmx_simd_set1_d(-0.431780540597889301512e-4);
    const gmx_simd_double_t CAP3      = gmx_simd_set1_d(-0.00578562306260059236059);
    const gmx_simd_double_t CAP2      = gmx_simd_set1_d(-0.028593586920219752446);
    const gmx_simd_double_t CAP1      = gmx_simd_set1_d(-0.315924962948621698209);
    const gmx_simd_double_t CAP0      = gmx_simd_set1_d(0.14952975608477029151);

    const gmx_simd_double_t CAQ5      = gmx_simd_set1_d(-0.374089300177174709737e-5);
    const gmx_simd_double_t CAQ4      = gmx_simd_set1_d(0.00015126584532155383535);
    const gmx_simd_double_t CAQ3      = gmx_simd_set1_d(0.00536692680669480725423);
    const gmx_simd_double_t CAQ2      = gmx_simd_set1_d(0.0668686825594046122636);
    const gmx_simd_double_t CAQ1      = gmx_simd_set1_d(0.402604990869284362773);
    /* CAQ0 == 1.0 */
    const gmx_simd_double_t CAoffset  = gmx_simd_set1_d(0.9788494110107421875);

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5] */
    const gmx_simd_double_t CBP6      = gmx_simd_set1_d(2.49650423685462752497647637088e-10);
    const gmx_simd_double_t CBP5      = gmx_simd_set1_d(0.00119770193298159629350136085658);
    const gmx_simd_double_t CBP4      = gmx_simd_set1_d(0.0164944422378370965881008942733);
    const gmx_simd_double_t CBP3      = gmx_simd_set1_d(0.0984581468691775932063932439252);
    const gmx_simd_double_t CBP2      = gmx_simd_set1_d(0.317364595806937763843589437418);
    const gmx_simd_double_t CBP1      = gmx_simd_set1_d(0.554167062641455850932670067075);
    const gmx_simd_double_t CBP0      = gmx_simd_set1_d(0.427583576155807163756925301060);
    const gmx_simd_double_t CBQ7      = gmx_simd_set1_d(0.00212288829699830145976198384930);
    const gmx_simd_double_t CBQ6      = gmx_simd_set1_d(0.0334810979522685300554606393425);
    const gmx_simd_double_t CBQ5      = gmx_simd_set1_d(0.2361713785181450957579508850717);
    const gmx_simd_double_t CBQ4      = gmx_simd_set1_d(0.955364736493055670530981883072);
    const gmx_simd_double_t CBQ3      = gmx_simd_set1_d(2.36815675631420037315349279199);
    const gmx_simd_double_t CBQ2      = gmx_simd_set1_d(3.55261649184083035537184223542);
    const gmx_simd_double_t CBQ1      = gmx_simd_set1_d(2.93501136050160872574376997993);
    /* CBQ0 == 1.0 */

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf] */
    const gmx_simd_double_t CCP6      = gmx_simd_set1_d(-2.8175401114513378771);
    const gmx_simd_double_t CCP5      = gmx_simd_set1_d(-3.22729451764143718517);
    const gmx_simd_double_t CCP4      = gmx_simd_set1_d(-2.5518551727311523996);
    const gmx_simd_double_t CCP3      = gmx_simd_set1_d(-0.687717681153649930619);
    const gmx_simd_double_t CCP2      = gmx_simd_set1_d(-0.212652252872804219852);
    const gmx_simd_double_t CCP1      = gmx_simd_set1_d(0.0175389834052493308818);
    const gmx_simd_double_t CCP0      = gmx_simd_set1_d(0.00628057170626964891937);

    const gmx_simd_double_t CCQ6      = gmx_simd_set1_d(5.48409182238641741584);
    const gmx_simd_double_t CCQ5      = gmx_simd_set1_d(13.5064170191802889145);
    const gmx_simd_double_t CCQ4      = gmx_simd_set1_d(22.9367376522880577224);
    const gmx_simd_double_t CCQ3      = gmx_simd_set1_d(15.930646027911794143);
    const gmx_simd_double_t CCQ2      = gmx_simd_set1_d(11.0567237927800161565);
    const gmx_simd_double_t CCQ1      = gmx_simd_set1_d(2.79257750980575282228);
    /* CCQ0 == 1.0 */
    const gmx_simd_double_t CCoffset  = gmx_simd_set1_d(0.5579090118408203125);

    const gmx_simd_double_t one       = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t two       = gmx_simd_set1_d(2.0);

    gmx_simd_double_t       xabs, x2, x4, t, t2, w, w2;
    gmx_simd_double_t       PolyAP0, PolyAP1, PolyAQ0, PolyAQ1;
    gmx_simd_double_t       PolyBP0, PolyBP1, PolyBQ0, PolyBQ1;
    gmx_simd_double_t       PolyCP0, PolyCP1, PolyCQ0, PolyCQ1;
    gmx_simd_double_t       res_erf, res_erfcB, res_erfcC, res_erfc, res;
    gmx_simd_double_t       expmx2;
    gmx_simd_dbool_t        mask, mask_erf;

    /* Calculate erf() */
    xabs     = gmx_simd_fabs_d(x);
    mask_erf = gmx_simd_cmplt_d(xabs, one);
    x2       = gmx_simd_mul_d(x, x);
    x4       = gmx_simd_mul_d(x2, x2);

    PolyAP0  = gmx_simd_mul_d(CAP4, x4);
    PolyAP1  = gmx_simd_mul_d(CAP3, x4);
    PolyAP0  = gmx_simd_add_d(PolyAP0, CAP2);
    PolyAP1  = gmx_simd_add_d(PolyAP1, CAP1);
    PolyAP0  = gmx_simd_mul_d(PolyAP0, x4);
    PolyAP1  = gmx_simd_mul_d(PolyAP1, x2);
    PolyAP0  = gmx_simd_add_d(PolyAP0, CAP0);
    PolyAP0  = gmx_simd_add_d(PolyAP0, PolyAP1);

    PolyAQ1  = gmx_simd_mul_d(CAQ5, x4);
    PolyAQ0  = gmx_simd_mul_d(CAQ4, x4);
    PolyAQ1  = gmx_simd_add_d(PolyAQ1, CAQ3);
    PolyAQ0  = gmx_simd_add_d(PolyAQ0, CAQ2);
    PolyAQ1  = gmx_simd_mul_d(PolyAQ1, x4);
    PolyAQ0  = gmx_simd_mul_d(PolyAQ0, x4);
    PolyAQ1  = gmx_simd_add_d(PolyAQ1, CAQ1);
    PolyAQ0  = gmx_simd_add_d(PolyAQ0, one);
    PolyAQ1  = gmx_simd_mul_d(PolyAQ1, x2);
    PolyAQ0  = gmx_simd_add_d(PolyAQ0, PolyAQ1);

    res_erf  = gmx_simd_mul_d(PolyAP0, gmx_simd_inv_maskfpe_d(PolyAQ0, mask_erf));
    res_erf  = gmx_simd_add_d(CAoffset, res_erf);
    res_erf  = gmx_simd_mul_d(x, res_erf);

    /* Calculate erfc() in range [1,4.5] */
    t       = gmx_simd_sub_d(xabs, one);
    t2      = gmx_simd_mul_d(t, t);

    PolyBP0  = gmx_simd_mul_d(CBP6, t2);
    PolyBP1  = gmx_simd_mul_d(CBP5, t2);
    PolyBP0  = gmx_simd_add_d(PolyBP0, CBP4);
    PolyBP1  = gmx_simd_add_d(PolyBP1, CBP3);
    PolyBP0  = gmx_simd_mul_d(PolyBP0, t2);
    PolyBP1  = gmx_simd_mul_d(PolyBP1, t2);
    PolyBP0  = gmx_simd_add_d(PolyBP0, CBP2);
    PolyBP1  = gmx_simd_add_d(PolyBP1, CBP1);
    PolyBP0  = gmx_simd_mul_d(PolyBP0, t2);
    PolyBP1  = gmx_simd_mul_d(PolyBP1, t);
    PolyBP0  = gmx_simd_add_d(PolyBP0, CBP0);
    PolyBP0  = gmx_simd_add_d(PolyBP0, PolyBP1);

    PolyBQ1 = gmx_simd_mul_d(CBQ7, t2);
    PolyBQ0 = gmx_simd_mul_d(CBQ6, t2);
    PolyBQ1 = gmx_simd_add_d(PolyBQ1, CBQ5);
    PolyBQ0 = gmx_simd_add_d(PolyBQ0, CBQ4);
    PolyBQ1 = gmx_simd_mul_d(PolyBQ1, t2);
    PolyBQ0 = gmx_simd_mul_d(PolyBQ0, t2);
    PolyBQ1 = gmx_simd_add_d(PolyBQ1, CBQ3);
    PolyBQ0 = gmx_simd_add_d(PolyBQ0, CBQ2);
    PolyBQ1 = gmx_simd_mul_d(PolyBQ1, t2);
    PolyBQ0 = gmx_simd_mul_d(PolyBQ0, t2);
    PolyBQ1 = gmx_simd_add_d(PolyBQ1, CBQ1);
    PolyBQ0 = gmx_simd_add_d(PolyBQ0, one);
    PolyBQ1 = gmx_simd_mul_d(PolyBQ1, t);
    PolyBQ0 = gmx_simd_add_d(PolyBQ0, PolyBQ1);

    res_erfcB = gmx_simd_mul_d(PolyBP0, gmx_simd_inv_notmaskfpe_d(PolyBQ0, mask_erf));

    res_erfcB = gmx_simd_mul_d(res_erfcB, xabs);

    /* Calculate erfc() in range [4.5,inf] */
    w       = gmx_simd_inv_notmaskfpe_d(xabs, mask_erf);
    w2      = gmx_simd_mul_d(w, w);

    PolyCP0  = gmx_simd_mul_d(CCP6, w2);
    PolyCP1  = gmx_simd_mul_d(CCP5, w2);
    PolyCP0  = gmx_simd_add_d(PolyCP0, CCP4);
    PolyCP1  = gmx_simd_add_d(PolyCP1, CCP3);
    PolyCP0  = gmx_simd_mul_d(PolyCP0, w2);
    PolyCP1  = gmx_simd_mul_d(PolyCP1, w2);
    PolyCP0  = gmx_simd_add_d(PolyCP0, CCP2);
    PolyCP1  = gmx_simd_add_d(PolyCP1, CCP1);
    PolyCP0  = gmx_simd_mul_d(PolyCP0, w2);
    PolyCP1  = gmx_simd_mul_d(PolyCP1, w);
    PolyCP0  = gmx_simd_add_d(PolyCP0, CCP0);
    PolyCP0  = gmx_simd_add_d(PolyCP0, PolyCP1);

    PolyCQ0  = gmx_simd_mul_d(CCQ6, w2);
    PolyCQ1  = gmx_simd_mul_d(CCQ5, w2);
    PolyCQ0  = gmx_simd_add_d(PolyCQ0, CCQ4);
    PolyCQ1  = gmx_simd_add_d(PolyCQ1, CCQ3);
    PolyCQ0  = gmx_simd_mul_d(PolyCQ0, w2);
    PolyCQ1  = gmx_simd_mul_d(PolyCQ1, w2);
    PolyCQ0  = gmx_simd_add_d(PolyCQ0, CCQ2);
    PolyCQ1  = gmx_simd_add_d(PolyCQ1, CCQ1);
    PolyCQ0  = gmx_simd_mul_d(PolyCQ0, w2);
    PolyCQ1  = gmx_simd_mul_d(PolyCQ1, w);
    PolyCQ0  = gmx_simd_add_d(PolyCQ0, one);
    PolyCQ0  = gmx_simd_add_d(PolyCQ0, PolyCQ1);

    expmx2   = gmx_simd_exp_d( gmx_simd_fneg_d(x2) );

    res_erfcC = gmx_simd_mul_d(PolyCP0, gmx_simd_inv_notmaskfpe_d(PolyCQ0, mask_erf));
    res_erfcC = gmx_simd_add_d(res_erfcC, CCoffset);
    res_erfcC = gmx_simd_mul_d(res_erfcC, w);

    mask     = gmx_simd_cmplt_d(gmx_simd_set1_d(4.5), xabs);
    res_erfc = gmx_simd_blendv_d(res_erfcB, res_erfcC, mask);

    res_erfc = gmx_simd_mul_d(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = gmx_simd_cmplt_d(x, gmx_simd_setzero_d());
    res_erfc = gmx_simd_blendv_d(res_erfc, gmx_simd_sub_d(two, res_erfc), mask);

    /* Select erf() or erfc() */
    res  = gmx_simd_blendv_d(res_erfc, gmx_simd_sub_d(one, res_erf), mask_erf);

    return res;
}

/*! \brief SIMD double sin \& cos.
 *
 * \copydetails gmx_simd_sincos_f
 */
static gmx_inline void gmx_simdcall
gmx_simd_sincos_d(gmx_simd_double_t x, gmx_simd_double_t *sinval, gmx_simd_double_t *cosval)
{
    /* Constants to subtract Pi/4*x from y while minimizing precision loss */
    const gmx_simd_double_t  argred0         = gmx_simd_set1_d(-2*0.78539816290140151978);
    const gmx_simd_double_t  argred1         = gmx_simd_set1_d(-2*4.9604678871439933374e-10);
    const gmx_simd_double_t  argred2         = gmx_simd_set1_d(-2*1.1258708853173288931e-18);
    const gmx_simd_double_t  argred3         = gmx_simd_set1_d(-2*1.7607799325916000908e-27);
    const gmx_simd_double_t  two_over_pi     = gmx_simd_set1_d(2.0/M_PI);
    const gmx_simd_double_t  const_sin5      = gmx_simd_set1_d( 1.58938307283228937328511e-10);
    const gmx_simd_double_t  const_sin4      = gmx_simd_set1_d(-2.50506943502539773349318e-08);
    const gmx_simd_double_t  const_sin3      = gmx_simd_set1_d( 2.75573131776846360512547e-06);
    const gmx_simd_double_t  const_sin2      = gmx_simd_set1_d(-0.000198412698278911770864914);
    const gmx_simd_double_t  const_sin1      = gmx_simd_set1_d( 0.0083333333333191845961746);
    const gmx_simd_double_t  const_sin0      = gmx_simd_set1_d(-0.166666666666666130709393);

    const gmx_simd_double_t  const_cos7      = gmx_simd_set1_d(-1.13615350239097429531523e-11);
    const gmx_simd_double_t  const_cos6      = gmx_simd_set1_d( 2.08757471207040055479366e-09);
    const gmx_simd_double_t  const_cos5      = gmx_simd_set1_d(-2.75573144028847567498567e-07);
    const gmx_simd_double_t  const_cos4      = gmx_simd_set1_d( 2.48015872890001867311915e-05);
    const gmx_simd_double_t  const_cos3      = gmx_simd_set1_d(-0.00138888888888714019282329);
    const gmx_simd_double_t  const_cos2      = gmx_simd_set1_d( 0.0416666666666665519592062);
    const gmx_simd_double_t  half            = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t  one             = gmx_simd_set1_d(1.0);
    gmx_simd_double_t        ssign, csign;
    gmx_simd_double_t        x2, y, z, psin, pcos, sss, ccc;
    gmx_simd_dbool_t         mask;
#if (defined GMX_SIMD_HAVE_DINT32) && (defined GMX_SIMD_HAVE_DINT32_ARITHMETICS) && (defined GMX_SIMD_HAVE_LOGICAL)
    const gmx_simd_dint32_t  ione            = gmx_simd_set1_di(1);
    const gmx_simd_dint32_t  itwo            = gmx_simd_set1_di(2);
    gmx_simd_dint32_t        iy;

    z       = gmx_simd_mul_d(x, two_over_pi);
    iy      = gmx_simd_cvt_d2i(z);
    y       = gmx_simd_round_d(z);

    mask    = gmx_simd_cvt_dib2db(gmx_simd_cmpeq_di(gmx_simd_and_di(iy, ione), gmx_simd_setzero_di()));
    ssign   = gmx_simd_blendzero_d(gmx_simd_set1_d(GMX_DOUBLE_NEGZERO), gmx_simd_cvt_dib2db(gmx_simd_cmpeq_di(gmx_simd_and_di(iy, itwo), itwo)));
    csign   = gmx_simd_blendzero_d(gmx_simd_set1_d(GMX_DOUBLE_NEGZERO), gmx_simd_cvt_dib2db(gmx_simd_cmpeq_di(gmx_simd_and_di(gmx_simd_add_di(iy, ione), itwo), itwo)));
#else
    const gmx_simd_double_t  quarter         = gmx_simd_set1_d(0.25);
    const gmx_simd_double_t  minusquarter    = gmx_simd_set1_d(-0.25);
    gmx_simd_double_t        q;
    gmx_simd_dbool_t         m1, m2, m3;

    /* The most obvious way to find the arguments quadrant in the unit circle
     * to calculate the sign is to use integer arithmetic, but that is not
     * present in all SIMD implementations. As an alternative, we have devised a
     * pure floating-point algorithm that uses truncation for argument reduction
     * so that we get a new value 0<=q<1 over the unit circle, and then
     * do floating-point comparisons with fractions. This is likely to be
     * slightly slower (~10%) due to the longer latencies of floating-point, so
     * we only use it when integer SIMD arithmetic is not present.
     */
    ssign   = x;
    x       = gmx_simd_fabs_d(x);
    /* It is critical that half-way cases are rounded down */
    z       = gmx_simd_fmadd_d(x, two_over_pi, half);
    y       = gmx_simd_trunc_d(z);
    q       = gmx_simd_mul_d(z, quarter);
    q       = gmx_simd_sub_d(q, gmx_simd_trunc_d(q));
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q       = gmx_simd_sub_d(q, half);
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#    ifdef GMX_SIMD_HAVE_LOGICAL
    ssign   = gmx_simd_and_d(ssign, gmx_simd_set1_d(GMX_DOUBLE_NEGZERO));
    csign   = gmx_simd_andnot_d(q, gmx_simd_set1_d(GMX_DOUBLE_NEGZERO));
    ssign   = gmx_simd_xor_d(ssign, csign);
#    else
    csign   = gmx_simd_xor_sign_d(gmx_simd_set1_d(-1.0), q);
    ssign   = gmx_simd_xor_sign_d(ssign, csign);    /* swap ssign if csign was set. */
#    endif
    /* Check if y had value 1 or 3 (remember we subtracted 0.5 from q) */
    m1      = gmx_simd_cmplt_d(q, minusquarter);
    m2      = gmx_simd_cmple_d(gmx_simd_setzero_d(), q);
    m3      = gmx_simd_cmplt_d(q, quarter);
    m2      = gmx_simd_and_db(m2, m3);
    mask    = gmx_simd_or_db(m1, m2);
    /* where mask is FALSE, set sign. */
    csign   = gmx_simd_xor_sign_d(csign, gmx_simd_blendv_d(gmx_simd_set1_d(-1.0), one, mask));
#endif
    x       = gmx_simd_fmadd_d(y, argred0, x);
    x       = gmx_simd_fmadd_d(y, argred1, x);
    x       = gmx_simd_fmadd_d(y, argred2, x);
    x       = gmx_simd_fmadd_d(y, argred3, x);
    x2      = gmx_simd_mul_d(x, x);

    psin    = gmx_simd_fmadd_d(const_sin5, x2, const_sin4);
    psin    = gmx_simd_fmadd_d(psin, x2, const_sin3);
    psin    = gmx_simd_fmadd_d(psin, x2, const_sin2);
    psin    = gmx_simd_fmadd_d(psin, x2, const_sin1);
    psin    = gmx_simd_fmadd_d(psin, x2, const_sin0);
    psin    = gmx_simd_fmadd_d(psin, gmx_simd_mul_d(x2, x), x);

    pcos    = gmx_simd_fmadd_d(const_cos7, x2, const_cos6);
    pcos    = gmx_simd_fmadd_d(pcos, x2, const_cos5);
    pcos    = gmx_simd_fmadd_d(pcos, x2, const_cos4);
    pcos    = gmx_simd_fmadd_d(pcos, x2, const_cos3);
    pcos    = gmx_simd_fmadd_d(pcos, x2, const_cos2);
    pcos    = gmx_simd_fmsub_d(pcos, x2, half);
    pcos    = gmx_simd_fmadd_d(pcos, x2, one);

    sss     = gmx_simd_blendv_d(pcos, psin, mask);
    ccc     = gmx_simd_blendv_d(psin, pcos, mask);
    /* See comment for GMX_SIMD_HAVE_LOGICAL section above. */
#ifdef GMX_SIMD_HAVE_LOGICAL
    *sinval = gmx_simd_xor_d(sss, ssign);
    *cosval = gmx_simd_xor_d(ccc, csign);
#else
    *sinval = gmx_simd_xor_sign_d(sss, ssign);
    *cosval = gmx_simd_xor_sign_d(ccc, csign);
#endif
}

/*! \brief SIMD double sin(x).
 *
 * \copydetails gmx_simd_sin_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_sin_d(gmx_simd_double_t x)
{
    gmx_simd_double_t s, c;
    gmx_simd_sincos_d(x, &s, &c);
    return s;
}

/*! \brief SIMD double cos(x).
 *
 * \copydetails gmx_simd_cos_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_cos_d(gmx_simd_double_t x)
{
    gmx_simd_double_t s, c;
    gmx_simd_sincos_d(x, &s, &c);
    return c;
}

/*! \brief SIMD double tan(x).
 *
 * \copydetails gmx_simd_tan_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_tan_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t  argred0         = gmx_simd_set1_d(-2*0.78539816290140151978);
    const gmx_simd_double_t  argred1         = gmx_simd_set1_d(-2*4.9604678871439933374e-10);
    const gmx_simd_double_t  argred2         = gmx_simd_set1_d(-2*1.1258708853173288931e-18);
    const gmx_simd_double_t  argred3         = gmx_simd_set1_d(-2*1.7607799325916000908e-27);
    const gmx_simd_double_t  two_over_pi     = gmx_simd_set1_d(2.0/M_PI);
    const gmx_simd_double_t  CT15            = gmx_simd_set1_d(1.01419718511083373224408e-05);
    const gmx_simd_double_t  CT14            = gmx_simd_set1_d(-2.59519791585924697698614e-05);
    const gmx_simd_double_t  CT13            = gmx_simd_set1_d(5.23388081915899855325186e-05);
    const gmx_simd_double_t  CT12            = gmx_simd_set1_d(-3.05033014433946488225616e-05);
    const gmx_simd_double_t  CT11            = gmx_simd_set1_d(7.14707504084242744267497e-05);
    const gmx_simd_double_t  CT10            = gmx_simd_set1_d(8.09674518280159187045078e-05);
    const gmx_simd_double_t  CT9             = gmx_simd_set1_d(0.000244884931879331847054404);
    const gmx_simd_double_t  CT8             = gmx_simd_set1_d(0.000588505168743587154904506);
    const gmx_simd_double_t  CT7             = gmx_simd_set1_d(0.00145612788922812427978848);
    const gmx_simd_double_t  CT6             = gmx_simd_set1_d(0.00359208743836906619142924);
    const gmx_simd_double_t  CT5             = gmx_simd_set1_d(0.00886323944362401618113356);
    const gmx_simd_double_t  CT4             = gmx_simd_set1_d(0.0218694882853846389592078);
    const gmx_simd_double_t  CT3             = gmx_simd_set1_d(0.0539682539781298417636002);
    const gmx_simd_double_t  CT2             = gmx_simd_set1_d(0.133333333333125941821962);
    const gmx_simd_double_t  CT1             = gmx_simd_set1_d(0.333333333333334980164153);

    gmx_simd_double_t        x2, p, y, z;
    gmx_simd_dbool_t         mask;

#if (defined GMX_SIMD_HAVE_DINT32) && (defined GMX_SIMD_HAVE_DINT32_ARITHMETICS) && (defined GMX_SIMD_HAVE_LOGICAL)
    gmx_simd_dint32_t  iy;
    gmx_simd_dint32_t  ione = gmx_simd_set1_di(1);

    z       = gmx_simd_mul_d(x, two_over_pi);
    iy      = gmx_simd_cvt_d2i(z);
    y       = gmx_simd_round_d(z);
    mask    = gmx_simd_cvt_dib2db(gmx_simd_cmpeq_di(gmx_simd_and_di(iy, ione), ione));

    x       = gmx_simd_fmadd_d(y, argred0, x);
    x       = gmx_simd_fmadd_d(y, argred1, x);
    x       = gmx_simd_fmadd_d(y, argred2, x);
    x       = gmx_simd_fmadd_d(y, argred3, x);
    x       = gmx_simd_xor_d(gmx_simd_blendzero_d(gmx_simd_set1_d(GMX_DOUBLE_NEGZERO), mask), x);
#else
    const gmx_simd_double_t  quarter         = gmx_simd_set1_d(0.25);
    const gmx_simd_double_t  half            = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t  threequarter    = gmx_simd_set1_d(0.75);
    gmx_simd_double_t        w, q;
    gmx_simd_dbool_t         m1, m2, m3;

    w       = gmx_simd_fabs_d(x);
    z       = gmx_simd_fmadd_d(w, two_over_pi, half);
    y       = gmx_simd_trunc_d(z);
    q       = gmx_simd_mul_d(z, quarter);
    q       = gmx_simd_sub_d(q, gmx_simd_trunc_d(q));
    m1      = gmx_simd_cmple_d(quarter, q);
    m2      = gmx_simd_cmplt_d(q, half);
    m3      = gmx_simd_cmple_d(threequarter, q);
    m1      = gmx_simd_and_db(m1, m2);
    mask    = gmx_simd_or_db(m1, m3);
    w       = gmx_simd_fmadd_d(y, argred0, w);
    w       = gmx_simd_fmadd_d(y, argred1, w);
    w       = gmx_simd_fmadd_d(y, argred2, w);
    w       = gmx_simd_fmadd_d(y, argred3, w);

    w       = gmx_simd_blendv_d(w, gmx_simd_fneg_d(w), mask);
    x       = gmx_simd_xor_sign_d(w, x);
#endif
    x2      = gmx_simd_mul_d(x, x);
    p       = gmx_simd_fmadd_d(CT15, x2, CT14);
    p       = gmx_simd_fmadd_d(p, x2, CT13);
    p       = gmx_simd_fmadd_d(p, x2, CT12);
    p       = gmx_simd_fmadd_d(p, x2, CT11);
    p       = gmx_simd_fmadd_d(p, x2, CT10);
    p       = gmx_simd_fmadd_d(p, x2, CT9);
    p       = gmx_simd_fmadd_d(p, x2, CT8);
    p       = gmx_simd_fmadd_d(p, x2, CT7);
    p       = gmx_simd_fmadd_d(p, x2, CT6);
    p       = gmx_simd_fmadd_d(p, x2, CT5);
    p       = gmx_simd_fmadd_d(p, x2, CT4);
    p       = gmx_simd_fmadd_d(p, x2, CT3);
    p       = gmx_simd_fmadd_d(p, x2, CT2);
    p       = gmx_simd_fmadd_d(p, x2, CT1);
    p       = gmx_simd_fmadd_d(x2, gmx_simd_mul_d(p, x), x);

    p       = gmx_simd_blendv_d( p, gmx_simd_inv_maskfpe_d(p, mask), mask);
    return p;
}

/*! \brief SIMD double asin(x).
 *
 * \copydetails gmx_simd_asin_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_asin_d(gmx_simd_double_t x)
{
    /* Same algorithm as cephes library */
    const gmx_simd_double_t limit1    = gmx_simd_set1_d(0.625);
    const gmx_simd_double_t limit2    = gmx_simd_set1_d(1e-8);
    const gmx_simd_double_t one       = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t quarterpi = gmx_simd_set1_d(M_PI/4.0);
    const gmx_simd_double_t morebits  = gmx_simd_set1_d(6.123233995736765886130e-17);

    const gmx_simd_double_t P5        = gmx_simd_set1_d(4.253011369004428248960e-3);
    const gmx_simd_double_t P4        = gmx_simd_set1_d(-6.019598008014123785661e-1);
    const gmx_simd_double_t P3        = gmx_simd_set1_d(5.444622390564711410273e0);
    const gmx_simd_double_t P2        = gmx_simd_set1_d(-1.626247967210700244449e1);
    const gmx_simd_double_t P1        = gmx_simd_set1_d(1.956261983317594739197e1);
    const gmx_simd_double_t P0        = gmx_simd_set1_d(-8.198089802484824371615e0);

    const gmx_simd_double_t Q4        = gmx_simd_set1_d(-1.474091372988853791896e1);
    const gmx_simd_double_t Q3        = gmx_simd_set1_d(7.049610280856842141659e1);
    const gmx_simd_double_t Q2        = gmx_simd_set1_d(-1.471791292232726029859e2);
    const gmx_simd_double_t Q1        = gmx_simd_set1_d(1.395105614657485689735e2);
    const gmx_simd_double_t Q0        = gmx_simd_set1_d(-4.918853881490881290097e1);

    const gmx_simd_double_t R4        = gmx_simd_set1_d(2.967721961301243206100e-3);
    const gmx_simd_double_t R3        = gmx_simd_set1_d(-5.634242780008963776856e-1);
    const gmx_simd_double_t R2        = gmx_simd_set1_d(6.968710824104713396794e0);
    const gmx_simd_double_t R1        = gmx_simd_set1_d(-2.556901049652824852289e1);
    const gmx_simd_double_t R0        = gmx_simd_set1_d(2.853665548261061424989e1);

    const gmx_simd_double_t S3        = gmx_simd_set1_d(-2.194779531642920639778e1);
    const gmx_simd_double_t S2        = gmx_simd_set1_d(1.470656354026814941758e2);
    const gmx_simd_double_t S1        = gmx_simd_set1_d(-3.838770957603691357202e2);
    const gmx_simd_double_t S0        = gmx_simd_set1_d(3.424398657913078477438e2);

    gmx_simd_double_t       xabs;
    gmx_simd_double_t       zz, ww, z, q, w, zz2, ww2;
    gmx_simd_double_t       PA, PB;
    gmx_simd_double_t       QA, QB;
    gmx_simd_double_t       RA, RB;
    gmx_simd_double_t       SA, SB;
    gmx_simd_double_t       nom, denom;
    gmx_simd_dbool_t        mask, mask2;

    xabs  = gmx_simd_fabs_d(x);

    mask  = gmx_simd_cmplt_d(limit1, xabs);

    zz    = gmx_simd_sub_d(one, xabs);
    ww    = gmx_simd_mul_d(xabs, xabs);
    zz2   = gmx_simd_mul_d(zz, zz);
    ww2   = gmx_simd_mul_d(ww, ww);

    /* R */
    RA    = gmx_simd_mul_d(R4, zz2);
    RB    = gmx_simd_mul_d(R3, zz2);
    RA    = gmx_simd_add_d(RA, R2);
    RB    = gmx_simd_add_d(RB, R1);
    RA    = gmx_simd_mul_d(RA, zz2);
    RB    = gmx_simd_mul_d(RB, zz);
    RA    = gmx_simd_add_d(RA, R0);
    RA    = gmx_simd_add_d(RA, RB);

    /* S, SA = zz2 */
    SB    = gmx_simd_mul_d(S3, zz2);
    SA    = gmx_simd_add_d(zz2, S2);
    SB    = gmx_simd_add_d(SB, S1);
    SA    = gmx_simd_mul_d(SA, zz2);
    SB    = gmx_simd_mul_d(SB, zz);
    SA    = gmx_simd_add_d(SA, S0);
    SA    = gmx_simd_add_d(SA, SB);

    /* P */
    PA    = gmx_simd_mul_d(P5, ww2);
    PB    = gmx_simd_mul_d(P4, ww2);
    PA    = gmx_simd_add_d(PA, P3);
    PB    = gmx_simd_add_d(PB, P2);
    PA    = gmx_simd_mul_d(PA, ww2);
    PB    = gmx_simd_mul_d(PB, ww2);
    PA    = gmx_simd_add_d(PA, P1);
    PB    = gmx_simd_add_d(PB, P0);
    PA    = gmx_simd_mul_d(PA, ww);
    PA    = gmx_simd_add_d(PA, PB);

    /* Q, QA = ww2 */
    QB    = gmx_simd_mul_d(Q4, ww2);
    QA    = gmx_simd_add_d(ww2, Q3);
    QB    = gmx_simd_add_d(QB, Q2);
    QA    = gmx_simd_mul_d(QA, ww2);
    QB    = gmx_simd_mul_d(QB, ww2);
    QA    = gmx_simd_add_d(QA, Q1);
    QB    = gmx_simd_add_d(QB, Q0);
    QA    = gmx_simd_mul_d(QA, ww);
    QA    = gmx_simd_add_d(QA, QB);

    RA    = gmx_simd_mul_d(RA, zz);
    PA    = gmx_simd_mul_d(PA, ww);

    nom   = gmx_simd_blendv_d( PA, RA, mask );
    denom = gmx_simd_blendv_d( QA, SA, mask );

    mask2 = gmx_simd_cmplt_d(limit2, xabs);
    q     = gmx_simd_mul_d( nom, gmx_simd_inv_maskfpe_d(denom, mask2) );

    zz    = gmx_simd_add_d(zz, zz);
    zz    = gmx_simd_sqrt_d(zz);
    z     = gmx_simd_sub_d(quarterpi, zz);
    zz    = gmx_simd_mul_d(zz, q);
    zz    = gmx_simd_sub_d(zz, morebits);
    z     = gmx_simd_sub_d(z, zz);
    z     = gmx_simd_add_d(z, quarterpi);

    w     = gmx_simd_mul_d(xabs, q);
    w     = gmx_simd_add_d(w, xabs);

    z     = gmx_simd_blendv_d( w, z, mask );

    z     = gmx_simd_blendv_d( xabs, z, mask2 );

    z = gmx_simd_xor_sign_d(z, x);

    return z;
}

/*! \brief SIMD double acos(x).
 *
 * \copydetails gmx_simd_acos_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_acos_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t one        = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t half       = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t quarterpi0 = gmx_simd_set1_d(7.85398163397448309616e-1);
    const gmx_simd_double_t quarterpi1 = gmx_simd_set1_d(6.123233995736765886130e-17);

    gmx_simd_dbool_t        mask1;
    gmx_simd_double_t       z, z1, z2;

    mask1 = gmx_simd_cmplt_d(half, x);
    z1    = gmx_simd_mul_d(half, gmx_simd_sub_d(one, x));
    z1    = gmx_simd_sqrt_d(z1);
    z     = gmx_simd_blendv_d( x, z1, mask1 );

    z     = gmx_simd_asin_d(z);

    z1    = gmx_simd_add_d(z, z);

    z2    = gmx_simd_sub_d(quarterpi0, z);
    z2    = gmx_simd_add_d(z2, quarterpi1);
    z2    = gmx_simd_add_d(z2, quarterpi0);

    z     = gmx_simd_blendv_d(z2, z1, mask1);

    return z;
}

/*! \brief SIMD double atan(x).
 *
 * \copydetails gmx_simd_atan_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_atan_d(gmx_simd_double_t x)
{
    /* Same algorithm as cephes library */
    const gmx_simd_double_t limit1    = gmx_simd_set1_d(0.66);
    const gmx_simd_double_t limit2    = gmx_simd_set1_d(2.41421356237309504880);
    const gmx_simd_double_t quarterpi = gmx_simd_set1_d(M_PI/4.0);
    const gmx_simd_double_t halfpi    = gmx_simd_set1_d(M_PI/2.0);
    const gmx_simd_double_t mone      = gmx_simd_set1_d(-1.0);
    const gmx_simd_double_t morebits1 = gmx_simd_set1_d(0.5*6.123233995736765886130E-17);
    const gmx_simd_double_t morebits2 = gmx_simd_set1_d(6.123233995736765886130E-17);

    const gmx_simd_double_t P4        = gmx_simd_set1_d(-8.750608600031904122785E-1);
    const gmx_simd_double_t P3        = gmx_simd_set1_d(-1.615753718733365076637E1);
    const gmx_simd_double_t P2        = gmx_simd_set1_d(-7.500855792314704667340E1);
    const gmx_simd_double_t P1        = gmx_simd_set1_d(-1.228866684490136173410E2);
    const gmx_simd_double_t P0        = gmx_simd_set1_d(-6.485021904942025371773E1);

    const gmx_simd_double_t Q4        = gmx_simd_set1_d(2.485846490142306297962E1);
    const gmx_simd_double_t Q3        = gmx_simd_set1_d(1.650270098316988542046E2);
    const gmx_simd_double_t Q2        = gmx_simd_set1_d(4.328810604912902668951E2);
    const gmx_simd_double_t Q1        = gmx_simd_set1_d(4.853903996359136964868E2);
    const gmx_simd_double_t Q0        = gmx_simd_set1_d(1.945506571482613964425E2);

    gmx_simd_double_t       y, xabs, t1, t2;
    gmx_simd_double_t       z, z2;
    gmx_simd_double_t       P_A, P_B, Q_A, Q_B;
    gmx_simd_dbool_t        mask1, mask2;

    xabs   = gmx_simd_fabs_d(x);

    mask1  = gmx_simd_cmplt_d(limit1, xabs);
    mask2  = gmx_simd_cmplt_d(limit2, xabs);

    t1     = gmx_simd_mul_d(gmx_simd_add_d(xabs, mone),
                            gmx_simd_inv_maskfpe_d(gmx_simd_sub_d(xabs, mone), mask1));
    t2     = gmx_simd_mul_d(mone, gmx_simd_inv_maskfpe_d(xabs, mask2));

    y      = gmx_simd_blendzero_d(quarterpi, mask1);
    y      = gmx_simd_blendv_d(y, halfpi, mask2);
    xabs   = gmx_simd_blendv_d(xabs, t1, mask1);
    xabs   = gmx_simd_blendv_d(xabs, t2, mask2);

    z      = gmx_simd_mul_d(xabs, xabs);
    z2     = gmx_simd_mul_d(z, z);

    P_A    = gmx_simd_mul_d(P4, z2);
    P_B    = gmx_simd_mul_d(P3, z2);
    P_A    = gmx_simd_add_d(P_A, P2);
    P_B    = gmx_simd_add_d(P_B, P1);
    P_A    = gmx_simd_mul_d(P_A, z2);
    P_B    = gmx_simd_mul_d(P_B, z);
    P_A    = gmx_simd_add_d(P_A, P0);
    P_A    = gmx_simd_add_d(P_A, P_B);

    /* Q_A = z2 */
    Q_B    = gmx_simd_mul_d(Q4, z2);
    Q_A    = gmx_simd_add_d(z2, Q3);
    Q_B    = gmx_simd_add_d(Q_B, Q2);
    Q_A    = gmx_simd_mul_d(Q_A, z2);
    Q_B    = gmx_simd_mul_d(Q_B, z2);
    Q_A    = gmx_simd_add_d(Q_A, Q1);
    Q_B    = gmx_simd_add_d(Q_B, Q0);
    Q_A    = gmx_simd_mul_d(Q_A, z);
    Q_A    = gmx_simd_add_d(Q_A, Q_B);

    z      = gmx_simd_mul_d(z, P_A);
    z      = gmx_simd_mul_d(z, gmx_simd_inv_d(Q_A));
    z      = gmx_simd_mul_d(z, xabs);
    z      = gmx_simd_add_d(z, xabs);

    t1     = gmx_simd_blendzero_d(morebits1, mask1);
    t1     = gmx_simd_blendv_d(t1, morebits2, mask2);

    z      = gmx_simd_add_d(z, t1);
    y      = gmx_simd_add_d(y, z);

    y      = gmx_simd_xor_sign_d(y, x);

    return y;
}

/*! \brief SIMD double atan2(y,x).
 *
 * \copydetails gmx_simd_atan2_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_atan2_d(gmx_simd_double_t y, gmx_simd_double_t x)
{
    const gmx_simd_double_t pi          = gmx_simd_set1_d(M_PI);
    const gmx_simd_double_t halfpi      = gmx_simd_set1_d(M_PI/2.0);
    gmx_simd_double_t       xinv, p, aoffset;
    gmx_simd_dbool_t        mask_x0, mask_y0, mask_xlt0, mask_ylt0;

    mask_x0   = gmx_simd_cmpeq_d(x, gmx_simd_setzero_d());
    mask_y0   = gmx_simd_cmpeq_d(y, gmx_simd_setzero_d());
    mask_xlt0 = gmx_simd_cmplt_d(x, gmx_simd_setzero_d());
    mask_ylt0 = gmx_simd_cmplt_d(y, gmx_simd_setzero_d());

    aoffset   = gmx_simd_blendzero_d(halfpi, mask_x0);
    aoffset   = gmx_simd_blendnotzero_d(aoffset, mask_y0);

    aoffset   = gmx_simd_blendv_d(aoffset, pi, mask_xlt0);
    aoffset   = gmx_simd_blendv_d(aoffset, gmx_simd_fneg_d(aoffset), mask_ylt0);

    xinv      = gmx_simd_blendnotzero_d(gmx_simd_inv_notmaskfpe_d(x, mask_x0), mask_x0);
    p         = gmx_simd_mul_d(y, xinv);
    p         = gmx_simd_atan_d(p);
    p         = gmx_simd_add_d(p, aoffset);

    return p;
}


/*! \brief Calculate the force correction due to PME analytically for SIMD double.
 *
 * \copydetails gmx_simd_pmecorrF_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_pmecorrF_d(gmx_simd_double_t z2)
{
    const gmx_simd_double_t  FN10     = gmx_simd_set1_d(-8.0072854618360083154e-14);
    const gmx_simd_double_t  FN9      = gmx_simd_set1_d(1.1859116242260148027e-11);
    const gmx_simd_double_t  FN8      = gmx_simd_set1_d(-8.1490406329798423616e-10);
    const gmx_simd_double_t  FN7      = gmx_simd_set1_d(3.4404793543907847655e-8);
    const gmx_simd_double_t  FN6      = gmx_simd_set1_d(-9.9471420832602741006e-7);
    const gmx_simd_double_t  FN5      = gmx_simd_set1_d(0.000020740315999115847456);
    const gmx_simd_double_t  FN4      = gmx_simd_set1_d(-0.00031991745139313364005);
    const gmx_simd_double_t  FN3      = gmx_simd_set1_d(0.0035074449373659008203);
    const gmx_simd_double_t  FN2      = gmx_simd_set1_d(-0.031750380176100813405);
    const gmx_simd_double_t  FN1      = gmx_simd_set1_d(0.13884101728898463426);
    const gmx_simd_double_t  FN0      = gmx_simd_set1_d(-0.75225277815249618847);

    const gmx_simd_double_t  FD5      = gmx_simd_set1_d(0.000016009278224355026701);
    const gmx_simd_double_t  FD4      = gmx_simd_set1_d(0.00051055686934806966046);
    const gmx_simd_double_t  FD3      = gmx_simd_set1_d(0.0081803507497974289008);
    const gmx_simd_double_t  FD2      = gmx_simd_set1_d(0.077181146026670287235);
    const gmx_simd_double_t  FD1      = gmx_simd_set1_d(0.41543303143712535988);
    const gmx_simd_double_t  FD0      = gmx_simd_set1_d(1.0);

    gmx_simd_double_t        z4;
    gmx_simd_double_t        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = gmx_simd_mul_d(z2, z2);

    polyFD1        = gmx_simd_fmadd_d(FD5, z4, FD3);
    polyFD1        = gmx_simd_fmadd_d(polyFD1, z4, FD1);
    polyFD1        = gmx_simd_mul_d(polyFD1, z2);
    polyFD0        = gmx_simd_fmadd_d(FD4, z4, FD2);
    polyFD0        = gmx_simd_fmadd_d(polyFD0, z4, FD0);
    polyFD0        = gmx_simd_add_d(polyFD0, polyFD1);

    polyFD0        = gmx_simd_inv_d(polyFD0);

    polyFN0        = gmx_simd_fmadd_d(FN10, z4, FN8);
    polyFN0        = gmx_simd_fmadd_d(polyFN0, z4, FN6);
    polyFN0        = gmx_simd_fmadd_d(polyFN0, z4, FN4);
    polyFN0        = gmx_simd_fmadd_d(polyFN0, z4, FN2);
    polyFN0        = gmx_simd_fmadd_d(polyFN0, z4, FN0);
    polyFN1        = gmx_simd_fmadd_d(FN9, z4, FN7);
    polyFN1        = gmx_simd_fmadd_d(polyFN1, z4, FN5);
    polyFN1        = gmx_simd_fmadd_d(polyFN1, z4, FN3);
    polyFN1        = gmx_simd_fmadd_d(polyFN1, z4, FN1);
    polyFN0        = gmx_simd_fmadd_d(polyFN1, z2, polyFN0);


    return gmx_simd_mul_d(polyFN0, polyFD0);
}



/*! \brief Calculate the potential correction due to PME analytically for SIMD double.
 *
 * \copydetails gmx_simd_pmecorrV_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_pmecorrV_d(gmx_simd_double_t z2)
{
    const gmx_simd_double_t  VN9      = gmx_simd_set1_d(-9.3723776169321855475e-13);
    const gmx_simd_double_t  VN8      = gmx_simd_set1_d(1.2280156762674215741e-10);
    const gmx_simd_double_t  VN7      = gmx_simd_set1_d(-7.3562157912251309487e-9);
    const gmx_simd_double_t  VN6      = gmx_simd_set1_d(2.6215886208032517509e-7);
    const gmx_simd_double_t  VN5      = gmx_simd_set1_d(-4.9532491651265819499e-6);
    const gmx_simd_double_t  VN4      = gmx_simd_set1_d(0.00025907400778966060389);
    const gmx_simd_double_t  VN3      = gmx_simd_set1_d(0.0010585044856156469792);
    const gmx_simd_double_t  VN2      = gmx_simd_set1_d(0.045247661136833092885);
    const gmx_simd_double_t  VN1      = gmx_simd_set1_d(0.11643931522926034421);
    const gmx_simd_double_t  VN0      = gmx_simd_set1_d(1.1283791671726767970);

    const gmx_simd_double_t  VD5      = gmx_simd_set1_d(0.000021784709867336150342);
    const gmx_simd_double_t  VD4      = gmx_simd_set1_d(0.00064293662010911388448);
    const gmx_simd_double_t  VD3      = gmx_simd_set1_d(0.0096311444822588683504);
    const gmx_simd_double_t  VD2      = gmx_simd_set1_d(0.085608012351550627051);
    const gmx_simd_double_t  VD1      = gmx_simd_set1_d(0.43652499166614811084);
    const gmx_simd_double_t  VD0      = gmx_simd_set1_d(1.0);

    gmx_simd_double_t        z4;
    gmx_simd_double_t        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = gmx_simd_mul_d(z2, z2);

    polyVD1        = gmx_simd_fmadd_d(VD5, z4, VD3);
    polyVD0        = gmx_simd_fmadd_d(VD4, z4, VD2);
    polyVD1        = gmx_simd_fmadd_d(polyVD1, z4, VD1);
    polyVD0        = gmx_simd_fmadd_d(polyVD0, z4, VD0);
    polyVD0        = gmx_simd_fmadd_d(polyVD1, z2, polyVD0);

    polyVD0        = gmx_simd_inv_d(polyVD0);

    polyVN1        = gmx_simd_fmadd_d(VN9, z4, VN7);
    polyVN0        = gmx_simd_fmadd_d(VN8, z4, VN6);
    polyVN1        = gmx_simd_fmadd_d(polyVN1, z4, VN5);
    polyVN0        = gmx_simd_fmadd_d(polyVN0, z4, VN4);
    polyVN1        = gmx_simd_fmadd_d(polyVN1, z4, VN3);
    polyVN0        = gmx_simd_fmadd_d(polyVN0, z4, VN2);
    polyVN1        = gmx_simd_fmadd_d(polyVN1, z4, VN1);
    polyVN0        = gmx_simd_fmadd_d(polyVN0, z4, VN0);
    polyVN0        = gmx_simd_fmadd_d(polyVN1, z2, polyVN0);

    return gmx_simd_mul_d(polyVN0, polyVD0);
}

/*! \} */


/*! \name SIMD math functions for double prec. data, single prec. accuracy
 *
 *  \note In some cases we do not need full double accuracy of individual
 *        SIMD math functions, although the data is stored in double precision
 *        SIMD registers. This might be the case for special algorithms, or
 *        if the architecture does not support single precision.
 *        Since the full double precision evaluation of math functions
 *        typically require much more expensive polynomial approximations
 *        these functions implement the algorithms used in the single precision
 *        SIMD math functions, but they operate on double precision
 *        SIMD variables.
 *
 *  \note You should normally not use these functions directly, but the
 *        real-precision wrappers instead. When Gromacs is compiled in single
 *        precision, those will be aliases to the normal single precision
 *        SIMD math functions.
 *  \{
 */

/*********************************************************************
 * SIMD MATH FUNCTIONS WITH DOUBLE PREC. DATA, SINGLE PREC. ACCURACY *
 *********************************************************************/

/*! \brief Calculate 1/sqrt(x) for SIMD double, but in single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_invsqrt_singleaccuracy_r.
 *
 *  \param x Argument that must be >0. This routine does not check arguments.
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_invsqrt_singleaccuracy_d(gmx_simd_double_t x)
{
    gmx_simd_double_t lu = gmx_simd_rsqrt_d(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rsqrt_iter_d(lu, x);
#endif
    return lu;
}

/*! \brief 1/sqrt(x) for masked entries of SIMD double, but in single accuracy.
 *
 * \copydetails gmx_simd_invsqrt_maskfpe_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_invsqrt_maskfpe_singleaccuracy_d(gmx_simd_double_t x, gmx_simd_dbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_invsqrt_singleaccuracy_d(x);
#else
    return gmx_simd_invsqrt_singleaccuracy_d(gmx_simd_blendv_d(gmx_simd_set1_d(1.0), x, m));
#endif
}

/*! \brief 1/sqrt(x) for non-masked entries of SIMD double, in single accuracy.
 *
 * \copydetails gmx_simd_invsqrt_notmaskfpe_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_invsqrt_notmaskfpe_singleaccuracy_d(gmx_simd_double_t x, gmx_simd_dbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_invsqrt_singleaccuracy_d(x);
#else
    return gmx_simd_invsqrt_singleaccuracy_d(gmx_simd_blendv_d(x, gmx_simd_set1_d(1.0), m));
#endif
}

/*! \brief Calculate 1/sqrt(x) for two SIMD doubles, but single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_invsqrt_pair_singleaccuracy_r.
 *
 * \param x0  First set of arguments, x0 must be positive - no argument checking.
 * \param x1  Second set of arguments, x1 must be positive - no argument checking.
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 *  In particular for double precision we can sometimes calculate square root
 *  pairs slightly faster by using single precision until the very last step.
 */
static gmx_inline void gmx_simdcall
gmx_simd_invsqrt_pair_singleaccuracy_d(gmx_simd_double_t x0,    gmx_simd_double_t x1,
                                       gmx_simd_double_t *out0, gmx_simd_double_t *out1)
{
#if (defined GMX_SIMD_HAVE_FLOAT) && (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH) && (GMX_SIMD_RSQRT_BITS < 22)
    gmx_simd_float_t  xf  = gmx_simd_cvt_dd2f(x0, x1);
    gmx_simd_float_t  luf = gmx_simd_rsqrt_f(xf);
    gmx_simd_double_t lu0, lu1;
    /* Intermediate target is single - mantissa+1 bits */
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = gmx_simd_rsqrt_iter_f(luf, xf);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = gmx_simd_rsqrt_iter_f(luf, xf);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = gmx_simd_rsqrt_iter_f(luf, xf);
#endif
    gmx_simd_cvt_f2dd(luf, &lu0, &lu1);
    /* We now have single-precision accuracy values in lu0/lu1 */
    *out0 = lu0;
    *out1 = lu1;
#else
    *out0 = gmx_simd_invsqrt_singleaccuracy_d(x0);
    *out1 = gmx_simd_invsqrt_singleaccuracy_d(x1);
#endif
}


/*! \brief Calculate 1/x for SIMD double, but in single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_inv_singleaccuracy_r.
 *
 *  \param x Argument that must be nonzero. This routine does not check arguments.
 *  \return 1/x. Result is undefined if your argument was invalid.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_inv_singleaccuracy_d(gmx_simd_double_t x)
{
    gmx_simd_double_t lu = gmx_simd_rcp_d(x);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rcp_iter_d(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rcp_iter_d(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd_rcp_iter_d(lu, x);
#endif
    return lu;
}

/*! \brief 1/x for masked entries of SIMD double, single accuracy.
 *
 * \copydetails gmx_simd_inv_maskfpe_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_inv_maskfpe_singleaccuracy_d(gmx_simd_double_t x, gmx_simd_dbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_inv_singleaccuracy_d(x);
#else
    return gmx_simd_inv_singleaccuracy_d(gmx_simd_blendv_d(gmx_simd_set1_d(1.0), x, m));
#endif
}

/*! \brief 1/x for non-masked entries of SIMD double, single accuracy.
 *
 * \copydetails gmx_simd_inv_notmaskfpe_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_inv_notmaskfpe_singleaccuracy_d(gmx_simd_double_t x, gmx_simd_dbool_t gmx_unused m)
{
#ifdef NDEBUG
    return gmx_simd_inv_singleaccuracy_d(x);
#else
    return gmx_simd_inv_singleaccuracy_d(gmx_simd_blendv_d(x, gmx_simd_set1_d(1.0), m));
#endif
}

/*! \brief Calculate sqrt(x) (correct for 0.0) for SIMD double, single accuracy.
 *
 * You should normally call the real-precision routine \ref gmx_simd_sqrt_r.
 *
 *  \param x Argument that must be >=0.
 *  \return sqrt(x). If x=0, the result will correctly be set to 0.
 *          The result is undefined if the input value is negative.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_sqrt_singleaccuracy_d(gmx_simd_double_t x)
{
    gmx_simd_dbool_t   mask;
    gmx_simd_double_t  res;

    mask = gmx_simd_cmpeq_d(x, gmx_simd_setzero_d());
    res  = gmx_simd_blendnotzero_d(gmx_simd_invsqrt_notmaskfpe_singleaccuracy_d(x, mask), mask);
    return gmx_simd_mul_d(res, x);
}

/*! \brief SIMD log(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_log_singleaccuracy_r.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 */
#ifndef gmx_simd_log_singleaccuracy_d
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_log_singleaccuracy_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t  half       = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t  one        = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t  sqrt2      = gmx_simd_set1_d(sqrt(2.0));
    const gmx_simd_double_t  corr       = gmx_simd_set1_d(0.693147180559945286226764);
    const gmx_simd_double_t  CL9        = gmx_simd_set1_d(0.2371599674224853515625);
    const gmx_simd_double_t  CL7        = gmx_simd_set1_d(0.285279005765914916992188);
    const gmx_simd_double_t  CL5        = gmx_simd_set1_d(0.400005519390106201171875);
    const gmx_simd_double_t  CL3        = gmx_simd_set1_d(0.666666567325592041015625);
    const gmx_simd_double_t  CL1        = gmx_simd_set1_d(2.0);
    gmx_simd_double_t        fexp, x2, p;
    gmx_simd_dbool_t         mask;

    fexp  = gmx_simd_get_exponent_d(x);
    x     = gmx_simd_get_mantissa_d(x);

    mask  = gmx_simd_cmplt_d(sqrt2, x);
    /* Adjust to non-IEEE format for x>sqrt(2): exponent += 1, mantissa *= 0.5 */
    fexp  = gmx_simd_add_d(fexp, gmx_simd_blendzero_d(one, mask));
    x     = gmx_simd_mul_d(x, gmx_simd_blendv_d(one, half, mask));

    x     = gmx_simd_mul_d( gmx_simd_sub_d(x, one), gmx_simd_inv_singleaccuracy_d( gmx_simd_add_d(x, one) ) );
    x2    = gmx_simd_mul_d(x, x);

    p     = gmx_simd_fmadd_d(CL9, x2, CL7);
    p     = gmx_simd_fmadd_d(p, x2, CL5);
    p     = gmx_simd_fmadd_d(p, x2, CL3);
    p     = gmx_simd_fmadd_d(p, x2, CL1);
    p     = gmx_simd_fmadd_d(p, x, gmx_simd_mul_d(corr, fexp));

    return p;
}
#endif

#ifndef gmx_simd_exp2_singleaccuracy_d
/*! \brief SIMD 2^x. Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_exp2_singleaccuracy_r.
 *
 * \param x Argument.
 * \result 2^x. Undefined if input argument caused overflow.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_exp2_singleaccuracy_d(gmx_simd_double_t x)
{
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const gmx_simd_double_t  arglimit = gmx_simd_set1_d(126.0);
    const gmx_simd_double_t  CC6      = gmx_simd_set1_d(0.0001534581200287996416911311);
    const gmx_simd_double_t  CC5      = gmx_simd_set1_d(0.001339993121934088894618990);
    const gmx_simd_double_t  CC4      = gmx_simd_set1_d(0.009618488957115180159497841);
    const gmx_simd_double_t  CC3      = gmx_simd_set1_d(0.05550328776964726865751735);
    const gmx_simd_double_t  CC2      = gmx_simd_set1_d(0.2402264689063408646490722);
    const gmx_simd_double_t  CC1      = gmx_simd_set1_d(0.6931472057372680777553816);
    const gmx_simd_double_t  one      = gmx_simd_set1_d(1.0);

    gmx_simd_double_t        fexppart;
    gmx_simd_double_t        intpart;
    gmx_simd_double_t        p;
    gmx_simd_dbool_t         valuemask;

    fexppart  = gmx_simd_set_exponent_d(x);
    intpart   = gmx_simd_round_d(x);
    valuemask = gmx_simd_cmple_d(gmx_simd_fabs_d(x), arglimit);
    fexppart  = gmx_simd_blendzero_d(fexppart, valuemask);
    x         = gmx_simd_sub_d(x, intpart);

    p         = gmx_simd_fmadd_d(CC6, x, CC5);
    p         = gmx_simd_fmadd_d(p, x, CC4);
    p         = gmx_simd_fmadd_d(p, x, CC3);
    p         = gmx_simd_fmadd_d(p, x, CC2);
    p         = gmx_simd_fmadd_d(p, x, CC1);
    p         = gmx_simd_fmadd_d(p, x, one);
    x         = gmx_simd_mul_d(p, fexppart);
    return x;
}
#endif

#ifndef gmx_simd_exp_singleaccuracy_d
/*! \brief SIMD exp(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_exp_singleaccuracy_r.
 *
 * \param x Argument.
 * \result exp(x). Undefined if input argument caused overflow.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_exp_singleaccuracy_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t  argscale     = gmx_simd_set1_d(1.44269504088896341);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const gmx_simd_double_t  arglimit     = gmx_simd_set1_d(126.0);
    const gmx_simd_double_t  invargscale  = gmx_simd_set1_d(0.69314718055994528623);
    const gmx_simd_double_t  CC4          = gmx_simd_set1_d(0.00136324646882712841033936);
    const gmx_simd_double_t  CC3          = gmx_simd_set1_d(0.00836596917361021041870117);
    const gmx_simd_double_t  CC2          = gmx_simd_set1_d(0.0416710823774337768554688);
    const gmx_simd_double_t  CC1          = gmx_simd_set1_d(0.166665524244308471679688);
    const gmx_simd_double_t  CC0          = gmx_simd_set1_d(0.499999850988388061523438);
    const gmx_simd_double_t  one          = gmx_simd_set1_d(1.0);
    gmx_simd_double_t        fexppart;
    gmx_simd_double_t        intpart;
    gmx_simd_double_t        y, p;
    gmx_simd_dbool_t         valuemask;

    y         = gmx_simd_mul_d(x, argscale);
    fexppart  = gmx_simd_set_exponent_d(y);  /* rounds to nearest int internally */
    intpart   = gmx_simd_round_d(y);         /* use same rounding algorithm here */
    valuemask = gmx_simd_cmple_d(gmx_simd_fabs_d(y), arglimit);
    fexppart  = gmx_simd_blendzero_d(fexppart, valuemask);

    /* Extended precision arithmetics not needed since
     * we have double precision and only need single accuracy.
     */
    x         = gmx_simd_fnmadd_d(invargscale, intpart, x);

    p         = gmx_simd_fmadd_d(CC4, x, CC3);
    p         = gmx_simd_fmadd_d(p, x, CC2);
    p         = gmx_simd_fmadd_d(p, x, CC1);
    p         = gmx_simd_fmadd_d(p, x, CC0);
    p         = gmx_simd_fmadd_d(gmx_simd_mul_d(x, x), p, x);
    p         = gmx_simd_add_d(p, one);
    x         = gmx_simd_mul_d(p, fexppart);
    return x;
}
#endif

/*! \brief SIMD erf(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_erf_singleaccuracy_r.
 *
 * \param x The value to calculate erf(x) for.
 * \result erf(x)
 *
 * This routine achieves very close to single precision, but we do not care about
 * the last bit or the subnormal result range.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_erf_singleaccuracy_d(gmx_simd_double_t x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const gmx_simd_double_t  CA6      = gmx_simd_set1_d(7.853861353153693e-5);
    const gmx_simd_double_t  CA5      = gmx_simd_set1_d(-8.010193625184903e-4);
    const gmx_simd_double_t  CA4      = gmx_simd_set1_d(5.188327685732524e-3);
    const gmx_simd_double_t  CA3      = gmx_simd_set1_d(-2.685381193529856e-2);
    const gmx_simd_double_t  CA2      = gmx_simd_set1_d(1.128358514861418e-1);
    const gmx_simd_double_t  CA1      = gmx_simd_set1_d(-3.761262582423300e-1);
    const gmx_simd_double_t  CA0      = gmx_simd_set1_d(1.128379165726710);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const gmx_simd_double_t  CB9      = gmx_simd_set1_d(-0.0018629930017603923);
    const gmx_simd_double_t  CB8      = gmx_simd_set1_d(0.003909821287598495);
    const gmx_simd_double_t  CB7      = gmx_simd_set1_d(-0.0052094582210355615);
    const gmx_simd_double_t  CB6      = gmx_simd_set1_d(0.005685614362160572);
    const gmx_simd_double_t  CB5      = gmx_simd_set1_d(-0.0025367682853477272);
    const gmx_simd_double_t  CB4      = gmx_simd_set1_d(-0.010199799682318782);
    const gmx_simd_double_t  CB3      = gmx_simd_set1_d(0.04369575504816542);
    const gmx_simd_double_t  CB2      = gmx_simd_set1_d(-0.11884063474674492);
    const gmx_simd_double_t  CB1      = gmx_simd_set1_d(0.2732120154030589);
    const gmx_simd_double_t  CB0      = gmx_simd_set1_d(0.42758357702025784);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const gmx_simd_double_t  CC10     = gmx_simd_set1_d(-0.0445555913112064);
    const gmx_simd_double_t  CC9      = gmx_simd_set1_d(0.21376355144663348);
    const gmx_simd_double_t  CC8      = gmx_simd_set1_d(-0.3473187200259257);
    const gmx_simd_double_t  CC7      = gmx_simd_set1_d(0.016690861551248114);
    const gmx_simd_double_t  CC6      = gmx_simd_set1_d(0.7560973182491192);
    const gmx_simd_double_t  CC5      = gmx_simd_set1_d(-1.2137903600145787);
    const gmx_simd_double_t  CC4      = gmx_simd_set1_d(0.8411872321232948);
    const gmx_simd_double_t  CC3      = gmx_simd_set1_d(-0.08670413896296343);
    const gmx_simd_double_t  CC2      = gmx_simd_set1_d(-0.27124782687240334);
    const gmx_simd_double_t  CC1      = gmx_simd_set1_d(-0.0007502488047806069);
    const gmx_simd_double_t  CC0      = gmx_simd_set1_d(0.5642114853803148);
    const gmx_simd_double_t  one      = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t  two      = gmx_simd_set1_d(2.0);

    gmx_simd_double_t        x2, x4, y;
    gmx_simd_double_t        t, t2, w, w2;
    gmx_simd_double_t        pA0, pA1, pB0, pB1, pC0, pC1;
    gmx_simd_double_t        expmx2;
    gmx_simd_double_t        res_erf, res_erfc, res;
    gmx_simd_dbool_t         mask, msk_erf;

    /* Calculate erf() */
    x2   = gmx_simd_mul_d(x, x);
    x4   = gmx_simd_mul_d(x2, x2);

    pA0  = gmx_simd_fmadd_d(CA6, x4, CA4);
    pA1  = gmx_simd_fmadd_d(CA5, x4, CA3);
    pA0  = gmx_simd_fmadd_d(pA0, x4, CA2);
    pA1  = gmx_simd_fmadd_d(pA1, x4, CA1);
    pA0  = gmx_simd_mul_d(pA0, x4);
    pA0  = gmx_simd_fmadd_d(pA1, x2, pA0);
    /* Constant term must come last for precision reasons */
    pA0  = gmx_simd_add_d(pA0, CA0);

    res_erf = gmx_simd_mul_d(x, pA0);

    /* Calculate erfc */
    y       = gmx_simd_fabs_d(x);
    msk_erf = gmx_simd_cmplt_d(y, gmx_simd_set1_d(0.75));
    t       = gmx_simd_inv_notmaskfpe_singleaccuracy_d(y, msk_erf);
    w       = gmx_simd_sub_d(t, one);
    t2      = gmx_simd_mul_d(t, t);
    w2      = gmx_simd_mul_d(w, w);

    expmx2  = gmx_simd_exp_singleaccuracy_d( gmx_simd_fneg_d( gmx_simd_mul_d(y, y)));

    pB1  = gmx_simd_fmadd_d(CB9, w2, CB7);
    pB0  = gmx_simd_fmadd_d(CB8, w2, CB6);
    pB1  = gmx_simd_fmadd_d(pB1, w2, CB5);
    pB0  = gmx_simd_fmadd_d(pB0, w2, CB4);
    pB1  = gmx_simd_fmadd_d(pB1, w2, CB3);
    pB0  = gmx_simd_fmadd_d(pB0, w2, CB2);
    pB1  = gmx_simd_fmadd_d(pB1, w2, CB1);
    pB0  = gmx_simd_fmadd_d(pB0, w2, CB0);
    pB0  = gmx_simd_fmadd_d(pB1, w, pB0);

    pC0  = gmx_simd_fmadd_d(CC10, t2, CC8);
    pC1  = gmx_simd_fmadd_d(CC9, t2, CC7);
    pC0  = gmx_simd_fmadd_d(pC0, t2, CC6);
    pC1  = gmx_simd_fmadd_d(pC1, t2, CC5);
    pC0  = gmx_simd_fmadd_d(pC0, t2, CC4);
    pC1  = gmx_simd_fmadd_d(pC1, t2, CC3);
    pC0  = gmx_simd_fmadd_d(pC0, t2, CC2);
    pC1  = gmx_simd_fmadd_d(pC1, t2, CC1);

    pC0  = gmx_simd_fmadd_d(pC0, t2, CC0);
    pC0  = gmx_simd_fmadd_d(pC1, t, pC0);
    pC0  = gmx_simd_mul_d(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = gmx_simd_cmplt_d(two, y);
    res_erfc = gmx_simd_blendv_d(pB0, pC0, mask);
    res_erfc = gmx_simd_mul_d(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = gmx_simd_cmplt_d(x, gmx_simd_setzero_d());
    res_erfc = gmx_simd_blendv_d(res_erfc, gmx_simd_sub_d(two, res_erfc), mask);

    /* Select erf() or erfc() */
    mask = gmx_simd_cmplt_d(y, gmx_simd_set1_d(0.75));
    res  = gmx_simd_blendv_d(gmx_simd_sub_d(one, res_erfc), res_erf, mask);

    return res;
}

/*! \brief SIMD erfc(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_erfc_singleaccuracy_r.
 *
 * \param x The value to calculate erfc(x) for.
 * \result erfc(x)
 *
 * This routine achieves singleprecision (bar the last bit) over most of the
 * input range, but for large arguments where the result is getting close
 * to the minimum representable numbers we accept slightly larger errors
 * (think results that are in the ballpark of 10^-30) since that is not
 * relevant for MD.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_erfc_singleaccuracy_d(gmx_simd_double_t x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const gmx_simd_double_t  CA6      = gmx_simd_set1_d(7.853861353153693e-5);
    const gmx_simd_double_t  CA5      = gmx_simd_set1_d(-8.010193625184903e-4);
    const gmx_simd_double_t  CA4      = gmx_simd_set1_d(5.188327685732524e-3);
    const gmx_simd_double_t  CA3      = gmx_simd_set1_d(-2.685381193529856e-2);
    const gmx_simd_double_t  CA2      = gmx_simd_set1_d(1.128358514861418e-1);
    const gmx_simd_double_t  CA1      = gmx_simd_set1_d(-3.761262582423300e-1);
    const gmx_simd_double_t  CA0      = gmx_simd_set1_d(1.128379165726710);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const gmx_simd_double_t  CB9      = gmx_simd_set1_d(-0.0018629930017603923);
    const gmx_simd_double_t  CB8      = gmx_simd_set1_d(0.003909821287598495);
    const gmx_simd_double_t  CB7      = gmx_simd_set1_d(-0.0052094582210355615);
    const gmx_simd_double_t  CB6      = gmx_simd_set1_d(0.005685614362160572);
    const gmx_simd_double_t  CB5      = gmx_simd_set1_d(-0.0025367682853477272);
    const gmx_simd_double_t  CB4      = gmx_simd_set1_d(-0.010199799682318782);
    const gmx_simd_double_t  CB3      = gmx_simd_set1_d(0.04369575504816542);
    const gmx_simd_double_t  CB2      = gmx_simd_set1_d(-0.11884063474674492);
    const gmx_simd_double_t  CB1      = gmx_simd_set1_d(0.2732120154030589);
    const gmx_simd_double_t  CB0      = gmx_simd_set1_d(0.42758357702025784);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const gmx_simd_double_t  CC10     = gmx_simd_set1_d(-0.0445555913112064);
    const gmx_simd_double_t  CC9      = gmx_simd_set1_d(0.21376355144663348);
    const gmx_simd_double_t  CC8      = gmx_simd_set1_d(-0.3473187200259257);
    const gmx_simd_double_t  CC7      = gmx_simd_set1_d(0.016690861551248114);
    const gmx_simd_double_t  CC6      = gmx_simd_set1_d(0.7560973182491192);
    const gmx_simd_double_t  CC5      = gmx_simd_set1_d(-1.2137903600145787);
    const gmx_simd_double_t  CC4      = gmx_simd_set1_d(0.8411872321232948);
    const gmx_simd_double_t  CC3      = gmx_simd_set1_d(-0.08670413896296343);
    const gmx_simd_double_t  CC2      = gmx_simd_set1_d(-0.27124782687240334);
    const gmx_simd_double_t  CC1      = gmx_simd_set1_d(-0.0007502488047806069);
    const gmx_simd_double_t  CC0      = gmx_simd_set1_d(0.5642114853803148);
    const gmx_simd_double_t  one      = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t  two      = gmx_simd_set1_d(2.0);

    gmx_simd_double_t        x2, x4, y;
    gmx_simd_double_t        t, t2, w, w2;
    gmx_simd_double_t        pA0, pA1, pB0, pB1, pC0, pC1;
    gmx_simd_double_t        expmx2;
    gmx_simd_double_t        res_erf, res_erfc, res;
    gmx_simd_dbool_t         mask, msk_erf;

    /* Calculate erf() */
    x2     = gmx_simd_mul_d(x, x);
    x4     = gmx_simd_mul_d(x2, x2);

    pA0  = gmx_simd_fmadd_d(CA6, x4, CA4);
    pA1  = gmx_simd_fmadd_d(CA5, x4, CA3);
    pA0  = gmx_simd_fmadd_d(pA0, x4, CA2);
    pA1  = gmx_simd_fmadd_d(pA1, x4, CA1);
    pA1  = gmx_simd_mul_d(pA1, x2);
    pA0  = gmx_simd_fmadd_d(pA0, x4, pA1);
    /* Constant term must come last for precision reasons */
    pA0  = gmx_simd_add_d(pA0, CA0);

    res_erf = gmx_simd_mul_d(x, pA0);

    /* Calculate erfc */
    y       = gmx_simd_fabs_d(x);
    msk_erf = gmx_simd_cmplt_d(y, gmx_simd_set1_d(0.75));
    t       = gmx_simd_inv_notmaskfpe_singleaccuracy_d(y, msk_erf);
    w       = gmx_simd_sub_d(t, one);
    t2      = gmx_simd_mul_d(t, t);
    w2      = gmx_simd_mul_d(w, w);

    expmx2  = gmx_simd_exp_singleaccuracy_d( gmx_simd_fneg_d( gmx_simd_mul_d(y, y) ) );

    pB1  = gmx_simd_fmadd_d(CB9, w2, CB7);
    pB0  = gmx_simd_fmadd_d(CB8, w2, CB6);
    pB1  = gmx_simd_fmadd_d(pB1, w2, CB5);
    pB0  = gmx_simd_fmadd_d(pB0, w2, CB4);
    pB1  = gmx_simd_fmadd_d(pB1, w2, CB3);
    pB0  = gmx_simd_fmadd_d(pB0, w2, CB2);
    pB1  = gmx_simd_fmadd_d(pB1, w2, CB1);
    pB0  = gmx_simd_fmadd_d(pB0, w2, CB0);
    pB0  = gmx_simd_fmadd_d(pB1, w, pB0);

    pC0  = gmx_simd_fmadd_d(CC10, t2, CC8);
    pC1  = gmx_simd_fmadd_d(CC9, t2, CC7);
    pC0  = gmx_simd_fmadd_d(pC0, t2, CC6);
    pC1  = gmx_simd_fmadd_d(pC1, t2, CC5);
    pC0  = gmx_simd_fmadd_d(pC0, t2, CC4);
    pC1  = gmx_simd_fmadd_d(pC1, t2, CC3);
    pC0  = gmx_simd_fmadd_d(pC0, t2, CC2);
    pC1  = gmx_simd_fmadd_d(pC1, t2, CC1);

    pC0  = gmx_simd_fmadd_d(pC0, t2, CC0);
    pC0  = gmx_simd_fmadd_d(pC1, t, pC0);
    pC0  = gmx_simd_mul_d(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = gmx_simd_cmplt_d(two, y);
    res_erfc = gmx_simd_blendv_d(pB0, pC0, mask);
    res_erfc = gmx_simd_mul_d(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = gmx_simd_cmplt_d(x, gmx_simd_setzero_d());
    res_erfc = gmx_simd_blendv_d(res_erfc, gmx_simd_sub_d(two, res_erfc), mask);

    /* Select erf() or erfc() */
    mask = gmx_simd_cmplt_d(y, gmx_simd_set1_d(0.75));
    res  = gmx_simd_blendv_d(res_erfc, gmx_simd_sub_d(one, res_erf), mask);

    return res;
}

/*! \brief SIMD sin \& cos. Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_sincos_singleaccuracy_r.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 */
static gmx_inline void gmx_simdcall
gmx_simd_sincos_singleaccuracy_d(gmx_simd_double_t x, gmx_simd_double_t *sinval, gmx_simd_double_t *cosval)
{
    /* Constants to subtract Pi/4*x from y while minimizing precision loss */
    const gmx_simd_double_t  argred0         = gmx_simd_set1_d(2*0.78539816290140151978);
    const gmx_simd_double_t  argred1         = gmx_simd_set1_d(2*4.9604678871439933374e-10);
    const gmx_simd_double_t  argred2         = gmx_simd_set1_d(2*1.1258708853173288931e-18);
    const gmx_simd_double_t  two_over_pi     = gmx_simd_set1_d(2.0/M_PI);
    const gmx_simd_double_t  const_sin2      = gmx_simd_set1_d(-1.9515295891e-4);
    const gmx_simd_double_t  const_sin1      = gmx_simd_set1_d( 8.3321608736e-3);
    const gmx_simd_double_t  const_sin0      = gmx_simd_set1_d(-1.6666654611e-1);
    const gmx_simd_double_t  const_cos2      = gmx_simd_set1_d( 2.443315711809948e-5);
    const gmx_simd_double_t  const_cos1      = gmx_simd_set1_d(-1.388731625493765e-3);
    const gmx_simd_double_t  const_cos0      = gmx_simd_set1_d( 4.166664568298827e-2);

    const gmx_simd_double_t  half            = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t  one             = gmx_simd_set1_d(1.0);
    gmx_simd_double_t        ssign, csign;
    gmx_simd_double_t        x2, y, z, psin, pcos, sss, ccc;
    gmx_simd_dbool_t         mask;
#if (defined GMX_SIMD_HAVE_FINT32) && (defined GMX_SIMD_HAVE_FINT32_ARITHMETICS) && (defined GMX_SIMD_HAVE_LOGICAL)
    const gmx_simd_dint32_t  ione            = gmx_simd_set1_di(1);
    const gmx_simd_dint32_t  itwo            = gmx_simd_set1_di(2);
    gmx_simd_dint32_t        iy;

    z       = gmx_simd_mul_d(x, two_over_pi);
    iy      = gmx_simd_cvt_d2i(z);
    y       = gmx_simd_round_d(z);

    mask    = gmx_simd_cvt_dib2db(gmx_simd_cmpeq_di(gmx_simd_and_di(iy, ione), gmx_simd_setzero_di()));
    ssign   = gmx_simd_blendzero_d(gmx_simd_set1_d(-0.0), gmx_simd_cvt_dib2db(gmx_simd_cmpeq_di(gmx_simd_and_di(iy, itwo), itwo)));
    csign   = gmx_simd_blendzero_d(gmx_simd_set1_d(-0.0), gmx_simd_cvt_dib2db(gmx_simd_cmpeq_di(gmx_simd_and_di(gmx_simd_add_di(iy, ione), itwo), itwo)));
#else
    const gmx_simd_double_t  quarter         = gmx_simd_set1_d(0.25);
    const gmx_simd_double_t  minusquarter    = gmx_simd_set1_d(-0.25);
    gmx_simd_double_t        q;
    gmx_simd_dbool_t         m1, m2, m3;

    /* The most obvious way to find the arguments quadrant in the unit circle
     * to calculate the sign is to use integer arithmetic, but that is not
     * present in all SIMD implementations. As an alternative, we have devised a
     * pure floating-point algorithm that uses truncation for argument reduction
     * so that we get a new value 0<=q<1 over the unit circle, and then
     * do floating-point comparisons with fractions. This is likely to be
     * slightly slower (~10%) due to the longer latencies of floating-point, so
     * we only use it when integer SIMD arithmetic is not present.
     */
    ssign   = x;
    x       = gmx_simd_fabs_d(x);
    /* It is critical that half-way cases are rounded down */
    z       = gmx_simd_fmadd_d(x, two_over_pi, half);
    y       = gmx_simd_trunc_d(z);
    q       = gmx_simd_mul_d(z, quarter);
    q       = gmx_simd_sub_d(q, gmx_simd_trunc_d(q));
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q       = gmx_simd_sub_d(q, half);
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#    ifdef GMX_SIMD_HAVE_LOGICAL
    ssign   = gmx_simd_and_d(ssign, gmx_simd_set1_d(-0.0));
    csign   = gmx_simd_andnot_d(q, gmx_simd_set1_d(-0.0));
    ssign   = gmx_simd_xor_d(ssign, csign);
#    else
    csign   = gmx_simd_xor_sign_d(gmx_simd_set1_d(-1.0), q);

    ssign   = gmx_simd_xor_sign_d(ssign, csign);    /* swap ssign if csign was set. */
#    endif
    /* Check if y had value 1 or 3 (remember we subtracted 0.5 from q) */
    m1      = gmx_simd_cmplt_d(q, minusquarter);
    m2      = gmx_simd_cmple_d(gmx_simd_setzero_d(), q);
    m3      = gmx_simd_cmplt_d(q, quarter);
    m2      = gmx_simd_and_db(m2, m3);
    mask    = gmx_simd_or_db(m1, m2);
    /* where mask is FALSE, set sign. */
    csign   = gmx_simd_xor_sign_d(csign, gmx_simd_blendv_d(gmx_simd_set1_d(-1.0), one, mask));
#endif
    x       = gmx_simd_fnmadd_d(y, argred0, x);
    x       = gmx_simd_fnmadd_d(y, argred1, x);
    x       = gmx_simd_fnmadd_d(y, argred2, x);
    x2      = gmx_simd_mul_d(x, x);

    psin    = gmx_simd_fmadd_d(const_sin2, x2, const_sin1);
    psin    = gmx_simd_fmadd_d(psin, x2, const_sin0);
    psin    = gmx_simd_fmadd_d(psin, gmx_simd_mul_d(x, x2), x);
    pcos    = gmx_simd_fmadd_d(const_cos2, x2, const_cos1);
    pcos    = gmx_simd_fmadd_d(pcos, x2, const_cos0);
    pcos    = gmx_simd_fmsub_d(pcos, x2, half);
    pcos    = gmx_simd_fmadd_d(pcos, x2, one);

    sss     = gmx_simd_blendv_d(pcos, psin, mask);
    ccc     = gmx_simd_blendv_d(psin, pcos, mask);
    /* See comment for GMX_SIMD_HAVE_LOGICAL section above. */
#ifdef GMX_SIMD_HAVE_LOGICAL
    *sinval = gmx_simd_xor_d(sss, ssign);
    *cosval = gmx_simd_xor_d(ccc, csign);
#else
    *sinval = gmx_simd_xor_sign_d(sss, ssign);
    *cosval = gmx_simd_xor_sign_d(ccc, csign);
#endif
}

/*! \brief SIMD sin(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_sin_singleaccuracy_r.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref gmx_simd_sincos_r and waste a factor 2 in performance.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_sin_singleaccuracy_d(gmx_simd_double_t x)
{
    gmx_simd_double_t s, c;
    gmx_simd_sincos_singleaccuracy_d(x, &s, &c);
    return s;
}

/*! \brief SIMD cos(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_cos_singleaccuracy_r.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref gmx_simd_sincos_r and waste a factor 2 in performance.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_cos_singleaccuracy_d(gmx_simd_double_t x)
{
    gmx_simd_double_t s, c;
    gmx_simd_sincos_singleaccuracy_d(x, &s, &c);
    return c;
}

/*! \brief SIMD tan(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_tan_singleaccuracy_r.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_tan_singleaccuracy_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t  argred0         = gmx_simd_set1_d(2*0.78539816290140151978);
    const gmx_simd_double_t  argred1         = gmx_simd_set1_d(2*4.9604678871439933374e-10);
    const gmx_simd_double_t  argred2         = gmx_simd_set1_d(2*1.1258708853173288931e-18);
    const gmx_simd_double_t  two_over_pi     = gmx_simd_set1_d(2.0/M_PI);
    const gmx_simd_double_t  CT6             = gmx_simd_set1_d(0.009498288995810566122993911);
    const gmx_simd_double_t  CT5             = gmx_simd_set1_d(0.002895755790837379295226923);
    const gmx_simd_double_t  CT4             = gmx_simd_set1_d(0.02460087336161924491836265);
    const gmx_simd_double_t  CT3             = gmx_simd_set1_d(0.05334912882656359828045988);
    const gmx_simd_double_t  CT2             = gmx_simd_set1_d(0.1333989091464957704418495);
    const gmx_simd_double_t  CT1             = gmx_simd_set1_d(0.3333307599244198227797507);

    gmx_simd_double_t        x2, p, y, z;
    gmx_simd_dbool_t         mask;

#if (defined GMX_SIMD_HAVE_FINT32) && (defined GMX_SIMD_HAVE_FINT32_ARITHMETICS) && (defined GMX_SIMD_HAVE_LOGICAL)
    gmx_simd_dint32_t  iy;
    gmx_simd_dint32_t  ione = gmx_simd_set1_di(1);

    z       = gmx_simd_mul_d(x, two_over_pi);
    iy      = gmx_simd_cvt_d2i(z);
    y       = gmx_simd_round_d(z);
    mask    = gmx_simd_cvt_dib2db(gmx_simd_cmpeq_di(gmx_simd_and_di(iy, ione), ione));

    x       = gmx_simd_fnmadd_d(y, argred0, x);
    x       = gmx_simd_fnmadd_d(y, argred1, x);
    x       = gmx_simd_fnmadd_d(y, argred2, x);
    x       = gmx_simd_xor_d(gmx_simd_blendzero_d(gmx_simd_set1_d(-0.0), mask), x);
#else
    const gmx_simd_double_t  quarter         = gmx_simd_set1_d(0.25);
    const gmx_simd_double_t  half            = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t  threequarter    = gmx_simd_set1_d(0.75);
    gmx_simd_double_t        w, q;
    gmx_simd_dbool_t         m1, m2, m3;

    w       = gmx_simd_fabs_d(x);
    z       = gmx_simd_fmadd_d(w, two_over_pi, half);
    y       = gmx_simd_trunc_d(z);
    q       = gmx_simd_mul_d(z, quarter);
    q       = gmx_simd_sub_d(q, gmx_simd_trunc_d(q));
    m1      = gmx_simd_cmple_d(quarter, q);
    m2      = gmx_simd_cmplt_d(q, half);
    m3      = gmx_simd_cmple_d(threequarter, q);
    m1      = gmx_simd_and_db(m1, m2);
    mask    = gmx_simd_or_db(m1, m3);
    w       = gmx_simd_fnmadd_d(y, argred0, w);
    w       = gmx_simd_fnmadd_d(y, argred1, w);
    w       = gmx_simd_fnmadd_d(y, argred2, w);

    w       = gmx_simd_blendv_d(w, gmx_simd_fneg_d(w), mask);
    x       = gmx_simd_xor_sign_d(w, x);
#endif
    x2      = gmx_simd_mul_d(x, x);
    p       = gmx_simd_fmadd_d(CT6, x2, CT5);
    p       = gmx_simd_fmadd_d(p, x2, CT4);
    p       = gmx_simd_fmadd_d(p, x2, CT3);
    p       = gmx_simd_fmadd_d(p, x2, CT2);
    p       = gmx_simd_fmadd_d(p, x2, CT1);
    p       = gmx_simd_fmadd_d(x2, gmx_simd_mul_d(p, x), x);

    p       = gmx_simd_blendv_d( p, gmx_simd_inv_maskfpe_singleaccuracy_d(p, mask), mask);
    return p;
}

/*! \brief SIMD asin(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_asin_singleaccuracy_r.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_asin_singleaccuracy_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t limitlow   = gmx_simd_set1_d(1e-4);
    const gmx_simd_double_t half       = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t one        = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t halfpi     = gmx_simd_set1_d(M_PI/2.0);
    const gmx_simd_double_t CC5        = gmx_simd_set1_d(4.2163199048E-2);
    const gmx_simd_double_t CC4        = gmx_simd_set1_d(2.4181311049E-2);
    const gmx_simd_double_t CC3        = gmx_simd_set1_d(4.5470025998E-2);
    const gmx_simd_double_t CC2        = gmx_simd_set1_d(7.4953002686E-2);
    const gmx_simd_double_t CC1        = gmx_simd_set1_d(1.6666752422E-1);
    gmx_simd_double_t       xabs;
    gmx_simd_double_t       z, z1, z2, q, q1, q2;
    gmx_simd_double_t       pA, pB;
    gmx_simd_dbool_t        mask, mask2;

    xabs  = gmx_simd_fabs_d(x);
    mask  = gmx_simd_cmplt_d(half, xabs);
    z1    = gmx_simd_mul_d(half, gmx_simd_sub_d(one, xabs));
    mask2 = gmx_simd_cmpeq_d(xabs, one);
    q1    = gmx_simd_mul_d(z1, gmx_simd_invsqrt_notmaskfpe_singleaccuracy_d(z1, mask2));
    q1    = gmx_simd_blendnotzero_d(q1, gmx_simd_cmpeq_d(xabs, one));
    q2    = xabs;
    z2    = gmx_simd_mul_d(q2, q2);
    z     = gmx_simd_blendv_d(z2, z1, mask);
    q     = gmx_simd_blendv_d(q2, q1, mask);

    z2    = gmx_simd_mul_d(z, z);
    pA    = gmx_simd_fmadd_d(CC5, z2, CC3);
    pB    = gmx_simd_fmadd_d(CC4, z2, CC2);
    pA    = gmx_simd_fmadd_d(pA, z2, CC1);
    pA    = gmx_simd_mul_d(pA, z);
    z     = gmx_simd_fmadd_d(pB, z2, pA);
    z     = gmx_simd_fmadd_d(z, q, q);
    q2    = gmx_simd_sub_d(halfpi, z);
    q2    = gmx_simd_sub_d(q2, z);
    z     = gmx_simd_blendv_d(z, q2, mask);

    mask  = gmx_simd_cmplt_d(limitlow, xabs);
    z     = gmx_simd_blendv_d( xabs, z, mask );
    z     = gmx_simd_xor_sign_d(z, x);

    return z;
}

/*! \brief SIMD acos(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_acos_singleaccuracy_r.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_acos_singleaccuracy_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t one       = gmx_simd_set1_d(1.0);
    const gmx_simd_double_t half      = gmx_simd_set1_d(0.5);
    const gmx_simd_double_t pi        = gmx_simd_set1_d(M_PI);
    const gmx_simd_double_t halfpi    = gmx_simd_set1_d(M_PI/2.0);
    gmx_simd_double_t       xabs;
    gmx_simd_double_t       z, z1, z2, z3;
    gmx_simd_dbool_t        mask1, mask2, mask3;

    xabs  = gmx_simd_fabs_d(x);
    mask1 = gmx_simd_cmplt_d(half, xabs);
    mask2 = gmx_simd_cmplt_d(gmx_simd_setzero_d(), x);

    z     = gmx_simd_mul_d(half, gmx_simd_sub_d(one, xabs));
    mask3 = gmx_simd_cmpeq_d(xabs, one);
    z     = gmx_simd_mul_d(z, gmx_simd_invsqrt_notmaskfpe_singleaccuracy_d(z, mask3));
    z     = gmx_simd_blendnotzero_d(z, gmx_simd_cmpeq_d(xabs, one));
    z     = gmx_simd_blendv_d(x, z, mask1);
    z     = gmx_simd_asin_singleaccuracy_d(z);

    z2    = gmx_simd_add_d(z, z);
    z1    = gmx_simd_sub_d(pi, z2);
    z3    = gmx_simd_sub_d(halfpi, z);
    z     = gmx_simd_blendv_d(z1, z2, mask2);
    z     = gmx_simd_blendv_d(z3, z, mask1);

    return z;
}

/*! \brief SIMD asin(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_atan_singleaccuracy_r.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x), same argument/value range as standard math library.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_atan_singleaccuracy_d(gmx_simd_double_t x)
{
    const gmx_simd_double_t halfpi    = gmx_simd_set1_d(M_PI/2);
    const gmx_simd_double_t CA17      = gmx_simd_set1_d(0.002823638962581753730774);
    const gmx_simd_double_t CA15      = gmx_simd_set1_d(-0.01595690287649631500244);
    const gmx_simd_double_t CA13      = gmx_simd_set1_d(0.04250498861074447631836);
    const gmx_simd_double_t CA11      = gmx_simd_set1_d(-0.07489009201526641845703);
    const gmx_simd_double_t CA9       = gmx_simd_set1_d(0.1063479334115982055664);
    const gmx_simd_double_t CA7       = gmx_simd_set1_d(-0.1420273631811141967773);
    const gmx_simd_double_t CA5       = gmx_simd_set1_d(0.1999269574880599975585);
    const gmx_simd_double_t CA3       = gmx_simd_set1_d(-0.3333310186862945556640);
    gmx_simd_double_t       x2, x3, x4, pA, pB;
    gmx_simd_dbool_t        mask, mask2;

    mask  = gmx_simd_cmplt_d(x, gmx_simd_setzero_d());
    x     = gmx_simd_fabs_d(x);
    mask2 = gmx_simd_cmplt_d(gmx_simd_set1_d(1.0), x);
    x     = gmx_simd_blendv_d(x, gmx_simd_inv_maskfpe_singleaccuracy_d(x, mask2), mask2);

    x2    = gmx_simd_mul_d(x, x);
    x3    = gmx_simd_mul_d(x2, x);
    x4    = gmx_simd_mul_d(x2, x2);
    pA    = gmx_simd_fmadd_d(CA17, x4, CA13);
    pB    = gmx_simd_fmadd_d(CA15, x4, CA11);
    pA    = gmx_simd_fmadd_d(pA, x4, CA9);
    pB    = gmx_simd_fmadd_d(pB, x4, CA7);
    pA    = gmx_simd_fmadd_d(pA, x4, CA5);
    pB    = gmx_simd_fmadd_d(pB, x4, CA3);
    pA    = gmx_simd_fmadd_d(pA, x2, pB);
    pA    = gmx_simd_fmadd_d(pA, x3, x);

    pA    = gmx_simd_blendv_d(pA, gmx_simd_sub_d(halfpi, pA), mask2);
    pA    = gmx_simd_blendv_d(pA, gmx_simd_fneg_d(pA), mask);

    return pA;
}

/*! \brief SIMD atan2(y,x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_atan2_singleaccuracy_r.
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x), same argument/value range as standard math library.
 *
 * \note This routine should provide correct results for all finite
 * non-zero or positive-zero arguments. However, negative zero arguments will
 * be treated as positive zero, which means the return value will deviate from
 * the standard math library atan2(y,x) for those cases. That should not be
 * of any concern in Gromacs, and in particular it will not affect calculations
 * of angles from vectors.
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_atan2_singleaccuracy_d(gmx_simd_double_t y, gmx_simd_double_t x)
{
    const gmx_simd_double_t pi          = gmx_simd_set1_d(M_PI);
    const gmx_simd_double_t halfpi      = gmx_simd_set1_d(M_PI/2.0);
    gmx_simd_double_t       xinv, p, aoffset;
    gmx_simd_dbool_t        mask_x0, mask_y0, mask_xlt0, mask_ylt0;

    mask_x0   = gmx_simd_cmpeq_d(x, gmx_simd_setzero_d());
    mask_y0   = gmx_simd_cmpeq_d(y, gmx_simd_setzero_d());
    mask_xlt0 = gmx_simd_cmplt_d(x, gmx_simd_setzero_d());
    mask_ylt0 = gmx_simd_cmplt_d(y, gmx_simd_setzero_d());

    aoffset   = gmx_simd_blendzero_d(halfpi, mask_x0);
    aoffset   = gmx_simd_blendnotzero_d(aoffset, mask_y0);

    aoffset   = gmx_simd_blendv_d(aoffset, pi, mask_xlt0);
    aoffset   = gmx_simd_blendv_d(aoffset, gmx_simd_fneg_d(aoffset), mask_ylt0);

    xinv      = gmx_simd_blendnotzero_d(gmx_simd_inv_notmaskfpe_singleaccuracy_d(x, mask_x0), mask_x0);
    p         = gmx_simd_mul_d(y, xinv);
    p         = gmx_simd_atan_singleaccuracy_d(p);
    p         = gmx_simd_add_d(p, aoffset);

    return p;
}

/*! \brief Analytical PME force correction, double SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_pmecorrF_singleaccuracy_r.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb force - see below for details.
 *
 * This routine is meant to enable analytical evaluation of the
 * direct-space PME electrostatic force to avoid tables.
 *
 * The direct-space potential should be \f$ \mbox{erfc}(\beta r)/r\f$, but there
 * are some problems evaluating that:
 *
 * First, the error function is difficult (read: expensive) to
 * approxmiate accurately for intermediate to large arguments, and
 * this happens already in ranges of \f$(\beta r)\f$ that occur in simulations.
 * Second, we now try to avoid calculating potentials in Gromacs but
 * use forces directly.
 *
 * We can simply things slight by noting that the PME part is really
 * a correction to the normal Coulomb force since \f$\mbox{erfc}(z)=1-\mbox{erf}(z)\f$, i.e.
 * \f[
 * V = \frac{1}{r} - \frac{\mbox{erf}(\beta r)}{r}
 * \f]
 * The first term we already have from the inverse square root, so
 * that we can leave out of this routine.
 *
 * For pme tolerances of 1e-3 to 1e-8 and cutoffs of 0.5nm to 1.8nm,
 * the argument \f$beta r\f$ will be in the range 0.15 to ~4. Use your
 * favorite plotting program to realize how well-behaved \f$\frac{\mbox{erf}(z)}{z}\f$ is
 * in this range!
 *
 * We approximate \f$f(z)=\mbox{erf}(z)/z\f$ with a rational minimax polynomial.
 * However, it turns out it is more efficient to approximate \f$f(z)/z\f$ and
 * then only use even powers. This is another minor optimization, since
 * we actually \a want \f$f(z)/z\f$, because it is going to be multiplied by
 * the vector between the two atoms to get the vectorial force. The
 * fastest flops are the ones we can avoid calculating!
 *
 * So, here's how it should be used:
 *
 * 1. Calculate \f$r^2\f$.
 * 2. Multiply by \f$\beta^2\f$, so you get \f$z^2=(\beta r)^2\f$.
 * 3. Evaluate this routine with \f$z^2\f$ as the argument.
 * 4. The return value is the expression:
 *
 * \f[
 *    \frac{2 \exp{-z^2}}{\sqrt{\pi} z^2}-\frac{\mbox{erf}(z)}{z^3}
 * \f]
 *
 * 5. Multiply the entire expression by \f$\beta^3\f$. This will get you
 *
 *  \f[
 *    \frac{2 \beta^3 \exp(-z^2)}{\sqrt{\pi} z^2} - \frac{\beta^3 \mbox{erf}(z)}{z^3}
 *  \f]
 *
 *    or, switching back to \f$r\f$ (since \f$z=r \beta\f$):
 *
 *  \f[
 *    \frac{2 \beta \exp(-r^2 \beta^2)}{\sqrt{\pi} r^2} - \frac{\mbox{erf}(r \beta)}{r^3}
 *  \f]
 *
 *    With a bit of math exercise you should be able to confirm that
 *    this is exactly
 *
 *  \f[
 *   \frac{\frac{d}{dr}\left( \frac{\mbox{erf}(\beta r)}{r} \right)}{r}
 *  \f]
 *
 * 6. Add the result to \f$r^{-3}\f$, multiply by the product of the charges,
 *    and you have your force (divided by \f$r\f$). A final multiplication
 *    with the vector connecting the two particles and you have your
 *    vectorial force to add to the particles.
 *
 * This approximation achieves an accuracy slightly lower than 1e-6; when
 * added to \f$1/r\f$ the error will be insignificant.
 *
 */
static gmx_simd_double_t gmx_simdcall
gmx_simd_pmecorrF_singleaccuracy_d(gmx_simd_double_t z2)
{
    const gmx_simd_double_t  FN6      = gmx_simd_set1_d(-1.7357322914161492954e-8);
    const gmx_simd_double_t  FN5      = gmx_simd_set1_d(1.4703624142580877519e-6);
    const gmx_simd_double_t  FN4      = gmx_simd_set1_d(-0.000053401640219807709149);
    const gmx_simd_double_t  FN3      = gmx_simd_set1_d(0.0010054721316683106153);
    const gmx_simd_double_t  FN2      = gmx_simd_set1_d(-0.019278317264888380590);
    const gmx_simd_double_t  FN1      = gmx_simd_set1_d(0.069670166153766424023);
    const gmx_simd_double_t  FN0      = gmx_simd_set1_d(-0.75225204789749321333);

    const gmx_simd_double_t  FD4      = gmx_simd_set1_d(0.0011193462567257629232);
    const gmx_simd_double_t  FD3      = gmx_simd_set1_d(0.014866955030185295499);
    const gmx_simd_double_t  FD2      = gmx_simd_set1_d(0.11583842382862377919);
    const gmx_simd_double_t  FD1      = gmx_simd_set1_d(0.50736591960530292870);
    const gmx_simd_double_t  FD0      = gmx_simd_set1_d(1.0);

    gmx_simd_double_t        z4;
    gmx_simd_double_t        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = gmx_simd_mul_d(z2, z2);

    polyFD0        = gmx_simd_fmadd_d(FD4, z4, FD2);
    polyFD1        = gmx_simd_fmadd_d(FD3, z4, FD1);
    polyFD0        = gmx_simd_fmadd_d(polyFD0, z4, FD0);
    polyFD0        = gmx_simd_fmadd_d(polyFD1, z2, polyFD0);

    polyFD0        = gmx_simd_inv_singleaccuracy_d(polyFD0);

    polyFN0        = gmx_simd_fmadd_d(FN6, z4, FN4);
    polyFN1        = gmx_simd_fmadd_d(FN5, z4, FN3);
    polyFN0        = gmx_simd_fmadd_d(polyFN0, z4, FN2);
    polyFN1        = gmx_simd_fmadd_d(polyFN1, z4, FN1);
    polyFN0        = gmx_simd_fmadd_d(polyFN0, z4, FN0);
    polyFN0        = gmx_simd_fmadd_d(polyFN1, z2, polyFN0);

    return gmx_simd_mul_d(polyFN0, polyFD0);
}



/*! \brief Analytical PME potential correction, double SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref gmx_simd_pmecorrV_singleaccuracy_r.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb potential - see below for details.
 *
 * See \ref gmx_simd_pmecorrF_f for details about the approximation.
 *
 * This routine calculates \f$\mbox{erf}(z)/z\f$, although you should provide \f$z^2\f$
 * as the input argument.
 *
 * Here's how it should be used:
 *
 * 1. Calculate \f$r^2\f$.
 * 2. Multiply by \f$\beta^2\f$, so you get \f$z^2=\beta^2*r^2\f$.
 * 3. Evaluate this routine with z^2 as the argument.
 * 4. The return value is the expression:
 *
 *  \f[
 *   \frac{\mbox{erf}(z)}{z}
 *  \f]
 *
 * 5. Multiply the entire expression by beta and switching back to \f$r\f$ (since \f$z=r \beta\f$):
 *
 *  \f[
 *    \frac{\mbox{erf}(r \beta)}{r}
 *  \f]
 *
 * 6. Subtract the result from \f$1/r\f$, multiply by the product of the charges,
 *    and you have your potential.
 *
 * This approximation achieves an accuracy slightly lower than 1e-6; when
 * added to \f$1/r\f$ the error will be insignificant.
 */
static gmx_simd_double_t gmx_simdcall
gmx_simd_pmecorrV_singleaccuracy_d(gmx_simd_double_t z2)
{
    const gmx_simd_double_t  VN6      = gmx_simd_set1_d(1.9296833005951166339e-8);
    const gmx_simd_double_t  VN5      = gmx_simd_set1_d(-1.4213390571557850962e-6);
    const gmx_simd_double_t  VN4      = gmx_simd_set1_d(0.000041603292906656984871);
    const gmx_simd_double_t  VN3      = gmx_simd_set1_d(-0.00013134036773265025626);
    const gmx_simd_double_t  VN2      = gmx_simd_set1_d(0.038657983986041781264);
    const gmx_simd_double_t  VN1      = gmx_simd_set1_d(0.11285044772717598220);
    const gmx_simd_double_t  VN0      = gmx_simd_set1_d(1.1283802385263030286);

    const gmx_simd_double_t  VD3      = gmx_simd_set1_d(0.0066752224023576045451);
    const gmx_simd_double_t  VD2      = gmx_simd_set1_d(0.078647795836373922256);
    const gmx_simd_double_t  VD1      = gmx_simd_set1_d(0.43336185284710920150);
    const gmx_simd_double_t  VD0      = gmx_simd_set1_d(1.0);

    gmx_simd_double_t        z4;
    gmx_simd_double_t        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = gmx_simd_mul_d(z2, z2);

    polyVD1        = gmx_simd_fmadd_d(VD3, z4, VD1);
    polyVD0        = gmx_simd_fmadd_d(VD2, z4, VD0);
    polyVD0        = gmx_simd_fmadd_d(polyVD1, z2, polyVD0);

    polyVD0        = gmx_simd_inv_singleaccuracy_d(polyVD0);

    polyVN0        = gmx_simd_fmadd_d(VN6, z4, VN4);
    polyVN1        = gmx_simd_fmadd_d(VN5, z4, VN3);
    polyVN0        = gmx_simd_fmadd_d(polyVN0, z4, VN2);
    polyVN1        = gmx_simd_fmadd_d(polyVN1, z4, VN1);
    polyVN0        = gmx_simd_fmadd_d(polyVN0, z4, VN0);
    polyVN0        = gmx_simd_fmadd_d(polyVN1, z2, polyVN0);

    return gmx_simd_mul_d(polyVN0, polyVD0);
}

#endif


/*! \name SIMD4 math functions
 *
 * \note Only a subset of the math functions are implemented for SIMD4.
 *  \{
 */


#ifdef GMX_SIMD4_HAVE_FLOAT

/*************************************************************************
 * SINGLE PRECISION SIMD4 MATH FUNCTIONS - JUST A SMALL SUBSET SUPPORTED *
 *************************************************************************/

/*! \brief SIMD4 utility function to sum a+b+c+d for SIMD4 floats.
 *
 * \copydetails gmx_simd_sum4_f
 */
static gmx_inline gmx_simd4_float_t gmx_simdcall
gmx_simd4_sum4_f(gmx_simd4_float_t a, gmx_simd4_float_t b,
                 gmx_simd4_float_t c, gmx_simd4_float_t d)
{
    return gmx_simd4_add_f(gmx_simd4_add_f(a, b), gmx_simd4_add_f(c, d));
}

/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD4 float.
 *
 * \copydetails gmx_simd_rsqrt_iter_f
 */
static gmx_inline gmx_simd4_float_t gmx_simdcall
gmx_simd4_rsqrt_iter_f(gmx_simd4_float_t lu, gmx_simd4_float_t x)
{
#    ifdef GMX_SIMD_HAVE_FMA
    return gmx_simd4_fmadd_f(gmx_simd4_fnmadd_f(x, gmx_simd4_mul_f(lu, lu), gmx_simd4_set1_f(1.0f)), gmx_simd4_mul_f(lu, gmx_simd4_set1_f(0.5f)), lu);
#    else
    return gmx_simd4_mul_f(gmx_simd4_set1_f(0.5f), gmx_simd4_mul_f(gmx_simd4_sub_f(gmx_simd4_set1_f(3.0f), gmx_simd4_mul_f(gmx_simd4_mul_f(lu, lu), x)), lu));
#    endif
}

/*! \brief Calculate 1/sqrt(x) for SIMD4 float.
 *
 * \copydetails gmx_simd_invsqrt_f
 */
static gmx_inline gmx_simd4_float_t gmx_simdcall
gmx_simd4_invsqrt_f(gmx_simd4_float_t x)
{
    gmx_simd4_float_t lu = gmx_simd4_rsqrt_f(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd4_rsqrt_iter_f(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd4_rsqrt_iter_f(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd4_rsqrt_iter_f(lu, x);
#endif
    return lu;
}

#endif /* GMX_SIMD4_HAVE_FLOAT */



#ifdef GMX_SIMD4_HAVE_DOUBLE
/*************************************************************************
 * DOUBLE PRECISION SIMD4 MATH FUNCTIONS - JUST A SMALL SUBSET SUPPORTED *
 *************************************************************************/


/*! \brief SIMD4 utility function to sum a+b+c+d for SIMD4 doubles.
 *
 * \copydetails gmx_simd_sum4_f
 */
static gmx_inline gmx_simd4_double_t gmx_simdcall
gmx_simd4_sum4_d(gmx_simd4_double_t a, gmx_simd4_double_t b,
                 gmx_simd4_double_t c, gmx_simd4_double_t d)
{
    return gmx_simd4_add_d(gmx_simd4_add_d(a, b), gmx_simd4_add_d(c, d));
}

/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD4 double.
 *
 * \copydetails gmx_simd_rsqrt_iter_f
 */
static gmx_inline gmx_simd4_double_t gmx_simdcall
gmx_simd4_rsqrt_iter_d(gmx_simd4_double_t lu, gmx_simd4_double_t x)
{
#ifdef GMX_SIMD_HAVE_FMA
    return gmx_simd4_fmadd_d(gmx_simd4_fnmadd_d(x, gmx_simd4_mul_d(lu, lu), gmx_simd4_set1_d(1.0)), gmx_simd4_mul_d(lu, gmx_simd4_set1_d(0.5)), lu);
#else
    return gmx_simd4_mul_d(gmx_simd4_set1_d(0.5), gmx_simd4_mul_d(gmx_simd4_sub_d(gmx_simd4_set1_d(3.0), gmx_simd4_mul_d(gmx_simd4_mul_d(lu, lu), x)), lu));
#endif
}

/*! \brief Calculate 1/sqrt(x) for SIMD4 double.
 *
 * \copydetails gmx_simd_invsqrt_f
 */
static gmx_inline gmx_simd4_double_t gmx_simdcall
gmx_simd4_invsqrt_d(gmx_simd4_double_t x)
{
    gmx_simd4_double_t lu = gmx_simd4_rsqrt_d(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd4_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd4_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd4_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = gmx_simd4_rsqrt_iter_d(lu, x);
#endif
    return lu;
}

/**********************************************************************
 * SIMD4 MATH FUNCTIONS WITH DOUBLE PREC. DATA, SINGLE PREC. ACCURACY *
 **********************************************************************/

/*! \brief Calculate 1/sqrt(x) for SIMD4 double, but in single accuracy.
 *
 * \copydetails gmx_simd_invsqrt_singleaccuracy_d
 */
static gmx_inline gmx_simd4_double_t gmx_simdcall
gmx_simd4_invsqrt_singleaccuracy_d(gmx_simd4_double_t x)
{
    gmx_simd4_double_t lu = gmx_simd4_rsqrt_d(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd4_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd4_rsqrt_iter_d(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = gmx_simd4_rsqrt_iter_d(lu, x);
#endif
    return lu;
}


#endif /* GMX_SIMD4_HAVE_DOUBLE */

/*! \} */


/* Set defines based on default Gromacs precision */
#ifdef GMX_DOUBLE
/* Documentation in single branch below */

#    define gmx_simd_sum4_r           gmx_simd_sum4_d
#    define gmx_simd_xor_sign_r       gmx_simd_xor_sign_d
#    define gmx_simd4_sum4_r          gmx_simd4_sum4_d

/* On hardware that only supports double precision SIMD it is possible to use
 * the faster _singleaccuracy_d routines everywhere by setting the requested SIMD
 * accuracy to single precision.
 */
#if (GMX_SIMD_ACCURACY_BITS_DOUBLE > GMX_SIMD_ACCURACY_BITS_SINGLE)

#        define gmx_simd_invsqrt_r        gmx_simd_invsqrt_d
#        define gmx_simd_invsqrt_pair_r   gmx_simd_invsqrt_pair_d
#        define gmx_simd_sqrt_r           gmx_simd_sqrt_d
#        define gmx_simd_inv_r            gmx_simd_inv_d
#        define gmx_simd_log_r            gmx_simd_log_d
#        define gmx_simd_exp2_r           gmx_simd_exp2_d
#        define gmx_simd_exp_r            gmx_simd_exp_d
#        define gmx_simd_erf_r            gmx_simd_erf_d
#        define gmx_simd_erfc_r           gmx_simd_erfc_d
#        define gmx_simd_sincos_r         gmx_simd_sincos_d
#        define gmx_simd_sin_r            gmx_simd_sin_d
#        define gmx_simd_cos_r            gmx_simd_cos_d
#        define gmx_simd_tan_r            gmx_simd_tan_d
#        define gmx_simd_asin_r           gmx_simd_asin_d
#        define gmx_simd_acos_r           gmx_simd_acos_d
#        define gmx_simd_atan_r           gmx_simd_atan_d
#        define gmx_simd_atan2_r          gmx_simd_atan2_d
#        define gmx_simd_pmecorrF_r       gmx_simd_pmecorrF_d
#        define gmx_simd_pmecorrV_r       gmx_simd_pmecorrV_d
#        define gmx_simd4_invsqrt_r       gmx_simd4_invsqrt_d

#else

#        define gmx_simd_invsqrt_r        gmx_simd_invsqrt_singleaccuracy_d
#        define gmx_simd_invsqrt_pair_r   gmx_simd_invsqrt_pair_singleaccuracy_d
#        define gmx_simd_sqrt_r           gmx_simd_sqrt_singleaccuracy_d
#        define gmx_simd_inv_r            gmx_simd_inv_singleaccuracy_d
#        define gmx_simd_log_r            gmx_simd_log_singleaccuracy_d
#        define gmx_simd_exp2_r           gmx_simd_exp2_singleaccuracy_d
#        define gmx_simd_exp_r            gmx_simd_exp_singleaccuracy_d
#        define gmx_simd_erf_r            gmx_simd_erf_singleaccuracy_d
#        define gmx_simd_erfc_r           gmx_simd_erfc_singleaccuracy_d
#        define gmx_simd_sincos_r         gmx_simd_sincos_singleaccuracy_d
#        define gmx_simd_sin_r            gmx_simd_sin_singleaccuracy_d
#        define gmx_simd_cos_r            gmx_simd_cos_singleaccuracy_d
#        define gmx_simd_tan_r            gmx_simd_tan_singleaccuracy_d
#        define gmx_simd_asin_r           gmx_simd_asin_singleaccuracy_d
#        define gmx_simd_acos_r           gmx_simd_acos_singleaccuracy_d
#        define gmx_simd_atan_r           gmx_simd_atan_singleaccuracy_d
#        define gmx_simd_atan2_r          gmx_simd_atan2_singleaccuracy_d
#        define gmx_simd_pmecorrF_r       gmx_simd_pmecorrF_singleaccuracy_d
#        define gmx_simd_pmecorrV_r       gmx_simd_pmecorrV_singleaccuracy_d
#        define gmx_simd4_invsqrt_r       gmx_simd4_invsqrt_singleaccuracy_d

#endif

#    define gmx_simd_invsqrt_singleaccuracy_r        gmx_simd_invsqrt_singleaccuracy_d
#    define gmx_simd_invsqrt_pair_singleaccuracy_r   gmx_simd_invsqrt_pair_singleaccuracy_d
#    define gmx_simd_sqrt_singleaccuracy_r           gmx_simd_sqrt_singleaccuracy_d
#    define gmx_simd_inv_singleaccuracy_r            gmx_simd_inv_singleaccuracy_d
#    define gmx_simd_log_singleaccuracy_r            gmx_simd_log_singleaccuracy_d
#    define gmx_simd_exp2_singleaccuracy_r           gmx_simd_exp2_singleaccuracy_d
#    define gmx_simd_exp_singleaccuracy_r            gmx_simd_exp_singleaccuracy_d
#    define gmx_simd_erf_singleaccuracy_r            gmx_simd_erf_singleaccuracy_d
#    define gmx_simd_erfc_singleaccuracy_r           gmx_simd_erfc_singleaccuracy_d
#    define gmx_simd_sincos_singleaccuracy_r         gmx_simd_sincos_singleaccuracy_d
#    define gmx_simd_sin_singleaccuracy_r            gmx_simd_sin_singleaccuracy_d
#    define gmx_simd_cos_singleaccuracy_r            gmx_simd_cos_singleaccuracy_d
#    define gmx_simd_tan_singleaccuracy_r            gmx_simd_tan_singleaccuracy_d
#    define gmx_simd_asin_singleaccuracy_r           gmx_simd_asin_singleaccuracy_d
#    define gmx_simd_acos_singleaccuracy_r           gmx_simd_acos_singleaccuracy_d
#    define gmx_simd_atan_singleaccuracy_r           gmx_simd_atan_singleaccuracy_d
#    define gmx_simd_atan2_singleaccuracy_r          gmx_simd_atan2_singleaccuracy_d
#    define gmx_simd_pmecorrF_singleaccuracy_r       gmx_simd_pmecorrF_singleaccuracy_d
#    define gmx_simd_pmecorrV_singleaccuracy_r       gmx_simd_pmecorrV_singleaccuracy_d
#    define gmx_simd4_invsqrt_singleaccuracy_r       gmx_simd4_invsqrt_singleaccuracy_d

#else /* GMX_DOUBLE */

/*! \name Real-precision SIMD math functions
 *
 *  These are the ones you should typically call in Gromacs.
 * \{
 */

/*! \brief SIMD utility function to sum a+b+c+d for SIMD reals.
 *
 * \copydetails gmx_simd_sum4_f
 */
#    define gmx_simd_sum4_r           gmx_simd_sum4_f

/*! \brief Return -a if b is negative, SIMD real.
 *
 * \copydetails gmx_simd_xor_sign_f
 */
#    define gmx_simd_xor_sign_r       gmx_simd_xor_sign_f

/*! \brief Calculate 1/sqrt(x) for SIMD real.
 *
 * \copydetails gmx_simd_invsqrt_f
 */
#    define gmx_simd_invsqrt_r        gmx_simd_invsqrt_f

/*! \brief Calculate 1/sqrt(x) for two SIMD reals.
 *
 * \copydetails gmx_simd_invsqrt_pair_f
 */
#    define gmx_simd_invsqrt_pair_r   gmx_simd_invsqrt_pair_f

/*! \brief Calculate sqrt(x) correctly for SIMD real, including argument 0.0.
 *
 * \copydetails gmx_simd_sqrt_f
 */
#    define gmx_simd_sqrt_r           gmx_simd_sqrt_f

/*! \brief Calculate 1/x for SIMD real.
 *
 * \copydetails gmx_simd_inv_f
 */
#    define gmx_simd_inv_r            gmx_simd_inv_f

/*! \brief SIMD real log(x). This is the natural logarithm.
 *
 * \copydetails gmx_simd_log_f
 */
#    define gmx_simd_log_r            gmx_simd_log_f

/*! \brief SIMD real 2^x.
 *
 * \copydetails gmx_simd_exp2_f
 */
#    define gmx_simd_exp2_r           gmx_simd_exp2_f

/*! \brief SIMD real e^x.
 *
 * \copydetails gmx_simd_exp_f
 */
#    define gmx_simd_exp_r            gmx_simd_exp_f

/*! \brief SIMD real erf(x).
 *
 * \copydetails gmx_simd_erf_f
 */
#    define gmx_simd_erf_r            gmx_simd_erf_f

/*! \brief SIMD real erfc(x).
 *
 * \copydetails gmx_simd_erfc_f
 */
#    define gmx_simd_erfc_r           gmx_simd_erfc_f

/*! \brief SIMD real sin \& cos.
 *
 * \copydetails gmx_simd_sincos_f
 */
#    define gmx_simd_sincos_r         gmx_simd_sincos_f

/*! \brief SIMD real sin(x).
 *
 * \copydetails gmx_simd_sin_f
 */
#    define gmx_simd_sin_r            gmx_simd_sin_f

/*! \brief SIMD real cos(x).
 *
 * \copydetails gmx_simd_cos_f
 */
#    define gmx_simd_cos_r            gmx_simd_cos_f

/*! \brief SIMD real tan(x).
 *
 * \copydetails gmx_simd_tan_f
 */
#    define gmx_simd_tan_r            gmx_simd_tan_f

/*! \brief SIMD real asin(x).
 *
 * \copydetails gmx_simd_asin_f
 */
#    define gmx_simd_asin_r           gmx_simd_asin_f

/*! \brief SIMD real acos(x).
 *
 * \copydetails gmx_simd_acos_f
 */
#    define gmx_simd_acos_r           gmx_simd_acos_f

/*! \brief SIMD real atan(x).
 *
 * \copydetails gmx_simd_atan_f
 */
#    define gmx_simd_atan_r           gmx_simd_atan_f

/*! \brief SIMD real atan2(y,x).
 *
 * \copydetails gmx_simd_atan2_f
 */
#    define gmx_simd_atan2_r          gmx_simd_atan2_f

/*! \brief SIMD Analytic PME force correction.
 *
 * \copydetails gmx_simd_pmecorrF_f
 */
#    define gmx_simd_pmecorrF_r       gmx_simd_pmecorrF_f

/*! \brief SIMD Analytic PME potential correction.
 *
 * \copydetails gmx_simd_pmecorrV_f
 */
#    define gmx_simd_pmecorrV_r       gmx_simd_pmecorrV_f

/*! \brief Calculate 1/sqrt(x) for SIMD, only targeting single accuracy.
 *
 * \copydetails gmx_simd_invsqrt_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_invsqrt_singleaccuracy_r        gmx_simd_invsqrt_f

/*! \brief Calculate 1/sqrt(x) for SIMD pair, only targeting single accuracy.
 *
 * \copydetails gmx_simd_invsqrt_pair_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_invsqrt_pair_singleaccuracy_r   gmx_simd_invsqrt_pair_f

/*! \brief Calculate sqrt(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_sqrt_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_sqrt_singleaccuracy_r           gmx_simd_sqrt_f

/*! \brief Calculate 1/x for SIMD real, only targeting single accuracy.
 *
 * \copydetails gmx_simd_inv_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_inv_singleaccuracy_r            gmx_simd_inv_f

/*! \brief SIMD real log(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_log_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_log_singleaccuracy_r            gmx_simd_log_f

/*! \brief SIMD real 2^x, only targeting single accuracy.
 *
 * \copydetails gmx_simd_exp2_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_exp2_singleaccuracy_r           gmx_simd_exp2_f

/*! \brief SIMD real e^x, only targeting single accuracy.
 *
 * \copydetails gmx_simd_exp_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_exp_singleaccuracy_r            gmx_simd_exp_f

/*! \brief SIMD real erf(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_erf_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_erf_singleaccuracy_r            gmx_simd_erf_f

/*! \brief SIMD real erfc(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_erfc_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_erfc_singleaccuracy_r           gmx_simd_erfc_f

/*! \brief SIMD real sin \& cos, only targeting single accuracy.
 *
 * \copydetails gmx_simd_sincos_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_sincos_singleaccuracy_r         gmx_simd_sincos_f

/*! \brief SIMD real sin(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_sin_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_sin_singleaccuracy_r            gmx_simd_sin_f

/*! \brief SIMD real cos(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_cos_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_cos_singleaccuracy_r            gmx_simd_cos_f

/*! \brief SIMD real tan(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_tan_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_tan_singleaccuracy_r            gmx_simd_tan_f

/*! \brief SIMD real asin(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_asin_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_asin_singleaccuracy_r           gmx_simd_asin_f

/*! \brief SIMD real acos(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_acos_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_acos_singleaccuracy_r           gmx_simd_acos_f

/*! \brief SIMD real atan(x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_atan_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_atan_singleaccuracy_r           gmx_simd_atan_f

/*! \brief SIMD real atan2(y,x), only targeting single accuracy.
 *
 * \copydetails gmx_simd_atan2_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_atan2_singleaccuracy_r          gmx_simd_atan2_f

/*! \brief SIMD Analytic PME force corr., only targeting single accuracy.
 *
 * \copydetails gmx_simd_pmecorrF_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_pmecorrF_singleaccuracy_r       gmx_simd_pmecorrF_f

/*! \brief SIMD Analytic PME potential corr., only targeting single accuracy.
 *
 * \copydetails gmx_simd_pmecorrV_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd_pmecorrV_singleaccuracy_r       gmx_simd_pmecorrV_f


/*! \}
 * \name SIMD4 math functions
 * \{
 */

/*! \brief SIMD4 utility function to sum a+b+c+d for SIMD4 reals.
 *
 * \copydetails gmx_simd_sum4_f
 */
#    define gmx_simd4_sum4_r          gmx_simd4_sum4_f

/*! \brief Calculate 1/sqrt(x) for SIMD4 real.
 *
 * \copydetails gmx_simd_invsqrt_f
 */
#    define gmx_simd4_invsqrt_r       gmx_simd4_invsqrt_f

/*! \brief 1/sqrt(x) for SIMD4 real. Single accuracy, even for double prec.
 *
 * \copydetails gmx_simd4_invsqrt_r
 *
 * \note This is a performance-targeted function that only achieves single
 *       precision accuracy, even when the SIMD data is double precision.
 */
#    define gmx_simd4_invsqrt_singleaccuracy_r       gmx_simd4_invsqrt_f

/*! \} */

#endif /* GMX_DOUBLE */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_SIMD_MATH_H_ */
