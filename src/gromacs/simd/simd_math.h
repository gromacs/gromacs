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
#ifndef GMX_SIMD_SIMD_MATH_H
#define GMX_SIMD_SIMD_MATH_H

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

#include <cmath>

#include "gromacs/math/utilities.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

#if GMX_SIMD

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name Implementation accuracy settings
 *  \{
 */

/*! \} */

#if GMX_SIMD_HAVE_FLOAT


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
 * You should normally call the real-precision routine \ref simdSum4.
 *
 * \param a term 1 (multiple values)
 * \param b term 2 (multiple values)
 * \param c term 3 (multiple values)
 * \param d term 4 (multiple values)
 * \return sum of terms 1-4 (multiple values)
 */
static inline SimdFloat gmx_simdcall
simdSum4F(SimdFloat a, SimdFloat b,
          SimdFloat c, SimdFloat d)
{
    return simdAddF(simdAddF(a, b), simdAddF(c, d));
}

/*! \brief Return -a if b is negative, SIMD float.
 *
 * You should normally call the real-precision routine \ref simdXorSign.
 *
 * \param a Values to set sign for
 * \param b Values used to set sign
 * \return if b is negative, the sign of a will be changed.
 *
 * This is equivalent to doing an xor operation on a with the sign bit of b,
 * with the exception that negative zero is not considered to be negative
 * on architectures where \ref GMX_SIMD_HAVE_LOGICAL is not set.
 */
static inline SimdFloat gmx_simdcall
simdXorSignF(SimdFloat a, SimdFloat b)
{
#if GMX_SIMD_HAVE_LOGICAL
    return simdXorF(a, simdAndF(simdSet1F(GMX_FLOAT_NEGZERO), b));
#else
    return simdBlendF(a, simdNegF(a), simdCmpLtF(b, simdSetZeroF()));
#endif
}

#ifndef simdRsqrtIterF
/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD float.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the inverse square root.
 *
 *  \param lu Approximation of 1/sqrt(x), typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/sqrt(x).
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static inline SimdFloat gmx_simdcall
simdRsqrtIterF(SimdFloat lu, SimdFloat x)
{
#    if GMX_SIMD_HAVE_FMA
    return simdFmaddF(simdFnmaddF(x, simdMulF(lu, lu), simdSet1F(1.0f)), simdMulF(lu, simdSet1F(0.5f)), lu);
#    else
    return simdMulF(simdSet1F(0.5f), simdMulF(simdSubF(simdSet1F(3.0f), simdMulF(simdMulF(lu, lu), x)), lu));
#    endif
}
#endif

/*! \brief Calculate 1/sqrt(x) for SIMD float.
 *
 * You should normally call the real-precision routine \ref simdInvsqrt.
 *
 *  \param x Argument that must be >0. This routine does not check arguments.
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline SimdFloat gmx_simdcall
simdInvsqrtF(SimdFloat x)
{
    SimdFloat lu = simdRsqrtF(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterF(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterF(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterF(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for two SIMD floats.
 *
 * You should normally call the real-precision routine \ref simdInvsqrtPair.
 *
 * \param x0  First set of arguments, x0 must be positive - no argument checking.
 * \param x1  Second set of arguments, x1 must be positive - no argument checking.
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 *  In particular for double precision we can sometimes calculate square root
 *  pairs slightly faster by using single precision until the very last step.
 */
static inline void gmx_simdcall
simdInvsqrtPairF(SimdFloat x0,    SimdFloat x1,
                 SimdFloat *out0, SimdFloat *out1)
{
    *out0 = simdInvsqrtF(x0);
    *out1 = simdInvsqrtF(x1);
}

#ifndef simdRcpIterF
/*! \brief Perform one Newton-Raphson iteration to improve 1/x for SIMD float.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the reciprocal.
 *
 *  \param lu Approximation of 1/x, typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/x.
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static inline SimdFloat gmx_simdcall
simdRcpIterF(SimdFloat lu, SimdFloat x)
{
    return simdMulF(lu, simdFnmaddF(lu, x, simdSet1F(2.0f)));
}
#endif

/*! \brief Calculate 1/x for SIMD float.
 *
 * You should normally call the real-precision routine \ref simdInv.
 *
 *  \param x Argument that must be nonzero. This routine does not check arguments.
 *  \return 1/x. Result is undefined if your argument was invalid.
 */
static inline SimdFloat gmx_simdcall
simdInvF(SimdFloat x)
{
    SimdFloat lu = simdRcpF(x);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterF(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterF(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterF(lu, x);
#endif
    return lu;
}


/*! \brief Calculate 1/sqrt(x) for masked entries of SIMD float.
 *
 *  This routine only evaluates 1/sqrt(x) for elements for which mask is true.
 *  Illegal values in the masked-out elements will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be >0 for masked-in entries
 *  \param m Mask
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 */
static inline SimdFloat
simdInvsqrtMaskF(SimdFloat x, SimdFBool m)
{
    SimdFloat lu = simdRsqrtMaskF(x, m);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterF(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterF(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterF(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/x for SIMD float, masked version.
 *
 * You should normally call the real-precision routine \ref gmx::simdInv.
 *
 *  \param x Argument that must be nonzero for non-masked entries.
 *  \param m Mask
 *  \return 1/x for elements where m is true, or 0.0 for masked-out entries.
 */
static inline SimdFloat gmx_simdcall
simdInvMaskF(SimdFloat x, SimdFBool m)
{
    SimdFloat lu = simdRcpMaskF(x, m);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterF(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterF(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterF(lu, x);
#endif
    return lu;
}

/*! \brief Calculate sqrt(x) correctly for SIMD floats, including argument 0.0.
 *
 * You should normally call the real-precision routine \ref simdSqrt.
 *
 *  \param x Argument that must be >=0.
 *  \return sqrt(x). If x=0, the result will correctly be set to 0.
 *          The result is undefined if the input value is negative.
 */
static inline SimdFloat gmx_simdcall
simdSqrtF(SimdFloat x)
{
    SimdFloat  res;

    res  = simdInvsqrtMaskF(x, simdCmpNzF(x));
    return simdMulF(res, x);
}

/*! \brief SIMD float log(x). This is the natural logarithm.
 *
 * You should normally call the real-precision routine \ref simdLog.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 */
#ifndef simdLogF
static inline SimdFloat gmx_simdcall
simdLogF(SimdFloat x)
{
    const SimdFloat  half       = simdSet1F(0.5f);
    const SimdFloat  one        = simdSet1F(1.0f);
    const SimdFloat  sqrt2      = simdSet1F(sqrt(2.0f));
    const SimdFloat  corr       = simdSet1F(0.693147180559945286226764f);
    const SimdFloat  CL9        = simdSet1F(0.2371599674224853515625f);
    const SimdFloat  CL7        = simdSet1F(0.285279005765914916992188f);
    const SimdFloat  CL5        = simdSet1F(0.400005519390106201171875f);
    const SimdFloat  CL3        = simdSet1F(0.666666567325592041015625f);
    const SimdFloat  CL1        = simdSet1F(2.0f);
    SimdFloat        fexp, x2, p;
    SimdFBool        mask;

    fexp  = simdGetExponentF(x);
    x     = simdGetMantissaF(x);

    mask  = simdCmpLtF(sqrt2, x);
    /* Adjust to non-IEEE format for x>sqrt(2): exponent += 1, mantissa *= 0.5 */
    fexp  = simdAddF(fexp, simdMaskF(one, mask));
    x     = simdMulF(x, simdBlendF(one, half, mask));

    x     = simdMulF( simdSubF(x, one), simdInvF( simdAddF(x, one) ) );
    x2    = simdMulF(x, x);

    p     = simdFmaddF(CL9, x2, CL7);
    p     = simdFmaddF(p, x2, CL5);
    p     = simdFmaddF(p, x2, CL3);
    p     = simdFmaddF(p, x2, CL1);
    p     = simdFmaddF(p, x, simdMulF(corr, fexp));

    return p;
}
#endif

#ifndef simdExp2F
/*! \brief SIMD float 2^x.
 *
 * You should normally call the real-precision routine \ref simdExp2.
 *
 * \param x Argument.
 * \result 2^x. Undefined if input argument caused overflow.
 */
static inline SimdFloat gmx_simdcall
simdExp2F(SimdFloat x)
{
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const SimdFloat  arglimit = simdSet1F(126.0f);
    const SimdFloat  CC6      = simdSet1F(0.0001534581200287996416911311);
    const SimdFloat  CC5      = simdSet1F(0.001339993121934088894618990);
    const SimdFloat  CC4      = simdSet1F(0.009618488957115180159497841);
    const SimdFloat  CC3      = simdSet1F(0.05550328776964726865751735);
    const SimdFloat  CC2      = simdSet1F(0.2402264689063408646490722);
    const SimdFloat  CC1      = simdSet1F(0.6931472057372680777553816);
    const SimdFloat  one      = simdSet1F(1.0f);

    SimdFloat        fexppart;
    SimdFloat        intpart;
    SimdFloat        p;
    SimdFBool        valuemask;

    fexppart  = simdSetExponentF(x);
    intpart   = simdRoundF(x);
    valuemask = simdCmpLeF(simdAbsF(x), arglimit);
    fexppart  = simdMaskF(fexppart, valuemask);
    x         = simdSubF(x, intpart);

    p         = simdFmaddF(CC6, x, CC5);
    p         = simdFmaddF(p, x, CC4);
    p         = simdFmaddF(p, x, CC3);
    p         = simdFmaddF(p, x, CC2);
    p         = simdFmaddF(p, x, CC1);
    p         = simdFmaddF(p, x, one);
    x         = simdMulF(p, fexppart);
    return x;
}
#endif

#ifndef simdExpF
/*! \brief SIMD float exp(x).
 *
 * You should normally call the real-precision routine \ref simdExp.
 *
 * In addition to scaling the argument for 2^x this routine correctly does
 * extended precision arithmetics to improve accuracy.
 *
 * \param x Argument.
 * \result exp(x). Undefined if input argument caused overflow,
 * which can happen if abs(x) \> 7e13.
 */
static inline SimdFloat gmx_simdcall
simdExpF(SimdFloat x)
{
    const SimdFloat  argscale     = simdSet1F(1.44269504088896341f);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const SimdFloat  arglimit     = simdSet1F(126.0f);
    const SimdFloat  invargscale0 = simdSet1F(-0.693145751953125f);
    const SimdFloat  invargscale1 = simdSet1F(-1.428606765330187045e-06f);
    const SimdFloat  CC4          = simdSet1F(0.00136324646882712841033936f);
    const SimdFloat  CC3          = simdSet1F(0.00836596917361021041870117f);
    const SimdFloat  CC2          = simdSet1F(0.0416710823774337768554688f);
    const SimdFloat  CC1          = simdSet1F(0.166665524244308471679688f);
    const SimdFloat  CC0          = simdSet1F(0.499999850988388061523438f);
    const SimdFloat  one          = simdSet1F(1.0f);
    SimdFloat        fexppart;
    SimdFloat        intpart;
    SimdFloat        y, p;
    SimdFBool        valuemask;

    y         = simdMulF(x, argscale);
    fexppart  = simdSetExponentF(y);  /* rounds to nearest int internally */
    intpart   = simdRoundF(y);        /* use same rounding algorithm here */
    valuemask = simdCmpLeF(simdAbsF(y), arglimit);
    fexppart  = simdMaskF(fexppart, valuemask);

    /* Extended precision arithmetics */
    x         = simdFmaddF(invargscale0, intpart, x);
    x         = simdFmaddF(invargscale1, intpart, x);

    p         = simdFmaddF(CC4, x, CC3);
    p         = simdFmaddF(p, x, CC2);
    p         = simdFmaddF(p, x, CC1);
    p         = simdFmaddF(p, x, CC0);
    p         = simdFmaddF(simdMulF(x, x), p, x);
    p         = simdAddF(p, one);
    x         = simdMulF(p, fexppart);
    return x;
}
#endif

/*! \brief SIMD float erf(x).
 *
 * You should normally call the real-precision routine \ref simdErf.
 *
 * \param x The value to calculate erf(x) for.
 * \result erf(x)
 *
 * This routine achieves very close to full precision, but we do not care about
 * the last bit or the subnormal result range.
 */
static inline SimdFloat gmx_simdcall
simdErfF(SimdFloat x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const SimdFloat  CA6      = simdSet1F(7.853861353153693e-5f);
    const SimdFloat  CA5      = simdSet1F(-8.010193625184903e-4f);
    const SimdFloat  CA4      = simdSet1F(5.188327685732524e-3f);
    const SimdFloat  CA3      = simdSet1F(-2.685381193529856e-2f);
    const SimdFloat  CA2      = simdSet1F(1.128358514861418e-1f);
    const SimdFloat  CA1      = simdSet1F(-3.761262582423300e-1f);
    const SimdFloat  CA0      = simdSet1F(1.128379165726710f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const SimdFloat  CB9      = simdSet1F(-0.0018629930017603923f);
    const SimdFloat  CB8      = simdSet1F(0.003909821287598495f);
    const SimdFloat  CB7      = simdSet1F(-0.0052094582210355615f);
    const SimdFloat  CB6      = simdSet1F(0.005685614362160572f);
    const SimdFloat  CB5      = simdSet1F(-0.0025367682853477272f);
    const SimdFloat  CB4      = simdSet1F(-0.010199799682318782f);
    const SimdFloat  CB3      = simdSet1F(0.04369575504816542f);
    const SimdFloat  CB2      = simdSet1F(-0.11884063474674492f);
    const SimdFloat  CB1      = simdSet1F(0.2732120154030589f);
    const SimdFloat  CB0      = simdSet1F(0.42758357702025784f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const SimdFloat  CC10     = simdSet1F(-0.0445555913112064f);
    const SimdFloat  CC9      = simdSet1F(0.21376355144663348f);
    const SimdFloat  CC8      = simdSet1F(-0.3473187200259257f);
    const SimdFloat  CC7      = simdSet1F(0.016690861551248114f);
    const SimdFloat  CC6      = simdSet1F(0.7560973182491192f);
    const SimdFloat  CC5      = simdSet1F(-1.2137903600145787f);
    const SimdFloat  CC4      = simdSet1F(0.8411872321232948f);
    const SimdFloat  CC3      = simdSet1F(-0.08670413896296343f);
    const SimdFloat  CC2      = simdSet1F(-0.27124782687240334f);
    const SimdFloat  CC1      = simdSet1F(-0.0007502488047806069f);
    const SimdFloat  CC0      = simdSet1F(0.5642114853803148f);
    const SimdFloat  one      = simdSet1F(1.0f);
    const SimdFloat  two      = simdSet1F(2.0f);

    SimdFloat        x2, x4, y;
    SimdFloat        t, t2, w, w2;
    SimdFloat        pA0, pA1, pB0, pB1, pC0, pC1;
    SimdFloat        expmx2;
    SimdFloat        res_erf, res_erfc, res;
    SimdFBool        mask, msk_erf;

    /* Calculate erf() */
    x2   = simdMulF(x, x);
    x4   = simdMulF(x2, x2);

    pA0  = simdFmaddF(CA6, x4, CA4);
    pA1  = simdFmaddF(CA5, x4, CA3);
    pA0  = simdFmaddF(pA0, x4, CA2);
    pA1  = simdFmaddF(pA1, x4, CA1);
    pA0  = simdMulF(pA0, x4);
    pA0  = simdFmaddF(pA1, x2, pA0);
    /* Constant term must come last for precision reasons */
    pA0  = simdAddF(pA0, CA0);

    res_erf = simdMulF(x, pA0);

    /* Calculate erfc */
    y       = simdAbsF(x);
    msk_erf = simdCmpLeF(simdSet1F(0.75f), y);
    t       = simdInvMaskF(y, msk_erf);
    w       = simdSubF(t, one);
    t2      = simdMulF(t, t);
    w2      = simdMulF(w, w);

    /* No need for a floating-point sieve here (as in erfc), since erf()
     * will never return values that are extremely small for large args.
     */
    expmx2  = simdExpF( simdNegF( simdMulF(y, y)));

    pB1  = simdFmaddF(CB9, w2, CB7);
    pB0  = simdFmaddF(CB8, w2, CB6);
    pB1  = simdFmaddF(pB1, w2, CB5);
    pB0  = simdFmaddF(pB0, w2, CB4);
    pB1  = simdFmaddF(pB1, w2, CB3);
    pB0  = simdFmaddF(pB0, w2, CB2);
    pB1  = simdFmaddF(pB1, w2, CB1);
    pB0  = simdFmaddF(pB0, w2, CB0);
    pB0  = simdFmaddF(pB1, w, pB0);

    pC0  = simdFmaddF(CC10, t2, CC8);
    pC1  = simdFmaddF(CC9, t2, CC7);
    pC0  = simdFmaddF(pC0, t2, CC6);
    pC1  = simdFmaddF(pC1, t2, CC5);
    pC0  = simdFmaddF(pC0, t2, CC4);
    pC1  = simdFmaddF(pC1, t2, CC3);
    pC0  = simdFmaddF(pC0, t2, CC2);
    pC1  = simdFmaddF(pC1, t2, CC1);

    pC0  = simdFmaddF(pC0, t2, CC0);
    pC0  = simdFmaddF(pC1, t, pC0);
    pC0  = simdMulF(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = simdCmpLtF(two, y);
    res_erfc = simdBlendF(pB0, pC0, mask);
    res_erfc = simdMulF(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = simdCmpLtF(x, simdSetZeroF());
    res_erfc = simdBlendF(res_erfc, simdSubF(two, res_erfc), mask);

    /* Select erf() or erfc() */
    res  = simdBlendF(res_erf, simdSubF(one, res_erfc), msk_erf);

    return res;
}

/*! \brief SIMD float erfc(x).
 *
 * You should normally call the real-precision routine \ref simdErfc.
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
static inline SimdFloat gmx_simdcall
simdErfcF(SimdFloat x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const SimdFloat  CA6      = simdSet1F(7.853861353153693e-5f);
    const SimdFloat  CA5      = simdSet1F(-8.010193625184903e-4f);
    const SimdFloat  CA4      = simdSet1F(5.188327685732524e-3f);
    const SimdFloat  CA3      = simdSet1F(-2.685381193529856e-2f);
    const SimdFloat  CA2      = simdSet1F(1.128358514861418e-1f);
    const SimdFloat  CA1      = simdSet1F(-3.761262582423300e-1f);
    const SimdFloat  CA0      = simdSet1F(1.128379165726710f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const SimdFloat  CB9      = simdSet1F(-0.0018629930017603923f);
    const SimdFloat  CB8      = simdSet1F(0.003909821287598495f);
    const SimdFloat  CB7      = simdSet1F(-0.0052094582210355615f);
    const SimdFloat  CB6      = simdSet1F(0.005685614362160572f);
    const SimdFloat  CB5      = simdSet1F(-0.0025367682853477272f);
    const SimdFloat  CB4      = simdSet1F(-0.010199799682318782f);
    const SimdFloat  CB3      = simdSet1F(0.04369575504816542f);
    const SimdFloat  CB2      = simdSet1F(-0.11884063474674492f);
    const SimdFloat  CB1      = simdSet1F(0.2732120154030589f);
    const SimdFloat  CB0      = simdSet1F(0.42758357702025784f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const SimdFloat  CC10     = simdSet1F(-0.0445555913112064f);
    const SimdFloat  CC9      = simdSet1F(0.21376355144663348f);
    const SimdFloat  CC8      = simdSet1F(-0.3473187200259257f);
    const SimdFloat  CC7      = simdSet1F(0.016690861551248114f);
    const SimdFloat  CC6      = simdSet1F(0.7560973182491192f);
    const SimdFloat  CC5      = simdSet1F(-1.2137903600145787f);
    const SimdFloat  CC4      = simdSet1F(0.8411872321232948f);
    const SimdFloat  CC3      = simdSet1F(-0.08670413896296343f);
    const SimdFloat  CC2      = simdSet1F(-0.27124782687240334f);
    const SimdFloat  CC1      = simdSet1F(-0.0007502488047806069f);
    const SimdFloat  CC0      = simdSet1F(0.5642114853803148f);
    /* Coefficients for expansion of exp(x) in [0,0.1] */
    /* CD0 and CD1 are both 1.0, so no need to declare them separately */
    const SimdFloat  CD2      = simdSet1F(0.5000066608081202f);
    const SimdFloat  CD3      = simdSet1F(0.1664795422874624f);
    const SimdFloat  CD4      = simdSet1F(0.04379839977652482f);
    const SimdFloat  one      = simdSet1F(1.0f);
    const SimdFloat  two      = simdSet1F(2.0f);

    /* We need to use a small trick here, since we cannot assume all SIMD
     * architectures support integers, and the flag we want (0xfffff000) would
     * evaluate to NaN (i.e., it cannot be expressed as a floating-point num).
     * Instead, we represent the flags 0xf0f0f000 and 0x0f0f0000 as valid
     * fp numbers, and perform a logical or. Since the expression is constant,
     * we can at least hope it is evaluated at compile-time.
     */
#if GMX_SIMD_HAVE_LOGICAL
    const SimdFloat         sieve    = simdOrF(simdSet1F(-5.965323564e+29f), simdSet1F(7.05044434e-30f));
#else
    const int               isieve   = 0xFFFFF000;
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  mem[GMX_SIMD_FLOAT_WIDTH];

    union {
        float f; int i;
    } conv;
    int                     i;
#endif

    SimdFloat        x2, x4, y;
    SimdFloat        q, z, t, t2, w, w2;
    SimdFloat        pA0, pA1, pB0, pB1, pC0, pC1;
    SimdFloat        expmx2, corr;
    SimdFloat        res_erf, res_erfc, res;
    SimdFBool        mask, msk_erf;

    /* Calculate erf() */
    x2     = simdMulF(x, x);
    x4     = simdMulF(x2, x2);

    pA0  = simdFmaddF(CA6, x4, CA4);
    pA1  = simdFmaddF(CA5, x4, CA3);
    pA0  = simdFmaddF(pA0, x4, CA2);
    pA1  = simdFmaddF(pA1, x4, CA1);
    pA1  = simdMulF(pA1, x2);
    pA0  = simdFmaddF(pA0, x4, pA1);
    /* Constant term must come last for precision reasons */
    pA0  = simdAddF(pA0, CA0);

    res_erf = simdMulF(x, pA0);

    /* Calculate erfc */
    y       = simdAbsF(x);
    msk_erf = simdCmpLeF(simdSet1F(0.75f), y);
    t       = simdInvMaskF(y, msk_erf);
    w       = simdSubF(t, one);
    t2      = simdMulF(t, t);
    w2      = simdMulF(w, w);
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
#if GMX_SIMD_HAVE_LOGICAL
    z       = simdAndF(y, sieve);
#else
    simdStoreF(mem, y);
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f  = mem[i];
        conv.i  = conv.i & isieve;
        mem[i]  = conv.f;
    }
    z = simdLoadF(mem);
#endif
    q       = simdMulF( simdSubF(z, y), simdAddF(z, y) );
    corr    = simdFmaddF(CD4, q, CD3);
    corr    = simdFmaddF(corr, q, CD2);
    corr    = simdFmaddF(corr, q, one);
    corr    = simdFmaddF(corr, q, one);

    expmx2  = simdExpF( simdNegF( simdMulF(z, z) ) );
    expmx2  = simdMulF(expmx2, corr);

    pB1  = simdFmaddF(CB9, w2, CB7);
    pB0  = simdFmaddF(CB8, w2, CB6);
    pB1  = simdFmaddF(pB1, w2, CB5);
    pB0  = simdFmaddF(pB0, w2, CB4);
    pB1  = simdFmaddF(pB1, w2, CB3);
    pB0  = simdFmaddF(pB0, w2, CB2);
    pB1  = simdFmaddF(pB1, w2, CB1);
    pB0  = simdFmaddF(pB0, w2, CB0);
    pB0  = simdFmaddF(pB1, w, pB0);

    pC0  = simdFmaddF(CC10, t2, CC8);
    pC1  = simdFmaddF(CC9, t2, CC7);
    pC0  = simdFmaddF(pC0, t2, CC6);
    pC1  = simdFmaddF(pC1, t2, CC5);
    pC0  = simdFmaddF(pC0, t2, CC4);
    pC1  = simdFmaddF(pC1, t2, CC3);
    pC0  = simdFmaddF(pC0, t2, CC2);
    pC1  = simdFmaddF(pC1, t2, CC1);

    pC0  = simdFmaddF(pC0, t2, CC0);
    pC0  = simdFmaddF(pC1, t, pC0);
    pC0  = simdMulF(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = simdCmpLtF(two, y);
    res_erfc = simdBlendF(pB0, pC0, mask);
    res_erfc = simdMulF(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = simdCmpLtF(x, simdSetZeroF());
    res_erfc = simdBlendF(res_erfc, simdSubF(two, res_erfc), mask);

    /* Select erf() or erfc() */
    res  = simdBlendF(simdSubF(one, res_erf), res_erfc, msk_erf);

    return res;
}

/*! \brief SIMD float sin \& cos.
 *
 * You should normally call the real-precision routine \ref simdSinCos.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 * This version achieves close to machine precision, but for very large
 * magnitudes of the argument we inherently begin to lose accuracy due to the
 * argument reduction, despite using extended precision arithmetics internally.
 */
static inline void gmx_simdcall
simdSinCosF(SimdFloat x, SimdFloat *sinval, SimdFloat *cosval)
{
    /* Constants to subtract Pi/4*x from y while minimizing precision loss */
    const SimdFloat  argred0         = simdSet1F(-1.5703125);
    const SimdFloat  argred1         = simdSet1F(-4.83751296997070312500e-04f);
    const SimdFloat  argred2         = simdSet1F(-7.54953362047672271729e-08f);
    const SimdFloat  argred3         = simdSet1F(-2.56334406825708960298e-12f);
    const SimdFloat  two_over_pi     = simdSet1F(2.0f/M_PI);
    const SimdFloat  const_sin2      = simdSet1F(-1.9515295891e-4f);
    const SimdFloat  const_sin1      = simdSet1F( 8.3321608736e-3f);
    const SimdFloat  const_sin0      = simdSet1F(-1.6666654611e-1f);
    const SimdFloat  const_cos2      = simdSet1F( 2.443315711809948e-5f);
    const SimdFloat  const_cos1      = simdSet1F(-1.388731625493765e-3f);
    const SimdFloat  const_cos0      = simdSet1F( 4.166664568298827e-2f);
    const SimdFloat  half            = simdSet1F(0.5f);
    const SimdFloat  one             = simdSet1F(1.0f);
    SimdFloat        ssign, csign;
    SimdFloat        x2, y, z, psin, pcos, sss, ccc;
    SimdFBool        mask;

#if GMX_SIMD_HAVE_FINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    const SimdFInt32 ione            = simdSet1FI(1);
    const SimdFInt32 itwo            = simdSet1FI(2);
    SimdFInt32       iy;

    z       = simdMulF(x, two_over_pi);
    iy      = simdCvtF2I(z);
    y       = simdRoundF(z);

    mask    = simdCvtFIB2FB(simdCmpEqFI(simdAndFI(iy, ione), simdSetZeroFI()));
    ssign   = simdMaskF(simdSet1F(GMX_FLOAT_NEGZERO), simdCvtFIB2FB(simdCmpEqFI(simdAndFI(iy, itwo), itwo)));
    csign   = simdMaskF(simdSet1F(GMX_FLOAT_NEGZERO), simdCvtFIB2FB(simdCmpEqFI(simdAndFI(simdAddFI(iy, ione), itwo), itwo)));
#else
    const SimdFloat  quarter         = simdSet1F(0.25f);
    const SimdFloat  minusquarter    = simdSet1F(-0.25f);
    SimdFloat        q;
    SimdFBool        m1, m2, m3;

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
    x       = simdAbsF(x);
    /* It is critical that half-way cases are rounded down */
    z       = simdFmaddF(x, two_over_pi, half);
    y       = simdTruncF(z);
    q       = simdMulF(z, quarter);
    q       = simdSubF(q, simdTruncF(q));
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q       = simdSubF(q, half);
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#    if GMX_SIMD_HAVE_LOGICAL
    ssign   = simdAndF(ssign, simdSet1F(GMX_FLOAT_NEGZERO));
    csign   = simdAndNotF(q, simdSet1F(GMX_FLOAT_NEGZERO));
    ssign   = simdXorF(ssign, csign);
#    else
    csign   = simdXorSignF(simdSet1F(-1.0f), q);

    ssign   = simdXorSignF(ssign, csign);    /* swap ssign if csign was set. */
#    endif
    /* Check if y had value 1 or 3 (remember we subtracted 0.5 from q) */
    m1      = simdCmpLtF(q, minusquarter);
    m2      = simdCmpLeF(simdSetZeroF(), q);
    m3      = simdCmpLtF(q, quarter);
    m2      = simdAndFB(m2, m3);
    mask    = simdOrFB(m1, m2);
    /* where mask is FALSE, set sign. */
    csign   = simdXorSignF(csign, simdBlendF(simdSet1F(-1.0f), one, mask));
#endif
    x       = simdFmaddF(y, argred0, x);
    x       = simdFmaddF(y, argred1, x);
    x       = simdFmaddF(y, argred2, x);
    x       = simdFmaddF(y, argred3, x);
    x2      = simdMulF(x, x);

    psin    = simdFmaddF(const_sin2, x2, const_sin1);
    psin    = simdFmaddF(psin, x2, const_sin0);
    psin    = simdFmaddF(psin, simdMulF(x, x2), x);
    pcos    = simdFmaddF(const_cos2, x2, const_cos1);
    pcos    = simdFmaddF(pcos, x2, const_cos0);
    pcos    = simdFmsubF(pcos, x2, half);
    pcos    = simdFmaddF(pcos, x2, one);

    sss     = simdBlendF(pcos, psin, mask);
    ccc     = simdBlendF(psin, pcos, mask);
    /* See comment for GMX_SIMD_HAVE_LOGICAL section above. */
#if GMX_SIMD_HAVE_LOGICAL
    *sinval = simdXorF(sss, ssign);
    *cosval = simdXorF(ccc, csign);
#else
    *sinval = simdXorSignF(sss, ssign);
    *cosval = simdXorSignF(ccc, csign);
#endif
}

/*! \brief SIMD float sin(x).
 *
 * You should normally call the real-precision routine \ref simdSin.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref simdSinCos and waste a factor 2 in performance.
 */
static inline SimdFloat gmx_simdcall
simdSinF(SimdFloat x)
{
    SimdFloat s, c;
    simdSinCosF(x, &s, &c);
    return s;
}

/*! \brief SIMD float cos(x).
 *
 * You should normally call the real-precision routine \ref simdCos.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref simdSinCos and waste a factor 2 in performance.
 */
static inline SimdFloat gmx_simdcall
simdCosF(SimdFloat x)
{
    SimdFloat s, c;
    simdSinCosF(x, &s, &c);
    return c;
}

/*! \brief SIMD float tan(x).
 *
 * You should normally call the real-precision routine \ref simdTan.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 */
static inline SimdFloat gmx_simdcall
simdTanF(SimdFloat x)
{
    const SimdFloat  argred0         = simdSet1F(-1.5703125);
    const SimdFloat  argred1         = simdSet1F(-4.83751296997070312500e-04f);
    const SimdFloat  argred2         = simdSet1F(-7.54953362047672271729e-08f);
    const SimdFloat  argred3         = simdSet1F(-2.56334406825708960298e-12f);
    const SimdFloat  two_over_pi     = simdSet1F(2.0f/M_PI);
    const SimdFloat  CT6             = simdSet1F(0.009498288995810566122993911);
    const SimdFloat  CT5             = simdSet1F(0.002895755790837379295226923);
    const SimdFloat  CT4             = simdSet1F(0.02460087336161924491836265);
    const SimdFloat  CT3             = simdSet1F(0.05334912882656359828045988);
    const SimdFloat  CT2             = simdSet1F(0.1333989091464957704418495);
    const SimdFloat  CT1             = simdSet1F(0.3333307599244198227797507);

    SimdFloat        x2, p, y, z;
    SimdFBool        mask;

#if GMX_SIMD_HAVE_FINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    SimdFInt32  iy;
    SimdFInt32  ione = simdSet1FI(1);

    z       = simdMulF(x, two_over_pi);
    iy      = simdCvtF2I(z);
    y       = simdRoundF(z);
    mask    = simdCvtFIB2FB(simdCmpEqFI(simdAndFI(iy, ione), ione));

    x       = simdFmaddF(y, argred0, x);
    x       = simdFmaddF(y, argred1, x);
    x       = simdFmaddF(y, argred2, x);
    x       = simdFmaddF(y, argred3, x);
    x       = simdXorF(simdMaskF(simdSet1F(GMX_FLOAT_NEGZERO), mask), x);
#else
    const SimdFloat  quarter         = simdSet1F(0.25f);
    const SimdFloat  half            = simdSet1F(0.5f);
    const SimdFloat  threequarter    = simdSet1F(0.75f);
    SimdFloat        w, q;
    SimdFBool        m1, m2, m3;

    w       = simdAbsF(x);
    z       = simdFmaddF(w, two_over_pi, half);
    y       = simdTruncF(z);
    q       = simdMulF(z, quarter);
    q       = simdSubF(q, simdTruncF(q));
    m1      = simdCmpLeF(quarter, q);
    m2      = simdCmpLtF(q, half);
    m3      = simdCmpLeF(threequarter, q);
    m1      = simdAndFB(m1, m2);
    mask    = simdOrFB(m1, m3);
    w       = simdFmaddF(y, argred0, w);
    w       = simdFmaddF(y, argred1, w);
    w       = simdFmaddF(y, argred2, w);
    w       = simdFmaddF(y, argred3, w);

    w       = simdBlendF(w, simdNegF(w), mask);
    x       = simdXorSignF(w, x);
#endif
    x2      = simdMulF(x, x);
    p       = simdFmaddF(CT6, x2, CT5);
    p       = simdFmaddF(p, x2, CT4);
    p       = simdFmaddF(p, x2, CT3);
    p       = simdFmaddF(p, x2, CT2);
    p       = simdFmaddF(p, x2, CT1);
    p       = simdFmaddF(x2, simdMulF(p, x), x);

    p       = simdBlendF( p, simdInvMaskF(p, mask), mask);
    return p;
}

/*! \brief SIMD float asin(x).
 *
 * You should normally call the real-precision routine \ref simdAsin.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 */
static inline SimdFloat gmx_simdcall
simdAsinF(SimdFloat x)
{
    const SimdFloat limitlow   = simdSet1F(1e-4f);
    const SimdFloat half       = simdSet1F(0.5f);
    const SimdFloat one        = simdSet1F(1.0f);
    const SimdFloat halfpi     = simdSet1F((float)M_PI/2.0f);
    const SimdFloat CC5        = simdSet1F(4.2163199048E-2f);
    const SimdFloat CC4        = simdSet1F(2.4181311049E-2f);
    const SimdFloat CC3        = simdSet1F(4.5470025998E-2f);
    const SimdFloat CC2        = simdSet1F(7.4953002686E-2f);
    const SimdFloat CC1        = simdSet1F(1.6666752422E-1f);
    SimdFloat       xabs;
    SimdFloat       z, z1, z2, q, q1, q2;
    SimdFloat       pA, pB;
    SimdFBool       mask, mask2;

    xabs  = simdAbsF(x);
    mask  = simdCmpLtF(half, xabs);
    z1    = simdMulF(half, simdSubF(one, xabs));
    mask2 = simdCmpLtF(xabs, one);
    q1    = simdMulF(z1, simdInvsqrtMaskF(z1, mask2));
    q2    = xabs;
    z2    = simdMulF(q2, q2);
    z     = simdBlendF(z2, z1, mask);
    q     = simdBlendF(q2, q1, mask);

    z2    = simdMulF(z, z);
    pA    = simdFmaddF(CC5, z2, CC3);
    pB    = simdFmaddF(CC4, z2, CC2);
    pA    = simdFmaddF(pA, z2, CC1);
    pA    = simdMulF(pA, z);
    z     = simdFmaddF(pB, z2, pA);
    z     = simdFmaddF(z, q, q);
    q2    = simdSubF(halfpi, z);
    q2    = simdSubF(q2, z);
    z     = simdBlendF(z, q2, mask);

    mask  = simdCmpLtF(limitlow, xabs);
    z     = simdBlendF( xabs, z, mask );
    z     = simdXorSignF(z, x);

    return z;
}

/*! \brief SIMD float acos(x).
 *
 * You should normally call the real-precision routine \ref simdAcos.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 */
static inline SimdFloat gmx_simdcall
simdAcosF(SimdFloat x)
{
    const SimdFloat one       = simdSet1F(1.0f);
    const SimdFloat half      = simdSet1F(0.5f);
    const SimdFloat pi        = simdSet1F((float)M_PI);
    const SimdFloat halfpi    = simdSet1F((float)M_PI/2.0f);
    SimdFloat       xabs;
    SimdFloat       z, z1, z2, z3;
    SimdFBool       mask1, mask2, mask3;

    xabs  = simdAbsF(x);
    mask1 = simdCmpLtF(half, xabs);
    mask2 = simdCmpLtF(simdSetZeroF(), x);

    z     = simdMulF(half, simdSubF(one, xabs));
    mask3 = simdCmpLtF(xabs, one);
    z     = simdMulF(z, simdInvsqrtMaskF(z, mask3));
    z     = simdBlendF(x, z, mask1);
    z     = simdAsinF(z);

    z2    = simdAddF(z, z);
    z1    = simdSubF(pi, z2);
    z3    = simdSubF(halfpi, z);
    z     = simdBlendF(z1, z2, mask2);
    z     = simdBlendF(z3, z, mask1);

    return z;
}

/*! \brief SIMD float asin(x).
 *
 * You should normally call the real-precision routine \ref simdAtan.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x), same argument/value range as standard math library.
 */
static inline SimdFloat gmx_simdcall
simdAtanF(SimdFloat x)
{
    const SimdFloat halfpi    = simdSet1F(M_PI/2);
    const SimdFloat CA17      = simdSet1F(0.002823638962581753730774f);
    const SimdFloat CA15      = simdSet1F(-0.01595690287649631500244f);
    const SimdFloat CA13      = simdSet1F(0.04250498861074447631836f);
    const SimdFloat CA11      = simdSet1F(-0.07489009201526641845703f);
    const SimdFloat CA9       = simdSet1F(0.1063479334115982055664f);
    const SimdFloat CA7       = simdSet1F(-0.1420273631811141967773f);
    const SimdFloat CA5       = simdSet1F(0.1999269574880599975585f);
    const SimdFloat CA3       = simdSet1F(-0.3333310186862945556640f);
    const SimdFloat one       = simdSet1F(1.0f);
    SimdFloat       x2, x3, x4, pA, pB;
    SimdFBool       mask, mask2;

    mask  = simdCmpLtF(x, simdSetZeroF());
    x     = simdAbsF(x);
    mask2 = simdCmpLtF(one, x);
    x     = simdBlendF(x, simdInvMaskF(x, mask2), mask2);

    x2    = simdMulF(x, x);
    x3    = simdMulF(x2, x);
    x4    = simdMulF(x2, x2);
    pA    = simdFmaddF(CA17, x4, CA13);
    pB    = simdFmaddF(CA15, x4, CA11);
    pA    = simdFmaddF(pA, x4, CA9);
    pB    = simdFmaddF(pB, x4, CA7);
    pA    = simdFmaddF(pA, x4, CA5);
    pB    = simdFmaddF(pB, x4, CA3);
    pA    = simdFmaddF(pA, x2, pB);
    pA    = simdFmaddF(pA, x3, x);

    pA    = simdBlendF(pA, simdSubF(halfpi, pA), mask2);
    pA    = simdBlendF(pA, simdNegF(pA), mask);

    return pA;
}

/*! \brief SIMD float atan2(y,x).
 *
 * You should normally call the real-precision routine \ref simdAtan2.
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
static inline SimdFloat gmx_simdcall
simdAtan2F(SimdFloat y, SimdFloat x)
{
    const SimdFloat pi          = simdSet1F(M_PI);
    const SimdFloat halfpi      = simdSet1F(M_PI/2.0);
    SimdFloat       xinv, p, aoffset;
    SimdFBool       mask_xnz, mask_ynz, mask_xlt0, mask_ylt0;

    mask_xnz  = simdCmpNzF(x);
    mask_ynz  = simdCmpNzF(y);
    mask_xlt0 = simdCmpLtF(x, simdSetZeroF());
    mask_ylt0 = simdCmpLtF(y, simdSetZeroF());

    aoffset   = simdMaskNotF(halfpi, mask_xnz);
    aoffset   = simdMaskF(aoffset, mask_ynz);

    aoffset   = simdBlendF(aoffset, pi, mask_xlt0);
    aoffset   = simdBlendF(aoffset, simdNegF(aoffset), mask_ylt0);

    xinv      = simdInvMaskF(x, mask_xnz);
    p         = simdMulF(y, xinv);
    p         = simdAtanF(p);
    p         = simdAddF(p, aoffset);

    return p;
}

/*! \brief Calculate the force correction due to PME analytically in SIMD float.
 *
 * You should normally call the real-precision routine \ref simdPmeCorrForce.
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
static inline SimdFloat gmx_simdcall
simdPmeCorrForceF(SimdFloat z2)
{
    const SimdFloat  FN6      = simdSet1F(-1.7357322914161492954e-8f);
    const SimdFloat  FN5      = simdSet1F(1.4703624142580877519e-6f);
    const SimdFloat  FN4      = simdSet1F(-0.000053401640219807709149f);
    const SimdFloat  FN3      = simdSet1F(0.0010054721316683106153f);
    const SimdFloat  FN2      = simdSet1F(-0.019278317264888380590f);
    const SimdFloat  FN1      = simdSet1F(0.069670166153766424023f);
    const SimdFloat  FN0      = simdSet1F(-0.75225204789749321333f);

    const SimdFloat  FD4      = simdSet1F(0.0011193462567257629232f);
    const SimdFloat  FD3      = simdSet1F(0.014866955030185295499f);
    const SimdFloat  FD2      = simdSet1F(0.11583842382862377919f);
    const SimdFloat  FD1      = simdSet1F(0.50736591960530292870f);
    const SimdFloat  FD0      = simdSet1F(1.0f);

    SimdFloat        z4;
    SimdFloat        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = simdMulF(z2, z2);

    polyFD0        = simdFmaddF(FD4, z4, FD2);
    polyFD1        = simdFmaddF(FD3, z4, FD1);
    polyFD0        = simdFmaddF(polyFD0, z4, FD0);
    polyFD0        = simdFmaddF(polyFD1, z2, polyFD0);

    polyFD0        = simdInvF(polyFD0);

    polyFN0        = simdFmaddF(FN6, z4, FN4);
    polyFN1        = simdFmaddF(FN5, z4, FN3);
    polyFN0        = simdFmaddF(polyFN0, z4, FN2);
    polyFN1        = simdFmaddF(polyFN1, z4, FN1);
    polyFN0        = simdFmaddF(polyFN0, z4, FN0);
    polyFN0        = simdFmaddF(polyFN1, z2, polyFN0);

    return simdMulF(polyFN0, polyFD0);
}



/*! \brief Calculate the potential correction due to PME analytically in SIMD float.
 *
 * You should normally call the real-precision routine \ref simdPmeCorrPotential.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb potential - see below for details.
 *
 * See \ref simdPmeCorrForceF for details about the approximation.
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
static inline SimdFloat gmx_simdcall
simdPmeCorrPotentialF(SimdFloat z2)
{
    const SimdFloat  VN6      = simdSet1F(1.9296833005951166339e-8f);
    const SimdFloat  VN5      = simdSet1F(-1.4213390571557850962e-6f);
    const SimdFloat  VN4      = simdSet1F(0.000041603292906656984871f);
    const SimdFloat  VN3      = simdSet1F(-0.00013134036773265025626f);
    const SimdFloat  VN2      = simdSet1F(0.038657983986041781264f);
    const SimdFloat  VN1      = simdSet1F(0.11285044772717598220f);
    const SimdFloat  VN0      = simdSet1F(1.1283802385263030286f);

    const SimdFloat  VD3      = simdSet1F(0.0066752224023576045451f);
    const SimdFloat  VD2      = simdSet1F(0.078647795836373922256f);
    const SimdFloat  VD1      = simdSet1F(0.43336185284710920150f);
    const SimdFloat  VD0      = simdSet1F(1.0f);

    SimdFloat        z4;
    SimdFloat        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = simdMulF(z2, z2);

    polyVD1        = simdFmaddF(VD3, z4, VD1);
    polyVD0        = simdFmaddF(VD2, z4, VD0);
    polyVD0        = simdFmaddF(polyVD1, z2, polyVD0);

    polyVD0        = simdInvF(polyVD0);

    polyVN0        = simdFmaddF(VN6, z4, VN4);
    polyVN1        = simdFmaddF(VN5, z4, VN3);
    polyVN0        = simdFmaddF(polyVN0, z4, VN2);
    polyVN1        = simdFmaddF(polyVN1, z4, VN1);
    polyVN0        = simdFmaddF(polyVN0, z4, VN0);
    polyVN0        = simdFmaddF(polyVN1, z2, polyVN0);

    return simdMulF(polyVN0, polyVD0);
}
#endif

/*! \} */

#if GMX_SIMD_HAVE_DOUBLE


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
 * \copydetails simdSum4F
 */
static inline SimdDouble gmx_simdcall
simdSum4D(SimdDouble a, SimdDouble b,
          SimdDouble c, SimdDouble d)
{
    return simdAddD(simdAddD(a, b), simdAddD(c, d));
}

/*! \brief Return -a if b is negative, SIMD double.
 *
 * You should normally call the real-precision routine \ref simdXorSign.
 *
 * \param a Values to set sign for
 * \param b Values used to set sign
 * \return if b is negative, the sign of a will be changed.
 *
 * This is equivalent to doing an xor operation on a with the sign bit of b,
 * with the exception that negative zero is not considered to be negative
 * on architectures where \ref GMX_SIMD_HAVE_LOGICAL is not set.
 */
static inline SimdDouble gmx_simdcall
simdXorSignD(SimdDouble a, SimdDouble b)
{
#if GMX_SIMD_HAVE_LOGICAL
    return simdXorD(a, simdAndD(simdSet1D(GMX_DOUBLE_NEGZERO), b));
#else
    return simdBlendD(a, simdNegD(a), simdCmpLtD(b, simdSetZeroD()));
#endif
}

#ifndef simdRsqrtIterD
/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD double.
 *
 * \copydetails simdRsqrtIterF
 */
static inline SimdDouble gmx_simdcall
simdRsqrtIterD(SimdDouble lu, SimdDouble x)
{
#if GMX_SIMD_HAVE_FMA
    return simdFmaddD(simdFnmaddD(x, simdMulD(lu, lu), simdSet1D(1.0)), simdMulD(lu, simdSet1D(0.5)), lu);
#else
    return simdMulD(simdSet1D(0.5), simdMulD(simdSubD(simdSet1D(3.0), simdMulD(simdMulD(lu, lu), x)), lu));
#endif
}
#endif

/*! \brief Calculate 1/sqrt(x) for SIMD double
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdDouble gmx_simdcall
simdInvsqrtD(SimdDouble x)
{
    SimdDouble lu = simdRsqrtD(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRsqrtIterD(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for two SIMD doubles.
 *
 * \copydetails simdInvsqrtPairF
 */
static inline void gmx_simdcall
simdInvsqrtPairD(SimdDouble x0,    SimdDouble x1,
                 SimdDouble *out0, SimdDouble *out1)
{
#if GMX_SIMD_HAVE_FLOAT && (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH) && (GMX_SIMD_RSQRT_BITS < 22)
    SimdFloat  xf  = simdCvtDD2F(x0, x1);
    SimdFloat  luf = simdRsqrtF(xf);
    SimdDouble lu0, lu1;
    /* Intermediate target is single - mantissa+1 bits */
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = simdRsqrtIterF(luf, xf);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = simdRsqrtIterF(luf, xf);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = simdRsqrtIterF(luf, xf);
#endif
    simdCvtF2DD(luf, &lu0, &lu1);
    /* Last iteration(s) performed in double - if we had 22 bits, this gets us to 44 (~1e-15) */
#if (GMX_SIMD_ACCURACY_BITS_SINGLE < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu0 = simdRsqrtIterD(lu0, x0);
    lu1 = simdRsqrtIterD(lu1, x1);
#endif
#if (GMX_SIMD_ACCURACY_BITS_SINGLE*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu0 = simdRsqrtIterD(lu0, x0);
    lu1 = simdRsqrtIterD(lu1, x1);
#endif
    *out0 = lu0;
    *out1 = lu1;
#else
    *out0 = simdInvsqrtD(x0);
    *out1 = simdInvsqrtD(x1);
#endif
}

#ifndef simdRcpIterD
/*! \brief Perform one Newton-Raphson iteration to improve 1/x for SIMD double.
 *
 * \copydetails simdRcpIterF
 */
static inline SimdDouble gmx_simdcall
simdRcpIterD(SimdDouble lu, SimdDouble x)
{
    return simdMulD(lu, simdFnmaddD(lu, x, simdSet1D(2.0)));
}
#endif

/*! \brief Calculate 1/x for SIMD double.
 *
 * \copydetails simdInvF
 */
static inline SimdDouble gmx_simdcall
simdInvD(SimdDouble x)
{
    SimdDouble lu = simdRcpD(x);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRcpIterD(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for masked entries of SIMD double.
 *
 *  This routine only evaluates 1/sqrt(x) for elements for which mask is true.
 *  Illegal values in the masked-out elements will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be >0 for masked-in entries
 *  \param m Mask
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 */
static inline SimdDouble
simdInvsqrtMaskD(SimdDouble x, SimdDBool m)
{
    SimdDouble lu = simdRsqrtMaskD(x, m);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRsqrtIterD(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/x for SIMD float, masked version.
 *
 * You should normally call the real-precision routine \ref gmx::simdInv.
 *
 *  \param x Argument that must be nonzero for non-masked entries.
 *  \param m Mask
 *  \return 1/x for elements where m is true, or 0.0 for masked-out entries.
 */
static inline SimdDouble gmx_simdcall
simdInvMaskD(SimdDouble x, SimdDBool m)
{
    SimdDouble lu = simdRcpMaskD(x, m);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simdRcpIterD(lu, x);
#endif
    return lu;
}


/*! \brief Calculate sqrt(x) correctly for SIMD doubles, including argument 0.0.
 *
 * \copydetails simdSqrtF
 */
static inline SimdDouble gmx_simdcall
simdSqrtD(SimdDouble x)
{
    SimdDouble  res;

    res  = simdInvsqrtMaskD(x, simdCmpNzD(x));
    return simdMulD(res, x);
}

/*! \brief SIMD double log(x). This is the natural logarithm.
 *
 * \copydetails simdLogF
 */
static inline SimdDouble gmx_simdcall
simdLogD(SimdDouble x)
{
    const SimdDouble  half       = simdSet1D(0.5);
    const SimdDouble  one        = simdSet1D(1.0);
    const SimdDouble  sqrt2      = simdSet1D(sqrt(2.0));
    const SimdDouble  corr       = simdSet1D(0.693147180559945286226764);
    const SimdDouble  CL15       = simdSet1D(0.148197055177935105296783);
    const SimdDouble  CL13       = simdSet1D(0.153108178020442575739679);
    const SimdDouble  CL11       = simdSet1D(0.181837339521549679055568);
    const SimdDouble  CL9        = simdSet1D(0.22222194152736701733275);
    const SimdDouble  CL7        = simdSet1D(0.285714288030134544449368);
    const SimdDouble  CL5        = simdSet1D(0.399999999989941956712869);
    const SimdDouble  CL3        = simdSet1D(0.666666666666685503450651);
    const SimdDouble  CL1        = simdSet1D(2.0);
    SimdDouble        fexp, x2, p;
    SimdDBool         mask;

    fexp  = simdGetExponentD(x);
    x     = simdGetMantissaD(x);

    mask  = simdCmpLtD(sqrt2, x);
    /* Adjust to non-IEEE format for x>sqrt(2): exponent += 1, mantissa *= 0.5 */
    fexp  = simdAddD(fexp, simdMaskD(one, mask));
    x     = simdMulD(x, simdBlendD(one, half, mask));

    x     = simdMulD( simdSubD(x, one), simdInvD( simdAddD(x, one) ) );
    x2    = simdMulD(x, x);

    p     = simdFmaddD(CL15, x2, CL13);
    p     = simdFmaddD(p, x2, CL11);
    p     = simdFmaddD(p, x2, CL9);
    p     = simdFmaddD(p, x2, CL7);
    p     = simdFmaddD(p, x2, CL5);
    p     = simdFmaddD(p, x2, CL3);
    p     = simdFmaddD(p, x2, CL1);
    p     = simdFmaddD(p, x, simdMulD(corr, fexp));

    return p;
}

/*! \brief SIMD double 2^x.
 *
 * \copydetails simdExp2F
 */
static inline SimdDouble gmx_simdcall
simdExp2D(SimdDouble x)
{
    const SimdDouble  arglimit      = simdSet1D(1022.0);
    const SimdDouble  CE11          = simdSet1D(4.435280790452730022081181e-10);
    const SimdDouble  CE10          = simdSet1D(7.074105630863314448024247e-09);
    const SimdDouble  CE9           = simdSet1D(1.017819803432096698472621e-07);
    const SimdDouble  CE8           = simdSet1D(1.321543308956718799557863e-06);
    const SimdDouble  CE7           = simdSet1D(0.00001525273348995851746990884);
    const SimdDouble  CE6           = simdSet1D(0.0001540353046251466849082632);
    const SimdDouble  CE5           = simdSet1D(0.001333355814678995257307880);
    const SimdDouble  CE4           = simdSet1D(0.009618129107588335039176502);
    const SimdDouble  CE3           = simdSet1D(0.05550410866481992147457793);
    const SimdDouble  CE2           = simdSet1D(0.2402265069591015620470894);
    const SimdDouble  CE1           = simdSet1D(0.6931471805599453304615075);
    const SimdDouble  one           = simdSet1D(1.0);
    SimdDouble        fexppart;
    SimdDouble        intpart;
    SimdDouble        p;
    SimdDBool         valuemask;

    fexppart  = simdSetExponentD(x);  /* rounds to nearest int internally */
    intpart   = simdRoundD(x);        /* use same rounding mode here */
    valuemask = simdCmpLeD(simdAbsD(x), arglimit);
    fexppart  = simdMaskD(fexppart, valuemask);
    x         = simdSubD(x, intpart);

    p         = simdFmaddD(CE11, x, CE10);
    p         = simdFmaddD(p, x, CE9);
    p         = simdFmaddD(p, x, CE8);
    p         = simdFmaddD(p, x, CE7);
    p         = simdFmaddD(p, x, CE6);
    p         = simdFmaddD(p, x, CE5);
    p         = simdFmaddD(p, x, CE4);
    p         = simdFmaddD(p, x, CE3);
    p         = simdFmaddD(p, x, CE2);
    p         = simdFmaddD(p, x, CE1);
    p         = simdFmaddD(p, x, one);
    x         = simdMulD(p, fexppart);
    return x;
}

/*! \brief SIMD double exp(x).
 *
 * \copydetails simdExpF
 */
static inline SimdDouble gmx_simdcall
simdExpD(SimdDouble x)
{
    const SimdDouble  argscale      = simdSet1D(1.44269504088896340735992468100);
    const SimdDouble  arglimit      = simdSet1D(1022.0);
    const SimdDouble  invargscale0  = simdSet1D(-0.69314718055966295651160180568695068359375);
    const SimdDouble  invargscale1  = simdSet1D(-2.8235290563031577122588448175013436025525412068e-13);
    const SimdDouble  CE12          = simdSet1D(2.078375306791423699350304e-09);
    const SimdDouble  CE11          = simdSet1D(2.518173854179933105218635e-08);
    const SimdDouble  CE10          = simdSet1D(2.755842049600488770111608e-07);
    const SimdDouble  CE9           = simdSet1D(2.755691815216689746619849e-06);
    const SimdDouble  CE8           = simdSet1D(2.480158383706245033920920e-05);
    const SimdDouble  CE7           = simdSet1D(0.0001984127043518048611841321);
    const SimdDouble  CE6           = simdSet1D(0.001388888889360258341755930);
    const SimdDouble  CE5           = simdSet1D(0.008333333332907368102819109);
    const SimdDouble  CE4           = simdSet1D(0.04166666666663836745814631);
    const SimdDouble  CE3           = simdSet1D(0.1666666666666796929434570);
    const SimdDouble  CE2           = simdSet1D(0.5);
    const SimdDouble  one           = simdSet1D(1.0);
    SimdDouble        fexppart;
    SimdDouble        intpart;
    SimdDouble        y, p;
    SimdDBool         valuemask;

    y         = simdMulD(x, argscale);
    fexppart  = simdSetExponentD(y);  /* rounds to nearest int internally */
    intpart   = simdRoundD(y);        /* use same rounding mode here */
    valuemask = simdCmpLeD(simdAbsD(y), arglimit);
    fexppart  = simdMaskD(fexppart, valuemask);

    /* Extended precision arithmetics */
    x         = simdFmaddD(invargscale0, intpart, x);
    x         = simdFmaddD(invargscale1, intpart, x);

    p         = simdFmaddD(CE12, x, CE11);
    p         = simdFmaddD(p, x, CE10);
    p         = simdFmaddD(p, x, CE9);
    p         = simdFmaddD(p, x, CE8);
    p         = simdFmaddD(p, x, CE7);
    p         = simdFmaddD(p, x, CE6);
    p         = simdFmaddD(p, x, CE5);
    p         = simdFmaddD(p, x, CE4);
    p         = simdFmaddD(p, x, CE3);
    p         = simdFmaddD(p, x, CE2);
    p         = simdFmaddD(p, simdMulD(x, x), simdAddD(x, one));
    x         = simdMulD(p, fexppart);
    return x;
}

/*! \brief SIMD double erf(x).
 *
 * \copydetails simdErfF
 */
static inline SimdDouble gmx_simdcall
simdErfD(SimdDouble x)
{
    /* Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75] */
    const SimdDouble CAP4      = simdSet1D(-0.431780540597889301512e-4);
    const SimdDouble CAP3      = simdSet1D(-0.00578562306260059236059);
    const SimdDouble CAP2      = simdSet1D(-0.028593586920219752446);
    const SimdDouble CAP1      = simdSet1D(-0.315924962948621698209);
    const SimdDouble CAP0      = simdSet1D(0.14952975608477029151);

    const SimdDouble CAQ5      = simdSet1D(-0.374089300177174709737e-5);
    const SimdDouble CAQ4      = simdSet1D(0.00015126584532155383535);
    const SimdDouble CAQ3      = simdSet1D(0.00536692680669480725423);
    const SimdDouble CAQ2      = simdSet1D(0.0668686825594046122636);
    const SimdDouble CAQ1      = simdSet1D(0.402604990869284362773);
    /* CAQ0 == 1.0 */
    const SimdDouble CAoffset  = simdSet1D(0.9788494110107421875);

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5] */
    const SimdDouble CBP6      = simdSet1D(2.49650423685462752497647637088e-10);
    const SimdDouble CBP5      = simdSet1D(0.00119770193298159629350136085658);
    const SimdDouble CBP4      = simdSet1D(0.0164944422378370965881008942733);
    const SimdDouble CBP3      = simdSet1D(0.0984581468691775932063932439252);
    const SimdDouble CBP2      = simdSet1D(0.317364595806937763843589437418);
    const SimdDouble CBP1      = simdSet1D(0.554167062641455850932670067075);
    const SimdDouble CBP0      = simdSet1D(0.427583576155807163756925301060);
    const SimdDouble CBQ7      = simdSet1D(0.00212288829699830145976198384930);
    const SimdDouble CBQ6      = simdSet1D(0.0334810979522685300554606393425);
    const SimdDouble CBQ5      = simdSet1D(0.2361713785181450957579508850717);
    const SimdDouble CBQ4      = simdSet1D(0.955364736493055670530981883072);
    const SimdDouble CBQ3      = simdSet1D(2.36815675631420037315349279199);
    const SimdDouble CBQ2      = simdSet1D(3.55261649184083035537184223542);
    const SimdDouble CBQ1      = simdSet1D(2.93501136050160872574376997993);
    /* CBQ0 == 1.0 */

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf] */
    const SimdDouble CCP6      = simdSet1D(-2.8175401114513378771);
    const SimdDouble CCP5      = simdSet1D(-3.22729451764143718517);
    const SimdDouble CCP4      = simdSet1D(-2.5518551727311523996);
    const SimdDouble CCP3      = simdSet1D(-0.687717681153649930619);
    const SimdDouble CCP2      = simdSet1D(-0.212652252872804219852);
    const SimdDouble CCP1      = simdSet1D(0.0175389834052493308818);
    const SimdDouble CCP0      = simdSet1D(0.00628057170626964891937);

    const SimdDouble CCQ6      = simdSet1D(5.48409182238641741584);
    const SimdDouble CCQ5      = simdSet1D(13.5064170191802889145);
    const SimdDouble CCQ4      = simdSet1D(22.9367376522880577224);
    const SimdDouble CCQ3      = simdSet1D(15.930646027911794143);
    const SimdDouble CCQ2      = simdSet1D(11.0567237927800161565);
    const SimdDouble CCQ1      = simdSet1D(2.79257750980575282228);
    /* CCQ0 == 1.0 */
    const SimdDouble CCoffset  = simdSet1D(0.5579090118408203125);

    const SimdDouble one       = simdSet1D(1.0);
    const SimdDouble two       = simdSet1D(2.0);

    SimdDouble       xabs, x2, x4, t, t2, w, w2;
    SimdDouble       PolyAP0, PolyAP1, PolyAQ0, PolyAQ1;
    SimdDouble       PolyBP0, PolyBP1, PolyBQ0, PolyBQ1;
    SimdDouble       PolyCP0, PolyCP1, PolyCQ0, PolyCQ1;
    SimdDouble       res_erf, res_erfcB, res_erfcC, res_erfc, res;
    SimdDouble       expmx2;
    SimdDBool        mask, mask_erf, notmask_erf;

    /* Calculate erf() */
    xabs        = simdAbsD(x);
    mask_erf    = simdCmpLtD(xabs, one);
    notmask_erf = simdCmpLeD(one, xabs);
    x2          = simdMulD(x, x);
    x4          = simdMulD(x2, x2);

    PolyAP0  = simdMulD(CAP4, x4);
    PolyAP1  = simdMulD(CAP3, x4);
    PolyAP0  = simdAddD(PolyAP0, CAP2);
    PolyAP1  = simdAddD(PolyAP1, CAP1);
    PolyAP0  = simdMulD(PolyAP0, x4);
    PolyAP1  = simdMulD(PolyAP1, x2);
    PolyAP0  = simdAddD(PolyAP0, CAP0);
    PolyAP0  = simdAddD(PolyAP0, PolyAP1);

    PolyAQ1  = simdMulD(CAQ5, x4);
    PolyAQ0  = simdMulD(CAQ4, x4);
    PolyAQ1  = simdAddD(PolyAQ1, CAQ3);
    PolyAQ0  = simdAddD(PolyAQ0, CAQ2);
    PolyAQ1  = simdMulD(PolyAQ1, x4);
    PolyAQ0  = simdMulD(PolyAQ0, x4);
    PolyAQ1  = simdAddD(PolyAQ1, CAQ1);
    PolyAQ0  = simdAddD(PolyAQ0, one);
    PolyAQ1  = simdMulD(PolyAQ1, x2);
    PolyAQ0  = simdAddD(PolyAQ0, PolyAQ1);

    res_erf  = simdMulD(PolyAP0, simdInvMaskD(PolyAQ0, mask_erf));
    res_erf  = simdAddD(CAoffset, res_erf);
    res_erf  = simdMulD(x, res_erf);

    /* Calculate erfc() in range [1,4.5] */
    t       = simdSubD(xabs, one);
    t2      = simdMulD(t, t);

    PolyBP0  = simdMulD(CBP6, t2);
    PolyBP1  = simdMulD(CBP5, t2);
    PolyBP0  = simdAddD(PolyBP0, CBP4);
    PolyBP1  = simdAddD(PolyBP1, CBP3);
    PolyBP0  = simdMulD(PolyBP0, t2);
    PolyBP1  = simdMulD(PolyBP1, t2);
    PolyBP0  = simdAddD(PolyBP0, CBP2);
    PolyBP1  = simdAddD(PolyBP1, CBP1);
    PolyBP0  = simdMulD(PolyBP0, t2);
    PolyBP1  = simdMulD(PolyBP1, t);
    PolyBP0  = simdAddD(PolyBP0, CBP0);
    PolyBP0  = simdAddD(PolyBP0, PolyBP1);

    PolyBQ1 = simdMulD(CBQ7, t2);
    PolyBQ0 = simdMulD(CBQ6, t2);
    PolyBQ1 = simdAddD(PolyBQ1, CBQ5);
    PolyBQ0 = simdAddD(PolyBQ0, CBQ4);
    PolyBQ1 = simdMulD(PolyBQ1, t2);
    PolyBQ0 = simdMulD(PolyBQ0, t2);
    PolyBQ1 = simdAddD(PolyBQ1, CBQ3);
    PolyBQ0 = simdAddD(PolyBQ0, CBQ2);
    PolyBQ1 = simdMulD(PolyBQ1, t2);
    PolyBQ0 = simdMulD(PolyBQ0, t2);
    PolyBQ1 = simdAddD(PolyBQ1, CBQ1);
    PolyBQ0 = simdAddD(PolyBQ0, one);
    PolyBQ1 = simdMulD(PolyBQ1, t);
    PolyBQ0 = simdAddD(PolyBQ0, PolyBQ1);

    /* The denominator polynomial can be zero outside the range */
    res_erfcB = simdMulD(PolyBP0, simdInvMaskD(PolyBQ0, notmask_erf));

    res_erfcB = simdMulD(res_erfcB, xabs);

    /* Calculate erfc() in range [4.5,inf] */
    w       = simdInvMaskD(xabs, notmask_erf);
    w2      = simdMulD(w, w);

    PolyCP0  = simdMulD(CCP6, w2);
    PolyCP1  = simdMulD(CCP5, w2);
    PolyCP0  = simdAddD(PolyCP0, CCP4);
    PolyCP1  = simdAddD(PolyCP1, CCP3);
    PolyCP0  = simdMulD(PolyCP0, w2);
    PolyCP1  = simdMulD(PolyCP1, w2);
    PolyCP0  = simdAddD(PolyCP0, CCP2);
    PolyCP1  = simdAddD(PolyCP1, CCP1);
    PolyCP0  = simdMulD(PolyCP0, w2);
    PolyCP1  = simdMulD(PolyCP1, w);
    PolyCP0  = simdAddD(PolyCP0, CCP0);
    PolyCP0  = simdAddD(PolyCP0, PolyCP1);

    PolyCQ0  = simdMulD(CCQ6, w2);
    PolyCQ1  = simdMulD(CCQ5, w2);
    PolyCQ0  = simdAddD(PolyCQ0, CCQ4);
    PolyCQ1  = simdAddD(PolyCQ1, CCQ3);
    PolyCQ0  = simdMulD(PolyCQ0, w2);
    PolyCQ1  = simdMulD(PolyCQ1, w2);
    PolyCQ0  = simdAddD(PolyCQ0, CCQ2);
    PolyCQ1  = simdAddD(PolyCQ1, CCQ1);
    PolyCQ0  = simdMulD(PolyCQ0, w2);
    PolyCQ1  = simdMulD(PolyCQ1, w);
    PolyCQ0  = simdAddD(PolyCQ0, one);
    PolyCQ0  = simdAddD(PolyCQ0, PolyCQ1);

    expmx2   = simdExpD( simdNegD(x2) );

    /* The denominator polynomial can be zero outside the range */
    res_erfcC = simdMulD(PolyCP0, simdInvMaskD(PolyCQ0, notmask_erf));
    res_erfcC = simdAddD(res_erfcC, CCoffset);
    res_erfcC = simdMulD(res_erfcC, w);

    mask     = simdCmpLtD(simdSet1D(4.5), xabs);
    res_erfc = simdBlendD(res_erfcB, res_erfcC, mask);

    res_erfc = simdMulD(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = simdCmpLtD(x, simdSetZeroD());
    res_erfc = simdBlendD(res_erfc, simdSubD(two, res_erfc), mask);

    /* Select erf() or erfc() */
    res  = simdBlendD(simdSubD(one, res_erfc), res_erf, mask_erf);

    return res;
}

/*! \brief SIMD double erfc(x).
 *
 * \copydetails simdErfcF
 */
static inline SimdDouble gmx_simdcall
simdErfcD(SimdDouble x)
{
    /* Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75] */
    const SimdDouble CAP4      = simdSet1D(-0.431780540597889301512e-4);
    const SimdDouble CAP3      = simdSet1D(-0.00578562306260059236059);
    const SimdDouble CAP2      = simdSet1D(-0.028593586920219752446);
    const SimdDouble CAP1      = simdSet1D(-0.315924962948621698209);
    const SimdDouble CAP0      = simdSet1D(0.14952975608477029151);

    const SimdDouble CAQ5      = simdSet1D(-0.374089300177174709737e-5);
    const SimdDouble CAQ4      = simdSet1D(0.00015126584532155383535);
    const SimdDouble CAQ3      = simdSet1D(0.00536692680669480725423);
    const SimdDouble CAQ2      = simdSet1D(0.0668686825594046122636);
    const SimdDouble CAQ1      = simdSet1D(0.402604990869284362773);
    /* CAQ0 == 1.0 */
    const SimdDouble CAoffset  = simdSet1D(0.9788494110107421875);

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5] */
    const SimdDouble CBP6      = simdSet1D(2.49650423685462752497647637088e-10);
    const SimdDouble CBP5      = simdSet1D(0.00119770193298159629350136085658);
    const SimdDouble CBP4      = simdSet1D(0.0164944422378370965881008942733);
    const SimdDouble CBP3      = simdSet1D(0.0984581468691775932063932439252);
    const SimdDouble CBP2      = simdSet1D(0.317364595806937763843589437418);
    const SimdDouble CBP1      = simdSet1D(0.554167062641455850932670067075);
    const SimdDouble CBP0      = simdSet1D(0.427583576155807163756925301060);
    const SimdDouble CBQ7      = simdSet1D(0.00212288829699830145976198384930);
    const SimdDouble CBQ6      = simdSet1D(0.0334810979522685300554606393425);
    const SimdDouble CBQ5      = simdSet1D(0.2361713785181450957579508850717);
    const SimdDouble CBQ4      = simdSet1D(0.955364736493055670530981883072);
    const SimdDouble CBQ3      = simdSet1D(2.36815675631420037315349279199);
    const SimdDouble CBQ2      = simdSet1D(3.55261649184083035537184223542);
    const SimdDouble CBQ1      = simdSet1D(2.93501136050160872574376997993);
    /* CBQ0 == 1.0 */

    /* Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf] */
    const SimdDouble CCP6      = simdSet1D(-2.8175401114513378771);
    const SimdDouble CCP5      = simdSet1D(-3.22729451764143718517);
    const SimdDouble CCP4      = simdSet1D(-2.5518551727311523996);
    const SimdDouble CCP3      = simdSet1D(-0.687717681153649930619);
    const SimdDouble CCP2      = simdSet1D(-0.212652252872804219852);
    const SimdDouble CCP1      = simdSet1D(0.0175389834052493308818);
    const SimdDouble CCP0      = simdSet1D(0.00628057170626964891937);

    const SimdDouble CCQ6      = simdSet1D(5.48409182238641741584);
    const SimdDouble CCQ5      = simdSet1D(13.5064170191802889145);
    const SimdDouble CCQ4      = simdSet1D(22.9367376522880577224);
    const SimdDouble CCQ3      = simdSet1D(15.930646027911794143);
    const SimdDouble CCQ2      = simdSet1D(11.0567237927800161565);
    const SimdDouble CCQ1      = simdSet1D(2.79257750980575282228);
    /* CCQ0 == 1.0 */
    const SimdDouble CCoffset  = simdSet1D(0.5579090118408203125);

    const SimdDouble one       = simdSet1D(1.0);
    const SimdDouble two       = simdSet1D(2.0);

    SimdDouble       xabs, x2, x4, t, t2, w, w2;
    SimdDouble       PolyAP0, PolyAP1, PolyAQ0, PolyAQ1;
    SimdDouble       PolyBP0, PolyBP1, PolyBQ0, PolyBQ1;
    SimdDouble       PolyCP0, PolyCP1, PolyCQ0, PolyCQ1;
    SimdDouble       res_erf, res_erfcB, res_erfcC, res_erfc, res;
    SimdDouble       expmx2;
    SimdDBool        mask, mask_erf, notmask_erf;

    /* Calculate erf() */
    xabs        = simdAbsD(x);
    mask_erf    = simdCmpLtD(xabs, one);
    notmask_erf = simdCmpLeD(one, xabs);
    x2          = simdMulD(x, x);
    x4          = simdMulD(x2, x2);

    PolyAP0  = simdMulD(CAP4, x4);
    PolyAP1  = simdMulD(CAP3, x4);
    PolyAP0  = simdAddD(PolyAP0, CAP2);
    PolyAP1  = simdAddD(PolyAP1, CAP1);
    PolyAP0  = simdMulD(PolyAP0, x4);
    PolyAP1  = simdMulD(PolyAP1, x2);
    PolyAP0  = simdAddD(PolyAP0, CAP0);
    PolyAP0  = simdAddD(PolyAP0, PolyAP1);

    PolyAQ1  = simdMulD(CAQ5, x4);
    PolyAQ0  = simdMulD(CAQ4, x4);
    PolyAQ1  = simdAddD(PolyAQ1, CAQ3);
    PolyAQ0  = simdAddD(PolyAQ0, CAQ2);
    PolyAQ1  = simdMulD(PolyAQ1, x4);
    PolyAQ0  = simdMulD(PolyAQ0, x4);
    PolyAQ1  = simdAddD(PolyAQ1, CAQ1);
    PolyAQ0  = simdAddD(PolyAQ0, one);
    PolyAQ1  = simdMulD(PolyAQ1, x2);
    PolyAQ0  = simdAddD(PolyAQ0, PolyAQ1);

    res_erf  = simdMulD(PolyAP0, simdInvMaskD(PolyAQ0, mask_erf));
    res_erf  = simdAddD(CAoffset, res_erf);
    res_erf  = simdMulD(x, res_erf);

    /* Calculate erfc() in range [1,4.5] */
    t       = simdSubD(xabs, one);
    t2      = simdMulD(t, t);

    PolyBP0  = simdMulD(CBP6, t2);
    PolyBP1  = simdMulD(CBP5, t2);
    PolyBP0  = simdAddD(PolyBP0, CBP4);
    PolyBP1  = simdAddD(PolyBP1, CBP3);
    PolyBP0  = simdMulD(PolyBP0, t2);
    PolyBP1  = simdMulD(PolyBP1, t2);
    PolyBP0  = simdAddD(PolyBP0, CBP2);
    PolyBP1  = simdAddD(PolyBP1, CBP1);
    PolyBP0  = simdMulD(PolyBP0, t2);
    PolyBP1  = simdMulD(PolyBP1, t);
    PolyBP0  = simdAddD(PolyBP0, CBP0);
    PolyBP0  = simdAddD(PolyBP0, PolyBP1);

    PolyBQ1 = simdMulD(CBQ7, t2);
    PolyBQ0 = simdMulD(CBQ6, t2);
    PolyBQ1 = simdAddD(PolyBQ1, CBQ5);
    PolyBQ0 = simdAddD(PolyBQ0, CBQ4);
    PolyBQ1 = simdMulD(PolyBQ1, t2);
    PolyBQ0 = simdMulD(PolyBQ0, t2);
    PolyBQ1 = simdAddD(PolyBQ1, CBQ3);
    PolyBQ0 = simdAddD(PolyBQ0, CBQ2);
    PolyBQ1 = simdMulD(PolyBQ1, t2);
    PolyBQ0 = simdMulD(PolyBQ0, t2);
    PolyBQ1 = simdAddD(PolyBQ1, CBQ1);
    PolyBQ0 = simdAddD(PolyBQ0, one);
    PolyBQ1 = simdMulD(PolyBQ1, t);
    PolyBQ0 = simdAddD(PolyBQ0, PolyBQ1);

    /* The denominator polynomial can be zero outside the range */
    res_erfcB = simdMulD(PolyBP0, simdInvMaskD(PolyBQ0, notmask_erf));

    res_erfcB = simdMulD(res_erfcB, xabs);

    /* Calculate erfc() in range [4.5,inf] */
    w       = simdInvMaskD(xabs, simdCmpNzD(xabs));
    w2      = simdMulD(w, w);

    PolyCP0  = simdMulD(CCP6, w2);
    PolyCP1  = simdMulD(CCP5, w2);
    PolyCP0  = simdAddD(PolyCP0, CCP4);
    PolyCP1  = simdAddD(PolyCP1, CCP3);
    PolyCP0  = simdMulD(PolyCP0, w2);
    PolyCP1  = simdMulD(PolyCP1, w2);
    PolyCP0  = simdAddD(PolyCP0, CCP2);
    PolyCP1  = simdAddD(PolyCP1, CCP1);
    PolyCP0  = simdMulD(PolyCP0, w2);
    PolyCP1  = simdMulD(PolyCP1, w);
    PolyCP0  = simdAddD(PolyCP0, CCP0);
    PolyCP0  = simdAddD(PolyCP0, PolyCP1);

    PolyCQ0  = simdMulD(CCQ6, w2);
    PolyCQ1  = simdMulD(CCQ5, w2);
    PolyCQ0  = simdAddD(PolyCQ0, CCQ4);
    PolyCQ1  = simdAddD(PolyCQ1, CCQ3);
    PolyCQ0  = simdMulD(PolyCQ0, w2);
    PolyCQ1  = simdMulD(PolyCQ1, w2);
    PolyCQ0  = simdAddD(PolyCQ0, CCQ2);
    PolyCQ1  = simdAddD(PolyCQ1, CCQ1);
    PolyCQ0  = simdMulD(PolyCQ0, w2);
    PolyCQ1  = simdMulD(PolyCQ1, w);
    PolyCQ0  = simdAddD(PolyCQ0, one);
    PolyCQ0  = simdAddD(PolyCQ0, PolyCQ1);

    expmx2   = simdExpD( simdNegD(x2) );

    /* The denominator polynomial can be zero outside the range */
    res_erfcC = simdMulD(PolyCP0, simdInvMaskD(PolyCQ0, notmask_erf));
    res_erfcC = simdAddD(res_erfcC, CCoffset);
    res_erfcC = simdMulD(res_erfcC, w);

    mask     = simdCmpLtD(simdSet1D(4.5), xabs);
    res_erfc = simdBlendD(res_erfcB, res_erfcC, mask);

    res_erfc = simdMulD(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = simdCmpLtD(x, simdSetZeroD());
    res_erfc = simdBlendD(res_erfc, simdSubD(two, res_erfc), mask);

    /* Select erf() or erfc() */
    res  = simdBlendD(res_erfc, simdSubD(one, res_erf), mask_erf);

    return res;
}

/*! \brief SIMD double sin \& cos.
 *
 * \copydetails simdSinCosF
 */
static inline void gmx_simdcall
simdSinCosD(SimdDouble x, SimdDouble *sinval, SimdDouble *cosval)
{
    /* Constants to subtract Pi/4*x from y while minimizing precision loss */
    const SimdDouble  argred0         = simdSet1D(-2*0.78539816290140151978);
    const SimdDouble  argred1         = simdSet1D(-2*4.9604678871439933374e-10);
    const SimdDouble  argred2         = simdSet1D(-2*1.1258708853173288931e-18);
    const SimdDouble  argred3         = simdSet1D(-2*1.7607799325916000908e-27);
    const SimdDouble  two_over_pi     = simdSet1D(2.0/M_PI);
    const SimdDouble  const_sin5      = simdSet1D( 1.58938307283228937328511e-10);
    const SimdDouble  const_sin4      = simdSet1D(-2.50506943502539773349318e-08);
    const SimdDouble  const_sin3      = simdSet1D( 2.75573131776846360512547e-06);
    const SimdDouble  const_sin2      = simdSet1D(-0.000198412698278911770864914);
    const SimdDouble  const_sin1      = simdSet1D( 0.0083333333333191845961746);
    const SimdDouble  const_sin0      = simdSet1D(-0.166666666666666130709393);

    const SimdDouble  const_cos7      = simdSet1D(-1.13615350239097429531523e-11);
    const SimdDouble  const_cos6      = simdSet1D( 2.08757471207040055479366e-09);
    const SimdDouble  const_cos5      = simdSet1D(-2.75573144028847567498567e-07);
    const SimdDouble  const_cos4      = simdSet1D( 2.48015872890001867311915e-05);
    const SimdDouble  const_cos3      = simdSet1D(-0.00138888888888714019282329);
    const SimdDouble  const_cos2      = simdSet1D( 0.0416666666666665519592062);
    const SimdDouble  half            = simdSet1D(0.5);
    const SimdDouble  one             = simdSet1D(1.0);
    SimdDouble        ssign, csign;
    SimdDouble        x2, y, z, psin, pcos, sss, ccc;
    SimdDBool         mask;
#if GMX_SIMD_HAVE_DINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    const SimdDInt32  ione            = simdSet1DI(1);
    const SimdDInt32  itwo            = simdSet1DI(2);
    SimdDInt32        iy;

    z       = simdMulD(x, two_over_pi);
    iy      = simdCvtD2I(z);
    y       = simdRoundD(z);

    mask    = simdCvtDIB2DB(simdCmpEqDI(simdAndDI(iy, ione), simdSetZeroDI()));
    ssign   = simdMaskD(simdSet1D(GMX_DOUBLE_NEGZERO), simdCvtDIB2DB(simdCmpEqDI(simdAndDI(iy, itwo), itwo)));
    csign   = simdMaskD(simdSet1D(GMX_DOUBLE_NEGZERO), simdCvtDIB2DB(simdCmpEqDI(simdAndDI(simdAddDI(iy, ione), itwo), itwo)));
#else
    const SimdDouble  quarter         = simdSet1D(0.25);
    const SimdDouble  minusquarter    = simdSet1D(-0.25);
    SimdDouble        q;
    SimdDBool         m1, m2, m3;

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
    x       = simdAbsD(x);
    /* It is critical that half-way cases are rounded down */
    z       = simdFmaddD(x, two_over_pi, half);
    y       = simdTruncD(z);
    q       = simdMulD(z, quarter);
    q       = simdSubD(q, simdTruncD(q));
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q       = simdSubD(q, half);
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#    if GMX_SIMD_HAVE_LOGICAL
    ssign   = simdAndD(ssign, simdSet1D(GMX_DOUBLE_NEGZERO));
    csign   = simdAndNotD(q, simdSet1D(GMX_DOUBLE_NEGZERO));
    ssign   = simdXorD(ssign, csign);
#    else
    csign   = simdXorSignD(simdSet1D(-1.0), q);
    ssign   = simdXorSignD(ssign, csign);    /* swap ssign if csign was set. */
#    endif
    /* Check if y had value 1 or 3 (remember we subtracted 0.5 from q) */
    m1      = simdCmpLtD(q, minusquarter);
    m2      = simdCmpLeD(simdSetZeroD(), q);
    m3      = simdCmpLtD(q, quarter);
    m2      = simdAndDB(m2, m3);
    mask    = simdOrDB(m1, m2);
    /* where mask is FALSE, set sign. */
    csign   = simdXorSignD(csign, simdBlendD(simdSet1D(-1.0), one, mask));
#endif
    x       = simdFmaddD(y, argred0, x);
    x       = simdFmaddD(y, argred1, x);
    x       = simdFmaddD(y, argred2, x);
    x       = simdFmaddD(y, argred3, x);
    x2      = simdMulD(x, x);

    psin    = simdFmaddD(const_sin5, x2, const_sin4);
    psin    = simdFmaddD(psin, x2, const_sin3);
    psin    = simdFmaddD(psin, x2, const_sin2);
    psin    = simdFmaddD(psin, x2, const_sin1);
    psin    = simdFmaddD(psin, x2, const_sin0);
    psin    = simdFmaddD(psin, simdMulD(x2, x), x);

    pcos    = simdFmaddD(const_cos7, x2, const_cos6);
    pcos    = simdFmaddD(pcos, x2, const_cos5);
    pcos    = simdFmaddD(pcos, x2, const_cos4);
    pcos    = simdFmaddD(pcos, x2, const_cos3);
    pcos    = simdFmaddD(pcos, x2, const_cos2);
    pcos    = simdFmsubD(pcos, x2, half);
    pcos    = simdFmaddD(pcos, x2, one);

    sss     = simdBlendD(pcos, psin, mask);
    ccc     = simdBlendD(psin, pcos, mask);
    /* See comment for GMX_SIMD_HAVE_LOGICAL section above. */
#if GMX_SIMD_HAVE_LOGICAL
    *sinval = simdXorD(sss, ssign);
    *cosval = simdXorD(ccc, csign);
#else
    *sinval = simdXorSignD(sss, ssign);
    *cosval = simdXorSignD(ccc, csign);
#endif
}

/*! \brief SIMD double sin(x).
 *
 * \copydetails simdSinF
 */
static inline SimdDouble gmx_simdcall
simdSinD(SimdDouble x)
{
    SimdDouble s, c;
    simdSinCosD(x, &s, &c);
    return s;
}

/*! \brief SIMD double cos(x).
 *
 * \copydetails simdCosF
 */
static inline SimdDouble gmx_simdcall
simdCosD(SimdDouble x)
{
    SimdDouble s, c;
    simdSinCosD(x, &s, &c);
    return c;
}

/*! \brief SIMD double tan(x).
 *
 * \copydetails simdTanF
 */
static inline SimdDouble gmx_simdcall
simdTanD(SimdDouble x)
{
    const SimdDouble  argred0         = simdSet1D(-2*0.78539816290140151978);
    const SimdDouble  argred1         = simdSet1D(-2*4.9604678871439933374e-10);
    const SimdDouble  argred2         = simdSet1D(-2*1.1258708853173288931e-18);
    const SimdDouble  argred3         = simdSet1D(-2*1.7607799325916000908e-27);
    const SimdDouble  two_over_pi     = simdSet1D(2.0/M_PI);
    const SimdDouble  CT15            = simdSet1D(1.01419718511083373224408e-05);
    const SimdDouble  CT14            = simdSet1D(-2.59519791585924697698614e-05);
    const SimdDouble  CT13            = simdSet1D(5.23388081915899855325186e-05);
    const SimdDouble  CT12            = simdSet1D(-3.05033014433946488225616e-05);
    const SimdDouble  CT11            = simdSet1D(7.14707504084242744267497e-05);
    const SimdDouble  CT10            = simdSet1D(8.09674518280159187045078e-05);
    const SimdDouble  CT9             = simdSet1D(0.000244884931879331847054404);
    const SimdDouble  CT8             = simdSet1D(0.000588505168743587154904506);
    const SimdDouble  CT7             = simdSet1D(0.00145612788922812427978848);
    const SimdDouble  CT6             = simdSet1D(0.00359208743836906619142924);
    const SimdDouble  CT5             = simdSet1D(0.00886323944362401618113356);
    const SimdDouble  CT4             = simdSet1D(0.0218694882853846389592078);
    const SimdDouble  CT3             = simdSet1D(0.0539682539781298417636002);
    const SimdDouble  CT2             = simdSet1D(0.133333333333125941821962);
    const SimdDouble  CT1             = simdSet1D(0.333333333333334980164153);

    SimdDouble        x2, p, y, z;
    SimdDBool         mask;

#if GMX_SIMD_HAVE_DINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    SimdDInt32  iy;
    SimdDInt32  ione = simdSet1DI(1);

    z       = simdMulD(x, two_over_pi);
    iy      = simdCvtD2I(z);
    y       = simdRoundD(z);
    mask    = simdCvtDIB2DB(simdCmpEqDI(simdAndDI(iy, ione), ione));

    x       = simdFmaddD(y, argred0, x);
    x       = simdFmaddD(y, argred1, x);
    x       = simdFmaddD(y, argred2, x);
    x       = simdFmaddD(y, argred3, x);
    x       = simdXorD(simdMaskD(simdSet1D(GMX_DOUBLE_NEGZERO), mask), x);
#else
    const SimdDouble  quarter         = simdSet1D(0.25);
    const SimdDouble  half            = simdSet1D(0.5);
    const SimdDouble  threequarter    = simdSet1D(0.75);
    SimdDouble        w, q;
    SimdDBool         m1, m2, m3;

    w       = simdAbsD(x);
    z       = simdFmaddD(w, two_over_pi, half);
    y       = simdTruncD(z);
    q       = simdMulD(z, quarter);
    q       = simdSubD(q, simdTruncD(q));
    m1      = simdCmpLeD(quarter, q);
    m2      = simdCmpLtD(q, half);
    m3      = simdCmpLeD(threequarter, q);
    m1      = simdAndDB(m1, m2);
    mask    = simdOrDB(m1, m3);
    w       = simdFmaddD(y, argred0, w);
    w       = simdFmaddD(y, argred1, w);
    w       = simdFmaddD(y, argred2, w);
    w       = simdFmaddD(y, argred3, w);

    w       = simdBlendD(w, simdNegD(w), mask);
    x       = simdXorSignD(w, x);
#endif
    x2      = simdMulD(x, x);
    p       = simdFmaddD(CT15, x2, CT14);
    p       = simdFmaddD(p, x2, CT13);
    p       = simdFmaddD(p, x2, CT12);
    p       = simdFmaddD(p, x2, CT11);
    p       = simdFmaddD(p, x2, CT10);
    p       = simdFmaddD(p, x2, CT9);
    p       = simdFmaddD(p, x2, CT8);
    p       = simdFmaddD(p, x2, CT7);
    p       = simdFmaddD(p, x2, CT6);
    p       = simdFmaddD(p, x2, CT5);
    p       = simdFmaddD(p, x2, CT4);
    p       = simdFmaddD(p, x2, CT3);
    p       = simdFmaddD(p, x2, CT2);
    p       = simdFmaddD(p, x2, CT1);
    p       = simdFmaddD(x2, simdMulD(p, x), x);

    p       = simdBlendD( p, simdInvMaskD(p, mask), mask);
    return p;
}

/*! \brief SIMD double asin(x).
 *
 * \copydetails simdAsinF
 */
static inline SimdDouble gmx_simdcall
simdAsinD(SimdDouble x)
{
    /* Same algorithm as cephes library */
    const SimdDouble limit1    = simdSet1D(0.625);
    const SimdDouble limit2    = simdSet1D(1e-8);
    const SimdDouble one       = simdSet1D(1.0);
    const SimdDouble quarterpi = simdSet1D(M_PI/4.0);
    const SimdDouble morebits  = simdSet1D(6.123233995736765886130e-17);

    const SimdDouble P5        = simdSet1D(4.253011369004428248960e-3);
    const SimdDouble P4        = simdSet1D(-6.019598008014123785661e-1);
    const SimdDouble P3        = simdSet1D(5.444622390564711410273e0);
    const SimdDouble P2        = simdSet1D(-1.626247967210700244449e1);
    const SimdDouble P1        = simdSet1D(1.956261983317594739197e1);
    const SimdDouble P0        = simdSet1D(-8.198089802484824371615e0);

    const SimdDouble Q4        = simdSet1D(-1.474091372988853791896e1);
    const SimdDouble Q3        = simdSet1D(7.049610280856842141659e1);
    const SimdDouble Q2        = simdSet1D(-1.471791292232726029859e2);
    const SimdDouble Q1        = simdSet1D(1.395105614657485689735e2);
    const SimdDouble Q0        = simdSet1D(-4.918853881490881290097e1);

    const SimdDouble R4        = simdSet1D(2.967721961301243206100e-3);
    const SimdDouble R3        = simdSet1D(-5.634242780008963776856e-1);
    const SimdDouble R2        = simdSet1D(6.968710824104713396794e0);
    const SimdDouble R1        = simdSet1D(-2.556901049652824852289e1);
    const SimdDouble R0        = simdSet1D(2.853665548261061424989e1);

    const SimdDouble S3        = simdSet1D(-2.194779531642920639778e1);
    const SimdDouble S2        = simdSet1D(1.470656354026814941758e2);
    const SimdDouble S1        = simdSet1D(-3.838770957603691357202e2);
    const SimdDouble S0        = simdSet1D(3.424398657913078477438e2);

    SimdDouble       xabs;
    SimdDouble       zz, ww, z, q, w, zz2, ww2;
    SimdDouble       PA, PB;
    SimdDouble       QA, QB;
    SimdDouble       RA, RB;
    SimdDouble       SA, SB;
    SimdDouble       nom, denom;
    SimdDBool        mask, mask2;

    xabs  = simdAbsD(x);

    mask  = simdCmpLtD(limit1, xabs);

    zz    = simdSubD(one, xabs);
    ww    = simdMulD(xabs, xabs);
    zz2   = simdMulD(zz, zz);
    ww2   = simdMulD(ww, ww);

    /* R */
    RA    = simdMulD(R4, zz2);
    RB    = simdMulD(R3, zz2);
    RA    = simdAddD(RA, R2);
    RB    = simdAddD(RB, R1);
    RA    = simdMulD(RA, zz2);
    RB    = simdMulD(RB, zz);
    RA    = simdAddD(RA, R0);
    RA    = simdAddD(RA, RB);

    /* S, SA = zz2 */
    SB    = simdMulD(S3, zz2);
    SA    = simdAddD(zz2, S2);
    SB    = simdAddD(SB, S1);
    SA    = simdMulD(SA, zz2);
    SB    = simdMulD(SB, zz);
    SA    = simdAddD(SA, S0);
    SA    = simdAddD(SA, SB);

    /* P */
    PA    = simdMulD(P5, ww2);
    PB    = simdMulD(P4, ww2);
    PA    = simdAddD(PA, P3);
    PB    = simdAddD(PB, P2);
    PA    = simdMulD(PA, ww2);
    PB    = simdMulD(PB, ww2);
    PA    = simdAddD(PA, P1);
    PB    = simdAddD(PB, P0);
    PA    = simdMulD(PA, ww);
    PA    = simdAddD(PA, PB);

    /* Q, QA = ww2 */
    QB    = simdMulD(Q4, ww2);
    QA    = simdAddD(ww2, Q3);
    QB    = simdAddD(QB, Q2);
    QA    = simdMulD(QA, ww2);
    QB    = simdMulD(QB, ww2);
    QA    = simdAddD(QA, Q1);
    QB    = simdAddD(QB, Q0);
    QA    = simdMulD(QA, ww);
    QA    = simdAddD(QA, QB);

    RA    = simdMulD(RA, zz);
    PA    = simdMulD(PA, ww);

    nom   = simdBlendD( PA, RA, mask );
    denom = simdBlendD( QA, SA, mask );

    mask2 = simdCmpLtD(limit2, xabs);
    q     = simdMulD( nom, simdInvMaskD(denom, mask2) );

    zz    = simdAddD(zz, zz);
    zz    = simdSqrtD(zz);
    z     = simdSubD(quarterpi, zz);
    zz    = simdMulD(zz, q);
    zz    = simdSubD(zz, morebits);
    z     = simdSubD(z, zz);
    z     = simdAddD(z, quarterpi);

    w     = simdMulD(xabs, q);
    w     = simdAddD(w, xabs);

    z     = simdBlendD( w, z, mask );

    z     = simdBlendD( xabs, z, mask2 );

    z = simdXorSignD(z, x);

    return z;
}

/*! \brief SIMD double acos(x).
 *
 * \copydetails simdAcosF
 */
static inline SimdDouble gmx_simdcall
simdAcosD(SimdDouble x)
{
    const SimdDouble one        = simdSet1D(1.0);
    const SimdDouble half       = simdSet1D(0.5);
    const SimdDouble quarterpi0 = simdSet1D(7.85398163397448309616e-1);
    const SimdDouble quarterpi1 = simdSet1D(6.123233995736765886130e-17);

    SimdDBool        mask1;
    SimdDouble       z, z1, z2;

    mask1 = simdCmpLtD(half, x);
    z1    = simdMulD(half, simdSubD(one, x));
    z1    = simdSqrtD(z1);
    z     = simdBlendD( x, z1, mask1 );

    z     = simdAsinD(z);

    z1    = simdAddD(z, z);

    z2    = simdSubD(quarterpi0, z);
    z2    = simdAddD(z2, quarterpi1);
    z2    = simdAddD(z2, quarterpi0);

    z     = simdBlendD(z2, z1, mask1);

    return z;
}

/*! \brief SIMD double atan(x).
 *
 * \copydetails simdAtanF
 */
static inline SimdDouble gmx_simdcall
simdAtanD(SimdDouble x)
{
    /* Same algorithm as cephes library */
    const SimdDouble limit1    = simdSet1D(0.66);
    const SimdDouble limit2    = simdSet1D(2.41421356237309504880);
    const SimdDouble quarterpi = simdSet1D(M_PI/4.0);
    const SimdDouble halfpi    = simdSet1D(M_PI/2.0);
    const SimdDouble mone      = simdSet1D(-1.0);
    const SimdDouble morebits1 = simdSet1D(0.5*6.123233995736765886130E-17);
    const SimdDouble morebits2 = simdSet1D(6.123233995736765886130E-17);

    const SimdDouble P4        = simdSet1D(-8.750608600031904122785E-1);
    const SimdDouble P3        = simdSet1D(-1.615753718733365076637E1);
    const SimdDouble P2        = simdSet1D(-7.500855792314704667340E1);
    const SimdDouble P1        = simdSet1D(-1.228866684490136173410E2);
    const SimdDouble P0        = simdSet1D(-6.485021904942025371773E1);

    const SimdDouble Q4        = simdSet1D(2.485846490142306297962E1);
    const SimdDouble Q3        = simdSet1D(1.650270098316988542046E2);
    const SimdDouble Q2        = simdSet1D(4.328810604912902668951E2);
    const SimdDouble Q1        = simdSet1D(4.853903996359136964868E2);
    const SimdDouble Q0        = simdSet1D(1.945506571482613964425E2);

    SimdDouble       y, xabs, t1, t2;
    SimdDouble       z, z2;
    SimdDouble       P_A, P_B, Q_A, Q_B;
    SimdDBool        mask1, mask2;

    xabs   = simdAbsD(x);

    mask1  = simdCmpLtD(limit1, xabs);
    mask2  = simdCmpLtD(limit2, xabs);

    t1     = simdMulD(simdAddD(xabs, mone),
                      simdInvMaskD(simdSubD(xabs, mone), mask1));
    t2     = simdMulD(mone, simdInvMaskD(xabs, mask2));

    y      = simdMaskD(quarterpi, mask1);
    y      = simdBlendD(y, halfpi, mask2);
    xabs   = simdBlendD(xabs, t1, mask1);
    xabs   = simdBlendD(xabs, t2, mask2);

    z      = simdMulD(xabs, xabs);
    z2     = simdMulD(z, z);

    P_A    = simdMulD(P4, z2);
    P_B    = simdMulD(P3, z2);
    P_A    = simdAddD(P_A, P2);
    P_B    = simdAddD(P_B, P1);
    P_A    = simdMulD(P_A, z2);
    P_B    = simdMulD(P_B, z);
    P_A    = simdAddD(P_A, P0);
    P_A    = simdAddD(P_A, P_B);

    /* Q_A = z2 */
    Q_B    = simdMulD(Q4, z2);
    Q_A    = simdAddD(z2, Q3);
    Q_B    = simdAddD(Q_B, Q2);
    Q_A    = simdMulD(Q_A, z2);
    Q_B    = simdMulD(Q_B, z2);
    Q_A    = simdAddD(Q_A, Q1);
    Q_B    = simdAddD(Q_B, Q0);
    Q_A    = simdMulD(Q_A, z);
    Q_A    = simdAddD(Q_A, Q_B);

    z      = simdMulD(z, P_A);
    z      = simdMulD(z, simdInvD(Q_A));
    z      = simdMulD(z, xabs);
    z      = simdAddD(z, xabs);

    t1     = simdMaskD(morebits1, mask1);
    t1     = simdBlendD(t1, morebits2, mask2);

    z      = simdAddD(z, t1);
    y      = simdAddD(y, z);

    y      = simdXorSignD(y, x);

    return y;
}

/*! \brief SIMD double atan2(y,x).
 *
 * \copydetails simdAtan2F
 */
static inline SimdDouble gmx_simdcall
simdAtan2D(SimdDouble y, SimdDouble x)
{
    const SimdDouble pi          = simdSet1D(M_PI);
    const SimdDouble halfpi      = simdSet1D(M_PI/2.0);
    SimdDouble       xinv, p, aoffset;
    SimdDBool        mask_xnz, mask_ynz, mask_xlt0, mask_ylt0;

    mask_xnz  = simdCmpNzD(x);
    mask_ynz  = simdCmpNzD(y);
    mask_xlt0 = simdCmpLtD(x, simdSetZeroD());
    mask_ylt0 = simdCmpLtD(y, simdSetZeroD());

    aoffset   = simdMaskNotD(halfpi, mask_xnz);
    aoffset   = simdMaskD(aoffset, mask_ynz);

    aoffset   = simdBlendD(aoffset, pi, mask_xlt0);
    aoffset   = simdBlendD(aoffset, simdNegD(aoffset), mask_ylt0);

    xinv      = simdInvMaskD(x, mask_xnz);
    p         = simdMulD(y, xinv);
    p         = simdAtanD(p);
    p         = simdAddD(p, aoffset);

    return p;
}


/*! \brief Calculate the force correction due to PME analytically for SIMD double.
 *
 * \copydetails simdPmeCorrForceF
 */
static inline SimdDouble gmx_simdcall
simdPmeCorrForceD(SimdDouble z2)
{
    const SimdDouble  FN10     = simdSet1D(-8.0072854618360083154e-14);
    const SimdDouble  FN9      = simdSet1D(1.1859116242260148027e-11);
    const SimdDouble  FN8      = simdSet1D(-8.1490406329798423616e-10);
    const SimdDouble  FN7      = simdSet1D(3.4404793543907847655e-8);
    const SimdDouble  FN6      = simdSet1D(-9.9471420832602741006e-7);
    const SimdDouble  FN5      = simdSet1D(0.000020740315999115847456);
    const SimdDouble  FN4      = simdSet1D(-0.00031991745139313364005);
    const SimdDouble  FN3      = simdSet1D(0.0035074449373659008203);
    const SimdDouble  FN2      = simdSet1D(-0.031750380176100813405);
    const SimdDouble  FN1      = simdSet1D(0.13884101728898463426);
    const SimdDouble  FN0      = simdSet1D(-0.75225277815249618847);

    const SimdDouble  FD5      = simdSet1D(0.000016009278224355026701);
    const SimdDouble  FD4      = simdSet1D(0.00051055686934806966046);
    const SimdDouble  FD3      = simdSet1D(0.0081803507497974289008);
    const SimdDouble  FD2      = simdSet1D(0.077181146026670287235);
    const SimdDouble  FD1      = simdSet1D(0.41543303143712535988);
    const SimdDouble  FD0      = simdSet1D(1.0);

    SimdDouble        z4;
    SimdDouble        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = simdMulD(z2, z2);

    polyFD1        = simdFmaddD(FD5, z4, FD3);
    polyFD1        = simdFmaddD(polyFD1, z4, FD1);
    polyFD1        = simdMulD(polyFD1, z2);
    polyFD0        = simdFmaddD(FD4, z4, FD2);
    polyFD0        = simdFmaddD(polyFD0, z4, FD0);
    polyFD0        = simdAddD(polyFD0, polyFD1);

    polyFD0        = simdInvD(polyFD0);

    polyFN0        = simdFmaddD(FN10, z4, FN8);
    polyFN0        = simdFmaddD(polyFN0, z4, FN6);
    polyFN0        = simdFmaddD(polyFN0, z4, FN4);
    polyFN0        = simdFmaddD(polyFN0, z4, FN2);
    polyFN0        = simdFmaddD(polyFN0, z4, FN0);
    polyFN1        = simdFmaddD(FN9, z4, FN7);
    polyFN1        = simdFmaddD(polyFN1, z4, FN5);
    polyFN1        = simdFmaddD(polyFN1, z4, FN3);
    polyFN1        = simdFmaddD(polyFN1, z4, FN1);
    polyFN0        = simdFmaddD(polyFN1, z2, polyFN0);


    return simdMulD(polyFN0, polyFD0);
}



/*! \brief Calculate the potential correction due to PME analytically for SIMD double.
 *
 * \copydetails simdPmeCorrPotentialF
 */
static inline SimdDouble gmx_simdcall
simdPmeCorrPotentialD(SimdDouble z2)
{
    const SimdDouble  VN9      = simdSet1D(-9.3723776169321855475e-13);
    const SimdDouble  VN8      = simdSet1D(1.2280156762674215741e-10);
    const SimdDouble  VN7      = simdSet1D(-7.3562157912251309487e-9);
    const SimdDouble  VN6      = simdSet1D(2.6215886208032517509e-7);
    const SimdDouble  VN5      = simdSet1D(-4.9532491651265819499e-6);
    const SimdDouble  VN4      = simdSet1D(0.00025907400778966060389);
    const SimdDouble  VN3      = simdSet1D(0.0010585044856156469792);
    const SimdDouble  VN2      = simdSet1D(0.045247661136833092885);
    const SimdDouble  VN1      = simdSet1D(0.11643931522926034421);
    const SimdDouble  VN0      = simdSet1D(1.1283791671726767970);

    const SimdDouble  VD5      = simdSet1D(0.000021784709867336150342);
    const SimdDouble  VD4      = simdSet1D(0.00064293662010911388448);
    const SimdDouble  VD3      = simdSet1D(0.0096311444822588683504);
    const SimdDouble  VD2      = simdSet1D(0.085608012351550627051);
    const SimdDouble  VD1      = simdSet1D(0.43652499166614811084);
    const SimdDouble  VD0      = simdSet1D(1.0);

    SimdDouble        z4;
    SimdDouble        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = simdMulD(z2, z2);

    polyVD1        = simdFmaddD(VD5, z4, VD3);
    polyVD0        = simdFmaddD(VD4, z4, VD2);
    polyVD1        = simdFmaddD(polyVD1, z4, VD1);
    polyVD0        = simdFmaddD(polyVD0, z4, VD0);
    polyVD0        = simdFmaddD(polyVD1, z2, polyVD0);

    polyVD0        = simdInvD(polyVD0);

    polyVN1        = simdFmaddD(VN9, z4, VN7);
    polyVN0        = simdFmaddD(VN8, z4, VN6);
    polyVN1        = simdFmaddD(polyVN1, z4, VN5);
    polyVN0        = simdFmaddD(polyVN0, z4, VN4);
    polyVN1        = simdFmaddD(polyVN1, z4, VN3);
    polyVN0        = simdFmaddD(polyVN0, z4, VN2);
    polyVN1        = simdFmaddD(polyVN1, z4, VN1);
    polyVN0        = simdFmaddD(polyVN0, z4, VN0);
    polyVN0        = simdFmaddD(polyVN1, z2, polyVN0);

    return simdMulD(polyVN0, polyVD0);
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
 * \ref simdInvsqrtSingleAccuracy.
 *
 *  \param x Argument that must be >0. This routine does not check arguments.
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline SimdDouble gmx_simdcall
simdInvsqrtSingleAccuracyD(SimdDouble x)
{
    SimdDouble lu = simdRsqrtD(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterD(lu, x);
#endif
    return lu;
}

/*! \brief 1/sqrt(x) for masked-in entries of SIMD double, but in single accuracy.
 *
 * \copydetails simdInvsqrtMaskF
 */
static inline SimdDouble
simdInvsqrtMaskSingleAccuracyD(SimdDouble x, SimdDBool m)
{
    SimdDouble lu = simdRsqrtMaskD(x, m);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRsqrtIterD(lu, x);
#endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for two SIMD doubles, but single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdInvsqrtPairSingleAccuracy.
 *
 * \param x0  First set of arguments, x0 must be positive - no argument checking.
 * \param x1  Second set of arguments, x1 must be positive - no argument checking.
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 *  In particular for double precision we can sometimes calculate square root
 *  pairs slightly faster by using single precision until the very last step.
 */
static inline void gmx_simdcall
simdInvsqrtPairSingleAccuracyD(SimdDouble x0,    SimdDouble x1,
                               SimdDouble *out0, SimdDouble *out1)
{
#if GMX_SIMD_HAVE_FLOAT && (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH) && (GMX_SIMD_RSQRT_BITS < 22)
    SimdFloat  xf  = simdCvtDD2F(x0, x1);
    SimdFloat  luf = simdRsqrtF(xf);
    SimdDouble lu0, lu1;
    /* Intermediate target is single - mantissa+1 bits */
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = simdRsqrtIterF(luf, xf);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = simdRsqrtIterF(luf, xf);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = simdRsqrtIterF(luf, xf);
#endif
    simdCvtF2DD(luf, &lu0, &lu1);
    /* We now have single-precision accuracy values in lu0/lu1 */
    *out0 = lu0;
    *out1 = lu1;
#else
    *out0 = simdInvsqrtSingleAccuracyD(x0);
    *out1 = simdInvsqrtSingleAccuracyD(x1);
#endif
}


/*! \brief Calculate 1/x for SIMD double, but in single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdInvSingleAccuracy.
 *
 *  \param x Argument that must be nonzero. This routine does not check arguments.
 *  \return 1/x. Result is undefined if your argument was invalid.
 */
static inline SimdDouble gmx_simdcall
simdInvSingleAccuracyD(SimdDouble x)
{
    SimdDouble lu = simdRcpD(x);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterD(lu, x);
#endif
    return lu;
}

/*! \brief 1/x for masked entries of SIMD double, single accuracy.
 *
 * \copydetails simdInvMaskF
 */
static inline SimdDouble
simdInvMaskSingleAccuracyD(SimdDouble x, SimdDBool m)
{
    SimdDouble lu = simdRcpMaskD(x, m);
#if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterD(lu, x);
#endif
#if (GMX_SIMD_RCP_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simdRcpIterD(lu, x);
#endif
    return lu;
}


/*! \brief Calculate sqrt(x) (correct for 0.0) for SIMD double, single accuracy.
 *
 * You should normally call the real-precision routine \ref simdSqrt.
 *
 *  \param x Argument that must be >=0.
 *  \return sqrt(x). If x=0, the result will correctly be set to 0.
 *          The result is undefined if the input value is negative.
 */
static inline SimdDouble gmx_simdcall
simdSqrtSingleAccuracyD(SimdDouble x)
{
    SimdDouble  res;

    res  = simdInvsqrtMaskSingleAccuracyD(x, simdCmpNzD(x));
    return simdMulD(res, x);
}

/*! \brief SIMD log(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdLogSingleAccuracy.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 */
#ifndef simdLogSingleAccuracyD
static inline SimdDouble gmx_simdcall
simdLogSingleAccuracyD(SimdDouble x)
{
    const SimdDouble  half       = simdSet1D(0.5);
    const SimdDouble  one        = simdSet1D(1.0);
    const SimdDouble  sqrt2      = simdSet1D(sqrt(2.0));
    const SimdDouble  corr       = simdSet1D(0.693147180559945286226764);
    const SimdDouble  CL9        = simdSet1D(0.2371599674224853515625);
    const SimdDouble  CL7        = simdSet1D(0.285279005765914916992188);
    const SimdDouble  CL5        = simdSet1D(0.400005519390106201171875);
    const SimdDouble  CL3        = simdSet1D(0.666666567325592041015625);
    const SimdDouble  CL1        = simdSet1D(2.0);
    SimdDouble        fexp, x2, p;
    SimdDBool         mask;

    fexp  = simdGetExponentD(x);
    x     = simdGetMantissaD(x);

    mask  = simdCmpLtD(sqrt2, x);
    /* Adjust to non-IEEE format for x>sqrt(2): exponent += 1, mantissa *= 0.5 */
    fexp  = simdAddD(fexp, simdMaskD(one, mask));
    x     = simdMulD(x, simdBlendD(one, half, mask));

    x     = simdMulD( simdSubD(x, one), simdInvSingleAccuracyD( simdAddD(x, one) ) );
    x2    = simdMulD(x, x);

    p     = simdFmaddD(CL9, x2, CL7);
    p     = simdFmaddD(p, x2, CL5);
    p     = simdFmaddD(p, x2, CL3);
    p     = simdFmaddD(p, x2, CL1);
    p     = simdFmaddD(p, x, simdMulD(corr, fexp));

    return p;
}
#endif

#ifndef simdExp2SingleAccuracyD
/*! \brief SIMD 2^x. Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdExp2SingleAccuracy.
 *
 * \param x Argument.
 * \result 2^x. Undefined if input argument caused overflow.
 */
static inline SimdDouble gmx_simdcall
simdExp2SingleAccuracyD(SimdDouble x)
{
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const SimdDouble  arglimit = simdSet1D(126.0);
    const SimdDouble  CC6      = simdSet1D(0.0001534581200287996416911311);
    const SimdDouble  CC5      = simdSet1D(0.001339993121934088894618990);
    const SimdDouble  CC4      = simdSet1D(0.009618488957115180159497841);
    const SimdDouble  CC3      = simdSet1D(0.05550328776964726865751735);
    const SimdDouble  CC2      = simdSet1D(0.2402264689063408646490722);
    const SimdDouble  CC1      = simdSet1D(0.6931472057372680777553816);
    const SimdDouble  one      = simdSet1D(1.0);

    SimdDouble        fexppart;
    SimdDouble        intpart;
    SimdDouble        p;
    SimdDBool         valuemask;

    fexppart  = simdSetExponentD(x);
    intpart   = simdRoundD(x);
    valuemask = simdCmpLeD(simdAbsD(x), arglimit);
    fexppart  = simdMaskD(fexppart, valuemask);
    x         = simdSubD(x, intpart);

    p         = simdFmaddD(CC6, x, CC5);
    p         = simdFmaddD(p, x, CC4);
    p         = simdFmaddD(p, x, CC3);
    p         = simdFmaddD(p, x, CC2);
    p         = simdFmaddD(p, x, CC1);
    p         = simdFmaddD(p, x, one);
    x         = simdMulD(p, fexppart);
    return x;
}
#endif

#ifndef simdExpSingleAccuracyD
/*! \brief SIMD exp(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdExpSingleAccuracy.
 *
 * \param x Argument.
 * \result exp(x). Undefined if input argument caused overflow.
 */
static inline SimdDouble gmx_simdcall
simdExpSingleAccuracyD(SimdDouble x)
{
    const SimdDouble  argscale     = simdSet1D(1.44269504088896341);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const SimdDouble  arglimit     = simdSet1D(126.0);
    const SimdDouble  invargscale  = simdSet1D(0.69314718055994528623);
    const SimdDouble  CC4          = simdSet1D(0.00136324646882712841033936);
    const SimdDouble  CC3          = simdSet1D(0.00836596917361021041870117);
    const SimdDouble  CC2          = simdSet1D(0.0416710823774337768554688);
    const SimdDouble  CC1          = simdSet1D(0.166665524244308471679688);
    const SimdDouble  CC0          = simdSet1D(0.499999850988388061523438);
    const SimdDouble  one          = simdSet1D(1.0);
    SimdDouble        fexppart;
    SimdDouble        intpart;
    SimdDouble        y, p;
    SimdDBool         valuemask;

    y         = simdMulD(x, argscale);
    fexppart  = simdSetExponentD(y);  /* rounds to nearest int internally */
    intpart   = simdRoundD(y);        /* use same rounding algorithm here */
    valuemask = simdCmpLeD(simdAbsD(y), arglimit);
    fexppart  = simdMaskD(fexppart, valuemask);

    /* Extended precision arithmetics not needed since
     * we have double precision and only need single accuracy.
     */
    x         = simdFnmaddD(invargscale, intpart, x);

    p         = simdFmaddD(CC4, x, CC3);
    p         = simdFmaddD(p, x, CC2);
    p         = simdFmaddD(p, x, CC1);
    p         = simdFmaddD(p, x, CC0);
    p         = simdFmaddD(simdMulD(x, x), p, x);
    p         = simdAddD(p, one);
    x         = simdMulD(p, fexppart);
    return x;
}
#endif

/*! \brief SIMD erf(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdErfSingleAccuracy.
 *
 * \param x The value to calculate erf(x) for.
 * \result erf(x)
 *
 * This routine achieves very close to single precision, but we do not care about
 * the last bit or the subnormal result range.
 */
static inline SimdDouble gmx_simdcall
simdErfSingleAccuracyD(SimdDouble x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const SimdDouble  CA6      = simdSet1D(7.853861353153693e-5);
    const SimdDouble  CA5      = simdSet1D(-8.010193625184903e-4);
    const SimdDouble  CA4      = simdSet1D(5.188327685732524e-3);
    const SimdDouble  CA3      = simdSet1D(-2.685381193529856e-2);
    const SimdDouble  CA2      = simdSet1D(1.128358514861418e-1);
    const SimdDouble  CA1      = simdSet1D(-3.761262582423300e-1);
    const SimdDouble  CA0      = simdSet1D(1.128379165726710);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const SimdDouble  CB9      = simdSet1D(-0.0018629930017603923);
    const SimdDouble  CB8      = simdSet1D(0.003909821287598495);
    const SimdDouble  CB7      = simdSet1D(-0.0052094582210355615);
    const SimdDouble  CB6      = simdSet1D(0.005685614362160572);
    const SimdDouble  CB5      = simdSet1D(-0.0025367682853477272);
    const SimdDouble  CB4      = simdSet1D(-0.010199799682318782);
    const SimdDouble  CB3      = simdSet1D(0.04369575504816542);
    const SimdDouble  CB2      = simdSet1D(-0.11884063474674492);
    const SimdDouble  CB1      = simdSet1D(0.2732120154030589);
    const SimdDouble  CB0      = simdSet1D(0.42758357702025784);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const SimdDouble  CC10     = simdSet1D(-0.0445555913112064);
    const SimdDouble  CC9      = simdSet1D(0.21376355144663348);
    const SimdDouble  CC8      = simdSet1D(-0.3473187200259257);
    const SimdDouble  CC7      = simdSet1D(0.016690861551248114);
    const SimdDouble  CC6      = simdSet1D(0.7560973182491192);
    const SimdDouble  CC5      = simdSet1D(-1.2137903600145787);
    const SimdDouble  CC4      = simdSet1D(0.8411872321232948);
    const SimdDouble  CC3      = simdSet1D(-0.08670413896296343);
    const SimdDouble  CC2      = simdSet1D(-0.27124782687240334);
    const SimdDouble  CC1      = simdSet1D(-0.0007502488047806069);
    const SimdDouble  CC0      = simdSet1D(0.5642114853803148);
    const SimdDouble  one      = simdSet1D(1.0);
    const SimdDouble  two      = simdSet1D(2.0);

    SimdDouble        x2, x4, y;
    SimdDouble        t, t2, w, w2;
    SimdDouble        pA0, pA1, pB0, pB1, pC0, pC1;
    SimdDouble        expmx2;
    SimdDouble        res_erf, res_erfc, res;
    SimdDBool         mask, msk_erf;

    /* Calculate erf() */
    x2   = simdMulD(x, x);
    x4   = simdMulD(x2, x2);

    pA0  = simdFmaddD(CA6, x4, CA4);
    pA1  = simdFmaddD(CA5, x4, CA3);
    pA0  = simdFmaddD(pA0, x4, CA2);
    pA1  = simdFmaddD(pA1, x4, CA1);
    pA0  = simdMulD(pA0, x4);
    pA0  = simdFmaddD(pA1, x2, pA0);
    /* Constant term must come last for precision reasons */
    pA0  = simdAddD(pA0, CA0);

    res_erf = simdMulD(x, pA0);

    /* Calculate erfc */
    y       = simdAbsD(x);
    msk_erf = simdCmpLeD(simdSet1D(0.75), y);
    t       = simdInvMaskSingleAccuracyD(y, msk_erf);
    w       = simdSubD(t, one);
    t2      = simdMulD(t, t);
    w2      = simdMulD(w, w);

    expmx2  = simdExpSingleAccuracyD( simdNegD( simdMulD(y, y)));

    pB1  = simdFmaddD(CB9, w2, CB7);
    pB0  = simdFmaddD(CB8, w2, CB6);
    pB1  = simdFmaddD(pB1, w2, CB5);
    pB0  = simdFmaddD(pB0, w2, CB4);
    pB1  = simdFmaddD(pB1, w2, CB3);
    pB0  = simdFmaddD(pB0, w2, CB2);
    pB1  = simdFmaddD(pB1, w2, CB1);
    pB0  = simdFmaddD(pB0, w2, CB0);
    pB0  = simdFmaddD(pB1, w, pB0);

    pC0  = simdFmaddD(CC10, t2, CC8);
    pC1  = simdFmaddD(CC9, t2, CC7);
    pC0  = simdFmaddD(pC0, t2, CC6);
    pC1  = simdFmaddD(pC1, t2, CC5);
    pC0  = simdFmaddD(pC0, t2, CC4);
    pC1  = simdFmaddD(pC1, t2, CC3);
    pC0  = simdFmaddD(pC0, t2, CC2);
    pC1  = simdFmaddD(pC1, t2, CC1);

    pC0  = simdFmaddD(pC0, t2, CC0);
    pC0  = simdFmaddD(pC1, t, pC0);
    pC0  = simdMulD(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = simdCmpLtD(two, y);
    res_erfc = simdBlendD(pB0, pC0, mask);
    res_erfc = simdMulD(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = simdCmpLtD(x, simdSetZeroD());
    res_erfc = simdBlendD(res_erfc, simdSubD(two, res_erfc), mask);

    /* Select erf() or erfc() */
    mask = simdCmpLtD(y, simdSet1D(0.75));
    res  = simdBlendD(simdSubD(one, res_erfc), res_erf, mask);

    return res;
}

/*! \brief SIMD erfc(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdErfcSingleAccuracy.
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
static inline SimdDouble gmx_simdcall
simdErfcSingleAccuracyD(SimdDouble x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const SimdDouble  CA6      = simdSet1D(7.853861353153693e-5);
    const SimdDouble  CA5      = simdSet1D(-8.010193625184903e-4);
    const SimdDouble  CA4      = simdSet1D(5.188327685732524e-3);
    const SimdDouble  CA3      = simdSet1D(-2.685381193529856e-2);
    const SimdDouble  CA2      = simdSet1D(1.128358514861418e-1);
    const SimdDouble  CA1      = simdSet1D(-3.761262582423300e-1);
    const SimdDouble  CA0      = simdSet1D(1.128379165726710);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const SimdDouble  CB9      = simdSet1D(-0.0018629930017603923);
    const SimdDouble  CB8      = simdSet1D(0.003909821287598495);
    const SimdDouble  CB7      = simdSet1D(-0.0052094582210355615);
    const SimdDouble  CB6      = simdSet1D(0.005685614362160572);
    const SimdDouble  CB5      = simdSet1D(-0.0025367682853477272);
    const SimdDouble  CB4      = simdSet1D(-0.010199799682318782);
    const SimdDouble  CB3      = simdSet1D(0.04369575504816542);
    const SimdDouble  CB2      = simdSet1D(-0.11884063474674492);
    const SimdDouble  CB1      = simdSet1D(0.2732120154030589);
    const SimdDouble  CB0      = simdSet1D(0.42758357702025784);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const SimdDouble  CC10     = simdSet1D(-0.0445555913112064);
    const SimdDouble  CC9      = simdSet1D(0.21376355144663348);
    const SimdDouble  CC8      = simdSet1D(-0.3473187200259257);
    const SimdDouble  CC7      = simdSet1D(0.016690861551248114);
    const SimdDouble  CC6      = simdSet1D(0.7560973182491192);
    const SimdDouble  CC5      = simdSet1D(-1.2137903600145787);
    const SimdDouble  CC4      = simdSet1D(0.8411872321232948);
    const SimdDouble  CC3      = simdSet1D(-0.08670413896296343);
    const SimdDouble  CC2      = simdSet1D(-0.27124782687240334);
    const SimdDouble  CC1      = simdSet1D(-0.0007502488047806069);
    const SimdDouble  CC0      = simdSet1D(0.5642114853803148);
    const SimdDouble  one      = simdSet1D(1.0);
    const SimdDouble  two      = simdSet1D(2.0);

    SimdDouble        x2, x4, y;
    SimdDouble        t, t2, w, w2;
    SimdDouble        pA0, pA1, pB0, pB1, pC0, pC1;
    SimdDouble        expmx2;
    SimdDouble        res_erf, res_erfc, res;
    SimdDBool         mask, msk_erf;

    /* Calculate erf() */
    x2     = simdMulD(x, x);
    x4     = simdMulD(x2, x2);

    pA0  = simdFmaddD(CA6, x4, CA4);
    pA1  = simdFmaddD(CA5, x4, CA3);
    pA0  = simdFmaddD(pA0, x4, CA2);
    pA1  = simdFmaddD(pA1, x4, CA1);
    pA1  = simdMulD(pA1, x2);
    pA0  = simdFmaddD(pA0, x4, pA1);
    /* Constant term must come last for precision reasons */
    pA0  = simdAddD(pA0, CA0);

    res_erf = simdMulD(x, pA0);

    /* Calculate erfc */
    y       = simdAbsD(x);
    msk_erf = simdCmpLeD(simdSet1D(0.75), y);
    t       = simdInvMaskSingleAccuracyD(y, msk_erf);
    w       = simdSubD(t, one);
    t2      = simdMulD(t, t);
    w2      = simdMulD(w, w);

    expmx2  = simdExpSingleAccuracyD( simdNegD( simdMulD(y, y) ) );

    pB1  = simdFmaddD(CB9, w2, CB7);
    pB0  = simdFmaddD(CB8, w2, CB6);
    pB1  = simdFmaddD(pB1, w2, CB5);
    pB0  = simdFmaddD(pB0, w2, CB4);
    pB1  = simdFmaddD(pB1, w2, CB3);
    pB0  = simdFmaddD(pB0, w2, CB2);
    pB1  = simdFmaddD(pB1, w2, CB1);
    pB0  = simdFmaddD(pB0, w2, CB0);
    pB0  = simdFmaddD(pB1, w, pB0);

    pC0  = simdFmaddD(CC10, t2, CC8);
    pC1  = simdFmaddD(CC9, t2, CC7);
    pC0  = simdFmaddD(pC0, t2, CC6);
    pC1  = simdFmaddD(pC1, t2, CC5);
    pC0  = simdFmaddD(pC0, t2, CC4);
    pC1  = simdFmaddD(pC1, t2, CC3);
    pC0  = simdFmaddD(pC0, t2, CC2);
    pC1  = simdFmaddD(pC1, t2, CC1);

    pC0  = simdFmaddD(pC0, t2, CC0);
    pC0  = simdFmaddD(pC1, t, pC0);
    pC0  = simdMulD(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = simdCmpLtD(two, y);
    res_erfc = simdBlendD(pB0, pC0, mask);
    res_erfc = simdMulD(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = simdCmpLtD(x, simdSetZeroD());
    res_erfc = simdBlendD(res_erfc, simdSubD(two, res_erfc), mask);

    /* Select erf() or erfc() */
    mask = simdCmpLtD(y, simdSet1D(0.75));
    res  = simdBlendD(res_erfc, simdSubD(one, res_erf), mask);

    return res;
}

/*! \brief SIMD sin \& cos. Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdSinCosSingleAccuracy.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 */
static inline void gmx_simdcall
simdSinCosSingleAccuracyD(SimdDouble x, SimdDouble *sinval, SimdDouble *cosval)
{
    /* Constants to subtract Pi/4*x from y while minimizing precision loss */
    const SimdDouble  argred0         = simdSet1D(2*0.78539816290140151978);
    const SimdDouble  argred1         = simdSet1D(2*4.9604678871439933374e-10);
    const SimdDouble  argred2         = simdSet1D(2*1.1258708853173288931e-18);
    const SimdDouble  two_over_pi     = simdSet1D(2.0/M_PI);
    const SimdDouble  const_sin2      = simdSet1D(-1.9515295891e-4);
    const SimdDouble  const_sin1      = simdSet1D( 8.3321608736e-3);
    const SimdDouble  const_sin0      = simdSet1D(-1.6666654611e-1);
    const SimdDouble  const_cos2      = simdSet1D( 2.443315711809948e-5);
    const SimdDouble  const_cos1      = simdSet1D(-1.388731625493765e-3);
    const SimdDouble  const_cos0      = simdSet1D( 4.166664568298827e-2);

    const SimdDouble  half            = simdSet1D(0.5);
    const SimdDouble  one             = simdSet1D(1.0);
    SimdDouble        ssign, csign;
    SimdDouble        x2, y, z, psin, pcos, sss, ccc;
    SimdDBool         mask;

#if GMX_SIMD_HAVE_FINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    const SimdDInt32  ione            = simdSet1DI(1);
    const SimdDInt32  itwo            = simdSet1DI(2);
    SimdDInt32        iy;

    z       = simdMulD(x, two_over_pi);
    iy      = simdCvtD2I(z);
    y       = simdRoundD(z);

    mask    = simdCvtDIB2DB(simdCmpEqDI(simdAndDI(iy, ione), simdSetZeroDI()));
    ssign   = simdMaskD(simdSet1D(-0.0), simdCvtDIB2DB(simdCmpEqDI(simdAndDI(iy, itwo), itwo)));
    csign   = simdMaskD(simdSet1D(-0.0), simdCvtDIB2DB(simdCmpEqDI(simdAndDI(simdAddDI(iy, ione), itwo), itwo)));
#else
    const SimdDouble  quarter         = simdSet1D(0.25);
    const SimdDouble  minusquarter    = simdSet1D(-0.25);
    SimdDouble        q;
    SimdDBool         m1, m2, m3;

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
    x       = simdAbsD(x);
    /* It is critical that half-way cases are rounded down */
    z       = simdFmaddD(x, two_over_pi, half);
    y       = simdTruncD(z);
    q       = simdMulD(z, quarter);
    q       = simdSubD(q, simdTruncD(q));
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q       = simdSubD(q, half);
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#    if GMX_SIMD_HAVE_LOGICAL
    ssign   = simdAndD(ssign, simdSet1D(-0.0));
    csign   = simdAndNotD(q, simdSet1D(-0.0));
    ssign   = simdXorD(ssign, csign);
#    else
    csign   = simdXorSignD(simdSet1D(-1.0), q);

    ssign   = simdXorSignD(ssign, csign);    /* swap ssign if csign was set. */
#    endif
    /* Check if y had value 1 or 3 (remember we subtracted 0.5 from q) */
    m1      = simdCmpLtD(q, minusquarter);
    m2      = simdCmpLeD(simdSetZeroD(), q);
    m3      = simdCmpLtD(q, quarter);
    m2      = simdAndDB(m2, m3);
    mask    = simdOrDB(m1, m2);
    /* where mask is FALSE, set sign. */
    csign   = simdXorSignD(csign, simdBlendD(simdSet1D(-1.0), one, mask));
#endif
    x       = simdFnmaddD(y, argred0, x);
    x       = simdFnmaddD(y, argred1, x);
    x       = simdFnmaddD(y, argred2, x);
    x2      = simdMulD(x, x);

    psin    = simdFmaddD(const_sin2, x2, const_sin1);
    psin    = simdFmaddD(psin, x2, const_sin0);
    psin    = simdFmaddD(psin, simdMulD(x, x2), x);
    pcos    = simdFmaddD(const_cos2, x2, const_cos1);
    pcos    = simdFmaddD(pcos, x2, const_cos0);
    pcos    = simdFmsubD(pcos, x2, half);
    pcos    = simdFmaddD(pcos, x2, one);

    sss     = simdBlendD(pcos, psin, mask);
    ccc     = simdBlendD(psin, pcos, mask);
    /* See comment for GMX_SIMD_HAVE_LOGICAL section above. */
#if GMX_SIMD_HAVE_LOGICAL
    *sinval = simdXorD(sss, ssign);
    *cosval = simdXorD(ccc, csign);
#else
    *sinval = simdXorSignD(sss, ssign);
    *cosval = simdXorSignD(ccc, csign);
#endif
}

/*! \brief SIMD sin(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdSinSingleAccuracy.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref simdSinCos and waste a factor 2 in performance.
 */
static inline SimdDouble gmx_simdcall
simdSinSingleAccuracyD(SimdDouble x)
{
    SimdDouble s, c;
    simdSinCosSingleAccuracyD(x, &s, &c);
    return s;
}

/*! \brief SIMD cos(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdCosSingleAccuracy.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref simdSinCos and waste a factor 2 in performance.
 */
static inline SimdDouble gmx_simdcall
simdCosSingleAccuracyD(SimdDouble x)
{
    SimdDouble s, c;
    simdSinCosSingleAccuracyD(x, &s, &c);
    return c;
}

/*! \brief SIMD tan(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdTanSingleAccuracy.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 */
static inline SimdDouble gmx_simdcall
simdTanSingleAccuracyD(SimdDouble x)
{
    const SimdDouble  argred0         = simdSet1D(2*0.78539816290140151978);
    const SimdDouble  argred1         = simdSet1D(2*4.9604678871439933374e-10);
    const SimdDouble  argred2         = simdSet1D(2*1.1258708853173288931e-18);
    const SimdDouble  two_over_pi     = simdSet1D(2.0/M_PI);
    const SimdDouble  CT6             = simdSet1D(0.009498288995810566122993911);
    const SimdDouble  CT5             = simdSet1D(0.002895755790837379295226923);
    const SimdDouble  CT4             = simdSet1D(0.02460087336161924491836265);
    const SimdDouble  CT3             = simdSet1D(0.05334912882656359828045988);
    const SimdDouble  CT2             = simdSet1D(0.1333989091464957704418495);
    const SimdDouble  CT1             = simdSet1D(0.3333307599244198227797507);

    SimdDouble        x2, p, y, z;
    SimdDBool         mask;

#if GMX_SIMD_HAVE_FINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    SimdDInt32  iy;
    SimdDInt32  ione = simdSet1DI(1);

    z       = simdMulD(x, two_over_pi);
    iy      = simdCvtD2I(z);
    y       = simdRoundD(z);
    mask    = simdCvtDIB2DB(simdCmpEqDI(simdAndDI(iy, ione), ione));

    x       = simdFnmaddD(y, argred0, x);
    x       = simdFnmaddD(y, argred1, x);
    x       = simdFnmaddD(y, argred2, x);
    x       = simdXorD(simdMaskD(simdSet1D(-0.0), mask), x);
#else
    const SimdDouble  quarter         = simdSet1D(0.25);
    const SimdDouble  half            = simdSet1D(0.5);
    const SimdDouble  threequarter    = simdSet1D(0.75);
    SimdDouble        w, q;
    SimdDBool         m1, m2, m3;

    w       = simdAbsD(x);
    z       = simdFmaddD(w, two_over_pi, half);
    y       = simdTruncD(z);
    q       = simdMulD(z, quarter);
    q       = simdSubD(q, simdTruncD(q));
    m1      = simdCmpLeD(quarter, q);
    m2      = simdCmpLtD(q, half);
    m3      = simdCmpLeD(threequarter, q);
    m1      = simdAndDB(m1, m2);
    mask    = simdOrDB(m1, m3);
    w       = simdFnmaddD(y, argred0, w);
    w       = simdFnmaddD(y, argred1, w);
    w       = simdFnmaddD(y, argred2, w);

    w       = simdBlendD(w, simdNegD(w), mask);
    x       = simdXorSignD(w, x);
#endif
    x2      = simdMulD(x, x);
    p       = simdFmaddD(CT6, x2, CT5);
    p       = simdFmaddD(p, x2, CT4);
    p       = simdFmaddD(p, x2, CT3);
    p       = simdFmaddD(p, x2, CT2);
    p       = simdFmaddD(p, x2, CT1);
    p       = simdFmaddD(x2, simdMulD(p, x), x);

    p       = simdBlendD( p, simdInvMaskSingleAccuracyD(p, mask), mask);
    return p;
}

/*! \brief SIMD asin(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdAsinSingleAccuracy.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 */
static inline SimdDouble gmx_simdcall
simdAsinSingleAccuracyD(SimdDouble x)
{
    const SimdDouble limitlow   = simdSet1D(1e-4);
    const SimdDouble half       = simdSet1D(0.5);
    const SimdDouble one        = simdSet1D(1.0);
    const SimdDouble halfpi     = simdSet1D(M_PI/2.0);
    const SimdDouble CC5        = simdSet1D(4.2163199048E-2);
    const SimdDouble CC4        = simdSet1D(2.4181311049E-2);
    const SimdDouble CC3        = simdSet1D(4.5470025998E-2);
    const SimdDouble CC2        = simdSet1D(7.4953002686E-2);
    const SimdDouble CC1        = simdSet1D(1.6666752422E-1);
    SimdDouble       xabs;
    SimdDouble       z, z1, z2, q, q1, q2;
    SimdDouble       pA, pB;
    SimdDBool        mask, mask2;

    xabs  = simdAbsD(x);
    mask  = simdCmpLtD(half, xabs);
    z1    = simdMulD(half, simdSubD(one, xabs));
    mask2 = simdCmpLtD(xabs, one);
    q1    = simdMulD(z1, simdInvsqrtMaskSingleAccuracyD(z1, mask2));
    q2    = xabs;
    z2    = simdMulD(q2, q2);
    z     = simdBlendD(z2, z1, mask);
    q     = simdBlendD(q2, q1, mask);

    z2    = simdMulD(z, z);
    pA    = simdFmaddD(CC5, z2, CC3);
    pB    = simdFmaddD(CC4, z2, CC2);
    pA    = simdFmaddD(pA, z2, CC1);
    pA    = simdMulD(pA, z);
    z     = simdFmaddD(pB, z2, pA);
    z     = simdFmaddD(z, q, q);
    q2    = simdSubD(halfpi, z);
    q2    = simdSubD(q2, z);
    z     = simdBlendD(z, q2, mask);

    mask  = simdCmpLtD(limitlow, xabs);
    z     = simdBlendD( xabs, z, mask );
    z     = simdXorSignD(z, x);

    return z;
}

/*! \brief SIMD acos(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdAcosSingleAccuracy.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 */
static inline SimdDouble gmx_simdcall
simdAcosSingleAccuracyD(SimdDouble x)
{
    const SimdDouble one       = simdSet1D(1.0);
    const SimdDouble half      = simdSet1D(0.5);
    const SimdDouble pi        = simdSet1D(M_PI);
    const SimdDouble halfpi    = simdSet1D(M_PI/2.0);
    SimdDouble       xabs;
    SimdDouble       z, z1, z2, z3;
    SimdDBool        mask1, mask2, mask3;

    xabs  = simdAbsD(x);
    mask1 = simdCmpLtD(half, xabs);
    mask2 = simdCmpLtD(simdSetZeroD(), x);

    z     = simdMulD(half, simdSubD(one, xabs));
    mask3 = simdCmpLtD(xabs, one);
    z     = simdMulD(z, simdInvsqrtMaskSingleAccuracyD(z, mask3));
    z     = simdBlendD(x, z, mask1);
    z     = simdAsinSingleAccuracyD(z);

    z2    = simdAddD(z, z);
    z1    = simdSubD(pi, z2);
    z3    = simdSubD(halfpi, z);
    z     = simdBlendD(z1, z2, mask2);
    z     = simdBlendD(z3, z, mask1);

    return z;
}

/*! \brief SIMD asin(x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdAtanSingleAccuracy.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x), same argument/value range as standard math library.
 */
static inline SimdDouble gmx_simdcall
simdAtanSingleAccuracyD(SimdDouble x)
{
    const SimdDouble halfpi    = simdSet1D(M_PI/2);
    const SimdDouble CA17      = simdSet1D(0.002823638962581753730774);
    const SimdDouble CA15      = simdSet1D(-0.01595690287649631500244);
    const SimdDouble CA13      = simdSet1D(0.04250498861074447631836);
    const SimdDouble CA11      = simdSet1D(-0.07489009201526641845703);
    const SimdDouble CA9       = simdSet1D(0.1063479334115982055664);
    const SimdDouble CA7       = simdSet1D(-0.1420273631811141967773);
    const SimdDouble CA5       = simdSet1D(0.1999269574880599975585);
    const SimdDouble CA3       = simdSet1D(-0.3333310186862945556640);
    SimdDouble       x2, x3, x4, pA, pB;
    SimdDBool        mask, mask2;

    mask  = simdCmpLtD(x, simdSetZeroD());
    x     = simdAbsD(x);
    mask2 = simdCmpLtD(simdSet1D(1.0), x);
    x     = simdBlendD(x, simdInvMaskSingleAccuracyD(x, mask2), mask2);

    x2    = simdMulD(x, x);
    x3    = simdMulD(x2, x);
    x4    = simdMulD(x2, x2);
    pA    = simdFmaddD(CA17, x4, CA13);
    pB    = simdFmaddD(CA15, x4, CA11);
    pA    = simdFmaddD(pA, x4, CA9);
    pB    = simdFmaddD(pB, x4, CA7);
    pA    = simdFmaddD(pA, x4, CA5);
    pB    = simdFmaddD(pB, x4, CA3);
    pA    = simdFmaddD(pA, x2, pB);
    pA    = simdFmaddD(pA, x3, x);

    pA    = simdBlendD(pA, simdSubD(halfpi, pA), mask2);
    pA    = simdBlendD(pA, simdNegD(pA), mask);

    return pA;
}

/*! \brief SIMD atan2(y,x). Double precision SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdAtan2SingleAccuracy.
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
static inline SimdDouble gmx_simdcall
simdAtan2SingleAccuracyD(SimdDouble y, SimdDouble x)
{
    const SimdDouble pi          = simdSet1D(M_PI);
    const SimdDouble halfpi      = simdSet1D(M_PI/2.0);
    SimdDouble       xinv, p, aoffset;
    SimdDBool        mask_xnz, mask_ynz, mask_xlt0, mask_ylt0;

    mask_xnz  = simdCmpNzD(x);
    mask_ynz  = simdCmpNzD(y);
    mask_xlt0 = simdCmpLtD(x, simdSetZeroD());
    mask_ylt0 = simdCmpLtD(y, simdSetZeroD());

    aoffset   = simdMaskNotD(halfpi, mask_xnz);
    aoffset   = simdMaskD(aoffset, mask_ynz);

    aoffset   = simdBlendD(aoffset, pi, mask_xlt0);
    aoffset   = simdBlendD(aoffset, simdNegD(aoffset), mask_ylt0);

    xinv      = simdInvMaskSingleAccuracyD(x, mask_xnz);
    p         = simdMulD(y, xinv);
    p         = simdAtanSingleAccuracyD(p);
    p         = simdAddD(p, aoffset);

    return p;
}

/*! \brief Analytical PME force correction, double SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdPmeCorrForceSingleAccuracy.
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
static SimdDouble gmx_simdcall
simdPmeCorrForceSingleAccuracyD(SimdDouble z2)
{
    const SimdDouble  FN6      = simdSet1D(-1.7357322914161492954e-8);
    const SimdDouble  FN5      = simdSet1D(1.4703624142580877519e-6);
    const SimdDouble  FN4      = simdSet1D(-0.000053401640219807709149);
    const SimdDouble  FN3      = simdSet1D(0.0010054721316683106153);
    const SimdDouble  FN2      = simdSet1D(-0.019278317264888380590);
    const SimdDouble  FN1      = simdSet1D(0.069670166153766424023);
    const SimdDouble  FN0      = simdSet1D(-0.75225204789749321333);

    const SimdDouble  FD4      = simdSet1D(0.0011193462567257629232);
    const SimdDouble  FD3      = simdSet1D(0.014866955030185295499);
    const SimdDouble  FD2      = simdSet1D(0.11583842382862377919);
    const SimdDouble  FD1      = simdSet1D(0.50736591960530292870);
    const SimdDouble  FD0      = simdSet1D(1.0);

    SimdDouble        z4;
    SimdDouble        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = simdMulD(z2, z2);

    polyFD0        = simdFmaddD(FD4, z4, FD2);
    polyFD1        = simdFmaddD(FD3, z4, FD1);
    polyFD0        = simdFmaddD(polyFD0, z4, FD0);
    polyFD0        = simdFmaddD(polyFD1, z2, polyFD0);

    polyFD0        = simdInvSingleAccuracyD(polyFD0);

    polyFN0        = simdFmaddD(FN6, z4, FN4);
    polyFN1        = simdFmaddD(FN5, z4, FN3);
    polyFN0        = simdFmaddD(polyFN0, z4, FN2);
    polyFN1        = simdFmaddD(polyFN1, z4, FN1);
    polyFN0        = simdFmaddD(polyFN0, z4, FN0);
    polyFN0        = simdFmaddD(polyFN1, z2, polyFN0);

    return simdMulD(polyFN0, polyFD0);
}



/*! \brief Analytical PME potential correction, double SIMD data, single accuracy.
 *
 * You should normally call the real-precision routine
 * \ref simdPmeCorrPotentialSingleAccuracy.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb potential - see below for details.
 *
 * See \ref simdPmeCorrForceF for details about the approximation.
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
static SimdDouble gmx_simdcall
simdPmeCorrPotentialSingleAccuracyD(SimdDouble z2)
{
    const SimdDouble  VN6      = simdSet1D(1.9296833005951166339e-8);
    const SimdDouble  VN5      = simdSet1D(-1.4213390571557850962e-6);
    const SimdDouble  VN4      = simdSet1D(0.000041603292906656984871);
    const SimdDouble  VN3      = simdSet1D(-0.00013134036773265025626);
    const SimdDouble  VN2      = simdSet1D(0.038657983986041781264);
    const SimdDouble  VN1      = simdSet1D(0.11285044772717598220);
    const SimdDouble  VN0      = simdSet1D(1.1283802385263030286);

    const SimdDouble  VD3      = simdSet1D(0.0066752224023576045451);
    const SimdDouble  VD2      = simdSet1D(0.078647795836373922256);
    const SimdDouble  VD1      = simdSet1D(0.43336185284710920150);
    const SimdDouble  VD0      = simdSet1D(1.0);

    SimdDouble        z4;
    SimdDouble        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = simdMulD(z2, z2);

    polyVD1        = simdFmaddD(VD3, z4, VD1);
    polyVD0        = simdFmaddD(VD2, z4, VD0);
    polyVD0        = simdFmaddD(polyVD1, z2, polyVD0);

    polyVD0        = simdInvSingleAccuracyD(polyVD0);

    polyVN0        = simdFmaddD(VN6, z4, VN4);
    polyVN1        = simdFmaddD(VN5, z4, VN3);
    polyVN0        = simdFmaddD(polyVN0, z4, VN2);
    polyVN1        = simdFmaddD(polyVN1, z4, VN1);
    polyVN0        = simdFmaddD(polyVN0, z4, VN0);
    polyVN0        = simdFmaddD(polyVN1, z2, polyVN0);

    return simdMulD(polyVN0, polyVD0);
}

#endif


/*! \name SIMD4 math functions
 *
 * \note Only a subset of the math functions are implemented for SIMD4.
 *  \{
 */


#if GMX_SIMD4_HAVE_FLOAT

/*************************************************************************
 * SINGLE PRECISION SIMD4 MATH FUNCTIONS - JUST A SMALL SUBSET SUPPORTED *
 *************************************************************************/

/*! \brief SIMD4 utility function to sum a+b+c+d for SIMD4 floats.
 *
 * \copydetails simdSum4F
 */
static inline Simd4Float gmx_simdcall
simd4Sum4F(Simd4Float a, Simd4Float b,
           Simd4Float c, Simd4Float d)
{
    return simd4AddF(simd4AddF(a, b), simd4AddF(c, d));
}

/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD4 float.
 *
 * \copydetails simdRsqrtIterF
 */
static inline Simd4Float gmx_simdcall
simd4RsqrtIterF(Simd4Float lu, Simd4Float x)
{
#    if GMX_SIMD_HAVE_FMA
    return simd4FmaddF(simd4FnmaddF(x, simd4MulF(lu, lu), simd4Set1F(1.0f)), simd4MulF(lu, simd4Set1F(0.5f)), lu);
#    else
    return simd4MulF(simd4Set1F(0.5f), simd4MulF(simd4SubF(simd4Set1F(3.0f), simd4MulF(simd4MulF(lu, lu), x)), lu));
#    endif
}

/*! \brief Calculate 1/sqrt(x) for SIMD4 float.
 *
 * \copydetails simdInvsqrtF
 */
static inline Simd4Float gmx_simdcall
simd4InvsqrtF(Simd4Float x)
{
    Simd4Float lu = simd4RsqrtF(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simd4RsqrtIterF(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simd4RsqrtIterF(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simd4RsqrtIterF(lu, x);
#endif
    return lu;
}


#endif /* GMX_SIMD4_HAVE_FLOAT */



#if GMX_SIMD4_HAVE_DOUBLE
/*************************************************************************
 * DOUBLE PRECISION SIMD4 MATH FUNCTIONS - JUST A SMALL SUBSET SUPPORTED *
 *************************************************************************/


/*! \brief SIMD4 utility function to sum a+b+c+d for SIMD4 doubles.
 *
 * \copydetails simdSum4F
 */
static inline Simd4Double gmx_simdcall
simd4Sum4D(Simd4Double a, Simd4Double b,
           Simd4Double c, Simd4Double d)
{
    return simd4AddD(simd4AddD(a, b), simd4AddD(c, d));
}

/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD4 double.
 *
 * \copydetails simdRsqrtIterF
 */
static inline Simd4Double gmx_simdcall
simd4RsqrtIterD(Simd4Double lu, Simd4Double x)
{
#if GMX_SIMD_HAVE_FMA
    return simd4FmaddD(simd4FnmaddD(x, simd4MulD(lu, lu), simd4Set1D(1.0)), simd4MulD(lu, simd4Set1D(0.5)), lu);
#else
    return simd4MulD(simd4Set1D(0.5), simd4MulD(simd4SubD(simd4Set1D(3.0), simd4MulD(simd4MulD(lu, lu), x)), lu));
#endif
}

/*! \brief Calculate 1/sqrt(x) for SIMD4 double.
 *
 * \copydetails simdInvsqrtF
 */
static inline Simd4Double gmx_simdcall
simd4InvsqrtD(Simd4Double x)
{
    Simd4Double lu = simd4RsqrtD(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simd4RsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simd4RsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simd4RsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = simd4RsqrtIterD(lu, x);
#endif
    return lu;
}


/**********************************************************************
 * SIMD4 MATH FUNCTIONS WITH DOUBLE PREC. DATA, SINGLE PREC. ACCURACY *
 **********************************************************************/

/*! \brief Calculate 1/sqrt(x) for SIMD4 double, but in single accuracy.
 *
 * \copydetails simdInvsqrtSingleAccuracyD
 */
static inline Simd4Double gmx_simdcall
simd4InvsqrtSingleAccuracyD(Simd4Double x)
{
    Simd4Double lu = simd4RsqrtD(x);
#if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simd4RsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simd4RsqrtIterD(lu, x);
#endif
#if (GMX_SIMD_RSQRT_BITS*4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = simd4RsqrtIterD(lu, x);
#endif
    return lu;
}



#endif /* GMX_SIMD4_HAVE_DOUBLE */

/*! \} */

/**********************************************************************
 *      WRAPPER FUNCTIONS FOR DEFAULT GROMACS "REAL" PRECISION        *
 **********************************************************************/

/*! \brief SIMD float utility to sum a+b+c+d.
 *
 * \copydetails simdSum4F
 */
static inline SimdReal gmx_simdcall
simdSum4(SimdReal a, SimdReal b,
         SimdReal c, SimdReal d)
{
#ifdef GMX_DOUBLE
    return simdSum4D(a, b, c, d);
#else
    return simdSum4F(a, b, c, d);
#endif
}

/*! \brief Return -a if b is negative, SIMD float.
 *
 * \copydetails simdXorSignF
 */
static inline SimdReal gmx_simdcall
simdXorSign(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdXorSignD(a, b);
#else
    return simdXorSignF(a, b);
#endif
}

#if GMX_SIMD4_HAVE_REAL
/*! \brief SIMD4 float utility to sum a+b+c+d.
 *
 * \copydetails simd4Sum4F
 */
static inline Simd4Real gmx_simdcall
simd4Sum4(Simd4Real a, Simd4Real b,
          Simd4Real c, Simd4Real d)
{
#ifdef GMX_DOUBLE
    return simd4Sum4D(a, b, c, d);
#else
    return simd4Sum4F(a, b, c, d);
#endif
}
#endif  // GMX_SIMD4_HAVE_REAL

/*! \brief Calculate 1/sqrt(x) for SIMD real.
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdReal gmx_simdcall
simdInvsqrt(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdInvsqrtSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdInvsqrtD(x);
#else
    return simdInvsqrtF(x);
#endif
}

/*! \brief Calculate 1/sqrt(x) for two SIMD reals.
 *
 * \copydetails simdInvsqrtPairF
 */
static inline void gmx_simdcall
simdInvsqrtPair(SimdReal x0,    SimdReal x1,
                SimdReal *out0, SimdReal *out1)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdInvsqrtPairSingleAccuracyD(x0, x1, out0, out1);
#elif defined GMX_DOUBLE
    return simdInvsqrtPairD(x0, x1, out0, out1);
#else
    return simdInvsqrtPairF(x0, x1, out0, out1);
#endif
}

/*! \brief Calculate 1/sqrt(x) for SIMD real.
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdReal gmx_simdcall
simdSqrt(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdSqrtSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdSqrtD(x);
#else
    return simdSqrtF(x);
#endif
}

/*! \brief Calculate 1/sqrt(x) for SIMD real.
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdReal gmx_simdcall
simdInv(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdInvSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdInvD(x);
#else
    return simdInvF(x);
#endif
}

/*! \brief Calculate 1/sqrt(x) for masked entries of SIMD reals.
 *
 *  \copydetails simdInvsqrtMaskF
 */
static inline SimdReal
simdInvsqrtMask(SimdReal x, SimdBool m)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdInvsqrtMaskSingleAccuracyD(x, m);
#elif defined GMX_DOUBLE
    return simdInvsqrtMaskD(x, m);
#else
    return simdInvsqrtMaskF(x, m);
#endif
}

/*! \brief Calculate 1/x for masked entries of SIMD reals.
 *
 *  \copydetails simdInvMaskF
 */
static inline SimdReal
simdInvMask(SimdReal x, SimdBool m)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdInvMaskSingleAccuracyD(x, m);
#elif defined GMX_DOUBLE
    return simdInvMaskD(x, m);
#else
    return simdInvMaskF(x, m);
#endif
}

/*! \brief SIMD real log(x). This is the natural logarithm.
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdReal gmx_simdcall
simdLog(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdLogSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdLogD(x);
#else
    return simdLogF(x);
#endif
}

/*! \brief SIMD real 2^x.
 *
 * \copydetails simdExp2F
 */
static inline SimdReal gmx_simdcall
simdExp2(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdExp2SingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdExp2D(x);
#else
    return simdExp2F(x);
#endif
}

/*! \brief SIMD real e^x.
 *
 * \copydetails simdExpF
 */
static inline SimdReal gmx_simdcall
simdExp(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdExpSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdExpD(x);
#else
    return simdExpF(x);
#endif
}

/*! \brief SIMD real erf(x).
 *
 * \copydetails simdErfF
 */
static inline SimdReal gmx_simdcall
simdErf(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdErfSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdErfD(x);
#else
    return simdErfF(x);
#endif
}

/*! \brief SIMD real erfc(x).
 *
 * \copydetails simdErfcF
 */
static inline SimdReal gmx_simdcall
simdErfc(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdErfcSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdErfcD(x);
#else
    return simdErfcF(x);
#endif
}

/*! \brief SIMD real sin \& cos.
 *
 * \copydetails simdSinCosF
 */
static inline void gmx_simdcall
simdSinCos(SimdReal x, SimdReal *sinval, SimdReal *cosval)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    simdSinCosSingleAccuracyD(x, sinval, cosval);
#elif defined GMX_DOUBLE
    simdSinCosD(x, sinval, cosval);
#else
    simdSinCosF(x, sinval, cosval);
#endif
}

/*! \brief SIMD real sin(x).
 *
 * \copydetails simdSinF
 */
static inline SimdReal gmx_simdcall
simdSin(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdSinSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdSinD(x);
#else
    return simdSinF(x);
#endif
}

/*! \brief SIMD real cos(x).
 *
 * \copydetails simdCosF
 */
static inline SimdReal gmx_simdcall
simdCos(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdCosSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdCosD(x);
#else
    return simdCosF(x);
#endif
}

/*! \brief SIMD real tan(x).
 *
 * \copydetails simdTanF
 */
static inline SimdReal gmx_simdcall
simdTan(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdTanSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdTanD(x);
#else
    return simdTanF(x);
#endif
}

/*! \brief SIMD real asin(x).
 *
 * \copydetails simdAsinF
 */
static inline SimdReal gmx_simdcall
simdAsin(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdAsinSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdAsinD(x);
#else
    return simdAsinF(x);
#endif
}

/*! \brief SIMD real acos(x).
 *
 * \copydetails simdAcosF
 */
static inline SimdReal gmx_simdcall
simdAcos(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdAcosSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdAcosD(x);
#else
    return simdAcosF(x);
#endif
}

/*! \brief SIMD real atan(x).
 *
 * \copydetails simdAtanF
 */
static inline SimdReal gmx_simdcall
simdAtan(SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdAtanSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simdAtanD(x);
#else
    return simdAtanF(x);
#endif
}

/*! \brief SIMD real atan2(y,x).
 *
 * \copydetails simdAtan2F
 */
static inline SimdReal gmx_simdcall
simdAtan2(SimdReal y, SimdReal x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdAtan2SingleAccuracyD(y, x);
#elif defined GMX_DOUBLE
    return simdAtan2D(y, x);
#else
    return simdAtan2F(y, x);
#endif
}

/*! \brief SIMD Analytic PME force correction.
 *
 * \copydetails simdPmeCorrForceF
 */
static inline SimdReal gmx_simdcall
simdPmeCorrForce(SimdReal z2)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdPmeCorrForceSingleAccuracyD(z2);
#elif defined GMX_DOUBLE
    return simdPmeCorrForceD(z2);
#else
    return simdPmeCorrForceF(z2);
#endif
}

/*! \brief SIMD Analytic PME potential correction.
 *
 * \copydetails simdPmeCorrPotentialF
 */
static inline SimdReal gmx_simdcall
simdPmeCorrPotential(SimdReal z2)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simdPmeCorrPotentialSingleAccuracyD(z2);
#elif defined GMX_DOUBLE
    return simdPmeCorrPotentialD(z2);
#else
    return simdPmeCorrPotentialF(z2);
#endif
}

#if GMX_SIMD4_HAVE_REAL
/*! \brief Calculate 1/sqrt(x) for SIMD4 real.
 *
 * \copydetails simdInvsqrtF
 */
static inline Simd4Real gmx_simdcall
simd4Invsqrt(Simd4Real x)
{
#if (defined GMX_DOUBLE) && (GMX_SIMD_ACCURACY_BITS_DOUBLE <= GMX_SIMD_ACCURACY_BITS_SINGLE)
    return simd4InvsqrtSingleAccuracyD(x);
#elif defined GMX_DOUBLE
    return simd4InvsqrtD(x);
#else
    return simd4InvsqrtF(x);
#endif
}
#endif

/*! \brief Calculate 1/sqrt(x) for SIMD real, only targeting single accuracy.
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdReal gmx_simdcall
simdInvsqrtSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdInvsqrtSingleAccuracyD(x);
#else
    return simdInvsqrtF(x);
#endif
}

/*! \brief Calculate 1/sqrt(x) for two SIMD reals, only targeting single accuracy.
 *
 * \copydetails simdInvsqrtPairF
 */
static inline void gmx_simdcall
simdInvsqrtPairSingleAccuracy(SimdReal x0,    SimdReal x1,
                              SimdReal *out0, SimdReal *out1)
{
#ifdef GMX_DOUBLE
    return simdInvsqrtPairSingleAccuracyD(x0, x1, out0, out1);
#else
    return simdInvsqrtPairF(x0, x1, out0, out1);
#endif
}

/*! \brief Calculate 1/sqrt(x) for SIMD real, only targeting single accuracy.
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdReal gmx_simdcall
simdSqrtSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdSqrtSingleAccuracyD(x);
#else
    return simdSqrtF(x);
#endif
}

/*! \brief Calculate 1/sqrt(x) for SIMD real, only targeting single accuracy.
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdReal gmx_simdcall
simdInvSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdInvSingleAccuracyD(x);
#else
    return simdInvF(x);
#endif
}

/*! \brief Calculate 1/sqrt(x) for masked SIMD reals, only targeting single accuracy.
 *
 *  \copydetails simdInvsqrtMaskF
 */
static inline SimdReal
simdInvsqrtMaskSingleAccuracy(SimdReal x, SimdBool m)
{
#ifdef GMX_DOUBLE
    return simdInvsqrtMaskSingleAccuracyD(x, m);
#else
    return simdInvsqrtMaskF(x, m);
#endif
}

/*! \brief Calculate 1/x for masked SIMD reals, only targeting single accuracy.
 *
 *  \copydetails simdInvMaskF
 */
static inline SimdReal
simdInvMaskSingleAccuracy(SimdReal x, SimdBool m)
{
#ifdef GMX_DOUBLE
    return simdInvMaskSingleAccuracyD(x, m);
#else
    return simdInvMaskF(x, m);
#endif
}

/*! \brief SIMD real log(x), only targeting single accuracy. This is the natural logarithm.
 *
 * \copydetails simdInvsqrtF
 */
static inline SimdReal gmx_simdcall
simdLogSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdLogSingleAccuracyD(x);
#else
    return simdLogF(x);
#endif
}

/*! \brief SIMD real 2^x, only targeting single accuracy.
 *
 * \copydetails simdExp2F
 */
static inline SimdReal gmx_simdcall
simdExp2SingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdExp2SingleAccuracyD(x);
#else
    return simdExp2F(x);
#endif
}

/*! \brief SIMD real e^x, only targeting single accuracy.
 *
 * \copydetails simdExpF
 */
static inline SimdReal gmx_simdcall
simdExpSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdExpSingleAccuracyD(x);
#else
    return simdExpF(x);
#endif
}

/*! \brief SIMD real erf(x), only targeting single accuracy.
 *
 * \copydetails simdErfF
 */
static inline SimdReal gmx_simdcall
simdErfSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdErfSingleAccuracyD(x);
#else
    return simdErfF(x);
#endif
}

/*! \brief SIMD real erfc(x), only targeting single accuracy.
 *
 * \copydetails simdErfcF
 */
static inline SimdReal gmx_simdcall
simdErfcSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdErfcSingleAccuracyD(x);
#else
    return simdErfcF(x);
#endif
}

/*! \brief SIMD real sin \& cos, only targeting single accuracy.
 *
 * \copydetails simdSinCosF
 */
static inline void gmx_simdcall
simdSinCosSingleAccuracy(SimdReal x, SimdReal *sinval, SimdReal *cosval)
{
#ifdef GMX_DOUBLE
    simdSinCosSingleAccuracyD(x, sinval, cosval);
#else
    simdSinCosF(x, sinval, cosval);
#endif
}

/*! \brief SIMD real sin(x), only targeting single accuracy.
 *
 * \copydetails simdSinF
 */
static inline SimdReal gmx_simdcall
simdSinSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdSinSingleAccuracyD(x);
#else
    return simdSinF(x);
#endif
}

/*! \brief SIMD real cos(x), only targeting single accuracy.
 *
 * \copydetails simdCosF
 */
static inline SimdReal gmx_simdcall
simdCosSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdCosSingleAccuracyD(x);
#else
    return simdCosF(x);
#endif
}

/*! \brief SIMD real tan(x), only targeting single accuracy.
 *
 * \copydetails simdTanF
 */
static inline SimdReal gmx_simdcall
simdTanSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdTanSingleAccuracyD(x);
#else
    return simdTanF(x);
#endif
}

/*! \brief SIMD real asin(x), only targeting single accuracy.
 *
 * \copydetails simdAsinF
 */
static inline SimdReal gmx_simdcall
simdAsinSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdAsinSingleAccuracyD(x);
#else
    return simdAsinF(x);
#endif
}

/*! \brief SIMD real acos(x), only targeting single accuracy.
 *
 * \copydetails simdAcosF
 */
static inline SimdReal gmx_simdcall
simdAcosSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdAcosSingleAccuracyD(x);
#else
    return simdAcosF(x);
#endif
}

/*! \brief SIMD real atan(x), only targeting single accuracy.
 *
 * \copydetails simdAtanF
 */
static inline SimdReal gmx_simdcall
simdAtanSingleAccuracy(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdAtanSingleAccuracyD(x);
#else
    return simdAtanF(x);
#endif
}

/*! \brief SIMD real atan2(y,x), only targeting single accuracy.
 *
 * \copydetails simdAtan2F
 */
static inline SimdReal gmx_simdcall
simdAtan2SingleAccuracy(SimdReal y, SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdAtan2SingleAccuracyD(y, x);
#else
    return simdAtan2F(y, x);
#endif
}

/*! \brief SIMD Analytic PME force correction, only targeting single accuracy.
 *
 * \copydetails simdPmeCorrForceF
 */
static inline SimdReal gmx_simdcall
simdPmeCorrForceSingleAccuracy(SimdReal z2)
{
#ifdef GMX_DOUBLE
    return simdPmeCorrForceSingleAccuracyD(z2);
#else
    return simdPmeCorrForceF(z2);
#endif
}

/*! \brief SIMD Analytic PME potential correction, only targeting single accuracy.
 *
 * \copydetails simdPmeCorrPotentialF
 */
static inline SimdReal gmx_simdcall
simdPmeCorrPotentialSingleAccuracy(SimdReal z2)
{
#ifdef GMX_DOUBLE
    return simdPmeCorrPotentialSingleAccuracyD(z2);
#else
    return simdPmeCorrPotentialF(z2);
#endif
}

#if GMX_SIMD4_HAVE_REAL
/*! \brief Calculate 1/sqrt(x) for SIMD4 real, only targeting single accuracy.
 *
 * \copydetails simdInvsqrtF
 */
static inline Simd4Real gmx_simdcall
simd4InvsqrtSingleAccuracy(Simd4Real x)
{
#ifdef GMX_DOUBLE
    return simd4InvsqrtSingleAccuracyD(x);
#else
    return simd4InvsqrtF(x);
#endif
}
#endif // GMX_SIMD4_HAVE_REAL

/*! \}   end of addtogroup module_simd */
/*! \endcond  end of condition libabl */

#endif /* GMX_SIMD */

}      // namespace gmx

#endif /* GMX_SIMD_SIMD_MATH_H */
