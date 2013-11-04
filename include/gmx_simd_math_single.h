/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
#ifndef _gmx_simd_math_single_h_
#define _gmx_simd_math_single_h_


/* 1.0/sqrt(x) */
static gmx_inline gmx_mm_pr
gmx_invsqrt_pr(gmx_mm_pr x)
{
    /* This is one of the few cases where FMA adds a FLOP, but ends up with
     * less instructions in total when FMA is available in hardware.
     * Usually we would not optimize this far, but invsqrt is used often.
     */
#ifdef GMX_SIMD_HAVE_FMA
    const gmx_mm_pr half  = gmx_set1_pr(0.5);
    const gmx_mm_pr one   = gmx_set1_pr(1.0);

    gmx_mm_pr       lu = gmx_rsqrt_pr(x);

    return gmx_madd_pr(gmx_nmsub_pr(x, gmx_mul_pr(lu, lu), one), gmx_mul_pr(lu, half), lu);
#else
    const gmx_mm_pr half  = gmx_set1_pr(0.5);
    const gmx_mm_pr three = gmx_set1_pr(3.0);

    gmx_mm_pr       lu = gmx_rsqrt_pr(x);
    
    return gmx_mul_pr(half, gmx_mul_pr(gmx_sub_pr(three, gmx_mul_pr(gmx_mul_pr(lu, lu), x)), lu));
#endif
}


/* 1.0/x */
static gmx_inline gmx_mm_pr
gmx_inv_pr(gmx_mm_pr x)
{
    const gmx_mm_pr two = gmx_set1_pr(2.0);

    gmx_mm_pr       lu = gmx_rcp_pr(x);

    return gmx_mul_pr(lu, gmx_nmsub_pr(lu, x, two));
}


/* Calculate the force correction due to PME analytically.
 *
 * This routine is meant to enable analytical evaluation of the
 * direct-space PME electrostatic force to avoid tables.
 *
 * The direct-space potential should be Erfc(beta*r)/r, but there
 * are some problems evaluating that:
 *
 * First, the error function is difficult (read: expensive) to
 * approxmiate accurately for intermediate to large arguments, and
 * this happens already in ranges of beta*r that occur in simulations.
 * Second, we now try to avoid calculating potentials in Gromacs but
 * use forces directly.
 *
 * We can simply things slight by noting that the PME part is really
 * a correction to the normal Coulomb force since Erfc(z)=1-Erf(z), i.e.
 *
 * V= 1/r - Erf(beta*r)/r
 *
 * The first term we already have from the inverse square root, so
 * that we can leave out of this routine.
 *
 * For pme tolerances of 1e-3 to 1e-8 and cutoffs of 0.5nm to 1.8nm,
 * the argument beta*r will be in the range 0.15 to ~4. Use your
 * favorite plotting program to realize how well-behaved Erf(z)/z is
 * in this range!
 *
 * We approximate f(z)=erf(z)/z with a rational minimax polynomial.
 * However, it turns out it is more efficient to approximate f(z)/z and
 * then only use even powers. This is another minor optimization, since
 * we actually WANT f(z)/z, because it is going to be multiplied by
 * the vector between the two atoms to get the vectorial force. The
 * fastest flops are the ones we can avoid calculating!
 *
 * So, here's how it should be used:
 *
 * 1. Calculate r^2.
 * 2. Multiply by beta^2, so you get z^2=beta^2*r^2.
 * 3. Evaluate this routine with z^2 as the argument.
 * 4. The return value is the expression:
 *
 *
 *       2*exp(-z^2)     erf(z)
 *       ------------ - --------
 *       sqrt(Pi)*z^2      z^3
 *
 * 5. Multiply the entire expression by beta^3. This will get you
 *
 *       beta^3*2*exp(-z^2)     beta^3*erf(z)
 *       ------------------  - ---------------
 *          sqrt(Pi)*z^2            z^3
 *
 *    or, switching back to r (z=r*beta):
 *
 *       2*beta*exp(-r^2*beta^2)   erf(r*beta)
 *       ----------------------- - -----------
 *            sqrt(Pi)*r^2            r^3
 *
 *
 *    With a bit of math exercise you should be able to confirm that
 *    this is exactly D[Erf[beta*r]/r,r] divided by r another time.
 *
 * 6. Add the result to 1/r^3, multiply by the product of the charges,
 *    and you have your force (divided by r). A final multiplication
 *    with the vector connecting the two particles and you have your
 *    vectorial force to add to the particles.
 *
 */
static gmx_mm_pr
gmx_pmecorrF_pr(gmx_mm_pr z2)
{
    const gmx_mm_pr  FN6      = gmx_set1_pr(-1.7357322914161492954e-8f);
    const gmx_mm_pr  FN5      = gmx_set1_pr(1.4703624142580877519e-6f);
    const gmx_mm_pr  FN4      = gmx_set1_pr(-0.000053401640219807709149f);
    const gmx_mm_pr  FN3      = gmx_set1_pr(0.0010054721316683106153f);
    const gmx_mm_pr  FN2      = gmx_set1_pr(-0.019278317264888380590f);
    const gmx_mm_pr  FN1      = gmx_set1_pr(0.069670166153766424023f);
    const gmx_mm_pr  FN0      = gmx_set1_pr(-0.75225204789749321333f);

    const gmx_mm_pr  FD4      = gmx_set1_pr(0.0011193462567257629232f);
    const gmx_mm_pr  FD3      = gmx_set1_pr(0.014866955030185295499f);
    const gmx_mm_pr  FD2      = gmx_set1_pr(0.11583842382862377919f);
    const gmx_mm_pr  FD1      = gmx_set1_pr(0.50736591960530292870f);
    const gmx_mm_pr  FD0      = gmx_set1_pr(1.0f);

    gmx_mm_pr        z4;
    gmx_mm_pr        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = gmx_mul_pr(z2, z2);

    polyFD0        = gmx_madd_pr(FD4, z4, FD2);
    polyFD1        = gmx_madd_pr(FD3, z4, FD1);
    polyFD0        = gmx_madd_pr(polyFD0, z4, FD0);
    polyFD0        = gmx_madd_pr(polyFD1, z2, polyFD0);

    polyFD0        = gmx_inv_pr(polyFD0);

    polyFN0        = gmx_madd_pr(FN6, z4, FN4);
    polyFN1        = gmx_madd_pr(FN5, z4, FN3);
    polyFN0        = gmx_madd_pr(polyFN0, z4, FN2);
    polyFN1        = gmx_madd_pr(polyFN1, z4, FN1);
    polyFN0        = gmx_madd_pr(polyFN0, z4, FN0);
    polyFN0        = gmx_madd_pr(polyFN1, z2, polyFN0);

    return gmx_mul_pr(polyFN0, polyFD0);
}


/* Calculate the potential correction due to PME analytically.
 *
 * See gmx_pmecorrF_pr() for details about the approximation.
 *
 * This routine calculates Erf(z)/z, although you should provide z^2
 * as the input argument.
 *
 * Here's how it should be used:
 *
 * 1. Calculate r^2.
 * 2. Multiply by beta^2, so you get z^2=beta^2*r^2.
 * 3. Evaluate this routine with z^2 as the argument.
 * 4. The return value is the expression:
 *
 *
 *        erf(z)
 *       --------
 *          z
 *
 * 5. Multiply the entire expression by beta and switching back to r (z=r*beta):
 *
 *       erf(r*beta)
 *       -----------
 *           r
 *
 * 6. Add the result to 1/r, multiply by the product of the charges,
 *    and you have your potential.
 */
static gmx_mm_pr
gmx_pmecorrV_pr(gmx_mm_pr z2)
{
    const gmx_mm_pr  VN6      = gmx_set1_pr(1.9296833005951166339e-8f);
    const gmx_mm_pr  VN5      = gmx_set1_pr(-1.4213390571557850962e-6f);
    const gmx_mm_pr  VN4      = gmx_set1_pr(0.000041603292906656984871f);
    const gmx_mm_pr  VN3      = gmx_set1_pr(-0.00013134036773265025626f);
    const gmx_mm_pr  VN2      = gmx_set1_pr(0.038657983986041781264f);
    const gmx_mm_pr  VN1      = gmx_set1_pr(0.11285044772717598220f);
    const gmx_mm_pr  VN0      = gmx_set1_pr(1.1283802385263030286f);

    const gmx_mm_pr  VD3      = gmx_set1_pr(0.0066752224023576045451f);
    const gmx_mm_pr  VD2      = gmx_set1_pr(0.078647795836373922256f);
    const gmx_mm_pr  VD1      = gmx_set1_pr(0.43336185284710920150f);
    const gmx_mm_pr  VD0      = gmx_set1_pr(1.0f);

    gmx_mm_pr        z4;
    gmx_mm_pr        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = gmx_mul_pr(z2, z2);

    polyVD1        = gmx_madd_pr(VD3, z4, VD1);
    polyVD0        = gmx_madd_pr(VD2, z4, VD0);
    polyVD0        = gmx_madd_pr(polyVD1, z2, polyVD0);

    polyVD0        = gmx_inv_pr(polyVD0);

    polyVN0        = gmx_madd_pr(VN6, z4, VN4);
    polyVN1        = gmx_madd_pr(VN5, z4, VN3);
    polyVN0        = gmx_madd_pr(polyVN0, z4, VN2);
    polyVN1        = gmx_madd_pr(polyVN1, z4, VN1);
    polyVN0        = gmx_madd_pr(polyVN0, z4, VN0);
    polyVN0        = gmx_madd_pr(polyVN1, z2, polyVN0);

    return gmx_mul_pr(polyVN0, polyVD0);
}


#endif /* _gmx_simd_math_single_h_ */
