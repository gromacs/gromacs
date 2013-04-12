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
#ifndef _gmx_simd_math_double_h_
#define _gmx_simd_math_double_h_


/* 1.0/sqrt(x) */
static gmx_inline gmx_mm_pr
gmx_invsqrt_pr(gmx_mm_pr x)
{
    const gmx_mm_pr half  = gmx_set1_pr(0.5);
    const gmx_mm_pr three = gmx_set1_pr(3.0);

    /* Lookup instruction only exists in single precision, convert back and forth... */
    gmx_mm_pr lu = gmx_rsqrt_pr(x);

    lu = gmx_mul_pr(gmx_mul_pr(half, lu), gmx_nmsub_pr(gmx_mul_pr(lu, lu), x, three));
    return gmx_mul_pr(gmx_mul_pr(half, lu), gmx_nmsub_pr(gmx_mul_pr(lu, lu), x, three));
}


/* 1.0/x */
static gmx_inline gmx_mm_pr
gmx_inv_pr(gmx_mm_pr x)
{
    const gmx_mm_pr two  = gmx_set1_pr(2.0);

    /* Lookup instruction only exists in single precision, convert back and forth... */
    gmx_mm_pr lu = gmx_rcp_pr(x);

    /* Perform two N-R steps for double precision */
    lu         = gmx_mul_pr(lu, gmx_nmsub_pr(lu, x, two));
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
    const gmx_mm_pr  FN10     = gmx_set1_pr(-8.0072854618360083154e-14);
    const gmx_mm_pr  FN9      = gmx_set1_pr(1.1859116242260148027e-11);
    const gmx_mm_pr  FN8      = gmx_set1_pr(-8.1490406329798423616e-10);
    const gmx_mm_pr  FN7      = gmx_set1_pr(3.4404793543907847655e-8);
    const gmx_mm_pr  FN6      = gmx_set1_pr(-9.9471420832602741006e-7);
    const gmx_mm_pr  FN5      = gmx_set1_pr(0.000020740315999115847456);
    const gmx_mm_pr  FN4      = gmx_set1_pr(-0.00031991745139313364005);
    const gmx_mm_pr  FN3      = gmx_set1_pr(0.0035074449373659008203);
    const gmx_mm_pr  FN2      = gmx_set1_pr(-0.031750380176100813405);
    const gmx_mm_pr  FN1      = gmx_set1_pr(0.13884101728898463426);
    const gmx_mm_pr  FN0      = gmx_set1_pr(-0.75225277815249618847);

    const gmx_mm_pr  FD5      = gmx_set1_pr(0.000016009278224355026701);
    const gmx_mm_pr  FD4      = gmx_set1_pr(0.00051055686934806966046);
    const gmx_mm_pr  FD3      = gmx_set1_pr(0.0081803507497974289008);
    const gmx_mm_pr  FD2      = gmx_set1_pr(0.077181146026670287235);
    const gmx_mm_pr  FD1      = gmx_set1_pr(0.41543303143712535988);
    const gmx_mm_pr  FD0      = gmx_set1_pr(1.0);

    gmx_mm_pr        z4;
    gmx_mm_pr        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = gmx_mul_pr(z2, z2);

    polyFD1        = gmx_madd_pr(FD5, z4, FD3);
    polyFD1        = gmx_madd_pr(polyFD1, z4, FD1);
    polyFD1        = gmx_mul_pr(polyFD1, z2);
    polyFD0        = gmx_madd_pr(FD4, z4, FD2);
    polyFD0        = gmx_madd_pr(polyFD0, z4, FD0);
    polyFD0        = gmx_add_pr(polyFD0, polyFD1);

    polyFD0        = gmx_inv_pr(polyFD0);

    polyFN0        = gmx_madd_pr(FN10, z4, FN8);
    polyFN0        = gmx_madd_pr(polyFN0, z4, FN6);
    polyFN0        = gmx_madd_pr(polyFN0, z4, FN4);
    polyFN0        = gmx_madd_pr(polyFN0, z4, FN2);
    polyFN0        = gmx_madd_pr(polyFN0, z4, FN0);
    polyFN1        = gmx_madd_pr(FN9, z4, FN7);
    polyFN1        = gmx_madd_pr(polyFN1, z4, FN5);
    polyFN1        = gmx_madd_pr(polyFN1, z4, FN3);
    polyFN1        = gmx_madd_pr(polyFN1, z4, FN1);
    polyFN0        = gmx_madd_pr(polyFN1, z2, polyFN0);

    return gmx_mul_pr(polyFN0, polyFD0);
}


/* Calculate the potential correction due to PME analytically.
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
 * 6. Subtract the result from 1/r, multiply by the product of the charges,
 *    and you have your potential.
 *
 */
static gmx_mm_pr
gmx_pmecorrV_pr(gmx_mm_pr z2)
{
    const gmx_mm_pr  VN9      = gmx_set1_pr(-9.3723776169321855475e-13);
    const gmx_mm_pr  VN8      = gmx_set1_pr(1.2280156762674215741e-10);
    const gmx_mm_pr  VN7      = gmx_set1_pr(-7.3562157912251309487e-9);
    const gmx_mm_pr  VN6      = gmx_set1_pr(2.6215886208032517509e-7);
    const gmx_mm_pr  VN5      = gmx_set1_pr(-4.9532491651265819499e-6);
    const gmx_mm_pr  VN4      = gmx_set1_pr(0.00025907400778966060389);
    const gmx_mm_pr  VN3      = gmx_set1_pr(0.0010585044856156469792);
    const gmx_mm_pr  VN2      = gmx_set1_pr(0.045247661136833092885);
    const gmx_mm_pr  VN1      = gmx_set1_pr(0.11643931522926034421);
    const gmx_mm_pr  VN0      = gmx_set1_pr(1.1283791671726767970);

    const gmx_mm_pr  VD5      = gmx_set1_pr(0.000021784709867336150342);
    const gmx_mm_pr  VD4      = gmx_set1_pr(0.00064293662010911388448);
    const gmx_mm_pr  VD3      = gmx_set1_pr(0.0096311444822588683504);
    const gmx_mm_pr  VD2      = gmx_set1_pr(0.085608012351550627051);
    const gmx_mm_pr  VD1      = gmx_set1_pr(0.43652499166614811084);
    const gmx_mm_pr  VD0      = gmx_set1_pr(1.0);

    gmx_mm_pr        z4;
    gmx_mm_pr        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = gmx_mul_pr(z2, z2);

    polyVD1        = gmx_madd_pr(VD5, z4, VD3);
    polyVD0        = gmx_madd_pr(VD4, z4, VD2);
    polyVD1        = gmx_madd_pr(polyVD1, z4, VD1);
    polyVD0        = gmx_madd_pr(polyVD0, z4, VD0);
    polyVD0        = gmx_madd_pr(polyVD1, z2, polyVD0);

    polyVD0        = gmx_inv_pr(polyVD0);

    polyVN1        = gmx_madd_pr(VN9, z4, VN7);
    polyVN0        = gmx_madd_pr(VN8, z4, VN6);
    polyVN1        = gmx_madd_pr(polyVN1, z4, VN5);
    polyVN0        = gmx_madd_pr(polyVN0, z4, VN4);
    polyVN1        = gmx_madd_pr(polyVN1, z4, VN3);
    polyVN0        = gmx_madd_pr(polyVN0, z4, VN2);
    polyVN1        = gmx_madd_pr(polyVN1, z4, VN1);
    polyVN0        = gmx_madd_pr(polyVN0, z4, VN0);
    polyVN0        = gmx_madd_pr(polyVN1, z2, polyVN0);

    return gmx_mul_pr(polyVN0, polyVD0);
}


#endif /*_gmx_simd_math_double_h_ */
