/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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
#ifndef _gmx_math_ibm_qpx_single_h_
#define _gmx_math_ibm_qpx_single_h_

#include <math.h>


#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288
#endif


/************************
 *                      *
 * Simple math routines *
 *                      *
 ************************/

/* 1.0/sqrt(x) */
static gmx_inline vector4double
gmx_ibm_qpx_invsqrt(vector4double x)
{
    const vector4double halfthree = vec_splats(1.5);

    /* Can we hope for 12 bits in practice? Otherwise we need more iterations... */
    vector4double lu = vec_rsqrte(x);

    /* halfx = 1.5 * x - x */
    vector4double halfx = vec_msub(halfthree,x,x);

    /* Each line does a single Newton-Raphson iterate as
    lu = lu * (1.5 - halfx * lu * lu) */
#ifdef GMX_DOUBLE
    lu = vec_mul(lu, vec_nmsub(halfx,vec_mul(lu,lu), halfthree));
#endif
    /* Need next line if 12 bits are not available in practice
    lu = vec_mul(lu, vec_nmsub(halfx,vec_mul(lu,lu), halfthree));
    */
    return vec_mul(lu, vec_nmsub(halfx,vec_mul(lu,lu), halfthree));
}


/* 1.0/x */
static gmx_inline vector4double
gmx_ibm_qpx_inv(vector4double x)
{
    const vector4double two = vec_splats(2.0);

    /* Can we hope for 12 bits in practice? Otherwise we need more iterations... */
    vector4double lu = vec_re(x);

#ifdef GMX_DOUBLE
    lu = vec_mul(lu,vec_nmsub(lu,x,two));
#endif
    return vec_mul(lu,vec_nmsub(lu,x,two));
}




/* Calculate the force correction due to PME analytically.
 *
 * See documentation in the SSE headers for details about algorithm and usage.
 */
static vector4double
gmx_ibm_qpx_pmecorrF(vector4double z2)
{
#ifdef GMX_DOUBLE
    const __m128d  FN10     = vec_splats(-8.0072854618360083154e-14);
    const __m128d  FN9      = vec_splats(1.1859116242260148027e-11);
    const __m128d  FN8      = vec_splats(-8.1490406329798423616e-10);
    const __m128d  FN7      = vec_splats(3.4404793543907847655e-8);
    const __m128d  FN6      = vec_splats(-9.9471420832602741006e-7);
    const __m128d  FN5      = vec_splats(0.000020740315999115847456);
    const __m128d  FN4      = vec_splats(-0.00031991745139313364005);
    const __m128d  FN3      = vec_splats(0.0035074449373659008203);
    const __m128d  FN2      = vec_splats(-0.031750380176100813405);
    const __m128d  FN1      = vec_splats(0.13884101728898463426);
    const __m128d  FN0      = vec_splats(-0.75225277815249618847);
    
    const __m128d  FD5      = vec_splats(0.000016009278224355026701);
    const __m128d  FD4      = vec_splats(0.00051055686934806966046);
    const __m128d  FD3      = vec_splats(0.0081803507497974289008);
    const __m128d  FD2      = vec_splats(0.077181146026670287235);
    const __m128d  FD1      = vec_splats(0.41543303143712535988);
    const __m128d  FD0      = vec_splats(1.0);
    
    __m128d z4;
    __m128d polyFN0,polyFN1,polyFD0,polyFD1;
    
    z4             = vec_mul(z2,z2);
    
    polyFD1        = vec_madd(FD5,z4,FD3);
    polyFD1        = vec_madd(polyFD1,z4,FD1);
    polyFD1        = vec_mul(polyFD1,z2);
    polyFD0        = vec_madd(FD4,z4,FD2);
    polyFD0        = vec_madd(polyFD0,z4,FD0);
    polyFD0        = vec_add(polyFD0,polyFD1);
    
    polyFD0        = gmx_ibm_qpx_inv(polyFD0);
    
    polyFN0        = vec_madd(FN10,z4,FN8);
    polyFN0        = vec_madd(polyFN0,z4,FN6);
    polyFN0        = vec_madd(polyFN0,z4,FN4);
    polyFN0        = vec_madd(polyFN0,z4,FN2);
    polyFN0        = vec_madd(polyFN0,z4,FN0);
    polyFN1        = vec_madd(FN9,z4,FN7);
    polyFN1        = vec_madd(polyFN1,z4,FN5);
    polyFN1        = vec_madd(polyFN1,z4,FN3);
    polyFN1        = vec_madd(polyFN1,z4,FN1);
    polyFN0        = vec_madd(polyFN1,z2,polyFN0);
#else
    const vector4double  FN6      = vec_splats(-1.7357322914161492954e-8);
    const vector4double  FN5      = vec_splats(1.4703624142580877519e-6);
    const vector4double  FN4      = vec_splats(-0.000053401640219807709149);
    const vector4double  FN3      = vec_splats(0.0010054721316683106153);
    const vector4double  FN2      = vec_splats(-0.019278317264888380590);
    const vector4double  FN1      = vec_splats(0.069670166153766424023);
    const vector4double  FN0      = vec_splats(-0.75225204789749321333);
    
    const vector4double  FD4      = vec_splats(0.0011193462567257629232);
    const vector4double  FD3      = vec_splats(0.014866955030185295499);
    const vector4double  FD2      = vec_splats(0.11583842382862377919);
    const vector4double  FD1      = vec_splats(0.50736591960530292870);
    const vector4double  FD0      = vec_splats(1.0);
    
    vector4double z4;
    vector4double polyFN0,polyFN1,polyFD0,polyFD1;
    
    z4             = vec_mul(z2,z2);
    
    polyFD0        = vec_madd(FD4,z4,FD2);
    polyFD1        = vec_madd(FD3,z4,FD1);
    polyFD0        = vec_madd(polyFD0,z4,FD0);
    polyFD0        = vec_madd(polyFD1,z2,polyFD0);
    
    polyFD0        = gmx_ibm_qpx_inv(polyFD0);
    
    polyFN0        = vec_madd(FN6,z4,FN4);
    polyFN1        = vec_madd(FN5,z4,FN3);
    polyFN0        = vec_madd(polyFN0,z4,FN2);
    polyFN1        = vec_madd(polyFN1,z4,FN1);
    polyFN0        = vec_madd(polyFN0,z4,FN0);
    polyFN0        = vec_madd(polyFN1,z2,polyFN0);
#endif
    return   vec_mul(polyFN0,polyFD0);
}




/* Calculate the potential correction due to PME analytically.
 *
 * See documentation in the SSE headers for details about algorithm and usage.
 */
static vector4double
gmx_ibm_qpx_pmecorrV(vector4double z2)
{
#ifdef GMX_DOUBLE
    const __m128d  VN9      = vec_splats(-9.3723776169321855475e-13);
    const __m128d  VN8      = vec_splats(1.2280156762674215741e-10);
    const __m128d  VN7      = vec_splats(-7.3562157912251309487e-9);
    const __m128d  VN6      = vec_splats(2.6215886208032517509e-7);
    const __m128d  VN5      = vec_splats(-4.9532491651265819499e-6);
    const __m128d  VN4      = vec_splats(0.00025907400778966060389);
    const __m128d  VN3      = vec_splats(0.0010585044856156469792);
    const __m128d  VN2      = vec_splats(0.045247661136833092885);
    const __m128d  VN1      = vec_splats(0.11643931522926034421);
    const __m128d  VN0      = vec_splats(1.1283791671726767970);
    
    const __m128d  VD5      = vec_splats(0.000021784709867336150342);
    const __m128d  VD4      = vec_splats(0.00064293662010911388448);
    const __m128d  VD3      = vec_splats(0.0096311444822588683504);
    const __m128d  VD2      = vec_splats(0.085608012351550627051);
    const __m128d  VD1      = vec_splats(0.43652499166614811084);
    const __m128d  VD0      = vec_splats(1.0);
    
    __m128d z4;
    __m128d polyVN0,polyVN1,polyVD0,polyVD1;
    
    z4             = vec_mul(z2,z2);
    
    polyVD1        = vec_madd(VD5,z4,VD3);
    polyVD0        = vec_madd(VD4,z4,VD2);
    polyVD1        = vec_madd(polyVD1,z4,VD1);
    polyVD0        = vec_madd(polyVD0,z4,VD0);
    polyVD0        = vec_madd(polyVD1,z2,polyVD0);
    
    polyVD0        = gmx_ibm_qpx_inv(polyVD0);
    
    polyVN1        = vec_madd(VN9,z4,VN7);
    polyVN0        = vec_madd(VN8,z4,VN6);
    polyVN1        = vec_madd(polyVN1,z4,VN5);
    polyVN0        = vec_madd(polyVN0,z4,VN4);
    polyVN1        = vec_madd(polyVN1,z4,VN3);
    polyVN0        = vec_madd(polyVN0,z4,VN2);
    polyVN1        = vec_madd(polyVN1,z4,VN1);
    polyVN0        = vec_madd(polyVN0,z4,VN0);
    polyVN0        = vec_madd(polyVN1,z2,polyVN0);
#else
    const vector4double  VN6      = vec_splats(1.9296833005951166339e-8);
    const vector4double  VN5      = vec_splats(-1.4213390571557850962e-6);
    const vector4double  VN4      = vec_splats(0.000041603292906656984871);
    const vector4double  VN3      = vec_splats(-0.00013134036773265025626);
    const vector4double  VN2      = vec_splats(0.038657983986041781264);
    const vector4double  VN1      = vec_splats(0.11285044772717598220);
    const vector4double  VN0      = vec_splats(1.1283802385263030286);
    
    const vector4double  VD3      = vec_splats(0.0066752224023576045451);
    const vector4double  VD2      = vec_splats(0.078647795836373922256);
    const vector4double  VD1      = vec_splats(0.43336185284710920150);
    const vector4double  VD0      = vec_splats(1.0);
    
    vector4double z4;
    vector4double polyVN0,polyVN1,polyVD0,polyVD1;
    
    z4             = vec_mul(z2,z2);
    
    polyVD1        = vec_madd(VD3,z4,VD1);
    polyVD0        = vec_madd(VD2,z4,VD0);
    polyVD0        = vec_madd(polyVD1,z2,polyVD0);
    
    polyVD0        = gmx_ibm_qpx_inv(polyVD0);
    
    polyVN0        = vec_madd(VN6,z4,VN4);
    polyVN1        = vec_madd(VN5,z4,VN3);    
    polyVN0        = vec_madd(polyVN0,z4,VN2);
    polyVN1        = vec_madd(polyVN1,z4,VN1);
    polyVN0        = vec_madd(polyVN0,z4,VN0);
    polyVN0        = vec_madd(polyVN1,z2,polyVN0);
#endif
    return   vec_mul(polyVN0,polyVD0);
}




#endif /* _gmx_math_ibm_qpx_single_h_ */
