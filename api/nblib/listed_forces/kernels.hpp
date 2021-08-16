/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*! \inpublicapi \file
 * \brief
 * Implements kernels for nblib supported bondtypes
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_LISTEDFORCES_KERNELS_HPP
#define NBLIB_LISTEDFORCES_KERNELS_HPP

#include <tuple>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "nblib/listed_forces/bondtypes.h"

#define NBLIB_ALWAYS_INLINE __attribute((always_inline))

namespace nblib
{

/*! \brief kernel to calculate the scalar part for simple harmonic bond forces
 *         for lambda = 0
 *
 * \param k spring constant
 * \param x0 equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template <class T>
inline std::tuple<T, T> harmonicScalarForce(T k, T x0, T x)
{
    T dx  = x - x0;
    T dx2 = dx * dx;

    T force = -k * dx;
    T epot = 0.5 * k * dx2;

    return std::make_tuple(force, epot);

    /* That was 6 flops */
}

/*! \brief kernel to calculate the scalar part for simple harmonic bond forces
 *         for non-zero lambda to interpolate between A and B states
 *
 * \param kA spring constant state A
 * \param kB spring constant state B
 * \param xA equilibrium distance state A
 * \param xB equilibrium distance state B
 * \param x  input bond length
 * \param lambda interpolation factor between A and B state
 *
 * \return tuple<force, potential energy, lambda-interpolated energy>
 */
template <class T>
inline std::tuple<T, T, T> harmonicScalarForce(T kA, T kB, T xA, T xB, T x, T lambda)
{
    // code unchanged relative to Gromacs

    T L1 = 1.0 - lambda;
    T kk = L1 * kA + lambda * kB;
    T x0 = L1 * xA + lambda * xB;

    T dx  = x - x0;
    T dx2 = dx * dx;

    T force     = -kk * dx;
    T epot      = 0.5 * kk * dx2;
    T dvdlambda = 0.5 * (kB - kA) * dx2 + (xA - xB) * kk * dx;

    return std::make_tuple(force, epot, dvdlambda);

    /* That was 19 flops */
}

//! abstraction layer for different 2-center bonds
template <class T>
inline auto bondKernel(T dr, const HarmonicBondType& bond)
{
    return harmonicScalarForce(bond.forceConstant(), bond.equilConstant(), dr);
}


/*! \brief kernel to calculate the scalar part for the forth power pontential bond forces
 *         for lambda = 0
 *
 * \param k spring constant
 * \param x0 squared equilibrium distance
 * \param x  squared input bond length
 *
 * \return tuple<force, potential energy>
 */
template <class T>
inline std::tuple<T, T> g96ScalarForce(T k, T x0, T x)
{
    T dx  = x - x0;
    T dx2 = dx * dx;

    T force = -k * dx;
    T epot = 0.5 * k * dx2;

    return std::make_tuple(force, epot);

    /* That was 6 flops */
}

/*! \brief kernel to calculate the scalar part for forth power pontential bond forces
 *         for non-zero lambda to interpolate between A and B states
 *
 * \param kA spring constant state A
 * \param kB spring constant state B
 * \param xA squared equilibrium distance state A
 * \param xB squared equilibrium distance state B
 * \param x  squared input bond length
 * \param lambda interpolation factor between A and B state
 *
 * \return tuple<force, potential energy, lambda-interpolated energy>
 */
template <class T>
inline std::tuple<T, T, T> g96ScalarForce(T kA, T kB, T xA, T xB, T x, T lambda)
{
    T L1 = 1.0 - lambda;
    T kk = L1 * kA + lambda * kB;
    T x0 = L1 * xA + lambda * xB;

    T dx  = x - x0;
    T dx2 = dx * dx;

    T force = -kk * dx;
    T epot = 0.5 * kk * dx2;
    // TODO: Check if this is 1/2 or 1/4
    T dvdlambda = 0.5 * (kB - kA) * dx2 + (xA - xB) * kk * dx;

    return std::make_tuple(force, epot, dvdlambda);

    /* That was 21 flops */
}

//! Abstraction layer for different 2-center bonds. Fourth power case
template <class T>
inline auto bondKernel(T dr, const G96BondType& bond)
{
    auto [force, ePot] = g96ScalarForce(bond.forceConstant(), bond.equilConstant(), dr*dr);
    force *= dr;
    ePot  *= 0.5;
    return std::make_tuple(force, ePot);
}


/*! \brief kernel to calculate the scalar part for the morse pontential bond forces
 *         for lambda = 0
 *
 * \param k force constant
 * \param beta beta exponent
 * \param x0 equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template <class T>
inline std::tuple<T, T> morseScalarForce(T k, T beta, T x0, T x)
{
    T exponent = std::exp(-beta * (x - x0));      /* 12 */
    T omexp = 1.0 - exponent;                     /*  1 */
    T kexp = k * omexp;                           /*  1 */

    T epot = kexp * omexp;                        /*  1 */
    T force = -2.0 * beta * exponent * kexp;      /*  4 */

    return std::make_tuple(force, epot);

    /* That was 20 flops */
}

/*! \brief kernel to calculate the scalar part for morse potential bond forces
 *         for non-zero lambda to interpolate between A and B states
 *
 * \param kA force constant state A
 * \param kB force constant state B
 * \param betaA beta exponent state A
 * \param betaB beta exponent state B
 * \param xA equilibrium distance state A
 * \param xB equilibrium distance state B
 * \param x  input bond length
 * \param lambda interpolation factor between A and B state
 *
 * \return tuple<force, potential energy, lambda-interpolated energy>
 */
template <class T>
inline std::tuple<T, T, T> morseScalarForce(T kA, T kB, T betaA, T betaB, T xA, T xB, T x, T lambda)
{
    T L1 = 1.0 - lambda;                          /* 1 */
    T x0 = L1 * xA + lambda * xB;                 /* 3 */
    T beta = L1 * betaA + lambda * betaB;         /* 3 */
    T k = L1 * kA + lambda * kB;                  /* 3 */

    T exponent = std::exp(-beta * (x - x0));      /* 12 */
    T omexp = 1.0 - exponent;                     /*  1 */
    T kexp = k * omexp;                           /*  1 */

    T epot = kexp * omexp;                        /*  1 */
    T force = -2.0 * beta * exponent * kexp;      /*  4 */

    T dvdlambda = (kB - kA) * omexp * omexp
                    - (2.0 - 2.0 * omexp) * omexp * k
                    * ((xB - xA) * beta - (betaB - betaA) * (x - x0)); /* 15 */

    return std::make_tuple(force, epot, dvdlambda);

    /* That was 44 flops */
}

//! Abstraction layer for different 2-center bonds. Morse case
template <class T>
inline auto bondKernel(T dr, const MorseBondType& bond)
{
    return morseScalarForce(bond.forceConstant(), bond.exponent(), bond.equilDistance(), dr);
}


/*! \brief kernel to calculate the scalar part for the FENE pontential bond forces
 *         for lambda = 0
 *
 * \param k spring constant
 * \param x0 equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template <class T>
inline std::tuple<T, T> FENEScalarForce(T k, T x0, T x)
{
    T x02 = x0 * x0;
    T x2 = x * x;

    T omx2_ox02 = 1.0 - (x2 / x02);

    T epot = -0.5 * k * x02 * std::log(omx2_ox02);
    T force = -k / omx2_ox02;

    return std::make_tuple(force, epot);

    /* That was 24 flops */
}

// TODO: Implement the free energy version of FENE (finitely extensible nonlinear elastic) bond types

//! Abstraction layer for different 2-center bonds. FENE case
template <class T>
inline auto bondKernel(T dr, const FENEBondType& bond)
{
    auto [force, ePot] = FENEScalarForce(bond.forceConstant(), bond.equilConstant(), dr);
    force *= dr;
    return std::make_tuple(force, ePot);
}


/*! \brief kernel to calculate the scalar part for cubic potential bond forces
 *         for lambda = 0
 *
 * \param kc cubic spring constant
 * \param kq quadratic spring constant
 * \param x0 equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template <class T>
inline std::tuple<T, T> cubicScalarForce(T kc, T kq, T x0, T x)
{
    T dx = x - x0;
    //T dx2 = dx * dx;

    T kdist  = kq * dx;
    T kdist2 = kdist * dx;

    T epot = kdist2 + (kc * kdist2 * dx);
    T force = -((2.0 * kdist) + (3.0 * kdist2 * kc));

    return std::make_tuple(force, epot);

    /* That was 16 flops */
}

// TODO: Implement the free energy version of Cubic bond types

template <class T>
inline auto bondKernel(T dr, const CubicBondType& bond)
{
    return cubicScalarForce(bond.cubicForceConstant(), bond.quadraticForceConstant(), bond.equilDistance(), dr);
}


/*! \brief kernel to calculate the scalar part for half attractive potential bond forces
 *         for lambda = 0
 *
 * \param k spring constant
 * \param x0 equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template <class T>
inline std::tuple<T, T> halfAttractiveScalarForce(T k, T x0, T x)
{
    T dx = x - x0;
    T dx2 = dx * dx;
    T dx3 = dx2 * dx;
    T dx4 = dx2 * dx2;

    T epot = -0.5 * k * dx4;
    T force = -2.0 * k * dx3;

    return std::make_tuple(force, epot);

    /* That was 10 flops */
}

/*! \brief kernel to calculate the scalar part for half attractive potential bond forces
 *         for non-zero lambda to interpolate between A and B states
 *
 * \param kA spring constant state A
 * \param kB spring constant state B
 * \param xA equilibrium distance state A
 * \param xB equilibrium distance state B
 * \param x  input bond length
 * \param lambda interpolation factor between A and B state
 *
 * \return tuple<force, potential energy, lambda-interpolated energy>
 */
template <class T>
inline std::tuple<T, T, T> halfAttractiveScalarForce(T kA, T kB, T xA, T xB, T x, T lambda)
{
    T L1 = 1.0 - lambda;
    T kk = L1 * kA + lambda * kB;
    T x0 = L1 * xA + lambda * xB;

    T dx  = x - x0;
    T dx2 = dx * dx;
    T dx3 = dx2 * dx;
    T dx4 = dx2 * dx2;

    T epot = -0.5 * kk * dx4;
    T force = -2.0 * kk * dx3;
    T dvdlambda = 0.5 * (kB - kA) * dx4 + (2.0 * (xA - xB) * kk * dx3);

    return std::make_tuple(force, epot, dvdlambda);

    /* That was 29 flops */
}

template <class T>
inline auto bondKernel(T dr, const HalfAttractiveQuarticBondType& bond)
{
    return halfAttractiveScalarForce(bond.forceConstant(), bond.equilConstant(), dr);
}


//! Three-center interaction type kernels


//! Harmonic Angle
template <class T>
inline auto threeCenterKernel(T dr, const HarmonicAngle& angle)
{
    return harmonicScalarForce(angle.forceConstant(), angle.equilConstant(), dr);
}


//! \brief Computes and returns the proper dihedral force
template <class T>
inline auto fourCenterKernel(T phi, const ProperDihedral& properDihedral)
{
    T deltaPhi = properDihedral.multiplicity() * phi - properDihedral.equilDistance();
    T force = -properDihedral.forceConstant() * properDihedral.multiplicity() * std::sin(deltaPhi);
    T ePot = properDihedral.forceConstant() * ( 1 + std::cos(deltaPhi) );
    return std::make_tuple(force, ePot);
}


//! \brief Ensure that a geometric quantity lies in (-pi, pi)
template<class T>
inline void makeAnglePeriodic(T& angle)
{
    if (angle >= M_PI)
    {
        angle -= 2 * M_PI;
    }
    else if (angle < -M_PI)
    {
        angle += 2 * M_PI;
    }
}

/*! \brief calculate the cosine of the angle between aInput and bInput
 *
 * \tparam T       float or double
 * \param aInput   aInput 3D vector
 * \param bInput   another 3D vector
 * \return         the cosine of the angle between aInput and bInput
 *
 *                       ax*bx + ay*by + az*bz
 * cos(aInput,bInput) = -----------------------, where aInput = (ax, ay, az)
 *                      ||aInput|| * ||bInput||
 */
template<class T>
inline T basicVectorCosAngle(gmx::BasicVector<T> aInput, gmx::BasicVector<T> bInput)
{
    gmx::BasicVector<double> a_double(aInput[0], aInput[1], aInput[2]);
    gmx::BasicVector<double> b_double(bInput[0], bInput[1], bInput[2]);

    double numerator     = dot(a_double, b_double);
    double denominatorSq = dot(a_double, a_double) * dot(b_double, b_double);

    T cosval = (denominatorSq > 0) ? static_cast<T>(numerator * gmx::invsqrt(denominatorSq)) : 1;
    cosval   = std::min(cosval, T(1.0));

    /* 25 TOTAL */
    return std::max(cosval, T(-1.0));
}

/*! \brief compute the angle between vectors a and b
 *
 * \tparam T    scalar type (float, double, or similar)
 * \param a     a 3D vector
 * \param b     another 3D vector
 * \return      the angle between a and b
 *
 * This routine calculates the angle between a & b without any loss of accuracy close to 0/PI.
 *
 * Note: This function is not (yet) implemented for the C++ replacement of the
 * deprecated rvec and dvec.
 */
template<class T>
inline T basicVectorAngle(gmx::BasicVector<T> a, gmx::BasicVector<T> b)
{
    gmx::BasicVector<T> w = cross(a, b);

    T wlen = norm(w);
    T s    = dot(a, b);

    return std::atan2(wlen, s);
}

/*! \brief Computes the dihedral phi angle and two cross products
 *
 * \tparam T        scalar type (float or double, or similar)
 * \param[in] dxIJ
 * \param[in] dxKJ
 * \param[in] dxKL
 * \param[out] m    output for \p dxIJ x \p dxKJ
 * \param[out] n    output for \p dxKJ x \p dxKL
 * \return          the angle between m and n
 */
template<class T>
inline T dihedralPhi(gmx::BasicVector<T> dxIJ,
                     gmx::BasicVector<T> dxKJ,
                     gmx::BasicVector<T> dxKL,
                     gmx::BasicVector<T>* m,
                     gmx::BasicVector<T>* n)
{
    *m = cross(dxIJ, dxKJ);
    *n = cross(dxKJ, dxKL);
    T phi  = basicVectorAngle(*m, *n);
    T ipr  = dot(dxIJ, *n);
    T sign = (ipr < 0.0) ? -1.0 : 1.0;
    phi    = sign * phi;
    return phi;
}

//! \brief Computes and returns the improper dihedral force
template <class T>
inline auto fourCenterKernel(T phi, const ImproperDihedral& improperDihedral)
{
    T deltaPhi = phi - improperDihedral.equilConstant();
    /* deltaPhi cannot be outside (-pi,pi) */
    makeAnglePeriodic(deltaPhi);
    const T force = -improperDihedral.forceConstant()  * deltaPhi;
    const T ePot = 0.5 * improperDihedral.forceConstant() * deltaPhi * deltaPhi;
    return std::make_tuple(force, ePot);
}

//! \brief Computes and returns the Ryckaert-Belleman dihedral force
template <class T>
inline auto fourCenterKernel(T phi, const RyckaertBellemanDihedral& ryckaertBellemanDihedral)
{
    /* Change to polymer convention */
    const T localPhi = (phi < 0) ? (phi += M_PI) : (phi -= M_PI);
    T cos_phi = std::cos(localPhi);
    T ePot = ryckaertBellemanDihedral[0];
    T force = 0;
    T cosineFactor = 1;

    for (int i = 1; i < int(ryckaertBellemanDihedral.size()); i++)
    {
        force += ryckaertBellemanDihedral[i] * cosineFactor * i;
        cosineFactor *= cos_phi;
        ePot += cosineFactor * ryckaertBellemanDihedral[i];
    }
    /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
    force = -force * std::sin(localPhi);
    return std::make_tuple(force, ePot);
}

//! Two-center category common

//! \brief add shift forces, if requested (static compiler decision)
template<class T, class ShiftForce>
inline void
addShiftForce(const gmx::BasicVector<T>& interactionForce, ShiftForce* shiftForce)
{
    *shiftForce += interactionForce;
}

//! \brief this will be called if shift forces are not computed (and removed by the compiler)
template<class T>
inline void addShiftForce([[maybe_unused]] const gmx::BasicVector<T>& fij,
                          [[maybe_unused]] std::nullptr_t*)
{
}

/*! \brief Spreads and accumulates the bonded forces to the two atoms and adds the virial contribution when needed
 *
 * \p shiftIndex is used as the periodic shift.
 */
template <class T, class ShiftForce>
inline void spreadTwoCenterForces(const T bondForce,
                                  const gmx::BasicVector<T>& dx,
                                  gmx::BasicVector<T>*       force_i,
                                  gmx::BasicVector<T>*       force_j,
                                  ShiftForce*                shf_ij,
                                  ShiftForce*                shf_c)
{
    gmx::BasicVector<T> fij = bondForce * dx;
    *force_i += fij;
    *force_j -= fij;

    addShiftForce(fij, shf_ij);
    addShiftForce(T(-1.0)*fij, shf_c);
    /* 15 Total */
}

//! Three-center category common

/*! \brief spread force to 3 centers based on scalar force and angle
 *
 * \tparam T
 * \param cos_theta
 * \param force
 * \param r_ij
 * \param r_kj
 * \param force_i
 * \param force_j
 * \param force_k
 */
template <class T, class ShiftForce>
inline void spreadThreeCenterForces(T                          cos_theta,
                                    T                          force,
                                    const gmx::BasicVector<T>& r_ij,
                                    const gmx::BasicVector<T>& r_kj,
                                    gmx::BasicVector<T>*       force_i,
                                    gmx::BasicVector<T>*       force_j,
                                    gmx::BasicVector<T>*       force_k,
                                    ShiftForce*                shf_ij,
                                    ShiftForce*                shf_kj,
                                    ShiftForce*                shf_c)
{
    T cos_theta2 = cos_theta * cos_theta;
    if (cos_theta2 < 1)                              /*   1		*/
    {
        T st    = force / std::sqrt(1 - cos_theta2); /*  12		*/
        T sth   = st * cos_theta;                    /*   1		*/
        T nrij2 = dot(r_ij, r_ij);                   /*   5		*/
        T nrkj2 = dot(r_kj, r_kj);                   /*   5		*/

        T cik = st / std::sqrt(nrij2 * nrkj2);       /*  11		*/
        T cii = sth / nrij2;                         /*   1		*/
        T ckk = sth / nrkj2;                         /*   1		*/

        /*  33 */
        gmx::BasicVector<T> f_i = cii * r_ij - cik * r_kj;
        gmx::BasicVector<T> f_k = ckk * r_kj - cik * r_ij;
        gmx::BasicVector<T> f_j = T(-1.0) * (f_i + f_k);
        *force_i += f_i;
        *force_j += f_j;
        *force_k += f_k;

        addShiftForce(f_i, shf_ij);
        addShiftForce(f_j, shf_c);
        addShiftForce(f_k, shf_kj);

    } /* 70 TOTAL	*/
}

//! Four-center category common
template <class T, class ShiftForce>
inline void spreadFourCenterForces(T force,
                                   const gmx::BasicVector<T>& dxIJ,
                                   const gmx::BasicVector<T>& dxJK,
                                   const gmx::BasicVector<T>& dxKL,
                                   const gmx::BasicVector<T>& m,
                                   const gmx::BasicVector<T>& n,
                                   gmx::BasicVector<T>* force_i,
                                   gmx::BasicVector<T>* force_j,
                                   gmx::BasicVector<T>* force_k,
                                   gmx::BasicVector<T>* force_l,
                                   ShiftForce* shf_ij,
                                   ShiftForce* shf_kj,
                                   ShiftForce* shf_lj,
                                   ShiftForce* shf_c)
{
    T norm2_m = dot(m, m);                             /* 5 */
    T norm2_n = dot(n, n);                             /* 5 */
    T norm2_jk = dot(dxJK, dxJK);                      /* 5 */
    T toler = norm2_jk * GMX_REAL_EPS;
    if ((norm2_m > toler) && (norm2_n > toler))
    {
        T rcp_norm2_jk = 1.0f / norm2_jk;              /* 1 */
        T norm_jk = std::sqrt(norm2_jk);               /* 10 */

        T a = -force * norm_jk / norm2_m;              /* 11 */
        gmx::BasicVector<T> f_i = a * m;               /* 3 */

        T b = force * norm_jk / norm2_n;               /* 11 */
        gmx::BasicVector<T> f_l = b * n;               /* 3 */

        T p = rcp_norm2_jk * dot(dxIJ, dxJK);          /* 6 */
        T q = rcp_norm2_jk * dot(dxKL, dxJK);          /* 6 */
        gmx::BasicVector<T> svec = p * f_i - q * f_l;  /* 9 */

        gmx::BasicVector<T> f_j = svec - f_i;          /* 3 */
        gmx::BasicVector<T> f_k = T(-1.0)*svec - f_l;  /* 6 */

        *force_i += f_i;                               /* 3 */
        *force_j += f_j;                               /* 3 */
        *force_k += f_k;                               /* 3 */
        *force_l += f_l;                               /* 3 */

        addShiftForce(f_i, shf_ij);
        addShiftForce(f_j, shf_c);
        addShiftForce(f_k, shf_kj);
        addShiftForce(f_l, shf_lj);
    }
}

} // namespace nblib

#undef NBLIB_ALWAYS_INLINE

#endif // NBLIB_LISTEDFORCES_KERNELS_HPP
