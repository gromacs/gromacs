/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#include <limits>

#include "nblib/listed_forces/bondtypes.h"
#include "nblib/listed_forces/selectors.hpp"
#include "nblib/util/array.hpp"
#include "nblib/util/functions.hpp"
#include "nblib/util/traits.hpp"
#include "nblib/util/tuple.hpp"

namespace nblib
{

template<typename T>
HOST_DEVICE_INLINE util::tuple<T, T> oneCenterKernel(T dr, const PositionRestraints& restraint, int dim)
{
    T force = -restraint.forceConstant(dim) * dr;
    T epot  = T(0.5) * restraint.forceConstant(dim) * dr * dr;

    return util::make_tuple(force, epot);
}

/*! \brief kernel to calculate the scalar part for simple harmonic bond forces
 *         for lambda = 0
 *
 * \param kA spring constant
 * \param xA equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy, (empty) dvdl>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
harmonicScalarForce(T kA, T /*kB*/, T xA, T /*xB*/, T x, NoFepLambdaType /*lambda*/)
{
    T dx  = x - xA;
    T dx2 = dx * dx;

    T force = -kA * dx;
    T epot  = T(0.5) * kA * dx2;

    return util::make_tuple(force, epot, T(0));
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
 * \return tuple<force, potential energy, dvdl>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T> harmonicScalarForce(T kA, T kB, T xA, T xB, T x, T lambda)
{
    T L1 = T(1.0) - lambda;
    T kk = L1 * kA + lambda * kB;
    T x0 = L1 * xA + lambda * xB;

    T dx  = x - x0;
    T dx2 = dx * dx;

    T force     = -kk * dx;
    T epot      = T(0.5) * kk * dx2;
    T dvdlambda = T(0.5) * (kB - kA) * dx2 + (xA - xB) * kk * dx;

    return util::make_tuple(force, epot, dvdlambda);
}

//! abstraction layer for different 2-center bonds
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
bondKernel(T dr, const HarmonicBondType& bondA, const HarmonicBondType& bondB, L lambda)
{
    return harmonicScalarForce(
            bondA.forceConstant(), bondB.forceConstant(), bondA.equilConstant(), bondB.equilConstant(), dr, lambda);
}


/*! \brief kernel to calculate the scalar part for the forth power pontential bond forces
 *         for lambda = 0
 *
 * \param kA spring constant
 * \param xA squared equilibrium distance
 * \param x  squared input bond length
 *
 * \return tuple<force, potential energy>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
g96ScalarForce(T kA, T /*kB*/, T xA, T /*xB*/, T x, NoFepLambdaType /*lambda*/)
{
    T dx  = x - xA;
    T dx2 = dx * dx;

    T force = -kA * dx;
    T epot  = T(0.5) * kA * dx2;

    return util::make_tuple(force, epot, T(0));
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
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T> g96ScalarForce(T kA, T kB, T xA, T xB, T x, T lambda)
{
    T L1 = T(1.0) - lambda;
    T kk = L1 * kA + lambda * kB;
    T x0 = L1 * xA + lambda * xB;

    T dx  = x - x0;
    T dx2 = dx * dx;

    T force = -kk * dx;
    T epot  = T(0.5) * kk * dx2;
    // TODO: Check if this is 1/2 or 1/4
    T dvdlambda = T(0.5) * (kB - kA) * dx2 + (xA - xB) * kk * dx;

    return util::make_tuple(force, epot, dvdlambda);
}

//! Abstraction layer for different 2-center bonds. G96 bonds.
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
bondKernel(T dr, const G96BondType& bondA, const G96BondType& bondB, L lambda)
{
    auto [force, ePot, dvdl] = g96ScalarForce(bondA.forceConstant(),
                                              bondB.forceConstant(),
                                              bondA.equilConstant(),
                                              bondB.forceConstant(),
                                              dr * dr,
                                              lambda);
    force *= dr;
    ePot *= T(0.5);
    return util::make_tuple(force, ePot, dvdl);
}


/*! \brief kernel to calculate the scalar part for the morse pontential bond forces
 *         for lambda = 0
 *
 * \param kA force constant
 * \param betaA beta exponent
 * \param xA equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   morseScalarForce(T kA, T /*kB*/, T betaA, T /*betaB*/, T xA, T /*xB*/, T x, NoFepLambdaType /*lambda*/)
{
    T exponent = std::exp(-betaA * (x - xA));
    T omexp    = T(1.0) - exponent;
    T kexp     = kA * omexp;

    T epot  = kexp * omexp;
    T force = T(-2.0) * betaA * exponent * kexp;

    return util::make_tuple(force, epot, T(0));
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
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
morseScalarForce(T kA, T kB, T betaA, T betaB, T xA, T xB, T x, T lambda)
{
    T L1   = T(1.0) - lambda;
    T x0   = L1 * xA + lambda * xB;
    T beta = L1 * betaA + lambda * betaB;
    T k    = L1 * kA + lambda * kB;

    T exponent = std::exp(-beta * (x - x0));
    T omexp    = T(1.0) - exponent;
    T kexp     = k * omexp;

    T epot  = kexp * omexp;
    T force = T(-2.0) * beta * exponent * kexp;

    T dvdlambda =
            (kB - kA) * omexp * omexp
            - (T(2.0) - T(2.0) * omexp) * omexp * k * ((xB - xA) * beta - (betaB - betaA) * (x - x0));

    return util::make_tuple(force, epot, dvdlambda);
}

//! Abstraction layer for different 2-center bonds. Morse case
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
bondKernel(T dr, const MorseBondType& bondA, const MorseBondType& bondB, L lambda)
{
    return morseScalarForce(bondA.forceConstant(),
                            bondB.forceConstant(),
                            bondA.exponent(),
                            bondB.exponent(),
                            bondA.equilDistance(),
                            bondB.equilDistance(),
                            dr,
                            lambda);
}


/*! \brief kernel to calculate the scalar part for the FENE pontential bond forces
 *         for lambda = 0
 *
 * \param kA spring constant
 * \param xA equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
FENEScalarForce(T kA, T /*kB*/, T xA, T /*xB*/, T x, NoFepLambdaType /*lambda*/)
{
    T x02 = xA * xA;
    T x2  = x * x;

    T omx2_ox02 = T(1.0) - (x2 / x02);

    T epot  = T(-0.5) * kA * x02 * std::log(omx2_ox02);
    T force = -kA / omx2_ox02;

    return util::make_tuple(force, epot, T(0));
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
        FENEScalarForce(T /*kA*/, T /*kB*/, T /*xA*/, T /*xB*/, T /*x*/, T /*lambda*/)
{
    printf("FENE FEP not implemented\n");
    return util::make_tuple(T(0), T(0), T(0));
}

//! Abstraction layer for different 2-center bonds. FENE case
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
bondKernel(T dr, const FENEBondType& bondA, const FENEBondType& bondB, L lambda)
{
    auto [force, ePot, dvdl] = FENEScalarForce(
            bondA.forceConstant(), bondB.forceConstant(), bondA.equilConstant(), bondB.equilConstant(), dr, lambda);
    force *= dr;
    return util::make_tuple(force, ePot, dvdl);
}


/*! \brief kernel to calculate the scalar part for cubic potential bond forces
 *         for lambda = 0
 *
 * \param kcA cubic spring constant
 * \param kqA quadratic spring constant
 * \param xA equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   cubicScalarForce(T kcA, T /*kcB*/, T kqA, T /*kqB*/, T xA, T /*xB*/, T x, NoFepLambdaType /*lambda*/)
{
    T dx = x - xA;

    T kdist  = kqA * dx;
    T kdist2 = kdist * dx;

    T epot  = kdist2 + (kcA * kdist2 * dx);
    T force = -((T(2.0) * kdist) + (T(3.0) * kdist2 * kcA));

    return util::make_tuple(force, epot, T(0));
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   cubicScalarForce(T /*kcA*/, T /*kcB*/, T /*kqA*/, T /*kqB*/, T /*xA*/, T /*xB*/, T /*x*/, T /*lambda*/)
{
    printf("Cubic FEP not implemented\n");
    return util::make_tuple(T(0), T(0), T(0));
}

template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
bondKernel(T dr, const CubicBondType& bondA, const CubicBondType& bondB, L lambda)
{
    return cubicScalarForce(bondA.cubicForceConstant(),
                            bondB.cubicForceConstant(),
                            bondA.quadraticForceConstant(),
                            bondB.quadraticForceConstant(),
                            bondA.equilDistance(),
                            bondB.equilDistance(),
                            dr,
                            lambda);
}


/*! \brief kernel to calculate the scalar part for half attractive potential bond forces
 *         for lambda = 0
 *
 * \param kA spring constant
 * \param xA equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
halfAttractiveScalarForce(T kA, T /*kb*/, T xA, T /*xB*/, T x, NoFepLambdaType /*lambda*/)
{
    T dx  = x - xA;
    T dx2 = dx * dx;
    T dx3 = dx2 * dx;
    T dx4 = dx2 * dx2;

    T epot  = T(-0.5) * kA * dx4;
    T force = T(-2.0) * kA * dx3;

    return util::make_tuple(force, epot, T(0));
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
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   halfAttractiveScalarForce(T kA, T kB, T xA, T xB, T x, T lambda)
{
    T L1 = T(1) - lambda;
    T kk = L1 * kA + lambda * kB;
    T x0 = L1 * xA + lambda * xB;

    T dx  = x - x0;
    T dx2 = dx * dx;
    T dx3 = dx2 * dx;
    T dx4 = dx2 * dx2;

    T epot      = T(-0.5) * kk * dx4;
    T force     = T(-2.0) * kk * dx3;
    T dvdlambda = T(0.5) * (kB - kA) * dx4 + (T(2.0) * (xA - xB) * kk * dx3);

    return util::make_tuple(force, epot, dvdlambda);
}

template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto bondKernel(T                                    dr,
                                                   const HalfAttractiveQuarticBondType& bondA,
                                                   const HalfAttractiveQuarticBondType& bondB,
                                                   L                                    lambda)
{
    return halfAttractiveScalarForce(
            bondA.forceConstant(), bondB.forceConstant(), bondA.equilConstant(), bondB.equilConstant(), dr, lambda);
}

/*! \brief kernel to calculate the scalar part for the 1-4 LJ non-bonded forces
 *
 * \param c6A  C6 parameter of LJ potential
 * \param c12A C12 parameter of LJ potential
 * \param qiqj product of the charges on the two particles
 * \param r    distance between the atoms
 *
 * \return tuple<force, ePot_VdW, ePot_Coulomb, dvdlambda>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T, T>
                                   pairLJScalarForce(C6 c6A, C6 /*c6B*/, C12 c12A, C12 /*c12B*/, T qiqj, T r, NoFepLambdaType /*lambda*/)
{
    T rinv  = T(1) / r;
    T rinv2 = rinv * rinv;
    T rinv6 = rinv2 * rinv2 * rinv2;

    T eVdw = rinv6 * (c12A * rinv6 - c6A);

    T c6_  = T(6) * c6A;
    T c12_ = T(12) * c12A;

    // T eCoul = ONE_4PI_EPS0 * qiqj * rinv;
    T eCoul = qiqj * rinv;

    T force = (rinv6 * (c12_ * rinv6 - c6_) + eCoul) * rinv;

    return util::make_tuple(force, eVdw, eCoul, T(0));
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T, T>
                                   pairLJScalarForce(C6 /*c6A*/, C6 /*c6B*/, C12 /*c12A*/, C12 /*c12B*/, T /*qiqj*/, T /*r*/, T /*lambda*/)
{
    printf("PairLJ FEP not implemented\n");
    return util::make_tuple(T(0), T(0), T(0), T(0));
}

//! Abstraction layer for different 2-center bonds. 1-4 LJ pair interactions case
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
bondKernel(T dr, T qiqj, const PairLJType& bondA, const PairLJType& bondB, L lambda)
{
    return pairLJScalarForce(bondA.c6(), bondB.c6(), bondA.c12(), bondB.c12(), qiqj, dr, lambda);
}

//! Abstraction layer for different 2-center bonds. Charged 1-4 LJ pair interactions case
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
bondKernel(T dr, T /*qiqj*/, const PairLJChargeType& bondA, const PairLJChargeType& bondB, L lambda)
{
    return pairLJScalarForce(
            bondA.c6(), bondB.c6(), bondA.c12(), bondB.c12(), bondA.qi() * bondA.qj() * bondA.ff(), dr, lambda);
}

template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
bondKernel(T dr, T qj, const SimplePolarization& polA, const SimplePolarization& polB, L lambda)
{
    T qj2  = qj * qj;
    T kshA = qj2 * polA.coefficient();
    T kshB = qj2 * polB.coefficient();

    return harmonicScalarForce(kshA, kshB, T(0.0), T(0.0), dr, lambda);
}
//! Three-center interaction type kernels

/*! \brief kernel to calculate the scalar part for linear angle forces
 *         for lambda = 0
 *
 * \param kA force constant
 * \param aA equilibrium angle
 * \param angle current angle vaule
 *
 * \return tuple<force, potential energy>
 */
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T, T>
linearAnglesScalarForce(T kA, T /*kb*/, T aA, T /*aB*/, T angle, NoFepLambdaType /*lambda*/)
{
    T b = T(1.0) - aA;

    T kdr  = kA * angle;
    T epot = T(0.5) * kdr * angle;

    T ci = aA * kA;
    T ck = b * kA;

    return util::make_tuple(ci, ck, epot, T(0));
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T, T>
        linearAnglesScalarForce(T /*kA*/, T /*kb*/, T /*aA*/, T /*aB*/, T /*angle*/, T /*lambda*/)
{
    printf("Linear angles FEP not implemented\n");
    return util::make_tuple(T(0), T(0), T(0), T(0));
}

template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
threeCenterKernel(T dr, const LinearAngle& angleA, const LinearAngle& angleB, L lambda)
{
    return linearAnglesScalarForce(angleA.forceConstant(),
                                   angleB.forceConstant(),
                                   angleA.equilConstant(),
                                   angleB.equilConstant(),
                                   dr,
                                   lambda);
}

//! Harmonic Angle
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
threeCenterKernel(T dr, const HarmonicAngle& angleA, const HarmonicAngle& angleB, L lambda)
{
    return harmonicScalarForce(angleA.forceConstant(),
                               angleB.forceConstant(),
                               angleA.equilConstant(),
                               angleB.equilConstant(),
                               dr,
                               lambda);
}

//! Cosine based (GROMOS-96) Angle
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
threeCenterKernel(T dr, const G96Angle& angleA, const G96Angle& angleB, L lambda)
{
    auto costheta = std::cos(dr);
    auto feTuple  = g96ScalarForce(angleA.forceConstant(),
                                  angleB.forceConstant(),
                                  angleA.equilConstant(),
                                  angleB.equilConstant(),
                                  costheta,
                                  lambda);

    // The above kernel call effectively computes the derivative of the potential with respect to
    // cos(theta). However, we need the derivative with respect to theta. We use this extra
    // -sin(theta) factor to account for this before the forces are spread between the particles.

    util::get<0>(feTuple) *= -std::sqrt(1 - costheta * costheta);
    return feTuple;
}

template<class T, class L, class... Ts>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto threeCenterKernel(T                                  dr,
                                                          const ThreeCenterAggregate<Ts...>& aggregateA,
                                                          const ThreeCenterAggregate<Ts...>& aggregateB,
                                                          L                                  lambda)
{
    return threeCenterKernel(dr, aggregateA.carrier(), aggregateB.carrier(), lambda);
}

template<class T, class L, class... Ts>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto fourCenterKernel(T                                 phi,
                                                         const FourCenterAggregate<Ts...>& aggregateA,
                                                         const FourCenterAggregate<Ts...>& aggregateB,
                                                         L                                 lambda)
{
    return fourCenterKernel(phi, aggregateA.carrier(), aggregateB.carrier(), lambda);
}


/*! \brief kernel to calculate the scalar part for cross bond-bond forces
 *         for lambda = 0
 *
 * \param k force constant
 * \param r0ij equilibrium distance between particles i & j
 * \param r0kj equilibrium distance between particles k & j
 * \param rij  input bond length between particles i & j
 * \param rkj  input bond length between particles k & j
 *
 * \return tuple<force scalar i, force scalar k, potential energy>
 */

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T> crossBondBondScalarForce(T k, T r0ij, T r0kj, T rij, T rkj)
{
    T si = rij - r0ij;
    T sk = rkj - r0kj;

    T epot = k * si * sk;

    T ci = -k * sk / rij;
    T ck = -k * si / rkj;

    return util::make_tuple(ci, ck, epot);

    /* That was 8 flops */
}

// FEP for cross-bond-bond not implemented in GROMACS

//! Cross bond-bond interaction
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto threeCenterKernel(T drij, T drkj, const CrossBondBond& crossBondBond)
{
    return crossBondBondScalarForce(crossBondBond.forceConstant(),
                                    crossBondBond.equilDistanceIJ(),
                                    crossBondBond.equilDistanceKJ(),
                                    drij,
                                    drkj);
}

/*! \brief kernel to calculate the scalar part for cross bond-angle forces
 *         for lambda = 0
 *
 * \param k force constant
 * \param r0ij equilibrium distance between particles i & j
 * \param r0kj equilibrium distance between particles k & j
 * \param r0ik equilibrium distance between particles i & k
 * \param rij  input bond length between particles i & j
 * \param rkj  input bond length between particles k & j
 * \param rik  input bond length between particles i & k
 *
 * \return tuple<atom i force, atom j force, atom k force, potential energy>
 */

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T, T>
crossBondAngleScalarForce(T k, T r0ij, T r0kj, T r0ik, T rij, T rkj, T rik)
{
    T sij = rij - r0ij;
    T skj = rkj - r0kj;
    T sik = rik - r0ik;

    T epot = k * sik * (sij + skj);

    T fi = -k * sik / rij;
    T fj = -k * sik / rkj;
    T fk = -k * (sij + skj) / rik;

    return util::make_tuple(fi, fj, fk, epot);

    /* That was 13 flops */
}

// FEP for cross-bond-angle not implemented in GROMACS

//! Cross bond-bond interaction
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto threeCenterKernel(T drij, T drkj, T drik, const CrossBondAngle& crossBondAngle)
{
    return crossBondAngleScalarForce(crossBondAngle.forceConstant(),
                                     crossBondAngle.equilDistanceIJ(),
                                     crossBondAngle.equilDistanceKJ(),
                                     crossBondAngle.equilDistanceIK(),
                                     drij,
                                     drkj,
                                     drik);
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto quarticAngleScalarForce(T                   dr,
                                                                const QuarticAngle& angle,
                                                                const QuarticAngle& /*angleB*/,
                                                                NoFepLambdaType /*lambda*/)
{
    T dt = dr - angle.equilConstant();

    T force  = 0;
    T energy = angle.forceConstant(0);
    T dtp    = 1.0;
    for (auto j = 1; j <= 4; j++)
    {
        T c = angle.forceConstant(j);
        force -= j * c * dtp;
        dtp *= dt;
        energy += c * dtp;
    }

    return util::make_tuple(force, energy, T(0));
}

template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto quarticAngleScalarForce(T /*dr*/,
                                                                const QuarticAngle& /*angle*/,
                                                                const QuarticAngle& /*angleB*/,
                                                                L /*lambda*/)
{
    printf("QuarticAngle FEP not implemented\n");
    return util::make_tuple(T(0), T(0), T(0));
}

template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
threeCenterKernel(T dr, const QuarticAngle& angleA, const QuarticAngle& angleB, L lambda)
{
    return quarticAngleScalarForce(dr, angleA, angleB, lambda);
}

//! \brief Restricted Angle potential. Returns scalar force and energy
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
threeCenterKernel(T theta, const RestrictedAngle& angleA, const RestrictedAngle& angleB, L lambda)
{
    T costheta               = std::cos(theta);
    auto [force, ePot, dvdl] = harmonicScalarForce(angleA.forceConstant(),
                                                   angleB.forceConstant(),
                                                   angleA.equilConstant(),
                                                   angleB.equilConstant(),
                                                   costheta,
                                                   lambda);

    // The above kernel call effectively computes the derivative of the potential with respect to
    // cos(theta). However, we need the derivative with respect to theta.
    // This introduces the extra (cos(theta)*cos(eqAngle) - 1)/(sin(theta)^3 factor for the force
    // The call also computes the potential energy without the sin(theta)^-2 factor

    T sintheta2 = T(1) - costheta * costheta;
    T sintheta  = std::sqrt(sintheta2);
    force *= (costheta * angleA.equilConstant() - T(1)) / (sintheta2 * sintheta);
    ePot /= sintheta2;

    // dvdl probably needs to be adjusted too
    dvdl = 0;

    return util::make_tuple(force, ePot, dvdl);
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   properDihedralForce(T cpA, T /*cpB*/, T phiA, T /*phiB*/, int multiplicity, T phi, NoFepLambdaType /*lambda*/)
{
    T deltaPhi = multiplicity * phi - phiA;
    T force    = -cpA * multiplicity * std::sin(deltaPhi);
    T ePot     = cpA * (T(1) + std::cos(deltaPhi));

    return util::make_tuple(force, ePot, T(0));
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   properDihedralForce(T /*cpA*/, T /*cpB*/, T /*phiA*/, T /*phiB*/, int /*mult*/, T /*phi*/, T /*lambda*/)
{
    printf("ProperDihedral FEP not implemented\n");
    return util::make_tuple(T(0), T(0), T(0));
}

//! \brief Computes and returns the proper dihedral force
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
fourCenterKernel(T phi, const ProperDihedral& pDihA, const ProperDihedral& pDihB, L lambda)
{
    return properDihedralForce(pDihA.forceConstant(),
                               pDihB.forceConstant(),
                               pDihA.equilDistance(),
                               pDihB.equilDistance(),
                               pDihA.multiplicity(),
                               phi,
                               lambda);
}

//! \brief Computes and returns the improper proper dihedral force
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto fourCenterKernel(T                             phi,
                                                         const ImproperProperDihedral& pDihA,
                                                         const ImproperProperDihedral& pDihB,
                                                         L                             lambda)
{
    return properDihedralForce(pDihA.forceConstant(),
                               pDihB.forceConstant(),
                               pDihA.equilDistance(),
                               pDihB.equilDistance(),
                               pDihA.multiplicity(),
                               phi,
                               lambda);
}


//! \brief Ensure that a geometric quantity lies in (-pi, pi)
template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE void makeAnglePeriodic(T& angle)
{
    if (angle >= T(M_PI))
    {
        angle -= T(2) * T(M_PI);
    }
    else if (angle < -T(M_PI))
    {
        angle += T(2) * T(M_PI);
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
template<class BasicVector>
HOST_DEVICE_FUN HOST_DEVICE_INLINE VectorValueType_t<BasicVector> basicVectorCosAngle(BasicVector aInput,
                                                                                      BasicVector bInput)
{
    using ValueType         = VectorValueType_t<BasicVector>;
    using BasicVectorDouble = typename SwapArg<BasicVector, double>::type;

    // Use a double vector here for increased precision
    BasicVectorDouble a_double{ aInput[0], aInput[1], aInput[2] };
    BasicVectorDouble b_double{ bInput[0], bInput[1], bInput[2] };

    double numerator     = dot(a_double, b_double);
    double denominatorSq = dot(a_double, a_double) * dot(b_double, b_double);

    ValueType cosval =
            (denominatorSq > 0) ? static_cast<ValueType>(numerator * util::invsqrt(denominatorSq)) : 1;
    cosval = util::min(cosval, ValueType(1.0));

    return util::max(cosval, ValueType(-1.0));
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
template<class BasicVector>
HOST_DEVICE_FUN HOST_DEVICE_INLINE VectorValueType_t<BasicVector> basicVectorAngle(BasicVector a,
                                                                                   BasicVector b)
{
    using ValueType = VectorValueType_t<BasicVector>;
    BasicVector w   = cross(a, b);

    ValueType wlen = norm(w);
    ValueType s    = dot(a, b);

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
template<class BasicVector>
HOST_DEVICE_FUN HOST_DEVICE_INLINE VectorValueType_t<BasicVector>
dihedralPhi(BasicVector dxIJ, BasicVector dxKJ, BasicVector dxKL, BasicVector* m, BasicVector* n)
{
    using ValueType = VectorValueType_t<BasicVector>;
    *m              = cross(dxIJ, dxKJ);
    *n              = cross(dxKJ, dxKL);
    ValueType phi   = basicVectorAngle(*m, *n);
    ValueType ipr   = dot(dxIJ, *n);
    ValueType sign  = (ipr < ValueType(0)) ? -1.0 : 1.0;
    phi             = sign * phi;
    return phi;
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
improperDihedralForce(T cpA, T /*cpB*/, T phiA, T /*phiB*/, T phi, NoFepLambdaType /*lambda*/)
{
    T deltaPhi = phi - phiA;
    /* deltaPhi cannot be outside (-pi,pi) */
    makeAnglePeriodic(deltaPhi);
    T force = cpA * deltaPhi;
    T ePot  = T(0.5) * cpA * deltaPhi * deltaPhi;
    return util::make_tuple(force, ePot, T(0));
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   improperDihedralForce(T /*cpA*/, T /*cpB*/, T /*phiA*/, T /*phiB*/, int /*mult*/, T /*phi*/, T /*lambda*/)
{
    printf("ImproperDihedral FEP not implemented\n");
    return util::make_tuple(T(0), T(0), T(0));
}

//! \brief Computes and returns the improper dihedral force
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto
fourCenterKernel(T phi, const ImproperDihedral& iDihA, const ImproperDihedral& iDihB, L lambda)
{
    return improperDihedralForce(
            iDihA.forceConstant(), iDihB.forceConstant(), iDihA.equilConstant(), iDihB.equilConstant(), phi, lambda);
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   ryckaertBellemanForce(T       phi,
                                                         const T parmA[RyckaertBellemanDihedral::RbNumParameters],
                                                         const T /*parmB*/[RyckaertBellemanDihedral::RbNumParameters],
                                                         NoFepLambdaType /*lambda*/)
{
    /* Change to polymer convention */
    T localPhi = (phi < T(0)) ? (phi += M_PI) : (phi -= M_PI);

    T cos_phi      = std::cos(localPhi);
    T ePot         = parmA[0];
    T force        = 0;
    T cosineFactor = 1;

    for (int i = 1; i < RyckaertBellemanDihedral::RbNumParameters; i++)
    {
        force += parmA[i] * cosineFactor * i;
        cosineFactor *= cos_phi;
        ePot += cosineFactor * parmA[i];
    }
    /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
    force = -force * std::sin(localPhi);

    return util::make_tuple(force, ePot, T(0));
}

template<class T>
HOST_DEVICE_FUN HOST_DEVICE_INLINE util::tuple<T, T, T>
                                   ryckaertBellemanForce(T /*phi*/,
                                                         const T /*parmA*/[RyckaertBellemanDihedral::RbNumParameters],
                                                         const T /*parmB*/[RyckaertBellemanDihedral::RbNumParameters],
                                                         T /*lambda*/)
{
    printf("Ryckaert Belleman FEP not implemented\n");
    return util::make_tuple(T(0), T(0), T(0));
}

//! \brief Computes and returns the Ryckaert-Belleman dihedral force
template<class T, class L>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto fourCenterKernel(T                               phi,
                                                         const RyckaertBellemanDihedral& rbDihA,
                                                         const RyckaertBellemanDihedral& rbDihB,
                                                         L                               lambda)
{
    return ryckaertBellemanForce(phi, rbDihA.parameters(), rbDihB.parameters(), lambda);
}

/*! \brief Spreads and accumulates the forces between two atoms and adds the virial contribution when needed
 *
 * \tparam T         The type of vector, e.g. float, double, etc
 * \param force      The Force
 * \param dx         Distance between centers
 * \param force_i    Force on center i
 * \param force_j    Force on center j
 */
template<class BasicVector>
HOST_DEVICE_FUN HOST_DEVICE_INLINE void spreadTwoCenterForces(const VectorValueType_t<BasicVector> force,
                                                              const BasicVector& dx,
                                                              BasicVector*       force_i,
                                                              BasicVector*       force_j)
{
    BasicVector fij = force * dx;
    *force_i += fij;
    *force_j -= fij;
}

//! Three-center category common

/*! \brief spread force to 3 centers based on scalar force and angle
 *
 * \tparam T         The type of vector, e.g. float, double, etc
 * \param cos_theta  Angle between two vectors
 * \param force      The Force
 * \param dxIJ       Distance between centers i and j
 * \param dxJK       Distance between centers j and k
 * \param force_i    Force on center i
 * \param force_j    Force on center j
 * \param force_k    Force on center k
 */
template<class BasicVector>
HOST_DEVICE_FUN HOST_DEVICE_INLINE void spreadThreeCenterForces(VectorValueType_t<BasicVector> cos_theta,
                                                                VectorValueType_t<BasicVector> force,
                                                                const BasicVector&             dxIJ,
                                                                const BasicVector&             dxKJ,
                                                                BasicVector* force_i,
                                                                BasicVector* force_j,
                                                                BasicVector* force_k)
{
    using ValueType      = VectorValueType_t<BasicVector>;
    ValueType cos_theta2 = cos_theta * cos_theta;
    if (cos_theta2 < ValueType(1))
    {
        ValueType st    = force * util::invsqrt(ValueType(1) - cos_theta2);
        ValueType sth   = st * cos_theta;
        ValueType nrij2 = dot(dxIJ, dxIJ);
        ValueType nrkj2 = dot(dxKJ, dxKJ);

        ValueType cik = st * util::invsqrt(nrij2 * nrkj2);
        ValueType cii = sth / nrij2;
        ValueType ckk = sth / nrkj2;

        BasicVector f_i = cii * dxIJ - cik * dxKJ;
        BasicVector f_k = ckk * dxKJ - cik * dxIJ;
        BasicVector f_j = -(f_i + f_k);
        *force_i += f_i;
        *force_j += f_j;
        *force_k += f_k;
    }
}

//! Four-center category common

/*! \brief spread force to 4 centers
 *
 * \tparam T         The type of vector, e.g. float, double, etc
 * \param dxIJ       Distance between centers i and j
 * \param dxKJ       Distance between centers j and k
 * \param dxKL       Distance between centers k and l
 * \param m          Cross product of \p dxIJ x \p dxKJ
 * \param m          Cross product of \p dxKJ x \p dxKL
 * \param force_i    Force on center i
 * \param force_j    Force on center j
 * \param force_k    Force on center k
 * \param force_k    Force on center l
 */
template<class BasicVector>
HOST_DEVICE_FUN HOST_DEVICE_INLINE void spreadFourCenterForces(VectorValueType_t<BasicVector> force,
                                                               const BasicVector&             dxIJ,
                                                               const BasicVector&             dxJK,
                                                               const BasicVector&             dxKL,
                                                               const BasicVector&             m,
                                                               const BasicVector&             n,
                                                               BasicVector* force_i,
                                                               BasicVector* force_j,
                                                               BasicVector* force_k,
                                                               BasicVector* force_l)
{
    using ValueType    = VectorValueType_t<BasicVector>;
    ValueType norm2_m  = dot(m, m);
    ValueType norm2_n  = dot(n, n);
    ValueType norm2_jk = dot(dxJK, dxJK);
    ValueType toler    = norm2_jk * MachineEpsilon<ValueType>::value;
    if ((norm2_m > toler) && (norm2_n > toler))
    {
        ValueType rcp_norm2_jk = ValueType(1.0) / norm2_jk;
        ValueType norm_jk      = std::sqrt(norm2_jk);

        ValueType   a   = -force * norm_jk / norm2_m;
        BasicVector f_i = a * m;

        ValueType   b   = force * norm_jk / norm2_n;
        BasicVector f_l = b * n;

        ValueType   p    = rcp_norm2_jk * dot(dxIJ, dxJK);
        ValueType   q    = rcp_norm2_jk * dot(dxKL, dxJK);
        BasicVector svec = p * f_i - q * f_l;

        BasicVector f_j = svec - f_i;
        BasicVector f_k = -svec - f_l;

        *force_i += f_i;
        *force_j += f_j;
        *force_k += f_k;
        *force_l += f_l;
    }
}

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_KERNELS_HPP
