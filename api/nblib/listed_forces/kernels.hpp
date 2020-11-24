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

#include "gromacs/math/vec.h"
#include "nblib/basicdefinitions.h"
#include "bondtypes.h"

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
    return harmonicScalarForce(bond.forceConstant(), bond.equilDistance(), dr);
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

    T force = -k * dx * x;
    T epot = 0.25 * k * dx2;

    return std::make_tuple(force, epot);

    /* That was 7 flops */
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

    T force = -kk * dx * x;
    T epot = 0.25 * kk * dx2;
    // TODO: Check if this is 1/2 or 1/4
    T dvdlambda = 0.5 * (kB - kA) * dx2 + (xA - xB) * kk * dx;

    return std::make_tuple(force, epot, dvdlambda);

    /* That was 21 flops */
}

//! Abstraction layer for different 2-center bonds. Forth power case
template <class T>
inline auto bondKernel(T dr, const G96BondType& bond)
{
    // NOTE: Not assuming GROMACS' convention of storing squared bond.equilDistance() for this type
    return g96ScalarForce(bond.forceConstant(), bond.equilDistance() * bond.equilDistance(), dr * dr);
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
    T force = -2.0 * beta * exponent * omexp;     /*  4 */

    return std::make_tuple(force, epot);

    /* That was 23 flops */
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
    T force = -2.0 * beta * exponent * omexp;     /*  4 */

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
    T force = -k * x / omx2_ox02;

    return std::make_tuple(force, epot);

    /* That was 24 flops */
}

// TODO: Implement the free energy version of FENE (finitely extensible nonlinear elastic) bond types

//! Abstraction layer for different 2-center bonds. FENE case
template <class T>
inline auto bondKernel(T dr, const FENEBondType& bond)
{
    return FENEScalarForce(bond.forceConstant(), bond.equilDistance(), dr);
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
    return halfAttractiveScalarForce(bond.forceConstant(), bond.equilDistance(), dr);
}


//! Three-center interaction type kernels

// linear angles go here
// quartic angles go here

//! Three-center interaction type dispatch

template <class T>
inline auto threeCenterKernel(T dr, const DefaultAngle& angle)
{
    return harmonicScalarForce(angle.forceConstant(), angle.equilDistance(), dr);
}


//! \brief Computes and returns the proper dihedral force
template <class T>
inline auto fourCenterKernel(T phi, const ProperDihedral& properDihedral)
{
    const T deltaPhi = properDihedral.multiplicity() * phi - properDihedral.equilDistance();
    const T force = -properDihedral.forceConstant() * properDihedral.multiplicity() * std::sin(deltaPhi);
    const T ePot = properDihedral.forceConstant() * ( 1 + std::cos(deltaPhi) );
    return std::make_tuple(force, ePot);
}


//! \brief Ensure that a geometric quantity lies in (-pi, pi)
static inline void makeAnglePeriodic(real& angle)
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

//! \brief Computes and returns a dihedral phi angle
static inline real dihedralPhi(rvec dxIJ, rvec dxKJ, rvec dxKL, rvec m, rvec n)
{
    cprod(dxIJ, dxKJ, m);
    cprod(dxKJ, dxKL, n);
    real phi  = gmx_angle(m, n);
    real ipr  = iprod(dxIJ, n);
    real sign = (ipr < 0.0) ? -1.0 : 1.0;
    phi       = sign * phi;
    return phi;
}

//! \brief Computes and returns the improper dihedral force
template <class T>
inline auto fourCenterKernel(T phi, const ImproperDihedral& improperDihedral)
{
    T deltaPhi = phi - improperDihedral.equilDistance();
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

/*! \brief Spreads and accumulates the bonded forces to the two atoms and adds the virial contribution when needed
 *
 * \p shiftIndex is used as the periodic shift.
 */
template <class T>
inline void spreadTwoCenterForces(const T bondForce,
                                  const gmx::RVec& dx,
                                  gmx::RVec* force_i,
                                  gmx::RVec* force_j)
{
    for (int m = 0; m < dimSize; m++) /*  15          */
    {
        const T fij = bondForce * dx[m];
        (*force_i)[m] += fij;
        (*force_j)[m] -= fij;
    }
}

//! Three-center category common

/*! \brief spread force to 3 centers based on scalar force and angle
 *
 * @tparam T
 * @param cos_theta
 * @param force
 * @param r_ij
 * @param r_kj
 * @param force_i
 * @param force_j
 * @param force_k
 */
template <class T>
inline void spreadThreeCenterForces(T cos_theta,
                             T force,
                             const gmx::RVec& r_ij,
                             const gmx::RVec& r_kj,
                             gmx::RVec* force_i,
                             gmx::RVec* force_j,
                             gmx::RVec* force_k)
{
    T cos_theta2 = cos_theta * cos_theta;
    if (cos_theta2 < 1)
    {
        T st    = force / std::sqrt(1 - cos_theta2); /*  12		*/
        T sth   = st * cos_theta;                      /*   1		*/
        T nrij2 = dot(r_ij, r_ij);                   /*   5		*/
        T nrkj2 = dot(r_kj, r_kj);                   /*   5		*/

        T nrij_1 = 1.0 / std::sqrt(nrij2); /*  10		*/
        T nrkj_1 = 1.0 / std::sqrt(nrkj2); /*  10		*/

        T cik = st * nrij_1 * nrkj_1;  /*   2		*/
        T cii = sth * nrij_1 * nrij_1; /*   2		*/
        T ckk = sth * nrkj_1 * nrkj_1; /*   2		*/

        gmx::RVec f_i{0, 0, 0};
        gmx::RVec f_j{0, 0, 0};
        gmx::RVec f_k{0, 0, 0};
        for (int m = 0; m < dimSize; m++)
        { /*  39		*/
            f_i[m] = -(cik * r_kj[m] - cii * r_ij[m]);
            f_k[m] = -(cik * r_ij[m] - ckk * r_kj[m]);
            f_j[m] = -f_i[m] - f_k[m];
            (*force_i)[m] += f_i[m];
            (*force_j)[m] += f_j[m];
            (*force_k)[m] += f_k[m];
        }
    } /* 161 TOTAL	*/
}

//! Four-center category common
template <class T>
inline void spreadFourCenterForces(T force, rvec dxIJ, rvec dxJK, rvec dxKL, rvec m, rvec n,
                            gmx::RVec* force_i,
                            gmx::RVec* force_j,
                            gmx::RVec* force_k,
                            gmx::RVec* force_l)
{
    rvec f_i, f_j, f_k, f_l;
    rvec uvec, vvec, svec;
    T iprm  = iprod(m, m);       /*  5    */
    T iprn  = iprod(n, n);       /*  5	*/
    T nrkj2 = iprod(dxJK, dxJK); /*  5	*/
    T toler = nrkj2 * GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        T nrkj_1 = gmx::invsqrt(nrkj2);  /* 10	*/
        T nrkj_2 = nrkj_1 * nrkj_1;      /*  1	*/
        T nrkj   = nrkj2 * nrkj_1;       /*  1	*/
        T a      = -force * nrkj / iprm; /* 11	*/
        svmul(a, m, f_i);              /*  3	*/
        T b = force * nrkj / iprn;       /* 11	*/
        svmul(b, n, f_l);              /*  3  */
        T p = iprod(dxIJ, dxJK);         /*  5	*/
        p *= nrkj_2;                   /*  1	*/
        T q = iprod(dxKL, dxJK);         /*  5	*/
        q *= nrkj_2;                   /*  1	*/
        svmul(p, f_i, uvec);           /*  3	*/
        svmul(q, f_l, vvec);           /*  3	*/
        rvec_sub(uvec, vvec, svec);    /*  3	*/
        rvec_sub(f_i, svec, f_j);      /*  3	*/
        rvec_add(f_l, svec, f_k);      /*  3	*/
        rvec_inc(force_i->as_vec(), f_i);           /*  3	*/
        rvec_dec(force_j->as_vec(), f_j);           /*  3	*/
        rvec_dec(force_k->as_vec(), f_k);           /*  3	*/
        rvec_inc(force_l->as_vec(), f_l);           /*  3	*/
    }
}

} // namespace nblib
#endif // NBLIB_LISTEDFORCES_KERNELS_HPP
