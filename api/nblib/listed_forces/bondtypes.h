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
 * Implements nblib supported bondtypes
 *
 * We choose to forward comparison operations to the
 * corresponding std::tuple comparison operations.
 * In order to do that without temporary copies,
 * we employ std::tie, which requires lvalues as input.
 * For this reason, bond type parameter getters are implemented
 * with a const lvalue reference return.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_LISTEDFORCES_BONDTYPES_H
#define NBLIB_LISTEDFORCES_BONDTYPES_H

#include <array>

#include "nblib/particletype.h"
#include "nblib/ppmap.h"
#include "nblib/util/user.h"

namespace nblib
{
using Name          = std::string;
using ForceConstant = real;
using EquilDistance = real;
using Exponent      = real;

using Degrees = StrongType<real, struct DegreeParameter>;
using Radians = StrongType<real, struct RadianParameter>;

/*! \brief Basic template for interactions with 2 parameters named forceConstant and equilDistance
 *
 * \tparam Phantom unused template parameter for type distinction
 *
 * Distinct bond types can be generated from this template with using declarations
 * and declared, but undefined structs. For example:
 * using HarmonicBondType = TwoParameterInteraction<struct HarmonicBondTypeParameter>;
 * Note that HarmonicBondTypeParameter does not have to be defined.
 */
template<class Phantom>
class TwoParameterInteraction
{
public:
    TwoParameterInteraction() = default;
    TwoParameterInteraction(ForceConstant f, EquilDistance d) : forceConstant_(f), equilDistance_(d)
    {
    }

    [[nodiscard]] const ForceConstant& forceConstant() const { return forceConstant_; }
    [[nodiscard]] const EquilDistance& equilDistance() const { return equilDistance_; }

private:
    ForceConstant forceConstant_;
    EquilDistance equilDistance_;
};

template<class Phantom>
inline bool operator<(const TwoParameterInteraction<Phantom>& a, const TwoParameterInteraction<Phantom>& b)
{
    return std::tie(a.forceConstant(), a.equilDistance())
           < std::tie(b.forceConstant(), b.equilDistance());
}

template<class Phantom>
inline bool operator==(const TwoParameterInteraction<Phantom>& a, const TwoParameterInteraction<Phantom>& b)
{
    return std::tie(a.forceConstant(), a.equilDistance())
           == std::tie(b.forceConstant(), b.equilDistance());
}

/*! \brief harmonic bond type
 *
 *  It represents the interaction of the form
 *  V(r; forceConstant, equilDistance) = 0.5 * forceConstant * (r - equilDistance)^2
 */
using HarmonicBondType = TwoParameterInteraction<struct HarmonicBondTypeParameter>;


/*! \brief GROMOS bond type
 *
 * It represents the interaction of the form
 * V(r; forceConstant, equilDistance) = 0.25 * forceConstant * (r^2 - equilDistance^2)^2
 */
using G96BondType = TwoParameterInteraction<struct G96BondTypeParameter>;


/*! \brief FENE bond type
 *
 * It represents the interaction of the form
 * V(r; forceConstant, equilDistance) = - 0.5 * forceConstant * equilDistance^2 * log( 1 - (r / equilDistance)^2)
 */
using FENEBondType = TwoParameterInteraction<struct FENEBondTypeParameter>;


/*! \brief Half-attractive quartic bond type
 *
 * It represents the interaction of the form
 * V(r; forceConstant, equilDistance) = 0.5 * forceConstant * (r - equilDistance)^4
 */
using HalfAttractiveQuarticBondType =
        TwoParameterInteraction<struct HalfAttractiveQuarticBondTypeParameter>;


/*! \brief Cubic bond type
 *
 * It represents the interaction of the form
 * V(r; quadraticForceConstant, cubicForceConstant, equilDistance) = quadraticForceConstant * (r -
 * equilDistance)^2 + quadraticForceConstant * cubicForceConstant * (r - equilDistance)
 */
struct CubicBondType
{
    CubicBondType() = default;
    CubicBondType(ForceConstant fq, ForceConstant fc, EquilDistance d) :
        quadraticForceConstant_(fq),
        cubicForceConstant_(fc),
        equilDistance_(d)
    {
    }

    [[nodiscard]] const ForceConstant& quadraticForceConstant() const
    {
        return quadraticForceConstant_;
    }
    [[nodiscard]] const ForceConstant& cubicForceConstant() const { return cubicForceConstant_; }
    [[nodiscard]] const EquilDistance& equilDistance() const { return equilDistance_; }

private:
    ForceConstant quadraticForceConstant_;
    ForceConstant cubicForceConstant_;
    EquilDistance equilDistance_;
};

inline bool operator<(const CubicBondType& a, const CubicBondType& b)
{
    return std::tie(a.quadraticForceConstant(), a.cubicForceConstant(), a.equilDistance())
           < std::tie(b.quadraticForceConstant(), b.cubicForceConstant(), b.equilDistance());
}

inline bool operator==(const CubicBondType& a, const CubicBondType& b)
{
    return std::tie(a.quadraticForceConstant(), a.cubicForceConstant(), a.equilDistance())
           == std::tie(b.quadraticForceConstant(), b.cubicForceConstant(), b.equilDistance());
}

/*! \brief Morse bond type
 *
 * It represents the interaction of the form
 * V(r; forceConstant, exponent, equilDistance) = forceConstant * ( 1 - exp( -exponent * (r - equilDistance))
 */
class MorseBondType
{
public:
    MorseBondType() = default;
    MorseBondType(ForceConstant f, Exponent e, EquilDistance d) :
        forceConstant_(f),
        exponent_(e),
        equilDistance_(d)
    {
    }

    [[nodiscard]] const ForceConstant& forceConstant() const { return forceConstant_; }
    [[nodiscard]] const Exponent&      exponent() const { return exponent_; }
    [[nodiscard]] const EquilDistance& equilDistance() const { return equilDistance_; }

private:
    ForceConstant forceConstant_;
    Exponent      exponent_;
    EquilDistance equilDistance_;
};

inline bool operator<(const MorseBondType& a, const MorseBondType& b)
{
    return std::tie(a.forceConstant(), a.exponent(), a.equilDistance())
           < std::tie(b.forceConstant(), b.exponent(), b.equilDistance());
}

inline bool operator==(const MorseBondType& a, const MorseBondType& b)
{
    return std::tie(a.forceConstant(), a.exponent(), a.equilDistance())
           == std::tie(b.forceConstant(), b.exponent(), b.equilDistance());
}


/*! \brief default angle type
 *
 * Note: the angle is always stored as radians internally
 */
struct DefaultAngle : public TwoParameterInteraction<struct DefaultAngleParameter>
{
    DefaultAngle() = default;
    //! \brief construct from angle given in radians
    DefaultAngle(Radians angle, ForceConstant f) :
        TwoParameterInteraction<struct DefaultAngleParameter>{ f, angle }
    {
    }

    //! \brief construct from angle given in degrees
    DefaultAngle(Degrees angle, ForceConstant f) :
        TwoParameterInteraction<struct DefaultAngleParameter>{ f, angle * DEG2RAD }
    {
    }
};

/*! \brief Proper Dihedral Implementation
 */
class ProperDihedral
{
public:
    using Multiplicity = int;

    ProperDihedral() = default;
    ProperDihedral(Radians phi, ForceConstant f, Multiplicity m) :
        phi_(phi),
        forceConstant_(f),
        multiplicity_(m)
    {
    }
    ProperDihedral(Degrees phi, ForceConstant f, Multiplicity m) :
        phi_(phi * DEG2RAD),
        forceConstant_(f),
        multiplicity_(m)
    {
    }

    [[nodiscard]] const EquilDistance& equilDistance() const { return phi_; }
    [[nodiscard]] const ForceConstant& forceConstant() const { return forceConstant_; }
    [[nodiscard]] const Multiplicity&  multiplicity() const { return multiplicity_; }

private:
    EquilDistance phi_;
    ForceConstant forceConstant_;
    Multiplicity  multiplicity_;
};

inline bool operator<(const ProperDihedral& a, const ProperDihedral& b)
{
    return std::tie(a.equilDistance(), a.forceConstant(), a.multiplicity())
           < std::tie(b.equilDistance(), b.forceConstant(), b.multiplicity());
}

inline bool operator==(const ProperDihedral& a, const ProperDihedral& b)
{
    return std::tie(a.equilDistance(), a.forceConstant(), a.multiplicity())
           == std::tie(b.equilDistance(), b.forceConstant(), b.multiplicity());
}


/*! \brief Improper Dihedral Implementation
 */
struct ImproperDihedral : public TwoParameterInteraction<struct ImproperDihdedralParameter>
{
    ImproperDihedral() = default;
    ImproperDihedral(Radians phi, ForceConstant f) :
        TwoParameterInteraction<struct ImproperDihdedralParameter>{ f, phi }
    {
    }
    ImproperDihedral(Degrees phi, ForceConstant f) :
        TwoParameterInteraction<struct ImproperDihdedralParameter>{ f, phi * DEG2RAD }
    {
    }
};

/*! \brief Ryckaert-Belleman Dihedral Implementation
 */
class RyckaertBellemanDihedral
{
public:
    RyckaertBellemanDihedral() = default;
    RyckaertBellemanDihedral(real p1, real p2, real p3, real p4, real p5, real p6) :
        parameters_{ p1, p2, p3, p4, p5, p6 }
    {
    }

    const real& operator[](std::size_t i) const { return parameters_[i]; }

    [[nodiscard]] const std::array<real, 6>& parameters() const { return parameters_; }

    [[nodiscard]] std::size_t size() const { return parameters_.size(); }

private:
    std::array<real, 6> parameters_;
};

inline bool operator<(const RyckaertBellemanDihedral& a, const RyckaertBellemanDihedral& b)
{
    return a.parameters() < b.parameters();
}

inline bool operator==(const RyckaertBellemanDihedral& a, const RyckaertBellemanDihedral& b)
{
    return a.parameters() == b.parameters();
}


/*! \brief Type for 5-center interaction (C-MAP)
 *
 *  Note: no kernels currently implemented
 */
class Default5Center
{
public:
    Default5Center() = default;
    Default5Center(Radians phi, Radians psi, ForceConstant fphi, ForceConstant fpsi) :
        phi_(phi),
        psi_(psi),
        fphi_(fphi),
        fpsi_(fpsi)
    {
    }

    [[nodiscard]] const Radians&       phi() const { return phi_; }
    [[nodiscard]] const Radians&       psi() const { return psi_; }
    [[nodiscard]] const ForceConstant& fphi() const { return fphi_; }
    [[nodiscard]] const ForceConstant& fpsi() const { return fpsi_; }

private:
    Radians       phi_, psi_;
    ForceConstant fphi_, fpsi_;
};

inline bool operator<(const Default5Center& a, const Default5Center& b)
{
    return std::tie(a.phi(), a.psi(), a.fphi(), a.fpsi())
           < std::tie(b.phi(), b.psi(), b.fphi(), b.fpsi());
}

inline bool operator==(const Default5Center& a, const Default5Center& b)
{
    return std::tie(a.phi(), a.psi(), a.fphi(), a.fpsi())
           == std::tie(b.phi(), b.psi(), b.fphi(), b.fpsi());
}


} // namespace nblib
#endif // NBLIB_LISTEDFORCES_BONDTYPES_H
