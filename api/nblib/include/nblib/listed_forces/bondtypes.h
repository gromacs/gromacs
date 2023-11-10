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
 * Implements nblib supported bondtypes
 *
 * We choose to forward comparison operations to the
 * corresponding std::tuple comparison operations.
 * In order to do that without temporary copies,
 * we employ util::tie, which requires lvalues as input.
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

#include "nblib/exception.h"
#include "nblib/particletype.h"
#include "nblib/util/array.hpp"
#include "nblib/util/tuple.hpp"
#include "nblib/util/util.hpp"

namespace nblib
{

using ForceConstant = real;
using EquilConstant = real;
using Exponent      = real;

using Degrees = StrongType<real, struct DegreeParameter>;
using Radians = StrongType<real, struct RadianParameter>;

/*! \brief Basic template for interactions with 2 parameters named forceConstant and equilConstant
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
    HOST_DEVICE_FUN
    TwoParameterInteraction(ForceConstant f, EquilConstant d) : forceConstant_(f), equilConstant_(d)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& forceConstant() const
    {
        return forceConstant_;
    }
    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilConstant() const
    {
        return equilConstant_;
    }

private:
    ForceConstant forceConstant_;
    EquilConstant equilConstant_;
};

template<class Phantom>
HOST_DEVICE_FUN inline bool operator<(const TwoParameterInteraction<Phantom>& a,
                                      const TwoParameterInteraction<Phantom>& b)
{
    return util::tie(a.forceConstant(), a.equilConstant())
           < util::tie(b.forceConstant(), b.equilConstant());
}

template<class Phantom>
HOST_DEVICE_FUN inline bool operator==(const TwoParameterInteraction<Phantom>& a,
                                       const TwoParameterInteraction<Phantom>& b)
{
    return util::tie(a.forceConstant(), a.equilConstant())
           == util::tie(b.forceConstant(), b.equilConstant());
}

/*! \brief harmonic bond type
 *
 *  It represents the interaction of the form
 *  V(r, forceConstant, equilConstant) = 0.5 * forceConstant * (r - equilConstant)^2
 */
using HarmonicBondType = TwoParameterInteraction<struct HarmonicBondTypeParameter>;


/*! \brief GROMOS bond type
 *
 * It represents the interaction of the form
 * V(r, forceConstant, equilConstant) = 0.25 * forceConstant * (r^2 - equilConstant^2)^2
 */
class G96BondType : public TwoParameterInteraction<struct G96BondTypeParameter>
{
public:
    G96BondType() = default;
    //! \brief Store square of equilibrium distance
    HOST_DEVICE_FUN G96BondType(ForceConstant f, EquilConstant equilConstant) :
        TwoParameterInteraction<struct G96BondTypeParameter>{ f, equilConstant * equilConstant }
    {
    }
};


/*! \brief FENE bond type
 *
 * It represents the interaction of the form
 * V(r, forceConstant, equilConstant) = - 0.5 * forceConstant * equilConstant^2 * log( 1 - (r / equilConstant)^2)
 */
using FENEBondType = TwoParameterInteraction<struct FENEBondTypeParameter>;


/*! \brief Half-attractive quartic bond type
 *
 * It represents the interaction of the form
 * V(r, forceConstant, equilConstant) = 0.5 * forceConstant * (r - equilConstant)^4
 */
using HalfAttractiveQuarticBondType =
        TwoParameterInteraction<struct HalfAttractiveQuarticBondTypeParameter>;


/*! \brief Cubic bond type
 *
 * It represents the interaction of the form
 * V(r, quadraticForceConstant, cubicForceConstant, equilConstant) = quadraticForceConstant * (r -
 * equilConstant)^2 + quadraticForceConstant * cubicForceConstant * (r - equilConstant)
 */
class CubicBondType
{
public:
    CubicBondType() = default;
    HOST_DEVICE_FUN CubicBondType(ForceConstant fq, ForceConstant fc, EquilConstant d) :
        quadraticForceConstant_(fq), cubicForceConstant_(fc), equilDistance_(d)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& quadraticForceConstant() const
    {
        return quadraticForceConstant_;
    }
    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& cubicForceConstant() const
    {
        return cubicForceConstant_;
    }
    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilDistance() const
    {
        return equilDistance_;
    }

private:
    ForceConstant quadraticForceConstant_;
    ForceConstant cubicForceConstant_;
    EquilConstant equilDistance_;
};

HOST_DEVICE_FUN inline bool operator<(const CubicBondType& a, const CubicBondType& b)
{
    return util::tie(a.quadraticForceConstant(), a.cubicForceConstant(), a.equilDistance())
           < util::tie(b.quadraticForceConstant(), b.cubicForceConstant(), b.equilDistance());
}

HOST_DEVICE_FUN inline bool operator==(const CubicBondType& a, const CubicBondType& b)
{
    return util::tie(a.quadraticForceConstant(), a.cubicForceConstant(), a.equilDistance())
           == util::tie(b.quadraticForceConstant(), b.cubicForceConstant(), b.equilDistance());
}

/*! \brief Morse bond type
 *
 * It represents the interaction of the form
 * V(r, forceConstant, exponent, equilConstant) = forceConstant * ( 1 - exp( -exponent * (r - equilConstant))
 */
class MorseBondType
{
public:
    MorseBondType() = default;
    HOST_DEVICE_FUN MorseBondType(ForceConstant f, Exponent e, EquilConstant d) :
        forceConstant_(f), exponent_(e), equilDistance_(d)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& forceConstant() const
    {
        return forceConstant_;
    }
    [[nodiscard]] HOST_DEVICE_FUN const Exponent& exponent() const { return exponent_; }
    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilDistance() const
    {
        return equilDistance_;
    }

private:
    ForceConstant forceConstant_;
    Exponent      exponent_;
    EquilConstant equilDistance_;
};

HOST_DEVICE_FUN inline bool operator<(const MorseBondType& a, const MorseBondType& b)
{
    return util::tie(a.forceConstant(), a.exponent(), a.equilDistance())
           < util::tie(b.forceConstant(), b.exponent(), b.equilDistance());
}

HOST_DEVICE_FUN inline bool operator==(const MorseBondType& a, const MorseBondType& b)
{
    return util::tie(a.forceConstant(), a.exponent(), a.equilDistance())
           == util::tie(b.forceConstant(), b.exponent(), b.equilDistance());
}

/*! \brief Non-Bonded Pair Interaction Type
 *
 * It represents the interaction of the form
 * of LJ interactions, but occur between atoms
 * of the same bonded chain
 * V(r, c_6, c_12) = c_12*(r^-12) - c_6*(r^-6)
 */
class PairLJType
{
public:
    PairLJType() = default;
    HOST_DEVICE_FUN PairLJType(C6 c6, C12 c12) : c6_(c6), c12_(c12) {}

    [[nodiscard]] HOST_DEVICE_FUN const C6& c6() const { return c6_; }
    [[nodiscard]] HOST_DEVICE_FUN const C12& c12() const { return c12_; }

private:
    C6  c6_;
    C12 c12_;
};

HOST_DEVICE_FUN inline bool operator<(const PairLJType& a, const PairLJType& b)
{
    return util::tie(a.c6(), a.c12()) < util::tie(b.c6(), b.c12());
}

HOST_DEVICE_FUN inline bool operator==(const PairLJType& a, const PairLJType& b)
{
    return util::tie(a.c6(), a.c12()) == util::tie(b.c6(), b.c12());
}

/*! \brief Non-Bonded Pair Interaction Type
 *
 * It represents the interaction of the form
 * of LJ interactions, but occur between atoms
 * of the same bonded chain
 * V(r, c_6, c_12) = c_12*(r^-12) - c_6*(r^-6)
 */
class PairLJChargeType
{
public:
    using FudgeFactor  = real;
    using Charge       = real;
    PairLJChargeType() = default;
    HOST_DEVICE_FUN PairLJChargeType(C6 c6, C12 c12, Charge qi, Charge qj, FudgeFactor ff) :
        c6_(c6), c12_(c12), qi_(qi), qj_(qj), ff_(ff)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const C6& c6() const { return c6_; }
    [[nodiscard]] HOST_DEVICE_FUN const C12& c12() const { return c12_; }
    [[nodiscard]] HOST_DEVICE_FUN const Charge& qi() const { return qi_; }
    [[nodiscard]] HOST_DEVICE_FUN const Charge& qj() const { return qj_; }
    [[nodiscard]] HOST_DEVICE_FUN const FudgeFactor& ff() const { return ff_; }

private:
    C6          c6_;
    C12         c12_;
    Charge      qi_;
    Charge      qj_;
    FudgeFactor ff_;
};

HOST_DEVICE_FUN inline bool operator<(const PairLJChargeType& a, const PairLJChargeType& b)
{
    return util::tie(a.c6(), a.c12(), a.qi(), a.qj(), a.ff())
           < util::tie(b.c6(), b.c12(), b.qi(), b.qj(), b.ff());
}

HOST_DEVICE_FUN inline bool operator==(const PairLJChargeType& a, const PairLJChargeType& b)
{
    return util::tie(a.c6(), a.c12(), a.qi(), a.qj(), a.ff())
           == util::tie(b.c6(), b.c12(), b.qi(), b.qj(), b.ff());
}

/*! \brief Simple polarization
 */
class SimplePolarization
{
public:
    SimplePolarization() = default;
    HOST_DEVICE_FUN explicit SimplePolarization(ForceConstant alpha) :
        alpha_(alpha), coefficient_(ONE_4PI_EPS0 / alpha)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& alpha() const { return alpha_; }
    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& coefficient() const { return coefficient_; }

private:
    ForceConstant alpha_;
    ForceConstant coefficient_;
};

HOST_DEVICE_FUN inline bool operator<(const SimplePolarization& a, const SimplePolarization& b)
{
    return util::tie(a.alpha(), a.coefficient()) < util::tie(b.alpha(), b.coefficient());
}

HOST_DEVICE_FUN inline bool operator==(const SimplePolarization& a, const SimplePolarization& b)
{
    return util::tie(a.alpha(), a.coefficient()) == util::tie(b.alpha(), b.coefficient());
}

/*! \brief Basic template for interactions with 2 parameters named forceConstant and equilAngle
 *
 * \tparam Phantom unused template parameter for type distinction
 *
 * Distinct angle types can be generated from this template with using declarations
 * and declared, but undefined structs. For example:
 * using HarmonicAngleType = AngleInteractionType<struct HarmonicAngleParameter>;
 * HarmonicAngleParameter does not have to be defined.
 *
 * Note: the angle is always stored as radians internally
 */
template<class Phantom>
class AngleInteractionType : public TwoParameterInteraction<Phantom>
{
public:
    AngleInteractionType() = default;
    //! \brief construct from angle given in radians
    HOST_DEVICE_FUN AngleInteractionType(ForceConstant f, Radians angle) :
        TwoParameterInteraction<Phantom>{ f, angle }
    {
    }

    //! \brief construct from angle given in degrees
    HOST_DEVICE_FUN AngleInteractionType(ForceConstant f, Degrees angle) :
        TwoParameterInteraction<Phantom>{ f, angle * DEG2RAD }
    {
    }
};

/*! \brief Harmonic angle type
 *
 * It represents the interaction of the form
 * V(theta, forceConstant, equilConstant) = 0.5 * forceConstant * (theta - equilConstant)^2
 */
using HarmonicAngle = AngleInteractionType<struct HarmonicAngleParameter>;

/*! \brief Urey-Bradly Interactions
 *
 * Dummy no op placeholder but implemented internally as a combination of a harmonic bond
 * and a harmonic angle
 */
class UreyBradley
{
};

/*! \brief linear angle type
 *
 * It represents the interaction of the form
 * V(theta, forceConstant, a) = 0.5 * forceConstant * (dr)^2
 * where dr = - a * r_ij - (1 - a) * r_kj
 */
using LinearAngle = TwoParameterInteraction<struct LinearAngleParameter>;

/*! \brief Basic template for angle types that use the cosines of their equilibrium angles
 *         in their potential expression
 */
template<class Phantom>
class CosineParamAngle : public TwoParameterInteraction<Phantom>
{
public:
    CosineParamAngle() = default;
    //! \brief construct from angle given in radians
    HOST_DEVICE_FUN CosineParamAngle(ForceConstant f, Radians angle) :
        TwoParameterInteraction<Phantom>{ f, std::cos(angle) }
    {
    }

    //! \brief construct from angle given in degrees
    HOST_DEVICE_FUN CosineParamAngle(ForceConstant f, Degrees angle) :
        TwoParameterInteraction<Phantom>{ f, std::cos(angle * DEG2RAD) }
    {
    }
};

/*! \brief G96 or Cosine-based angle type
 *
 * This represents the interaction of the form
 * V(cos(theta), forceConstant, cos(equilConstant)) = 0.5 * forceConstant * (cos(theta) - cos(equilConstant))^2
 */
using G96Angle = CosineParamAngle<struct G96AngleParameter>;

/*! \brief Restricted angle type
 *
 * This represents the interaction of the form
 * V(cos(theta), forceConstant, cos(equilConstant)) =
 *     0.5 * forceConstant * (cos(theta) - cos(equilConstant))^2 / (sin(theta))^2
 */
using RestrictedAngle = CosineParamAngle<struct RestrictedAngleParameter>;

/*! \brief Quartic angle type
 *
 * It represents the interaction of the form of a fourth order polynomial
 * V(theta, forceConstant, equilConstant) = sum[i = 0 -> 4](forceConstant_i * (theta - equilConstant)^i
 */
class QuarticAngle
{
public:
    QuarticAngle() = default;
    //! \brief construct from given angle in radians
    HOST_DEVICE_FUN QuarticAngle(ForceConstant f0,
                                 ForceConstant f1,
                                 ForceConstant f2,
                                 ForceConstant f3,
                                 ForceConstant f4,
                                 Radians       angle) :
        forceConstants_{ f0, f1, f2, f3, f4 }, equilConstant_(angle)
    {
    }

    //! \brief construct from given angle in degrees
    HOST_DEVICE_FUN QuarticAngle(ForceConstant f0,
                                 ForceConstant f1,
                                 ForceConstant f2,
                                 ForceConstant f3,
                                 ForceConstant f4,
                                 Degrees       angle) :
        forceConstants_{ f0, f1, f2, f3, f4 }, equilConstant_(angle * DEG2RAD)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& forceConstant(int order) const
    {
        assert(order < 5);
        return forceConstants_[order];
    }

    [[nodiscard]] HOST_DEVICE_FUN const Radians& equilConstant() const { return equilConstant_; }

private:
    ForceConstant forceConstants_[5];
    Radians       equilConstant_;
};

HOST_DEVICE_FUN inline bool operator<(const QuarticAngle& a, const QuarticAngle& b)
{
    return util::tie(a.forceConstant(0),
                     a.forceConstant(1),
                     a.forceConstant(2),
                     a.forceConstant(3),
                     a.forceConstant(4),
                     a.equilConstant())
           < util::tie(b.forceConstant(0),
                       b.forceConstant(1),
                       b.forceConstant(2),
                       b.forceConstant(3),
                       b.forceConstant(4),
                       b.equilConstant());
}

HOST_DEVICE_FUN inline bool operator==(const QuarticAngle& a, const QuarticAngle& b)
{
    return util::tie(a.forceConstant(0),
                     a.forceConstant(1),
                     a.forceConstant(2),
                     a.forceConstant(3),
                     a.forceConstant(4),
                     a.equilConstant())
           == util::tie(b.forceConstant(0),
                        b.forceConstant(1),
                        b.forceConstant(2),
                        b.forceConstant(3),
                        b.forceConstant(4),
                        b.equilConstant());
}


/*! \brief Cross bond-bond interaction type
 */
class CrossBondBond
{
public:
    CrossBondBond() = default;
    HOST_DEVICE_FUN CrossBondBond(ForceConstant f, EquilConstant r0ij, EquilConstant r0kj) :
        forceConstant_(f), r0ij_(r0ij), r0kj_(r0kj)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& forceConstant() const
    {
        return forceConstant_;
    }
    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilDistanceIJ() const { return r0ij_; }
    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilDistanceKJ() const { return r0kj_; }

private:
    ForceConstant forceConstant_;
    EquilConstant r0ij_;
    EquilConstant r0kj_;
};

HOST_DEVICE_FUN inline bool operator<(const CrossBondBond& a, const CrossBondBond& b)
{
    return util::tie(a.forceConstant(), a.equilDistanceIJ(), a.equilDistanceKJ())
           < util::tie(b.forceConstant(), b.equilDistanceIJ(), b.equilDistanceKJ());
}

HOST_DEVICE_FUN inline bool operator==(const CrossBondBond& a, const CrossBondBond& b)
{
    return util::tie(a.forceConstant(), a.equilDistanceIJ(), a.equilDistanceKJ())
           == util::tie(b.forceConstant(), b.equilDistanceIJ(), b.equilDistanceKJ());
}

/*! \brief Cross bond-angle interaction type
 */
class CrossBondAngle
{
public:
    CrossBondAngle() = default;
    HOST_DEVICE_FUN CrossBondAngle(ForceConstant f, EquilConstant r0ij, EquilConstant r0kj, EquilConstant r0ik) :
        forceConstant_(f), r0ij_(r0ij), r0kj_(r0kj), r0ik_(r0ik)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& forceConstant() const
    {
        return forceConstant_;
    }
    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilDistanceIJ() const { return r0ij_; }
    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilDistanceKJ() const { return r0kj_; }
    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilDistanceIK() const { return r0ik_; }

private:
    ForceConstant forceConstant_;
    EquilConstant r0ij_;
    EquilConstant r0kj_;
    EquilConstant r0ik_;
};

HOST_DEVICE_FUN inline bool operator<(const CrossBondAngle& a, const CrossBondAngle& b)
{
    return util::tie(a.forceConstant(), a.equilDistanceIJ(), a.equilDistanceKJ(), a.equilDistanceIK())
           < util::tie(b.forceConstant(), b.equilDistanceIJ(), b.equilDistanceKJ(), b.equilDistanceIK());
}

HOST_DEVICE_FUN inline bool operator==(const CrossBondAngle& a, const CrossBondAngle& b)
{
    return util::tie(a.forceConstant(), a.equilDistanceIJ(), a.equilDistanceKJ(), a.equilDistanceIK())
           == util::tie(b.forceConstant(), b.equilDistanceIJ(), b.equilDistanceKJ(), b.equilDistanceIK());
}

template<class Phantom>
class DihedralInteraction
{
public:
    using Multiplicity = int;

    DihedralInteraction() = default;
    HOST_DEVICE_FUN DihedralInteraction(Radians phi, ForceConstant f, Multiplicity m) :
        phi_(phi), forceConstant_(f), multiplicity_(m)
    {
    }
    HOST_DEVICE_FUN DihedralInteraction(Degrees phi, ForceConstant f, Multiplicity m) :
        phi_(phi * DEG2RAD), forceConstant_(f), multiplicity_(m)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const EquilConstant& equilDistance() const { return phi_; }
    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& forceConstant() const
    {
        return forceConstant_;
    }
    [[nodiscard]] HOST_DEVICE_FUN const Multiplicity& multiplicity() const { return multiplicity_; }

private:
    EquilConstant phi_;
    ForceConstant forceConstant_;
    Multiplicity  multiplicity_;
};

template<class Phantom>
HOST_DEVICE_FUN inline bool operator<(const DihedralInteraction<Phantom>& a,
                                      const DihedralInteraction<Phantom>& b)
{
    return util::tie(a.equilDistance(), a.forceConstant(), a.multiplicity())
           < util::tie(b.equilDistance(), b.forceConstant(), b.multiplicity());
}

template<class Phantom>
HOST_DEVICE_FUN inline bool operator==(const DihedralInteraction<Phantom>& a,
                                       const DihedralInteraction<Phantom>& b)
{
    return util::tie(a.equilDistance(), a.forceConstant(), a.multiplicity())
           == util::tie(b.equilDistance(), b.forceConstant(), b.multiplicity());
}

/*! \brief Proper Dihedral Implementation
 */
using ProperDihedral = DihedralInteraction<struct ProperDihedralParam>;

/*! \brief Proper Dihedral Implementation
 */
using ImproperProperDihedral = DihedralInteraction<struct ImproperProperDihedralParam>;

/*! \brief Improper Dihedral Implementation
 */
class ImproperDihedral : public TwoParameterInteraction<struct ImproperDihdedralParameter>
{
public:
    ImproperDihedral() = default;
    HOST_DEVICE_FUN ImproperDihedral(Radians phi, ForceConstant f) :
        TwoParameterInteraction<struct ImproperDihdedralParameter>{ f, phi }
    {
    }
    HOST_DEVICE_FUN ImproperDihedral(Degrees phi, ForceConstant f) :
        TwoParameterInteraction<struct ImproperDihdedralParameter>{ f, phi * DEG2RAD }
    {
    }
};

/*! \brief Ryckaert-Belleman Dihedral Implementation
 */
class RyckaertBellemanDihedral
{
public:
    static constexpr int RbNumParameters = 6;

    RyckaertBellemanDihedral() = default;
    HOST_DEVICE_FUN RyckaertBellemanDihedral(real p1, real p2, real p3, real p4, real p5, real p6) :
        parameters_{ p1, p2, p3, p4, p5, p6 }
    {
    }

    HOST_DEVICE_FUN const real& operator[](std::size_t i) const { return parameters_[i]; }

    [[nodiscard]] HOST_DEVICE_FUN const real* parameters() const { return parameters_; }

private:
    real parameters_[RbNumParameters];
};

HOST_DEVICE_FUN inline bool operator<(const RyckaertBellemanDihedral& a, const RyckaertBellemanDihedral& b)
{
    return util::tie(a[0], a[1], a[2], a[3], a[4], a[5]) < util::tie(b[0], b[1], b[2], b[3], b[4], b[5]);
}

HOST_DEVICE_FUN inline bool operator==(const RyckaertBellemanDihedral& a, const RyckaertBellemanDihedral& b)
{
    return util::tie(a[0], a[1], a[2], a[3], a[4], a[5]) == util::tie(b[0], b[1], b[2], b[3], b[4], b[5]);
}


/*! \brief Type for 5-center interaction (C-MAP)
 *
 *  Note: no kernels currently implemented
 */
class Default5Center
{
public:
    Default5Center() = default;
    HOST_DEVICE_FUN Default5Center(Radians phi, Radians psi, ForceConstant fphi, ForceConstant fpsi) :
        phi_(phi), psi_(psi), fphi_(fphi), fpsi_(fpsi)
    {
    }

    [[nodiscard]] HOST_DEVICE_FUN const Radians& phi() const { return phi_; }
    [[nodiscard]] HOST_DEVICE_FUN const Radians& psi() const { return psi_; }
    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& fphi() const { return fphi_; }
    [[nodiscard]] HOST_DEVICE_FUN const ForceConstant& fpsi() const { return fpsi_; }

private:
    Radians       phi_, psi_;
    ForceConstant fphi_, fpsi_;
};

HOST_DEVICE_FUN inline bool operator<(const Default5Center& a, const Default5Center& b)
{
    return util::tie(a.phi(), a.psi(), a.fphi(), a.fpsi())
           < util::tie(b.phi(), b.psi(), b.fphi(), b.fpsi());
}

HOST_DEVICE_FUN inline bool operator==(const Default5Center& a, const Default5Center& b)
{
    return util::tie(a.phi(), a.psi(), a.fphi(), a.fpsi())
           == util::tie(b.phi(), b.psi(), b.fphi(), b.fpsi());
}

class PositionRestraints
{
public:
    PositionRestraints() = default;
    HOST_DEVICE_FUN PositionRestraints(EquilConstant p0x,
                                       EquilConstant p0y,
                                       EquilConstant p0z,
                                       ForceConstant fcx,
                                       ForceConstant fcy,
                                       ForceConstant fcz) :
        position0_({ p0x, p0y, p0z }), forceConstant_({ fcx, fcy, fcz })
    {
    }

    HOST_DEVICE_FUN const EquilConstant& position0(int index) const
    {
        assert(index < 3);
        return position0_[index];
    }

    HOST_DEVICE_FUN const ForceConstant& forceConstant(int index) const
    {
        assert(index < 3);
        return forceConstant_[index];
    }

private:
    util::array<EquilConstant, dimSize> position0_;
    util::array<ForceConstant, dimSize> forceConstant_;

    friend HOST_DEVICE_FUN inline bool operator<(const PositionRestraints& a, const PositionRestraints& b)
    {
        return util::tie(a.position0_, a.forceConstant_) < util::tie(b.position0_, b.forceConstant_);
    }

    friend HOST_DEVICE_FUN inline bool operator==(const PositionRestraints& a, const PositionRestraints& b)
    {
        return util::tie(a.position0_, a.forceConstant_) == util::tie(b.position0_, b.forceConstant_);
    }
};

/*! \brief Type for two-three-center aggregates
 */
template<class TwoCenterType, class ThreeCenterType>
class ThreeCenterAggregate
{
public:
    using CarrierType            = ThreeCenterType;
    using TwoCenterAggregateType = TwoCenterType;

    ThreeCenterAggregate() = default;
    HOST_DEVICE_FUN ThreeCenterAggregate(const TwoCenterType& twoC, const ThreeCenterType& threeC) :
        twoC_(twoC), threeC_(threeC)
    {
    }

    HOST_DEVICE_FUN TwoCenterType& twoCenter() { return twoC_; }
    HOST_DEVICE_FUN const TwoCenterType& twoCenter() const { return twoC_; }

    HOST_DEVICE_FUN ThreeCenterType& carrier() { return threeC_; }
    HOST_DEVICE_FUN const ThreeCenterType& carrier() const { return threeC_; }

    enum Cargo : int
    {
        bond_ij  = 1 << 1,
        bond_jk  = 1 << 2,
        bond_jl  = 1 << 3,
        has_bond = bond_ij + bond_jk + bond_jl
    };

    int manifest = 0;

private:
    TwoCenterType   twoC_;
    ThreeCenterType threeC_;

    using Aggregate = ThreeCenterAggregate<TwoCenterType, ThreeCenterType>;

    friend HOST_DEVICE_FUN inline bool operator<(const Aggregate& a, const Aggregate& b)
    {
        return util::tie(a.twoC_, a.threeC_) < util::tie(b.twoC_, b.threeC_);
    }

    friend HOST_DEVICE_FUN inline bool operator==(const Aggregate& a, const Aggregate& b)
    {
        return util::tie(a.twoC_, a.threeC_) == util::tie(b.twoC_, b.threeC_);
    }
};

/*! \brief Type for two-three-four-center aggregates
 */
template<class TwoCenterType, class ThreeCenterType, class PairType, class FourCenterType>
class FourCenterAggregate
{
public:
    using CarrierType              = FourCenterType;
    using ThreeCenterAggregateType = ThreeCenterType;
    using TwoCenterAggregateType   = TwoCenterType;
    using PairAggregateType        = PairType;

    FourCenterAggregate() = default;

    HOST_DEVICE_FUN TwoCenterType& twoCenter() { return twoC_; }
    HOST_DEVICE_FUN const TwoCenterType& twoCenter() const { return twoC_; }

    HOST_DEVICE_FUN ThreeCenterType& threeCenter() { return threeC_; }
    HOST_DEVICE_FUN const ThreeCenterType& threeCenter() const { return threeC_; }

    HOST_DEVICE_FUN PairType& pair() { return pair_; }
    HOST_DEVICE_FUN const PairType& pair() const { return pair_; }

    HOST_DEVICE_FUN FourCenterType& carrier() { return fourC_; }
    HOST_DEVICE_FUN const FourCenterType& carrier() const { return fourC_; }

    enum Cargo : int
    {
        angle_j   = 1 << 1,
        angle_k   = 1 << 2,
        has_angle = angle_j + angle_k,
        bond_ij   = 1 << 3,
        bond_jk   = 1 << 4,
        bond_jl   = 1 << 5,
        has_bond  = bond_ij + bond_jk + bond_jl,
        pair_14   = 1 << 6
    };

    int manifest = 0;

private:
    TwoCenterType   twoC_;
    ThreeCenterType threeC_;
    PairType        pair_;
    FourCenterType  fourC_;

    using Aggregate = FourCenterAggregate<TwoCenterType, ThreeCenterType, PairType, FourCenterType>;

    friend HOST_DEVICE_FUN inline bool operator<(const Aggregate& a, const Aggregate& b)
    {
        return util::tie(a.twoC_, a.threeC_, a.pair_, a.fourC_, a.manifest)
               < util::tie(b.twoC_, b.threeC_, b.pair_, b.fourC_, b.manifest);
    }

    friend HOST_DEVICE_FUN inline bool operator==(const Aggregate& a, const Aggregate& b)
    {
        return util::tie(a.twoC_, a.threeC_, a.pair_, a.fourC_, a.manifest)
               == util::tie(b.twoC_, b.threeC_, b.pair_, b.fourC_, b.manifest);
    }
};


} // namespace nblib

#endif // NBLIB_LISTEDFORCES_BONDTYPES_H
