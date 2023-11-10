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
/*! \internal \file
 * \brief
 * Implements the data flow from ListedInteractionData and coordinates
 * down to the individual interaction type kernels
 *
 * Intended for internal use inside ListedCalculator only.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#ifndef NBLIB_LISTEDFORCES_DATAFLOW_HPP
#define NBLIB_LISTEDFORCES_DATAFLOW_HPP

#include "gromacs/utility/arrayref.h"

#include "nblib/listed_forces/definitions.h"
#include "nblib/listed_forces/kernels.hpp"
#include "nblib/listed_forces/traits.h"
#include "nblib/pbc.hpp"
#include "nblib/util/tuple.hpp"
#include "nblib/util/util.hpp"

namespace nblib
{

template<class TwoCenterType, class StackVector, class Lambda>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto computeTwoCenter(const TwoCenterType& parameterA,
                                                         const TwoCenterType& parameterB,
                                                         const StackVector&   dx,
                                                         Lambda               lambda,
                                                         StackVector*         fi,
                                                         StackVector*         fj)
{
    using ValueType = VectorValueType_t<StackVector>;

    ValueType dr2              = dot(dx, dx);
    ValueType dr               = std::sqrt(dr2);
    auto [force, energy, dvdl] = bondKernel(dr, parameterA, parameterB, lambda);

    // avoid division by 0
    if (dr2 != ValueType(0.0))
    {
        force /= dr;
        spreadTwoCenterForces(force, dx, fi, fj);
    }

    return util::make_tuple(energy, dvdl);
}

template<class TwoCenterType, class StackVector, class Lambda>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto computeTwoCenter(const TwoCenterType&           parametersA,
                                                         const TwoCenterType&           parametersB,
                                                         const StackVector&             dx,
                                                         VectorValueType_t<StackVector> q,
                                                         Lambda                         lambda,
                                                         StackVector*                   fi,
                                                         StackVector*                   fj)
{
    using ValueType = VectorValueType_t<StackVector>;

    ValueType dr2 = dot(dx, dx);
    ValueType dr  = std::sqrt(dr2);

    auto      forceAndEnergy = bondKernel(dr, q, parametersA, parametersB, lambda);
    ValueType force          = util::get<0>(forceAndEnergy);

    // avoid division by 0
    if (dr2 != ValueType(0.0))
    {
        force /= dr;
        spreadTwoCenterForces(force, dx, fi, fj);
    }

    return util::discardFirstElement(forceAndEnergy);
}


/*! \brief calculate charged pair interactions with XYZQ inputs
 *
 * @tparam TwoCenterType         interaction type to compute
 * @tparam MemVector             4-wide xyzq vector, DeviceTagged when called from GPU code
 * @tparam Lambda                real if computing FEP, NoFepLambdaType otherwise
 * @tparam Buffer                pointer to 3-wide vector buffer for force output
 * @tparam ShiftForce            pointer to 3-wide vectors if computing shift forces, dummy type otherwise
 * @tparam Pbc                   something that provides pbc_dx_aiuc
 * @param[in]    index           two interaction particle indices and the parameter index
 * @param[in]    bondInstancesA  pointer to interaction parameters for system state A
 * @param[in]    bondInstancesB  pointer to interaction parameters for system state B
 * @param[in]    xyzq            pointer particle input coordinates and charges
 * @param[in]    lambda          interpolation between A and B states
 * @param[inout] forces          force output to add to
 * @param[inout] shiftForces     shift force output to add to
 * @param[in]    pbc             pbc object to use for calculating distances
 * @return                       packed interaction energies (ePot, eVdW, eCoul, dVdl)
 */
template<class TwoCenterType, class MemVector, class Lambda, class Buffer, class ShiftForce, class Pbc, std::enable_if_t<HasCharge<TwoCenterType>{}, int> = 0>
HOST_DEVICE_FUN HOST_INLINE auto dispatchInteraction(IndexArray<3>        index,
                                                     const TwoCenterType* bondInstancesA,
                                                     const TwoCenterType* bondInstancesB,
                                                     const MemVector*     xyzq,
                                                     Lambda               lambda,
                                                     Buffer*              forces,
                                                     ShiftForce*          shiftForces,
                                                     const Pbc&           pbc)
{
    using ValueType = VectorValueType_t<MemVector>;
    using Vec       = StackVec3<ValueType>;
    KernelEnergy<ValueType> energy;

    int  i      = index[0];
    int  j      = index[1];
    auto xyzq_i = xyzq[i];
    auto xyzq_j = xyzq[j];

    TwoCenterType bondA = bondInstancesA[index[2]];
    // conditional load of B parameters only if Lambda is not NoFepLambdaType
    TwoCenterType bondB = loadInteractionParameters<Lambda>(bondInstancesB, index[2]);

    Vec fi{ 0, 0, 0 }, fj{ 0, 0, 0 };

    Vec dx;
    int sIdx = pbc.dxAiuc(xyzq_i, xyzq_j, dx);

    ValueType charge = (Contains<TwoCenterType, PairListedTypes>{}) ? xyzq_i[3] * xyzq_j[3] : xyzq_j[3];

    auto interactionEnergies = computeTwoCenter(bondA, bondB, dx, charge, lambda, &fi, &fj);
    addEnergy<TwoCenterType>(&energy, interactionEnergies);

    addForce(&(*forces)[i], fi);
    addForce(&(*forces)[j], fj);

    if (sIdx != gmx::c_centralShiftIndex)
    {
        addForce(shiftForces + sIdx, fi);
        addForce(shiftForces + gmx::c_centralShiftIndex, fj);
    }

    return energy;
}

/*! \brief calculate two-center interactions
 *
 * \tparam Force buffer type
 * \tparam TwoCenterType The bond type to compute; used for type deduction
 * \tparam Cartesian vector type
 * \tparam PBC type
 * \param[in] index The atom and parameter indices used for computing the interaction
 * \param[in] bondInstances The full type-specific interaction list
 * \param[in] x The coordinates
 * \param[in/out] forces The forces
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template<class TwoCenterType, class MemVector, class Lambda, class Buffer, class ShiftForce, class Pbc, std::enable_if_t<!HasCharge<TwoCenterType>{}, int> = 0>
HOST_DEVICE_FUN HOST_INLINE auto dispatchInteraction(IndexArray<3>        index,
                                                     const TwoCenterType* bondInstancesA,
                                                     const TwoCenterType* bondInstancesB,
                                                     const MemVector*     x,
                                                     Lambda               lambda,
                                                     Buffer*              forces,
                                                     ShiftForce*          shiftForces,
                                                     const Pbc&           pbc)
{
    using ValueType = VectorValueType_t<MemVector>;
    using Vec       = StackVec3<ValueType>;
    KernelEnergy<ValueType> energy;

    int i  = index[0];
    int j  = index[1];
    Vec xi = loadVec(x[i]);
    Vec xj = loadVec(x[j]);

    TwoCenterType bondA = bondInstancesA[index[2]];
    // conditional load of B parameters only if Lambda is not NoFepLambdaType
    TwoCenterType bondB = loadInteractionParameters<Lambda>(bondInstancesB, index[2]);

    Vec fi{ 0, 0, 0 }, fj{ 0, 0, 0 };

    Vec dx;
    int sIdx = pbc.dxAiuc(xi, xj, dx);

    auto [ePot, dvdl] = computeTwoCenter(bondA, bondB, dx, lambda, &fi, &fj);

    energy.carrier()              = ePot;
    energy.freeEnergyDerivative() = dvdl;

    addForce(&(*forces)[i], fi);
    addForce(&(*forces)[j], fj);

    if (sIdx != gmx::c_centralShiftIndex)
    {
        addForce(shiftForces + sIdx, fi);
        addForce(shiftForces + gmx::c_centralShiftIndex, fj);
    }

    return energy;
}


template<class ThreeCenterType, class StackVector, class Lambda, std::enable_if_t<HasTwoCenterAggregate<ThreeCenterType>::value, int> = 0>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto addTwoCenterAggregate(const ThreeCenterType& parametersA,
                                                              const ThreeCenterType& parametersB,
                                                              const StackVector&     rij,
                                                              const StackVector&     rkj,
                                                              Lambda                 lambda,
                                                              StackVector*           fi,
                                                              StackVector*           fj,
                                                              StackVector*           fk)
{
    using ValueType = VectorValueType_t<StackVector>;

    if (parametersA.manifest & ThreeCenterType::bond_ij)
    {
        // i-j bond
        return computeTwoCenter(parametersA.twoCenter(), parametersB.twoCenter(), rij, lambda, fi, fj);
    }
    if (parametersA.manifest & ThreeCenterType::bond_jk)
    {
        // j-k bond
        return computeTwoCenter(parametersA.twoCenter(), parametersB.twoCenter(), rkj, lambda, fk, fj);
    }

    // aggregate is empty
    return util::make_tuple(ValueType(0), ValueType(0));
};

template<class ThreeCenterType, class StackVector, class Lambda, std::enable_if_t<!HasTwoCenterAggregate<ThreeCenterType>::value, int> = 0>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto addTwoCenterAggregate(const ThreeCenterType& /* parametersA */,
                                                              const ThreeCenterType& /* parametersB */,
                                                              const StackVector& /* rij */,
                                                              const StackVector& /* rkj */,
                                                              Lambda /*lambda*/,
                                                              StackVector* /* fi */,
                                                              StackVector* /* fj */,
                                                              StackVector* /* fk */)
{
    using ValueType = VectorValueType_t<StackVector>;
    return util::make_tuple(ValueType(0), ValueType(0));
}

template<class ThreeCenterType, class StackVector, class Lambda>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto computeThreeCenter(const ThreeCenterType& parametersA,
                                                           const ThreeCenterType& parametersB,
                                                           const StackVector&     rij,
                                                           const StackVector&     rkj,
                                                           const StackVector& /* rik */,
                                                           Lambda       lambda,
                                                           StackVector* fi,
                                                           StackVector* fj,
                                                           StackVector* fk)
{
    using ValueType = VectorValueType_t<StackVector>;
    // calculate 3-center common quantities: angle between x1-x2 and x2-x3
    ValueType costh = basicVectorCosAngle(rij, rkj);
    ValueType theta = std::acos(costh);

    // call type-specific angle kernel, e.g. harmonic, restricted, quartic, etc.
    auto [force, energy, dvdl] = threeCenterKernel(theta, parametersA, parametersB, lambda);

    spreadThreeCenterForces(costh, force, rij, rkj, fi, fj, fk);

    return util::make_tuple(energy, dvdl);
}

template<class StackVector, class Lambda>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto computeThreeCenter(const LinearAngle& parametersA,
                                                           const LinearAngle& parametersB,
                                                           const StackVector& rij,
                                                           const StackVector& rkj,
                                                           const StackVector& /* rik */,
                                                           Lambda       lambda,
                                                           StackVector* fi,
                                                           StackVector* fj,
                                                           StackVector* fk)
{
    using ValueType = VectorValueType_t<StackVector>;

    ValueType b      = parametersA.equilConstant() - ValueType(1);
    auto      dr_vec = b * rkj - parametersA.equilConstant() * rij;
    ValueType dr     = norm(dr_vec);

    auto [ci, ck, energy, dvdl] = threeCenterKernel(dr, parametersA, parametersB, lambda);

    StackVector fi_loc = ci * dr_vec;
    StackVector fk_loc = ck * dr_vec;
    StackVector fj_loc = -(fi_loc + fk_loc);
    *fi += fi_loc;
    *fj += fj_loc;
    *fk += fk_loc;

    return util::make_tuple(energy, dvdl);
}

template<class StackVector, class Lambda>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto computeThreeCenter(const CrossBondBond& parametersA,
                                                           const CrossBondBond& /*parametersB*/,
                                                           const StackVector& rij,
                                                           const StackVector& rkj,
                                                           const StackVector& /* rik */,
                                                           Lambda /*lambda*/,
                                                           StackVector* fi,
                                                           StackVector* fj,
                                                           StackVector* fk)
{
    using ValueType       = VectorValueType_t<StackVector>;
    auto [ci, ck, energy] = threeCenterKernel(norm(rij), norm(rkj), parametersA);

    StackVector fi_loc = ci * rij;
    StackVector fk_loc = ck * rkj;
    StackVector fj_loc = -(fi_loc + fk_loc);
    *fi += fi_loc;
    *fj += fj_loc;
    *fk += fk_loc;

    return util::make_tuple(energy, ValueType(0));
}

template<class StackVector, class Lambda>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto computeThreeCenter(const CrossBondAngle& parameters,
                                                           const CrossBondAngle& /*parametersB*/,
                                                           const StackVector& rij,
                                                           const StackVector& rkj,
                                                           const StackVector& rik,
                                                           Lambda /*lambda*/,
                                                           StackVector* fi,
                                                           StackVector* fj,
                                                           StackVector* fk)
{
    using ValueType           = VectorValueType_t<StackVector>;
    auto [ci, cj, ck, energy] = threeCenterKernel(norm(rij), norm(rkj), norm(rik), parameters);

    StackVector fi_loc = ci * rij + ck * rik;
    StackVector fk_loc = cj * rkj - ck * rik;
    StackVector fj_loc = -(fi_loc + fk_loc);
    *fi += fi_loc;
    *fj += fj_loc;
    *fk += fk_loc;

    return util::make_tuple(energy, ValueType(0));
}

/*! \brief Calculate three-center interactions
 *
 * \tparam Force buffer type
 * \tparam Three centre interaction parameters
 * \tparam Cartesian vector type
 * \tparam PBC type
 * \param[in] index
 * \param[in] Bond parameters
 * \param[in] x coordinate array
 * \param[in/out] Force buffer
 * \param[in] PBC
 * \return Computed kernel energies
 */
template<class Buffer, class ThreeCenterType, class MemVector, class Lambda, class ShiftForce, class Pbc>
HOST_DEVICE_FUN HOST_INLINE auto dispatchInteraction(IndexArray<4>          index,
                                                     const ThreeCenterType* parametersA,
                                                     const ThreeCenterType* parametersB,
                                                     const MemVector*       x,
                                                     Lambda                 lambda,
                                                     Buffer*                forces,
                                                     ShiftForce*            shiftForces,
                                                     const Pbc&             pbc)
{
    using ValueType = VectorValueType_t<MemVector>;
    using Vec       = StackVec3<ValueType>;

    KernelEnergy<ValueType> energy;

    int i  = index[0];
    int j  = index[1];
    int k  = index[2];
    Vec xi = loadVec(x[i]);
    Vec xj = loadVec(x[j]);
    Vec xk = loadVec(x[k]);

    ThreeCenterType parmA = parametersA[index[3]];
    // conditional load of B parameters only if Lambda is not NoFepLambdaType
    ThreeCenterType parmB = loadInteractionParameters<Lambda>(parametersB, index[3]);

    Vec fi{ 0, 0, 0 }, fj{ 0, 0, 0 }, fk{ 0, 0, 0 };

    Vec rij, rkj, rik;
    int sIdx_ij = pbc.dxAiuc(xi, xj, rij);
    int sIdx_kj = pbc.dxAiuc(xk, xj, rkj);
    pbc.dxAiuc(xi, xk, rik);

    auto [e3c, dvdl3c] = computeThreeCenter(parmA, parmB, rij, rkj, rik, lambda, &fi, &fj, &fk);
    auto [e2c, dvdl2c] = addTwoCenterAggregate(parmA, parmB, rij, rkj, lambda, &fi, &fj, &fk);

    energy.carrier()              = e3c;
    energy.twoCenterAggregate()   = e2c;
    energy.freeEnergyDerivative() = dvdl3c + dvdl2c;

    addForce(&(*forces)[i], fi);
    addForce(&(*forces)[j], fj);
    addForce(&(*forces)[k], fk);

    if (sIdx_ij != gmx::c_centralShiftIndex || sIdx_kj != gmx::c_centralShiftIndex)
    {
        addForce(shiftForces + sIdx_ij, fi);
        addForce(shiftForces + gmx::c_centralShiftIndex, fj);
        addForce(shiftForces + sIdx_kj, fk);
    }

    return energy;
}

template<class FourCenterType, class StackVector, class Lambda, std::enable_if_t<HasThreeCenterAggregate<FourCenterType>::value, int> = 0>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto addThreeCenterAggregate(const FourCenterType& parameterA,
                                                                const FourCenterType& parameterB,
                                                                const StackVector&    rij,
                                                                const StackVector&    rkj,
                                                                const StackVector&    rkl,
                                                                Lambda                lambda,
                                                                StackVector*          fi,
                                                                StackVector*          fj,
                                                                StackVector*          fk,
                                                                StackVector*          fl)
{
    using ValueType = VectorValueType_t<StackVector>;
    if (parameterA.manifest & FourCenterType::angle_j)
    {
        return computeThreeCenter(
                parameterA.threeCenter(), parameterB.threeCenter(), rij, rkj, StackVector{}, lambda, fi, fj, fk);
    }
    if (parameterA.manifest & FourCenterType::angle_k)
    {
        return computeThreeCenter(
                parameterA.threeCenter(), parameterB.threeCenter(), -rkj, -rkl, StackVector{}, lambda, fj, fk, fl);
    }

    return util::make_tuple(ValueType(0), ValueType(0));
}

template<class FourCenterType, class StackVector, class Lambda, std::enable_if_t<!HasThreeCenterAggregate<FourCenterType>::value, int> = 0>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto addThreeCenterAggregate(const FourCenterType& /* parameterA*/,
                                                                const FourCenterType& /*parameterB*/,
                                                                const StackVector& /* rij */,
                                                                const StackVector& /* rkj */,
                                                                const StackVector& /* rkl */,
                                                                Lambda /*lambda*/,
                                                                StackVector* /* fi */,
                                                                StackVector* /* fj */,
                                                                StackVector* /* fk */,
                                                                StackVector* /* fl */)
{
    using ValueType = VectorValueType_t<StackVector>;
    return util::make_tuple(ValueType(0), ValueType(0));
}

template<class FourCenterType, class StackVector, class Lambda, std::enable_if_t<HasTwoCenterAggregate<FourCenterType>::value, int> = 0>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto addTwoCenterAggregate(const FourCenterType& parameterA,
                                                              const FourCenterType& parameterB,
                                                              const StackVector&    rij,
                                                              const StackVector&    rkj,
                                                              const StackVector&    rlj,
                                                              Lambda                lambda,
                                                              StackVector*          fi,
                                                              StackVector*          fj,
                                                              StackVector*          fk,
                                                              StackVector*          fl)
{
    using ValueType = VectorValueType_t<StackVector>;
    if (parameterA.manifest & FourCenterType::bond_ij)
    {
        return computeTwoCenter(parameterA.twoCenter(), parameterB.twoCenter(), rij, lambda, fi, fj);
    }
    if (parameterA.manifest & FourCenterType::bond_jk)
    {
        return computeTwoCenter(parameterA.twoCenter(), parameterB.twoCenter(), rkj, lambda, fk, fj);
    }
    if (parameterA.manifest & FourCenterType::bond_jl)
    {
        return computeTwoCenter(parameterA.twoCenter(), parameterB.twoCenter(), rlj, lambda, fl, fj);
    }

    return util::make_tuple(ValueType(0), ValueType(0));
}

template<class FourCenterType, class StackVector, class Lambda, std::enable_if_t<!HasTwoCenterAggregate<FourCenterType>::value, int> = 0>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto addTwoCenterAggregate(const FourCenterType& /* parameterA*/,
                                                              const FourCenterType& /*parameterB*/,
                                                              const StackVector& /* rij */,
                                                              const StackVector& /* rkj */,
                                                              const StackVector& /* rlj */,
                                                              Lambda /*lambda*/,
                                                              StackVector* /* fi */,
                                                              StackVector* /* fj */,
                                                              StackVector* /* fk */,
                                                              StackVector* /* fl */)
{
    using ValueType = VectorValueType_t<StackVector>;
    return util::make_tuple(ValueType(0), ValueType(0));
}

template<class FourCenterType, class StackVector, class ForceVector, class Lambda, class Pbc, std::enable_if_t<HasPairAggregate<FourCenterType>{}, int> = 0>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto addPairAggregate(const FourCenterType&          parameterA,
                                                         const FourCenterType&          parameterB,
                                                         const StackVector&             xi,
                                                         const StackVector&             xl,
                                                         VectorValueType_t<StackVector> qi,
                                                         VectorValueType_t<StackVector> ql,
                                                         Lambda                         lambda,
                                                         ForceVector*                   fi,
                                                         ForceVector*                   fl,
                                                         const Pbc&                     pbc)
{
    using ValueType = VectorValueType_t<StackVector>;

    if (parameterA.manifest & FourCenterType::pair_14)
    {
        util::array<ValueType, 3> dxIL;
        pbc.dxAiuc(xi, xl, dxIL);
        ValueType charge = (StackVector{}.size() == 4) ? xi[3] * xl[3] : qi * ql;
        return computeTwoCenter(parameterA.pair(), parameterB.pair(), dxIL, charge, lambda, fi, fl);
    }
    return util::make_tuple(ValueType(0), ValueType(0), ValueType(0));
}

template<class FourCenterType, class StackVector, class ForceVector, class Lambda, class Pbc, std::enable_if_t<!HasPairAggregate<FourCenterType>{}, int> = 0>
HOST_DEVICE_FUN HOST_DEVICE_INLINE auto addPairAggregate(const FourCenterType& /*parameterA*/,
                                                         const FourCenterType& /*parameterB*/,
                                                         const StackVector& /*xi*/,
                                                         const StackVector& /*xl*/,
                                                         VectorValueType_t<StackVector> /*qi*/,
                                                         VectorValueType_t<StackVector> /*ql*/,
                                                         Lambda /*lambda*/,
                                                         ForceVector* /*fi*/,
                                                         ForceVector* /*fl*/,
                                                         const Pbc& /*pbc*/)
{
    using ValueType = VectorValueType_t<StackVector>;

    return util::make_tuple(ValueType(0), ValueType(0), ValueType(0));
}

/*! \brief Calculate four-center interactions
 *
 * \tparam Force buffer type
 * \tparam FourCenterType The bond type to compute; used for type deduction
 * \tparam Cartesian vector type
 * \tparam PBC type
 * \param[in] index The atom and parameter indices used for computing the interaction
 * \param[in] parametersA The full type-specific interaction list
 * \param[in] x The coordinates
 * \param[in/out] forces The forces
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template<class Buffer, class FourCenterType, class MemVector, class Lambda, class ShiftForce, class Pbc>
HOST_DEVICE_FUN HOST_INLINE auto dispatchInteraction(IndexArray<5>         index,
                                                     const FourCenterType* parametersA,
                                                     const FourCenterType* parametersB,
                                                     const MemVector*      x,
                                                     Lambda                lambda,
                                                     Buffer*               forces,
                                                     ShiftForce*           shiftForces,
                                                     const Pbc&            pbc)
{
    using ValueType = VectorValueType_t<MemVector>;
    using Vec       = StackVec3<ValueType>;
    KernelEnergy<ValueType> energy;

    int i = index[0];
    int j = index[1];
    int k = index[2];
    int l = index[3];

    auto xi = x[i];
    auto xj = x[j];
    auto xk = x[k];
    auto xl = x[l];

    Vec fi{ 0, 0, 0 }, fj{ 0, 0, 0 }, fk{ 0, 0, 0 }, fl{ 0, 0, 0 };

    Vec dxIJ, dxKJ, dxKL, dxLJ;
    int sIdx_ij = pbc.dxAiuc(xi, xj, dxIJ);
    int sIdx_kj = pbc.dxAiuc(xk, xj, dxKJ);
    pbc.dxAiuc(xk, xl, dxKL);

    int sIdx_lj = pbc.dxAiuc(xl, xj, dxLJ);

    FourCenterType paramsA = parametersA[index[4]];
    // conditional load of B parameters only if Lambda is not NoFepLambdaType
    FourCenterType paramsB = loadInteractionParameters<Lambda>(parametersB, index[4]);

    Vec       m, n;
    ValueType phi = dihedralPhi(dxIJ, dxKJ, dxKL, &m, &n);

    auto [force, kernelEnergy, dvdl] = fourCenterKernel(phi, paramsA, paramsB, lambda);

    auto [e3c, dvdl3c] =
            addThreeCenterAggregate(paramsA, paramsB, dxIJ, dxKJ, dxKL, lambda, &fi, &fj, &fk, &fl);
    auto [e2c, dvdl2c] =
            addTwoCenterAggregate(paramsA, paramsB, dxIJ, dxKJ, dxLJ, lambda, &fi, &fj, &fk, &fl);
    auto [eVdw, eCoul, dvdlPair] =
            addPairAggregate(paramsA, paramsB, xi, xl, 0, 0, lambda, &fi, &fl, pbc);

    energy.carrier()              = kernelEnergy;
    energy.twoCenterAggregate()   = e2c;
    energy.threeCenterAggregate() = e3c;
    energy.eVdw()                 = eVdw;
    energy.eCoul()                = eCoul;
    energy.freeEnergyDerivative() = dvdl + dvdl3c + dvdl2c + dvdlPair;

    spreadFourCenterForces(force, dxIJ, dxKJ, dxKL, m, n, &fi, &fj, &fk, &fl);

    addForce(&(*forces)[i], fi);
    addForce(&(*forces)[j], fj);
    addForce(&(*forces)[k], fk);
    addForce(&(*forces)[l], fl);

    if (sIdx_ij != gmx::c_centralShiftIndex || sIdx_kj != gmx::c_centralShiftIndex
        || sIdx_lj != gmx::c_centralShiftIndex)
    {
        addForce(shiftForces + sIdx_ij, fi);
        addForce(shiftForces + gmx::c_centralShiftIndex, fj);
        addForce(shiftForces + sIdx_kj, fk);
        addForce(shiftForces + sIdx_lj, fl);
    }

    return energy;
}

/*! \brief Calculate five-center interactions
 *
 * \tparam Force buffer type
 * \tparam FiveCenterType The bond type to compute; used for type deduction
 * \tparam Cartesian vector type
 * \tparam PBC type
 * \param[in] index The atom and parameter indices used for computing the interaction
 * \param[in] parametersA The full type-specific interaction list
 * \param[in] x The coordinates
 * \param[in/out] forces The forces
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template<class Buffer, class FiveCenterType, class MemVector, class Lambda, class ShiftForce, class Pbc>
HOST_DEVICE_FUN HOST_INLINE auto dispatchInteraction(IndexArray<6>         index,
                                                     const FiveCenterType* parametersA,
                                                     const FiveCenterType* parametersB,
                                                     const MemVector*      x,
                                                     Lambda                lambda,
                                                     Buffer*               forces,
                                                     ShiftForce*           shiftForces,
                                                     const Pbc&            pbc)
{
    using ValueType = VectorValueType_t<MemVector>;
    using Vec       = StackVec3<ValueType>;
    KernelEnergy<ValueType> energy;

    int i = index[0];
    int j = index[1];
    int k = index[2];
    int l = index[3];
    int m = index[4];

    auto xi = x[i];
    auto xj = x[j];
    auto xk = x[k];
    auto xl = x[l];
    auto xm = x[m];

    Vec dxIJ, dxJK, dxKL, dxLM;
    pbc.dxAiuc(xi, xj, dxIJ);
    pbc.dxAiuc(xj, xk, dxJK);
    pbc.dxAiuc(xk, xl, dxKL);
    pbc.dxAiuc(xl, xm, dxLM);

    FiveCenterType paramsA = parametersA[index[5]];
    // conditional load of B parameters only if Lambda is not NoFepLambdaType
    FiveCenterType paramsB = loadInteractionParameters<Lambda>(parametersB, index[5]);

    // this dispatch function is not in use yet, because CMap is not yet implemented
    // we don't want to add [[maybe_unused]] in the signature since the params will
    // be used once CMap is implemented, and we also don't want compiler warnings,
    // so we cast to void.
    (void)paramsA;
    (void)paramsB;
    (void)lambda;
    (void)forces;
    (void)shiftForces;

    return energy;
}


/*! \brief implement a loop over bonds for a given BondType and Kernel
 *  corresponds to e.g. the "bonds" function at Gromacs:bonded.cpp@450
 *
 * \param[in] indices interaction atom pair indices + bond parameter index
 * \param[in] interactionParameters bond/interaction parameters
 * \param[in] x coordinate input
 * \param[in/out] forces The forces
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template<class Index, class InteractionType, class MemVector, class Buffer, class Lambda, class ShiftForce, class Pbc>
auto computeForces(gmx::ArrayRef<const Index>           indices,
                   gmx::ArrayRef<const InteractionType> parametersA,
                   gmx::ArrayRef<const InteractionType> parametersB,
                   gmx::ArrayRef<const MemVector>       x,
                   Lambda                               lambda,
                   Buffer*                              forces,
                   gmx::ArrayRef<ShiftForce>            shiftForces,
                   const Pbc&                           pbc)
{
    KernelEnergy<VectorValueType_t<MemVector>> energy;

    for (const auto& index : indices)
    {
        energy += dispatchInteraction(
                index, parametersA.data(), parametersB.data(), x.data(), lambda, forces, shiftForces.data(), pbc);
    }

    return energy;
}

//! \brief convenience overload without shift forces
template<class Index, class InteractionType, class MemVector, class Buffer, class Pbc>
auto computeForces(gmx::ArrayRef<const Index>           indices,
                   gmx::ArrayRef<const InteractionType> parameters,
                   gmx::ArrayRef<const MemVector>       x,
                   Buffer*                              forces,
                   const Pbc&                           pbc)
{
    return computeForces(
            indices, parameters, parameters, x, NoFepLambdaType{}, forces, gmx::ArrayRef<std::nullptr_t>{}, pbc);
}

/*! \brief implement a loop over bond types and accumulate their force contributions
 *
 * \param[in] interactions interaction pairs and bond parameters
 * \param[in] x coordinate input
 * \param[in/out] forces output force buffer
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template<class MemVector, class Buffer, class ShiftForce, class Pbc>
auto reduceListedForces(const ListedInteractionData&   interactions,
                        gmx::ArrayRef<const MemVector> x,
                        Buffer*                        forces,
                        gmx::ArrayRef<ShiftForce>      shiftForces,
                        const Pbc&                     pbc)
{
    using ValueType = VectorValueType_t<MemVector>;
    ListedEnergies energies{ 0 };
    std::fill(energies.begin(), energies.end(), 0);

    // calculate one bond type
    auto computeForceType = [forces, x, shiftForces, &energies, &pbc](const auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;

        gmx::ArrayRef<const InteractionIndex<InteractionType>> indices(interactionElement.indices);
        gmx::ArrayRef<const InteractionType> parametersA(interactionElement.parametersA);
        gmx::ArrayRef<const InteractionType> parametersB(interactionElement.parametersB);

        KernelEnergy<ValueType> energy = computeForces(
                indices, parametersA, parametersB, x, NoFepLambdaType{}, forces, shiftForces, pbc);

        energies[CarrierIndex<InteractionType>{}] += energy.carrier();
        energies[TwoCenterAggregateIndex<InteractionType>{}] += energy.twoCenterAggregate();
        energies[ThreeCenterAggregateIndex<InteractionType>{}] += energy.threeCenterAggregate();
        energies[FepIndex{}] += energy.freeEnergyDerivative();
    };

    auto computeIndices = subsetIndices(BasicListedTypes{}, AllListedTypes{});
    // calculate all bond types, returns a tuple with the energies for each type
    for_each_tuple(computeForceType, tieElements(interactions, computeIndices));

    return energies;
}


//! \brief convenience overload without shift forces
template<class MemVector, class Buffer, class Pbc>
auto reduceListedForces(const ListedInteractionData&   interactions,
                        gmx::ArrayRef<const MemVector> x,
                        Buffer*                        forces,
                        const Pbc&                     pbc)
{
    return reduceListedForces(interactions, x, forces, gmx::ArrayRef<std::nullptr_t>{}, pbc);
}

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_DATAFLOW_HPP
