/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

#include <tuple>
#include <vector>

#include "nblib/listed_forces/traits.h"
#include "nblib/listed_forces/kernels.hpp"
#include "nblib/util/util.hpp"
#include "nblib/pbc.hpp"
#include "nblib/vector.h"
#include "gromacs/utility/arrayref.h"

#define NBLIB_ALWAYS_INLINE __attribute((always_inline))

namespace nblib
{

//! \brief returns the address of an element in the shiftForces buffer
template<class ShiftForce>
inline ShiftForce* accessShiftForces(int shiftIndex, gmx::ArrayRef<ShiftForce> shiftForces)
{
    return shiftForces.data() + shiftIndex;
}

//! \brief dummy op in case shift forces are not computed (will be removed by the compiler)
inline std::nullptr_t accessShiftForces(int /* shiftIndex */,
                                        gmx::ArrayRef<std::nullptr_t> /* shiftForces */)
{
    return nullptr;
}

template<class TwoCenterType, class BasicVector, class ShiftForce>
inline NBLIB_ALWAYS_INLINE
auto computeTwoCenter(const TwoCenterType& parameters, const BasicVector& dx, BasicVector* fi, BasicVector* fj,
                      ShiftForce* sh_f, ShiftForce* sh_fc)
{
    using ValueType = BasicVectorValueType_t<BasicVector>;

    ValueType dr2 = dot(dx, dx);
    ValueType dr  = std::sqrt(dr2);
    auto [force, energy] = bondKernel(dr, parameters);

    // avoid division by 0
    if (dr2 != 0.0)
    {
        force /= dr;
        spreadTwoCenterForces(force, dx, fi, fj, sh_f, sh_fc);
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
template <class Buffer, class TwoCenterType, class BasicVector, class ShiftForce, class Pbc,
          std::enable_if_t<Contains<TwoCenterType, SupportedTwoCenterTypes>{}>* = nullptr>
inline NBLIB_ALWAYS_INLINE
auto dispatchInteraction(InteractionIndex<TwoCenterType> index,
                         gmx::ArrayRef<const TwoCenterType> bondInstances,
                         gmx::ArrayRef<const BasicVector> x,
                         Buffer* forces,
                         gmx::ArrayRef<ShiftForce> shiftForces,
                         const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<BasicVector>> energy;

    int i = std::get<0>(index);
    int j = std::get<1>(index);
    BasicVector   xi   = x[i];
    BasicVector   xj   = x[j];
    TwoCenterType bond = bondInstances[std::get<2>(index)];

    BasicVector dx;
    // calculate xi - xj modulo pbc
    int sIdx = pbc.dxAiuc(xi, xj, dx);

    ShiftForce* sh_f  = accessShiftForces(sIdx, shiftForces);
    ShiftForce* sh_fc = accessShiftForces(gmx::c_centralShiftIndex, shiftForces);

    energy.carrier() = computeTwoCenter(bond, dx, &(*forces)[i], &(*forces)[j], sh_f, sh_fc);
    return energy;
}


template<class ThreeCenterType, class BasicVector, class ShiftForce>
inline NBLIB_ALWAYS_INLINE
std::enable_if_t<HasTwoCenterAggregate<ThreeCenterType>::value, BasicVectorValueType_t<BasicVector>>
addTwoCenterAggregate(const ThreeCenterType& parameters, const BasicVector& rij, const BasicVector& rkj,
                      BasicVector* fi, BasicVector* fj, BasicVector* fk,
                      ShiftForce* shf_ij, ShiftForce* shf_kj, ShiftForce* shf_c)
{
if (parameters.manifest == ThreeCenterType::Cargo::ij)
    {
        // i-j bond
        return computeTwoCenter(parameters.twoCenter(), rij, fi, fj, shf_ij, shf_c);
    }
    if (parameters.manifest == ThreeCenterType::Cargo::jk)
    {
        // j-k bond
        return computeTwoCenter(parameters.twoCenter(), rkj, fk, fj, shf_kj, shf_c);
    }

    // aggregate is empty
    return 0.0;
};

template<class ThreeCenterType, class BasicVector, class ShiftForce>
inline NBLIB_ALWAYS_INLINE
std::enable_if_t<!HasTwoCenterAggregate<ThreeCenterType>::value, BasicVectorValueType_t<BasicVector>>
addTwoCenterAggregate(const ThreeCenterType& /* parameters */,
                      const BasicVector& /* rij */,
                      const BasicVector& /* rkj */,
                      BasicVector* /* fi */,
                      BasicVector* /* fj */,
                      BasicVector* /* fk */,
                      ShiftForce* /* shf_ij */,
                      ShiftForce* /* shf_kj */,
                      ShiftForce* /* shf_c */)
{
    return 0.0;
};

template<class ThreeCenterType, class BasicVector, class ShiftForce>
inline NBLIB_ALWAYS_INLINE
auto computeThreeCenter(const ThreeCenterType& parameters, const BasicVector& rij, const BasicVector& rkj,
                        const BasicVector& /* rik */, BasicVector* fi, BasicVector* fj, BasicVector* fk,
                        ShiftForce* shf_ij, ShiftForce* shf_kj, ShiftForce* shf_c)
{
    using ValueType = BasicVectorValueType_t<BasicVector>;
    // calculate 3-center common quantities: angle between x1-x2 and x2-x3
    // Todo: after sufficient evaluation, switch over to atan2 based algorithm
    ValueType costh = basicVectorCosAngle(rij, rkj); /* 25 */
    ValueType theta = std::acos(costh);    /* 10 */

    // call type-specific angle kernel, e.g. harmonic, restricted, quartic, etc.
    auto [force, energy] = threeCenterKernel(theta, parameters);

    spreadThreeCenterForces(costh, force, rij, rkj, fi, fj, fk, shf_ij, shf_kj, shf_c);

    return energy;
}

template<class BasicVector, class ShiftForce>
inline NBLIB_ALWAYS_INLINE
auto computeThreeCenter(const LinearAngle& parameters, const BasicVector& rij, const BasicVector& rkj,
                        const BasicVector& /* rik */, BasicVector* fi, BasicVector* fj, BasicVector* fk,
                        ShiftForce* shf_ij, ShiftForce* shf_kj, ShiftForce* shf_c)
{
    using ValueType = BasicVectorValueType_t<BasicVector>;

    ValueType b  = parameters.equilConstant() - 1;
    auto dr_vec  = b * rkj - parameters.equilConstant() * rij;
    ValueType dr = norm(dr_vec);

    auto [ci, ck, energy] = threeCenterKernel(dr, parameters);

    BasicVector fi_loc = ci * dr_vec;
    BasicVector fk_loc = ck * dr_vec;
    BasicVector fj_loc = ValueType(-1.0) * (fi_loc + fk_loc);
    *fi += fi_loc;
    *fj += fj_loc;
    *fk += fk_loc;

    addShiftForce(fi_loc, shf_ij);
    addShiftForce(fj_loc, shf_c);
    addShiftForce(fk_loc, shf_kj);

    return energy;
}

template<class BasicVector, class ShiftForce>
inline NBLIB_ALWAYS_INLINE
auto computeThreeCenter(const CrossBondBond& parameters, const BasicVector& rij, const BasicVector& rkj,
                        const BasicVector& /* rik */, BasicVector* fi, BasicVector* fj, BasicVector* fk,
                        ShiftForce* shf_ij, ShiftForce* shf_kj, ShiftForce* shf_c)
{
    using ValueType = BasicVectorValueType_t<BasicVector>;
    // 28 flops from the norm() calls
    auto [ci, ck, energy] = threeCenterKernel(norm(rij), norm(rkj), parameters);

    BasicVector fi_loc = ci * rij;
    BasicVector fk_loc = ck * rkj;
    BasicVector fj_loc = ValueType(-1.0) * (fi_loc + fk_loc);
    *fi += fi_loc;
    *fj += fj_loc;
    *fk += fk_loc;

    addShiftForce(fi_loc, shf_ij);
    addShiftForce(fj_loc, shf_c);
    addShiftForce(fk_loc, shf_kj);

    return energy;
}

template<class BasicVector, class ShiftForce>
inline NBLIB_ALWAYS_INLINE
auto computeThreeCenter(const CrossBondAngle& parameters, const BasicVector& rij, const BasicVector& rkj,
                        const BasicVector& rik, BasicVector* fi, BasicVector* fj, BasicVector* fk,
                        ShiftForce* shf_ij, ShiftForce* shf_kj, ShiftForce* shf_c)
{
    using ValueType = BasicVectorValueType_t<BasicVector>;
    // 42 flops from the norm() calls
    auto [ci, cj, ck, energy] = threeCenterKernel(norm(rij), norm(rkj), norm(rik), parameters);

    BasicVector fi_loc = ci * rij + ck * rik;
    BasicVector fk_loc = cj * rkj - ck * rik;
    BasicVector fj_loc = ValueType(-1.0) * (fi_loc + fk_loc);
    *fi += fi_loc;
    *fj += fj_loc;
    *fk += fk_loc;

    addShiftForce(fi_loc, shf_ij);
    addShiftForce(fj_loc, shf_c);
    addShiftForce(fk_loc, shf_kj);

    return energy;
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
template <class Buffer, class ThreeCenterType, class BasicVector, class ShiftForce, class Pbc,
          std::enable_if_t<Contains<ThreeCenterType, SupportedThreeCenterTypes>{}>* = nullptr>
inline NBLIB_ALWAYS_INLINE
auto dispatchInteraction(InteractionIndex<ThreeCenterType> index,
                         gmx::ArrayRef<const ThreeCenterType> parameters,
                         gmx::ArrayRef<const BasicVector> x,
                         Buffer* forces,
                         gmx::ArrayRef<ShiftForce> shiftForces,
                         const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<BasicVector>> energy;

    //! fetch input data: position vectors x1-x3 and interaction parameters
    int i = std::get<0>(index);
    int j = std::get<1>(index);
    int k = std::get<2>(index);
    BasicVector xi = x[i];
    BasicVector xj = x[j];
    BasicVector xk = x[k];
    ThreeCenterType threeCenterParameters = parameters[std::get<3>(index)];

    BasicVector fi{0,0,0}, fj{0,0,0}, fk{0,0,0};

    BasicVector rij, rkj, rik;
    int sIdx_ij = pbc.dxAiuc(xi, xj, rij); /* 3 */
    int sIdx_kj = pbc.dxAiuc(xk, xj, rkj); /* 3 */
    pbc.dxAiuc(xi, xk, rik); /* 3 */

    ShiftForce* shf_ij = accessShiftForces(sIdx_ij, shiftForces);
    ShiftForce* shf_kj = accessShiftForces(sIdx_kj, shiftForces);
    ShiftForce* shf_c  = accessShiftForces(gmx::c_centralShiftIndex, shiftForces);

    energy.carrier()            = computeThreeCenter(threeCenterParameters, rij, rkj, rik, &fi, &fj, &fk, shf_ij, shf_kj, shf_c);
    energy.twoCenterAggregate() = addTwoCenterAggregate(threeCenterParameters, rij, rkj, &fi, &fj, &fk, shf_ij, shf_kj, shf_c);

    (*forces)[i] += fi;
    (*forces)[j] += fj;
    (*forces)[k] += fk;

    return energy;
}

template<class FourCenterType, class BasicVector>
inline NBLIB_ALWAYS_INLINE
std::enable_if_t<HasThreeCenterAggregate<FourCenterType>::value, BasicVectorValueType_t<BasicVector>>
addThreeCenterAggregate(const FourCenterType& parameters,
                        const BasicVector& rij, const BasicVector& rkj, const BasicVector& rkl,
                        BasicVector* fi, BasicVector* fj, BasicVector* fk, BasicVector* fl)
{
    using ValueType = BasicVectorValueType_t<BasicVector>;
    if (parameters.manifest == FourCenterType::Cargo::j)
    {
        return computeThreeCenter(parameters.threeCenter(), rij, rkj, fi, fj, fk);
    }
    if (parameters.manifest == FourCenterType::Cargo::k)
    {
        return computeThreeCenter(parameters.threeCenter(), ValueType(-1.0)*rkj, ValueType(-1.0)*rkl, fj, fk, fl);
        //return computeThreeCenter(parameters.threeCenter(), rkj, rkl, fj, fk, fl);
    }

    // aggregate is empty
    return 0.0;
};

template<class FourCenterType, class BasicVector>
inline NBLIB_ALWAYS_INLINE
std::enable_if_t<!HasThreeCenterAggregate<FourCenterType>::value, BasicVectorValueType_t<BasicVector>>
addThreeCenterAggregate(const FourCenterType& /* parameters*/,
                        const BasicVector& /* rij */,
                        const BasicVector& /* rkj */,
                        const BasicVector& /* rkl */,
                        BasicVector* /* fi */,
                        BasicVector* /* fj */,
                        BasicVector* /* fk */,
                        BasicVector* /* fl */)
{
    return 0.0;
};

/*! \brief Calculate four-center interactions
 *
 * \tparam Force buffer type
 * \tparam FourCenterType The bond type to compute; used for type deduction
 * \tparam Cartesian vector type
 * \tparam PBC type
 * \param[in] index The atom and parameter indices used for computing the interaction
 * \param[in] parameters The full type-specific interaction list
 * \param[in] x The coordinates
 * \param[in/out] forces The forces
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template <class Buffer, class FourCenterType, class BasicVector, class ShiftForce, class Pbc,
          std::enable_if_t<Contains<FourCenterType, SupportedFourCenterTypes>{}>* = nullptr>
inline NBLIB_ALWAYS_INLINE
auto dispatchInteraction(InteractionIndex<FourCenterType> index,
                         gmx::ArrayRef<const FourCenterType> parameters,
                         gmx::ArrayRef<const BasicVector> x,
                         Buffer* forces,
                         gmx::ArrayRef<ShiftForce> shiftForces,
                         const Pbc& pbc)
{
    using RealScalar = BasicVectorValueType_t<BasicVector>;
    KernelEnergy<RealScalar> energy;

    int i = std::get<0>(index);
    int j = std::get<1>(index);
    int k = std::get<2>(index);
    int l = std::get<3>(index);

    BasicVector xi = x[i];
    BasicVector xj = x[j];
    BasicVector xk = x[k];
    BasicVector xl = x[l];

    BasicVector fi{0,0,0}, fj{0,0,0}, fk{0,0,0}, fl{0,0,0};

    BasicVector dxIJ, dxKJ, dxKL, dxLJ;
    int sIdx_ij = pbc.dxAiuc(xi, xj, dxIJ);
    int sIdx_kj = pbc.dxAiuc(xk, xj, dxKJ);
    pbc.dxAiuc(xk, xl, dxKL);

    int sIdx_lj = pbc.dxAiuc(xl, xj, dxLJ);

    ShiftForce* shf_ij = accessShiftForces(sIdx_ij, shiftForces);
    ShiftForce* shf_kj = accessShiftForces(sIdx_kj, shiftForces);
    ShiftForce* shf_lj = accessShiftForces(sIdx_lj, shiftForces);
    ShiftForce* shf_c  = accessShiftForces(gmx::c_centralShiftIndex, shiftForces);

    FourCenterType fourCenterTypeParams = parameters[std::get<4>(index)];

    BasicVector m, n;
    RealScalar phi = dihedralPhi(dxIJ, dxKJ, dxKL, &m, &n);

    auto [force, kernelEnergy] = fourCenterKernel(phi, fourCenterTypeParams);

    energy.carrier()              = kernelEnergy;
    energy.threeCenterAggregate() = addThreeCenterAggregate(fourCenterTypeParams, dxIJ, dxKJ, dxKL, &fi, &fj, &fk, &fl);

    spreadFourCenterForces(force, dxIJ, dxKJ, dxKL, m, n, &fi, &fj, &fk, &fl, shf_ij, shf_kj, shf_lj, shf_c);

    (*forces)[i] += fi;
    (*forces)[j] += fj;
    (*forces)[k] += fk;
    (*forces)[l] += fl;

    return energy;
}

/*! \brief Calculate five-center interactions
 *
 * \tparam Force buffer type
 * \tparam FiveCenterType The bond type to compute; used for type deduction
 * \tparam Cartesian vector type
 * \tparam PBC type
 * \param[in] index The atom and parameter indices used for computing the interaction
 * \param[in] parameters The full type-specific interaction list
 * \param[in] x The coordinates
 * \param[in/out] forces The forces
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template <class Buffer, class FiveCenterType, class BasicVector, class ShiftForce, class Pbc,
          std::enable_if_t<Contains<FiveCenterType, SupportedFiveCenterTypes>{}>* = nullptr>
inline NBLIB_ALWAYS_INLINE
auto dispatchInteraction(InteractionIndex<FiveCenterType> index,
                         gmx::ArrayRef<const FiveCenterType> parameters,
                         gmx::ArrayRef<const BasicVector> x,
                         Buffer* forces,
                         [[maybe_unused]] gmx::ArrayRef<ShiftForce> shiftForces,
                         const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<BasicVector>> energy;

    int i = std::get<0>(index);
    int j = std::get<1>(index);
    int k = std::get<2>(index);
    int l = std::get<3>(index);
    int m = std::get<4>(index);

    BasicVector xi = x[i];
    BasicVector xj = x[j];
    BasicVector xk = x[k];
    BasicVector xl = x[l];
    BasicVector xm = x[m];

    BasicVector dxIJ, dxJK, dxKL, dxLM;
    pbc.dxAiuc(xi, xj, dxIJ);
    pbc.dxAiuc(xj, xk, dxJK);
    pbc.dxAiuc(xk, xl, dxKL);
    pbc.dxAiuc(xl, xm, dxLM);

    FiveCenterType fiveCenterTypeParams = parameters[std::get<5>(index)];

    // this dispatch function is not in use yet, because CMap is not yet implemented
    // we don't want to add [[maybe_unused]] in the signature since the params will
    // be used once CMap is implemented, and we also don't want compiler warnings,
    // so we cast to void.
    (void)fiveCenterTypeParams;
    (void)forces;

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
template <class Index, class InteractionType, class Buffer, class ShiftForce, class Pbc>
auto computeForces(gmx::ArrayRef<const Index> indices,
                   gmx::ArrayRef<const InteractionType> parameters,
                   gmx::ArrayRef<const Vec3> x,
                   Buffer* forces,
                   gmx::ArrayRef<ShiftForce> shiftForces,
                   const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<Vec3>> energy;

    for (const auto& index : indices)
    {
        energy += dispatchInteraction(index, parameters, x, forces, shiftForces, pbc);
    }

    return energy;
}

//! \brief convenience overload without shift forces
template <class Index, class InteractionType, class Buffer, class Pbc>
auto computeForces(gmx::ArrayRef<const Index> indices,
                   gmx::ArrayRef<const InteractionType> parameters,
                   gmx::ArrayRef<const Vec3> x,
                   Buffer* forces,
                   const Pbc& pbc)
{
    return computeForces(indices, parameters, x, forces, gmx::ArrayRef<std::nullptr_t>{}, pbc);
}

/*! \brief implement a loop over bond types and accumulate their force contributions
 *
 * \param[in] interactions interaction pairs and bond parameters
 * \param[in] x coordinate input
 * \param[in/out] forces output force buffer
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template<class Buffer, class ShiftForce, class Pbc>
auto reduceListedForces(const ListedInteractionData& interactions,
                        gmx::ArrayRef<const Vec3> x,
                        Buffer* forces,
                        gmx::ArrayRef<ShiftForce> shiftForces,
                        const Pbc& pbc)
{
    using ValueType = BasicVectorValueType_t<Vec3>;
    std::array<ValueType, std::tuple_size<ListedInteractionData>::value> energies{0};
    energies.fill(0);

    // calculate one bond type
    auto computeForceType = [forces, x, shiftForces, &energies, &pbc](const auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;

        gmx::ArrayRef<const InteractionIndex<InteractionType>> indices(interactionElement.indices);
        gmx::ArrayRef<const InteractionType>                   parameters(interactionElement.parameters);

        KernelEnergy<ValueType> energy = computeForces(indices, parameters, x, forces, shiftForces, pbc);

        energies[CarrierIndex<InteractionType>{}] += energy.carrier();
        energies[TwoCenterAggregateIndex<InteractionType>{}] += energy.twoCenterAggregate();
        energies[ThreeCenterAggregateIndex<InteractionType>{}] += energy.threeCenterAggregate();
        // energies[FepIndex{}] += energy.freeEnergyDerivative();
    };

    // calculate all bond types, returns a tuple with the energies for each type
    for_each_tuple(computeForceType, interactions);

    return energies;
}

//! \brief convenience overload without shift forces
template<class Buffer, class Pbc>
auto reduceListedForces(const ListedInteractionData& interactions,
                        gmx::ArrayRef<const Vec3> x,
                        Buffer* forces,
                        const Pbc& pbc)
{
    return reduceListedForces(interactions, x, forces, gmx::ArrayRef<std::nullptr_t>{}, pbc);
}

} // namespace nblib

#undef NBLIB_ALWAYS_INLINE

#endif // NBLIB_LISTEDFORCES_DATAFLOW_HPP
