/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
#include "gromacs/math/vec.h"
#include "gromacs/utility/arrayref.h"

#define NBLIB_ALWAYS_INLINE __attribute((always_inline))

namespace nblib
{

template<class TwoCenterType, class BasicVector>
inline NBLIB_ALWAYS_INLINE
auto computeTwoCenter(const TwoCenterType& parameters, const BasicVector& dx, BasicVector* fi, BasicVector* fj)
{
    using ValueType = BasicVectorValueType_t<BasicVector>;

    ValueType dr2 = dot(dx, dx);
    ValueType dr  = std::sqrt(dr2);
    auto [force, energy] = bondKernel(dr, parameters);

    // avoid division by 0
    if (dr2 != 0.0)
    {
        force /= dr;
        spreadTwoCenterForces(force, dx, fi, fj);
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
template <class Buffer, class TwoCenterType, class BasicVector, class Pbc,
          std::enable_if_t<Contains<TwoCenterType, SupportedTwoCenterTypes>{}>* = nullptr>
inline NBLIB_ALWAYS_INLINE
auto dispatchInteraction(const InteractionIndex<TwoCenterType>& index,
                         const std::vector<TwoCenterType>& bondInstances,
                         gmx::ArrayRef<const BasicVector> x,
                         Buffer* forces,
                         const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<BasicVector>> energy;

    int i = std::get<0>(index);
    int j = std::get<1>(index);
    const gmx::RVec& x1 = x[i];
    const gmx::RVec& x2 = x[j];
    const TwoCenterType& bond = bondInstances[std::get<2>(index)];

    gmx::RVec dx;
    // calculate x1 - x2 modulo pbc
    pbc.dxAiuc(x1, x2, dx);

    energy.carrier() = computeTwoCenter(bond, dx, &(*forces)[i], &(*forces)[j]);
    return energy;
}


template<class ThreeCenterType, class BasicVector>
inline NBLIB_ALWAYS_INLINE
std::enable_if_t<HasTwoCenterAggregate<ThreeCenterType>::value, BasicVectorValueType_t<BasicVector>>
addTwoCenterAggregate(const ThreeCenterType& parameters, const BasicVector& rij, const BasicVector& rkj,
                      BasicVector* fi, BasicVector* fj, BasicVector* fk)
{
    if (parameters.manifest == ThreeCenterType::Cargo::ij)
    {
        // i-j bond
        return computeTwoCenter(parameters.twoCenter(), rij, fi, fj);
    }
    if (parameters.manifest == ThreeCenterType::Cargo::jk)
    {
        // j-k bond
        return computeTwoCenter(parameters.twoCenter(), rkj, fk, fj);
    }

    // aggregate is empty
    return 0.0;
};

template<class ThreeCenterType, class BasicVector>
inline NBLIB_ALWAYS_INLINE
std::enable_if_t<!HasTwoCenterAggregate<ThreeCenterType>::value, BasicVectorValueType_t<BasicVector>>
addTwoCenterAggregate([[maybe_unused]] const ThreeCenterType& parameters,
                      [[maybe_unused]] const BasicVector& rij,
                      [[maybe_unused]] const BasicVector& rkj,
                      [[maybe_unused]] BasicVector* fi,
                      [[maybe_unused]] BasicVector* fj,
                      [[maybe_unused]] BasicVector* fk)
{
    return 0.0;
};

template<class ThreeCenterType, class BasicVector>
inline NBLIB_ALWAYS_INLINE
auto computeThreeCenter(const ThreeCenterType& parameters, const BasicVector& rij, const BasicVector& rkj,
                        BasicVector* fi, BasicVector* fj, BasicVector* fk)
{
    using ValueType = BasicVectorValueType_t<BasicVector>;
    // calculate 3-center common quantities: angle between x1-x2 and x2-x3
    // Todo: after sufficient evaluation, switch over to atan2 based algorithm
    ValueType costh = cos_angle(rij, rkj); /* 25             */
    ValueType theta = std::acos(costh);    /* 10             */

    // call type-specific angle kernel, e.g. harmonic, linear, quartic,  etc.
    auto [force, energy] = threeCenterKernel(theta, parameters);

    spreadThreeCenterForces(costh, force, rij, rkj, fi, fj, fk);

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
template <class Buffer, class ThreeCenterType, class BasicVector, class Pbc,
          std::enable_if_t<Contains<ThreeCenterType, SupportedThreeCenterTypes>{}>* = nullptr>
inline NBLIB_ALWAYS_INLINE
auto dispatchInteraction(const InteractionIndex<ThreeCenterType>& index,
                         const std::vector<ThreeCenterType>& parameters,
                         gmx::ArrayRef<const BasicVector> x,
                         Buffer* forces,
                         const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<BasicVector>> energy;

    //! fetch input data: position vectors x1-x3 and interaction parameters
    int i = std::get<0>(index);
    int j = std::get<1>(index);
    int k = std::get<2>(index);
    const gmx::RVec& xi = x[i];
    const gmx::RVec& xj = x[j];
    const gmx::RVec& xk = x[k];
    const ThreeCenterType& threeCenterParameters = parameters[std::get<3>(index)];

    gmx::RVec fi{0,0,0}, fj{0,0,0}, fk{0,0,0};

    gmx::RVec rij, rkj;
    pbc.dxAiuc(xi, xj, rij); /*  3              */
    pbc.dxAiuc(xk, xj, rkj); /*  3              */

    energy.carrier()            = computeThreeCenter(threeCenterParameters, rij, rkj, &fi, &fj, &fk);
    energy.twoCenterAggregate() = addTwoCenterAggregate(threeCenterParameters, rij, rkj, &fi, &fj, &fk);

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
addThreeCenterAggregate([[maybe_unused]] const FourCenterType& parameters,
                        [[maybe_unused]] const BasicVector& rij,
                        [[maybe_unused]] const BasicVector& rkj,
                        [[maybe_unused]] const BasicVector& rkl,
                        [[maybe_unused]] BasicVector* fi,
                        [[maybe_unused]] BasicVector* fj,
                        [[maybe_unused]] BasicVector* fk,
                        [[maybe_unused]] BasicVector* fl)
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
template <class Buffer, class FourCenterType, class BasicVector, class Pbc,
          std::enable_if_t<Contains<FourCenterType, SupportedFourCenterTypes>{}>* = nullptr>
inline NBLIB_ALWAYS_INLINE
auto dispatchInteraction(const InteractionIndex<FourCenterType>& index,
                         const std::vector<FourCenterType>& parameters,
                         gmx::ArrayRef<const BasicVector> x,
                         Buffer* forces,
                         const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<BasicVector>> energy;

    int i = std::get<0>(index);
    int j = std::get<1>(index);
    int k = std::get<2>(index);
    int l = std::get<3>(index);

    const gmx::RVec& xi = x[i];
    const gmx::RVec& xj = x[j];
    const gmx::RVec& xk = x[k];
    const gmx::RVec& xl = x[l];

    gmx::RVec fi{0,0,0}, fj{0,0,0}, fk{0,0,0}, fl{0,0,0};

    gmx::RVec dxIJ, dxKJ, dxKL;
    pbc.dxAiuc(xi, xj, dxIJ);
    pbc.dxAiuc(xk, xj, dxKJ);
    pbc.dxAiuc(xk, xl, dxKL);

    const FourCenterType& fourCenterTypeParams = parameters[std::get<4>(index)];

    rvec m, n;
    real phi = dihedralPhi(dxIJ, dxKJ, dxKL, m, n);

    auto [force, kernelEnergy] = fourCenterKernel(phi, fourCenterTypeParams);

    energy.carrier()              = kernelEnergy;
    energy.threeCenterAggregate() = addThreeCenterAggregate(fourCenterTypeParams, dxIJ, dxKJ, dxKL, &fi, &fj, &fk, &fl);

    spreadFourCenterForces(force, dxIJ, dxKJ, dxKL, m, n, &fi, &fj, &fk, &fl);

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
template <class Buffer, class FiveCenterType, class BasicVector, class Pbc,
          std::enable_if_t<Contains<FiveCenterType, SupportedFiveCenterTypes>{}>* = nullptr>
inline NBLIB_ALWAYS_INLINE
auto dispatchInteraction(const InteractionIndex<FiveCenterType>& index,
                         const std::vector<FiveCenterType>& parameters,
                         gmx::ArrayRef<const BasicVector> x,
                         Buffer* forces,
                         const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<BasicVector>> energy;

    int i = std::get<0>(index);
    int j = std::get<1>(index);
    int k = std::get<2>(index);
    int l = std::get<3>(index);
    int m = std::get<4>(index);

    const gmx::RVec& xi = x[i];
    const gmx::RVec& xj = x[j];
    const gmx::RVec& xk = x[k];
    const gmx::RVec& xl = x[l];
    const gmx::RVec& xm = x[m];

    gmx::RVec dxIJ, dxJK, dxKL, dxLM;
    pbc.dxAiuc(xi, xj, dxIJ);
    pbc.dxAiuc(xj, xk, dxJK);
    pbc.dxAiuc(xk, xl, dxKL);
    pbc.dxAiuc(xl, xm, dxLM);

    const FiveCenterType& fiveCenterTypeParams = parameters[std::get<5>(index)];

    // this dispatch function is not in use yet, because CMap is not yet implemented
    // we don't want to add [[maybe_unused]] in the signature
    // and we also don't want compiler warnings, so we cast to void
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
template <class Index, class InteractionType, class Buffer, class Pbc>
auto computeForces(const std::vector<Index>& indices,
                   const std::vector<InteractionType>& interactionParameters,
                   gmx::ArrayRef<const Vec3> x,
                   Buffer* forces,
                   const Pbc& pbc)
{
    KernelEnergy<BasicVectorValueType_t<Vec3>> energy;

    for (const auto& index : indices)
    {
        energy += dispatchInteraction(index, interactionParameters, x, forces, pbc);
    }

    return energy;
}

/*! \brief implement a loop over bond types and accumulate their force contributions
 *
 * \param[in] interactions interaction pairs and bond parameters
 * \param[in] x coordinate input
 * \param[in/out] forces output force buffer
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template<class Buffer, class Pbc>
auto reduceListedForces(const ListedInteractionData& interactions,
                        gmx::ArrayRef<const Vec3> x,
                        Buffer* forces,
                        const Pbc& pbc)
{
    using ValueType = BasicVectorValueType_t<Vec3>;
    std::array<ValueType, std::tuple_size<ListedInteractionData>::value> energies{0};
    energies.fill(0);

    // calculate one bond type
    auto computeForceType = [forces, &x, &energies, &pbc](const auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;

        KernelEnergy<ValueType> energy = computeForces(interactionElement.indices, interactionElement.parameters, x, forces, pbc);

        energies[CarrierIndex<InteractionType>{}] += energy.carrier();
        energies[TwoCenterAggregateIndex<InteractionType>{}] += energy.twoCenterAggregate();
        energies[ThreeCenterAggregateIndex<InteractionType>{}] += energy.threeCenterAggregate();
        // energies[FepIndex{}] += energy.freeEnergyDerivative();
    };

    // calculate all bond types, returns a tuple with the energies for each type
    for_each_tuple(computeForceType, interactions);

    return energies;
}

} // namespace nblib

#undef NBLIB_ALWAYS_INLINE

#endif // NBLIB_LISTEDFORCES_DATAFLOW_HPP
