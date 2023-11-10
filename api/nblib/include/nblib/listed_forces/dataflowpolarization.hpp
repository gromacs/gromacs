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

#ifndef NBLIB_LISTEDFORCES_DATAFLOWPOLARIZATION_HPP
#define NBLIB_LISTEDFORCES_DATAFLOWPOLARIZATION_HPP

#include "nblib/util/util.hpp"

namespace nblib
{

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
template<class Buffer, class ChargedType, class MemVector, class ShiftForce, class Pbc, class Lambda>
HOST_DEVICE_INLINE auto dispatchChargedInteraction(IndexArray<3>      index,
                                                   const ChargedType* bondInstancesA,
                                                   const ChargedType* bondInstancesB,
                                                   const MemVector*   x,
                                                   const real*        q,
                                                   Lambda             lambda,
                                                   Buffer*            forces,
                                                   ShiftForce*        shiftForces,
                                                   const Pbc&         pbc)
{
    using ValueType = VectorValueType_t<MemVector>;
    using Vec       = StackVec3<ValueType>;
    KernelEnergy<ValueType> energy;

    int i  = index[0];
    int j  = index[1];
    Vec xi = loadVec(x[i]);
    Vec xj = loadVec(x[j]);

    ChargedType bondA = bondInstancesA[index[2]];
    // conditional load of B parameters only if Lambda is not NoFepLambdaType
    ChargedType bondB = loadInteractionParameters<Lambda>(bondInstancesB, index[2]);

    Vec fi{ 0, 0, 0 }, fj{ 0, 0, 0 };

    Vec dx;
    // calculate xi - xj modulo pbc
    int sIdx = pbc.dxAiuc(xi, xj, dx);

    ValueType charge = q[j];
    if constexpr (Contains<ChargedType, PairListedTypes>{})
    {
        charge *= q[i];
    }

    auto interactionEnergies = computeTwoCenter(bondA, bondB, dx, charge, lambda, &fi, &fj);
    addEnergy<ChargedType>(&energy, interactionEnergies);

    addForce(&(*forces)[i], fi);
    addForce(&(*forces)[j], fj);

    addForce(shiftForces + sIdx, fi);
    addForce(shiftForces + gmx::c_centralShiftIndex, fj);

    return energy;
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
HOST_INLINE auto dispatchChargedInteraction(IndexArray<5>                       index,
                                            const FourCenterType*               parametersA,
                                            const FourCenterType*               parametersB,
                                            const MemVector*                    x,
                                            const VectorValueType_t<MemVector>* q,
                                            Lambda                              lambda,
                                            Buffer*                             forces,
                                            ShiftForce*                         shiftForces,
                                            const Pbc&                          pbc)
{
    using ValueType = VectorValueType_t<MemVector>;
    using Vec       = StackVec3<ValueType>;
    KernelEnergy<ValueType> energy;

    int i = index[0];
    int j = index[1];
    int k = index[2];
    int l = index[3];

    Vec xi = loadVec(x[i]);
    Vec xj = loadVec(x[j]);
    Vec xk = loadVec(x[k]);
    Vec xl = loadVec(x[l]);

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
            addPairAggregate(paramsA, paramsB, xi, xl, q[i], q[l], lambda, &fi, &fl, pbc);

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
template<class Index, class InteractionType, class MemVector, class Buffer, class ShiftForce, class Pbc, class Lambda>
auto computePolarizationForces(gmx::ArrayRef<const Index>           indices,
                               gmx::ArrayRef<const InteractionType> parametersA,
                               gmx::ArrayRef<const InteractionType> parametersB,
                               gmx::ArrayRef<const MemVector>       x,
                               gmx::ArrayRef<const real>            q,
                               Lambda                               lambda,
                               Buffer*                              forces,
                               gmx::ArrayRef<ShiftForce>            shiftForces,
                               const Pbc&                           pbc)
{
    KernelEnergy<VectorValueType_t<MemVector>> energy;
    for (const auto& index : indices)
    {
        energy += dispatchChargedInteraction(index,
                                             parametersA.data(),
                                             parametersB.data(),
                                             x.data(),
                                             q.data(),
                                             lambda,
                                             forces,
                                             shiftForces.data(),
                                             pbc);
    }
    return energy;
}

//! \brief convenience overload without FEP
template<class Index, class InteractionType, class MemVector, class Buffer, class ShiftForce, class Pbc>
auto computePolarizationForces(gmx::ArrayRef<const Index>           indices,
                               gmx::ArrayRef<const InteractionType> parameters,
                               gmx::ArrayRef<const MemVector>       x,
                               gmx::ArrayRef<const real>            q,
                               Buffer*                              forces,
                               gmx::ArrayRef<ShiftForce>            shiftForces,
                               const Pbc&                           pbc)
{
    return computePolarizationForces(
            indices, parameters, parameters, x, q, NoFepLambdaType{}, forces, shiftForces, pbc);
}

/*! \brief implement a loop over polarization types and accumulate their force contributions
 *
 * \param[in] interactions interaction pairs and bond parameters
 * \param[in] x coordinate input
 * \param[in] q charge input
 * \param[in/out] forces output force buffer
 * \param[in] pbc Object used for computing distances accounting for PBC's
 * \return Computed kernel energies
 */
template<class Buffer, class MemVector, class ShiftForce, class Pbc>
auto reducePolarization(const ListedInteractionData&   interactions,
                        gmx::ArrayRef<const MemVector> x,
                        gmx::ArrayRef<const real>      q,
                        Buffer*                        forces,
                        gmx::ArrayRef<ShiftForce>      shiftForces,
                        const Pbc&                     pbc)
{
    using ValueType = VectorValueType_t<MemVector>;

    if (x.size() != q.size())
    {
        throw InputException("Charges array size mismatch");
    }

    ListedEnergies energies{ 0 };
    std::fill(energies.begin(), energies.end(), 0);

    // calculate one bond type
    auto computeForceType = [forces, x, q, shiftForces, &energies, &pbc](const auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;

        gmx::ArrayRef<const InteractionIndex<InteractionType>> indices(interactionElement.indices);
        gmx::ArrayRef<const InteractionType> parameters(interactionElement.parametersA);

        KernelEnergy<ValueType> energy =
                computePolarizationForces(indices, parameters, x, q, forces, shiftForces, pbc);

        energies[CarrierIndex<InteractionType>{}] += energy.carrier();
        energies[TwoCenterAggregateIndex<InteractionType>{}] += energy.twoCenterAggregate();
        energies[ThreeCenterAggregateIndex<InteractionType>{}] += energy.threeCenterAggregate();
        energies[VdwIndex{}] += energy.eVdw();
        energies[CoulombIndex{}] += energy.eCoul();
        energies[FepIndex{}] += energy.freeEnergyDerivative();
    };

    auto computeIndices = subsetIndices(ChargedListedTypes{}, AllListedTypes{});
    // calculate all bond types, returns a tuple with the energies for each type
    for_each_tuple(computeForceType, tieElements(interactions, computeIndices));

    return energies;
}

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_DATAFLOWPOLARIZATION_HPP
