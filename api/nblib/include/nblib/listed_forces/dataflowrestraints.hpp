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

#ifndef NBLIB_LISTEDFORCES_DATAFLOWRESTRAINTS_HPP
#define NBLIB_LISTEDFORCES_DATAFLOWRESTRAINTS_HPP

#include <vector>

#include "gromacs/utility/arrayref.h"

#include "nblib/listed_forces/dataflow.hpp"
#include "nblib/listed_forces/kernels.hpp"
#include "nblib/listed_forces/positionrestraints.hpp"
#include "nblib/listed_forces/traits.h"
#include "nblib/pbc.hpp"
#include "nblib/util/util.hpp"

namespace nblib
{

template<class OneCenterType, class StackVector, class Lambda, class Virial>
HOST_DEVICE_INLINE auto computeOneCenter(const OneCenterType& parameterA,
                                         const OneCenterType& /* parameterB */,
                                         const StackVector& dx,
                                         const StackVector& rdist,
                                         const Lambda /* lambda */,
                                         StackVector* fi,
                                         Virial*      virial)
{
    using ValueType    = VectorValueType_t<StackVector>;
    ValueType   energy = 0;
    StackVector vir    = { 0, 0, 0 };
    StackVector force  = { 0, 0, 0 };

    for (int m = 0; m < 3; m++)
    {
        auto [forceM, epot] = oneCenterKernel(dx[m], parameterA, m);
        force[m]            = forceM;

        energy += epot;
        vir[m] -= ValueType(0.5) * (dx[m] + rdist[m]) * forceM;
    }

    virial[0] = vir[0]; // XX
    virial[4] = vir[1]; // YY
    virial[8] = vir[2]; // ZZ

    *fi += force;

    return energy;
}

template<class OneCenterType, class MemVector, class Lambda, class Buffer, class Virial, class Pbc>
HOST_DEVICE_INLINE auto dispatchRestraints(IndexArray<2>                           index,
                                           const OneCenterType*                    bondInstancesA,
                                           const OneCenterType*                    bondInstancesB,
                                           const MemVector*                        x,
                                           Lambda                                  lambda,
                                           Buffer*                                 forces,
                                           Virial*                                 virials,
                                           const Pbc&                              pbc,
                                           const Box&                              box,
                                           PbcType                                 pbcType,
                                           RefCoordScaling                         refcoord_scaling,
                                           StackVec3<VectorValueType_t<MemVector>> com)
{
    using ValueType = VectorValueType_t<MemVector>;
    using Vec       = StackVec3<ValueType>;
    KernelEnergy<ValueType> energy;

    int i  = index[0];
    Vec xi = loadVec(x[i]);

    OneCenterType bondA = bondInstancesA[index[1]];
    // conditional load of B parameters only if Lambda is not NoFepLambdaType
    OneCenterType bondB = loadInteractionParameters<Lambda>(bondInstancesB, index[1]);

    Vec fi{ 0, 0, 0 };

    Vec com_sc       = posresCOM(box, refcoord_scaling, pbcType, com);
    auto [dx, rdist] = posres_dx(xi, com, com_sc, pbc, box, refcoord_scaling, numPbcDimensions(pbcType));

    energy.carrier() = computeOneCenter(bondA, bondB, dx, rdist, lambda, &fi, virials);

    addForce(&(*forces)[i], fi);

    return energy;
}

template<class Index, class InteractionType, class Buffer, class Virial, class MemVector, class Pbc>
auto computeForcesRestraints(gmx::ArrayRef<const Index>              indices,
                             gmx::ArrayRef<const InteractionType>    parametersA,
                             gmx::ArrayRef<const InteractionType>    parametersB,
                             gmx::ArrayRef<const MemVector>          x,
                             Buffer*                                 forces,
                             gmx::ArrayRef<Virial>                   virials,
                             const Pbc&                              pbc,
                             const Box&                              box,
                             PbcType                                 pbcType,
                             RefCoordScaling                         refcoord_scaling,
                             StackVec3<VectorValueType_t<MemVector>> com)
{
    KernelEnergy<VectorValueType_t<MemVector>> energy;

    for (const auto& index : indices)
    {
        energy += dispatchRestraints(index,
                                     parametersA.data(),
                                     parametersB.data(),
                                     x.data(),
                                     NoFepLambdaType{},
                                     forces,
                                     virials.data(),
                                     pbc,
                                     box,
                                     pbcType,
                                     refcoord_scaling,
                                     com);
    }

    return energy;
}

template<class Buffer, class Virial, class Pbc, class MemVector, class StackVector>
auto reduceRestraints(const ListedInteractionData&   interactions,
                      gmx::ArrayRef<const MemVector> x,
                      Buffer*                        forces,
                      gmx::ArrayRef<Virial>          virials,
                      const Pbc&                     pbc,
                      const Box&                     box,
                      PbcType                        pbcType,
                      RefCoordScaling                refcoord_scaling,
                      StackVector                    com)
{
    using ValueType = VectorValueType_t<MemVector>;

    if (virials.size() != 9)
    {
        throw InputException("Virial forces array size mismatch");
    }

    ListedEnergies energies{ 0 };
    std::fill(energies.begin(), energies.end(), 0);

    // calculate one bond type
    auto computeForceType = [forces, x, virials, &energies, &pbc, box, refcoord_scaling, pbcType, com](
                                    const auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;
        gmx::ArrayRef<const InteractionIndex<InteractionType>> indices(interactionElement.indices);
        gmx::ArrayRef<const InteractionType> parametersA(interactionElement.parametersA);
        gmx::ArrayRef<const InteractionType> parametersB(interactionElement.parametersB);

        KernelEnergy<ValueType> energy =
                computeForcesRestraints(indices,
                                        parametersA,
                                        parametersB,
                                        x,
                                        forces,
                                        virials,
                                        pbc,
                                        box,
                                        pbcType,
                                        refcoord_scaling,
                                        StackVec3<ValueType>{ com[0], com[1], com[2] });
        energies[FindIndex<InteractionType, AllListedTypes>{}] += energy.carrier();
        energies[FepIndex{}] += energy.freeEnergyDerivative();
    };

    auto computeIndices = subsetIndices(RestraintTypes{}, AllListedTypes{});
    // calculate all bond types, returns a tuple with the energies for each type
    for_each_tuple(computeForceType, tieElements(interactions, computeIndices));

    return energies;
}


} // namespace nblib

#endif // NBLIB_LISTEDFORCES_DATAFLOWRESTRAINTS_HPP
