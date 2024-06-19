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
 * This implements a fixture for calling calc_listed in gromacs.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxcalculator.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <tuple>
#include <type_traits>

#include "listed_forces/conversionscommon.h"

#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/booltype.h"

#include "nblib/exception.h"
#include "nblib/listed_forces/bondtypes.h"
#include "nblib/util/traits.hpp"
#include "nblib/util/util.hpp"

namespace nblib
{

ListedGmxCalculator::ListedGmxCalculator(const ListedInteractionData& interactions,
                                         int                          nP,
                                         int                          nThr,
                                         const Box&                   box) :
    numParticles(nP),
    numThreads(nThr),
    box_(box),
    shiftBuffer(gmx::c_numShiftVectors, gmx::RVec{ 0, 0, 0 }),
    forceBuffer(2 * numParticles, gmx::RVec{ 0, 0, 0 }),
    shiftProxy(gmx::ArrayRefWithPadding<gmx::RVec>(forceBuffer.data(),
                                                   forceBuffer.data() + numParticles,
                                                   forceBuffer.data() + 2 * numParticles),
               true,
               shiftBuffer),
    virialProxy(forceBuffer, true),
    forceOutputs(shiftProxy, true, virialProxy),
    enerd(1, nullptr),
    lambdaBuffer(42) // values unused; just initialized with something larger than the number of enum types in FreeEnergyPerturbationCouplingType
{
    std::tie(idef, ffparams) = convertToGmxInteractions(interactions);
    idef->ilsort             = ilsortNO_FE;

    ListedForces::InteractionSelection interactionSelection = ListedForces::interactionSelectionAll();

    // no distance restraints
    disres_.nres   = 0;
    fcdata_.disres = &disres_;

    gmxListedForces_ =
            std::make_unique<ListedForces>(*ffparams, 1, numThreads, interactionSelection, nullptr);
    gmxListedForces_->setup(*idef, nP, false);

    wcycle = wallcycle_init(nullptr, 0, &cr);
    set_pbc(&pbc, PbcType::Xyz, box_.legacyMatrix());

    stepWork.computeForces       = true;
    stepWork.computeListedForces = true;
    stepWork.computeDhdl         = false;
    stepWork.computeEnergy       = true;

    fr.natoms_force = numParticles;

    mdatoms_.nr = nP;
}

void ListedGmxCalculator::compute(gmx::ArrayRef<const gmx::RVec>     x,
                                  gmx::ArrayRef<gmx::RVec>           forces,
                                  gmx::ArrayRef<gmx::RVec>           shiftForces,
                                  ListedForceCalculator::EnergyType& energies,
                                  bool                               usePbc)
{
    if (forces.size() != x.size() || forces.size() >= forceBuffer.size())
    {
        throw InputException("Provided force and/or coordinate buffers inconsistent");
    }

    // gmx listed forces only uses PBC if this flag is true
    fr.bMolPBC = usePbc;

    // the forceOutput proxies also have computeVirial flags, but they are not accessed
    // by gmx::ListedForces::calculate(), only stepWork matters
    if (shiftForces.size() == gmx::c_numShiftVectors)
    {
        stepWork.computeVirial = true;
    }
    else
    {
        stepWork.computeVirial = false;
    }

    // zero output data structure
    energies.fill(0);
    std::fill(std::begin(enerd.term), std::end(enerd.term), 0.0);
    std::fill(forceBuffer.begin(), forceBuffer.end(), gmx::RVec{ 0, 0, 0 });
    std::fill(shiftBuffer.begin(), shiftBuffer.end(), gmx::RVec{ 0, 0, 0 });

    gmxListedForces_->calculate(wcycle.get(),
                                box_.legacyMatrix(),
                                &cr,
                                nullptr,
                                { x.data(), x.data() + x.size(), x.data() + x.size() },
                                gmx::ArrayRef<const gmx::RVec>{},
                                &fcdata_,
                                nullptr,
                                &forceOutputs,
                                &fr,
                                &pbc,
                                &enerd,
                                &nrnb,
                                lambdaBuffer,
                                mdatoms_.chargeA,
                                mdatoms_.chargeB,
                                makeConstArrayRef(mdatoms_.bPerturbed),
                                mdatoms_.cENER,
                                mdatoms_.nPerturbed,
                                nullptr,
                                stepWork);

    auto transferEnergy = [&energies, this](auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;
        if constexpr (ListedTypeIsImplemented<InteractionType>{})
        {
            constexpr int index       = FindIndex<InteractionType, ListedInteractionData>::value;
            auto          gmxListedID = FindIndex<InteractionType, GmxToNblibMapping>::value;
            energies[index]           = this->enerd.term[gmxListedID];
        }
    };
    for_each_tuple(transferEnergy, ListedInteractionData{});

    // add forces to output force buffers
    for (int pIndex = 0; pIndex < forces.ssize(); pIndex++)
    {
        forces[pIndex] += forceBuffer[pIndex];
    }

    // add shift forces to the output buffer
    for (int i = 0; i < shiftForces.ssize(); ++i)
    {
        shiftForces[i] += shiftBuffer[i];
    }
}

void ListedGmxCalculator::compute(gmx::ArrayRef<const gmx::RVec>     x,
                                  gmx::ArrayRef<gmx::RVec>           forces,
                                  ListedForceCalculator::EnergyType& energies,
                                  bool                               usePbc)
{
    compute(x, forces, gmx::ArrayRef<gmx::RVec>{}, energies, usePbc);
}

const InteractionDefinitions& ListedGmxCalculator::getIdef() const
{
    return *idef;
}


} // namespace nblib
