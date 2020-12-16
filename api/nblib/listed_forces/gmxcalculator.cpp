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

/*! \internal \file
 * \brief
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/exception.h"
#include "nblib/listed_forces/gmxcalculator.h"
#include "gromacs/listed_forces/listed_forces.cpp"
#include "gromacs/listed_forces/manage_threading.h"
#include "gromacs/timing/wallcycle.h"
#include "conversions.hpp"

namespace nblib
{

ListedGmxCalculator::ListedGmxCalculator(const ListedInteractionData& interactions,
                                         int                          nP,
                                         int                          nThr,
                                         const Box&                   box) :
    numParticles(nP),
    numThreads(nThr),
    bondedThreading(numThreads, 1, nullptr),
    shiftBuffer(SHIFTS), // this is gromacs setup code so here use SHIFTS instead of numShiftVectors
    forceBuffer(2 * numParticles, gmx::RVec{ 0, 0, 0 }),
    shiftProxy(gmx::ArrayRefWithPadding<gmx::RVec>(&forceBuffer[0],
                                                   &forceBuffer[numParticles],
                                                   &forceBuffer[2 * numParticles]),
               false,
               shiftBuffer),
    virialProxy(forceBuffer, false),
    forceOutputs(shiftProxy, false, virialProxy),
    enerd(1, 0),
    lambdaBuffer(numParticles) // just something big enough
{
    std::tie(idef, ffparams) = createFFparams(interactions);
    idef->ilsort             = ilsortNO_FE;

    setup_bonded_threading(&bondedThreading, numParticles, false, *idef);

    wcycle = wallcycle_init(nullptr, 0, &cr);
    set_pbc(&pbc, PbcType::Xyz, box.legacyMatrix());

    stepWork.computeDhdl   = false;
    stepWork.computeVirial = false;
    stepWork.computeEnergy = true;

    fr.natoms_force = numParticles;
}

void ListedGmxCalculator::compute(const std::vector<gmx::RVec>&      x,
                                  std::vector<gmx::RVec>&            forces,
                                  ListedForceCalculator::EnergyType& energies)
{
    if (forces.size() != x.size() || forces.size() >= forceBuffer.size())
    {
        throw InputException("Provided force and/or coordinate buffers inconsistent");
    }

    const rvec* xdata = &(x[0].as_vec());

    energies.fill(0);
    std::fill(enerd.term, enerd.term + F_NRE, 0.0);

    calc_listed(wcycle,
                *idef,
                &bondedThreading,
                xdata,
                &forceOutputs,
                &fr,
                &pbc,
                &enerd,
                &nrnb,
                lambdaBuffer.data(),
                nullptr,
                nullptr,
                nullptr,
                stepWork);

    auto transferEnergy = [&energies, this](auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;
        if constexpr (ListedTypeIsImplemented<InteractionType>{})
        {
            constexpr int index = FindIndex<InteractionType, ListedInteractionData>::value;
            energies[index]     = this->enerd.term[ListedIndex<InteractionType>::value];
        }
    };
    for_each_tuple(transferEnergy, ListedInteractionData{});

    // add forces to output force buffers
    for (int pIndex = 0; pIndex < forces.size(); pIndex++)
    {
        forces[pIndex] += forceBuffer[pIndex];
    }
}

const InteractionDefinitions& ListedGmxCalculator::getIdef() const
{
    return *idef;
}


} // namespace nblib
