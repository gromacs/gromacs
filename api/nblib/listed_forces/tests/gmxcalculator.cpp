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

#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/timing/wallcycle.h"

#include "nblib/exception.h"
#include "nblib/listed_forces/conversionscommon.h"

namespace nblib
{

ListedGmxCalculator::ListedGmxCalculator(const InteractionDefinitions& idefs,
                                         const gmx_ffparams_t&         ffparams,
                                         size_t                        nP,
                                         int                           nThr,
                                         const Box&                    box) :
    idef_(std::make_unique<InteractionDefinitions>(idefs)),
    ffparams_(std::make_unique<gmx_ffparams_t>(ffparams)),
    numParticles(nP),
    numThreads(nThr),
    box_(box),
    shiftBuffer(gmx::c_numShiftVectors, gmx::RVec{ 0, 0, 0 }),
    forceBuffer(2 * numParticles, gmx::RVec{ 0, 0, 0 }),
    charges_(nP, 0.0),
    shiftProxy(gmx::ArrayRefWithPadding<gmx::RVec>(&forceBuffer[0],
                                                   &forceBuffer[numParticles],
                                                   &forceBuffer[2 * numParticles]),
               true,
               shiftBuffer),
    virialProxy(forceBuffer, true),
    forceOutputs(shiftProxy, true, virialProxy),
    enerd(1, nullptr),
    lambdaBuffer(42) // values unused; just initialized with something larger than the number of enum types in FreeEnergyPerturbationCouplingType
{
    idef_->ilsort = ilsortNO_FE;

    ListedForces::InteractionSelection interactionSelection = ListedForces::interactionSelectionAll();

    // no distance restraints
    disres_.nres   = 0;
    fcdata_.disres = &disres_;

    gmxListedForces_ =
            std::make_unique<ListedForces>(*ffparams_, 1, numThreads, interactionSelection, nullptr);
    gmxListedForces_->setup(*idef_, numParticles, false);

    wcycle = wallcycle_init(nullptr, 0, &cr);
    set_pbc(&pbc, PbcType::Xyz, box_.legacyMatrix());

    stepWork.computeForces       = true;
    stepWork.computeListedForces = true;
    stepWork.computeDhdl         = false;
    stepWork.computeEnergy       = true;

    fr.natoms_force = numParticles;
    mdatoms_        = { 0 };
    mdatoms_.nr     = nP;
    energyGroup_.resize(nP, 0);
    mdatoms_.cENER = energyGroup_;

    interaction_const_t interactionConst;
    interactionConst.vdwtype      = VanDerWaalsType::Cut;
    interactionConst.vdw_modifier = InteractionModifiers::PotShift;
    interactionConst.rvdw         = 1.0;
    fr.ic = std::make_unique<interaction_const_t>(std::move(interactionConst));

    real tableRange = 2.9;
    fr.pairsTable   = make_tables(nullptr, fr.ic.get(), nullptr, tableRange, GMX_MAKETABLES_14ONLY);

    mdatoms_.nPerturbed = 0;
}

void ListedGmxCalculator::compute(gmx::ArrayRef<const gmx::RVec> x,
                                  gmx::ArrayRef<const real>      q,
                                  gmx::ArrayRef<gmx::RVec>       forces,
                                  gmx::ArrayRef<gmx::RVec>       shiftForces,
                                  gmx::ArrayRef<real>            virials,
                                  ListedEnergies&                energies,
                                  gmx::RVec                      com)
{
    if (forces.size() != x.size() || forces.size() >= forceBuffer.size())
    {
        throw InputException("Provided force and/or coordinate buffers inconsistent");
    }

    // gmx listed forces only uses PBC if this flag is true
    fr.bMolPBC     = true;
    fr.rc_scaling  = RefCoordScaling::No;
    fr.pbcType     = PbcType::Xyz;
    fr.posres_com  = com;
    fr.posres_comB = com;
    fr.fudgeQQ     = 1.0;

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
    std::fill(energies.begin(), energies.end(), 0);
    enerd.grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJ14][0]      = 0;
    enerd.grpp.energyGroupPairTerms[NonBondedEnergyTerms::Coulomb14][0] = 0;
    std::fill(enerd.term.data(), enerd.term.data() + F_NRE, 0.0);
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
                                q,
                                q,
                                makeConstArrayRef(mdatoms_.bPerturbed),
                                mdatoms_.cENER,
                                mdatoms_.nPerturbed,
                                nullptr,
                                stepWork);

    auto transferEnergy = [&energies, this](auto& interactionElement) {
        using InteractionType = std::decay_t<decltype(interactionElement)>;
        if constexpr (ListedTypeIsImplemented<InteractionType>{})
        {
            int index = FindIndex<InteractionType, AllListedTypes>{};

            auto gmxListedID = FindIndex<InteractionType, GmxToNblibMapping>::value;
            energies[index]  = this->enerd.term[gmxListedID];
        }
    };
    for_each_tuple(transferEnergy, Reduce<std::tuple, AllListedTypes>{});
    // We enforce that there is only one energy group in this calculator, so we don't need to accumulate
    energies[VdwIndex{}] += enerd.grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJ14][0];
    energies[CoulombIndex{}] += enerd.grpp.energyGroupPairTerms[NonBondedEnergyTerms::Coulomb14][0];

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

    if (!virials.empty())
    {
        // add virial forces to the output buffer
        auto virialMatrix = forceOutputs.forceWithVirial().getVirial();
        // Only the diagonal elememts are populated by position restraints
        virials[0] += virialMatrix[0][0];
        virials[4] += virialMatrix[1][1];
        virials[8] += virialMatrix[2][2];
    }
}

void ListedGmxCalculator::compute(gmx::ArrayRef<const gmx::RVec> x,
                                  gmx::ArrayRef<const real>      charges,
                                  gmx::ArrayRef<gmx::RVec>       forces,
                                  gmx::ArrayRef<gmx::RVec>       shiftForces,
                                  ListedEnergies&                energies)
{
    compute(x, charges, forces, shiftForces, gmx::ArrayRef<real>{}, energies, {});
}

void ListedGmxCalculator::compute(gmx::ArrayRef<const gmx::RVec> x,
                                  gmx::ArrayRef<gmx::RVec>       forces,
                                  ListedEnergies&                energies)
{
    compute(x,
            gmx::ArrayRef<real>{},
            forces,
            gmx::ArrayRef<gmx::RVec>{},
            gmx::ArrayRef<real>{},
            energies,
            { 0, 0, 0 });
}

const InteractionDefinitions& ListedGmxCalculator::getIdef() const
{
    return *idef_;
}

} // namespace nblib
