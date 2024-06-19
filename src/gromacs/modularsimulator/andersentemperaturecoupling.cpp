/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief Defines Andersen temperature coupling for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "andersentemperaturecoupling.h"

#include <cmath>

#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "gromacs/compat/pointers.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdrun/isimulator.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/modularsimulator/modularsimulatorinterfaces.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

#include "compositesimulatorelement.h"
#include "constraintelement.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"

namespace gmx
{
class EnergyData;
class FreeEnergyPerturbationData;
enum class ReferenceTemperatureChangeAlgorithm;

AndersenTemperatureCoupling::AndersenTemperatureCoupling(double               simulationTimestep,
                                                         bool                 doMassive,
                                                         int64_t              seed,
                                                         ArrayRef<const real> referenceTemperature,
                                                         ArrayRef<const real> couplingTime,
                                                         StatePropagatorData* statePropagatorData,
                                                         const MDAtoms*       mdAtoms,
                                                         const t_commrec*     cr) :
    doMassive_(doMassive),
    randomizationRate_(simulationTimestep / couplingTime[0]),
    couplingFrequency_(doMassive ? roundToInt(1. / randomizationRate_) : 1),
    seed_(seed),
    referenceTemperature_(referenceTemperature),
    couplingTime_(couplingTime),
    statePropagatorData_(statePropagatorData),
    mdAtoms_(mdAtoms->mdatoms()),
    cr_(cr)
{
}

void AndersenTemperatureCoupling::scheduleTask(Step                       step,
                                               Time gmx_unused            time,
                                               const RegisterRunFunction& registerRunFunction)
{
    if (do_per_step(step, couplingFrequency_))
    {
        registerRunFunction([this, step]() { apply(step); });
    }
}

void AndersenTemperatureCoupling::apply(Step step)
{
    ThreeFry2x64<0>                       rng(seed_, RandomDomain::Thermostat);
    UniformRealDistribution<real>         uniformDist;
    TabulatedNormalDistribution<real, 14> normalDist;

    const bool atomOrderingIsDD = haveDDAtomOrdering(*cr_);

    auto velocities = statePropagatorData_->velocitiesView().unpaddedArrayRef();

    for (int atomIdx = 0; atomIdx < mdAtoms_->homenr; ++atomIdx)
    {
        const int temperatureGroup = !mdAtoms_->cTC.empty() ? mdAtoms_->cTC[atomIdx] : 0;
        if (referenceTemperature_[temperatureGroup] <= 0 || couplingTime_[temperatureGroup] <= 0)
        {
            continue;
        }

        const int globalAtomIdx = atomOrderingIsDD ? cr_->dd->globalAtomIndices[atomIdx] : atomIdx;
        rng.restart(step, globalAtomIdx);

        // For massive Andersen, this function is only called periodically, but we apply each time
        // Otherwise, this function is called every step, but we randomize atoms probabilistically
        if (!doMassive_)
        {
            uniformDist.reset();
        }
        if (doMassive_ || (uniformDist(rng) < randomizationRate_))
        {
            const real scalingFactor = std::sqrt(c_boltz * referenceTemperature_[temperatureGroup]
                                                 * mdAtoms_->invmass[atomIdx]);
            normalDist.reset();
            for (int d = 0; d < DIM; d++)
            {
                velocities[atomIdx][d] = scalingFactor * normalDist(rng);
            }
        }
    }
}

int AndersenTemperatureCoupling::frequency() const
{
    return couplingFrequency_;
}

void AndersenTemperatureCoupling::updateReferenceTemperature(ArrayRef<const real> gmx_unused temperatures,
                                                             ReferenceTemperatureChangeAlgorithm gmx_unused algorithm)
{
    // Currently, we don't know about any temperature change algorithms, so we assert this never gets called
    GMX_ASSERT(false, "AndersenTemperatureCoupling: Unknown ReferenceTemperatureChangeAlgorithm.");
}

void               AndersenTemperatureCoupling::elementSetup() {}
ISimulatorElement* AndersenTemperatureCoupling::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData*                    statePropagatorData,
        EnergyData*                             energyData,
        FreeEnergyPerturbationData*             freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer gmx_unused* observablesReducer)
{
    GMX_RELEASE_ASSERT(legacySimulatorData->inputRec_->etc == TemperatureCoupling::Andersen
                               || legacySimulatorData->inputRec_->etc == TemperatureCoupling::AndersenMassive,
                       "Expected the thermostat type to be andersen or andersen-massive.");
    auto andersenThermostat = std::make_unique<AndersenTemperatureCoupling>(
            legacySimulatorData->inputRec_->delta_t,
            legacySimulatorData->inputRec_->etc == TemperatureCoupling::AndersenMassive,
            legacySimulatorData->inputRec_->andersen_seed,
            constArrayRefFromArray(legacySimulatorData->inputRec_->opts.ref_t,
                                   legacySimulatorData->inputRec_->opts.ngtc),
            constArrayRefFromArray(legacySimulatorData->inputRec_->opts.tau_t,
                                   legacySimulatorData->inputRec_->opts.ngtc),
            statePropagatorData,
            legacySimulatorData->mdAtoms_,
            legacySimulatorData->cr_);
    auto* andersenThermostatPtr = andersenThermostat.get();
    builderHelper->registerReferenceTemperatureUpdate(
            [andersenThermostatPtr](ArrayRef<const real>                temperatures,
                                    ReferenceTemperatureChangeAlgorithm algorithm) {
                andersenThermostatPtr->updateReferenceTemperature(temperatures, algorithm);
            });

    // T-coupling frequency will be composite element frequency
    const auto frequency = andersenThermostat->frequency();
    // Set up call list for composite element
    std::vector<compat::not_null<ISimulatorElement*>> elementCallList = { compat::make_not_null(
            andersenThermostat.get()) };
    // Set up element list for composite element
    std::vector<std::unique_ptr<gmx::ISimulatorElement>> elements;
    elements.emplace_back(std::move(andersenThermostat));

    // If there are constraints, add constraint element after Andersen element
    if (legacySimulatorData->constr_)
    {
        // This is excluded in preprocessing -
        // asserted here to make sure things don't get out of sync
        GMX_RELEASE_ASSERT(
                legacySimulatorData->inputRec_->etc == TemperatureCoupling::AndersenMassive,
                "Per-particle Andersen thermostat is not implemented for systems with constrains.");
        // Build constraint element
        auto constraintElement = std::make_unique<ConstraintsElement<ConstraintVariable::Velocities>>(
                legacySimulatorData->constr_,
                statePropagatorData,
                energyData,
                freeEnergyPerturbationData,
                MAIN(legacySimulatorData->cr_),
                legacySimulatorData->fpLog_,
                legacySimulatorData->inputRec_,
                legacySimulatorData->mdAtoms_->mdatoms());
        // Add call to composite element call list
        elementCallList.emplace_back(compat::make_not_null(constraintElement.get()));
        // Move ownership of constraint element to composite element
        elements.emplace_back(std::move(constraintElement));
    }

    // Store composite element in builder helper and return pointer
    return builderHelper->storeElement(std::make_unique<CompositeSimulatorElement>(
            std::move(elementCallList), std::move(elements), frequency));
}

} // namespace gmx
