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
 * \brief Defines the force element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "forceelement.h"

#include "gromacs/domdec/mdsetup.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdrun/shellfc.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/taskassignment/include/gromacs/taskassignment/decidesimulationworkload.h"
#include "gromacs/topology/topology.h"

#include "energydata.h"
#include "freeenergyperturbationdata.h"
#include "modularsimulator.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"

struct gmx_edsam;
struct gmx_enfrot;
struct gmx_multisim_t;
class history_t;
struct MDModulesNotifiers;

namespace gmx
{
ForceElement::ForceElement(StatePropagatorData*        statePropagatorData,
                           EnergyData*                 energyData,
                           FreeEnergyPerturbationData* freeEnergyPerturbationData,
                           bool                        isVerbose,
                           bool                        isDynamicBox,
                           FILE*                       fplog,
                           const t_commrec*            cr,
                           const t_inputrec*           inputrec,
                           const MDModulesNotifiers&   mdModulesNotifiers,
                           const MDAtoms*              mdAtoms,
                           t_nrnb*                     nrnb,
                           t_forcerec*                 fr,
                           gmx_wallcycle*              wcycle,
                           MdrunScheduleWorkload*      runScheduleWork,
                           VirtualSitesHandler*        vsite,
                           ImdSession*                 imdSession,
                           pull_t*                     pull_work,
                           Constraints*                constr,
                           const gmx_mtop_t&           globalTopology,
                           gmx_enfrot*                 enforcedRotation) :
    shellfc_(init_shell_flexcon(fplog,
                                globalTopology,
                                constr ? constr->numFlexibleConstraints() : 0,
                                inputrec->nstcalcenergy,
                                haveDDAtomOrdering(*cr),
                                runScheduleWork->simulationWork.useGpuPme)),
    doShellFC_(shellfc_ != nullptr),
    nextNSStep_(-1),
    nextEnergyCalculationStep_(-1),
    nextVirialCalculationStep_(-1),
    nextFreeEnergyCalculationStep_(-1),
    statePropagatorData_(statePropagatorData),
    energyData_(energyData),
    freeEnergyPerturbationData_(freeEnergyPerturbationData),
    localTopology_(nullptr),
    isDynamicBox_(isDynamicBox),
    isVerbose_(isVerbose),
    nShellRelaxationSteps_(0),
    ddBalanceRegionHandler_(cr),
    longRangeNonbondeds_(std::make_unique<CpuPpLongRangeNonbondeds>(fr->n_tpi,
                                                                    fr->ic->ewaldcoeff_q,
                                                                    fr->ic->epsilon_r,
                                                                    fr->qsum,
                                                                    fr->ic->eeltype,
                                                                    fr->ic->vdwtype,
                                                                    *inputrec,
                                                                    nrnb,
                                                                    wcycle,
                                                                    fplog)),
    lambda_(),
    fplog_(fplog),
    cr_(cr),
    inputrec_(inputrec),
    mdModulesNotifiers_(mdModulesNotifiers),
    mdAtoms_(mdAtoms),
    nrnb_(nrnb),
    wcycle_(wcycle),
    fr_(fr),
    vsite_(vsite),
    imdSession_(imdSession),
    pull_work_(pull_work),
    runScheduleWork_(runScheduleWork),
    constr_(constr),
    enforcedRotation_(enforcedRotation)
{
    std::fill(lambda_.begin(), lambda_.end(), 0);

    if (doShellFC_ && !haveDDAtomOrdering(*cr))
    {
        // This was done in mdAlgorithmsSetupAtomData(), but shellfc
        // won't be available outside this element.
        make_local_shells(cr, *mdAtoms->mdatoms(), shellfc_);
    }
}

ForceElement::~ForceElement() = default;

void ForceElement::scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction)
{
    unsigned int flags =
            (GMX_FORCE_STATECHANGED | GMX_FORCE_ALLFORCES | (isDynamicBox_ ? GMX_FORCE_DYNAMICBOX : 0)
             | (doShellFC_ && isVerbose_ ? GMX_FORCE_ENERGY : 0)
             | (nextVirialCalculationStep_ == step ? GMX_FORCE_VIRIAL : 0)
             | (nextEnergyCalculationStep_ == step ? GMX_FORCE_ENERGY : 0)
             | (nextFreeEnergyCalculationStep_ == step ? GMX_FORCE_DHDL : 0)
             | (nextNSStep_ == step ? GMX_FORCE_NS : 0));

    registerRunFunction([this, step, time, flags]() {
        if (doShellFC_)
        {
            run<true>(step, time, flags);
        }
        else
        {
            run<false>(step, time, flags);
        }
    });
}

void ForceElement::elementSetup()
{
    GMX_ASSERT(localTopology_, "Setup called before local topology was set.");
}

template<bool doShellFC>
void ForceElement::run(Step step, Time time, unsigned int flags)
{
    // Disabled functionality
    gmx_multisim_t* ms = nullptr;


    if (!haveDDAtomOrdering(*cr_) && (flags & GMX_FORCE_NS) && inputrecDynamicBox(inputrec_))
    {
        // TODO: Correcting the box is done in DomDecHelper (if using DD) or here (non-DD simulations).
        //       Think about unifying this responsibility, could this be done in one place?
        auto* box = statePropagatorData_->box();
        correct_box(fplog_, step, box);
    }

    if (flags & GMX_FORCE_NS)
    {
        if (fr_->listedForcesGpu)
        {
            fr_->listedForcesGpu->updateHaveInteractions(localTopology_->idef);
        }
        gmx_edsam* ed                = nullptr; // disabled
        runScheduleWork_->domainWork = setupDomainLifetimeWorkload(
                *inputrec_, *fr_, pull_work_, ed, *mdAtoms_->mdatoms(), runScheduleWork_->simulationWork);
    }

    runScheduleWork_->stepWork = setupStepWorkload(
            flags, inputrec_->mtsLevels, step, runScheduleWork_->domainWork, runScheduleWork_->simulationWork);

    /* The coordinates (x) are shifted (to get whole molecules)
     * in do_force.
     * This is parallelized as well, and does communication too.
     * Check comments in sim_util.c
     */
    auto        x      = statePropagatorData_->positionsView();
    auto&       forces = statePropagatorData_->forcesView();
    const auto* box    = statePropagatorData_->constBox();
    history_t*  hist   = nullptr; // disabled

    tensor force_vir = { { 0 } };
    // TODO: Make lambda const (needs some adjustments in lower force routines)
    ArrayRef<real> lambda =
            freeEnergyPerturbationData_ ? freeEnergyPerturbationData_->lambdaView() : lambda_;

    longRangeNonbondeds_->updateAfterPartition(*mdAtoms_->mdatoms());

    if (doShellFC)
    {
        auto v = statePropagatorData_->velocitiesView();

        relax_shell_flexcon(fplog_,
                            cr_,
                            ms,
                            isVerbose_,
                            enforcedRotation_,
                            step,
                            inputrec_,
                            mdModulesNotifiers_,
                            imdSession_,
                            pull_work_,
                            step == nextNSStep_,
                            localTopology_,
                            constr_,
                            energyData_->enerdata(),
                            statePropagatorData_->localNumAtoms(),
                            x,
                            v,
                            box,
                            lambda,
                            hist,
                            &forces,
                            force_vir,
                            *mdAtoms_->mdatoms(),
                            longRangeNonbondeds_.get(),
                            nrnb_,
                            wcycle_,
                            shellfc_,
                            fr_,
                            *runScheduleWork_,
                            time,
                            energyData_->muTot(),
                            vsite_,
                            ddBalanceRegionHandler_);
        nShellRelaxationSteps_++;
    }
    else
    {
        // Disabled functionality
        Awh*       awh = nullptr;
        gmx_edsam* ed  = nullptr;

        auto v = statePropagatorData_->velocitiesView();

        do_force(fplog_,
                 cr_,
                 ms,
                 *inputrec_,
                 mdModulesNotifiers_,
                 awh,
                 enforcedRotation_,
                 imdSession_,
                 pull_work_,
                 step,
                 nrnb_,
                 wcycle_,
                 localTopology_,
                 box,
                 x,
                 v.unpaddedArrayRef(),
                 hist,
                 &forces,
                 force_vir,
                 mdAtoms_->mdatoms(),
                 energyData_->enerdata(),
                 lambda,
                 fr_,
                 *runScheduleWork_,
                 vsite_,
                 energyData_->muTot(),
                 time,
                 ed,
                 longRangeNonbondeds_.get(),
                 ddBalanceRegionHandler_);
    }
    energyData_->addToForceVirial(force_vir, step);
}

void ForceElement::elementTeardown()
{
    if (doShellFC_)
    {
        done_shellfc(fplog_, shellfc_, nShellRelaxationSteps_);
    }
}

void ForceElement::setTopology(const gmx_localtop_t* top)
{
    localTopology_ = top;
}

std::optional<SignallerCallback> ForceElement::registerNSCallback()
{
    return [this](Step step, Time gmx_unused time) { this->nextNSStep_ = step; };
}

std::optional<SignallerCallback> ForceElement::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::EnergyCalculationStep)
    {
        return [this](Step step, Time /*unused*/) { nextEnergyCalculationStep_ = step; };
    }
    if (event == EnergySignallerEvent::VirialCalculationStep)
    {
        return [this](Step step, Time /*unused*/) { nextVirialCalculationStep_ = step; };
    }
    if (event == EnergySignallerEvent::FreeEnergyCalculationStep)
    {
        return [this](Step step, Time /*unused*/) { nextFreeEnergyCalculationStep_ = step; };
    }
    return std::nullopt;
}

DomDecCallback ForceElement::registerDomDecCallback()
{
    return [this]() { longRangeNonbondeds_->updateAfterPartition(*mdAtoms_->mdatoms()); };
}

ISimulatorElement*
ForceElement::getElementPointerImpl(LegacySimulatorData*                    legacySimulatorData,
                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                    StatePropagatorData*                    statePropagatorData,
                                    EnergyData*                             energyData,
                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                    GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
                                    ObservablesReducer* /*observablesReducer*/)
{
    const bool isVerbose    = legacySimulatorData->mdrunOptions_.verbose;
    const bool isDynamicBox = inputrecDynamicBox(legacySimulatorData->inputRec_);
    return builderHelper->storeElement(
            std::make_unique<ForceElement>(statePropagatorData,
                                           energyData,
                                           freeEnergyPerturbationData,
                                           isVerbose,
                                           isDynamicBox,
                                           legacySimulatorData->fpLog_,
                                           legacySimulatorData->cr_,
                                           legacySimulatorData->inputRec_,
                                           legacySimulatorData->mdModulesNotifiers_,
                                           legacySimulatorData->mdAtoms_,
                                           legacySimulatorData->nrnb_,
                                           legacySimulatorData->fr_,
                                           legacySimulatorData->wallCycleCounters_,
                                           legacySimulatorData->runScheduleWork_,
                                           legacySimulatorData->virtualSites_,
                                           legacySimulatorData->imdSession_,
                                           legacySimulatorData->pullWork_,
                                           legacySimulatorData->constr_,
                                           legacySimulatorData->topGlobal_,
                                           legacySimulatorData->enforcedRotation_));
}
} // namespace gmx
