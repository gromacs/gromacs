/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief Defines the shell / flex constraints element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "shellfcelement.h"

#include "gromacs/domdec/mdsetup.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdrun/shellfc.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"

#include "energyelement.h"
#include "freeenergyperturbationelement.h"
#include "statepropagatordata.h"

struct gmx_edsam;
struct gmx_enfrot;
struct gmx_multisim_t;
class history_t;
struct t_graph;

namespace gmx
{
bool ShellFCElement::doShellsOrFlexConstraints(const gmx_mtop_t& mtop, int nflexcon)
{
    if (nflexcon != 0)
    {
        return true;
    }
    std::array<int, eptNR> n = gmx_mtop_particletype_count(mtop);
    return n[eptShell] != 0;
}

ShellFCElement::ShellFCElement(StatePropagatorData*           statePropagatorData,
                               EnergyElement*                 energyElement,
                               FreeEnergyPerturbationElement* freeEnergyPerturbationElement,
                               bool                           isVerbose,
                               bool                           isDynamicBox,
                               FILE*                          fplog,
                               const t_commrec*               cr,
                               const t_inputrec*              inputrec,
                               const MDAtoms*                 mdAtoms,
                               t_nrnb*                        nrnb,
                               t_forcerec*                    fr,
                               t_fcdata*                      fcd,
                               gmx_wallcycle*                 wcycle,
                               MdrunScheduleWorkload*         runScheduleWork,
                               gmx_vsite_t*                   vsite,
                               ImdSession*                    imdSession,
                               pull_t*                        pull_work,
                               Constraints*                   constr,
                               const gmx_mtop_t*              globalTopology,
                               gmx_enfrot*                    enforcedRotation) :
    nextNSStep_(-1),
    nextEnergyCalculationStep_(-1),
    nextVirialCalculationStep_(-1),
    nextFreeEnergyCalculationStep_(-1),
    statePropagatorData_(statePropagatorData),
    energyElement_(energyElement),
    freeEnergyPerturbationElement_(freeEnergyPerturbationElement),
    localTopology_(nullptr),
    isDynamicBox_(isDynamicBox),
    isVerbose_(isVerbose),
    nSteps_(0),
    ddBalanceRegionHandler_(cr),
    fplog_(fplog),
    cr_(cr),
    inputrec_(inputrec),
    mdAtoms_(mdAtoms),
    nrnb_(nrnb),
    wcycle_(wcycle),
    fr_(fr),
    vsite_(vsite),
    imdSession_(imdSession),
    pull_work_(pull_work),
    fcd_(fcd),
    runScheduleWork_(runScheduleWork),
    constr_(constr),
    enforcedRotation_(enforcedRotation)
{
    lambda_.fill(0);

    shellfc_ = init_shell_flexcon(fplog, globalTopology, constr_ ? constr_->numFlexibleConstraints() : 0,
                                  inputrec->nstcalcenergy, DOMAINDECOMP(cr));

    GMX_ASSERT(shellfc_, "ShellFCElement built, but init_shell_flexcon returned a nullptr");

    if (!DOMAINDECOMP(cr))
    {
        // This was done in mdAlgorithmsSetupAtomData(), but shellfc
        // won't be available outside this element.
        make_local_shells(cr, mdAtoms->mdatoms(), shellfc_);
    }
}

void ShellFCElement::scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction)
{
    unsigned int flags =
            (GMX_FORCE_STATECHANGED | GMX_FORCE_ALLFORCES | (isDynamicBox_ ? GMX_FORCE_DYNAMICBOX : 0)
             | (nextVirialCalculationStep_ == step ? GMX_FORCE_VIRIAL : 0)
             | (nextEnergyCalculationStep_ == step ? GMX_FORCE_ENERGY : 0)
             | (nextFreeEnergyCalculationStep_ == step ? GMX_FORCE_DHDL : 0));

    const bool isNSStep = (step == nextNSStep_);
    (*registerRunFunction)(std::make_unique<SimulatorRunFunction>(
            [this, step, time, flags, isNSStep]() { run(step, time, isNSStep, flags); }));
    nSteps_++;
}

void ShellFCElement::elementSetup()
{
    GMX_ASSERT(localTopology_, "Setup called before local topology was set.");
}

void ShellFCElement::run(Step step, Time time, bool isNSStep, unsigned int flags)
{
    // Disabled functionality
    gmx_multisim_t* ms    = nullptr;
    t_graph*        graph = nullptr;

    if (!DOMAINDECOMP(cr_) && isNSStep && inputrecDynamicBox(inputrec_))
    {
        // TODO: Correcting the box is done in DomDecHelper (if using DD) or here (non-DD simulations).
        //       Think about unifying this responsibility, could this be done in one place?
        auto box = statePropagatorData_->box();
        correct_box(fplog_, step, box, graph);
    }

    auto       x      = statePropagatorData_->positionsView();
    auto       v      = statePropagatorData_->velocitiesView();
    auto       forces = statePropagatorData_->forcesView();
    auto       box    = statePropagatorData_->constBox();
    history_t* hist   = nullptr; // disabled

    tensor force_vir = { { 0 } };
    // TODO: Make lambda const (needs some adjustments in lower force routines)
    ArrayRef<real> lambda =
            freeEnergyPerturbationElement_ ? freeEnergyPerturbationElement_->lambdaView() : lambda_;
    relax_shell_flexcon(fplog_, cr_, ms, isVerbose_, enforcedRotation_, step, inputrec_, imdSession_,
                        pull_work_, isNSStep, static_cast<int>(flags), localTopology_, constr_,
                        energyElement_->enerdata(), fcd_, statePropagatorData_->localNumAtoms(), x,
                        v, box, lambda, hist, forces, force_vir, mdAtoms_->mdatoms(), nrnb_,
                        wcycle_, graph, shellfc_, fr_, runScheduleWork_, time,
                        energyElement_->muTot(), vsite_, ddBalanceRegionHandler_);
    energyElement_->addToForceVirial(force_vir, step);
}

void ShellFCElement::elementTeardown()
{
    done_shellfc(fplog_, shellfc_, nSteps_);
}

void ShellFCElement::setTopology(const gmx_localtop_t* top)
{
    localTopology_ = top;
}

SignallerCallbackPtr ShellFCElement::registerNSCallback()
{
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->nextNSStep_ = step; });
}

SignallerCallbackPtr ShellFCElement::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::EnergyCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time /*unused*/) { nextEnergyCalculationStep_ = step; });
    }
    if (event == EnergySignallerEvent::VirialCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time /*unused*/) { nextVirialCalculationStep_ = step; });
    }
    if (event == EnergySignallerEvent::FreeEnergyCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time /*unused*/) { nextFreeEnergyCalculationStep_ = step; });
    }
    return nullptr;
}
} // namespace gmx
