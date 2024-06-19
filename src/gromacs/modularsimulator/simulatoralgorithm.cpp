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
 * \brief Defines the modular simulator algorithm
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "simulatoralgorithm.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <iterator>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdrun/replicaexchange.h"
#include "gromacs/mdrun/shellfc.h"
#include "gromacs/mdrunutility/freeenergy.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/modularsimulator/modularsimulatorinterfaces.h"
#include "gromacs/modularsimulator/signallers.h"
#include "gromacs/modularsimulator/topologyholder.h"
#include "gromacs/modularsimulator/trajectoryelement.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

#include "checkpointhelper.h"
#include "domdechelper.h"
#include "energydata.h"
#include "firstorderpressurecoupling.h"
#include "freeenergyperturbationdata.h"
#include "modularsimulator.h"
#include "pmeloadbalancehelper.h"
#include "propagator.h"
#include "referencetemperaturemanager.h"
#include "statepropagatordata.h"

struct gmx_walltime_accounting;

namespace gmx
{
ModularSimulatorAlgorithm::ModularSimulatorAlgorithm(std::string              topologyName,
                                                     FILE*                    fplog,
                                                     t_commrec*               cr,
                                                     const MDLogger&          mdlog,
                                                     const MdrunOptions&      mdrunOptions,
                                                     const t_inputrec*        inputrec,
                                                     t_nrnb*                  nrnb,
                                                     gmx_wallcycle*           wcycle,
                                                     t_forcerec*              fr,
                                                     gmx_walltime_accounting* walltime_accounting) :
    taskIterator_(taskQueue_.end()),
    statePropagatorData_(nullptr),
    energyData_(nullptr),
    freeEnergyPerturbationData_(nullptr),
    step_(-1),
    runFinished_(false),
    topologyName_(std::move(topologyName)),
    fpLog_(fplog),
    cr_(cr),
    mdLog_(mdlog),
    mdrunOptions_(mdrunOptions),
    inputRec_(inputrec),
    nrnb_(nrnb),
    wallCycle_(wcycle),
    fr_(fr),
    wallTimeAccounting_(walltime_accounting)
{
    signalHelper_ = std::make_unique<SignalHelper>();
}

void ModularSimulatorAlgorithm::setup()
{
    simulatorSetup();
    for (auto& signaller : signallerList_)
    {
        signaller->setup();
    }
    if (domDecHelper_)
    {
        domDecHelper_->setup();
    }

    for (auto& element : elementSetupTeardownList_)
    {
        element->elementSetup();
    }
    statePropagatorData_->setup();
    if (pmeLoadBalanceHelper_)
    {
        // State must have been initialized so pmeLoadBalanceHelper_ gets a valid box
        pmeLoadBalanceHelper_->setup();
    }
}

const SimulatorRunFunction* ModularSimulatorAlgorithm::getNextTask()
{
    if (!taskQueue_.empty())
    {
        ++taskIterator_;
    }
    if (taskIterator_ == taskQueue_.end())
    {
        if (runFinished_)
        {
            return nullptr;
        }
        updateTaskQueue();
        taskIterator_ = taskQueue_.begin();
    }
    return &*taskIterator_;
}

void ModularSimulatorAlgorithm::updateTaskQueue()
{
    // For now, we'll just clean the task queue and then re-populate
    // TODO: If tasks are periodic around updates of the task queue,
    //       we should reuse it instead
    taskQueue_.clear();
    populateTaskQueue();
}

void ModularSimulatorAlgorithm::teardown()
{
    for (auto& element : elementSetupTeardownList_)
    {
        element->elementTeardown();
    }
    energyData_->teardown();
    if (pmeLoadBalanceHelper_)
    {
        pmeLoadBalanceHelper_->teardown();
    }
    simulatorTeardown();
}

void ModularSimulatorAlgorithm::simulatorSetup()
{
    if (!mdrunOptions_.writeConfout)
    {
        // This is on by default, and the main known use case for
        // turning it off is for convenience in benchmarking, which is
        // something that should not show up in the general user
        // interface.
        GMX_LOG(mdLog_.info)
                .asParagraph()
                .appendText(
                        "The -noconfout functionality is deprecated, and "
                        "may be removed in a future version.");
    }

    if (MAIN(cr_))
    {
        char        sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
        std::string timeString;
        fprintf(stderr, "starting mdrun '%s'\n", topologyName_.c_str());
        if (inputRec_->nsteps >= 0)
        {
            timeString = formatString("%8.1f",
                                      static_cast<double>(inputRec_->init_step + inputRec_->nsteps)
                                              * inputRec_->delta_t);
        }
        else
        {
            timeString = "infinite";
        }
        if (inputRec_->init_step > 0)
        {
            fprintf(stderr,
                    "%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                    gmx_step_str(inputRec_->init_step + inputRec_->nsteps, sbuf),
                    timeString.c_str(),
                    gmx_step_str(inputRec_->init_step, sbuf2),
                    inputRec_->init_step * inputRec_->delta_t);
        }
        else
        {
            fprintf(stderr, "%s steps, %s ps.\n", gmx_step_str(inputRec_->nsteps, sbuf), timeString.c_str());
        }
        fprintf(fpLog_, "\n");
    }

    walltime_accounting_start_time(wallTimeAccounting_);
    wallcycle_start(wallCycle_, WallCycleCounter::Run);
    print_start(fpLog_, cr_, wallTimeAccounting_, "mdrun");

    step_ = inputRec_->init_step;
}

void ModularSimulatorAlgorithm::simulatorTeardown()
{

    // Stop measuring walltime
    walltime_accounting_end_time(wallTimeAccounting_);

    if (!thisRankHasDuty(cr_, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr_);
    }

    walltime_accounting_set_nsteps_done(wallTimeAccounting_, step_ - inputRec_->init_step);
}

void ModularSimulatorAlgorithm::preStep(Step step, Time gmx_unused time, bool isNeighborSearchingStep)
{
    if (stopHandler_->stoppingAfterCurrentStep(isNeighborSearchingStep) && step != signalHelper_->lastStep_)
    {
        /*
         * Stop handler wants to stop after the current step, which was
         * not known when building the current task queue. This happens
         * e.g. when a stop is signalled by OS. We therefore want to purge
         * the task queue now, and re-schedule this step as last step.
         */
        // clear task queue
        taskQueue_.clear();
        // rewind step
        step_ = step;
        return;
    }

    resetHandler_->setSignal(wallTimeAccounting_);
    // This is a hack to avoid having to rewrite StopHandler to be a NeighborSearchSignaller
    // and accept the step as input. Eventually, we want to do that, but currently this would
    // require introducing NeighborSearchSignaller in the legacy do_md or a lot of code
    // duplication.
    stophandlerIsNSStep_    = isNeighborSearchingStep;
    stophandlerCurrentStep_ = step;
    stopHandler_->setSignal();

    wallcycle_start(wallCycle_, WallCycleCounter::Step);
}

void ModularSimulatorAlgorithm::postStep(Step step, Time gmx_unused time)
{
    // Output stuff
    if (MAIN(cr_))
    {
        if (do_per_step(step, inputRec_->nstlog))
        {
            if (fflush(fpLog_) != 0)
            {
                gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
            }
        }
    }
    const bool do_verbose = mdrunOptions_.verbose
                            && (step % mdrunOptions_.verboseStepPrintInterval == 0
                                || step == inputRec_->init_step || step == signalHelper_->lastStep_);
    // Print the remaining wall clock time for the run
    if (MAIN(cr_) && (do_verbose || gmx_got_usr_signal())
        && !(pmeLoadBalanceHelper_ && pmeLoadBalanceHelper_->pmePrinting()))
    {
        print_time(stderr, wallTimeAccounting_, step, inputRec_, cr_);
    }

    double cycles = wallcycle_stop(wallCycle_, WallCycleCounter::Step);
    if (haveDDAtomOrdering(*cr_) && wallCycle_)
    {
        dd_cycles_add(cr_->dd, static_cast<float>(cycles), ddCyclStep);
    }

    resetHandler_->resetCounters(step,
                                 step - inputRec_->init_step,
                                 mdLog_,
                                 fpLog_,
                                 cr_,
                                 fr_->nbv.get(),
                                 nrnb_,
                                 fr_->pmedata,
                                 pmeLoadBalanceHelper_ ? pmeLoadBalanceHelper_->loadBalancingObject() : nullptr,
                                 wallCycle_,
                                 wallTimeAccounting_);
}

void ModularSimulatorAlgorithm::populateTaskQueue()
{
    /*
     * The registerRunFunction emplaces functions to the task queue.
     * All elements are owned by the ModularSimulatorAlgorithm, as is the task queue.
     * Elements can hence register lambdas capturing their `this` pointers without expecting
     * life time issues, as the task queue and the elements are in the same scope.
     */
    auto registerRunFunction = [this](SimulatorRunFunction function) {
        taskQueue_.emplace_back(std::move(function));
    };

    Time startTime = inputRec_->init_t;
    Time timeStep  = inputRec_->delta_t;
    Time time      = startTime + step_ * timeStep;

    // Run an initial call to the signallers
    for (auto& signaller : signallerList_)
    {
        signaller->signal(step_, time);
    }

    if (checkpointHelper_)
    {
        checkpointHelper_->run(step_, time);
    }

    if (pmeLoadBalanceHelper_)
    {
        pmeLoadBalanceHelper_->run(step_, time);
    }
    if (domDecHelper_)
    {
        domDecHelper_->run(step_, time);
    }

    do
    {
        // local variables for lambda capturing
        const int  step     = step_;
        const bool isNSStep = step == signalHelper_->nextNSStep_;

        // register pre-step (task queue is local, so no problem with `this`)
        registerRunFunction([this, step, time, isNSStep]() { preStep(step, time, isNSStep); });
        // register pre step functions
        for (const auto& schedulingFunction : preStepScheduling_)
        {
            schedulingFunction(step_, time, registerRunFunction);
        }
        // register elements for step
        for (auto& element : elementCallList_)
        {
            element->scheduleTask(step_, time, registerRunFunction);
        }
        // register post step functions
        for (const auto& schedulingFunction : postStepScheduling_)
        {
            schedulingFunction(step_, time, registerRunFunction);
        }
        // register post-step (task queue is local, so no problem with `this`)
        registerRunFunction([this, step, time]() { postStep(step, time); });

        // prepare next step
        step_++;
        time = startTime + step_ * timeStep;
        for (auto& signaller : signallerList_)
        {
            signaller->signal(step_, time);
        }
    } while (step_ != signalHelper_->nextNSStep_ && step_ <= signalHelper_->lastStep_);

    runFinished_ = (step_ > signalHelper_->lastStep_);

    if (runFinished_)
    {
        // task queue is local, so no problem with `this`
        registerRunFunction([this]() { teardown(); });
    }
}

ModularSimulatorAlgorithmBuilder::ModularSimulatorAlgorithmBuilder(
        compat::not_null<LegacySimulatorData*>    legacySimulatorData,
        std::unique_ptr<ReadCheckpointDataHolder> checkpointDataHolder) :
    legacySimulatorData_(legacySimulatorData),
    signals_(std::make_unique<SimulationSignals>()),
    elementAdditionHelper_(this),
    globalCommunicationHelper_(computeGlobalCommunicationPeriod(legacySimulatorData->mdLog_,
                                                                legacySimulatorData->inputRec_,
                                                                legacySimulatorData->cr_),
                               signals_.get()),
    observablesReducer_(legacySimulatorData->observablesReducerBuilder_->build()),
    checkpointHelperBuilder_(std::move(checkpointDataHolder),
                             legacySimulatorData->startingBehavior_,
                             legacySimulatorData->cr_)
{
    if (legacySimulatorData->inputRec_->efep != FreeEnergyPerturbationType::No)
    {
        freeEnergyPerturbationData_ =
                std::make_unique<FreeEnergyPerturbationData>(legacySimulatorData->fpLog_,
                                                             *legacySimulatorData->inputRec_,
                                                             legacySimulatorData->mdAtoms_,
                                                             legacySimulatorData->ekind_);
        registerExistingElement(freeEnergyPerturbationData_->element());
    }

    statePropagatorData_ = std::make_unique<StatePropagatorData>(
            legacySimulatorData->topGlobal_.natoms,
            legacySimulatorData->fpLog_,
            legacySimulatorData->cr_,
            legacySimulatorData->stateGlobal_,
            legacySimulatorData->state_,
            legacySimulatorData->fr_->nbv->useGpu(),
            legacySimulatorData->fr_->bMolPBC,
            legacySimulatorData->mdrunOptions_.writeConfout,
            opt2fn("-c", legacySimulatorData->nFile_, legacySimulatorData->fnm_),
            legacySimulatorData->inputRec_,
            legacySimulatorData->mdAtoms_->mdatoms(),
            legacySimulatorData->topGlobal_);
    registerExistingElement(statePropagatorData_->element());

    // Multi sim is turned off
    const bool simulationsShareHamiltonian = false;

    energyData_ = std::make_unique<EnergyData>(statePropagatorData_.get(),
                                               freeEnergyPerturbationData_.get(),
                                               legacySimulatorData->topGlobal_,
                                               legacySimulatorData->inputRec_,
                                               legacySimulatorData->mdAtoms_,
                                               legacySimulatorData->enerd_,
                                               legacySimulatorData->ekind_,
                                               legacySimulatorData->constr_,
                                               legacySimulatorData->fpLog_,
                                               legacySimulatorData->fr_->fcdata.get(),
                                               legacySimulatorData->mdModulesNotifiers_,
                                               MAIN(legacySimulatorData->cr_),
                                               legacySimulatorData->observablesHistory_,
                                               legacySimulatorData->startingBehavior_,
                                               simulationsShareHamiltonian,
                                               legacySimulatorData->pullWork_);
    registerExistingElement(energyData_->element());

    storeSimulationData("ReferenceTemperatureManager",
                        ReferenceTemperatureManager(legacySimulatorData->ekind_));
    auto* referenceTemperatureManager =
            simulationData<ReferenceTemperatureManager>("ReferenceTemperatureManager").value();

    // State propagator data is scaling velocities if reference temperature is updated
    auto* statePropagatorDataPtr = statePropagatorData_.get();
    referenceTemperatureManager->registerUpdateCallback(
            [statePropagatorDataPtr](ArrayRef<const real>                temperatures,
                                     ReferenceTemperatureChangeAlgorithm algorithm) {
                statePropagatorDataPtr->updateReferenceTemperature(temperatures, algorithm);
            });
}

ModularSimulatorAlgorithm ModularSimulatorAlgorithmBuilder::build()
{
    if (algorithmHasBeenBuilt_)
    {
        GMX_THROW(SimulationAlgorithmSetupError(
                "Tried to build ModularSimulationAlgorithm more than once."));
    }
    algorithmHasBeenBuilt_ = true;

    // Connect propagators with thermostat / barostat
    for (const auto& registrationFunction : pressureTemperatureControlRegistrationFunctions_)
    {
        for (const auto& connection : propagatorConnections_)
        {
            registrationFunction(connection);
        }
    }

    ModularSimulatorAlgorithm algorithm(*(legacySimulatorData_->topGlobal_.name),
                                        legacySimulatorData_->fpLog_,
                                        legacySimulatorData_->cr_,
                                        legacySimulatorData_->mdLog_,
                                        legacySimulatorData_->mdrunOptions_,
                                        legacySimulatorData_->inputRec_,
                                        legacySimulatorData_->nrnb_,
                                        legacySimulatorData_->wallCycleCounters_,
                                        legacySimulatorData_->fr_,
                                        legacySimulatorData_->wallTimeAccounting_);
    registerWithInfrastructureAndSignallers(algorithm.signalHelper_.get());
    algorithm.statePropagatorData_        = std::move(statePropagatorData_);
    algorithm.energyData_                 = std::move(energyData_);
    algorithm.freeEnergyPerturbationData_ = std::move(freeEnergyPerturbationData_);
    algorithm.signals_                    = std::move(signals_);
    algorithm.simulationData_             = std::move(simulationData_);

    // Multi sim is turned off
    const bool simulationsShareState = false;

    // Build stop handler
    algorithm.stopHandler_ = legacySimulatorData_->stopHandlerBuilder_->getStopHandlerMD(
            compat::not_null<SimulationSignal*>(
                    &(*globalCommunicationHelper_.simulationSignals())[eglsSTOPCOND]),
            simulationsShareState,
            MAIN(legacySimulatorData_->cr_),
            legacySimulatorData_->inputRec_->nstlist,
            legacySimulatorData_->mdrunOptions_.reproducible,
            globalCommunicationHelper_.nstglobalcomm(),
            legacySimulatorData_->mdrunOptions_.maximumHoursToRun,
            legacySimulatorData_->inputRec_->nstlist == 0,
            legacySimulatorData_->fpLog_,
            algorithm.stophandlerCurrentStep_,
            algorithm.stophandlerIsNSStep_,
            legacySimulatorData_->wallTimeAccounting_);

    // Build reset handler
    const bool simulationsShareResetCounters = false;
    algorithm.resetHandler_                  = std::make_unique<ResetHandler>(
            compat::make_not_null<SimulationSignal*>(
                    &(*globalCommunicationHelper_.simulationSignals())[eglsRESETCOUNTERS]),
            simulationsShareResetCounters,
            legacySimulatorData_->inputRec_->nsteps,
            MAIN(legacySimulatorData_->cr_),
            legacySimulatorData_->mdrunOptions_.timingOptions.resetHalfway,
            legacySimulatorData_->mdrunOptions_.maximumHoursToRun,
            legacySimulatorData_->mdLog_,
            legacySimulatorData_->wallCycleCounters_,
            legacySimulatorData_->wallTimeAccounting_);

    // Build topology holder
    algorithm.topologyHolder_ = topologyHolderBuilder_.build(legacySimulatorData_->topGlobal_,
                                                             legacySimulatorData_->top_,
                                                             legacySimulatorData_->cr_,
                                                             legacySimulatorData_->inputRec_,
                                                             legacySimulatorData_->fr_,
                                                             legacySimulatorData_->mdAtoms_,
                                                             legacySimulatorData_->constr_,
                                                             legacySimulatorData_->virtualSites_);
    registerWithInfrastructureAndSignallers(algorithm.topologyHolder_.get());

    // Build PME load balance helper
    if (PmeLoadBalanceHelper::doPmeLoadBalancing(legacySimulatorData_->mdrunOptions_,
                                                 legacySimulatorData_->inputRec_,
                                                 legacySimulatorData_->fr_,
                                                 legacySimulatorData_->runScheduleWork_->simulationWork))
    {
        algorithm.pmeLoadBalanceHelper_ =
                std::make_unique<PmeLoadBalanceHelper>(legacySimulatorData_->mdrunOptions_.verbose,
                                                       algorithm.statePropagatorData_.get(),
                                                       legacySimulatorData_->fpLog_,
                                                       legacySimulatorData_->cr_,
                                                       legacySimulatorData_->mdLog_,
                                                       legacySimulatorData_->inputRec_,
                                                       legacySimulatorData_->wallCycleCounters_,
                                                       legacySimulatorData_->fr_);
        registerWithInfrastructureAndSignallers(algorithm.pmeLoadBalanceHelper_.get());
    }

    // Build trajectory element
    auto trajectoryElement = trajectoryElementBuilder_.build(legacySimulatorData_->fpLog_,
                                                             legacySimulatorData_->nFile_,
                                                             legacySimulatorData_->fnm_,
                                                             legacySimulatorData_->mdrunOptions_,
                                                             legacySimulatorData_->cr_,
                                                             legacySimulatorData_->outputProvider_,
                                                             legacySimulatorData_->mdModulesNotifiers_,
                                                             legacySimulatorData_->inputRec_,
                                                             legacySimulatorData_->topGlobal_,
                                                             legacySimulatorData_->oenv_,
                                                             legacySimulatorData_->wallCycleCounters_,
                                                             legacySimulatorData_->startingBehavior_,
                                                             simulationsShareState);
    registerWithInfrastructureAndSignallers(trajectoryElement.get());

    // Build domdec helper (free energy element is a client, so keep this after it is built)
    if (haveDDAtomOrdering(*legacySimulatorData_->cr_))
    {
        algorithm.domDecHelper_ =
                domDecHelperBuilder_.build(legacySimulatorData_->mdrunOptions_.verbose,
                                           legacySimulatorData_->mdrunOptions_.verboseStepPrintInterval,
                                           algorithm.statePropagatorData_.get(),
                                           algorithm.topologyHolder_.get(),
                                           legacySimulatorData_->fpLog_,
                                           legacySimulatorData_->cr_,
                                           legacySimulatorData_->mdLog_,
                                           legacySimulatorData_->constr_,
                                           legacySimulatorData_->inputRec_,
                                           legacySimulatorData_->mdModulesNotifiers_,
                                           legacySimulatorData_->mdAtoms_,
                                           legacySimulatorData_->nrnb_,
                                           legacySimulatorData_->wallCycleCounters_,
                                           legacySimulatorData_->fr_,
                                           legacySimulatorData_->virtualSites_,
                                           legacySimulatorData_->imdSession_,
                                           legacySimulatorData_->pullWork_);
        registerWithInfrastructureAndSignallers(algorithm.domDecHelper_.get());
    }
    // Build checkpoint helper (do this last so everyone else can be a checkpoint client!)
    {
        checkpointHelperBuilder_.setCheckpointHandler(std::make_unique<CheckpointHandler>(
                compat::make_not_null<SimulationSignal*>(&(*algorithm.signals_)[eglsCHKPT]),
                simulationsShareState,
                legacySimulatorData_->inputRec_->nstlist == 0,
                MAIN(legacySimulatorData_->cr_),
                legacySimulatorData_->mdrunOptions_.writeConfout,
                legacySimulatorData_->mdrunOptions_.checkpointOptions.period));
        algorithm.checkpointHelper_ =
                checkpointHelperBuilder_.build(legacySimulatorData_->inputRec_->init_step,
                                               trajectoryElement.get(),
                                               legacySimulatorData_->fpLog_,
                                               legacySimulatorData_->cr_,
                                               legacySimulatorData_->observablesHistory_,
                                               legacySimulatorData_->wallTimeAccounting_,
                                               legacySimulatorData_->stateGlobal_,
                                               legacySimulatorData_->mdrunOptions_.writeConfout);
        registerWithInfrastructureAndSignallers(algorithm.checkpointHelper_.get());
    }

    // Build signallers
    {
        /* Signallers need to be called in an exact order. Some signallers are clients
         * of other signallers, which requires the clients signallers to be called
         * _after_ any signaller they are registered to - otherwise, they couldn't
         * adapt their behavior to the information they got signalled.
         *
         * Signallers being clients of other signallers require registration.
         * That registration happens during construction, which in turn means that
         * we want to construct the signallers in the reverse order of their later
         * call order.
         *
         * For the above reasons, the `addSignaller` lambda defined below emplaces
         * added signallers at the beginning of the signaller list, which will yield
         * a signaller list which is inverse to the build order (and hence equal to
         * the intended call order).
         */
        auto addSignaller = [this, &algorithm](auto signaller) {
            registerWithInfrastructureAndSignallers(signaller.get());
            algorithm.signallerList_.emplace(algorithm.signallerList_.begin(), std::move(signaller));
        };
        const auto* inputrec   = legacySimulatorData_->inputRec_;
        auto        virialMode = EnergySignallerVirialMode::Off;
        if (inputrec->pressureCouplingOptions.epc != PressureCoupling::No)
        {
            if (EI_VV(inputrec->eI))
            {
                virialMode = EnergySignallerVirialMode::OnStepAndNext;
            }
            else
            {
                virialMode = EnergySignallerVirialMode::OnStep;
            }
        }
        addSignaller(energySignallerBuilder_.build(
                inputrec->nstcalcenergy,
                computeFepPeriod(*inputrec, legacySimulatorData_->replExParams_),
                inputrec->pressureCouplingOptions.nstpcouple,
                virialMode));
        addSignaller(trajectorySignallerBuilder_.build(inputrec->nstxout,
                                                       inputrec->nstvout,
                                                       inputrec->nstfout,
                                                       inputrec->nstxout_compressed,
                                                       trajectoryElement->tngBoxOut(),
                                                       trajectoryElement->tngLambdaOut(),
                                                       trajectoryElement->tngBoxOutCompressed(),
                                                       trajectoryElement->tngLambdaOutCompressed(),
                                                       inputrec->nstenergy));
        addSignaller(loggingSignallerBuilder_.build(
                inputrec->nstlog, inputrec->init_step, legacySimulatorData_->startingBehavior_));
        addSignaller(lastStepSignallerBuilder_.build(
                inputrec->nsteps, inputrec->init_step, algorithm.stopHandler_.get()));
        addSignaller(neighborSearchSignallerBuilder_.build(
                inputrec->nstlist, inputrec->init_step, inputrec->init_t));
    }

    // Move setup / teardown list
    algorithm.elementSetupTeardownList_ = std::move(setupAndTeardownList_);
    // Move pre- / post-step scheduling lists
    algorithm.preStepScheduling_  = std::move(preStepScheduling_);
    algorithm.postStepScheduling_ = std::move(postStepScheduling_);

    // Create element list
    // Checkpoint helper needs to be in the call list (as first element!) to react to last step
    algorithm.elementCallList_.emplace_back(algorithm.checkpointHelper_.get());
    // Next, update the free energy lambda vector if needed
    if (algorithm.freeEnergyPerturbationData_)
    {
        algorithm.elementCallList_.emplace_back(algorithm.freeEnergyPerturbationData_->element());
    }
    // Then, move the built algorithm
    algorithm.elementsOwnershipList_.insert(algorithm.elementsOwnershipList_.end(),
                                            std::make_move_iterator(elements_.begin()),
                                            std::make_move_iterator(elements_.end()));
    algorithm.elementCallList_.insert(algorithm.elementCallList_.end(),
                                      std::make_move_iterator(callList_.begin()),
                                      std::make_move_iterator(callList_.end()));
    // Finally, all trajectory writing is happening after the step
    // (relevant data was stored by elements through energy signaller)
    algorithm.elementsOwnershipList_.emplace_back(std::move(trajectoryElement));
    algorithm.elementCallList_.emplace_back(algorithm.elementsOwnershipList_.back().get());
    algorithm.elementSetupTeardownList_.emplace_back(algorithm.elementsOwnershipList_.back().get());

    algorithm.setup();
    return algorithm;
}

bool ModularSimulatorAlgorithmBuilder::elementExists(const ISimulatorElement* element) const
{
    // Check whether element exists in element list
    if (std::any_of(elements_.begin(), elements_.end(), [element](auto& existingElement) {
            return element == existingElement.get();
        }))
    {
        return true;
    }
    // Check whether element exists in other places controlled by *this
    return ((statePropagatorData_ && statePropagatorData_->element() == element)
            || (energyData_ && energyData_->element() == element)
            || (freeEnergyPerturbationData_ && freeEnergyPerturbationData_->element() == element));
}

std::optional<SignallerCallback> ModularSimulatorAlgorithm::SignalHelper::registerLastStepCallback()
{
    return [this](Step step, Time gmx_unused time) { this->lastStep_ = step; };
}

std::optional<SignallerCallback> ModularSimulatorAlgorithm::SignalHelper::registerNSCallback()
{
    return [this](Step step, Time gmx_unused time) { this->nextNSStep_ = step; };
}

GlobalCommunicationHelper::GlobalCommunicationHelper(int nstglobalcomm, SimulationSignals* simulationSignals) :
    nstglobalcomm_(nstglobalcomm), simulationSignals_(simulationSignals)
{
}

int GlobalCommunicationHelper::nstglobalcomm() const
{
    return nstglobalcomm_;
}

SimulationSignals* GlobalCommunicationHelper::simulationSignals()
{
    return simulationSignals_;
}

ModularSimulatorAlgorithmBuilderHelper::ModularSimulatorAlgorithmBuilderHelper(
        ModularSimulatorAlgorithmBuilder* builder) :
    builder_(builder)
{
}

bool ModularSimulatorAlgorithmBuilderHelper::elementIsStored(const ISimulatorElement* element) const
{
    return builder_->elementExists(element);
}

[[maybe_unused]] void ModularSimulatorAlgorithmBuilderHelper::registerPreStepScheduling(SchedulingFunction schedulingFunction)
{
    builder_->preStepScheduling_.emplace_back(std::move(schedulingFunction));
}

[[maybe_unused]] void ModularSimulatorAlgorithmBuilderHelper::registerPostStepScheduling(SchedulingFunction schedulingFunction)
{
    builder_->postStepScheduling_.emplace_back(std::move(schedulingFunction));
}

std::optional<std::any> ModularSimulatorAlgorithmBuilderHelper::builderData(const std::string& key) const
{
    const auto iter = builder_->builderData_.find(key);
    if (iter == builder_->builderData_.end())
    {
        return std::nullopt;
    }
    else
    {
        return iter->second;
    }
}

void ModularSimulatorAlgorithmBuilderHelper::registerTemperaturePressureControl(
        std::function<void(const PropagatorConnection&)> registrationFunction)
{
    builder_->pressureTemperatureControlRegistrationFunctions_.emplace_back(std::move(registrationFunction));
}

void ModularSimulatorAlgorithmBuilderHelper::registerPropagator(PropagatorConnection connectionData)
{
    builder_->propagatorConnections_.emplace_back(std::move(connectionData));
}

void ModularSimulatorAlgorithmBuilderHelper::registerReferenceTemperatureUpdate(
        ReferenceTemperatureCallback referenceTemperatureCallback)
{
    auto* referenceTemperatureManager =
            simulationData<ReferenceTemperatureManager>("ReferenceTemperatureManager").value();
    referenceTemperatureManager->registerUpdateCallback(std::move(referenceTemperatureCallback));
}

ReferenceTemperatureCallback ModularSimulatorAlgorithmBuilderHelper::changeReferenceTemperatureCallback()
{
    // Capture is safe because SimulatorAlgorithm will manage life time of both the
    // recipient of the callback and the reference temperature manager
    auto* referenceTemperatureManager =
            simulationData<ReferenceTemperatureManager>("ReferenceTemperatureManager").value();
    return [referenceTemperatureManager](ArrayRef<const real>                temperatures,
                                         ReferenceTemperatureChangeAlgorithm algorithm) {
        referenceTemperatureManager->setReferenceTemperature(temperatures, algorithm);
    };
}

} // namespace gmx
