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
 * \brief Defines the modular simulator algorithm
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "simulatoralgorithm.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdrun/replicaexchange.h"
#include "gromacs/mdrun/shellfc.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"

#include "domdechelper.h"
#include "freeenergyperturbationelement.h"
#include "modularsimulator.h"
#include "parrinellorahmanbarostat.h"
#include "signallers.h"
#include "trajectoryelement.h"
#include "vrescalethermostat.h"

namespace gmx
{
ModularSimulatorAlgorithm::ModularSimulatorAlgorithm(std::string              topologyName,
                                                     FILE*                    fplog,
                                                     t_commrec*               cr,
                                                     const MDLogger&          mdlog,
                                                     const MdrunOptions&      mdrunOptions,
                                                     t_inputrec*              inputrec,
                                                     t_nrnb*                  nrnb,
                                                     gmx_wallcycle*           wcycle,
                                                     t_forcerec*              fr,
                                                     gmx_walltime_accounting* walltime_accounting) :
    taskIterator_(taskQueue_.end()),
    step_(-1),
    runFinished_(false),
    topologyName_(std::move(topologyName)),
    fplog(fplog),
    cr(cr),
    mdlog(mdlog),
    mdrunOptions(mdrunOptions),
    inputrec(inputrec),
    nrnb(nrnb),
    wcycle(wcycle),
    fr(fr),
    walltime_accounting(walltime_accounting)
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

    for (auto& element : elementsOwnershipList_)
    {
        element->elementSetup();
    }
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
        taskIterator_++;
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
    return taskIterator_->get();
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
    for (auto& element : elementsOwnershipList_)
    {
        element->elementTeardown();
    }
    if (pmeLoadBalanceHelper_)
    {
        pmeLoadBalanceHelper_->teardown();
    }
    simulatorTeardown();
}

void ModularSimulatorAlgorithm::simulatorSetup()
{
    if (!mdrunOptions.writeConfout)
    {
        // This is on by default, and the main known use case for
        // turning it off is for convenience in benchmarking, which is
        // something that should not show up in the general user
        // interface.
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendText(
                        "The -noconfout functionality is deprecated, and "
                        "may be removed in a future version.");
    }

    if (MASTER(cr))
    {
        char        sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
        std::string timeString;
        fprintf(stderr, "starting mdrun '%s'\n", topologyName_.c_str());
        if (inputrec->nsteps >= 0)
        {
            timeString = formatString("%8.1f", static_cast<double>(inputrec->init_step + inputrec->nsteps)
                                                       * inputrec->delta_t);
        }
        else
        {
            timeString = "infinite";
        }
        if (inputrec->init_step > 0)
        {
            fprintf(stderr, "%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                    gmx_step_str(inputrec->init_step + inputrec->nsteps, sbuf), timeString.c_str(),
                    gmx_step_str(inputrec->init_step, sbuf2), inputrec->init_step * inputrec->delta_t);
        }
        else
        {
            fprintf(stderr, "%s steps, %s ps.\n", gmx_step_str(inputrec->nsteps, sbuf),
                    timeString.c_str());
        }
        fprintf(fplog, "\n");
    }

    walltime_accounting_start_time(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, "mdrun");

    step_ = inputrec->init_step;
}

void ModularSimulatorAlgorithm::simulatorTeardown()
{

    // Stop measuring walltime
    walltime_accounting_end_time(walltime_accounting);

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    walltime_accounting_set_nsteps_done(walltime_accounting, step_ - inputrec->init_step);
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

    resetHandler_->setSignal(walltime_accounting);
    // This is a hack to avoid having to rewrite StopHandler to be a NeighborSearchSignaller
    // and accept the step as input. Eventually, we want to do that, but currently this would
    // require introducing NeighborSearchSignaller in the legacy do_md or a lot of code
    // duplication.
    stophandlerIsNSStep_    = isNeighborSearchingStep;
    stophandlerCurrentStep_ = step;
    stopHandler_->setSignal();

    wallcycle_start(wcycle, ewcSTEP);
}

void ModularSimulatorAlgorithm::postStep(Step step, Time gmx_unused time)
{
    // Output stuff
    if (MASTER(cr))
    {
        if (do_per_step(step, inputrec->nstlog))
        {
            if (fflush(fplog) != 0)
            {
                gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
            }
        }
    }
    const bool do_verbose = mdrunOptions.verbose
                            && (step % mdrunOptions.verboseStepPrintInterval == 0
                                || step == inputrec->init_step || step == signalHelper_->lastStep_);
    // Print the remaining wall clock time for the run
    if (MASTER(cr) && (do_verbose || gmx_got_usr_signal())
        && !(pmeLoadBalanceHelper_ && pmeLoadBalanceHelper_->pmePrinting()))
    {
        print_time(stderr, walltime_accounting, step, inputrec, cr);
    }

    double cycles = wallcycle_stop(wcycle, ewcSTEP);
    if (DOMAINDECOMP(cr) && wcycle)
    {
        dd_cycles_add(cr->dd, static_cast<float>(cycles), ddCyclStep);
    }

    resetHandler_->resetCounters(
            step, step - inputrec->init_step, mdlog, fplog, cr, fr->nbv.get(), nrnb, fr->pmedata,
            pmeLoadBalanceHelper_ ? pmeLoadBalanceHelper_->loadBalancingObject() : nullptr, wcycle,
            walltime_accounting);
}

void ModularSimulatorAlgorithm::populateTaskQueue()
{
    /*
     * The registerRunFunction emplaces functions to the task queue.
     * All elements are owned by the ModularSimulatorAlgorithm, as is the task queue.
     * Elements can hence register lambdas capturing their `this` pointers without expecting
     * life time issues, as the task queue and the elements are in the same scope.
     */
    auto registerRunFunction = std::make_unique<RegisterRunFunction>(
            [this](SimulatorRunFunctionPtr ptr) { taskQueue_.emplace_back(std::move(ptr)); });

    Time startTime = inputrec->init_t;
    Time timeStep  = inputrec->delta_t;
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
        (*registerRunFunction)(std::make_unique<SimulatorRunFunction>(
                [this, step, time, isNSStep]() { preStep(step, time, isNSStep); }));
        // register elements for step
        for (auto& element : elementCallList_)
        {
            element->scheduleTask(step_, time, registerRunFunction);
        }
        // register post-step (task queue is local, so no problem with `this`)
        (*registerRunFunction)(
                std::make_unique<SimulatorRunFunction>([this, step, time]() { postStep(step, time); }));

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
        (*registerRunFunction)(std::make_unique<SimulatorRunFunction>([this]() { teardown(); }));
    }
}

ModularSimulatorAlgorithm ModularSimulatorAlgorithmBuilder::constructElementsAndSignallers()
{
    ModularSimulatorAlgorithm algorithm(*(top_global->name), fplog, cr, mdlog, mdrunOptions,
                                        inputrec, nrnb, wcycle, fr, walltime_accounting);
    /* When restarting from a checkpoint, it can be appropriate to
     * initialize ekind from quantities in the checkpoint. Otherwise,
     * compute_globals must initialize ekind before the simulation
     * starts/restarts. However, only the master rank knows what was
     * found in the checkpoint file, so we have to communicate in
     * order to coordinate the restart.
     *
     * TODO (modular) This should become obsolete when checkpoint reading
     *      happens within the modular simulator framework: The energy
     *      element should read its data from the checkpoint file pointer,
     *      and signal to the compute globals element if it needs anything
     *      reduced.
     *
     * TODO (legacy) Consider removing this communication if/when checkpoint
     *      reading directly follows .tpr reading, because all ranks can
     *      agree on hasReadEkinState at that time.
     */
    bool hasReadEkinState = MASTER(cr) ? state_global->ekinstate.hasReadEkinState : false;
    if (PAR(cr))
    {
        gmx_bcast(sizeof(hasReadEkinState), &hasReadEkinState, cr->mpi_comm_mygroup);
    }
    if (hasReadEkinState)
    {
        restore_ekinstate_from_state(cr, ekind, &state_global->ekinstate);
    }

    /*
     * Build data structures
     */
    algorithm.topologyHolder_ =
            std::make_unique<TopologyHolder>(*top_global, cr, inputrec, fr, mdAtoms, constr, vsite);

    std::unique_ptr<FreeEnergyPerturbationElement> freeEnergyPerturbationElement    = nullptr;
    FreeEnergyPerturbationElement*                 freeEnergyPerturbationElementPtr = nullptr;
    if (inputrec->efep != efepNO)
    {
        freeEnergyPerturbationElement =
                std::make_unique<FreeEnergyPerturbationElement>(fplog, inputrec, mdAtoms);
        freeEnergyPerturbationElementPtr = freeEnergyPerturbationElement.get();
    }

    auto statePropagatorData = std::make_unique<StatePropagatorData>(
            top_global->natoms, fplog, cr, state_global, inputrec->nstxout, inputrec->nstvout,
            inputrec->nstfout, inputrec->nstxout_compressed, fr->nbv->useGpu(),
            freeEnergyPerturbationElementPtr, algorithm.topologyHolder_.get(), fr->bMolPBC,
            mdrunOptions.writeConfout, opt2fn("-c", nfile, fnm), inputrec, mdAtoms->mdatoms());
    auto statePropagatorDataPtr = compat::make_not_null(statePropagatorData.get());

    auto energyElement = std::make_unique<EnergyElement>(
            statePropagatorDataPtr, freeEnergyPerturbationElementPtr, top_global, inputrec, mdAtoms,
            enerd, ekind, constr, fplog, &fr->listedForces->fcdata(), mdModulesNotifier, MASTER(cr),
            observablesHistory, startingBehavior);
    auto energyElementPtr = compat::make_not_null(energyElement.get());

    /*
     * Build stop handler
     */
    const bool simulationsShareState = false;
    algorithm.stopHandler_           = stopHandlerBuilder->getStopHandlerMD(
            compat::not_null<SimulationSignal*>(&algorithm.signals_[eglsSTOPCOND]),
            simulationsShareState, MASTER(cr), inputrec->nstlist, mdrunOptions.reproducible,
            nstglobalcomm_, mdrunOptions.maximumHoursToRun, inputrec->nstlist == 0, fplog,
            algorithm.stophandlerCurrentStep_, algorithm.stophandlerIsNSStep_, walltime_accounting);

    /*
     * Create simulator builders
     */
    SignallerBuilder<NeighborSearchSignaller> neighborSearchSignallerBuilder;
    SignallerBuilder<LastStepSignaller>       lastStepSignallerBuilder;
    SignallerBuilder<LoggingSignaller>        loggingSignallerBuilder;
    SignallerBuilder<EnergySignaller>         energySignallerBuilder;
    SignallerBuilder<TrajectorySignaller>     trajectorySignallerBuilder;
    TrajectoryElementBuilder                  trajectoryElementBuilder;

    /*
     * Register data structures to signallers
     */
    trajectoryElementBuilder.registerWriterClient(statePropagatorDataPtr);
    trajectorySignallerBuilder.registerSignallerClient(statePropagatorDataPtr);
    lastStepSignallerBuilder.registerSignallerClient(statePropagatorDataPtr);

    trajectoryElementBuilder.registerWriterClient(energyElementPtr);
    trajectorySignallerBuilder.registerSignallerClient(energyElementPtr);
    energySignallerBuilder.registerSignallerClient(energyElementPtr);

    // Register the simulator itself to the neighbor search / last step signaller
    neighborSearchSignallerBuilder.registerSignallerClient(
            compat::make_not_null(algorithm.signalHelper_.get()));
    lastStepSignallerBuilder.registerSignallerClient(compat::make_not_null(algorithm.signalHelper_.get()));

    /*
     * Build integrator - this takes care of force calculation, propagation,
     * constraining, and of the place the statePropagatorData and the energy element
     * have a full timestep state.
     */
    // TODO: Make a CheckpointHelperBuilder
    std::vector<ICheckpointHelperClient*> checkpointClients = { statePropagatorDataPtr, energyElementPtr,
                                                                freeEnergyPerturbationElementPtr };
    CheckBondedInteractionsCallbackPtr checkBondedInteractionsCallback = nullptr;
    auto                               integrator                      = buildIntegrator(
            &neighborSearchSignallerBuilder, &energySignallerBuilder, &loggingSignallerBuilder,
            &trajectorySignallerBuilder, &checkpointClients, &checkBondedInteractionsCallback,
            statePropagatorDataPtr, energyElementPtr, freeEnergyPerturbationElementPtr,
            hasReadEkinState, algorithm.topologyHolder_.get(), &algorithm.signals_);

    /*
     * Build infrastructure elements
     */

    if (PmeLoadBalanceHelper::doPmeLoadBalancing(mdrunOptions, inputrec, fr))
    {
        algorithm.pmeLoadBalanceHelper_ = std::make_unique<PmeLoadBalanceHelper>(
                mdrunOptions.verbose, statePropagatorDataPtr, fplog, cr, mdlog, inputrec, wcycle, fr);
        neighborSearchSignallerBuilder.registerSignallerClient(
                compat::make_not_null(algorithm.pmeLoadBalanceHelper_.get()));
    }

    if (DOMAINDECOMP(cr))
    {
        GMX_ASSERT(checkBondedInteractionsCallback,
                   "Domain decomposition needs a callback for check the number of bonded "
                   "interactions.");
        algorithm.domDecHelper_ = std::make_unique<DomDecHelper>(
                mdrunOptions.verbose, mdrunOptions.verboseStepPrintInterval, statePropagatorDataPtr,
                algorithm.topologyHolder_.get(), std::move(checkBondedInteractionsCallback),
                nstglobalcomm_, fplog, cr, mdlog, constr, inputrec, mdAtoms, nrnb, wcycle, fr,
                vsite, imdSession, pull_work);
        neighborSearchSignallerBuilder.registerSignallerClient(
                compat::make_not_null(algorithm.domDecHelper_.get()));
    }

    const bool simulationsShareResetCounters = false;
    algorithm.resetHandler_                  = std::make_unique<ResetHandler>(
            compat::make_not_null<SimulationSignal*>(&algorithm.signals_[eglsRESETCOUNTERS]),
            simulationsShareResetCounters, inputrec->nsteps, MASTER(cr),
            mdrunOptions.timingOptions.resetHalfway, mdrunOptions.maximumHoursToRun, mdlog, wcycle,
            walltime_accounting);

    /*
     * Build signaller list
     *
     * Note that as signallers depend on each others, the order of calling the signallers
     * matters. It is the responsibility of this builder to ensure that the order is
     * maintained.
     */
    auto energySignaller = energySignallerBuilder.build(
            inputrec->nstcalcenergy, inputrec->fepvals->nstdhdl, inputrec->nstpcouple);
    trajectorySignallerBuilder.registerSignallerClient(compat::make_not_null(energySignaller.get()));
    loggingSignallerBuilder.registerSignallerClient(compat::make_not_null(energySignaller.get()));
    auto trajectoryElement = trajectoryElementBuilder.build(
            fplog, nfile, fnm, mdrunOptions, cr, outputProvider, mdModulesNotifier, inputrec,
            top_global, oenv, wcycle, startingBehavior, simulationsShareState);
    loggingSignallerBuilder.registerSignallerClient(compat::make_not_null(trajectoryElement.get()));
    trajectorySignallerBuilder.registerSignallerClient(compat::make_not_null(trajectoryElement.get()));
    auto trajectorySignaller = trajectorySignallerBuilder.build(
            inputrec->nstxout, inputrec->nstvout, inputrec->nstfout, inputrec->nstxout_compressed,
            trajectoryElement->tngBoxOut(), trajectoryElement->tngLambdaOut(),
            trajectoryElement->tngBoxOutCompressed(), trajectoryElement->tngLambdaOutCompressed(),
            inputrec->nstenergy);

    // Add checkpoint helper here since we need a pointer to the trajectory element and
    // need to register it with the lastStepSignallerBuilder
    auto checkpointHandler = std::make_unique<CheckpointHandler>(
            compat::make_not_null<SimulationSignal*>(&algorithm.signals_[eglsCHKPT]),
            simulationsShareState, inputrec->nstlist == 0, MASTER(cr), mdrunOptions.writeConfout,
            mdrunOptions.checkpointOptions.period);
    algorithm.checkpointHelper_ = std::make_unique<CheckpointHelper>(
            std::move(checkpointClients), std::move(checkpointHandler), inputrec->init_step,
            trajectoryElement.get(), top_global->natoms, fplog, cr, observablesHistory,
            walltime_accounting, state_global, mdrunOptions.writeConfout);
    lastStepSignallerBuilder.registerSignallerClient(
            compat::make_not_null(algorithm.checkpointHelper_.get()));

    lastStepSignallerBuilder.registerSignallerClient(compat::make_not_null(trajectorySignaller.get()));
    auto loggingSignaller =
            loggingSignallerBuilder.build(inputrec->nstlog, inputrec->init_step, inputrec->init_t);
    lastStepSignallerBuilder.registerSignallerClient(compat::make_not_null(loggingSignaller.get()));
    auto lastStepSignaller = lastStepSignallerBuilder.build(inputrec->nsteps, inputrec->init_step,
                                                            algorithm.stopHandler_.get());
    neighborSearchSignallerBuilder.registerSignallerClient(compat::make_not_null(lastStepSignaller.get()));
    auto neighborSearchSignaller = neighborSearchSignallerBuilder.build(
            inputrec->nstlist, inputrec->init_step, inputrec->init_t);

    algorithm.signallerList_.emplace_back(std::move(neighborSearchSignaller));
    algorithm.signallerList_.emplace_back(std::move(lastStepSignaller));
    algorithm.signallerList_.emplace_back(std::move(loggingSignaller));
    algorithm.signallerList_.emplace_back(std::move(trajectorySignaller));
    algorithm.signallerList_.emplace_back(std::move(energySignaller));

    /*
     * Build the element list
     *
     * This is the actual sequence of (non-infrastructure) elements to be run.
     * For NVE, the trajectory element is used outside of the integrator
     * (composite) element, as well as the checkpoint helper. The checkpoint
     * helper should be on top of the loop, and is only part of the simulator
     * call list to be able to react to the last step being signalled.
     */
    addToCallList(algorithm.checkpointHelper_, algorithm.elementCallList_);
    if (freeEnergyPerturbationElement)
    {
        addToCallListAndMove(std::move(freeEnergyPerturbationElement), algorithm.elementCallList_,
                             algorithm.elementsOwnershipList_);
    }
    addToCallListAndMove(std::move(integrator), algorithm.elementCallList_, algorithm.elementsOwnershipList_);
    addToCallListAndMove(std::move(trajectoryElement), algorithm.elementCallList_,
                         algorithm.elementsOwnershipList_);
    // for vv, we need to setup statePropagatorData after the compute
    // globals so that we reset the right velocities
    // TODO: Avoid this by getting rid of the need of resetting velocities in vv
    algorithm.elementsOwnershipList_.emplace_back(std::move(statePropagatorData));
    algorithm.elementsOwnershipList_.emplace_back(std::move(energyElement));

    return algorithm;
}

ModularSimulatorAlgorithmBuilder::ModularSimulatorAlgorithmBuilder(compat::not_null<ModularSimulator*> simulator) :
    nstglobalcomm_(computeGlobalCommunicationPeriod(simulator->mdlog, simulator->inputrec, simulator->cr)),
    fplog(simulator->fplog),
    cr(simulator->cr),
    ms(simulator->ms),
    mdlog(simulator->mdlog),
    nfile(simulator->nfile),
    fnm(simulator->fnm),
    oenv(simulator->oenv),
    mdrunOptions(simulator->mdrunOptions),
    startingBehavior(simulator->startingBehavior),
    vsite(simulator->vsite),
    constr(simulator->constr),
    enforcedRotation(simulator->enforcedRotation),
    deform(simulator->deform),
    outputProvider(simulator->outputProvider),
    mdModulesNotifier(simulator->mdModulesNotifier),
    inputrec(simulator->inputrec),
    imdSession(simulator->imdSession),
    pull_work(simulator->pull_work),
    swap(simulator->swap),
    top_global(simulator->top_global),
    state_global(simulator->state_global),
    observablesHistory(simulator->observablesHistory),
    mdAtoms(simulator->mdAtoms),
    nrnb(simulator->nrnb),
    wcycle(simulator->wcycle),
    fr(simulator->fr),
    enerd(simulator->enerd),
    ekind(simulator->ekind),
    runScheduleWork(simulator->runScheduleWork),
    replExParams(simulator->replExParams),
    membed(simulator->membed),
    walltime_accounting(simulator->walltime_accounting),
    stopHandlerBuilder(simulator->stopHandlerBuilder.get()),
    doRerun(simulator->doRerun)
{
}

ModularSimulatorAlgorithm ModularSimulatorAlgorithmBuilder::build()
{
    auto algorithm = constructElementsAndSignallers();
    algorithm.setup();
    return algorithm;
}

SignallerCallbackPtr ModularSimulatorAlgorithm::SignalHelper::registerLastStepCallback()
{
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->lastStep_ = step; });
}

SignallerCallbackPtr ModularSimulatorAlgorithm::SignalHelper::registerNSCallback()
{
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->nextNSStep_ = step; });
}
} // namespace gmx
