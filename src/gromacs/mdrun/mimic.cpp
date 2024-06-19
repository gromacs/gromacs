/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 *
 * \brief Declares the loop for MiMiC QM/MM
 *
 * \author Viacheslav Bolnykh <v.bolnykh@hpc-leap.eu>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <vector>

#include "gromacs/applied_forces/awh/awh.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/pointers.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/localtopologychecker.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/ebin.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/expanded.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/freeenergyparameters.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/membed.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mimic/communicator.h"
#include "gromacs/mimic/utilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/taskassignment/include/gromacs/taskassignment/decidesimulationworkload.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

#include "legacysimulator.h"
#include "replicaexchange.h"
#include "shellfc.h"

struct gmx_edsam;
struct gmx_mdoutf;
struct gmx_shellfc_t;

using gmx::SimulationSignaller;

void gmx::LegacySimulator::do_mimic()
{
    const t_inputrec* ir = inputRec_;
    double            t;
    bool              isLastStep               = false;
    bool              doFreeEnergyPerturbation = false;
    unsigned int      force_flags;
    tensor            force_vir, shake_vir, total_vir, pres;
    rvec              mu_tot;
    ForceBuffers      f;
    gmx_global_stat_t gstat;
    gmx_shellfc_t*    shellfc;

    double cycles;

    SimulationSignals signals;
    // Most global communnication stages don't propagate mdrun
    // signals, and will use this object to achieve that.
    SimulationSignaller nullSignaller(nullptr, nullptr, nullptr, false, false);

    if (ir->bExpanded)
    {
        gmx_fatal(FARGS, "Expanded ensemble not supported by MiMiC.");
    }
    if (ir->bSimTemp)
    {
        gmx_fatal(FARGS, "Simulated tempering not supported by MiMiC.");
    }
    if (ir->bDoAwh)
    {
        gmx_fatal(FARGS, "AWH not supported by MiMiC.");
    }
    if (replExParams_.exchangeInterval > 0)
    {
        gmx_fatal(FARGS, "Replica exchange not supported by MiMiC.");
    }
    if (opt2bSet("-ei", nFile_, fnm_) || observablesHistory_->edsamHistory != nullptr)
    {
        gmx_fatal(FARGS, "Essential dynamics not supported by MiMiC.");
    }
    if (ir->bIMD)
    {
        gmx_fatal(FARGS, "Interactive MD not supported by MiMiC.");
    }
    if (isMultiSim(ms_))
    {
        gmx_fatal(FARGS, "Multiple simulations not supported by MiMiC.");
    }
    if (std::any_of(ir->opts.annealing, ir->opts.annealing + ir->opts.ngtc, [](SimulatedAnnealing i) {
            return i != SimulatedAnnealing::No;
        }))
    {
        gmx_fatal(FARGS, "Simulated annealing not supported by MiMiC.");
    }

    /* Settings for rerun */
    {
        // TODO: Avoid changing inputrec (#3854)
        auto* nonConstInputrec               = const_cast<t_inputrec*>(inputRec_);
        nonConstInputrec->nstlist            = 1;
        nonConstInputrec->nstcalcenergy      = 1;
        nonConstInputrec->nstxout_compressed = 0;
    }
    int        nstglobalcomm = 1;
    const bool bNS           = true;

    ObservablesReducer observablesReducer = observablesReducerBuilder_->build();

    if (MAIN(cr_))
    {
        MimicCommunicator::init();
        auto* nonConstGlobalTopology = const_cast<gmx_mtop_t*>(&topGlobal_);
        MimicCommunicator::sendInitData(nonConstGlobalTopology, stateGlobal_->x);
        // TODO: Avoid changing inputrec (#3854)
        auto* nonConstInputrec   = const_cast<t_inputrec*>(inputRec_);
        nonConstInputrec->nsteps = MimicCommunicator::getStepNumber();
    }
    if (haveDDAtomOrdering(*cr_))
    {
        // TODO: Avoid changing inputrec (#3854)
        auto* nonConstInputrec = const_cast<t_inputrec*>(inputRec_);
        gmx_bcast(sizeof(ir->nsteps), &nonConstInputrec->nsteps, cr_->mpi_comm_mygroup);
    }

    const SimulationGroups* groups = &topGlobal_.groups;
    {
        auto* nonConstGlobalTopology                         = const_cast<gmx_mtop_t*>(&topGlobal_);
        nonConstGlobalTopology->intermolecularExclusionGroup = genQmmmIndices(topGlobal_);
    }

    initialize_lambdas(fpLog_,
                       ir->efep,
                       ir->bSimTemp,
                       *ir->fepvals,
                       ir->simtempvals->temperatures,
                       ekind_,
                       MAIN(cr_),
                       &stateGlobal_->fep_state,
                       stateGlobal_->lambda);

    const bool        simulationsShareState = false;
    gmx_mdoutf*       outf                  = init_mdoutf(fpLog_,
                                   nFile_,
                                   fnm_,
                                   mdrunOptions_,
                                   cr_,
                                   outputProvider_,
                                   mdModulesNotifiers_,
                                   ir,
                                   topGlobal_,
                                   oenv_,
                                   wallCycleCounters_,
                                   StartingBehavior::NewSimulation,
                                   simulationsShareState,
                                   ms_);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf),
                                   topGlobal_,
                                   *ir,
                                   pullWork_,
                                   mdoutf_get_fp_dhdl(outf),
                                   true,
                                   StartingBehavior::NewSimulation,
                                   simulationsShareState,
                                   mdModulesNotifiers_);

    gstat = global_stat_init(ir);

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fpLog_,
                                 topGlobal_,
                                 constr_ ? constr_->numFlexibleConstraints() : 0,
                                 ir->nstcalcenergy,
                                 haveDDAtomOrdering(*cr_),
                                 runScheduleWork_->simulationWork.useGpuPme);

    {
        double io = compute_io(ir, topGlobal_.natoms, *groups, energyOutput.numEnergyTerms(), 1);
        if ((io > 2000) && MAIN(cr_))
        {
            fprintf(stderr, "\nWARNING: This run will generate roughly %.0f Mb of data\n\n", io);
        }
    }

    if (haveDDAtomOrdering(*cr_))
    {
        // Local state only becomes valid now.
        dd_init_local_state(*cr_->dd, stateGlobal_, state_);

        /* Distribute the charge groups over the nodes from the main node */
        dd_partition_system(fpLog_,
                            mdLog_,
                            ir->init_step,
                            cr_,
                            TRUE,
                            stateGlobal_,
                            topGlobal_,
                            *ir,
                            mdModulesNotifiers_,
                            imdSession_,
                            pullWork_,
                            state_,
                            &f,
                            mdAtoms_,
                            top_,
                            fr_,
                            virtualSites_,
                            constr_,
                            nrnb_,
                            nullptr,
                            FALSE);
    }
    else
    {
        mdAlgorithmsSetupAtomData(
                cr_, *ir, topGlobal_, top_, fr_, &f, mdAtoms_, constr_, virtualSites_, shellfc);
    }

    auto* mdatoms = mdAtoms_->mdatoms();

    // NOTE: The global state is no longer used at this point.
    // But state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdatoms, state_->lambda[FreeEnergyPerturbationCouplingType::Mass]);
    fr_->longRangeNonbondeds->updateAfterPartition(*mdatoms);

    if (ir->efep != FreeEnergyPerturbationType::No && ir->fepvals->nstdhdl != 0)
    {
        doFreeEnergyPerturbation = true;
    }

    int64_t step     = ir->init_step;
    int64_t step_rel = 0;

    {
        int    cglo_flags   = CGLO_GSTAT;
        bool   bSumEkinhOld = false;
        t_vcm* vcm          = nullptr;
        compute_globals(gstat,
                        cr_,
                        ir,
                        fr_,
                        ekind_,
                        makeConstArrayRef(state_->x),
                        makeConstArrayRef(state_->v),
                        state_->box,
                        mdatoms,
                        nrnb_,
                        vcm,
                        nullptr,
                        enerd_,
                        force_vir,
                        shake_vir,
                        total_vir,
                        pres,
                        &nullSignaller,
                        state_->box,
                        &bSumEkinhOld,
                        cglo_flags,
                        step,
                        &observablesReducer);
        // Clean up after pre-step use of compute_globals()
        observablesReducer.markAsReadyToReduce();
    }

    if (MAIN(cr_))
    {
        fprintf(stderr, "starting MiMiC MD run '%s'\n\n", *(topGlobal_.name));
        if (mdrunOptions_.verbose)
        {
            fprintf(stderr,
                    "Calculated time to finish depends on nsteps from "
                    "run input file,\nwhich may not correspond to the time "
                    "needed to process input trajectory.\n\n");
        }
        fprintf(fpLog_, "\n");
    }

    walltime_accounting_start_time(wallTimeAccounting_);
    wallcycle_start(wallCycleCounters_, WallCycleCounter::Run);
    print_start(fpLog_, cr_, wallTimeAccounting_, "mdrun");

    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    if (constr_)
    {
        GMX_LOG(mdLog_.info)
                .asParagraph()
                .appendText(
                        "Simulations has constraints. Constraints will "
                        "be handled by CPMD.");
    }

    GMX_LOG(mdLog_.info)
            .asParagraph()
            .appendText(
                    "MiMiC does not report kinetic energy, total energy, temperature, virial and "
                    "pressure.");

    auto stopHandler = stopHandlerBuilder_->getStopHandlerMD(
            compat::not_null<SimulationSignal*>(&signals[eglsSTOPCOND]),
            false,
            MAIN(cr_),
            ir->nstlist,
            mdrunOptions_.reproducible,
            nstglobalcomm,
            mdrunOptions_.maximumHoursToRun,
            ir->nstlist == 0,
            fpLog_,
            step,
            bNS,
            wallTimeAccounting_);

    // we don't do counter resetting in rerun - finish will always be valid
    walltime_accounting_set_valid_finish(wallTimeAccounting_);

    const DDBalanceRegionHandler ddBalanceRegionHandler(cr_);

    /* and stop now if we should */
    isLastStep = (isLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!isLastStep)
    {
        isLastStep = (isLastStep || (ir->nsteps >= 0 && step_rel == ir->nsteps));
        wallcycle_start(wallCycleCounters_, WallCycleCounter::Step);

        t = step;

        if (MAIN(cr_))
        {
            MimicCommunicator::getCoords(stateGlobal_->x, stateGlobal_->numAtoms());
        }

        if (ir->efep != FreeEnergyPerturbationType::No)
        {
            state_->lambda = currentLambdas(step, *(ir->fepvals), stateGlobal_->fep_state);
        }

        if (MAIN(cr_))
        {
            const bool constructVsites =
                    ((virtualSites_ != nullptr) && mdrunOptions_.rerunConstructVsites);
            if (constructVsites && haveDDAtomOrdering(*cr_))
            {
                gmx_fatal(FARGS,
                          "Vsite recalculation with -rerun is not implemented with domain "
                          "decomposition, "
                          "use a single rank");
            }
            if (constructVsites)
            {
                wallcycle_start(wallCycleCounters_, WallCycleCounter::VsiteConstr);
                virtualSites_->construct(
                        state_->x, state_->v, state_->box, VSiteOperation::PositionsAndVelocities);
                wallcycle_stop(wallCycleCounters_, WallCycleCounter::VsiteConstr);
            }
        }

        if (haveDDAtomOrdering(*cr_))
        {
            /* Repartition the domain decomposition */
            const bool bMainState = true;
            dd_partition_system(fpLog_,
                                mdLog_,
                                step,
                                cr_,
                                bMainState,
                                stateGlobal_,
                                topGlobal_,
                                *ir,
                                mdModulesNotifiers_,
                                imdSession_,
                                pullWork_,
                                state_,
                                &f,
                                mdAtoms_,
                                top_,
                                fr_,
                                virtualSites_,
                                constr_,
                                nrnb_,
                                wallCycleCounters_,
                                mdrunOptions_.verbose);
        }

        if (MAIN(cr_))
        {
            EnergyOutput::printHeader(fpLog_, step, t); /* can we improve the information printed here? */
        }

        if (ir->efep != FreeEnergyPerturbationType::No)
        {
            update_mdatoms(mdatoms, state_->lambda[FreeEnergyPerturbationCouplingType::Mass]);
        }

        fr_->longRangeNonbondeds->updateAfterPartition(*mdatoms);

        force_flags = (GMX_FORCE_STATECHANGED | GMX_FORCE_DYNAMICBOX | GMX_FORCE_ALLFORCES
                       | GMX_FORCE_VIRIAL | // TODO: Get rid of this once #2649 is solved
                       GMX_FORCE_ENERGY | (doFreeEnergyPerturbation ? GMX_FORCE_DHDL : 0));

        const int shellfcFlags     = force_flags | (mdrunOptions_.verbose ? GMX_FORCE_ENERGY : 0);
        const int legacyForceFlags = ((shellfc) ? shellfcFlags : force_flags) | GMX_FORCE_NS;

        gmx_edsam* const ed = nullptr;

        if (bNS)
        {
            if (fr_->listedForcesGpu)
            {
                fr_->listedForcesGpu->updateHaveInteractions(top_->idef);
            }
            runScheduleWork_->domainWork = setupDomainLifetimeWorkload(
                    *ir, *fr_, pullWork_, ed, *mdatoms, runScheduleWork_->simulationWork);
        }


        runScheduleWork_->stepWork = setupStepWorkload(legacyForceFlags,
                                                       ir->mtsLevels,
                                                       step,
                                                       runScheduleWork_->domainWork,
                                                       runScheduleWork_->simulationWork);

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            relax_shell_flexcon(fpLog_,
                                cr_,
                                ms_,
                                mdrunOptions_.verbose,
                                enforcedRotation_,
                                step,
                                ir,
                                mdModulesNotifiers_,
                                imdSession_,
                                pullWork_,
                                bNS,
                                top_,
                                constr_,
                                enerd_,
                                state_->numAtoms(),
                                state_->x.arrayRefWithPadding(),
                                state_->v.arrayRefWithPadding(),
                                state_->box,
                                state_->lambda,
                                &state_->hist,
                                &f.view(),
                                force_vir,
                                *mdatoms,
                                fr_->longRangeNonbondeds.get(),
                                nrnb_,
                                wallCycleCounters_,
                                shellfc,
                                fr_,
                                *runScheduleWork_,
                                t,
                                mu_tot,
                                virtualSites_,
                                ddBalanceRegionHandler);
        }
        else
        {
            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            Awh* awh = nullptr;

            do_force(fpLog_,
                     cr_,
                     ms_,
                     *ir,
                     mdModulesNotifiers_,
                     awh,
                     enforcedRotation_,
                     imdSession_,
                     pullWork_,
                     step,
                     nrnb_,
                     wallCycleCounters_,
                     top_,
                     state_->box,
                     state_->x.arrayRefWithPadding(),
                     state_->v.arrayRefWithPadding().unpaddedArrayRef(),
                     &state_->hist,
                     &f.view(),
                     force_vir,
                     mdatoms,
                     enerd_,
                     state_->lambda,
                     fr_,
                     *runScheduleWork_,
                     virtualSites_,
                     mu_tot,
                     t,
                     ed,
                     fr_->longRangeNonbondeds.get(),
                     ddBalanceRegionHandler);
        }

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t.
         */
        {
            const bool isCheckpointingStep = false;
            const bool doRerun             = false;
            const bool bSumEkinhOld        = false;
            do_md_trajectory_writing(fpLog_,
                                     cr_,
                                     nFile_,
                                     fnm_,
                                     step,
                                     step_rel,
                                     t,
                                     ir,
                                     state_,
                                     stateGlobal_,
                                     observablesHistory_,
                                     topGlobal_,
                                     fr_,
                                     outf,
                                     energyOutput,
                                     ekind_,
                                     f.view().force(),
                                     isCheckpointingStep,
                                     doRerun,
                                     isLastStep,
                                     mdrunOptions_.writeConfout,
                                     bSumEkinhOld ? EkindataState::UsedNeedToReduce
                                                  : EkindataState::UsedDoNotNeedToReduce);
        }

        stopHandler->setSignal();

        {
            const bool          doInterSimSignal = false;
            const bool          doIntraSimSignal = true;
            bool                bSumEkinhOld     = false;
            t_vcm*              vcm              = nullptr;
            SimulationSignaller signaller(&signals, cr_, ms_, doInterSimSignal, doIntraSimSignal);

            int cglo_flags = CGLO_GSTAT | CGLO_ENERGY;
            compute_globals(gstat,
                            cr_,
                            ir,
                            fr_,
                            ekind_,
                            makeConstArrayRef(state_->x),
                            makeConstArrayRef(state_->v),
                            state_->box,
                            mdatoms,
                            nrnb_,
                            vcm,
                            wallCycleCounters_,
                            enerd_,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            &signaller,
                            state_->box,
                            &bSumEkinhOld,
                            cglo_flags,
                            step,
                            &observablesReducer);
        }

        {
            gmx::HostVector<gmx::RVec>     fglobal(topGlobal_.natoms);
            gmx::ArrayRef<gmx::RVec>       ftemp;
            gmx::ArrayRef<const gmx::RVec> flocal = f.view().force();
            if (haveDDAtomOrdering(*cr_))
            {
                ftemp = gmx::makeArrayRef(fglobal);
                dd_collect_vec(
                        cr_->dd, state_->ddp_count, state_->ddp_count_cg_gl, state_->cg_gl, flocal, ftemp);
            }
            else
            {
                ftemp = f.view().force();
            }

            if (MAIN(cr_))
            {
                MimicCommunicator::sendEnergies(enerd_->term[F_EPOT]);
                MimicCommunicator::sendForces(ftemp, stateGlobal_->numAtoms());
            }
        }


        /* Note: this is OK, but there are some numerical precision issues with using the convergence of
           the virial that should probably be addressed eventually. state->veta has better properies,
           but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
           generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

        if (ir->efep != FreeEnergyPerturbationType::No)
        {
            /* Sum up the foreign energy and dhdl terms for md and sd.
               Currently done every step so that dhdl is correct in the .edr */
            accumulateKineticLambdaComponents(enerd_, state_->lambda, *ir->fepvals);
        }

        /* Output stuff */
        if (MAIN(cr_))
        {
            const bool bCalcEnerStep = true;
            energyOutput.addDataAtEnergyStep(doFreeEnergyPerturbation,
                                             bCalcEnerStep,
                                             t,
                                             mdatoms->tmass,
                                             enerd_,
                                             ir->fepvals.get(),
                                             state_->box,
                                             PTCouplingArrays({ state_->boxv,
                                                                state_->nosehoover_xi,
                                                                state_->nosehoover_vxi,
                                                                state_->nhpres_xi,
                                                                state_->nhpres_vxi }),
                                             state_->fep_state,
                                             total_vir,
                                             pres,
                                             ekind_,
                                             mu_tot,
                                             constr_);

            const bool do_ene = true;
            const bool do_log = true;
            Awh*       awh    = nullptr;
            const bool do_dr  = ir->nstdisreout != 0;
            const bool do_or  = ir->nstorireout != 0;

            EnergyOutput::printAnnealingTemperatures(do_log ? fpLog_ : nullptr, *groups, ir->opts, *ekind_);
            energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                               do_ene,
                                               do_dr,
                                               do_or,
                                               do_log ? fpLog_ : nullptr,
                                               step,
                                               t,
                                               fr_->fcdata.get(),
                                               awh);

            if (do_per_step(step, ir->nstlog))
            {
                if (fflush(fpLog_) != 0)
                {
                    gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }

        /* Print the remaining wall clock time for the run */
        if (isMainSimMainRank(ms_, MAIN(cr_)) && (mdrunOptions_.verbose || gmx_got_usr_signal()))
        {
            if (shellfc)
            {
                fprintf(stderr, "\n");
            }
            print_time(stderr, wallTimeAccounting_, step, ir, cr_);
        }

        cycles = wallcycle_stop(wallCycleCounters_, WallCycleCounter::Step);
        if (haveDDAtomOrdering(*cr_) && wallCycleCounters_)
        {
            dd_cycles_add(cr_->dd, cycles, ddCyclStep);
        }

        /* increase the MD step number */
        step++;
        step_rel++;
        observablesReducer.markAsReadyToReduce();
    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end_time(wallTimeAccounting_);

    if (MAIN(cr_))
    {
        MimicCommunicator::finalize();
    }

    if (!thisRankHasDuty(cr_, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr_);
    }

    done_mdoutf(outf);

    done_shellfc(fpLog_, shellfc, step_rel);

    walltime_accounting_set_nsteps_done(wallTimeAccounting_, step_rel);
}
