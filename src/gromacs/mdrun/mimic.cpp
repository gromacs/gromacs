/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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
#include <memory>

#include "gromacs/applied_forces/awh/awh.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/listed_forces/listed_forces.h"
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
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mimic/communicator.h"
#include "gromacs/mimic/utilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

#include "legacysimulator.h"
#include "replicaexchange.h"
#include "shellfc.h"

using gmx::SimulationSignaller;

void gmx::LegacySimulator::do_mimic()
{
    t_inputrec*       ir = inputrec;
    int64_t           step, step_rel;
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

    /* Domain decomposition could incorrectly miss a bonded
       interaction, but checking for that requires a global
       communication stage, which does not otherwise happen in DD
       code. So we do that alongside the first global energy reduction
       after a new DD is made. These variables handle whether the
       check happens, and the result it returns. */
    bool shouldCheckNumberOfBondedInteractions = false;
    int  totalNumberOfBondedInteractions       = -1;

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
    if (replExParams.exchangeInterval > 0)
    {
        gmx_fatal(FARGS, "Replica exchange not supported by MiMiC.");
    }
    if (opt2bSet("-ei", nfile, fnm) || observablesHistory->edsamHistory != nullptr)
    {
        gmx_fatal(FARGS, "Essential dynamics not supported by MiMiC.");
    }
    if (ir->bIMD)
    {
        gmx_fatal(FARGS, "Interactive MD not supported by MiMiC.");
    }
    if (isMultiSim(ms))
    {
        gmx_fatal(FARGS, "Multiple simulations not supported by MiMiC.");
    }
    if (std::any_of(ir->opts.annealing, ir->opts.annealing + ir->opts.ngtc,
                    [](int i) { return i != eannNO; }))
    {
        gmx_fatal(FARGS, "Simulated annealing not supported by MiMiC.");
    }

    /* Settings for rerun */
    ir->nstlist              = 1;
    ir->nstcalcenergy        = 1;
    int        nstglobalcomm = 1;
    const bool bNS           = true;

    if (MASTER(cr))
    {
        MimicCommunicator::init();
        auto nonConstGlobalTopology = const_cast<gmx_mtop_t*>(top_global);
        MimicCommunicator::sendInitData(nonConstGlobalTopology, state_global->x);
        ir->nsteps = MimicCommunicator::getStepNumber();
    }

    ir->nstxout_compressed         = 0;
    const SimulationGroups* groups = &top_global->groups;
    {
        auto nonConstGlobalTopology                          = const_cast<gmx_mtop_t*>(top_global);
        nonConstGlobalTopology->intermolecularExclusionGroup = genQmmmIndices(*top_global);
    }

    initialize_lambdas(fplog, *ir, MASTER(cr), &state_global->fep_state, state_global->lambda);

    const bool        simulationsShareState = false;
    gmx_mdoutf*       outf = init_mdoutf(fplog, nfile, fnm, mdrunOptions, cr, outputProvider,
                                   mdModulesNotifier, ir, top_global, oenv, wcycle,
                                   StartingBehavior::NewSimulation, simulationsShareState, ms);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf), top_global, ir, pull_work,
                                   mdoutf_get_fp_dhdl(outf), true, StartingBehavior::NewSimulation,
                                   simulationsShareState, mdModulesNotifier);

    gstat = global_stat_init(ir);

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog, top_global, constr ? constr->numFlexibleConstraints() : 0,
                                 ir->nstcalcenergy, DOMAINDECOMP(cr),
                                 runScheduleWork->simulationWork.useGpuPme);

    {
        double io = compute_io(ir, top_global->natoms, *groups, energyOutput.numEnergyTerms(), 1);
        if ((io > 2000) && MASTER(cr))
        {
            fprintf(stderr, "\nWARNING: This run will generate roughly %.0f Mb of data\n\n", io);
        }
    }

    // Local state only becomes valid now.
    std::unique_ptr<t_state> stateInstance;
    t_state*                 state;

    gmx_localtop_t top(top_global->ffparams);

    if (DOMAINDECOMP(cr))
    {
        stateInstance = std::make_unique<t_state>();
        state         = stateInstance.get();
        dd_init_local_state(cr->dd, state_global, state);

        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, mdlog, ir->init_step, cr, TRUE, 1, state_global, *top_global, ir,
                            imdSession, pull_work, state, &f, mdAtoms, &top, fr, vsite, constr,
                            nrnb, nullptr, FALSE);
        shouldCheckNumberOfBondedInteractions = true;
        gmx_bcast(sizeof(ir->nsteps), &ir->nsteps, cr->mpi_comm_mygroup);
    }
    else
    {
        state_change_natoms(state_global, state_global->natoms);
        /* Copy the pointer to the global state */
        state = state_global;

        mdAlgorithmsSetupAtomData(cr, ir, *top_global, &top, fr, &f, mdAtoms, constr, vsite, shellfc);
    }

    auto mdatoms = mdAtoms->mdatoms();

    // NOTE: The global state is no longer used at this point.
    // But state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdatoms, state->lambda[efptMASS]);

    if (ir->efep != efepNO && ir->fepvals->nstdhdl != 0)
    {
        doFreeEnergyPerturbation = true;
    }

    {
        int cglo_flags =
                (CGLO_GSTAT
                 | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0));
        bool   bSumEkinhOld = false;
        t_vcm* vcm          = nullptr;
        compute_globals(gstat, cr, ir, fr, ekind, makeConstArrayRef(state->x),
                        makeConstArrayRef(state->v), state->box, mdatoms, nrnb, vcm, nullptr, enerd,
                        force_vir, shake_vir, total_vir, pres, constr, &nullSignaller, state->box,
                        &totalNumberOfBondedInteractions, &bSumEkinhOld, cglo_flags);
    }
    checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions, top_global, &top,
                                    makeConstArrayRef(state->x), state->box,
                                    &shouldCheckNumberOfBondedInteractions);

    if (MASTER(cr))
    {
        fprintf(stderr, "starting MiMiC MD run '%s'\n\n", *(top_global->name));
        if (mdrunOptions.verbose)
        {
            fprintf(stderr,
                    "Calculated time to finish depends on nsteps from "
                    "run input file,\nwhich may not correspond to the time "
                    "needed to process input trajectory.\n\n");
        }
        fprintf(fplog, "\n");
    }

    walltime_accounting_start_time(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, "mdrun");

    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    if (constr)
    {
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendText(
                        "Simulations has constraints. Constraints will "
                        "be handled by CPMD.");
    }

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(
                    "MiMiC does not report kinetic energy, total energy, temperature, virial and "
                    "pressure.");

    step     = ir->init_step;
    step_rel = 0;

    auto stopHandler = stopHandlerBuilder->getStopHandlerMD(
            compat::not_null<SimulationSignal*>(&signals[eglsSTOPCOND]), false, MASTER(cr),
            ir->nstlist, mdrunOptions.reproducible, nstglobalcomm, mdrunOptions.maximumHoursToRun,
            ir->nstlist == 0, fplog, step, bNS, walltime_accounting);

    // we don't do counter resetting in rerun - finish will always be valid
    walltime_accounting_set_valid_finish(walltime_accounting);

    const DDBalanceRegionHandler ddBalanceRegionHandler(cr);

    /* and stop now if we should */
    isLastStep = (isLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!isLastStep)
    {
        isLastStep = (isLastStep || (ir->nsteps >= 0 && step_rel == ir->nsteps));
        wallcycle_start(wcycle, ewcSTEP);

        t = step;

        if (MASTER(cr))
        {
            MimicCommunicator::getCoords(&state_global->x, state_global->natoms);
        }

        if (ir->efep != efepNO)
        {
            state->lambda = currentLambdas(step, *(ir->fepvals), state_global->fep_state);
        }

        if (MASTER(cr))
        {
            const bool constructVsites = ((vsite != nullptr) && mdrunOptions.rerunConstructVsites);
            if (constructVsites && DOMAINDECOMP(cr))
            {
                gmx_fatal(FARGS,
                          "Vsite recalculation with -rerun is not implemented with domain "
                          "decomposition, "
                          "use a single rank");
            }
            if (constructVsites)
            {
                wallcycle_start(wcycle, ewcVSITECONSTR);
                vsite->construct(state->x, ir->delta_t, state->v, state->box);
                wallcycle_stop(wcycle, ewcVSITECONSTR);
            }
        }

        if (DOMAINDECOMP(cr))
        {
            /* Repartition the domain decomposition */
            const bool bMasterState = true;
            dd_partition_system(fplog, mdlog, step, cr, bMasterState, nstglobalcomm, state_global,
                                *top_global, ir, imdSession, pull_work, state, &f, mdAtoms, &top,
                                fr, vsite, constr, nrnb, wcycle, mdrunOptions.verbose);
            shouldCheckNumberOfBondedInteractions = true;
        }

        if (MASTER(cr))
        {
            EnergyOutput::printHeader(fplog, step, t); /* can we improve the information printed here? */
        }

        if (ir->efep != efepNO)
        {
            update_mdatoms(mdatoms, state->lambda[efptMASS]);
        }

        force_flags = (GMX_FORCE_STATECHANGED | GMX_FORCE_DYNAMICBOX | GMX_FORCE_ALLFORCES
                       | GMX_FORCE_VIRIAL | // TODO: Get rid of this once #2649 is solved
                       GMX_FORCE_ENERGY | (doFreeEnergyPerturbation ? GMX_FORCE_DHDL : 0));

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            relax_shell_flexcon(fplog, cr, ms, mdrunOptions.verbose, enforcedRotation, step, ir,
                                imdSession, pull_work, bNS, force_flags, &top, constr, enerd,
                                state->natoms, state->x.arrayRefWithPadding(),
                                state->v.arrayRefWithPadding(), state->box, state->lambda,
                                &state->hist, &f.view(), force_vir, mdatoms, nrnb, wcycle, shellfc,
                                fr, runScheduleWork, t, mu_tot, vsite, ddBalanceRegionHandler);
        }
        else
        {
            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            Awh*       awh = nullptr;
            gmx_edsam* ed  = nullptr;
            do_force(fplog, cr, ms, ir, awh, enforcedRotation, imdSession, pull_work, step, nrnb,
                     wcycle, &top, state->box, state->x.arrayRefWithPadding(), &state->hist,
                     &f.view(), force_vir, mdatoms, enerd, state->lambda, fr, runScheduleWork,
                     vsite, mu_tot, t, ed, GMX_FORCE_NS | force_flags, ddBalanceRegionHandler);
        }

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t.
         */
        {
            const bool isCheckpointingStep = false;
            const bool doRerun             = false;
            const bool bSumEkinhOld        = false;
            do_md_trajectory_writing(fplog, cr, nfile, fnm, step, step_rel, t, ir, state,
                                     state_global, observablesHistory, top_global, fr, outf,
                                     energyOutput, ekind, f.view().force(), isCheckpointingStep,
                                     doRerun, isLastStep, mdrunOptions.writeConfout, bSumEkinhOld);
        }

        stopHandler->setSignal();

        {
            const bool          doInterSimSignal = false;
            const bool          doIntraSimSignal = true;
            bool                bSumEkinhOld     = false;
            t_vcm*              vcm              = nullptr;
            SimulationSignaller signaller(&signals, cr, ms, doInterSimSignal, doIntraSimSignal);

            compute_globals(gstat, cr, ir, fr, ekind, makeConstArrayRef(state->x),
                            makeConstArrayRef(state->v), state->box, mdatoms, nrnb, vcm, wcycle,
                            enerd, nullptr, nullptr, nullptr, nullptr, constr, &signaller,
                            state->box, &totalNumberOfBondedInteractions, &bSumEkinhOld,
                            CGLO_GSTAT | CGLO_ENERGY
                                    | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS
                                                                             : 0));
            checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions, top_global,
                                            &top, makeConstArrayRef(state->x), state->box,
                                            &shouldCheckNumberOfBondedInteractions);
        }

        {
            gmx::HostVector<gmx::RVec>     fglobal(top_global->natoms);
            gmx::ArrayRef<gmx::RVec>       ftemp;
            gmx::ArrayRef<const gmx::RVec> flocal = f.view().force();
            if (DOMAINDECOMP(cr))
            {
                ftemp = gmx::makeArrayRef(fglobal);
                dd_collect_vec(cr->dd, state->ddp_count, state->ddp_count_cg_gl, state->cg_gl,
                               flocal, ftemp);
            }
            else
            {
                ftemp = f.view().force();
            }

            if (MASTER(cr))
            {
                MimicCommunicator::sendEnergies(enerd->term[F_EPOT]);
                MimicCommunicator::sendForces(ftemp, state_global->natoms);
            }
        }


        /* Note: this is OK, but there are some numerical precision issues with using the convergence of
           the virial that should probably be addressed eventually. state->veta has better properies,
           but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
           generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

        if (ir->efep != efepNO)
        {
            /* Sum up the foreign energy and dhdl terms for md and sd.
               Currently done every step so that dhdl is correct in the .edr */
            accumulateKineticLambdaComponents(enerd, state->lambda, *ir->fepvals);
        }

        /* Output stuff */
        if (MASTER(cr))
        {
            const bool bCalcEnerStep = true;
            energyOutput.addDataAtEnergyStep(
                    doFreeEnergyPerturbation, bCalcEnerStep, t, mdatoms->tmass, enerd, ir->fepvals,
                    ir->expandedvals, state->box,
                    PTCouplingArrays({ state->boxv, state->nosehoover_xi, state->nosehoover_vxi,
                                       state->nhpres_xi, state->nhpres_vxi }),
                    state->fep_state, shake_vir, force_vir, total_vir, pres, ekind, mu_tot, constr);

            const bool do_ene = true;
            const bool do_log = true;
            Awh*       awh    = nullptr;
            const bool do_dr  = ir->nstdisreout != 0;
            const bool do_or  = ir->nstorireout != 0;

            EnergyOutput::printAnnealingTemperatures(do_log ? fplog : nullptr, groups, &(ir->opts));
            energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), do_ene, do_dr, do_or,
                                               do_log ? fplog : nullptr, step, t, fr->fcdata.get(), awh);

            if (do_per_step(step, ir->nstlog))
            {
                if (fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }

        /* Print the remaining wall clock time for the run */
        if (isMasterSimMasterRank(ms, MASTER(cr)) && (mdrunOptions.verbose || gmx_got_usr_signal()))
        {
            if (shellfc)
            {
                fprintf(stderr, "\n");
            }
            print_time(stderr, walltime_accounting, step, ir, cr);
        }

        cycles = wallcycle_stop(wcycle, ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        /* increase the MD step number */
        step++;
        step_rel++;
    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end_time(walltime_accounting);

    if (MASTER(cr))
    {
        MimicCommunicator::finalize();
    }

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    done_mdoutf(outf);

    done_shellfc(fplog, shellfc, step_rel);

    walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);
}
