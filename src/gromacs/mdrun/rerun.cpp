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
 * \brief Implements the loop for simulation reruns
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
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
#include "gromacs/domdec/localtopologychecker.h"
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
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mimic/utilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/output.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

#include "legacysimulator.h"
#include "replicaexchange.h"
#include "shellfc.h"

using gmx::SimulationSignaller;
using gmx::VirtualSitesHandler;

/*! \brief Copy the state from \p rerunFrame to \p globalState and, if requested, construct vsites
 *
 * \param[in]     rerunFrame      The trajectory frame to compute energy/forces for
 * \param[in,out] globalState     The global state container
 * \param[in]     constructVsites When true, vsite coordinates are constructed
 * \param[in]     vsite           Vsite setup, can be nullptr when \p constructVsites = false
 */
static void prepareRerunState(const t_trxframe&          rerunFrame,
                              t_state*                   globalState,
                              bool                       constructVsites,
                              const VirtualSitesHandler* vsite)
{
    auto x      = makeArrayRef(globalState->x);
    auto rerunX = arrayRefFromArray(reinterpret_cast<gmx::RVec*>(rerunFrame.x), globalState->natoms);
    std::copy(rerunX.begin(), rerunX.end(), x.begin());
    copy_mat(rerunFrame.box, globalState->box);

    if (constructVsites)
    {
        GMX_ASSERT(vsite, "Need valid vsite for constructing vsites");

        vsite->construct(globalState->x, globalState->v, globalState->box, gmx::VSiteOperation::PositionsAndVelocities);
    }
}

void gmx::LegacySimulator::do_rerun()
{
    // TODO Historically, the EM and MD "integrators" used different
    // names for the t_inputrec *parameter, but these must have the
    // same name, now that it's a member of a struct. We use this ir
    // alias to avoid a large ripple of nearly useless changes.
    // t_inputrec is being replaced by IMdpOptionsProvider, so this
    // will go away eventually.
    const t_inputrec* ir = inputrec;
    double            t;
    bool              isLastStep               = false;
    bool              doFreeEnergyPerturbation = false;
    unsigned int      force_flags;
    tensor            force_vir, shake_vir, total_vir, pres;
    t_trxstatus*      status = nullptr;
    rvec              mu_tot;
    t_trxframe        rerun_fr;
    ForceBuffers      f;
    gmx_global_stat_t gstat;
    gmx_shellfc_t*    shellfc;

    double cycles;

    SimulationSignals signals;
    // Most global communnication stages don't propagate mdrun
    // signals, and will use this object to achieve that.
    SimulationSignaller nullSignaller(nullptr, nullptr, nullptr, false, false);

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(
                    "Note that it is planned that the command gmx mdrun -rerun will "
                    "be available in a different form in a future version of GROMACS, "
                    "e.g. gmx rerun -f.");

    if (ir->efep != FreeEnergyPerturbationType::No
        && (mdAtoms->mdatoms()->nMassPerturbed > 0 || (constr && constr->havePerturbedConstraints())))
    {
        gmx_fatal(FARGS,
                  "Perturbed masses or constraints are not supported by rerun. "
                  "Either make a .tpr without mass and constraint perturbation, "
                  "or use GROMACS 2018.4, 2018.5 or later 2018 version.");
    }
    if (ir->bExpanded)
    {
        gmx_fatal(FARGS, "Expanded ensemble not supported by rerun.");
    }
    if (ir->bSimTemp)
    {
        gmx_fatal(FARGS, "Simulated tempering not supported by rerun.");
    }
    if (ir->bDoAwh)
    {
        gmx_fatal(FARGS, "AWH not supported by rerun.");
    }
    if (replExParams.exchangeInterval > 0)
    {
        gmx_fatal(FARGS, "Replica exchange not supported by rerun.");
    }
    if (opt2bSet("-ei", nfile, fnm) || observablesHistory->edsamHistory != nullptr)
    {
        gmx_fatal(FARGS, "Essential dynamics not supported by rerun.");
    }
    if (ir->bIMD)
    {
        gmx_fatal(FARGS, "Interactive MD not supported by rerun.");
    }
    if (isMultiSim(ms))
    {
        gmx_fatal(FARGS, "Multiple simulations not supported by rerun.");
    }
    if (std::any_of(ir->opts.annealing, ir->opts.annealing + ir->opts.ngtc, [](SimulatedAnnealing i) {
            return i != SimulatedAnnealing::No;
        }))
    {
        gmx_fatal(FARGS, "Simulated annealing not supported by rerun.");
    }

    /* Rerun can't work if an output file name is the same as the input file name.
     * If this is the case, the user will get an error telling them what the issue is.
     */
    if (strcmp(opt2fn("-rerun", nfile, fnm), opt2fn("-o", nfile, fnm)) == 0
        || strcmp(opt2fn("-rerun", nfile, fnm), opt2fn("-x", nfile, fnm)) == 0)
    {
        gmx_fatal(FARGS,
                  "When using mdrun -rerun, the name of the input trajectory file "
                  "%s cannot be identical to the name of an output file (whether "
                  "given explicitly with -o or -x, or by default)",
                  opt2fn("-rerun", nfile, fnm));
    }

    /* Settings for rerun */
    {
        // TODO: Avoid changing inputrec (#3854)
        auto* nonConstInputrec               = const_cast<t_inputrec*>(inputrec);
        nonConstInputrec->nstlist            = 1;
        nonConstInputrec->nstcalcenergy      = 1;
        nonConstInputrec->nstxout_compressed = 0;
    }
    int        nstglobalcomm = 1;
    const bool bNS           = true;

    ObservablesReducer observablesReducer = observablesReducerBuilder->build();

    const SimulationGroups* groups = &top_global.groups;
    if (ir->eI == IntegrationAlgorithm::Mimic)
    {
        auto* nonConstGlobalTopology                         = const_cast<gmx_mtop_t*>(&top_global);
        nonConstGlobalTopology->intermolecularExclusionGroup = genQmmmIndices(top_global);
    }
    int*                fep_state = MAIN(cr) ? &state_global->fep_state : nullptr;
    gmx::ArrayRef<real> lambda    = MAIN(cr) ? state_global->lambda : gmx::ArrayRef<real>();
    initialize_lambdas(fplog,
                       ir->efep,
                       ir->bSimTemp,
                       *ir->fepvals,
                       ir->simtempvals->temperatures,
                       gmx::arrayRefFromArray(ir->opts.ref_t, ir->opts.ngtc),
                       MAIN(cr),
                       fep_state,
                       lambda);
    const bool        simulationsShareState = false;
    gmx_mdoutf*       outf                  = init_mdoutf(fplog,
                                   nfile,
                                   fnm,
                                   mdrunOptions,
                                   cr,
                                   outputProvider,
                                   mdModulesNotifiers,
                                   ir,
                                   top_global,
                                   oenv,
                                   wcycle,
                                   StartingBehavior::NewSimulation,
                                   simulationsShareState,
                                   ms);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf),
                                   top_global,
                                   *ir,
                                   pull_work,
                                   mdoutf_get_fp_dhdl(outf),
                                   true,
                                   StartingBehavior::NewSimulation,
                                   simulationsShareState,
                                   mdModulesNotifiers);

    gstat = global_stat_init(ir);

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global,
                                 constr ? constr->numFlexibleConstraints() : 0,
                                 ir->nstcalcenergy,
                                 haveDDAtomOrdering(*cr),
                                 runScheduleWork->simulationWork.useGpuPme);

    if (haveDDAtomOrdering(*cr))
    {
        // Local state only becomes valid now.
        dd_init_local_state(*cr->dd, state_global, state);

        /* Distribute the charge groups over the nodes from the main node */
        dd_partition_system(fplog,
                            mdlog,
                            ir->init_step,
                            cr,
                            TRUE,
                            1,
                            state_global,
                            top_global,
                            *ir,
                            imdSession,
                            pull_work,
                            state,
                            &f,
                            mdAtoms,
                            top,
                            fr,
                            vsite,
                            constr,
                            nrnb,
                            nullptr,
                            FALSE);
    }
    else
    {
        state_change_natoms(state_global, state_global->natoms);
        /* Copy the pointer to the global state */
        state = state_global;

        mdAlgorithmsSetupAtomData(cr, *ir, top_global, top, fr, &f, mdAtoms, constr, vsite, shellfc);
    }

    auto* mdatoms = mdAtoms->mdatoms();
    fr->longRangeNonbondeds->updateAfterPartition(*mdatoms);

    // NOTE: The global state is no longer used at this point.
    // But state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdatoms, state->lambda[FreeEnergyPerturbationCouplingType::Mass]);

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
                        cr,
                        ir,
                        fr,
                        ekind,
                        makeConstArrayRef(state->x),
                        makeConstArrayRef(state->v),
                        state->box,
                        mdatoms,
                        nrnb,
                        vcm,
                        nullptr,
                        enerd,
                        force_vir,
                        shake_vir,
                        total_vir,
                        pres,
                        &nullSignaller,
                        state->box,
                        &bSumEkinhOld,
                        cglo_flags,
                        step,
                        &observablesReducer);
        // Clean up after pre-step use of compute_globals()
        observablesReducer.markAsReadyToReduce();
    }

    if (MAIN(cr))
    {
        fprintf(stderr,
                "starting md rerun '%s', reading coordinates from"
                " input trajectory '%s'\n\n",
                *(top_global.name),
                opt2fn("-rerun", nfile, fnm));
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
    wallcycle_start(wcycle, WallCycleCounter::Run);
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
                .appendText("Simulations has constraints. Rerun does not recalculate constraints.");
    }

    rerun_fr.natoms = 0;
    if (MAIN(cr))
    {
        isLastStep = !read_first_frame(oenv, &status, opt2fn("-rerun", nfile, fnm), &rerun_fr, TRX_NEED_X);
        if (rerun_fr.natoms != top_global.natoms)
        {
            gmx_fatal(FARGS,
                      "Number of atoms in trajectory (%d) does not match the "
                      "run input file (%d)\n",
                      rerun_fr.natoms,
                      top_global.natoms);
        }

        if (ir->pbcType != PbcType::No)
        {
            if (!rerun_fr.bBox)
            {
                gmx_fatal(FARGS,
                          "Rerun trajectory frame step %" PRId64
                          " time %f "
                          "does not contain a box, while pbc is used",
                          rerun_fr.step,
                          rerun_fr.time);
            }
            if (max_cutoff2(ir->pbcType, rerun_fr.box) < gmx::square(fr->rlist))
            {
                gmx_fatal(FARGS,
                          "Rerun trajectory frame step %" PRId64
                          " time %f "
                          "has too small box dimensions",
                          rerun_fr.step,
                          rerun_fr.time);
            }
        }
    }

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(
                    "Rerun does not report kinetic energy, total energy, temperature, virial and "
                    "pressure.");

    if (PAR(cr))
    {
        rerun_parallel_comm(cr, &rerun_fr, &isLastStep);
    }

    if (ir->pbcType != PbcType::No)
    {
        /* Set the shift vectors.
         * Necessary here when have a static box different from the tpr box.
         */
        calc_shifts(rerun_fr.box, fr->shift_vec);
    }

    auto stopHandler = stopHandlerBuilder->getStopHandlerMD(
            compat::not_null<SimulationSignal*>(&signals[eglsSTOPCOND]),
            false,
            MAIN(cr),
            ir->nstlist,
            mdrunOptions.reproducible,
            nstglobalcomm,
            mdrunOptions.maximumHoursToRun,
            ir->nstlist == 0,
            fplog,
            step,
            bNS,
            walltime_accounting);

    // we don't do counter resetting in rerun - finish will always be valid
    walltime_accounting_set_valid_finish(walltime_accounting);

    const DDBalanceRegionHandler ddBalanceRegionHandler(cr);

    /* and stop now if we should */
    isLastStep = (isLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!isLastStep)
    {
        wallcycle_start(wcycle, WallCycleCounter::Step);

        if (rerun_fr.bStep)
        {
            step     = rerun_fr.step;
            step_rel = step - ir->init_step;
        }
        if (rerun_fr.bTime)
        {
            t = rerun_fr.time;
        }
        else
        {
            t = step;
        }

        if (ir->efep != FreeEnergyPerturbationType::No && MAIN(cr))
        {
            if (rerun_fr.bLambda)
            {
                ir->fepvals->init_lambda = rerun_fr.lambda;
            }
            else
            {
                if (rerun_fr.bFepState)
                {
                    state->fep_state = rerun_fr.fep_state;
                }
            }

            state_global->lambda = currentLambdas(step, *(ir->fepvals), state->fep_state);
        }

        if (MAIN(cr))
        {
            const bool constructVsites = ((vsite != nullptr) && mdrunOptions.rerunConstructVsites);
            if (constructVsites && haveDDAtomOrdering(*cr))
            {
                gmx_fatal(FARGS,
                          "Vsite recalculation with -rerun is not implemented with domain "
                          "decomposition, "
                          "use a single rank");
            }
            prepareRerunState(rerun_fr, state_global, constructVsites, vsite);
        }

        isLastStep = isLastStep || stopHandler->stoppingAfterCurrentStep(bNS);

        if (haveDDAtomOrdering(*cr))
        {
            /* Repartition the domain decomposition */
            const bool bMainState = true;
            dd_partition_system(fplog,
                                mdlog,
                                step,
                                cr,
                                bMainState,
                                nstglobalcomm,
                                state_global,
                                top_global,
                                *ir,
                                imdSession,
                                pull_work,
                                state,
                                &f,
                                mdAtoms,
                                top,
                                fr,
                                vsite,
                                constr,
                                nrnb,
                                wcycle,
                                mdrunOptions.verbose);
        }

        if (MAIN(cr))
        {
            EnergyOutput::printHeader(fplog, step, t); /* can we improve the information printed here? */
        }

        if (ir->efep != FreeEnergyPerturbationType::No)
        {
            update_mdatoms(mdatoms, state->lambda[FreeEnergyPerturbationCouplingType::Mass]);
        }

        fr->longRangeNonbondeds->updateAfterPartition(*mdatoms);

        force_flags = (GMX_FORCE_STATECHANGED | GMX_FORCE_DYNAMICBOX | GMX_FORCE_ALLFORCES
                       | GMX_FORCE_VIRIAL | // TODO: Get rid of this once #2649 and #3400 are solved
                       GMX_FORCE_ENERGY | (doFreeEnergyPerturbation ? GMX_FORCE_DHDL : 0));

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            relax_shell_flexcon(fplog,
                                cr,
                                ms,
                                mdrunOptions.verbose,
                                enforcedRotation,
                                step,
                                ir,
                                imdSession,
                                pull_work,
                                bNS,
                                force_flags,
                                top,
                                constr,
                                enerd,
                                state->natoms,
                                state->x.arrayRefWithPadding(),
                                state->v.arrayRefWithPadding(),
                                state->box,
                                state->lambda,
                                &state->hist,
                                &f.view(),
                                force_vir,
                                *mdatoms,
                                fr->longRangeNonbondeds.get(),
                                nrnb,
                                wcycle,
                                shellfc,
                                fr,
                                runScheduleWork,
                                t,
                                mu_tot,
                                vsite,
                                ddBalanceRegionHandler);
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
            try
            {
                do_force(fplog,
                         cr,
                         ms,
                         *ir,
                         awh,
                         enforcedRotation,
                         imdSession,
                         pull_work,
                         step,
                         nrnb,
                         wcycle,
                         top,
                         state->box,
                         state->x.arrayRefWithPadding(),
                         &state->hist,
                         &f.view(),
                         force_vir,
                         mdatoms,
                         enerd,
                         state->lambda,
                         fr,
                         runScheduleWork,
                         vsite,
                         mu_tot,
                         t,
                         ed,
                         fr->longRangeNonbondeds.get(),
                         GMX_FORCE_NS | force_flags,
                         ddBalanceRegionHandler);
            }
            catch (const gmx::InternalError&)
            {
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendText(
                                "Continuing with next frame after catching invalid force in "
                                "previous frame");
            };
        }

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t.
         */
        {
            const bool isCheckpointingStep = false;
            const bool doRerun             = true;
            do_md_trajectory_writing(fplog,
                                     cr,
                                     nfile,
                                     fnm,
                                     step,
                                     step_rel,
                                     t,
                                     ir,
                                     state,
                                     state_global,
                                     observablesHistory,
                                     top_global,
                                     fr,
                                     outf,
                                     energyOutput,
                                     ekind,
                                     f.view().force(),
                                     isCheckpointingStep,
                                     doRerun,
                                     isLastStep,
                                     mdrunOptions.writeConfout,
                                     EkindataState::NotUsed);
        }

        stopHandler->setSignal();

        {
            const bool          doInterSimSignal = false;
            const bool          doIntraSimSignal = true;
            bool                bSumEkinhOld     = false;
            t_vcm*              vcm              = nullptr;
            SimulationSignaller signaller(&signals, cr, ms, doInterSimSignal, doIntraSimSignal);

            int cglo_flags = CGLO_GSTAT | CGLO_ENERGY;
            compute_globals(gstat,
                            cr,
                            ir,
                            fr,
                            ekind,
                            makeConstArrayRef(state->x),
                            makeConstArrayRef(state->v),
                            state->box,
                            mdatoms,
                            nrnb,
                            vcm,
                            wcycle,
                            enerd,
                            force_vir,
                            shake_vir,
                            total_vir,
                            pres,
                            &signaller,
                            state->box,
                            &bSumEkinhOld,
                            cglo_flags,
                            step,
                            &observablesReducer);
            // Clean up after pre-step use of compute_globals()
            observablesReducer.markAsReadyToReduce();
        }

        /* Note: this is OK, but there are some numerical precision issues with using the convergence of
           the virial that should probably be addressed eventually. state->veta has better properies,
           but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
           generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

        /* Output stuff */
        if (MAIN(cr))
        {
            const bool bCalcEnerStep = true;
            energyOutput.addDataAtEnergyStep(doFreeEnergyPerturbation,
                                             bCalcEnerStep,
                                             t,
                                             mdatoms->tmass,
                                             enerd,
                                             ir->fepvals.get(),
                                             state->box,
                                             PTCouplingArrays({ state->boxv,
                                                                state->nosehoover_xi,
                                                                state->nosehoover_vxi,
                                                                state->nhpres_xi,
                                                                state->nhpres_vxi }),
                                             state->fep_state,
                                             total_vir,
                                             pres,
                                             ekind,
                                             mu_tot,
                                             constr);

            const bool do_ene = true;
            const bool do_log = true;
            Awh*       awh    = nullptr;
            const bool do_dr  = ir->nstdisreout != 0;
            const bool do_or  = ir->nstorireout != 0;

            EnergyOutput::printAnnealingTemperatures(do_log ? fplog : nullptr, groups, &(ir->opts));
            energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                               do_ene,
                                               do_dr,
                                               do_or,
                                               do_log ? fplog : nullptr,
                                               step,
                                               t,
                                               fr->fcdata.get(),
                                               awh);

            if (ir->bPull)
            {
                pull_print_output(pull_work, step, t);
            }

            if (do_per_step(step, ir->nstlog))
            {
                if (fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }

        /* Print the remaining wall clock time for the run */
        if (isMainSimMainRank(ms, MAIN(cr)) && (mdrunOptions.verbose || gmx_got_usr_signal()))
        {
            if (shellfc)
            {
                fprintf(stderr, "\n");
            }
            print_time(stderr, walltime_accounting, step, ir, cr);
        }

        /* Ion/water position swapping.
         * Not done in last step since trajectory writing happens before this call
         * in the MD loop and exchanges would be lost anyway. */
        if ((ir->eSwapCoords != SwapType::No) && (step > 0) && !isLastStep
            && do_per_step(step, ir->swap->nstswap))
        {
            const bool doRerun = true;
            do_swapcoords(cr,
                          step,
                          t,
                          ir,
                          swap,
                          wcycle,
                          rerun_fr.x,
                          rerun_fr.box,
                          MAIN(cr) && mdrunOptions.verbose,
                          doRerun);
        }

        if (MAIN(cr))
        {
            /* read next frame from input trajectory */
            isLastStep = !read_next_frame(oenv, status, &rerun_fr);
        }

        if (PAR(cr))
        {
            rerun_parallel_comm(cr, &rerun_fr, &isLastStep);
        }

        cycles = wallcycle_stop(wcycle, WallCycleCounter::Step);
        if (haveDDAtomOrdering(*cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        if (!rerun_fr.bStep)
        {
            /* increase the MD step number */
            step++;
            step_rel++;
        }
        observablesReducer.markAsReadyToReduce();
    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end_time(walltime_accounting);

    if (MAIN(cr))
    {
        close_trx(status);
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
