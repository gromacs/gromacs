/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011-2019,2020,2021, by the GROMACS development team, led by
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
 * \brief Implements the integrator for normal molecular dynamics simulations
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <memory>
#include <numeric>

#include "gromacs/applied_forces/awh/awh.h"
#include "gromacs/applied_forces/awh/read_params.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/domdec/localtopologychecker.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/coupling.h"
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
#include "gromacs/mdlib/update_constrain_gpu.h"
#include "gromacs/mdlib/update_vv.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrunutility/freeenergy.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/fcdata.h"
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
#include "gromacs/mdtypes/pullhistory.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/modularsimulator/energydata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
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
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "legacysimulator.h"
#include "replicaexchange.h"
#include "shellfc.h"

using gmx::SimulationSignaller;

void gmx::LegacySimulator::do_md()
{
    // TODO Historically, the EM and MD "integrators" used different
    // names for the t_inputrec *parameter, but these must have the
    // same name, now that it's a member of a struct. We use this ir
    // alias to avoid a large ripple of nearly useless changes.
    // t_inputrec is being replaced by IMdpOptionsProvider, so this
    // will go away eventually.
    const t_inputrec* ir = inputrec;

    double       t, t0 = ir->init_t;
    gmx_bool     bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
    gmx_bool     bNS = FALSE, bNStList, bStopCM, bFirstStep, bInitStep, bLastStep = FALSE;
    gmx_bool     bDoExpanded = FALSE;
    gmx_bool     do_ene, do_log, do_verbose;
    gmx_bool     bMasterState;
    unsigned int force_flags;
    tensor force_vir = { { 0 } }, shake_vir = { { 0 } }, total_vir = { { 0 } }, pres = { { 0 } };
    int    i, m;
    rvec   mu_tot;
    matrix pressureCouplingMu, M;
    gmx_repl_ex_t     repl_ex = nullptr;
    gmx_global_stat_t gstat;
    gmx_shellfc_t*    shellfc;
    gmx_bool          bSumEkinhOld, bDoReplEx, bExchanged, bNeedRepartition;
    gmx_bool          bTrotter;
    real              dvdl_constr;
    std::vector<RVec> cbuf;
    matrix            lastbox;
    int               lamnew = 0;
    /* for FEP */
    double    cycles;
    real      saved_conserved_quantity = 0;
    real      last_ekin                = 0;
    t_extmass MassQ;
    char      sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];

    /* PME load balancing data for GPU kernels */
    gmx_bool bPMETune         = FALSE;
    gmx_bool bPMETunePrinting = FALSE;

    bool bInteractiveMDstep = false;

    SimulationSignals signals;
    // Most global communnication stages don't propagate mdrun
    // signals, and will use this object to achieve that.
    SimulationSignaller nullSignaller(nullptr, nullptr, nullptr, false, false);

    if (!mdrunOptions.writeConfout)
    {
        // This is on by default, and the main known use case for
        // turning it off is for convenience in benchmarking, which is
        // something that should not show up in the general user
        // interface.
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendText(
                        "The -noconfout functionality is deprecated, and may be removed in a "
                        "future version.");
    }

    /* md-vv uses averaged full step velocities for T-control
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bTrotter = (EI_VV(ir->eI)
                && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir) || inputrecNvtTrotter(ir)));

    const bool bRerunMD = false;

    int nstglobalcomm = computeGlobalCommunicationPeriod(mdlog, ir, cr);
    bGStatEveryStep   = (nstglobalcomm == 1);

    const SimulationGroups* groups = &top_global.groups;

    std::unique_ptr<EssentialDynamics> ed = nullptr;
    if (opt2bSet("-ei", nfile, fnm))
    {
        /* Initialize essential dynamics sampling */
        ed = init_edsam(mdlog,
                        opt2fn_null("-ei", nfile, fnm),
                        opt2fn("-eo", nfile, fnm),
                        top_global,
                        *ir,
                        cr,
                        constr,
                        state_global,
                        observablesHistory,
                        oenv,
                        startingBehavior);
    }
    else if (observablesHistory->edsamHistory)
    {
        gmx_fatal(FARGS,
                  "The checkpoint is from a run with essential dynamics sampling, "
                  "but the current run did not specify the -ei option. "
                  "Either specify the -ei option to mdrun, or do not use this checkpoint file.");
    }

    int*                fep_state = MASTER(cr) ? &state_global->fep_state : nullptr;
    gmx::ArrayRef<real> lambda    = MASTER(cr) ? state_global->lambda : gmx::ArrayRef<real>();
    initialize_lambdas(fplog,
                       ir->efep,
                       ir->bSimTemp,
                       *ir->fepvals,
                       ir->simtempvals->temperatures,
                       gmx::arrayRefFromArray(ir->opts.ref_t, ir->opts.ngtc),
                       MASTER(cr),
                       fep_state,
                       lambda);
    Update upd(*ir, deform);
    bool   doSimulatedAnnealing = false;
    {
        // TODO: Avoid changing inputrec (#3854)
        // Simulated annealing updates the reference temperature.
        auto* nonConstInputrec = const_cast<t_inputrec*>(inputrec);
        doSimulatedAnnealing   = initSimulatedAnnealing(nonConstInputrec, &upd);
    }
    const bool useReplicaExchange = (replExParams.exchangeInterval > 0);

    t_fcdata& fcdata = *fr->fcdata;

    bool simulationsShareState = false;
    int  nstSignalComm         = nstglobalcomm;
    {
        // TODO This implementation of ensemble orientation restraints is nasty because
        // a user can't just do multi-sim with single-sim orientation restraints.
        bool usingEnsembleRestraints = (fcdata.disres->nsystems > 1) || ((ms != nullptr) && fcdata.orires);
        bool awhUsesMultiSim = (ir->bDoAwh && ir->awhParams->shareBiasMultisim() && (ms != nullptr));

        // Replica exchange, ensemble restraints and AWH need all
        // simulations to remain synchronized, so they need
        // checkpoints and stop conditions to act on the same step, so
        // the propagation of such signals must take place between
        // simulations, not just within simulations.
        // TODO: Make algorithm initializers set these flags.
        simulationsShareState = useReplicaExchange || usingEnsembleRestraints || awhUsesMultiSim;

        if (simulationsShareState)
        {
            // Inter-simulation signal communication does not need to happen
            // often, so we use a minimum of 200 steps to reduce overhead.
            const int c_minimumInterSimulationSignallingInterval = 200;
            nstSignalComm = ((c_minimumInterSimulationSignallingInterval + nstglobalcomm - 1) / nstglobalcomm)
                            * nstglobalcomm;
        }
    }

    if (startingBehavior != StartingBehavior::RestartWithAppending)
    {
        pleaseCiteCouplingAlgorithms(fplog, *ir);
    }
    gmx_mdoutf*       outf = init_mdoutf(fplog,
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
                                   startingBehavior,
                                   simulationsShareState,
                                   ms);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf),
                                   top_global,
                                   *ir,
                                   pull_work,
                                   mdoutf_get_fp_dhdl(outf),
                                   false,
                                   startingBehavior,
                                   simulationsShareState,
                                   mdModulesNotifiers);

    gstat = global_stat_init(ir);

    const auto& simulationWork     = runScheduleWork->simulationWork;
    const bool  useGpuForPme       = simulationWork.useGpuPme;
    const bool  useGpuForNonbonded = simulationWork.useGpuNonbonded;
    const bool  useGpuForUpdate    = simulationWork.useGpuUpdate;

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global,
                                 constr ? constr->numFlexibleConstraints() : 0,
                                 ir->nstcalcenergy,
                                 haveDDAtomOrdering(*cr),
                                 useGpuForPme);

    {
        double io = compute_io(ir, top_global.natoms, *groups, energyOutput.numEnergyTerms(), 1);
        if ((io > 2000) && MASTER(cr))
        {
            fprintf(stderr, "\nWARNING: This run will generate roughly %.0f Mb of data\n\n", io);
        }
    }

    ObservablesReducer observablesReducer = observablesReducerBuilder->build();

    ForceBuffers     f(simulationWork.useMts,
                   (simulationWork.useGpuFBufferOps || useGpuForUpdate) ? PinningPolicy::PinnedIfSupported
                                                                            : PinningPolicy::CannotBePinned);
    const t_mdatoms* md = mdAtoms->mdatoms();
    if (haveDDAtomOrdering(*cr))
    {
        // Local state only becomes valid now.
        dd_init_local_state(*cr->dd, state_global, state);

        /* Distribute the charge groups over the nodes from the master node */
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
        upd.updateAfterPartition(state->natoms,
                                 md->cFREEZE ? gmx::arrayRefFromArray(md->cFREEZE, md->nr)
                                             : gmx::ArrayRef<const unsigned short>(),
                                 md->cTC ? gmx::arrayRefFromArray(md->cTC, md->nr)
                                         : gmx::ArrayRef<const unsigned short>(),
                                 md->cACC ? gmx::arrayRefFromArray(md->cACC, md->nr)
                                          : gmx::ArrayRef<const unsigned short>());
        fr->longRangeNonbondeds->updateAfterPartition(*md);
    }
    else
    {
        state_change_natoms(state_global, state_global->natoms);

        /* Generate and initialize new topology */
        mdAlgorithmsSetupAtomData(cr, *ir, top_global, top, fr, &f, mdAtoms, constr, vsite, shellfc);

        upd.updateAfterPartition(state->natoms,
                                 md->cFREEZE ? gmx::arrayRefFromArray(md->cFREEZE, md->nr)
                                             : gmx::ArrayRef<const unsigned short>(),
                                 md->cTC ? gmx::arrayRefFromArray(md->cTC, md->nr)
                                         : gmx::ArrayRef<const unsigned short>(),
                                 md->cACC ? gmx::arrayRefFromArray(md->cACC, md->nr)
                                          : gmx::ArrayRef<const unsigned short>());
        fr->longRangeNonbondeds->updateAfterPartition(*md);
    }

    std::unique_ptr<UpdateConstrainGpu> integrator;

    StatePropagatorDataGpu* stateGpu = fr->stateGpu;

    // TODO: the assertions below should be handled by UpdateConstraintsBuilder.
    if (useGpuForUpdate)
    {
        GMX_RELEASE_ASSERT(!haveDDAtomOrdering(*cr) || ddUsesUpdateGroups(*cr->dd)
                                   || constr == nullptr || constr->numConstraintsTotal() == 0,
                           "Constraints in domain decomposition are only supported with update "
                           "groups if using GPU update.\n");
        GMX_RELEASE_ASSERT(ir->eConstrAlg != ConstraintAlgorithm::Shake || constr == nullptr
                                   || constr->numConstraintsTotal() == 0,
                           "SHAKE is not supported with GPU update.");
        GMX_RELEASE_ASSERT(useGpuForPme || (useGpuForNonbonded && simulationWork.useGpuXBufferOps),
                           "Either PME or short-ranged non-bonded interaction tasks must run on "
                           "the GPU to use GPU update.\n");
        GMX_RELEASE_ASSERT(ir->eI == IntegrationAlgorithm::MD,
                           "Only the md integrator is supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(
                ir->etc != TemperatureCoupling::NoseHoover,
                "Nose-Hoover temperature coupling is not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(
                ir->epc == PressureCoupling::No || ir->epc == PressureCoupling::ParrinelloRahman
                        || ir->epc == PressureCoupling::Berendsen || ir->epc == PressureCoupling::CRescale,
                "Only Parrinello-Rahman, Berendsen, and C-rescale pressure coupling are supported "
                "with the GPU update.\n");
        GMX_RELEASE_ASSERT(!md->haveVsites,
                           "Virtual sites are not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(ed == nullptr,
                           "Essential dynamics is not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(!ir->bPull || !pull_have_constraint(*ir->pull),
                           "Constraints pulling is not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(fcdata.orires == nullptr,
                           "Orientation restraints are not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(
                ir->efep == FreeEnergyPerturbationType::No
                        || (!haveFepPerturbedMasses(top_global) && !havePerturbedConstraints(top_global)),
                "Free energy perturbation of masses and constraints are not supported with the GPU "
                "update.");

        if (constr != nullptr && constr->numConstraintsTotal() > 0)
        {
            GMX_LOG(mdlog.info)
                    .asParagraph()
                    .appendText("Updating coordinates and applying constraints on the GPU.");
        }
        else
        {
            GMX_LOG(mdlog.info).asParagraph().appendText("Updating coordinates on the GPU.");
        }
        GMX_RELEASE_ASSERT(fr->deviceStreamManager != nullptr,
                           "Device stream manager should be initialized in order to use GPU "
                           "update-constraints.");
        GMX_RELEASE_ASSERT(
                fr->deviceStreamManager->streamIsValid(gmx::DeviceStreamType::UpdateAndConstraints),
                "Update stream should be initialized in order to use GPU "
                "update-constraints.");
        integrator = std::make_unique<UpdateConstrainGpu>(
                *ir,
                top_global,
                ekind->ngtc,
                fr->deviceStreamManager->context(),
                fr->deviceStreamManager->stream(gmx::DeviceStreamType::UpdateAndConstraints),
                wcycle);

        stateGpu->setXUpdatedOnDeviceEvent(integrator->xUpdatedOnDeviceEvent());

        integrator->setPbc(PbcType::Xyz, state->box);
    }

    if (useGpuForPme || simulationWork.useGpuXBufferOps || useGpuForUpdate)
    {
        changePinningPolicy(&state->x, PinningPolicy::PinnedIfSupported);
    }
    if (useGpuForUpdate)
    {
        changePinningPolicy(&state->v, PinningPolicy::PinnedIfSupported);
    }

    // NOTE: The global state is no longer used at this point.
    // But state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdAtoms->mdatoms(), state->lambda[FreeEnergyPerturbationCouplingType::Mass]);

    if (ir->bExpanded)
    {
        /* Check nstexpanded here, because the grompp check was broken */
        if (ir->expandedvals->nstexpanded % ir->nstcalcenergy != 0)
        {
            gmx_fatal(FARGS,
                      "With expanded ensemble, nstexpanded should be a multiple of nstcalcenergy");
        }
        init_expanded_ensemble(startingBehavior != StartingBehavior::NewSimulation, ir, state->dfhist);
    }

    if (MASTER(cr))
    {
        EnergyData::initializeEnergyHistory(startingBehavior, observablesHistory, &energyOutput);
    }

    preparePrevStepPullCom(ir,
                           pull_work,
                           gmx::arrayRefFromArray(md->massT, md->nr),
                           state,
                           state_global,
                           cr,
                           startingBehavior != StartingBehavior::NewSimulation);

    // TODO: Remove this by converting AWH into a ForceProvider
    auto awh = prepareAwhModule(fplog,
                                *ir,
                                state_global,
                                cr,
                                ms,
                                startingBehavior != StartingBehavior::NewSimulation,
                                shellfc != nullptr,
                                opt2fn("-awh", nfile, fnm),
                                pull_work);

    if (useReplicaExchange && MASTER(cr))
    {
        repl_ex = init_replica_exchange(fplog, ms, top_global.natoms, ir, replExParams);
    }
    /* PME tuning is only supported in the Verlet scheme, with PME for
     * Coulomb. It is not supported with only LJ PME. */
    bPMETune = (mdrunOptions.tunePme && EEL_PME(fr->ic->eeltype) && !mdrunOptions.reproducible
                && ir->cutoff_scheme != CutoffScheme::Group);

    pme_load_balancing_t* pme_loadbal = nullptr;
    if (bPMETune)
    {
        pme_loadbal_init(
                &pme_loadbal, cr, mdlog, *ir, state->box, *fr->ic, *fr->nbv, fr->pmedata, fr->nbv->useGpu());
    }

    if (!ir->bContinuation)
    {
        if (state->flags & enumValueToBitMask(StateEntry::V))
        {
            auto v = makeArrayRef(state->v);
            /* Set the velocities of vsites, shells and frozen atoms to zero */
            for (i = 0; i < md->homenr; i++)
            {
                if (md->ptype[i] == ParticleType::Shell)
                {
                    clear_rvec(v[i]);
                }
                else if (md->cFREEZE)
                {
                    for (m = 0; m < DIM; m++)
                    {
                        if (ir->opts.nFreeze[md->cFREEZE[i]][m])
                        {
                            v[i][m] = 0;
                        }
                    }
                }
            }
        }

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog,
                               constr,
                               ir,
                               md->nr,
                               md->homenr,
                               state->x.arrayRefWithPadding(),
                               state->v.arrayRefWithPadding(),
                               state->box,
                               state->lambda[FreeEnergyPerturbationCouplingType::Bonded]);
        }
    }

    const int nstfep = computeFepPeriod(*ir, replExParams);

    /* Be REALLY careful about what flags you set here. You CANNOT assume
     * this is the first step, since we might be restarting from a checkpoint,
     * and in that case we should not do any modifications to the state.
     */
    bStopCM = (ir->comm_mode != ComRemovalAlgorithm::No && !ir->bContinuation);

    // When restarting from a checkpoint, it can be appropriate to
    // initialize ekind from quantities in the checkpoint. Otherwise,
    // compute_globals must initialize ekind before the simulation
    // starts/restarts. However, only the master rank knows what was
    // found in the checkpoint file, so we have to communicate in
    // order to coordinate the restart.
    //
    // TODO Consider removing this communication if/when checkpoint
    // reading directly follows .tpr reading, because all ranks can
    // agree on hasReadEkinState at that time.
    bool hasReadEkinState = MASTER(cr) ? state_global->ekinstate.hasReadEkinState : false;
    if (PAR(cr))
    {
        gmx_bcast(sizeof(hasReadEkinState), &hasReadEkinState, cr->mpi_comm_mygroup);
    }
    if (hasReadEkinState)
    {
        restore_ekinstate_from_state(cr, ekind, &state_global->ekinstate);
    }

    unsigned int cglo_flags =
            (CGLO_TEMPERATURE | CGLO_GSTAT | (EI_VV(ir->eI) ? CGLO_PRESSURE : 0)
             | (EI_VV(ir->eI) ? CGLO_CONSTRAINT : 0) | (hasReadEkinState ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;

    t_vcm vcm(top_global.groups, *ir);
    reportComRemovalInfo(fplog, vcm);

    int64_t step     = ir->init_step;
    int64_t step_rel = 0;

    /* To minimize communication, compute_globals computes the COM velocity
     * and the kinetic energy for the velocities without COM motion removed.
     * Thus to get the kinetic energy without the COM contribution, we need
     * to call compute_globals twice.
     */
    for (int cgloIteration = 0; cgloIteration < (bStopCM ? 2 : 1); cgloIteration++)
    {
        unsigned int cglo_flags_iteration = cglo_flags;
        if (bStopCM && cgloIteration == 0)
        {
            cglo_flags_iteration |= CGLO_STOPCM;
            cglo_flags_iteration &= ~CGLO_TEMPERATURE;
        }
        compute_globals(gstat,
                        cr,
                        ir,
                        fr,
                        ekind,
                        makeConstArrayRef(state->x),
                        makeConstArrayRef(state->v),
                        state->box,
                        md,
                        nrnb,
                        &vcm,
                        nullptr,
                        enerd,
                        force_vir,
                        shake_vir,
                        total_vir,
                        pres,
                        &nullSignaller,
                        state->box,
                        &bSumEkinhOld,
                        cglo_flags_iteration,
                        step,
                        &observablesReducer);
        // Clean up after pre-step use of compute_globals()
        observablesReducer.markAsReadyToReduce();

        if (cglo_flags_iteration & CGLO_STOPCM)
        {
            /* At initialization, do not pass x with acceleration-correction mode
             * to avoid (incorrect) correction of the initial coordinates.
             */
            auto x = (vcm.mode == ComRemovalAlgorithm::LinearAccelerationCorrection)
                             ? ArrayRef<RVec>()
                             : makeArrayRef(state->x);
            process_and_stopcm_grp(fplog, &vcm, *md, x, makeArrayRef(state->v));
            inc_nrnb(nrnb, eNR_STOPCM, md->homenr);
        }
    }
    if (ir->eI == IntegrationAlgorithm::VVAK)
    {
        /* a second call to get the half step temperature initialized as well */
        /* we do the same call as above, but turn the pressure off -- internally to
           compute_globals, this is recognized as a velocity verlet half-step
           kinetic energy calculation.  This minimized excess variables, but
           perhaps loses some logic?*/

        compute_globals(gstat,
                        cr,
                        ir,
                        fr,
                        ekind,
                        makeConstArrayRef(state->x),
                        makeConstArrayRef(state->v),
                        state->box,
                        md,
                        nrnb,
                        &vcm,
                        nullptr,
                        enerd,
                        force_vir,
                        shake_vir,
                        total_vir,
                        pres,
                        &nullSignaller,
                        state->box,
                        &bSumEkinhOld,
                        cglo_flags & ~CGLO_PRESSURE,
                        step,
                        &observablesReducer);
        // Clean up after pre-step use of compute_globals()
        observablesReducer.markAsReadyToReduce();
    }

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (startingBehavior == StartingBehavior::NewSimulation)
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants
       for non-trotter temperature control */
    auto trotter_seq = init_npt_vars(ir, state, &MassQ, bTrotter);

    if (MASTER(cr))
    {
        if (!ir->bContinuation)
        {
            if (constr && ir->eConstrAlg == ConstraintAlgorithm::Lincs)
            {
                fprintf(fplog,
                        "RMS relative constraint deviation after constraining: %.2e\n",
                        constr->rmsd());
            }
            if (EI_STATE_VELOCITY(ir->eI))
            {
                real temp = enerd->term[F_TEMP];
                if (ir->eI != IntegrationAlgorithm::VV)
                {
                    /* Result of Ekin averaged over velocities of -half
                     * and +half step, while we only have -half step here.
                     */
                    temp *= 2;
                }
                fprintf(fplog, "Initial temperature: %g K\n", temp);
            }
        }

        char tbuf[20];
        fprintf(stderr, "starting mdrun '%s'\n", *(top_global.name));
        if (ir->nsteps >= 0)
        {
            sprintf(tbuf, "%8.1f", (ir->init_step + ir->nsteps) * ir->delta_t);
        }
        else
        {
            sprintf(tbuf, "%s", "infinite");
        }
        if (ir->init_step > 0)
        {
            fprintf(stderr,
                    "%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                    gmx_step_str(ir->init_step + ir->nsteps, sbuf),
                    tbuf,
                    gmx_step_str(ir->init_step, sbuf2),
                    ir->init_step * ir->delta_t);
        }
        else
        {
            fprintf(stderr, "%s steps, %s ps.\n", gmx_step_str(ir->nsteps, sbuf), tbuf);
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

    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bInitStep        = startingBehavior == StartingBehavior::NewSimulation || EI_VV(ir->eI);
    bSumEkinhOld     = FALSE;
    bExchanged       = FALSE;
    bNeedRepartition = FALSE;

    auto stopHandler = stopHandlerBuilder->getStopHandlerMD(
            compat::not_null<SimulationSignal*>(&signals[eglsSTOPCOND]),
            simulationsShareState,
            MASTER(cr),
            ir->nstlist,
            mdrunOptions.reproducible,
            nstSignalComm,
            mdrunOptions.maximumHoursToRun,
            ir->nstlist == 0,
            fplog,
            step,
            bNS,
            walltime_accounting);

    auto checkpointHandler = std::make_unique<CheckpointHandler>(
            compat::make_not_null<SimulationSignal*>(&signals[eglsCHKPT]),
            simulationsShareState,
            ir->nstlist == 0,
            MASTER(cr),
            mdrunOptions.writeConfout,
            mdrunOptions.checkpointOptions.period);

    const bool resetCountersIsLocal = true;
    auto       resetHandler         = std::make_unique<ResetHandler>(
            compat::make_not_null<SimulationSignal*>(&signals[eglsRESETCOUNTERS]),
            !resetCountersIsLocal,
            ir->nsteps,
            MASTER(cr),
            mdrunOptions.timingOptions.resetHalfway,
            mdrunOptions.maximumHoursToRun,
            mdlog,
            wcycle,
            walltime_accounting);

    const DDBalanceRegionHandler ddBalanceRegionHandler(cr);

    if (MASTER(cr) && isMultiSim(ms) && !useReplicaExchange)
    {
        logInitialMultisimStatus(ms, cr, mdlog, simulationsShareState, ir->nsteps, ir->init_step);
    }

    /* and stop now if we should */
    bLastStep = (bLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!bLastStep)
    {

        /* Determine if this is a neighbor search step */
        bNStList = (ir->nstlist > 0 && step % ir->nstlist == 0);

        if (bPMETune && bNStList)
        {
            // This has to be here because PME load balancing is called so early.
            // TODO: Move to after all booleans are defined.
            if (useGpuForUpdate && !bFirstStep)
            {
                stateGpu->copyCoordinatesFromGpu(state->x, AtomLocality::Local);
                stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
            }
            /* PME grid + cut-off optimization with GPUs or PME nodes */
            pme_loadbal_do(pme_loadbal,
                           cr,
                           (mdrunOptions.verbose && MASTER(cr)) ? stderr : nullptr,
                           fplog,
                           mdlog,
                           *ir,
                           fr,
                           state->box,
                           state->x,
                           wcycle,
                           step,
                           step_rel,
                           &bPMETunePrinting,
                           simulationWork.useGpuPmePpCommunication);
        }

        wallcycle_start(wcycle, WallCycleCounter::Step);

        bLastStep = (step_rel == ir->nsteps);
        t         = t0 + step * ir->delta_t;

        // TODO Refactor this, so that nstfep does not need a default value of zero
        if (ir->efep != FreeEnergyPerturbationType::No || ir->bSimTemp)
        {
            /* find and set the current lambdas */
            state->lambda = currentLambdas(step, *(ir->fepvals), state->fep_state);

            bDoExpanded = (do_per_step(step, ir->expandedvals->nstexpanded) && (ir->bExpanded)
                           && (!bFirstStep));
        }

        bDoReplEx = (useReplicaExchange && (step > 0) && !bLastStep
                     && do_per_step(step, replExParams.exchangeInterval));

        if (doSimulatedAnnealing)
        {
            // TODO: Avoid changing inputrec (#3854)
            // Simulated annealing updates the reference temperature.
            auto* nonConstInputrec = const_cast<t_inputrec*>(inputrec);
            update_annealing_target_temp(nonConstInputrec, t, &upd);
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ComRemovalAlgorithm::No && do_per_step(step, ir->nstcomm));

        /* Determine whether or not to do Neighbour Searching */
        bNS = (bFirstStep || bNStList || bExchanged || bNeedRepartition);

        /* Note that the stopHandler will cause termination at nstglobalcomm
         * steps. Since this concides with nstcalcenergy, nsttcouple and/or
         * nstpcouple steps, we have computed the half-step kinetic energy
         * of the previous step and can always output energies at the last step.
         */
        bLastStep = bLastStep || stopHandler->stoppingAfterCurrentStep(bNS);

        /* do_log triggers energy and virial calculation. Because this leads
         * to different code paths, forces can be different. Thus for exact
         * continuation we should avoid extra log output.
         * Note that the || bLastStep can result in non-exact continuation
         * beyond the last step. But we don't consider that to be an issue.
         */
        do_log     = (do_per_step(step, ir->nstlog)
                  || (bFirstStep && startingBehavior == StartingBehavior::NewSimulation) || bLastStep);
        do_verbose = mdrunOptions.verbose
                     && (step % mdrunOptions.verboseStepPrintInterval == 0 || bFirstStep || bLastStep);

        // On search steps, when doing the update on the GPU, copy
        // the coordinates and velocities to the host unless they are
        // already there (ie on the first step and after replica
        // exchange).
        if (useGpuForUpdate && bNS && !bFirstStep && !bExchanged)
        {
            stateGpu->copyVelocitiesFromGpu(state->v, AtomLocality::Local);
            stateGpu->copyCoordinatesFromGpu(state->x, AtomLocality::Local);
            stateGpu->waitVelocitiesReadyOnHost(AtomLocality::Local);
            stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
        }

        // We only need to calculate virtual velocities if we are writing them in the current step
        const bool needVirtualVelocitiesThisStep =
                (vsite != nullptr)
                && (do_per_step(step, ir->nstvout) || checkpointHandler->isCheckpointingStep());

        if (vsite != nullptr)
        {
            // Virtual sites need to be updated before domain decomposition and forces are calculated
            wallcycle_start(wcycle, WallCycleCounter::VsiteConstr);
            // md-vv calculates virtual velocities once it has full-step real velocities
            vsite->construct(state->x,
                             state->v,
                             state->box,
                             (!EI_VV(inputrec->eI) && needVirtualVelocitiesThisStep)
                                     ? VSiteOperation::PositionsAndVelocities
                                     : VSiteOperation::Positions);
            wallcycle_stop(wcycle, WallCycleCounter::VsiteConstr);
        }

        if (bNS && !(bFirstStep && ir->bContinuation))
        {
            bMasterState = FALSE;
            /* Correct the new box if it is too skewed */
            if (inputrecDynamicBox(ir))
            {
                if (correct_box(fplog, step, state->box))
                {
                    bMasterState = TRUE;
                }
            }
            // If update is offloaded, and the box was changed either
            // above or in a replica exchange on the previous step,
            // the GPU Update object should be informed
            if (useGpuForUpdate && (bMasterState || bExchanged))
            {
                integrator->setPbc(PbcType::Xyz, state->box);
            }
            if (haveDDAtomOrdering(*cr) && bMasterState)
            {
                dd_collect_state(cr->dd, state, state_global);
            }

            if (haveDDAtomOrdering(*cr))
            {
                /* Repartition the domain decomposition */
                dd_partition_system(fplog,
                                    mdlog,
                                    step,
                                    cr,
                                    bMasterState,
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
                                    do_verbose && !bPMETunePrinting);
                upd.updateAfterPartition(state->natoms,
                                         md->cFREEZE ? gmx::arrayRefFromArray(md->cFREEZE, md->nr)
                                                     : gmx::ArrayRef<const unsigned short>(),
                                         md->cTC ? gmx::arrayRefFromArray(md->cTC, md->nr)
                                                 : gmx::ArrayRef<const unsigned short>(),
                                         md->cACC ? gmx::arrayRefFromArray(md->cACC, md->nr)
                                                  : gmx::ArrayRef<const unsigned short>());
                fr->longRangeNonbondeds->updateAfterPartition(*md);
            }
        }

        // Allocate or re-size GPU halo exchange object, if necessary
        if (bNS && simulationWork.havePpDomainDecomposition && simulationWork.useGpuHaloExchange)
        {
            GMX_RELEASE_ASSERT(fr->deviceStreamManager != nullptr,
                               "GPU device manager has to be initialized to use GPU "
                               "version of halo exchange.");
            constructGpuHaloExchange(mdlog, *cr, *fr->deviceStreamManager, wcycle);
        }

        if (MASTER(cr) && do_log)
        {
            gmx::EnergyOutput::printHeader(
                    fplog, step, t); /* can we improve the information printed here? */
        }

        if (ir->efep != FreeEnergyPerturbationType::No)
        {
            update_mdatoms(mdAtoms->mdatoms(), state->lambda[FreeEnergyPerturbationCouplingType::Mass]);
        }

        if (bExchanged)
        {
            /* We need the kinetic energy at minus the half step for determining
             * the full step kinetic energy and possibly for T-coupling.*/
            /* This may not be quite working correctly yet . . . . */
            int cglo_flags = CGLO_GSTAT | CGLO_TEMPERATURE;
            compute_globals(gstat,
                            cr,
                            ir,
                            fr,
                            ekind,
                            makeConstArrayRef(state->x),
                            makeConstArrayRef(state->v),
                            state->box,
                            md,
                            nrnb,
                            &vcm,
                            wcycle,
                            enerd,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            &nullSignaller,
                            state->box,
                            &bSumEkinhOld,
                            cglo_flags,
                            step,
                            &observablesReducer);
        }
        clear_mat(force_vir);

        checkpointHandler->decideIfCheckpointingThisStep(bNS, bFirstStep, bLastStep);

        /* Determine the energy and pressure:
         * at nstcalcenergy steps and at energy output steps (set below).
         */
        if (EI_VV(ir->eI) && (!bInitStep))
        {
            bCalcEnerStep = do_per_step(step, ir->nstcalcenergy);
            bCalcVir      = bCalcEnerStep
                       || (ir->epc != PressureCoupling::No
                           && (do_per_step(step, ir->nstpcouple) || do_per_step(step - 1, ir->nstpcouple)));
        }
        else
        {
            bCalcEnerStep = do_per_step(step, ir->nstcalcenergy);
            bCalcVir      = bCalcEnerStep
                       || (ir->epc != PressureCoupling::No && do_per_step(step, ir->nstpcouple));
        }
        bCalcEner = bCalcEnerStep;

        do_ene = (do_per_step(step, ir->nstenergy) || bLastStep);

        if (do_ene || do_log || bDoReplEx)
        {
            bCalcVir  = TRUE;
            bCalcEner = TRUE;
        }

        // bCalcEner is only here for when the last step is not a mulitple of nstfep
        const bool computeDHDL = ((ir->efep != FreeEnergyPerturbationType::No || ir->bSimTemp)
                                  && (do_per_step(step, nstfep) || bCalcEner));

        /* Do we need global communication ? */
        bGStat = (bCalcVir || bCalcEner || bStopCM || do_per_step(step, nstglobalcomm)
                  || (EI_VV(ir->eI) && inputrecNvtTrotter(ir) && do_per_step(step - 1, nstglobalcomm)));

        force_flags = (GMX_FORCE_STATECHANGED | ((inputrecDynamicBox(ir)) ? GMX_FORCE_DYNAMICBOX : 0)
                       | GMX_FORCE_ALLFORCES | (bCalcVir ? GMX_FORCE_VIRIAL : 0)
                       | (bCalcEner ? GMX_FORCE_ENERGY : 0) | (computeDHDL ? GMX_FORCE_DHDL : 0));
        if (simulationWork.useMts && !do_per_step(step, ir->nstfout))
        {
            // TODO: merge this with stepWork.useOnlyMtsCombinedForceBuffer
            force_flags |= GMX_FORCE_DO_NOT_NEED_NORMAL_FORCE;
        }

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
                                *md,
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
            /* The AWH history need to be saved _before_ doing force calculations where the AWH bias
               is updated (or the AWH update will be performed twice for one step when continuing).
               It would be best to call this update function from do_md_trajectory_writing but that
               would occur after do_force. One would have to divide the update_awh function into one
               function applying the AWH force and one doing the AWH bias update. The update AWH
               bias function could then be called after do_md_trajectory_writing (then containing
               update_awh_history). The checkpointing will in the future probably moved to the start
               of the md loop which will rid of this issue. */
            if (awh && checkpointHandler->isCheckpointingStep() && MASTER(cr))
            {
                awh->updateHistory(state_global->awhHistory.get());
            }

            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            do_force(fplog,
                     cr,
                     ms,
                     *ir,
                     awh.get(),
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
                     md,
                     enerd,
                     state->lambda,
                     fr,
                     runScheduleWork,
                     vsite,
                     mu_tot,
                     t,
                     ed ? ed->getLegacyED() : nullptr,
                     fr->longRangeNonbondeds.get(),
                     (bNS ? GMX_FORCE_NS : 0) | force_flags,
                     ddBalanceRegionHandler);
        }

        // VV integrators do not need the following velocity half step
        // if it is the first step after starting from a checkpoint.
        // That is, the half step is needed on all other steps, and
        // also the first step when starting from a .tpr file.
        if (EI_VV(ir->eI))
        {
            integrateVVFirstStep(step,
                                 bFirstStep,
                                 bInitStep,
                                 startingBehavior,
                                 nstglobalcomm,
                                 ir,
                                 fr,
                                 cr,
                                 state,
                                 mdAtoms->mdatoms(),
                                 &fcdata,
                                 &MassQ,
                                 &vcm,
                                 enerd,
                                 &observablesReducer,
                                 ekind,
                                 gstat,
                                 &last_ekin,
                                 bCalcVir,
                                 total_vir,
                                 shake_vir,
                                 force_vir,
                                 pres,
                                 M,
                                 do_log,
                                 do_ene,
                                 bCalcEner,
                                 bGStat,
                                 bStopCM,
                                 bTrotter,
                                 bExchanged,
                                 &bSumEkinhOld,
                                 &saved_conserved_quantity,
                                 &f,
                                 &upd,
                                 constr,
                                 &nullSignaller,
                                 trotter_seq,
                                 nrnb,
                                 fplog,
                                 wcycle);
            if (vsite != nullptr && needVirtualVelocitiesThisStep)
            {
                // Positions were calculated earlier
                wallcycle_start(wcycle, WallCycleCounter::VsiteConstr);
                vsite->construct(state->x, state->v, state->box, VSiteOperation::Velocities);
                wallcycle_stop(wcycle, WallCycleCounter::VsiteConstr);
            }
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */
        if (bDoExpanded)
        {
            /* perform extended ensemble sampling in lambda - we don't
               actually move to the new state before outputting
               statistics, but if performing simulated tempering, we
               do update the velocities and the tau_t. */
            // TODO: Avoid changing inputrec (#3854)
            // Simulated tempering updates the reference temperature.
            // Expanded ensemble without simulated tempering does not change the inputrec.
            auto* nonConstInputrec = const_cast<t_inputrec*>(inputrec);
            lamnew                 = ExpandedEnsembleDynamics(fplog,
                                              nonConstInputrec,
                                              enerd,
                                              state,
                                              &MassQ,
                                              state->fep_state,
                                              state->dfhist,
                                              step,
                                              state->v.rvec_array(),
                                              md->homenr,
                                              md->cTC ? gmx::arrayRefFromArray(md->cTC, md->nr)
                                                                      : gmx::ArrayRef<const unsigned short>());
            /* history is maintained in state->dfhist, but state_global is what is sent to trajectory and log output */
            if (MASTER(cr))
            {
                copy_df_history(state_global->dfhist, state->dfhist);
            }
        }

        // Copy coordinate from the GPU for the output/checkpointing if the update is offloaded and
        // coordinates have not already been copied for i) search or ii) CPU force tasks.
        if (useGpuForUpdate && !bNS && !runScheduleWork->domainWork.haveCpuLocalForceWork
            && (do_per_step(step, ir->nstxout) || do_per_step(step, ir->nstxout_compressed)
                || checkpointHandler->isCheckpointingStep()))
        {
            stateGpu->copyCoordinatesFromGpu(state->x, AtomLocality::Local);
            stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
        }
        // Copy velocities if needed for the output/checkpointing.
        // NOTE: Copy on the search steps is done at the beginning of the step.
        if (useGpuForUpdate && !bNS
            && (do_per_step(step, ir->nstvout) || checkpointHandler->isCheckpointingStep()))
        {
            stateGpu->copyVelocitiesFromGpu(state->v, AtomLocality::Local);
            stateGpu->waitVelocitiesReadyOnHost(AtomLocality::Local);
        }
        // Copy forces for the output if the forces were reduced on the GPU (not the case on virial steps)
        // and update is offloaded hence forces are kept on the GPU for update and have not been
        // already transferred in do_force().
        // TODO: There should be an improved, explicit mechanism that ensures this copy is only executed
        //       when the forces are ready on the GPU -- the same synchronizer should be used as the one
        //       prior to GPU update.
        // TODO: When the output flags will be included in step workload, this copy can be combined with the
        //       copy call in do_force(...).
        // NOTE: The forces should not be copied here if the vsites are present, since they were modified
        //       on host after the D2H copy in do_force(...).
        if (runScheduleWork->stepWork.useGpuFBufferOps && (simulationWork.useGpuUpdate && !vsite)
            && do_per_step(step, ir->nstfout))
        {
            stateGpu->copyForcesFromGpu(f.view().force(), AtomLocality::Local);
            stateGpu->waitForcesReadyOnHost(AtomLocality::Local);
        }
        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         */
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
                                 checkpointHandler->isCheckpointingStep(),
                                 bRerunMD,
                                 bLastStep,
                                 mdrunOptions.writeConfout,
                                 bSumEkinhOld);
        /* Check if IMD step and do IMD communication, if bIMD is TRUE. */
        bInteractiveMDstep = imdSession->run(step, bNS, state->box, state->x, t);

        /* kludge -- virial is lost with restart for MTTK NPT control. Must reload (saved earlier). */
        if (startingBehavior != StartingBehavior::NewSimulation && bFirstStep
            && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir)))
        {
            copy_mat(state->svir_prev, shake_vir);
            copy_mat(state->fvir_prev, force_vir);
        }

        stopHandler->setSignal();
        resetHandler->setSignal(walltime_accounting);

        if (bGStat || !PAR(cr))
        {
            /* In parallel we only have to check for checkpointing in steps
             * where we do global communication,
             *  otherwise the other nodes don't know.
             */
            checkpointHandler->setSignal(walltime_accounting);
        }

        /* #########   START SECOND UPDATE STEP ################# */

        /* at the start of step, randomize or scale the velocities ((if vv. Restriction of Andersen
           controlled in preprocessing */

        if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
        {
            gmx_bool bIfRandomize;
            bIfRandomize = update_randomize_velocities(ir,
                                                       step,
                                                       cr,
                                                       md->homenr,
                                                       md->cTC ? gmx::arrayRefFromArray(md->cTC, md->nr)
                                                               : gmx::ArrayRef<const unsigned short>(),
                                                       gmx::arrayRefFromArray(md->invmass, md->nr),
                                                       state->v,
                                                       &upd,
                                                       constr);
            /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
            if (constr && bIfRandomize)
            {
                constrain_velocities(constr, do_log, do_ene, step, state, nullptr, false, nullptr);
            }
        }
        /* Box is changed in update() when we do pressure coupling,
         * but we should still use the old box for energy corrections and when
         * writing it to the energy file, so it matches the trajectory files for
         * the same timestep above. Make a copy in a separate array.
         */
        copy_mat(state->box, lastbox);

        dvdl_constr = 0;

        if (!useGpuForUpdate)
        {
            wallcycle_start(wcycle, WallCycleCounter::Update);
        }
        /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
        if (bTrotter)
        {
            trotter_update(ir,
                           step,
                           ekind,
                           enerd,
                           state,
                           total_vir,
                           md->homenr,
                           md->cTC ? gmx::arrayRefFromArray(md->cTC, md->nr)
                                   : gmx::ArrayRef<const unsigned short>(),
                           gmx::arrayRefFromArray(md->invmass, md->nr),
                           &MassQ,
                           trotter_seq,
                           TrotterSequence::Three);
            /* We can only do Berendsen coupling after we have summed
             * the kinetic energy or virial. Since the happens
             * in global_state after update, we should only do it at
             * step % nstlist = 1 with bGStatEveryStep=FALSE.
             */
        }
        else
        {
            update_tcouple(step,
                           ir,
                           state,
                           ekind,
                           &MassQ,
                           md->homenr,
                           md->cTC ? gmx::arrayRefFromArray(md->cTC, md->nr)
                                   : gmx::ArrayRef<const unsigned short>());
            update_pcouple_before_coordinates(fplog, step, ir, state, pressureCouplingMu, M, bInitStep);
        }

        /* With leap-frog type integrators we compute the kinetic energy
         * at a whole time step as the average of the half-time step kinetic
         * energies of two subsequent steps. Therefore we need to compute the
         * half step kinetic energy also if we need energies at the next step.
         */
        const bool needHalfStepKineticEnergy =
                (!EI_VV(ir->eI) && (do_per_step(step + 1, nstglobalcomm) || step_rel + 1 == ir->nsteps));

        // Parrinello-Rahman requires the pressure to be availible before the update to compute
        // the velocity scaling matrix. Hence, it runs one step after the nstpcouple step.
        const bool doParrinelloRahman = (ir->epc == PressureCoupling::ParrinelloRahman
                                         && do_per_step(step + ir->nstpcouple - 1, ir->nstpcouple));

        if (EI_VV(ir->eI))
        {
            GMX_ASSERT(!useGpuForUpdate, "GPU update is not supported with VVAK integrator.");

            integrateVVSecondStep(step,
                                  ir,
                                  fr,
                                  cr,
                                  state,
                                  mdAtoms->mdatoms(),
                                  &fcdata,
                                  &MassQ,
                                  &vcm,
                                  pull_work,
                                  enerd,
                                  &observablesReducer,
                                  ekind,
                                  gstat,
                                  &dvdl_constr,
                                  bCalcVir,
                                  total_vir,
                                  shake_vir,
                                  force_vir,
                                  pres,
                                  M,
                                  lastbox,
                                  do_log,
                                  do_ene,
                                  bGStat,
                                  &bSumEkinhOld,
                                  &f,
                                  &cbuf,
                                  &upd,
                                  constr,
                                  &nullSignaller,
                                  trotter_seq,
                                  nrnb,
                                  wcycle);
        }
        else
        {
            if (useGpuForUpdate)
            {
                // On search steps, update handles to device vectors
                if (bNS && (bFirstStep || haveDDAtomOrdering(*cr) || bExchanged))
                {
                    integrator->set(stateGpu->getCoordinates(),
                                    stateGpu->getVelocities(),
                                    stateGpu->getForces(),
                                    top->idef,
                                    *md);

                    // Copy data to the GPU after buffers might have being reinitialized
                    /* The velocity copy is redundant if we had Center-of-Mass motion removed on
                     * the previous step. We don't check that now. */
                    stateGpu->copyVelocitiesToGpu(state->v, AtomLocality::Local);
                    if (bExchanged
                        || (!runScheduleWork->stepWork.haveGpuPmeOnThisRank
                            && !runScheduleWork->stepWork.useGpuXBufferOps))
                    {
                        stateGpu->copyCoordinatesToGpu(state->x, AtomLocality::Local);
                        // Coordinates are later used by the integrator running in the same stream.
                        stateGpu->consumeCoordinatesCopiedToDeviceEvent(AtomLocality::Local);
                    }
                }

                if ((simulationWork.useGpuPme && simulationWork.useCpuPmePpCommunication)
                    || (!runScheduleWork->stepWork.useGpuFBufferOps))
                {
                    // The PME forces were recieved to the host, and reduced on the CPU with the
                    // rest of the forces computed on the GPU, so the final forces have to be copied
                    // back to the GPU. Or the buffer ops were not offloaded this step, so the
                    // forces are on the host and have to be copied
                    stateGpu->copyForcesToGpu(f.view().force(), AtomLocality::Local);
                }
                const bool doTemperatureScaling =
                        (ir->etc != TemperatureCoupling::No
                         && do_per_step(step + ir->nsttcouple - 1, ir->nsttcouple));

                // This applies Leap-Frog, LINCS and SETTLE in succession
                integrator->integrate(stateGpu->getLocalForcesReadyOnDeviceEvent(
                                              runScheduleWork->stepWork, runScheduleWork->simulationWork),
                                      ir->delta_t,
                                      true,
                                      bCalcVir,
                                      shake_vir,
                                      doTemperatureScaling,
                                      ekind->tcstat,
                                      doParrinelloRahman,
                                      ir->nstpcouple * ir->delta_t,
                                      M);
            }
            else
            {
                /* With multiple time stepping we need to do an additional normal
                 * update step to obtain the virial, as the actual MTS integration
                 * using an acceleration where the slow forces are multiplied by mtsFactor.
                 * Using that acceleration would result in a virial with the slow
                 * force contribution would be a factor mtsFactor too large.
                 */
                if (simulationWork.useMts && bCalcVir && constr != nullptr)
                {
                    upd.update_for_constraint_virial(*ir,
                                                     md->homenr,
                                                     md->havePartiallyFrozenAtoms,
                                                     gmx::arrayRefFromArray(md->invmass, md->nr),
                                                     gmx::arrayRefFromArray(md->invMassPerDim, md->nr),
                                                     *state,
                                                     f.view().forceWithPadding(),
                                                     *ekind);

                    constrain_coordinates(constr,
                                          do_log,
                                          do_ene,
                                          step,
                                          state,
                                          upd.xp()->arrayRefWithPadding(),
                                          &dvdl_constr,
                                          bCalcVir,
                                          shake_vir);
                }

                ArrayRefWithPadding<const RVec> forceCombined =
                        (simulationWork.useMts && step % ir->mtsLevels[1].stepFactor == 0)
                                ? f.view().forceMtsCombinedWithPadding()
                                : f.view().forceWithPadding();
                upd.update_coords(*ir,
                                  step,
                                  md->homenr,
                                  md->havePartiallyFrozenAtoms,
                                  gmx::arrayRefFromArray(md->ptype, md->nr),
                                  gmx::arrayRefFromArray(md->invmass, md->nr),
                                  gmx::arrayRefFromArray(md->invMassPerDim, md->nr),
                                  state,
                                  forceCombined,
                                  &fcdata,
                                  ekind,
                                  M,
                                  etrtPOSITION,
                                  cr,
                                  constr != nullptr);

                wallcycle_stop(wcycle, WallCycleCounter::Update);

                constrain_coordinates(constr,
                                      do_log,
                                      do_ene,
                                      step,
                                      state,
                                      upd.xp()->arrayRefWithPadding(),
                                      &dvdl_constr,
                                      bCalcVir && !simulationWork.useMts,
                                      shake_vir);

                upd.update_sd_second_half(*ir,
                                          step,
                                          &dvdl_constr,
                                          md->homenr,
                                          gmx::arrayRefFromArray(md->ptype, md->nr),
                                          gmx::arrayRefFromArray(md->invmass, md->nr),
                                          state,
                                          cr,
                                          nrnb,
                                          wcycle,
                                          constr,
                                          do_log,
                                          do_ene);
                upd.finish_update(
                        *ir, md->havePartiallyFrozenAtoms, md->homenr, state, wcycle, constr != nullptr);
            }

            if (ir->bPull && ir->pull->bSetPbcRefToPrevStepCOM)
            {
                updatePrevStepPullCom(pull_work, state->pull_com_prev_step);
            }

            enerd->term[F_DVDL_CONSTR] += dvdl_constr;
        }

        /* ############## IF NOT VV, Calculate globals HERE  ############ */
        /* With Leap-Frog we can skip compute_globals at
         * non-communication steps, but we need to calculate
         * the kinetic energy one step before communication.
         */
        {
            // Organize to do inter-simulation signalling on steps if
            // and when algorithms require it.
            const bool doInterSimSignal = (simulationsShareState && do_per_step(step, nstSignalComm));

            if (useGpuForUpdate)
            {
                const bool coordinatesRequiredForStopCM =
                        bStopCM && (bGStat || needHalfStepKineticEnergy || doInterSimSignal)
                        && !EI_VV(ir->eI);

                // Copy coordinates when needed to stop the CM motion or for replica exchange
                if (coordinatesRequiredForStopCM || bDoReplEx)
                {
                    stateGpu->copyCoordinatesFromGpu(state->x, AtomLocality::Local);
                    stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
                }

                // Copy velocities back to the host if:
                // - Globals are computed this step (includes the energy output steps).
                // - Temperature is needed for the next step.
                // - This is a replica exchange step (even though we will only need
                //     the velocities if an exchange succeeds)
                if (bGStat || needHalfStepKineticEnergy || bDoReplEx)
                {
                    stateGpu->copyVelocitiesFromGpu(state->v, AtomLocality::Local);
                    stateGpu->waitVelocitiesReadyOnHost(AtomLocality::Local);
                }
            }

            if (bGStat || needHalfStepKineticEnergy || doInterSimSignal)
            {
                // Since we're already communicating at this step, we
                // can propagate intra-simulation signals. Note that
                // check_nstglobalcomm has the responsibility for
                // choosing the value of nstglobalcomm that is one way
                // bGStat becomes true, so we can't get into a
                // situation where e.g. checkpointing can't be
                // signalled.
                bool                doIntraSimSignal = true;
                SimulationSignaller signaller(&signals, cr, ms, doInterSimSignal, doIntraSimSignal);

                compute_globals(gstat,
                                cr,
                                ir,
                                fr,
                                ekind,
                                makeConstArrayRef(state->x),
                                makeConstArrayRef(state->v),
                                state->box,
                                md,
                                nrnb,
                                &vcm,
                                wcycle,
                                enerd,
                                force_vir,
                                shake_vir,
                                total_vir,
                                pres,
                                &signaller,
                                lastbox,
                                &bSumEkinhOld,
                                (bGStat ? CGLO_GSTAT : 0) | (!EI_VV(ir->eI) && bCalcEner ? CGLO_ENERGY : 0)
                                        | (!EI_VV(ir->eI) && bStopCM ? CGLO_STOPCM : 0)
                                        | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                                        | (!EI_VV(ir->eI) ? CGLO_PRESSURE : 0) | CGLO_CONSTRAINT,
                                step,
                                &observablesReducer);
                if (!EI_VV(ir->eI) && bStopCM)
                {
                    process_and_stopcm_grp(
                            fplog, &vcm, *md, makeArrayRef(state->x), makeArrayRef(state->v));
                    inc_nrnb(nrnb, eNR_STOPCM, md->homenr);

                    // TODO: The special case of removing CM motion should be dealt more gracefully
                    if (useGpuForUpdate)
                    {
                        // Issue #3988, #4106.
                        stateGpu->resetCoordinatesCopiedToDeviceEvent(AtomLocality::Local);
                        stateGpu->copyCoordinatesToGpu(state->x, AtomLocality::Local);
                        // Here we block until the H2D copy completes because event sync with the
                        // force kernels that use the coordinates on the next steps is not implemented
                        // (not because of a race on state->x being modified on the CPU while H2D is in progress).
                        stateGpu->waitCoordinatesCopiedToDevice(AtomLocality::Local);
                        // If the COM removal changed the velocities on the CPU, this has to be accounted for.
                        if (vcm.mode != ComRemovalAlgorithm::No)
                        {
                            stateGpu->copyVelocitiesToGpu(state->v, AtomLocality::Local);
                        }
                    }
                }
            }
        }

        /* #############  END CALC EKIN AND PRESSURE ################# */

        /* Note: this is OK, but there are some numerical precision issues with using the convergence of
           the virial that should probably be addressed eventually. state->veta has better properies,
           but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
           generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

        if (ir->efep != FreeEnergyPerturbationType::No && !EI_VV(ir->eI))
        {
            /* Sum up the foreign energy and dK/dl terms for md and sd.
               Currently done every step so that dH/dl is correct in the .edr */
            accumulateKineticLambdaComponents(enerd, state->lambda, *ir->fepvals);
        }

        bool scaleCoordinates = !useGpuForUpdate || bDoReplEx;
        update_pcouple_after_coordinates(fplog,
                                         step,
                                         ir,
                                         md->homenr,
                                         md->cFREEZE ? gmx::arrayRefFromArray(md->cFREEZE, md->nr)
                                                     : gmx::ArrayRef<const unsigned short>(),
                                         pres,
                                         force_vir,
                                         shake_vir,
                                         pressureCouplingMu,
                                         state,
                                         nrnb,
                                         upd.deform(),
                                         scaleCoordinates);

        const bool doBerendsenPressureCoupling = (inputrec->epc == PressureCoupling::Berendsen
                                                  && do_per_step(step, inputrec->nstpcouple));
        const bool doCRescalePressureCoupling  = (inputrec->epc == PressureCoupling::CRescale
                                                 && do_per_step(step, inputrec->nstpcouple));
        if (useGpuForUpdate
            && (doBerendsenPressureCoupling || doCRescalePressureCoupling || doParrinelloRahman))
        {
            integrator->scaleCoordinates(pressureCouplingMu);
            if (doCRescalePressureCoupling)
            {
                matrix pressureCouplingInvMu;
                gmx::invertBoxMatrix(pressureCouplingMu, pressureCouplingInvMu);
                integrator->scaleVelocities(pressureCouplingInvMu);
            }
            integrator->setPbc(PbcType::Xyz, state->box);
        }

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

        /* The coordinates (x) were unshifted in update */
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            bSumEkinhOld = TRUE;
        }

        if (bCalcEner)
        {
            /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

            /* use the directly determined last velocity, not actually the averaged half steps */
            if (bTrotter && ir->eI == IntegrationAlgorithm::VV)
            {
                enerd->term[F_EKIN] = last_ekin;
            }
            enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

            if (integratorHasConservedEnergyQuantity(ir))
            {
                if (EI_VV(ir->eI))
                {
                    enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + saved_conserved_quantity;
                }
                else
                {
                    enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + NPT_energy(ir, state, &MassQ);
                }
            }
            /* #########  END PREPARING EDR OUTPUT  ###########  */
        }

        /* Output stuff */
        if (MASTER(cr))
        {
            if (fplog && do_log && bDoExpanded)
            {
                /* only needed if doing expanded ensemble */
                PrintFreeEnergyInfoToFile(fplog,
                                          ir->fepvals.get(),
                                          ir->expandedvals.get(),
                                          ir->bSimTemp ? ir->simtempvals.get() : nullptr,
                                          state_global->dfhist,
                                          state->fep_state,
                                          ir->nstlog,
                                          step);
            }
            if (bCalcEner)
            {
                const bool outputDHDL = (computeDHDL && do_per_step(step, ir->fepvals->nstdhdl));

                energyOutput.addDataAtEnergyStep(outputDHDL,
                                                 bCalcEnerStep,
                                                 t,
                                                 md->tmass,
                                                 enerd,
                                                 ir->fepvals.get(),
                                                 lastbox,
                                                 PTCouplingArrays{ state->boxv,
                                                                   state->nosehoover_xi,
                                                                   state->nosehoover_vxi,
                                                                   state->nhpres_xi,
                                                                   state->nhpres_vxi },
                                                 state->fep_state,
                                                 total_vir,
                                                 pres,
                                                 ekind,
                                                 mu_tot,
                                                 constr);
            }
            else
            {
                energyOutput.recordNonEnergyStep();
            }

            gmx_bool do_dr = do_per_step(step, ir->nstdisreout);
            gmx_bool do_or = do_per_step(step, ir->nstorireout);

            if (doSimulatedAnnealing)
            {
                gmx::EnergyOutput::printAnnealingTemperatures(
                        do_log ? fplog : nullptr, groups, &(ir->opts));
            }
            if (do_log || do_ene || do_dr || do_or)
            {
                energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                                   do_ene,
                                                   do_dr,
                                                   do_or,
                                                   do_log ? fplog : nullptr,
                                                   step,
                                                   t,
                                                   fr->fcdata.get(),
                                                   awh.get());
            }
            if (do_log && ir->bDoAwh && awh->hasFepLambdaDimension())
            {
                const bool isInitialOutput = false;
                printLambdaStateToLog(fplog, state->lambda, isInitialOutput);
            }

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
        if (bDoExpanded)
        {
            /* Have to do this part _after_ outputting the logfile and the edr file */
            /* Gets written into the state at the beginning of next loop*/
            state->fep_state = lamnew;
        }
        else if (ir->bDoAwh && awh->needForeignEnergyDifferences(step))
        {
            state->fep_state = awh->fepLambdaState();
        }
        /* Print the remaining wall clock time for the run */
        if (isMasterSimMasterRank(ms, MASTER(cr)) && (do_verbose || gmx_got_usr_signal()) && !bPMETunePrinting)
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
        bNeedRepartition = FALSE;
        if ((ir->eSwapCoords != SwapType::No) && (step > 0) && !bLastStep
            && do_per_step(step, ir->swap->nstswap))
        {
            bNeedRepartition = do_swapcoords(cr,
                                             step,
                                             t,
                                             ir,
                                             swap,
                                             wcycle,
                                             as_rvec_array(state->x.data()),
                                             state->box,
                                             MASTER(cr) && mdrunOptions.verbose,
                                             bRerunMD);

            if (bNeedRepartition && haveDDAtomOrdering(*cr))
            {
                dd_collect_state(cr->dd, state, state_global);
            }
        }

        /* Replica exchange */
        bExchanged = FALSE;
        if (bDoReplEx)
        {
            bExchanged = replica_exchange(fplog, cr, ms, repl_ex, state_global, enerd, state, step, t);
        }

        if ((bExchanged || bNeedRepartition) && haveDDAtomOrdering(*cr))
        {
            dd_partition_system(fplog,
                                mdlog,
                                step,
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
                                wcycle,
                                FALSE);
            upd.updateAfterPartition(state->natoms,
                                     md->cFREEZE ? gmx::arrayRefFromArray(md->cFREEZE, md->nr)
                                                 : gmx::ArrayRef<const unsigned short>(),
                                     md->cTC ? gmx::arrayRefFromArray(md->cTC, md->nr)
                                             : gmx::ArrayRef<const unsigned short>(),
                                     md->cACC ? gmx::arrayRefFromArray(md->cACC, md->nr)
                                              : gmx::ArrayRef<const unsigned short>());
            fr->longRangeNonbondeds->updateAfterPartition(*md);
        }

        bFirstStep = FALSE;
        bInitStep  = FALSE;

        /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
        /* With all integrators, except VV, we need to retain the pressure
         * at the current step for coupling at the next step.
         */
        if ((state->flags & enumValueToBitMask(StateEntry::PressurePrevious))
            && (bGStatEveryStep || (ir->nstpcouple > 0 && step % ir->nstpcouple == 0)))
        {
            /* Store the pressure in t_state for pressure coupling
             * at the next MD step.
             */
            copy_mat(pres, state->pres_prev);
        }

        /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

        if ((membed != nullptr) && (!bLastStep))
        {
            rescale_membed(step_rel, membed, as_rvec_array(state_global->x.data()));
        }

        cycles = wallcycle_stop(wcycle, WallCycleCounter::Step);
        if (haveDDAtomOrdering(*cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        /* increase the MD step number */
        step++;
        step_rel++;
        observablesReducer.markAsReadyToReduce();

#if GMX_FAHCORE
        if (MASTER(cr))
        {
            fcReportProgress(ir->nsteps + ir->init_step, step);
        }
#endif

        resetHandler->resetCounters(
                step, step_rel, mdlog, fplog, cr, fr->nbv.get(), nrnb, fr->pmedata, pme_loadbal, wcycle, walltime_accounting);

        /* If bIMD is TRUE, the master updates the IMD energy record and sends positions to VMD client */
        imdSession->updateEnergyRecordAndSendPositionsAndEnergies(bInteractiveMDstep, step, bCalcEner);
    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end_time(walltime_accounting);

    if (simulationWork.haveSeparatePmeRank)
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0)
        {
            energyOutput.printEnergyConservation(fplog, ir->simulation_part, EI_MD(ir->eI));

            gmx::EnergyOutput::printAnnealingTemperatures(fplog, groups, &(ir->opts));
            energyOutput.printAverages(fplog, groups);
        }
    }
    done_mdoutf(outf);

    if (bPMETune)
    {
        pme_loadbal_done(pme_loadbal, fplog, mdlog, fr->nbv->useGpu());
    }

    done_shellfc(fplog, shellfc, step_rel);

    if (useReplicaExchange && MASTER(cr))
    {
        print_replica_exchange_statistics(fplog, repl_ex);
    }

    walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);

    global_stat_destroy(gstat);
}
