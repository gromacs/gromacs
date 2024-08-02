/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include <array>
#include <filesystem>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/applied_forces/awh/awh.h"
#include "gromacs/applied_forces/awh/read_params.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/pointers.h"
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
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/matrix.h"
#include "gromacs/math/paddedvector.h"
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
#include "gromacs/mdlib/mdgraph_gpu.h"
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
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/mdtypes/pull_params.h"
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
#include "gromacs/taskassignment/include/gromacs/taskassignment/decidesimulationworkload.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "legacysimulator.h"
#include "replicaexchange.h"
#include "shellfc.h"

struct gmx_mdoutf;
struct gmx_shellfc_t;
struct pme_load_balancing_t;

using gmx::SimulationSignaller;

void gmx::LegacySimulator::do_md()
{
    // TODO Historically, the EM and MD "integrators" used different
    // names for the t_inputrec *parameter, but these must have the
    // same name, now that it's a member of a struct. We use this ir
    // alias to avoid a large ripple of nearly useless changes.
    // t_inputrec is being replaced by IMdpOptionsProvider, so this
    // will go away eventually.
    const t_inputrec* ir = inputRec_;

    double       t, t0 = ir->init_t;
    gmx_bool     bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
    gmx_bool     bNS = FALSE, bNStList, bStopCM, bFirstStep, bInitStep, bLastStep = FALSE;
    gmx_bool     bDoExpanded = FALSE;
    gmx_bool     do_ene, do_log, do_verbose;
    gmx_bool     bMainState;
    unsigned int force_flags;
    tensor    force_vir = { { 0 } }, shake_vir = { { 0 } }, total_vir = { { 0 } }, pres = { { 0 } };
    int       i, m;
    rvec      mu_tot;
    Matrix3x3 pressureCouplingMu{ { 0. } }, parrinelloRahmanM{ { 0. } };
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

    if (!mdrunOptions_.writeConfout)
    {
        // This is on by default, and the main known use case for
        // turning it off is for convenience in benchmarking, which is
        // something that should not show up in the general user
        // interface.
        GMX_LOG(mdLog_.info)
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

    int nstglobalcomm = computeGlobalCommunicationPeriod(mdLog_, ir, cr_);
    bGStatEveryStep   = (nstglobalcomm == 1);

    const SimulationGroups* groups = &topGlobal_.groups;

    std::unique_ptr<EssentialDynamics> ed = nullptr;
    if (opt2bSet("-ei", nFile_, fnm_))
    {
        /* Initialize essential dynamics sampling */
        ed = init_edsam(mdLog_,
                        opt2fn_null("-ei", nFile_, fnm_),
                        opt2fn("-eo", nFile_, fnm_),
                        topGlobal_,
                        *ir,
                        cr_,
                        constr_,
                        stateGlobal_,
                        observablesHistory_,
                        oenv_,
                        startingBehavior_);
    }
    else if (observablesHistory_->edsamHistory)
    {
        gmx_fatal(FARGS,
                  "The checkpoint is from a run with essential dynamics sampling, "
                  "but the current run did not specify the -ei option. "
                  "Either specify the -ei option to mdrun, or do not use this checkpoint file.");
    }

    int*                fep_state = MAIN(cr_) ? &stateGlobal_->fep_state : nullptr;
    gmx::ArrayRef<real> lambda    = MAIN(cr_) ? stateGlobal_->lambda : gmx::ArrayRef<real>();
    initialize_lambdas(
            fpLog_, ir->efep, ir->bSimTemp, *ir->fepvals, ir->simtempvals->temperatures, ekind_, MAIN(cr_), fep_state, lambda);
    Update upd(*ir, *ekind_, deform_);

    // Simulated annealing updates the reference temperature.
    const bool doSimulatedAnnealing = initSimulatedAnnealing(*ir, ekind_, &upd);

    const bool useReplicaExchange = (replExParams_.exchangeInterval > 0);

    t_fcdata& fcdata = *fr_->fcdata;

    bool simulationsShareState       = false;
    bool simulationsShareHamiltonian = false;
    int  nstSignalComm               = nstglobalcomm;
    {
        // TODO This implementation of ensemble orientation restraints is nasty because
        // a user can't just do multi-sim with single-sim orientation restraints.
        bool usingEnsembleRestraints =
                (fcdata.disres->nsystems > 1) || ((ms_ != nullptr) && fcdata.orires);
        bool awhUsesMultiSim = (ir->bDoAwh && ir->awhParams->shareBiasMultisim() && (ms_ != nullptr));

        // Replica exchange, ensemble restraints and AWH need all
        // simulations to remain synchronized, so they need
        // checkpoints and stop conditions to act on the same step, so
        // the propagation of such signals must take place between
        // simulations, not just within simulations.
        // TODO: Make algorithm initializers set these flags.
        simulationsShareState = useReplicaExchange || usingEnsembleRestraints || awhUsesMultiSim;

        // With AWH with bias sharing each simulation uses an non-shared, but identical, Hamiltonian
        simulationsShareHamiltonian = useReplicaExchange || usingEnsembleRestraints;

        if (simulationsShareState)
        {
            // Inter-simulation signal communication does not need to happen
            // often, so we use a minimum of 200 steps to reduce overhead.
            const int c_minimumInterSimulationSignallingInterval = 200;
            nstSignalComm = gmx::divideRoundUp(c_minimumInterSimulationSignallingInterval, nstglobalcomm)
                            * nstglobalcomm;
        }
    }

    if (startingBehavior_ != StartingBehavior::RestartWithAppending)
    {
        pleaseCiteCouplingAlgorithms(fpLog_, *ir);
    }
    gmx_mdoutf*       outf = init_mdoutf(fpLog_,
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
                                   startingBehavior_,
                                   simulationsShareState,
                                   ms_);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf),
                                   topGlobal_,
                                   *ir,
                                   pullWork_,
                                   mdoutf_get_fp_dhdl(outf),
                                   false,
                                   startingBehavior_,
                                   simulationsShareHamiltonian,
                                   mdModulesNotifiers_);

    gstat = global_stat_init(ir);

    const auto& simulationWork     = runScheduleWork_->simulationWork;
    const bool  useGpuForPme       = simulationWork.useGpuPme;
    const bool  useGpuForNonbonded = simulationWork.useGpuNonbonded;
    const bool  useGpuForUpdate    = simulationWork.useGpuUpdate;

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fpLog_,
                                 topGlobal_,
                                 constr_ ? constr_->numFlexibleConstraints() : 0,
                                 ir->nstcalcenergy,
                                 haveDDAtomOrdering(*cr_),
                                 useGpuForPme);

    {
        double io = compute_io(ir, topGlobal_.natoms, *groups, energyOutput.numEnergyTerms(), 1);
        if ((io > 2000) && MAIN(cr_))
        {
            fprintf(stderr, "\nWARNING: This run will generate roughly %.0f Mb of data\n\n", io);
        }
    }

    ObservablesReducer observablesReducer = observablesReducerBuilder_->build();

    ForceBuffers     f(simulationWork.useMts,
                   (simulationWork.useGpuFBufferOpsWhenAllowed || useGpuForUpdate)
                               ? PinningPolicy::PinnedIfSupported
                               : PinningPolicy::CannotBePinned);
    const t_mdatoms* md = mdAtoms_->mdatoms();
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
        upd.updateAfterPartition(state_->numAtoms(), md->cFREEZE, md->cTC, md->cACC);
        fr_->longRangeNonbondeds->updateAfterPartition(*md);
    }
    else
    {
        /* Generate and initialize new topology */
        mdAlgorithmsSetupAtomData(
                cr_, *ir, topGlobal_, top_, fr_, &f, mdAtoms_, constr_, virtualSites_, shellfc);

        upd.updateAfterPartition(state_->numAtoms(), md->cFREEZE, md->cTC, md->cACC);
        fr_->longRangeNonbondeds->updateAfterPartition(*md);
    }

    // Now that the state is valid we can set up Parrinello-Rahman
    init_parrinellorahman(ir->pressureCouplingOptions,
                          ir->deform,
                          ir->delta_t * ir->pressureCouplingOptions.nstpcouple,
                          state_->box,
                          state_->box_rel,
                          state_->boxv,
                          &parrinelloRahmanM,
                          &pressureCouplingMu);

    std::unique_ptr<UpdateConstrainGpu> integrator;

    StatePropagatorDataGpu* stateGpu = fr_->stateGpu;

    // TODO: the assertions below should be handled by UpdateConstraintsBuilder.
    if (useGpuForUpdate)
    {
        GMX_RELEASE_ASSERT(!haveDDAtomOrdering(*cr_) || ddUsesUpdateGroups(*cr_->dd)
                                   || constr_ == nullptr || constr_->numConstraintsTotal() == 0,
                           "Constraints in domain decomposition are only supported with update "
                           "groups if using GPU update.\n");
        GMX_RELEASE_ASSERT(ir->eConstrAlg != ConstraintAlgorithm::Shake || constr_ == nullptr
                                   || constr_->numConstraintsTotal() == 0,
                           "SHAKE is not supported with GPU update.");
        GMX_RELEASE_ASSERT(useGpuForPme || (useGpuForNonbonded && simulationWork.useGpuXBufferOpsWhenAllowed),
                           "Either PME or short-ranged non-bonded interaction tasks must run on "
                           "the GPU to use GPU update.\n");
        GMX_RELEASE_ASSERT(ir->eI == IntegrationAlgorithm::MD,
                           "Only the md integrator is supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(
                ir->etc != TemperatureCoupling::NoseHoover,
                "Nose-Hoover temperature coupling is not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(
                ir->pressureCouplingOptions.epc == PressureCoupling::No
                        || ir->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman
                        || ir->pressureCouplingOptions.epc == PressureCoupling::Berendsen
                        || ir->pressureCouplingOptions.epc == PressureCoupling::CRescale,
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
                        || (!haveFepPerturbedMasses(topGlobal_) && !havePerturbedConstraints(topGlobal_)),
                "Free energy perturbation of masses and constraints are not supported with the GPU "
                "update.");

        if (constr_ != nullptr && constr_->numConstraintsTotal() > 0)
        {
            GMX_LOG(mdLog_.info)
                    .asParagraph()
                    .appendText("Updating coordinates and applying constraints on the GPU.");
        }
        else
        {
            GMX_LOG(mdLog_.info).asParagraph().appendText("Updating coordinates on the GPU.");
        }
        GMX_RELEASE_ASSERT(fr_->deviceStreamManager != nullptr,
                           "Device stream manager should be initialized in order to use GPU "
                           "update-constraints.");
        GMX_RELEASE_ASSERT(
                fr_->deviceStreamManager->streamIsValid(gmx::DeviceStreamType::UpdateAndConstraints),
                "Update stream should be initialized in order to use GPU "
                "update-constraints.");
        integrator = std::make_unique<UpdateConstrainGpu>(
                *ir,
                topGlobal_,
                ekind_->numTemperatureCouplingGroups(),
                fr_->deviceStreamManager->context(),
                fr_->deviceStreamManager->stream(gmx::DeviceStreamType::UpdateAndConstraints),
                wallCycleCounters_);

        stateGpu->setXUpdatedOnDeviceEvent(integrator->xUpdatedOnDeviceEvent());

        integrator->setPbc(PbcType::Xyz, state_->box);
    }

    if (useGpuForPme || simulationWork.useGpuXBufferOpsWhenAllowed || useGpuForUpdate)
    {
        changePinningPolicy(&state_->x, PinningPolicy::PinnedIfSupported);
    }
    if (useGpuForUpdate)
    {
        changePinningPolicy(&state_->v, PinningPolicy::PinnedIfSupported);
    }

    // NOTE: The global state is no longer used at this point.
    // But state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdAtoms_->mdatoms(), state_->lambda[FreeEnergyPerturbationCouplingType::Mass]);

    if (ir->bExpanded)
    {
        /* Check nstexpanded here, because the grompp check was broken */
        if (ir->expandedvals->nstexpanded % ir->nstcalcenergy != 0)
        {
            gmx_fatal(FARGS,
                      "With expanded ensemble, nstexpanded should be a multiple of nstcalcenergy");
        }
        init_expanded_ensemble(startingBehavior_ != StartingBehavior::NewSimulation, ir, state_->dfhist);
    }

    if (MAIN(cr_))
    {
        EnergyData::initializeEnergyHistory(startingBehavior_, observablesHistory_, &energyOutput);
    }

    preparePrevStepPullCom(
            ir, pullWork_, md->massT, state_, stateGlobal_, cr_, startingBehavior_ != StartingBehavior::NewSimulation);

    // TODO: Remove this by converting AWH into a ForceProvider
    auto awh = prepareAwhModule(fpLog_,
                                *ir,
                                stateGlobal_,
                                cr_,
                                ms_,
                                startingBehavior_ != StartingBehavior::NewSimulation,
                                shellfc != nullptr,
                                opt2fn("-awh", nFile_, fnm_),
                                pullWork_);

    if (useReplicaExchange && MAIN(cr_))
    {
        repl_ex = init_replica_exchange(fpLog_, ms_, topGlobal_.natoms, ir, replExParams_);
    }
    /* PME tuning is only supported in the Verlet scheme, with PME for
     * Coulomb. It is not supported with only LJ PME.
     * Disable PME tuning with GPU PME decomposition */
    bPMETune = (mdrunOptions_.tunePme && usingPme(fr_->ic->eeltype) && !mdrunOptions_.reproducible
                && ir->cutoff_scheme != CutoffScheme::Group && !simulationWork.useGpuPmeDecomposition);

    pme_load_balancing_t* pme_loadbal = nullptr;
    if (bPMETune)
    {
        pme_loadbal_init(
                &pme_loadbal, cr_, mdLog_, *ir, state_->box, *fr_->ic, *fr_->nbv, fr_->pmedata, fr_->nbv->useGpu());
    }

    if (!ir->bContinuation)
    {
        if (state_->hasEntry(StateEntry::V))
        {
            auto v = makeArrayRef(state_->v);
            /* Set the velocities of vsites, shells and frozen atoms to zero */
            for (i = 0; i < md->homenr; i++)
            {
                if (md->ptype[i] == ParticleType::Shell)
                {
                    clear_rvec(v[i]);
                }
                else if (!md->cFREEZE.empty())
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

        if (constr_)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fpLog_,
                               constr_,
                               *ir,
                               md->nr,
                               md->homenr,
                               state_->x.arrayRefWithPadding(),
                               state_->v.arrayRefWithPadding(),
                               state_->box,
                               state_->lambda[FreeEnergyPerturbationCouplingType::Bonded]);
        }
    }

    const int nstfep = computeFepPeriod(*ir, replExParams_);

    /* Be REALLY careful about what flags you set here. You CANNOT assume
     * this is the first step, since we might be restarting from a checkpoint,
     * and in that case we should not do any modifications to the state.
     */
    bStopCM = (ir->comm_mode != ComRemovalAlgorithm::No && !ir->bContinuation);

    // When restarting from a checkpoint, it can be appropriate to
    // initialize ekind from quantities in the checkpoint. Otherwise,
    // compute_globals must initialize ekind before the simulation
    // starts/restarts. However, only the main rank knows what was
    // found in the checkpoint file, so we have to communicate in
    // order to coordinate the restart.
    //
    // TODO Consider removing this communication if/when checkpoint
    // reading directly follows .tpr reading, because all ranks can
    // agree on hasReadEkinState at that time.
    bool hasReadEkinState = MAIN(cr_) ? stateGlobal_->ekinstate.hasReadEkinState : false;
    if (PAR(cr_))
    {
        gmx_bcast(sizeof(hasReadEkinState), &hasReadEkinState, cr_->mpi_comm_mygroup);
    }
    if (hasReadEkinState)
    {
        restore_ekinstate_from_state(cr_, ekind_, MAIN(cr_) ? &stateGlobal_->ekinstate : nullptr);
    }

    unsigned int cglo_flags =
            (CGLO_TEMPERATURE | CGLO_GSTAT | (EI_VV(ir->eI) ? CGLO_PRESSURE : 0)
             | (EI_VV(ir->eI) ? CGLO_CONSTRAINT : 0) | (hasReadEkinState ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;

    t_vcm vcm(topGlobal_.groups, *ir);
    reportComRemovalInfo(fpLog_, vcm);

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
                        cr_,
                        ir,
                        fr_,
                        ekind_,
                        makeConstArrayRef(state_->x),
                        makeConstArrayRef(state_->v),
                        state_->box,
                        md,
                        nrnb_,
                        &vcm,
                        nullptr,
                        enerd_,
                        force_vir,
                        shake_vir,
                        total_vir,
                        pres,
                        &nullSignaller,
                        state_->box,
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
                             : makeArrayRef(state_->x);
            process_and_stopcm_grp(fpLog_, &vcm, *md, x, makeArrayRef(state_->v));
            inc_nrnb(nrnb_, eNR_STOPCM, md->homenr);
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
                        cr_,
                        ir,
                        fr_,
                        ekind_,
                        makeConstArrayRef(state_->x),
                        makeConstArrayRef(state_->v),
                        state_->box,
                        md,
                        nrnb_,
                        &vcm,
                        nullptr,
                        enerd_,
                        force_vir,
                        shake_vir,
                        total_vir,
                        pres,
                        &nullSignaller,
                        state_->box,
                        &bSumEkinhOld,
                        cglo_flags & ~CGLO_PRESSURE,
                        step,
                        &observablesReducer);
        // Clean up after pre-step use of compute_globals()
        observablesReducer.markAsReadyToReduce();
    }

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (startingBehavior_ == StartingBehavior::NewSimulation)
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind_->tcstat[i].ekinh, ekind_->tcstat[i].ekinh_old);
        }
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants
       for non-trotter temperature control */
    auto trotter_seq = init_npt_vars(ir, *ekind_, state_, &MassQ, bTrotter);

    if (MAIN(cr_))
    {
        if (!ir->bContinuation)
        {
            if (constr_ && ir->eConstrAlg == ConstraintAlgorithm::Lincs)
            {
                fprintf(fpLog_,
                        "RMS relative constraint deviation after constraining: %.2e\n",
                        constr_->rmsd());
            }
            if (EI_STATE_VELOCITY(ir->eI))
            {
                real temp = enerd_->term[F_TEMP];
                if (ir->eI != IntegrationAlgorithm::VV)
                {
                    /* Result of Ekin averaged over velocities of -half
                     * and +half step, while we only have -half step here.
                     */
                    temp *= 2;
                }
                fprintf(fpLog_, "Initial temperature: %g K\n", temp);
            }
        }

        char tbuf[20];
        fprintf(stderr, "starting mdrun '%s'\n", *(topGlobal_.name));
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

    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bInitStep        = startingBehavior_ == StartingBehavior::NewSimulation || EI_VV(ir->eI);
    bSumEkinhOld     = FALSE;
    bExchanged       = FALSE;
    bNeedRepartition = FALSE;

    auto stopHandler = stopHandlerBuilder_->getStopHandlerMD(
            compat::not_null<SimulationSignal*>(&signals[eglsSTOPCOND]),
            simulationsShareState,
            MAIN(cr_),
            ir->nstlist,
            mdrunOptions_.reproducible,
            nstSignalComm,
            mdrunOptions_.maximumHoursToRun,
            ir->nstlist == 0,
            fpLog_,
            step,
            bNS,
            wallTimeAccounting_);

    real checkpointPeriod = mdrunOptions_.checkpointOptions.period;
    if (ir->bExpanded)
    {
        GMX_LOG(mdLog_.info)
                .asParagraph()
                .appendText(
                        "Expanded ensemble with the legacy simulator does not always "
                        "checkpoint correctly, so checkpointing is disabled. You will "
                        "not be able to do a checkpoint restart of this simulation. "
                        "If you use the modular simulator (e.g. by choosing md-vv integrator) "
                        "then checkpointing is enabled. See "
                        "https://gitlab.com/gromacs/gromacs/-/issues/4629 for details.");
        // Use a negative period to disable checkpointing.
        checkpointPeriod = -1;
    }
    auto checkpointHandler = std::make_unique<CheckpointHandler>(
            compat::make_not_null<SimulationSignal*>(&signals[eglsCHKPT]),
            simulationsShareState,
            ir->nstlist == 0,
            MAIN(cr_),
            mdrunOptions_.writeConfout,
            checkpointPeriod);

    const bool resetCountersIsLocal = true;
    auto       resetHandler         = std::make_unique<ResetHandler>(
            compat::make_not_null<SimulationSignal*>(&signals[eglsRESETCOUNTERS]),
            !resetCountersIsLocal,
            ir->nsteps,
            MAIN(cr_),
            mdrunOptions_.timingOptions.resetHalfway,
            mdrunOptions_.maximumHoursToRun,
            mdLog_,
            wallCycleCounters_,
            wallTimeAccounting_);

    const DDBalanceRegionHandler ddBalanceRegionHandler(cr_);

    if (MAIN(cr_) && isMultiSim(ms_) && !useReplicaExchange)
    {
        logInitialMultisimStatus(ms_, cr_, mdLog_, simulationsShareState, ir->nsteps, ir->init_step);
    }

    bool usedMdGpuGraphLastStep = false;
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
                stateGpu->copyCoordinatesFromGpu(state_->x, AtomLocality::Local);
                stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
            }
            /* PME grid + cut-off optimization with GPUs or PME nodes */
            pme_loadbal_do(pme_loadbal,
                           cr_,
                           (mdrunOptions_.verbose && MAIN(cr_)) ? stderr : nullptr,
                           fpLog_,
                           mdLog_,
                           *ir,
                           fr_,
                           state_->box,
                           state_->x,
                           wallCycleCounters_,
                           step,
                           step_rel,
                           &bPMETunePrinting,
                           simulationWork.useGpuPmePpCommunication);
        }

        wallcycle_start(wallCycleCounters_, WallCycleCounter::Step);

        bLastStep = (step_rel == ir->nsteps);
        t         = t0 + step * ir->delta_t;

        // TODO Refactor this, so that nstfep does not need a default value of zero
        if (ir->efep != FreeEnergyPerturbationType::No || ir->bSimTemp)
        {
            /* find and set the current lambdas */
            state_->lambda = currentLambdas(step, *(ir->fepvals), state_->fep_state);

            bDoExpanded = (do_per_step(step, ir->expandedvals->nstexpanded) && (ir->bExpanded)
                           && (!bFirstStep));
        }

        bDoReplEx = (useReplicaExchange && (step > 0) && !bLastStep
                     && do_per_step(step, replExParams_.exchangeInterval));

        if (doSimulatedAnnealing)
        {
            // Simulated annealing updates the reference temperature.
            update_annealing_target_temp(*ir, t, ekind_, &upd);
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
                  || (bFirstStep && startingBehavior_ == StartingBehavior::NewSimulation) || bLastStep);
        do_verbose = mdrunOptions_.verbose
                     && (step % mdrunOptions_.verboseStepPrintInterval == 0 || bFirstStep || bLastStep);

        // On search steps, when doing the update on the GPU, copy
        // the coordinates and velocities to the host unless they are
        // already there (ie on the first step and after replica
        // exchange).
        if (useGpuForUpdate && bNS && !bFirstStep && !bExchanged)
        {
            if (usedMdGpuGraphLastStep)
            {
                // Wait on coordinates produced from GPU graph
                stateGpu->waitCoordinatesUpdatedOnDevice();
            }
            stateGpu->copyVelocitiesFromGpu(state_->v, AtomLocality::Local);
            stateGpu->copyCoordinatesFromGpu(state_->x, AtomLocality::Local);
            stateGpu->waitVelocitiesReadyOnHost(AtomLocality::Local);
            stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
        }

        // We need to calculate virtual velocities if we are writing them in the current step.
        // They also need to be periodically updated. Every 1000 steps is arbitrary, but a reasonable number.
        // The reason why the velocities need to be updated regularly is that the virtual site coordinates
        // are updated using these velocities during integration. Those coordinates are used for, e.g., domain
        // decomposition. Before computing any forces the positions of the virtual sites are recalculated.
        // This fixes a bug, #4879, which was introduced in MR !979.
        const int  c_virtualSiteVelocityUpdateInterval = 1000;
        const bool needVirtualVelocitiesThisStep =
                (virtualSites_ != nullptr)
                && (do_per_step(step, ir->nstvout) || checkpointHandler->isCheckpointingStep()
                    || do_per_step(step, c_virtualSiteVelocityUpdateInterval));

        if (virtualSites_ != nullptr)
        {
            // Virtual sites need to be updated before domain decomposition and forces are calculated
            wallcycle_start(wallCycleCounters_, WallCycleCounter::VsiteConstr);
            // md-vv calculates virtual velocities once it has full-step real velocities
            virtualSites_->construct(state_->x,
                                     state_->v,
                                     state_->box,
                                     (!EI_VV(inputRec_->eI) && needVirtualVelocitiesThisStep)
                                             ? VSiteOperation::PositionsAndVelocities
                                             : VSiteOperation::Positions);
            wallcycle_stop(wallCycleCounters_, WallCycleCounter::VsiteConstr);
        }

        if (bNS && !(bFirstStep && ir->bContinuation))
        {
            bMainState = FALSE;
            /* Correct the new box if it is too skewed */
            if (inputrecDynamicBox(ir))
            {
                if (correct_box(fpLog_, step, state_->box))
                {
                    bMainState = TRUE;
                }
            }
            // If update is offloaded, and the box was changed either
            // above or in a replica exchange on the previous step,
            // the GPU Update object should be informed
            if (useGpuForUpdate && (bMainState || bExchanged))
            {
                integrator->setPbc(PbcType::Xyz, state_->box);
            }
            if (haveDDAtomOrdering(*cr_) && bMainState)
            {
                dd_collect_state(cr_->dd, state_, stateGlobal_);
            }

            if (haveDDAtomOrdering(*cr_))
            {
                /* Repartition the domain decomposition */
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
                                    do_verbose && !bPMETunePrinting);
                upd.updateAfterPartition(state_->numAtoms(), md->cFREEZE, md->cTC, md->cACC);
                fr_->longRangeNonbondeds->updateAfterPartition(*md);
            }
        }

        // Allocate or re-size GPU halo exchange object, if necessary
        if (bNS && simulationWork.havePpDomainDecomposition && simulationWork.useGpuHaloExchange)
        {
            GMX_RELEASE_ASSERT(fr_->deviceStreamManager != nullptr,
                               "GPU device manager has to be initialized to use GPU "
                               "version of halo exchange.");
            constructGpuHaloExchange(*cr_, *fr_->deviceStreamManager, wallCycleCounters_);
        }

        if (MAIN(cr_) && do_log)
        {
            gmx::EnergyOutput::printHeader(
                    fpLog_, step, t); /* can we improve the information printed here? */
        }

        if (ir->efep != FreeEnergyPerturbationType::No)
        {
            update_mdatoms(mdAtoms_->mdatoms(), state_->lambda[FreeEnergyPerturbationCouplingType::Mass]);
        }

        if (bExchanged)
        {
            /* We need the kinetic energy at minus the half step for determining
             * the full step kinetic energy and possibly for T-coupling.*/
            /* This may not be quite working correctly yet . . . . */
            int cglo_flags = CGLO_GSTAT | CGLO_TEMPERATURE;
            compute_globals(gstat,
                            cr_,
                            ir,
                            fr_,
                            ekind_,
                            makeConstArrayRef(state_->x),
                            makeConstArrayRef(state_->v),
                            state_->box,
                            md,
                            nrnb_,
                            &vcm,
                            wallCycleCounters_,
                            enerd_,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            &nullSignaller,
                            state_->box,
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
                       || (ir->pressureCouplingOptions.epc != PressureCoupling::No
                           && (do_per_step(step, ir->pressureCouplingOptions.nstpcouple)
                               || do_per_step(step - 1, ir->pressureCouplingOptions.nstpcouple)));
        }
        else
        {
            bCalcEnerStep = do_per_step(step, ir->nstcalcenergy);
            bCalcVir      = bCalcEnerStep
                       || (ir->pressureCouplingOptions.epc != PressureCoupling::No
                           && do_per_step(step, ir->pressureCouplingOptions.nstpcouple));
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

        if (bNS)
        {
            if (fr_->listedForcesGpu)
            {
                fr_->listedForcesGpu->updateHaveInteractions(top_->idef);
            }
            runScheduleWork_->domainWork = setupDomainLifetimeWorkload(
                    *ir, *fr_, pullWork_, ed ? ed->getLegacyED() : nullptr, *md, simulationWork);
        }

        const int shellfcFlags = force_flags | (mdrunOptions_.verbose ? GMX_FORCE_ENERGY : 0);
        const int legacyForceFlags = ((shellfc) ? shellfcFlags : force_flags) | (bNS ? GMX_FORCE_NS : 0);

        runScheduleWork_->stepWork = setupStepWorkload(
                legacyForceFlags, ir->mtsLevels, step, runScheduleWork_->domainWork, simulationWork);

        const bool doTemperatureScaling = (ir->etc != TemperatureCoupling::No
                                           && do_per_step(step + ir->nsttcouple - 1, ir->nsttcouple));

        /* With leap-frog type integrators we compute the kinetic energy
         * at a whole time step as the average of the half-time step kinetic
         * energies of two subsequent steps. Therefore we need to compute the
         * half step kinetic energy also if we need energies at the next step.
         */
        const bool needHalfStepKineticEnergy =
                (!EI_VV(ir->eI) && (do_per_step(step + 1, nstglobalcomm) || step_rel + 1 == ir->nsteps));

        // Parrinello-Rahman requires the pressure to be availible before the update to compute
        // the velocity scaling matrix. Hence, it runs one step after the nstpcouple step.
        const bool doParrinelloRahman =
                (ir->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman
                 && do_per_step(step + ir->pressureCouplingOptions.nstpcouple - 1,
                                ir->pressureCouplingOptions.nstpcouple));

        MdGpuGraph* mdGraph = simulationWork.useMdGpuGraph ? fr_->mdGraph[step % 2].get() : nullptr;

        if (simulationWork.useMdGpuGraph)
        {
            // Reset graph on search step (due to changing neighbour list etc)
            // or virial step (due to changing shifts and box).
            if (bNS || bCalcVir)
            {
                fr_->mdGraph[MdGraphEvenOrOddStep::EvenStep]->reset();
                fr_->mdGraph[MdGraphEvenOrOddStep::OddStep]->reset();
            }
            else
            {
                mdGraph->setUsedGraphLastStep(usedMdGpuGraphLastStep);
                bool canUseMdGpuGraphThisStep =
                        !bNS && !bCalcVir && !doTemperatureScaling && !doParrinelloRahman && !bGStat
                        && !needHalfStepKineticEnergy && !do_per_step(step, ir->nstxout)
                        && !do_per_step(step, ir->nstxout_compressed)
                        && !do_per_step(step, ir->nstvout) && !do_per_step(step, ir->nstfout)
                        && !checkpointHandler->isCheckpointingStep();
                if (mdGraph->captureThisStep(canUseMdGpuGraphThisStep))
                {
                    mdGraph->startRecord(stateGpu->getCoordinatesReadyOnDeviceEvent(
                            AtomLocality::Local, simulationWork, runScheduleWork_->stepWork));
                }
            }
        }
        if (!simulationWork.useMdGpuGraph || mdGraph->graphIsCapturingThisStep()
            || !mdGraph->useGraphThisStep())
        {

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
                                    *md,
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
                /* The AWH history need to be saved _before_ doing force calculations where the AWH bias
                   is updated (or the AWH update will be performed twice for one step when continuing).
                   It would be best to call this update function from do_md_trajectory_writing but that
                   would occur after do_force. One would have to divide the update_awh function into one
                   function applying the AWH force and one doing the AWH bias update. The update AWH
                   bias function could then be called after do_md_trajectory_writing (then containing
                   update_awh_history). The checkpointing will in the future probably moved to the start
                   of the md loop which will rid of this issue. */
                if (awh && checkpointHandler->isCheckpointingStep() && MAIN(cr_))
                {
                    awh->updateHistory(stateGlobal_->awhHistory.get());
                }

                /* The coordinates (x) are shifted (to get whole molecules)
                 * in do_force.
                 * This is parallellized as well, and does communication too.
                 * Check comments in sim_util.c
                 */
                do_force(fpLog_,
                         cr_,
                         ms_,
                         *ir,
                         mdModulesNotifiers_,
                         awh.get(),
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
                         md,
                         enerd_,
                         state_->lambda,
                         fr_,
                         *runScheduleWork_,
                         virtualSites_,
                         mu_tot,
                         t,
                         ed ? ed->getLegacyED() : nullptr,
                         fr_->longRangeNonbondeds.get(),
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
                                     startingBehavior_,
                                     nstglobalcomm,
                                     ir,
                                     fr_,
                                     cr_,
                                     state_,
                                     mdAtoms_->mdatoms(),
                                     &fcdata,
                                     &MassQ,
                                     &vcm,
                                     enerd_,
                                     &observablesReducer,
                                     ekind_,
                                     gstat,
                                     &last_ekin,
                                     bCalcVir,
                                     total_vir,
                                     shake_vir,
                                     force_vir,
                                     pres,
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
                                     constr_,
                                     &nullSignaller,
                                     trotter_seq,
                                     nrnb_,
                                     fpLog_,
                                     wallCycleCounters_);
                if (virtualSites_ != nullptr && needVirtualVelocitiesThisStep)
                {
                    // Positions were calculated earlier
                    wallcycle_start(wallCycleCounters_, WallCycleCounter::VsiteConstr);
                    virtualSites_->construct(state_->x, state_->v, state_->box, VSiteOperation::Velocities);
                    wallcycle_stop(wallCycleCounters_, WallCycleCounter::VsiteConstr);
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
                lamnew = ExpandedEnsembleDynamics(fpLog_,
                                                  *inputRec_,
                                                  *enerd_,
                                                  ekind_,
                                                  state_,
                                                  &MassQ,
                                                  state_->fep_state,
                                                  state_->dfhist,
                                                  step,
                                                  state_->v.rvec_array(),
                                                  md->homenr,
                                                  md->cTC);
                /* history is maintained in state->dfhist, but state_global is what is sent to trajectory and log output */
                if (MAIN(cr_))
                {
                    copy_df_history(stateGlobal_->dfhist, state_->dfhist);
                }
            }

            // Copy coordinate from the GPU for the output/checkpointing if the update is offloaded
            // and coordinates have not already been copied for i) search or ii) CPU force tasks.
            if (useGpuForUpdate && !bNS && !runScheduleWork_->domainWork.haveCpuLocalForceWork
                && (do_per_step(step, ir->nstxout) || do_per_step(step, ir->nstxout_compressed)
                    || checkpointHandler->isCheckpointingStep()))
            {
                stateGpu->copyCoordinatesFromGpu(state_->x, AtomLocality::Local);
                stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
            }
            // Copy velocities if needed for the output/checkpointing.
            // NOTE: Copy on the search steps is done at the beginning of the step.
            if (useGpuForUpdate && !bNS
                && (do_per_step(step, ir->nstvout) || checkpointHandler->isCheckpointingStep()))
            {
                stateGpu->copyVelocitiesFromGpu(state_->v, AtomLocality::Local);
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
            if (runScheduleWork_->stepWork.useGpuFBufferOps
                && (simulationWork.useGpuUpdate && !virtualSites_) && do_per_step(step, ir->nstfout))
            {
                stateGpu->copyForcesFromGpu(f.view().force(), AtomLocality::Local);
                stateGpu->waitForcesReadyOnHost(AtomLocality::Local);
            }
            /* Now we have the energies and forces corresponding to the
             * coordinates at time t. We must output all of this before
             * the update.
             */
            const EkindataState ekindataState =
                    bGStat ? (bSumEkinhOld ? EkindataState::UsedNeedToReduce
                                           : EkindataState::UsedDoNotNeedToReduce)
                           : EkindataState::NotUsed;
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
                                     checkpointHandler->isCheckpointingStep(),
                                     bRerunMD,
                                     bLastStep,
                                     mdrunOptions_.writeConfout,
                                     ekindataState);
            /* Check if IMD step and do IMD communication, if bIMD is TRUE. */
            bInteractiveMDstep = imdSession_->run(step, bNS, state_->box, state_->x, t);

            /* kludge -- virial is lost with restart for MTTK NPT control. Must reload (saved earlier). */
            if (startingBehavior_ != StartingBehavior::NewSimulation && bFirstStep
                && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir)))
            {
                copy_mat(state_->svir_prev, shake_vir);
                copy_mat(state_->fvir_prev, force_vir);
            }

            stopHandler->setSignal();
            resetHandler->setSignal(wallTimeAccounting_);

            if (bGStat || !PAR(cr_))
            {
                /* In parallel we only have to check for checkpointing in steps
                 * where we do global communication,
                 *  otherwise the other nodes don't know.
                 */
                checkpointHandler->setSignal(wallTimeAccounting_);
            }

            /* #########   START SECOND UPDATE STEP ################# */

            /* at the start of step, randomize or scale the velocities ((if vv. Restriction of
               Andersen controlled in preprocessing */

            if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
            {
                gmx_bool bIfRandomize;
                bIfRandomize = update_randomize_velocities(
                        ir, step, cr_, md->homenr, md->cTC, md->invmass, state_->v, &upd, constr_);
                /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
                if (constr_ && bIfRandomize)
                {
                    constrain_velocities(constr_, do_log || do_ene, step, state_, nullptr, false, nullptr);
                }
            }
            /* Box is changed in update() when we do pressure coupling,
             * but we should still use the old box for energy corrections and when
             * writing it to the energy file, so it matches the trajectory files for
             * the same timestep above. Make a copy in a separate array.
             */
            copy_mat(state_->box, lastbox);

            dvdl_constr = 0;

            if (!useGpuForUpdate)
            {
                wallcycle_start(wallCycleCounters_, WallCycleCounter::Update);
            }
            /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
            if (bTrotter)
            {
                trotter_update(ir,
                               step,
                               ekind_,
                               state_,
                               total_vir,
                               md->homenr,
                               md->cTC,
                               md->invmass,
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
                update_tcouple(step, ir, state_, ekind_, &MassQ, md->homenr, md->cTC);
                update_pcouple_before_coordinates(mdLog_,
                                                  step,
                                                  ir->pressureCouplingOptions,
                                                  ir->deform,
                                                  ir->delta_t,
                                                  state_,
                                                  &pressureCouplingMu,
                                                  &parrinelloRahmanM);
            }

            if (EI_VV(ir->eI))
            {
                GMX_ASSERT(!useGpuForUpdate, "GPU update is not supported with VVAK integrator.");

                integrateVVSecondStep(step,
                                      ir,
                                      fr_,
                                      cr_,
                                      state_,
                                      mdAtoms_->mdatoms(),
                                      &fcdata,
                                      &MassQ,
                                      &vcm,
                                      pullWork_,
                                      enerd_,
                                      &observablesReducer,
                                      ekind_,
                                      gstat,
                                      &dvdl_constr,
                                      bCalcVir,
                                      total_vir,
                                      shake_vir,
                                      force_vir,
                                      pres,
                                      lastbox,
                                      do_log,
                                      do_ene,
                                      bGStat,
                                      &bSumEkinhOld,
                                      &f,
                                      &cbuf,
                                      &upd,
                                      constr_,
                                      &nullSignaller,
                                      trotter_seq,
                                      nrnb_,
                                      wallCycleCounters_);
            }
            else
            {
                if (useGpuForUpdate)
                {
                    // On search steps, update handles to device vectors
                    // TODO: this condition has redundant / unnecessary clauses
                    if (bNS && (bFirstStep || haveDDAtomOrdering(*cr_) || bExchanged))
                    {
                        integrator->set(stateGpu->getCoordinates(),
                                        stateGpu->getVelocities(),
                                        stateGpu->getForces(),
                                        top_->idef,
                                        *md);

                        // Copy data to the GPU after buffers might have been reinitialized
                        /* The velocity copy is redundant if we had Center-of-Mass motion removed on
                         * the previous step. We don't check that now. */
                        stateGpu->copyVelocitiesToGpu(state_->v, AtomLocality::Local);
                    }

                    // Copy x to the GPU unless we have already transferred in do_force().
                    // We transfer in do_force() if a GPU force task requires x (PME or x buffer ops).
                    if (!(runScheduleWork_->stepWork.haveGpuPmeOnThisRank
                          || runScheduleWork_->stepWork.useGpuXBufferOps))
                    {
                        stateGpu->copyCoordinatesToGpu(state_->x, AtomLocality::Local);
                        // Coordinates are later used by the integrator running in the same stream.
                        stateGpu->consumeCoordinatesCopiedToDeviceEvent(AtomLocality::Local);
                    }

                    if ((simulationWork.useGpuPme && simulationWork.useCpuPmePpCommunication)
                        || (!runScheduleWork_->stepWork.useGpuFBufferOps))
                    {
                        // The PME forces were recieved to the host, and reduced on the CPU with the
                        // rest of the forces computed on the GPU, so the final forces have to be
                        // copied back to the GPU. Or the buffer ops were not offloaded this step,
                        // so the forces are on the host and have to be copied
                        stateGpu->copyForcesToGpu(f.view().force(), AtomLocality::Local);
                    }
                    const bool doTemperatureScaling =
                            (ir->etc != TemperatureCoupling::No
                             && do_per_step(step + ir->nsttcouple - 1, ir->nsttcouple));

                    // This applies Leap-Frog, LINCS and SETTLE in succession
                    integrator->integrate(
                            stateGpu->getLocalForcesReadyOnDeviceEvent(
                                    runScheduleWork_->stepWork, runScheduleWork_->simulationWork),
                            ir->delta_t,
                            true,
                            bCalcVir,
                            shake_vir,
                            doTemperatureScaling,
                            ekind_->tcstat,
                            doParrinelloRahman,
                            ir->pressureCouplingOptions.nstpcouple * ir->delta_t,
                            parrinelloRahmanM);
                }
                else
                {
                    /* With multiple time stepping we need to do an additional normal
                     * update step to obtain the virial and dH/dl, as the actual MTS integration
                     * using an acceleration where the slow forces are multiplied by mtsFactor.
                     * Using that acceleration would result in a virial with the slow
                     * force contribution would be a factor mtsFactor too large.
                     */
                    const bool separateVirialConstraining =
                            (simulationWork.useMts && (bCalcVir || computeDHDL) && constr_ != nullptr);
                    if (separateVirialConstraining)
                    {
                        upd.update_for_constraint_virial(*ir,
                                                         md->homenr,
                                                         md->havePartiallyFrozenAtoms,
                                                         md->invmass,
                                                         md->invMassPerDim,
                                                         *state_,
                                                         f.view().forceWithPadding(),
                                                         *ekind_);

                        // Call apply() directly so we can avoid constraining the velocities
                        constr_->apply(false,
                                       step,
                                       1,
                                       1.0,
                                       state_->x.arrayRefWithPadding(),
                                       upd.xp()->arrayRefWithPadding(),
                                       {},
                                       state_->box,
                                       state_->lambda[FreeEnergyPerturbationCouplingType::Bonded],
                                       &dvdl_constr,
                                       {},
                                       bCalcVir,
                                       shake_vir,
                                       ConstraintVariable::Positions);
                    }

                    ArrayRefWithPadding<const RVec> forceCombined =
                            (simulationWork.useMts && step % ir->mtsLevels[1].stepFactor == 0)
                                    ? f.view().forceMtsCombinedWithPadding()
                                    : f.view().forceWithPadding();
                    upd.update_coords(*ir,
                                      step,
                                      md->homenr,
                                      md->havePartiallyFrozenAtoms,
                                      md->ptype,
                                      md->invmass,
                                      md->invMassPerDim,
                                      state_,
                                      forceCombined,
                                      &fcdata,
                                      ekind_,
                                      parrinelloRahmanM,
                                      etrtPOSITION,
                                      cr_,
                                      constr_ != nullptr);

                    wallcycle_stop(wallCycleCounters_, WallCycleCounter::Update);

                    constrain_coordinates(constr_,
                                          do_log || do_ene,
                                          step,
                                          state_,
                                          upd.xp()->arrayRefWithPadding(),
                                          separateVirialConstraining ? nullptr : &dvdl_constr,
                                          bCalcVir && !separateVirialConstraining,
                                          shake_vir);

                    upd.update_sd_second_half(*ir,
                                              step,
                                              &dvdl_constr,
                                              md->homenr,
                                              md->ptype,
                                              md->invmass,
                                              state_,
                                              cr_,
                                              nrnb_,
                                              wallCycleCounters_,
                                              constr_,
                                              do_log,
                                              do_ene);
                    upd.finish_update(*ir,
                                      md->havePartiallyFrozenAtoms,
                                      md->homenr,
                                      state_,
                                      wallCycleCounters_,
                                      constr_ != nullptr);
                }

                if (ir->bPull && ir->pull->bSetPbcRefToPrevStepCOM)
                {
                    updatePrevStepPullCom(pullWork_, state_->pull_com_prev_step);
                }

                enerd_->term[F_DVDL_CONSTR] += dvdl_constr;
            }
        }

        if (simulationWork.useMdGpuGraph)
        {
            GMX_ASSERT((mdGraph != nullptr), "MD GPU graph does not exist.");
            if (mdGraph->graphIsCapturingThisStep())
            {
                mdGraph->endRecord();
                // Force graph reinstantiation (instead of graph exec
                // update): with PME tuning, since the GPU kernels
                // chosen by the FFT library can vary with grid size;
                // or with an odd nstlist, since the odd/even step
                // pruning pattern will change
                bool forceGraphReinstantiation =
                        pme_loadbal_is_active(pme_loadbal) || ((ir->nstlist % 2) == 1);
                mdGraph->createExecutableGraph(forceGraphReinstantiation);
            }
            if (mdGraph->useGraphThisStep())
            {
                mdGraph->launchGraphMdStep(integrator->xUpdatedOnDeviceEvent());
            }
            if (bNS)
            {
                // TODO: merge disableForDomainIfAnyPpRankHasCpuForces() back into reset() when
                // domainWork initialization is moved out of do_force().
                fr_->mdGraph[MdGraphEvenOrOddStep::EvenStep]->disableForDomainIfAnyPpRankHasCpuForces(
                        runScheduleWork_->domainWork.haveCpuLocalForceWork);
                fr_->mdGraph[MdGraphEvenOrOddStep::OddStep]->disableForDomainIfAnyPpRankHasCpuForces(
                        runScheduleWork_->domainWork.haveCpuLocalForceWork);
            }
            usedMdGpuGraphLastStep = mdGraph->useGraphThisStep();
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
                    stateGpu->copyCoordinatesFromGpu(state_->x, AtomLocality::Local);
                    stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
                }

                // Copy velocities back to the host if:
                // - Globals are computed this step (includes the energy output steps).
                // - Temperature is needed for the next step.
                // - This is a replica exchange step (even though we will only need
                //     the velocities if an exchange succeeds)
                if (bGStat || needHalfStepKineticEnergy || bDoReplEx)
                {
                    stateGpu->copyVelocitiesFromGpu(state_->v, AtomLocality::Local);
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
                SimulationSignaller signaller(&signals, cr_, ms_, doInterSimSignal, doIntraSimSignal);

                compute_globals(gstat,
                                cr_,
                                ir,
                                fr_,
                                ekind_,
                                makeConstArrayRef(state_->x),
                                makeConstArrayRef(state_->v),
                                state_->box,
                                md,
                                nrnb_,
                                &vcm,
                                wallCycleCounters_,
                                enerd_,
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
                            fpLog_, &vcm, *md, makeArrayRef(state_->x), makeArrayRef(state_->v));
                    inc_nrnb(nrnb_, eNR_STOPCM, md->homenr);

                    // TODO: The special case of removing CM motion should be dealt more gracefully
                    if (useGpuForUpdate)
                    {
                        // Issue #3988, #4106.
                        stateGpu->resetCoordinatesCopiedToDeviceEvent(AtomLocality::Local);
                        stateGpu->copyCoordinatesToGpu(state_->x, AtomLocality::Local);
                        // Here we block until the H2D copy completes because event sync with the
                        // force kernels that use the coordinates on the next steps is not implemented
                        // (not because of a race on state->x being modified on the CPU while H2D is in progress).
                        stateGpu->waitCoordinatesCopiedToDevice(AtomLocality::Local);
                        // If the COM removal changed the velocities on the CPU, this has to be accounted for.
                        if (vcm.mode != ComRemovalAlgorithm::No)
                        {
                            stateGpu->copyVelocitiesToGpu(state_->v, AtomLocality::Local);
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
            accumulateKineticLambdaComponents(enerd_, state_->lambda, *ir->fepvals);
        }

        const real currentSystemRefT =
                (haveEnsembleTemperature(*ir) ? ekind_->currentEnsembleTemperature() : 0.0_real);
        const bool scaleCoordinates = !useGpuForUpdate || bDoReplEx;
        update_pcouple_after_coordinates(fpLog_,
                                         step,
                                         ir->pressureCouplingOptions,
                                         ir->ld_seed,
                                         currentSystemRefT,
                                         ir->opts.nFreeze,
                                         ir->deform,
                                         ir->delta_t,
                                         md->homenr,
                                         md->cFREEZE,
                                         pres,
                                         force_vir,
                                         shake_vir,
                                         &pressureCouplingMu,
                                         state_,
                                         nrnb_,
                                         upd.deform(),
                                         scaleCoordinates);

        const bool doBerendsenPressureCoupling =
                (inputRec_->pressureCouplingOptions.epc == PressureCoupling::Berendsen
                 && do_per_step(step, inputRec_->pressureCouplingOptions.nstpcouple));
        const bool doCRescalePressureCoupling =
                (inputRec_->pressureCouplingOptions.epc == PressureCoupling::CRescale
                 && do_per_step(step, inputRec_->pressureCouplingOptions.nstpcouple));
        if (useGpuForUpdate
            && (doBerendsenPressureCoupling || doCRescalePressureCoupling || doParrinelloRahman))
        {
            integrator->scaleCoordinates(pressureCouplingMu);
            if (doCRescalePressureCoupling)
            {
                integrator->scaleVelocities(invertBoxMatrix(pressureCouplingMu));
            }
            integrator->setPbc(PbcType::Xyz, state_->box);
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
                enerd_->term[F_EKIN] = last_ekin;
            }
            enerd_->term[F_ETOT] = enerd_->term[F_EPOT] + enerd_->term[F_EKIN];

            if (integratorHasConservedEnergyQuantity(ir))
            {
                if (EI_VV(ir->eI))
                {
                    enerd_->term[F_ECONSERVED] = enerd_->term[F_ETOT] + saved_conserved_quantity;
                }
                else
                {
                    enerd_->term[F_ECONSERVED] =
                            enerd_->term[F_ETOT]
                            + NPT_energy(ir->pressureCouplingOptions,
                                         ir->etc,
                                         gmx::constArrayRefFromArray(ir->opts.nrdf, ir->opts.ngtc),
                                         *ekind_,
                                         inputrecNvtTrotter(ir) || inputrecNptTrotter(ir),
                                         state_,
                                         &MassQ);
                }
            }
            /* #########  END PREPARING EDR OUTPUT  ###########  */
        }

        /* Output stuff */
        if (MAIN(cr_))
        {
            if (fpLog_ && do_log && bDoExpanded)
            {
                /* only needed if doing expanded ensemble */
                PrintFreeEnergyInfoToFile(fpLog_,
                                          ir->fepvals.get(),
                                          ir->expandedvals.get(),
                                          ir->bSimTemp ? ir->simtempvals.get() : nullptr,
                                          stateGlobal_->dfhist,
                                          state_->fep_state,
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
                                                 enerd_,
                                                 ir->fepvals.get(),
                                                 lastbox,
                                                 PTCouplingArrays{ state_->boxv,
                                                                   state_->nosehoover_xi,
                                                                   state_->nosehoover_vxi,
                                                                   state_->nhpres_xi,
                                                                   state_->nhpres_vxi },
                                                 state_->fep_state,
                                                 total_vir,
                                                 pres,
                                                 ekind_,
                                                 mu_tot,
                                                 constr_);
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
                        do_log ? fpLog_ : nullptr, *groups, ir->opts, *ekind_);
            }
            if (do_log || do_ene || do_dr || do_or)
            {
                energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                                   do_ene,
                                                   do_dr,
                                                   do_or,
                                                   do_log ? fpLog_ : nullptr,
                                                   step,
                                                   t,
                                                   fr_->fcdata.get(),
                                                   awh.get());
            }
            if (do_log && ((ir->bDoAwh && awh->hasFepLambdaDimension()) || ir->fepvals->delta_lambda != 0))
            {
                const bool isInitialOutput = false;
                printLambdaStateToLog(fpLog_, state_->lambda, isInitialOutput);
            }

            if (ir->bPull)
            {
                pull_print_output(pullWork_, step, t);
            }

            if (do_per_step(step, ir->nstlog))
            {
                if (fflush(fpLog_) != 0)
                {
                    gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }
        if (bDoExpanded)
        {
            /* Have to do this part _after_ outputting the logfile and the edr file */
            /* Gets written into the state at the beginning of next loop*/
            state_->fep_state = lamnew;
        }
        else if (ir->bDoAwh && awh->needForeignEnergyDifferences(step))
        {
            state_->fep_state = awh->fepLambdaState();
        }
        /* Print the remaining wall clock time for the run */
        if (isMainSimMainRank(ms_, MAIN(cr_)) && (do_verbose || gmx_got_usr_signal()) && !bPMETunePrinting)
        {
            if (shellfc)
            {
                fprintf(stderr, "\n");
            }
            print_time(stderr, wallTimeAccounting_, step, ir, cr_);
        }

        /* Ion/water position swapping.
         * Not done in last step since trajectory writing happens before this call
         * in the MD loop and exchanges would be lost anyway. */
        bNeedRepartition = FALSE;
        if ((ir->eSwapCoords != SwapType::No) && (step > 0) && !bLastStep
            && do_per_step(step, ir->swap->nstswap))
        {
            bNeedRepartition = do_swapcoords(cr_,
                                             step,
                                             t,
                                             ir,
                                             swap_,
                                             wallCycleCounters_,
                                             as_rvec_array(state_->x.data()),
                                             state_->box,
                                             MAIN(cr_) && mdrunOptions_.verbose,
                                             bRerunMD);

            if (bNeedRepartition && haveDDAtomOrdering(*cr_))
            {
                dd_collect_state(cr_->dd, state_, stateGlobal_);
            }
        }

        /* Replica exchange */
        bExchanged = FALSE;
        if (bDoReplEx)
        {
            bExchanged =
                    replica_exchange(fpLog_, cr_, ms_, repl_ex, stateGlobal_, enerd_, state_, step, t);
        }

        if ((bExchanged || bNeedRepartition) && haveDDAtomOrdering(*cr_))
        {
            dd_partition_system(fpLog_,
                                mdLog_,
                                step,
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
                                wallCycleCounters_,
                                FALSE);
            upd.updateAfterPartition(state_->numAtoms(), md->cFREEZE, md->cTC, md->cACC);
            fr_->longRangeNonbondeds->updateAfterPartition(*md);
        }

        bFirstStep = FALSE;
        bInitStep  = FALSE;

        /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
        /* With all integrators, except VV, we need to retain the pressure
         * at the current step for coupling at the next step.
         */
        if (state_->hasEntry(StateEntry::PressurePrevious)
            && (bGStatEveryStep
                || (ir->pressureCouplingOptions.nstpcouple > 0
                    && step % ir->pressureCouplingOptions.nstpcouple == 0)))
        {
            /* Store the pressure in t_state for pressure coupling
             * at the next MD step.
             */
            copy_mat(pres, state_->pres_prev);
        }

        /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

        if ((membed_ != nullptr) && (!bLastStep))
        {
            rescale_membed(step_rel, membed_, as_rvec_array(stateGlobal_->x.data()));
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

#if GMX_FAHCORE
        if (MAIN(cr))
        {
            fcReportProgress(ir->nsteps + ir->init_step, step);
        }
#endif

        resetHandler->resetCounters(
                step, step_rel, mdLog_, fpLog_, cr_, fr_->nbv.get(), nrnb_, fr_->pmedata, pme_loadbal, wallCycleCounters_, wallTimeAccounting_);

        /* If bIMD is TRUE, the main updates the IMD energy record and sends positions to VMD client */
        imdSession_->updateEnergyRecordAndSendPositionsAndEnergies(bInteractiveMDstep, step, bCalcEner);

        // any run that uses GPUs must be at least offloading nonbondeds
        const bool usingGpu = simulationWork.useGpuNonbonded;
        if (usingGpu)
        {
            // ensure that GPU errors do not propagate between MD steps
            checkPendingDeviceErrorBetweenSteps();
        }
    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end_time(wallTimeAccounting_);

    if (simulationWork.haveSeparatePmeRank)
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr_);
    }

    if (MAIN(cr_))
    {
        if (ir->nstcalcenergy > 0)
        {
            energyOutput.printEnergyConservation(fpLog_, ir->simulation_part, EI_MD(ir->eI));

            gmx::EnergyOutput::printAnnealingTemperatures(fpLog_, *groups, ir->opts, *ekind_);
            energyOutput.printAverages(fpLog_, groups);
        }
    }
    done_mdoutf(outf);

    if (bPMETune)
    {
        pme_loadbal_done(pme_loadbal, fpLog_, mdLog_, fr_->nbv->useGpu());
    }

    done_shellfc(fpLog_, shellfc, step_rel);

    if (useReplicaExchange && MAIN(cr_))
    {
        print_replica_exchange_statistics(fpLog_, repl_ex);
    }

    walltime_accounting_set_nsteps_done(wallTimeAccounting_, step_rel);

    global_stat_destroy(gstat);
}
