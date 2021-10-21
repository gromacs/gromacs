/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
 * \brief Defines the microstate for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "energydata.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/pullhistory.h"
#include "gromacs/topology/topology.h"

#include "freeenergyperturbationdata.h"
#include "modularsimulator.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"

struct pull_t;
class t_state;

namespace gmx
{
class Awh;

EnergyData::EnergyData(StatePropagatorData*        statePropagatorData,
                       FreeEnergyPerturbationData* freeEnergyPerturbationData,
                       const gmx_mtop_t&           globalTopology,
                       const t_inputrec*           inputrec,
                       const MDAtoms*              mdAtoms,
                       gmx_enerdata_t*             enerd,
                       gmx_ekindata_t*             ekind,
                       const Constraints*          constr,
                       FILE*                       fplog,
                       t_fcdata*                   fcd,
                       const MDModulesNotifiers&   mdModulesNotifiers,
                       bool                        isMasterRank,
                       ObservablesHistory*         observablesHistory,
                       StartingBehavior            startingBehavior,
                       bool                        simulationsShareState,
                       pull_t*                     pullWork) :
    element_(std::make_unique<Element>(this, isMasterRank)),
    isMasterRank_(isMasterRank),
    forceVirialStep_(-1),
    shakeVirialStep_(-1),
    totalVirialStep_(-1),
    pressureStep_(-1),
    needToSumEkinhOld_(false),
    hasReadEkinFromCheckpoint_(false),
    startingBehavior_(startingBehavior),
    statePropagatorData_(statePropagatorData),
    freeEnergyPerturbationData_(freeEnergyPerturbationData),
    inputrec_(inputrec),
    top_global_(globalTopology),
    mdAtoms_(mdAtoms),
    enerd_(enerd),
    ekind_(ekind),
    constr_(constr),
    fplog_(fplog),
    fcd_(fcd),
    mdModulesNotifiers_(mdModulesNotifiers),
    groups_(&globalTopology.groups),
    observablesHistory_(observablesHistory),
    simulationsShareState_(simulationsShareState),
    pullWork_(pullWork)
{
    clear_mat(forceVirial_);
    clear_mat(shakeVirial_);
    clear_mat(totalVirial_);
    clear_mat(pressure_);
    clear_rvec(muTot_);

    init_ekinstate(&ekinstate_, inputrec_);
    observablesHistory_->energyHistory = std::make_unique<energyhistory_t>();
}

void EnergyData::Element::scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction)
{
    if (!isMasterRank_)
    {
        return;
    }
    auto writeEnergy                 = energyWritingStep_ == step;
    auto isEnergyCalculationStep     = energyCalculationStep_ == step;
    auto isFreeEnergyCalculationStep = freeEnergyCalculationStep_ == step;
    if (isEnergyCalculationStep || writeEnergy)
    {
        registerRunFunction([this, step, time, isEnergyCalculationStep, isFreeEnergyCalculationStep]() {
            energyData_->doStep(step, time, isEnergyCalculationStep, isFreeEnergyCalculationStep);
        });
    }
    else
    {
        registerRunFunction([this]() { energyData_->energyOutput_->recordNonEnergyStep(); });
    }
}

void EnergyData::teardown()
{
    if (inputrec_->nstcalcenergy > 0 && isMasterRank_)
    {
        energyOutput_->printEnergyConservation(fplog_, inputrec_->simulation_part, EI_MD(inputrec_->eI));
        energyOutput_->printAverages(fplog_, groups_);
    }
}

void EnergyData::Element::trajectoryWriterSetup(gmx_mdoutf* outf)
{
    energyData_->setup(outf);
}

void EnergyData::setup(gmx_mdoutf* outf)
{
    energyOutput_ = std::make_unique<EnergyOutput>(mdoutf_get_fp_ene(outf),
                                                   top_global_,
                                                   *inputrec_,
                                                   pullWork_,
                                                   mdoutf_get_fp_dhdl(outf),
                                                   false,
                                                   startingBehavior_,
                                                   simulationsShareState_,
                                                   mdModulesNotifiers_);

    if (!isMasterRank_)
    {
        return;
    }

    initializeEnergyHistory(startingBehavior_, observablesHistory_, energyOutput_.get());

    // TODO: This probably doesn't really belong here...
    //       but we have all we need in this element,
    //       so we'll leave it here for now!
    double io = compute_io(inputrec_, top_global_.natoms, *groups_, energyOutput_->numEnergyTerms(), 1);
    if ((io > 2000) && isMasterRank_)
    {
        fprintf(stderr, "\nWARNING: This run will generate roughly %.0f Mb of data\n\n", io);
    }
    if (!inputrec_->bContinuation)
    {
        real temp = enerd_->term[F_TEMP];
        if (inputrec_->eI != IntegrationAlgorithm::VV)
        {
            /* Result of Ekin averaged over velocities of -half
             * and +half step, while we only have -half step here.
             */
            temp *= 2;
        }
        fprintf(fplog_, "Initial temperature: %g K\n", temp);
    }
}

std::optional<ITrajectoryWriterCallback> EnergyData::Element::registerTrajectoryWriterCallback(TrajectoryEvent event)
{
    if (event == TrajectoryEvent::EnergyWritingStep && isMasterRank_)
    {
        return [this](gmx_mdoutf* mdoutf, Step step, Time time, bool writeTrajectory, bool writeLog) {
            energyData_->write(mdoutf, step, time, writeTrajectory, writeLog);
        };
    }
    return std::nullopt;
}

std::optional<SignallerCallback> EnergyData::Element::registerTrajectorySignallerCallback(gmx::TrajectoryEvent event)
{
    if (event == TrajectoryEvent::EnergyWritingStep && isMasterRank_)
    {
        return [this](Step step, Time /*unused*/) { energyWritingStep_ = step; };
    }
    return std::nullopt;
}

std::optional<SignallerCallback> EnergyData::Element::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::EnergyCalculationStep && isMasterRank_)
    {
        return [this](Step step, Time /*unused*/) { energyCalculationStep_ = step; };
    }
    if (event == EnergySignallerEvent::FreeEnergyCalculationStep && isMasterRank_)
    {
        return [this](Step step, Time /*unused*/) { freeEnergyCalculationStep_ = step; };
    }
    return std::nullopt;
}

void EnergyData::doStep(Step step, Time time, bool isEnergyCalculationStep, bool isFreeEnergyCalculationStep)
{
    enerd_->term[F_ETOT] = enerd_->term[F_EPOT] + enerd_->term[F_EKIN];
    if (freeEnergyPerturbationData_)
    {
        accumulateKineticLambdaComponents(
                enerd_, freeEnergyPerturbationData_->constLambdaView(), *inputrec_->fepvals);
    }
    if (integratorHasConservedEnergyQuantity(inputrec_))
    {
        enerd_->term[F_ECONSERVED] = enerd_->term[F_ETOT];
        for (const auto& energyContibution : conservedEnergyContributions_)
        {
            enerd_->term[F_ECONSERVED] += energyContibution(step, time);
        }
    }
    matrix nullMatrix = {};
    energyOutput_->addDataAtEnergyStep(
            isFreeEnergyCalculationStep,
            isEnergyCalculationStep,
            time,
            mdAtoms_->mdatoms()->tmass,
            enerd_,
            inputrec_->fepvals.get(),
            statePropagatorData_->constPreviousBox(),
            PTCouplingArrays({ parrinelloRahmanBoxVelocities_ ? parrinelloRahmanBoxVelocities_() : nullMatrix,
                               {},
                               {},
                               {},
                               {} }),
            freeEnergyPerturbationData_ ? freeEnergyPerturbationData_->currentFEPState() : 0,
            totalVirial_,
            pressure_,
            ekind_,
            muTot_,
            constr_);
}

void EnergyData::write(gmx_mdoutf* outf, Step step, Time time, bool writeTrajectory, bool writeLog)
{
    if (writeLog)
    {
        energyOutput_->printHeader(fplog_, step, time);
    }

    bool do_dr = do_per_step(step, inputrec_->nstdisreout);
    bool do_or = do_per_step(step, inputrec_->nstorireout);

    // energyOutput_->printAnnealingTemperatures(writeLog ? fplog_ : nullptr, groups_, &(inputrec_->opts));
    Awh* awh = nullptr;
    energyOutput_->printStepToEnergyFile(
            mdoutf_get_fp_ene(outf), writeTrajectory, do_dr, do_or, writeLog ? fplog_ : nullptr, step, time, fcd_, awh);
}

void EnergyData::addToForceVirial(const tensor virial, Step step)
{
    if (step > forceVirialStep_)
    {
        forceVirialStep_ = step;
        clear_mat(forceVirial_);
    }
    m_add(forceVirial_, virial, forceVirial_);
}

void EnergyData::addToConstraintVirial(const tensor virial, Step step)
{
    if (step > shakeVirialStep_)
    {
        shakeVirialStep_ = step;
        clear_mat(shakeVirial_);
    }
    m_add(shakeVirial_, virial, shakeVirial_);
}

rvec* EnergyData::forceVirial(Step gmx_unused step)
{
    if (step > forceVirialStep_)
    {
        forceVirialStep_ = step;
        clear_mat(forceVirial_);
    }
    GMX_ASSERT(step >= forceVirialStep_ || forceVirialStep_ == -1,
               "Asked for force virial of previous step.");
    return forceVirial_;
}

rvec* EnergyData::constraintVirial(Step gmx_unused step)
{
    if (step > shakeVirialStep_)
    {
        shakeVirialStep_ = step;
        clear_mat(shakeVirial_);
    }
    GMX_ASSERT(step >= shakeVirialStep_ || shakeVirialStep_ == -1,
               "Asked for constraint virial of previous step.");
    return shakeVirial_;
}

rvec* EnergyData::totalVirial(Step gmx_unused step)
{
    if (step > totalVirialStep_)
    {
        totalVirialStep_ = step;
        clear_mat(totalVirial_);
    }
    GMX_ASSERT(step >= totalVirialStep_ || totalVirialStep_ == -1,
               "Asked for total virial of previous step.");
    return totalVirial_;
}

rvec* EnergyData::pressure(Step gmx_unused step)
{
    if (step > pressureStep_)
    {
        pressureStep_ = step;
        clear_mat(pressure_);
    }
    GMX_ASSERT(step >= pressureStep_ || pressureStep_ == -1, "Asked for pressure of previous step.");
    return pressure_;
}

real* EnergyData::muTot()
{
    return muTot_;
}

gmx_enerdata_t* EnergyData::enerdata()
{
    return enerd_;
}

const gmx_enerdata_t* EnergyData::enerdata() const
{
    return enerd_;
}

gmx_ekindata_t* EnergyData::ekindata()
{
    return ekind_;
}

bool* EnergyData::needToSumEkinhOld()
{
    return &needToSumEkinhOld_;
}

bool EnergyData::hasReadEkinFromCheckpoint() const
{
    return hasReadEkinFromCheckpoint_;
}

namespace
{
/*!
 * \brief Enum describing the contents EnergyData::Element writes to modular checkpoint
 *
 * When changing the checkpoint content, add a new element just above Count, and adjust the
 * checkpoint functionality.
 */
enum class CheckpointVersion
{
    Base, //!< First version of modular checkpointing
    Count //!< Number of entries. Add new versions right above this!
};
constexpr auto c_currentVersion = CheckpointVersion(int(CheckpointVersion::Count) - 1);
} // namespace

template<CheckpointDataOperation operation>
void EnergyData::Element::doCheckpointData(CheckpointData<operation>* checkpointData)
{
    checkpointVersion(checkpointData, "EnergyData version", c_currentVersion);

    energyData_->observablesHistory_->energyHistory->doCheckpoint<operation>(
            checkpointData->subCheckpointData("energy history"));
    energyData_->ekinstate_.doCheckpoint<operation>(checkpointData->subCheckpointData("ekinstate"));
}

void EnergyData::Element::saveCheckpointState(std::optional<WriteCheckpointData> checkpointData,
                                              const t_commrec*                   cr)
{
    if (MASTER(cr))
    {
        if (energyData_->needToSumEkinhOld_)
        {
            energyData_->ekinstate_.bUpToDate = false;
        }
        else
        {
            update_ekinstate(&energyData_->ekinstate_, energyData_->ekind_);
            energyData_->ekinstate_.bUpToDate = true;
        }
        energyData_->energyOutput_->fillEnergyHistory(
                energyData_->observablesHistory_->energyHistory.get());
        doCheckpointData<CheckpointDataOperation::Write>(&checkpointData.value());
    }
}

void EnergyData::Element::restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData,
                                                 const t_commrec*                  cr)
{
    if (MASTER(cr))
    {
        doCheckpointData<CheckpointDataOperation::Read>(&checkpointData.value());
    }
    energyData_->hasReadEkinFromCheckpoint_ = MASTER(cr) ? energyData_->ekinstate_.bUpToDate : false;
    if (PAR(cr))
    {
        gmx_bcast(sizeof(hasReadEkinFromCheckpoint_),
                  &energyData_->hasReadEkinFromCheckpoint_,
                  cr->mpi_comm_mygroup);
    }
    if (energyData_->hasReadEkinFromCheckpoint_)
    {
        // this takes care of broadcasting from master to agents
        restore_ekinstate_from_state(cr, energyData_->ekind_, &energyData_->ekinstate_);
    }
}

const std::string& EnergyData::Element::clientID()
{
    return identifier_;
}

void EnergyData::initializeEnergyHistory(StartingBehavior    startingBehavior,
                                         ObservablesHistory* observablesHistory,
                                         EnergyOutput*       energyOutput)
{
    if (startingBehavior != StartingBehavior::NewSimulation)
    {
        /* Restore from energy history if appending to output files */
        if (startingBehavior == StartingBehavior::RestartWithAppending)
        {
            /* If no history is available (because a checkpoint is from before
             * it was written) make a new one later, otherwise restore it.
             */
            if (observablesHistory->energyHistory)
            {
                energyOutput->restoreFromEnergyHistory(*observablesHistory->energyHistory);
            }
        }
        else if (observablesHistory->energyHistory)
        {
            /* We might have read an energy history from checkpoint.
             * As we are not appending, we want to restart the statistics.
             * Free the allocated memory and reset the counts.
             */
            observablesHistory->energyHistory = {};
            /* We might have read a pull history from checkpoint.
             * We will still want to keep the statistics, so that the files
             * can be joined and still be meaningful.
             * This means that observablesHistory_->pullHistory
             * should not be reset.
             */
        }
    }
    if (!observablesHistory->energyHistory)
    {
        observablesHistory->energyHistory = std::make_unique<energyhistory_t>();
    }
    if (!observablesHistory->pullHistory)
    {
        observablesHistory->pullHistory = std::make_unique<PullHistory>();
    }
    /* Set the initial energy history */
    energyOutput->fillEnergyHistory(observablesHistory->energyHistory.get());
}

void EnergyData::addConservedEnergyContribution(EnergyContribution&& energyContribution)
{
    conservedEnergyContributions_.emplace_back(std::move(energyContribution));
}

void EnergyData::setParrinelloRahmanBoxVelocities(std::function<const rvec*()>&& parrinelloRahmanBoxVelocities)
{
    GMX_RELEASE_ASSERT(!parrinelloRahmanBoxVelocities_,
                       "Received a second callback to the Parrinello-Rahman velocities");
    parrinelloRahmanBoxVelocities_ = parrinelloRahmanBoxVelocities;
}

void EnergyData::updateKineticEnergy()
{
    // The legacy sum_ekin function does not offer named types, so define variables for readability
    // dEkin/dlambda is not handled here
    real* dEkinDLambda = nullptr;
    // Whether we use the full step kinetic energy (vs the average of half step KEs)
    const bool useFullStepKineticEnergy = (inputrec_->eI == IntegrationAlgorithm::VV);
    /* Whether we're ignoring the NHC scaling factor, only used if useFullStepKineticEnergy
     * is true. (This parameter is confusing, as it is named `bScaleEkin`, but prompts the
     * function to ignore scaling. There is no use case within modular simulator to ignore
     * these, so we set this to false.) */
    const bool ignoreScalingFactor = false;

    enerd_->term[F_TEMP] = sum_ekin(
            &(inputrec_->opts), ekind_, dEkinDLambda, useFullStepKineticEnergy, ignoreScalingFactor);
    enerd_->term[F_EKIN] = trace(ekind_->ekin);
}

EnergyData::Element* EnergyData::element()
{
    return element_.get();
}

EnergyData::Element::Element(EnergyData* energyData, bool isMasterRank) :
    energyData_(energyData),
    isMasterRank_(isMasterRank),
    energyWritingStep_(-1),
    energyCalculationStep_(-1),
    freeEnergyCalculationStep_(-1)
{
}

ISimulatorElement* EnergyData::Element::getElementPointerImpl(
        LegacySimulatorData gmx_unused*        legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper gmx_unused* builderHelper,
        StatePropagatorData gmx_unused* statePropagatorData,
        EnergyData*                     energyData,
        FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer* /*observablesReducer*/)
{
    return energyData->element();
}

} // namespace gmx
