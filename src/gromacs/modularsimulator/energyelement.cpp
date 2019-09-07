/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal
 * \brief Defines the microstate for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "energyelement.h"

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"

#include "statepropagatordata.h"

struct pull_t;
class t_state;

namespace gmx
{
class Awh;

EnergyElement::EnergyElement(
        StatePropagatorData     *statePropagatorData,
        const gmx_mtop_t        *globalTopology,
        const t_inputrec        *inputrec,
        const MDAtoms           *mdAtoms,
        gmx_enerdata_t          *enerd,
        gmx_ekindata_t          *ekind,
        const Constraints       *constr,
        FILE                    *fplog,
        t_fcdata                *fcd,
        const MdModulesNotifier &mdModulesNotifier,
        bool                     isMaster) :
    isMaster_(isMaster),
    energyWritingStep_(-1),
    energyCalculationStep_(-1),
    freeEnergyCalculationStep_(-1),
    logWritingStep_(-1),
    forceVirialStep_(-1),
    shakeVirialStep_(-1),
#ifndef NDEBUG
    totalVirialStep_(-1),
    pressureStep_(-1),
#endif
    statePropagatorData_(statePropagatorData),
    inputrec_(inputrec),
    top_global_(globalTopology),
    mdAtoms_(mdAtoms),
    enerd_(enerd),
    ekind_(ekind),
    constr_(constr),
    fplog_(fplog),
    fcd_(fcd),
    mdModulesNotifier_(mdModulesNotifier),
    groups_(&globalTopology->groups)
{
    clear_mat(forceVirial_);
    clear_mat(shakeVirial_);
    clear_mat(totalVirial_);
    clear_mat(pressure_);
    clear_rvec(muTot_);

    // TODO: If energy history is re-introduced, it should probably be initialized here.
}

void EnergyElement::scheduleTask(
        Step step, Time time, const RegisterRunFunctionPtr &registerRunFunction)
{
    if (!isMaster_)
    {
        return;
    }
    auto writeEnergy                 = energyWritingStep_ == step;
    auto isEnergyCalculationStep     = energyCalculationStep_ == step;
    auto isFreeEnergyCalculationStep = freeEnergyCalculationStep_ == step;
    if (isEnergyCalculationStep || writeEnergy)
    {
        (*registerRunFunction)(
                std::make_unique<SimulatorRunFunction>(
                        [this, time, isEnergyCalculationStep, isFreeEnergyCalculationStep]() {
                            doStep(time, isEnergyCalculationStep, isFreeEnergyCalculationStep);
                        }));
    }
    else
    {
        (*registerRunFunction)(
                std::make_unique<SimulatorRunFunction>(
                        [this]() { energyOutput_->recordNonEnergyStep(); }));
    }
}

void EnergyElement::elementTeardown()
{
    if (inputrec_->nstcalcenergy > 0 && isMaster_)
    {
        energyOutput_->printAverages(fplog_, groups_);
    }
}

void EnergyElement::trajectoryWriterSetup(gmx_mdoutf *outf)
{
    pull_t *pull_work = nullptr;
    energyOutput_ = std::make_unique<EnergyOutput>(
                mdoutf_get_fp_ene(outf), top_global_,
                inputrec_, pull_work, mdoutf_get_fp_dhdl(outf), false,
                mdModulesNotifier_);

    if (!isMaster_)
    {
        return;
    }

    // TODO: This probably doesn't really belong here...
    //       but we have all we need in this element,
    //       so we'll leave it here for now!
    double io = compute_io(inputrec_, top_global_->natoms, *groups_,
                           energyOutput_->numEnergyTerms(), 1);
    if ((io > 2000) && isMaster_)
    {
        fprintf(stderr,
                "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                io);
    }
    if (!inputrec_->bContinuation)
    {
        real temp = enerd_->term[F_TEMP];
        if (inputrec_->eI != eiVV)
        {
            /* Result of Ekin averaged over velocities of -half
             * and +half step, while we only have -half step here.
             */
            temp *= 2;
        }
        fprintf(fplog_, "Initial temperature: %g K\n", temp);
    }
}

ITrajectoryWriterCallbackPtr EnergyElement::registerTrajectoryWriterCallback(TrajectoryEvent event)
{
    if (event == TrajectoryEvent::energyWritingStep && isMaster_)
    {
        return std::make_unique<ITrajectoryWriterCallback>(
                [this](gmx_mdoutf *mdoutf, Step step, Time time)
                {write(mdoutf, step, time); });
    }
    return nullptr;
}

SignallerCallbackPtr EnergyElement::registerTrajectorySignallerCallback(gmx::TrajectoryEvent event)
{
    if (event == TrajectoryEvent::energyWritingStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time){energyWritingStep_ = step; });
    }
    return nullptr;
}

SignallerCallbackPtr EnergyElement::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::energyCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time) {energyCalculationStep_ = step; });
    }
    if (event == EnergySignallerEvent::freeEnergyCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time){freeEnergyCalculationStep_ = step; });
    }
    return nullptr;

}

SignallerCallbackPtr EnergyElement::registerLoggingCallback()
{
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time){logWritingStep_ = step; });
}

void EnergyElement::doStep(
        Time time,
        bool isEnergyCalculationStep,
        bool isFreeEnergyCalculationStep)
{
    enerd_->term[F_ETOT] = enerd_->term[F_EPOT] + enerd_->term[F_EKIN];
    // All the state is used for is turned off for now (free energy, temperature / pressure coupling)
    t_state *localState = nullptr;
    energyOutput_->addDataAtEnergyStep(
            isFreeEnergyCalculationStep, isEnergyCalculationStep,
            time, mdAtoms_->mdatoms()->tmass, enerd_, localState,
            inputrec_->fepvals, inputrec_->expandedvals,
            statePropagatorData_->constPreviousBox(),
            shakeVirial_, forceVirial_, totalVirial_, pressure_,
            ekind_, muTot_, constr_);
}

void EnergyElement::write(gmx_mdoutf *outf, Step step, Time time)
{
    // This gets called by the trajectory writer, which means that
    // we only write to log when writing energy.
    // TODO: Think if we need to change that, i.e. writing to log without
    //       writing to file.

    auto writeLog = logWritingStep_ == step;
    if (writeLog)
    {
        energyOutput_->printHeader(fplog_, step, time);
    }

    bool do_dr  = do_per_step(step, inputrec_->nstdisreout);
    bool do_or  = do_per_step(step, inputrec_->nstorireout);

    // energyOutput_->printAnnealingTemperatures(writeLog ? fplog_ : nullptr, groups_, &(inputrec_->opts));
    Awh *awh = nullptr;
    energyOutput_->printStepToEnergyFile(
            mdoutf_get_fp_ene(outf), true, do_dr, do_or,
            writeLog ? fplog_ : nullptr, step, time, fcd_, awh);
}

void EnergyElement::addToForceVirial(const tensor virial, Step step)
{
    if (step > forceVirialStep_)
    {
        forceVirialStep_ = step;
        clear_mat(forceVirial_);
    }
    m_add(forceVirial_, virial, forceVirial_);
}

void EnergyElement::addToConstraintVirial(const tensor virial, Step step)
{
    if (step > shakeVirialStep_)
    {
        shakeVirialStep_ = step;
        clear_mat(shakeVirial_);
    }
    m_add(shakeVirial_, virial, shakeVirial_);
}

rvec* EnergyElement::forceVirial(Step gmx_unused step)
{
    GMX_ASSERT(step >= forceVirialStep_ || forceVirialStep_ == -1,
               "Asked for force virial of previous step.");
    return forceVirial_;
}

rvec* EnergyElement::constraintVirial(Step gmx_unused step)
{
    GMX_ASSERT(step >= shakeVirialStep_ || shakeVirialStep_ == -1,
               "Asked for constraint virial of previous step.");
    return shakeVirial_;
}

rvec* EnergyElement::totalVirial(Step gmx_unused step)
{
    GMX_ASSERT(step >= totalVirialStep_ || totalVirialStep_ == -1,
               "Asked for total virial of previous step.");
    return totalVirial_;
}

rvec* EnergyElement::pressure(Step gmx_unused step)
{
    GMX_ASSERT(step >= pressureStep_ || pressureStep_ == -1,
               "Asked for pressure of previous step.");
    return pressure_;
}

real* EnergyElement::muTot()
{
    return muTot_;
}

gmx_enerdata_t* EnergyElement::enerdata()
{
    return enerd_;
}

gmx_ekindata_t* EnergyElement::ekindata()
{
    return ekind_;
}

} // namespace gmx
