/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Defines the constraint element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "constraintelement.h"

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/fatalerror.h"

#include "energydata.h"
#include "freeenergyperturbationdata.h"
#include "modularsimulator.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"

namespace gmx
{
template<ConstraintVariable variable>
ConstraintsElement<variable>::ConstraintsElement(Constraints*                constr,
                                                 StatePropagatorData*        statePropagatorData,
                                                 EnergyData*                 energyData,
                                                 FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                 bool                        isMain,
                                                 FILE*                       fplog,
                                                 const t_inputrec*           inputrec,
                                                 const t_mdatoms*            mdAtoms) :
    nextVirialCalculationStep_(-1),
    nextEnergyWritingStep_(-1),
    nextLogWritingStep_(-1),
    isMainRank_(isMain),
    statePropagatorData_(statePropagatorData),
    energyData_(energyData),
    freeEnergyPerturbationData_(freeEnergyPerturbationData),
    constr_(constr),
    fplog_(fplog),
    inputrec_(inputrec),
    mdAtoms_(mdAtoms)
{
    GMX_ASSERT(constr_, "Constraint element created but constr == nullptr");
}

template<ConstraintVariable variable>
void ConstraintsElement<variable>::elementSetup()
{
    if (!inputrec_->bContinuation
        && ((variable == ConstraintVariable::Positions && inputrec_->eI == IntegrationAlgorithm::MD)
            || (variable == ConstraintVariable::Velocities && inputrec_->eI == IntegrationAlgorithm::VV)))
    {
        const real lambdaBonded =
                freeEnergyPerturbationData_
                        ? freeEnergyPerturbationData_->constLambdaView()[static_cast<int>(
                                FreeEnergyPerturbationCouplingType::Bonded)]
                        : 0;
        // Constrain the initial coordinates and velocities
        do_constrain_first(fplog_,
                           constr_,
                           *inputrec_,
                           statePropagatorData_->totalNumAtoms(),
                           statePropagatorData_->localNumAtoms(),
                           statePropagatorData_->positionsView(),
                           statePropagatorData_->velocitiesView(),
                           statePropagatorData_->box(),
                           lambdaBonded);

        if (isMainRank_)
        {
            if (inputrec_->eConstrAlg == ConstraintAlgorithm::Lincs)
            {
                fprintf(fplog_,
                        "RMS relative constraint deviation after constraining: %.2e\n",
                        constr_->rmsd());
            }
        }
    }
}

template<ConstraintVariable variable>
void ConstraintsElement<variable>::scheduleTask(Step                       step,
                                                Time gmx_unused            time,
                                                const RegisterRunFunction& registerRunFunction)
{
    bool calculateVirial = (step == nextVirialCalculationStep_);
    bool writeLog        = (step == nextLogWritingStep_);
    bool writeEnergy     = (step == nextEnergyWritingStep_);

    // register constraining
    registerRunFunction([this, step, calculateVirial, writeLog, writeEnergy]() {
        apply(step, calculateVirial, writeLog, writeEnergy);
    });
}

template<ConstraintVariable variable>
void ConstraintsElement<variable>::apply(Step step, bool calculateVirial, bool writeLog, bool writeEnergy)
{
    tensor vir_con;

    ArrayRefWithPadding<RVec> x;
    ArrayRefWithPadding<RVec> xprime;
    ArrayRef<RVec>            min_proj;
    ArrayRefWithPadding<RVec> v;

    const real lambdaBonded =
            freeEnergyPerturbationData_
                    ? freeEnergyPerturbationData_
                              ->constLambdaView()[static_cast<int>(FreeEnergyPerturbationCouplingType::Bonded)]
                    : 0;
    real dvdlambda = 0;

    switch (variable)
    {
        case ConstraintVariable::Positions:
            x      = statePropagatorData_->previousPositionsView();
            xprime = statePropagatorData_->positionsView();
            v      = statePropagatorData_->velocitiesView();
            break;
        case ConstraintVariable::Velocities:
            x        = statePropagatorData_->positionsView();
            xprime   = statePropagatorData_->velocitiesView();
            min_proj = statePropagatorData_->velocitiesView().unpaddedArrayRef();
            break;
        default: gmx_fatal(FARGS, "Constraint algorithm not implemented for modular simulator.");
    }

    constr_->apply(writeLog || writeEnergy,
                   step,
                   1,
                   1.0,
                   x,
                   xprime,
                   min_proj,
                   statePropagatorData_->box(),
                   lambdaBonded,
                   &dvdlambda,
                   v,
                   calculateVirial,
                   vir_con,
                   variable);

    if (calculateVirial)
    {
        if (inputrec_->eI == IntegrationAlgorithm::VV)
        {
            // For some reason, the shake virial in VV is reset twice a step.
            // Energy element will only do this once per step.
            // TODO: Investigate this
            clear_mat(energyData_->constraintVirial(step));
        }
        energyData_->addToConstraintVirial(vir_con, step);
    }

    /* The factor of 2 correction is necessary because half of the constraint
     * force is removed in the VV step. This factor is either exact or a very
     * good approximation, statistically insignificant in any real free energy
     * calculation. Any possible error is not a simulation propagation error,
     * but a potential reporting error in the data that goes to dh/dlambda.
     * Cf. Issue #1255
     */
    const real c_dvdlConstraintCorrectionFactor = EI_VV(inputrec_->eI) ? 2.0 : 1.0;
    energyData_->enerdata()->term[F_DVDL_CONSTR] += c_dvdlConstraintCorrectionFactor * dvdlambda;
}

template<ConstraintVariable variable>
std::optional<SignallerCallback> ConstraintsElement<variable>::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::VirialCalculationStep)
    {
        return [this](Step step, Time /*unused*/) { nextVirialCalculationStep_ = step; };
    }
    return std::nullopt;
}

template<ConstraintVariable variable>
std::optional<SignallerCallback> ConstraintsElement<variable>::registerTrajectorySignallerCallback(TrajectoryEvent event)
{
    if (event == TrajectoryEvent::EnergyWritingStep)
    {
        return [this](Step step, Time /*unused*/) { nextEnergyWritingStep_ = step; };
    }
    return std::nullopt;
}

template<ConstraintVariable variable>
std::optional<SignallerCallback> ConstraintsElement<variable>::registerLoggingCallback()
{
    return [this](Step step, Time /*unused*/) { nextLogWritingStep_ = step; };
}

template<ConstraintVariable variable>
ISimulatorElement* ConstraintsElement<variable>::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData*                    statePropagatorData,
        EnergyData*                             energyData,
        FreeEnergyPerturbationData*             freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer* /*observablesReducer*/)
{
    return builderHelper->storeElement(std::make_unique<ConstraintsElement<variable>>(
            legacySimulatorData->constr_,
            statePropagatorData,
            energyData,
            freeEnergyPerturbationData,
            MAIN(legacySimulatorData->cr_),
            legacySimulatorData->fpLog_,
            legacySimulatorData->inputRec_,
            legacySimulatorData->mdAtoms_->mdatoms()));
}

// Explicit template initializations
template class ConstraintsElement<ConstraintVariable::Positions>;
template class ConstraintsElement<ConstraintVariable::Velocities>;

} // namespace gmx
