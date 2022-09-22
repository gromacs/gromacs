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
 * \brief Defines the Parrinello-Rahman barostat for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "parrinellorahmanbarostat.h"

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/boxutilities.h"

#include "energydata.h"
#include "modularsimulator.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"

namespace gmx
{

ParrinelloRahmanBarostat::ParrinelloRahmanBarostat(int                  nstpcouple,
                                                   int                  offset,
                                                   real                 couplingTimePeriod,
                                                   Step                 initStep,
                                                   StatePropagatorData* statePropagatorData,
                                                   EnergyData*          energyData,
                                                   const MDLogger&      mdlog,
                                                   const t_inputrec*    inputrec,
                                                   const MDAtoms*       mdAtoms) :
    nstpcouple_(nstpcouple),
    offset_(offset),
    couplingTimePeriod_(couplingTimePeriod),
    initStep_(initStep),
    mu_{ { 0 } },
    boxRel_{ { 0 } },
    boxVelocity_{ { 0 } },
    statePropagatorData_(statePropagatorData),
    energyData_(energyData),
    nextEnergyCalculationStep_(-1),
    mdlog_(mdlog),
    inputrec_(inputrec),
    mdAtoms_(mdAtoms)
{
    energyData->setParrinelloRahmanBoxVelocities([this]() { return boxVelocity_; });
    energyData->addConservedEnergyContribution([this](Step gmx_used_in_debug step, Time /*unused*/) {
        GMX_ASSERT(conservedEnergyContributionStep_ == step,
                   "Parrinello-Rahman conserved energy step mismatch.");
        return conservedEnergyContribution_;
    });
}

void ParrinelloRahmanBarostat::connectWithMatchingPropagator(const PropagatorConnection& connectionData,
                                                             const PropagatorTag& propagatorTag)
{
    if (connectionData.tag == propagatorTag)
    {
        GMX_RELEASE_ASSERT(connectionData.hasParrinelloRahmanScaling(),
                           "Connection data lacks Parrinello-Rahman scaling");
        scalingTensor_      = connectionData.getViewOnPRScalingMatrix();
        propagatorCallback_ = connectionData.getPRScalingCallback();
    }
}

void ParrinelloRahmanBarostat::scheduleTask(Step                       step,
                                            Time gmx_unused            time,
                                            const RegisterRunFunction& registerRunFunction)
{
    const bool scaleOnNextStep = do_per_step(step + nstpcouple_ + offset_ + 1, nstpcouple_);
    const bool scaleOnThisStep = do_per_step(step + nstpcouple_ + offset_, nstpcouple_);
    const bool contributeEnergyThisStep = (step == nextEnergyCalculationStep_);

    if (contributeEnergyThisStep)
    {
        // For compatibility with legacy md, we store this before integrating the box velocities
        registerRunFunction([this, step]() {
            conservedEnergyContribution_     = conservedEnergyContribution();
            conservedEnergyContributionStep_ = step;
        });
    }
    if (scaleOnThisStep)
    {
        registerRunFunction([this]() { scaleBoxAndPositions(); });
    }
    if (scaleOnNextStep)
    {
        registerRunFunction([this, step]() { integrateBoxVelocityEquations(step); });
        // let propagator know that it will have to scale on next step
        propagatorCallback_(step + 1);
    }
}

void ParrinelloRahmanBarostat::integrateBoxVelocityEquations(Step step)
{
    const auto* box = statePropagatorData_->constBox();
    parrinellorahman_pcoupl(mdlog_,
                            step,
                            inputrec_->pressureCouplingOptions,
                            inputrec_->deform,
                            couplingTimePeriod_,
                            energyData_->pressure(step),
                            box,
                            boxRel_,
                            boxVelocity_,
                            scalingTensor_,
                            &mu_);
    // multiply matrix by the coupling time step to avoid having the propagator needing to know about that
    (*scalingTensor_) = (*scalingTensor_) * couplingTimePeriod_;
}

void ParrinelloRahmanBarostat::scaleBoxAndPositions()
{
    // Propagate the box by the box velocities
    auto* box = statePropagatorData_->box();
    for (int i = 0; i < DIM; i++)
    {
        for (int m = 0; m <= i; m++)
        {
            box[i][m] += couplingTimePeriod_ * boxVelocity_[i][m];
        }
    }
    preserveBoxShape(inputrec_->pressureCouplingOptions, inputrec_->deform, boxRel_, box);

    // Scale the coordinates
    const int      start   = 0;
    const int      homenr  = mdAtoms_->mdatoms()->homenr;
    ArrayRef<RVec> x       = statePropagatorData_->positionsView().paddedArrayRef();
    ivec*          nFreeze = inputrec_->opts.nFreeze;
    for (int n = start; n < start + homenr; n++)
    {
        if (mdAtoms_->mdatoms()->cFREEZE.empty())
        {
            x[n] = multiplyVectorByTransposeOfBoxMatrix(mu_, x[n]);
        }
        else
        {
            int g = mdAtoms_->mdatoms()->cFREEZE[n];
            if (!nFreeze[g][XX])
            {
                x[n][XX] = mu_(XX, XX) * x[n][XX] + mu_(YY, XX) * x[n][YY] + mu_(ZZ, XX) * x[n][ZZ];
            }
            if (!nFreeze[g][YY])
            {
                x[n][YY] = mu_(YY, YY) * x[n][YY] + mu_(ZZ, YY) * x[n][ZZ];
            }
            if (!nFreeze[g][ZZ])
            {
                x[n][ZZ] = mu_(ZZ, ZZ) * x[n][ZZ];
            }
        }
    }
}

void ParrinelloRahmanBarostat::elementSetup()
{
    if (!propagatorCallback_ || scalingTensor_ == nullptr || scalingTensor_->asConstView().rank() == 0)
    {
        throw MissingElementConnectionError(
                "Parrinello-Rahman barostat was not connected to a propagator.\n"
                "Connection to a propagator element is needed to scale the velocities.\n"
                "Use connectWithMatchingPropagator(...) before building the "
                "ModularSimulatorAlgorithm "
                "object.");
    }

    if (shouldPreserveBoxShape(inputrec_->pressureCouplingOptions, inputrec_->deform))
    {
        auto*     box = statePropagatorData_->box();
        const int ndim =
                inputrec_->pressureCouplingOptions.epct == PressureCouplingType::SemiIsotropic ? 2 : 3;
        do_box_rel(ndim, inputrec_->deform, boxRel_, box, true);
    }

    const bool scaleOnInitStep = do_per_step(initStep_ + nstpcouple_ + offset_, nstpcouple_);
    if (scaleOnInitStep)
    {
        // If we need to scale on the first step, we need to set the scaling matrix using the
        // current box velocity. If this is a fresh start, we will hence not move the box (this does
        // currently never happen as the offset is set to -1 in all cases). If this is a restart, we
        // will use the saved box velocity which we would have updated right before checkpointing.
        const auto* box = statePropagatorData_->constBox();
        init_parrinellorahman(inputrec_->pressureCouplingOptions,
                              inputrec_->deform,
                              couplingTimePeriod_,
                              box,
                              boxRel_,
                              boxVelocity_,
                              scalingTensor_,
                              &mu_);
        // multiply matrix by the coupling time step to avoid having the propagator needing to know about that
        (*scalingTensor_) = (*scalingTensor_) * couplingTimePeriod_;

        propagatorCallback_(initStep_);
    }
}

const rvec* ParrinelloRahmanBarostat::boxVelocities() const
{
    return boxVelocity_;
}

real ParrinelloRahmanBarostat::conservedEnergyContribution() const
{
    real        energy       = 0;
    const auto* box          = statePropagatorData_->constBox();
    real        maxBoxLength = std::max({ box[XX][XX], box[YY][YY], box[ZZ][ZZ] });
    real        volume       = det(box);

    // contribution from the pressure momenta
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            real invMass = c_presfac
                           * (4 * M_PI * M_PI * inputrec_->pressureCouplingOptions.compress[i][j])
                           / (3 * inputrec_->pressureCouplingOptions.tau_p
                              * inputrec_->pressureCouplingOptions.tau_p * maxBoxLength);
            if (invMass > 0)
            {
                energy += 0.5 * boxVelocity_[i][j] * boxVelocity_[i][j] / invMass;
            }
        }
    }

    /* Contribution from the PV term.
     * Note that with non-zero off-diagonal reference pressures,
     * i.e. applied shear stresses, there are additional terms.
     * We don't support this here, since that requires keeping
     * track of unwrapped box diagonal elements. This case is
     * excluded in integratorHasConservedEnergyQuantity().
     */
    energy += volume * ::trace(inputrec_->pressureCouplingOptions.ref_p) / (DIM * c_presfac);

    return energy;
}

namespace
{
/*!
 * \brief Enum describing the contents ParrinelloRahmanBarostat writes to modular checkpoint
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
void ParrinelloRahmanBarostat::doCheckpointData(CheckpointData<operation>* checkpointData)
{
    checkpointVersion(checkpointData, "ParrinelloRahmanBarostat version", c_currentVersion);

    checkpointData->tensor("box velocity", boxVelocity_);
    checkpointData->tensor("relative box vector", boxRel_);
}

void ParrinelloRahmanBarostat::saveCheckpointState(std::optional<WriteCheckpointData> checkpointData,
                                                   const t_commrec*                   cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Write>(&checkpointData.value());
    }
}

void ParrinelloRahmanBarostat::restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData,
                                                      const t_commrec*                  cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Read>(&checkpointData.value());
    }
    if (haveDDAtomOrdering(*cr))
    {
        dd_bcast(cr->dd, sizeof(boxVelocity_), boxVelocity_);
        dd_bcast(cr->dd, sizeof(boxRel_), boxRel_);
    }
}

const std::string& ParrinelloRahmanBarostat::clientID()
{
    return identifier_;
}

std::optional<SignallerCallback> ParrinelloRahmanBarostat::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::EnergyCalculationStep)
    {
        return [this](Step step, Time /*unused*/) { nextEnergyCalculationStep_ = step; };
    }
    return std::nullopt;
}

ISimulatorElement* ParrinelloRahmanBarostat::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData*                    statePropagatorData,
        EnergyData*                             energyData,
        FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer* /*observablesReducer*/,
        Offset               offset,
        const PropagatorTag& propagatorTag)
{
    auto* element  = builderHelper->storeElement(std::make_unique<ParrinelloRahmanBarostat>(
            legacySimulatorData->inputrec->pressureCouplingOptions.nstpcouple,
            offset,
            legacySimulatorData->inputrec->delta_t
                    * legacySimulatorData->inputrec->pressureCouplingOptions.nstpcouple,
            legacySimulatorData->inputrec->init_step,
            statePropagatorData,
            energyData,
            legacySimulatorData->mdlog,
            legacySimulatorData->inputrec,
            legacySimulatorData->mdAtoms));
    auto* barostat = static_cast<ParrinelloRahmanBarostat*>(element);
    builderHelper->registerTemperaturePressureControl(
            [barostat, propagatorTag](const PropagatorConnection& connection) {
                barostat->connectWithMatchingPropagator(connection, propagatorTag);
            });
    return element;
}

} // namespace gmx
