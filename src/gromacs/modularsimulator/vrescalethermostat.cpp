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
/*! \internal \file
 * \brief Defines the v-rescale thermostat for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "vrescalethermostat.h"

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

VRescaleThermostat::VRescaleThermostat(int                   nstcouple,
                                       int                   offset,
                                       bool                  useFullStepKE,
                                       int64_t               seed,
                                       int                   numTemperatureGroups,
                                       double                couplingTimeStep,
                                       const real*           referenceTemperature,
                                       const real*           couplingTime,
                                       const real*           numDegreesOfFreedom,
                                       EnergyElement*        energyElement,
                                       ArrayRef<real>        lambdaView,
                                       PropagatorCallbackPtr propagatorCallback,
                                       const t_state*        globalState,
                                       t_commrec*            cr,
                                       bool                  isRestart) :
    nstcouple_(nstcouple),
    offset_(offset),
    useFullStepKE_(useFullStepKE),
    seed_(seed),
    numTemperatureGroups_(numTemperatureGroups),
    couplingTimeStep_(couplingTimeStep),
    referenceTemperature_(referenceTemperature, referenceTemperature + numTemperatureGroups),
    couplingTime_(couplingTime, couplingTime + numTemperatureGroups),
    numDegreesOfFreedom_(numDegreesOfFreedom, numDegreesOfFreedom + numTemperatureGroups),
    thermostatIntegral_(numTemperatureGroups, 0.0),
    energyElement_(energyElement),
    lambda_(lambdaView),
    propagatorCallback_(std::move(propagatorCallback))
{
    // TODO: This is only needed to restore the thermostatIntegral_ from cpt. Remove this when
    //       switching to purely client-based checkpointing.
    if (isRestart)
    {
        if (MASTER(cr))
        {
            for (unsigned long i = 0; i < thermostatIntegral_.size(); ++i)
            {
                thermostatIntegral_[i] = globalState->therm_integral[i];
            }
        }
        if (DOMAINDECOMP(cr))
        {
            dd_bcast(cr->dd, int(thermostatIntegral_.size() * sizeof(double)), thermostatIntegral_.data());
        }
    }
}

void VRescaleThermostat::scheduleTask(Step step, Time gmx_unused time, const RegisterRunFunctionPtr& registerRunFunction)
{
    /* The thermostat will need a valid kinetic energy when it is running.
     * Currently, computeGlobalCommunicationPeriod() is making sure this
     * happens on time.
     * TODO: Once we're switching to a new global communication scheme, we
     *       will want the thermostat to signal that global reduction
     *       of the kinetic energy is needed.
     *
     */
    if (do_per_step(step + nstcouple_ + offset_, nstcouple_))
    {
        // do T-coupling this step
        (*registerRunFunction)(
                std::make_unique<SimulatorRunFunction>([this, step]() { setLambda(step); }));

        // Let propagator know that we want to do T-coupling
        (*propagatorCallback_)(step);
    }
}

void VRescaleThermostat::setLambda(Step step)
{
    real currentKineticEnergy, referenceKineticEnergy, newKineticEnergy;

    auto ekind = energyElement_->ekindata();

    for (int i = 0; (i < numTemperatureGroups_); i++)
    {
        if (useFullStepKE_)
        {
            currentKineticEnergy = trace(ekind->tcstat[i].ekinf);
        }
        else
        {
            currentKineticEnergy = trace(ekind->tcstat[i].ekinh);
        }

        if (couplingTime_[i] >= 0 && numDegreesOfFreedom_[i] > 0 && currentKineticEnergy > 0)
        {
            referenceKineticEnergy = 0.5 * referenceTemperature_[i] * BOLTZ * numDegreesOfFreedom_[i];

            newKineticEnergy = vrescale_resamplekin(currentKineticEnergy, referenceKineticEnergy,
                                                    numDegreesOfFreedom_[i],
                                                    couplingTime_[i] / couplingTimeStep_, step, seed_);

            // Analytically newKineticEnergy >= 0, but we check for rounding errors
            if (newKineticEnergy <= 0)
            {
                lambda_[i] = 0.0;
            }
            else
            {
                lambda_[i] = std::sqrt(newKineticEnergy / currentKineticEnergy);
            }

            thermostatIntegral_[i] -= newKineticEnergy - currentKineticEnergy;

            if (debug)
            {
                fprintf(debug, "TC: group %d: Ekr %g, Ek %g, Ek_new %g, Lambda: %g\n", i,
                        referenceKineticEnergy, currentKineticEnergy, newKineticEnergy, lambda_[i]);
            }
        }
        else
        {
            lambda_[i] = 1.0;
        }
    }
}

void VRescaleThermostat::writeCheckpoint(t_state* localState, t_state gmx_unused* globalState)
{
    localState->therm_integral = thermostatIntegral_;
    localState->flags |= (1U << estTHERM_INT);
}

const std::vector<double>& VRescaleThermostat::thermostatIntegral() const
{
    return thermostatIntegral_;
}

} // namespace gmx
