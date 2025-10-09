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
 * \brief Declares the PME load balancing helper for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "pmeloadbalancehelper.h"

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/nbnxm.h"

#include "statepropagatordata.h"

namespace gmx
{
bool PmeLoadBalanceHelper::doPmeLoadBalancing(const MdrunOptions&       mdrunOptions,
                                              const t_forcerec*         fr,
                                              const SimulationWorkload& simWorkload)
{
    return (mdrunOptions.tunePme
            && pmeTuningIsSupported(fr->ic->coulomb.type, mdrunOptions.reproducible, simWorkload));
}

PmeLoadBalanceHelper::PmeLoadBalanceHelper(bool                      isVerbose,
                                           StatePropagatorData*      statePropagatorData,
                                           gmx_domdec_t*             dd,
                                           const MDLogger&           mdlog,
                                           const t_inputrec*         inputrec,
                                           gmx_wallcycle*            wcycle,
                                           t_forcerec*               fr,
                                           const SimulationWorkload& simWorkload) :
    pme_loadbal_(dd, mdlog, *inputrec, statePropagatorData->constBox(), *fr->ic, *fr->nbv, fr->pmedata, simWorkload),
    nextNSStep_(-1),
    isVerbose_(isVerbose),
    statePropagatorData_(statePropagatorData),
    dd_(dd),
    inputrec_(inputrec),
    wcycle_(wcycle),
    fr_(fr)
{
}

void PmeLoadBalanceHelper::setup()
{
    const auto* box = statePropagatorData_->constBox();
    GMX_RELEASE_ASSERT(box[0][0] != 0 && box[1][1] != 0 && box[2][2] != 0,
                       "PmeLoadBalanceHelper cannot be initialized with zero box.");
}

void PmeLoadBalanceHelper::run(gmx::Step step, gmx::Time gmx_unused time)
{
    if (step != nextNSStep_ || step == inputrec_->init_step)
    {
        return;
    }

    // PME grid + cut-off optimization with GPUs or PME nodes
    // TODO pass SimulationWork object into this function, such that last argument can be set as
    // simulationWork.useGpuPmePpCommunication as is done in main MD loop.
    pme_loadbal_.addCycles((isVerbose_ && (dd_ == nullptr || DDMAIN(*dd_))) ? stderr : nullptr,
                           fr_,
                           statePropagatorData_->constBox(),
                           statePropagatorData_->constPositionsView().paddedArrayRef(),
                           wcycle_,
                           step,
                           step - inputrec_->init_step);
}

void PmeLoadBalanceHelper::teardown()
{
    pme_loadbal_.printSettings();
}

bool PmeLoadBalanceHelper::pmePrinting() const
{
    return pme_loadbal_.isPrintingLoad();
}

const PmeLoadBalancing& PmeLoadBalanceHelper::loadBalancingObject()
{
    return pme_loadbal_;
}

std::optional<SignallerCallback> PmeLoadBalanceHelper::registerNSCallback()
{
    return [this](Step step, Time gmx_unused time) { nextNSStep_ = step; };
}

} // namespace gmx
