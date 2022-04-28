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

#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/mdtypes/commrec.h"
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
                                              const t_inputrec*         inputrec,
                                              const t_forcerec*         fr,
                                              const SimulationWorkload& simWorkload)
{
    return (mdrunOptions.tunePme && usingPme(fr->ic->eeltype) && !mdrunOptions.reproducible
            && inputrec->cutoff_scheme != CutoffScheme::Group && !simWorkload.useGpuPmeDecomposition);
}

PmeLoadBalanceHelper::PmeLoadBalanceHelper(bool                 isVerbose,
                                           StatePropagatorData* statePropagatorData,
                                           FILE*                fplog,
                                           t_commrec*           cr,
                                           const MDLogger&      mdlog,
                                           const t_inputrec*    inputrec,
                                           gmx_wallcycle*       wcycle,
                                           t_forcerec*          fr) :
    pme_loadbal_(nullptr),
    nextNSStep_(-1),
    isVerbose_(isVerbose),
    bPMETunePrinting_(false),
    statePropagatorData_(statePropagatorData),
    fplog_(fplog),
    cr_(cr),
    mdlog_(mdlog),
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
    pme_loadbal_init(
            &pme_loadbal_, cr_, mdlog_, *inputrec_, box, *fr_->ic, *fr_->nbv, fr_->pmedata, fr_->nbv->useGpu());
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
    pme_loadbal_do(pme_loadbal_,
                   cr_,
                   (isVerbose_ && MASTER(cr_)) ? stderr : nullptr,
                   fplog_,
                   mdlog_,
                   *inputrec_,
                   fr_,
                   statePropagatorData_->constBox(),
                   statePropagatorData_->constPositionsView().paddedArrayRef(),
                   wcycle_,
                   step,
                   step - inputrec_->init_step,
                   &bPMETunePrinting_,
                   false);
}

void PmeLoadBalanceHelper::teardown()
{
    pme_loadbal_done(pme_loadbal_, fplog_, mdlog_, fr_->nbv->useGpu());
}

bool PmeLoadBalanceHelper::pmePrinting() const
{
    return bPMETunePrinting_;
}

const pme_load_balancing_t* PmeLoadBalanceHelper::loadBalancingObject()
{
    return pme_loadbal_;
}

std::optional<SignallerCallback> PmeLoadBalanceHelper::registerNSCallback()
{
    return [this](Step step, Time gmx_unused time) { nextNSStep_ = step; };
}

} // namespace gmx
