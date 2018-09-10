/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \brief
 * Defines the checkpoint handler class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "checkpointhandler.h"

#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"

using namespace gmx;

CheckpointHandler::CheckpointHandler(
        const Integrator *integrator,
        MDState *mdState,
        const SimulationSetup *setup) :
    doSet_(false),
    doCheckpointing_(false),
    checkpointThisStep_(false),
    nNextCheckpoint_(1),
    checkedStep_(-1),
    isParallel_(PAR(integrator->cr)),
    writeConfout_(bool(integrator->mdrunOptions.writeConfout)),
    doNsEveryStep_(integrator->inputrec->nstlist == 0),
    initialStep_(integrator->inputrec->init_step),
    checkPointingPeriod_(integrator->mdrunOptions.checkpointOptions.period)
{
    if (setup->simulationsShareState)
    {
        mdState->signals.write()->at(eglsCHKPT).isLocal = false;
    }

    if (!(integrator->mdrunOptions.rerun || integrator->mdrunOptions.checkpointOptions.period < 0))
    {
        if (MASTER(integrator->cr))
        {
            doSet_ = true;
        }
        doCheckpointing_ = true;
    }
}

void CheckpointHandler::setSignalImpl_(
        const Integrator *integrator,
        MDState *mdState,
        const SimulationSetup *setup)
{
    /* In parallel we only have to check for checkpointing in steps
     * where we do global communication,
     *  otherwise the other nodes don't know.
     */
    auto signal = &mdState->signals.write()->at(eglsCHKPT);
    const double secondsSinceStart = walltime_accounting_get_time_since_start(integrator->walltime_accounting);
    if ((setup->bGStat || !isParallel_) &&
        (checkPointingPeriod_ == 0 || secondsSinceStart >= nNextCheckpoint_ * checkPointingPeriod_ * 60.0) &&
        signal->set == 0)
    {
        signal->sig = 1;
    }
}

void CheckpointHandler::doCheckpointImpl_(
        const Integrator *integrator,
        MDState *mdState,
        const SimulationSetup *setup)
{
    auto signal = &mdState->signals.write()->at(eglsCHKPT);
    checkpointThisStep_ = (((signal->set != 0 && (setup->bNS || doNsEveryStep_)) ||
                            (setup->bLastStep && writeConfout_)) &&
                           setup->step > initialStep_);
    if (checkpointThisStep_)
    {
        signal->set = 0;
        nNextCheckpoint_++;
    }
    checkedStep_ = setup->step;
}
