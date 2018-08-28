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
        gmx::SimulationSignal     *sig,
        bool                       needSync,
        const t_inputrec          *ir,
        const t_commrec           *cr,
        const MdrunOptions        &mdrunOptions,
        const gmx_bool            &bNS,
        const gmx_bool            &bLastStep,
        const int64_t             &step,
        const gmx_bool            &bGStat,
        gmx_walltime_accounting_t  walltime_accounting) :
    signal(sig),
    isParallel(PAR(cr)),
    writeConfout(bool(mdrunOptions.writeConfout)),
    nstlist(ir->nstlist),
    init_step(ir->init_step),
    cpt_period(mdrunOptions.checkpointOptions.period),
    bGStat(bGStat),
    bNS(bNS),
    bLastStep(bLastStep),
    step(step),
    walltime_accounting(walltime_accounting)

{
    if (needSync)
    {
        signal->isLocal = false;
    }

    if (!(mdrunOptions.rerun || mdrunOptions.checkpointOptions.period < 0))
    {
        if (MASTER(cr))
        {
            doSet = true;
        }
        doHandle = true;
    }
}

void CheckpointHandler::setSignalImpl()
{
    /* In parallel we only have to check for checkpointing in steps
     * where we do global communication,
     *  otherwise the other nodes don't know.
     */
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
    if ((bGStat || !isParallel) &&
        (cpt_period == 0 || secondsSinceStart >= nchkpt * cpt_period * 60.0) &&
        signal->set == 0)
    {
        signal->sig = 1;
        nchkpt++;
    }
}

void CheckpointHandler::handleSignalImpl()
{
    /* We write a checkpoint at this MD step when:
     * either at an NS step when we signalled through gs,
     * or at the last step (but not when we do not want confout),
     * but never at the first step or with rerun.
     */
    checkpointThisStep = (((signal->set && (bNS || nstlist == 0)) ||
                           (bLastStep && writeConfout)) &&
                          step > init_step);
    if (checkpointThisStep)
    {
        signal->set = 0;
    }
}
