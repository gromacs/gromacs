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
 * Defines the stop helper class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "stophandler.h"

#include "config.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/cstringutil.h"

using namespace gmx;

StopHandler::StopHandler(
        const std::vector < std::unique_ptr < gmx::IStopCondition>> &stopConditions,
        gmx::SimulationSignal                                       *sig,
        const t_inputrec                                            *ir,
        const gmx_bool                                              &bNS,
        bool                                                         needSync) :
    stopConditions(stopConditions),
    nstlist(ir->nstlist),
    signal(sig),
    bNS(bNS)
{
    if (needSync)
    {
        signal->isLocal = false;
    }
}

StopConditionSignal::StopConditionSignal(
        const t_inputrec* const   ir,
        const MdrunOptions       &mdrunOptions,
        int                       nstSignalComm,
        FILE                     *fplog) :
    reproducible(bool(mdrunOptions.reproducible)),
    nstSignalComm(nstSignalComm),
    nstlist(ir->nstlist),
    fplog(fplog)
{}

signed char StopConditionSignal::getSignal()
{
    signed char signal = 0;

    /* Check whether everything is still alright */
    if (static_cast<int>(gmx_get_stop_condition()) > handled_stop_condition)
    {
        int nsteps_stop   = -1;

        /* this just makes signals[].sig compatible with the hack
           of sending signals around by MPI_Reduce together with
           other floats */
        if ((gmx_get_stop_condition() == gmx_stop_cond_next_ns) ||
            (reproducible && gmx_get_stop_condition() == gmx_stop_cond_next))
        {
            /* We need at least two global communication steps to pass
             * around the signal. We stop at a pair-list creation step
             * to allow for exact continuation, when possible.
             */
            signal      = 1;
            nsteps_stop = std::max(nstlist, 2 * nstSignalComm);
        }
        else if (gmx_get_stop_condition() == gmx_stop_cond_next)
        {
            /* Stop directly after the next global communication step.
             * This breaks exact continuation.
             */
            signal      = -1;
            nsteps_stop = nstSignalComm + 1;
        }
        if (fplog)
        {
            fprintf(fplog,
                    "\n\nReceived the %s signal, stopping within %d steps\n\n",
                    gmx_get_signal_name(), nsteps_stop);
            fflush(fplog);
        }
        fprintf(stderr,
                "\n\nReceived the %s signal, stopping within %d steps\n\n",
                gmx_get_signal_name(), nsteps_stop);
        fflush(stderr);
        handled_stop_condition = static_cast<int>(gmx_get_stop_condition());
    }

    return signal;
}

StopConditionTime::StopConditionTime(
        const t_inputrec* const    ir,
        const MdrunOptions        &mdrunOptions,
        FILE                      *fplog,
        const gmx_bool            &bNS,
        const int64_t             &step,
        int                        nstSignalComm,
        gmx_walltime_accounting_t  walltime_accounting) :
    maximumHoursToRun(mdrunOptions.maximumHoursToRun),
    nstlist(ir->nstlist),
    nstSignalComm(nstSignalComm),
    bNS(bNS),
    step(step),
    fplog(fplog),
    walltime_accounting(walltime_accounting)
{}

signed char StopConditionTime::getSignal()
{
    if (signalSent)
    {
        // We only want to send it once, but might be called again before run is terminated
        return 0;
    }
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
    if ((bNS || nstlist <= 0) &&
        secondsSinceStart > maximumHoursToRun * 60.0 * 60.0 * 0.99)
    {
        /* Signal to terminate the run */
        char sbuf[STEPSTRSIZE];
        int  nsteps_stop = std::max(nstlist, 2 * nstSignalComm);
        if (fplog)
        {
            fprintf(fplog, "\nStep %s: Run time exceeded %.3f hours, "
                    "will terminate the run within %d steps\n",
                    gmx_step_str(step, sbuf), maximumHoursToRun * 0.99, nsteps_stop);
        }
        fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, "
                "will terminate the run within %d steps\n",
                gmx_step_str(step, sbuf), maximumHoursToRun * 0.99, nsteps_stop);
        signalSent = true;
        return 1;
    }
    return 0;
}

void StopHandlerHelper::registerStopCondition(std::unique_ptr<IStopCondition> stopCondition)
{
    stopConditions.emplace_back(std::move(stopCondition));
};

StopHandler StopHandlerHelper::getStopHandlerMD(
        gmx::SimulationSignal    *sig,
        bool                      needSync,
        const t_inputrec         *ir,
        const t_commrec          *cr,
        const MdrunOptions       &mdrunOptions,
        int                       nstSignalComm,
        FILE                     *fplog,
        const int64_t            &step,
        const gmx_bool           &bNS,
        gmx_walltime_accounting_t walltime_accounting)
{
#if GMX_THREAD_MPI
    if (MASTER(cr))
    {
        registerStopCondition(compat::make_unique<StopConditionSignal>(
                                      ir, mdrunOptions, nstSignalComm, fplog));
    }
#else
    registerStopCondition(compat::make_unique<StopConditionSignal>(
                                  ir, mdrunOptions, nstSignalComm, fplog));
#endif

    if (MASTER(cr) && mdrunOptions.maximumHoursToRun > 0)
    {
        registerStopCondition(compat::make_unique<StopConditionTime>(
                                      ir, mdrunOptions, fplog, bNS, step, nstSignalComm, walltime_accounting));
    }

    return {stopConditions, sig, ir, bNS, needSync};
}
