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
 * \brief Defines StopHandler, a helper class and two stop conditions.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "stophandler.h"

#include "config.h"

#include <algorithm>

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/cstringutil.h"

using namespace gmx;

StopHandler::StopHandler(
        const std::vector < std::shared_ptr < std::function<signed char()>>> &stopConditions,
        gmx::SimulationSignal                                                *signal,
        const t_inputrec                                                     *ir,
        bool                                                                  needSync,
        std::function<void()>                                                 destructorCallback) :
    stopConditions_(stopConditions),
    doNsEveryStep_(ir->nstlist == 0),
    destructorCallback_(std::move(destructorCallback))
{
    if (needSync)
    {
        signal->isLocal = false;
    }
}

StopConditionSignal::StopConditionSignal(
        const t_inputrec* const ir,
        const MdrunOptions     &mdrunOptions,
        int                     nstSignalComm) :
    handled_stop_condition(gmx_stop_cond_none),
    reproducible(bool(mdrunOptions.reproducible)),
    nstSignalComm(nstSignalComm),
    nstlist(ir->nstlist)
{}

signed char StopConditionSignal::getSignal(FILE *fplog)
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
        const t_inputrec*   ir,
        const MdrunOptions &mdrunOptions,
        int                 nstSignalComm) :
    maximumHoursToRun(mdrunOptions.maximumHoursToRun),
    nstlist(ir->nstlist),
    nstSignalComm(nstSignalComm)
{}

signed char StopConditionTime::getSignal(
        bool                      bNS,
        int64_t                   step,
        FILE                     *fplog,
        gmx_walltime_accounting_t walltime_accounting)
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

void StopHandlerHelper::registerStopCondition(std::shared_ptr < std::function < signed char()>> stopCondition)
{
    stopConditions_.emplace_back(stopCondition);
};

std::unique_ptr<StopHandler> StopHandlerHelper::getStopHandlerMD(
        gmx::SimulationSignal    *signal,
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
        stopConditionSignal_ = compat::make_unique<StopConditionSignal>(
                    ir, mdrunOptions, nstSignalComm);
        stopConditionSignalFunction_ = std::make_shared < std::function < signed char()>>(
                [this, fplog](){return this->stopConditionSignal_->getSignal(fplog); });
        registerStopCondition(stopConditionSignalFunction_);
    }
#else
    stopConditionSignal_ = compat::make_unique<StopConditionSignal>(
                ir, mdrunOptions, nstSignalComm);
    stopConditionSignalFunction_ = std::make_shared < std::function < signed char()>>(
            [this, fplog](){return this->stopConditionSignal_->getSignal(fplog); });
    registerStopCondition(stopConditionSignalFunction_);
#endif

    if (MASTER(cr) && mdrunOptions.maximumHoursToRun > 0)
    {
        stopConditionTime_ = compat::make_unique<StopConditionTime>(
                    ir, mdrunOptions, nstSignalComm);
        stopConditionTimeFunction_ = std::make_shared < std::function < signed char()>>(
                [this, bNS, &step, fplog, walltime_accounting]()
                {return this->stopConditionTime_->getSignal(bNS, step, fplog, walltime_accounting); });
        registerStopCondition(stopConditionTimeFunction_);

    }

    return compat::make_unique<StopHandler>(
            stopConditions_, signal, ir, needSync,
            [this](){this->deleteStopHandlerMD();});
}

void StopHandlerHelper::deleteStopHandlerMD()
{
    stopConditions_.erase(
            std::remove(stopConditions_.begin(), stopConditions_.end(), stopConditionSignalFunction_),
            stopConditions_.end());
    stopConditionSignalFunction_.reset();
    stopConditionSignal_.reset();
    stopConditions_.erase(
            std::remove(stopConditions_.begin(), stopConditions_.end(), stopConditionTimeFunction_),
            stopConditions_.end());
    stopConditionTimeFunction_.reset();
    stopConditionTime_.reset();
}
