/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Defines StopHandler, a helper class and two stop conditions.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "stophandler.h"

#include "config.h"

#include <algorithm>
#include <memory>
#include <utility>

#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/cstringutil.h"

namespace gmx
{

StopHandler::StopHandler(compat::not_null<SimulationSignal*>      signal,
                         bool                                     simulationShareState,
                         std::vector<std::function<StopSignal()>> stopConditions,
                         bool                                     neverUpdateNeighborList) :
    signal_(*signal),
    stopConditions_(std::move(stopConditions)),
    neverUpdateNeighborlist_(neverUpdateNeighborList)
{
    if (simulationShareState)
    {
        signal_.isLocal = false;
    }
}

StopConditionSignal::StopConditionSignal(int nstList, bool makeBinaryReproducibleSimulation, int nstSignalComm) :
    handledStopCondition_(StopCondition::None),
    makeBinaryReproducibleSimulation_(makeBinaryReproducibleSimulation),
    nstSignalComm_(nstSignalComm),
    nstList_(nstList)
{
}

StopSignal StopConditionSignal::getSignal(FILE* fplog)
{
    StopSignal signal = StopSignal::noSignal;

    /* Check whether everything is still alright */
    if (gmx_get_stop_condition() > handledStopCondition_)
    {
        int nsteps_stop = -1;

        /* this just makes signals[].sig compatible with the hack
           of sending signals around by MPI_Reduce together with
           other floats */
        if ((gmx_get_stop_condition() == StopCondition::NextNS)
            || (makeBinaryReproducibleSimulation_ && gmx_get_stop_condition() == StopCondition::Next))
        {
            /* We need at least two global communication steps to pass
             * around the signal. We stop at a pair-list creation step
             * to allow for exact continuation, when possible.
             */
            signal      = StopSignal::stopAtNextNSStep;
            nsteps_stop = std::max(nstList_, 2 * nstSignalComm_);
        }
        else if (gmx_get_stop_condition() == StopCondition::Next)
        {
            /* Stop directly after the next global communication step.
             * This breaks exact continuation.
             */
            signal      = StopSignal::stopImmediately;
            nsteps_stop = nstSignalComm_ + 1;
        }
        if (fplog)
        {
            fprintf(fplog,
                    "\n\nReceived the %s signal, stopping within %d steps\n\n",
                    gmx_get_signal_name(),
                    nsteps_stop);
            fflush(fplog);
        }
        fprintf(stderr,
                "\n\nReceived the %s signal, stopping within %d steps\n\n",
                gmx_get_signal_name(),
                nsteps_stop);
        fflush(stderr);
        handledStopCondition_ = gmx_get_stop_condition();
    }

    return signal;
}

StopConditionTime::StopConditionTime(int nstList, real maximumHoursToRun, int nstSignalComm) :
    signalSent_(false),
    maximumHoursToRun_(maximumHoursToRun),
    nstList_(nstList),
    nstSignalComm_(nstSignalComm),
    neverUpdateNeighborlist_(nstList <= 0)
{
}

StopSignal StopConditionTime::getSignal(bool bNS, int64_t step, FILE* fplog, gmx_walltime_accounting_t walltime_accounting)
{
    if (signalSent_)
    {
        // We only want to send it once, but might be called again before run is terminated
        return StopSignal::noSignal;
    }
    if ((bNS || neverUpdateNeighborlist_)
        && walltime_accounting_get_time_since_start(walltime_accounting)
                   > maximumHoursToRun_ * 60.0 * 60.0 * 0.99)
    {
        /* Signal to terminate the run */
        char sbuf[STEPSTRSIZE];
        int  nsteps_stop = std::max(nstList_, 2 * nstSignalComm_);
        if (fplog)
        {
            fprintf(fplog,
                    "\nStep %s: Run time exceeded %.3f hours, "
                    "will terminate the run within %d steps\n",
                    gmx_step_str(step, sbuf),
                    maximumHoursToRun_ * 0.99,
                    nsteps_stop);
        }
        fprintf(stderr,
                "\nStep %s: Run time exceeded %.3f hours, "
                "will terminate the run within %d steps\n",
                gmx_step_str(step, sbuf),
                maximumHoursToRun_ * 0.99,
                nsteps_stop);
        signalSent_ = true;
        return StopSignal::stopAtNextNSStep;
    }
    return StopSignal::noSignal;
}

void StopHandlerBuilder::registerStopCondition(std::function<StopSignal()> stopCondition)
{
    stopConditions_.emplace_back(std::move(stopCondition));
};

std::unique_ptr<StopHandler> StopHandlerBuilder::getStopHandlerMD(compat::not_null<SimulationSignal*> signal,
                                                                  bool simulationShareState,
                                                                  bool isMain,
                                                                  int  nstList,
                                                                  bool makeBinaryReproducibleSimulation,
                                                                  int   nstSignalComm,
                                                                  real  maximumHoursToRun,
                                                                  bool  neverUpdateNeighborList,
                                                                  FILE* fplog,
                                                                  const int64_t&  step,
                                                                  const gmx_bool& bNS,
                                                                  gmx_walltime_accounting_t walltime_accounting)
{
    if (!GMX_THREAD_MPI || isMain)
    {
        // Using shared ptr because move-only callable not supported by std::function.
        // Would require replacement such as fu2::function or cxx_function.
        auto stopConditionSignal = std::make_shared<StopConditionSignal>(
                nstList, makeBinaryReproducibleSimulation, nstSignalComm);
        registerStopCondition(
                [stopConditionSignal, fplog]() { return stopConditionSignal->getSignal(fplog); });
    }

    if (isMain && maximumHoursToRun > 0)
    {
        auto stopConditionTime =
                std::make_shared<StopConditionTime>(nstList, maximumHoursToRun, nstSignalComm);
        registerStopCondition([stopConditionTime, &bNS, &step, fplog, walltime_accounting]() {
            return stopConditionTime->getSignal(bNS, step, fplog, walltime_accounting);
        });
    }

    return std::make_unique<StopHandler>(
            signal, simulationShareState, stopConditions_, neverUpdateNeighborList);
}

} // namespace gmx
