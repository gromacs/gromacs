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
/*! \libinternal \file
 * \brief Declares StopHandler, a helper class and two stop conditions.
 *
 * These classes encapsulate the setting and handling of stop signals.
 *
 * StopHandler lives during the lifetime of do_md. It checks via registered stop
 * conditions whether the simulation should be stopped at the next possible step or
 * at the next possible neighbor-searching step. It communicates this via signal to
 * all ranks and communicates this to do_md via stoppingAfterCurrentStep().
 *
 * StopHandlerBuilder is owned by the runner, and allows to register stop conditions
 * at a higher level, outside of do_md. Within do_md, it is creating a StopHandler
 * object by binding local data and passing a reference to the stop conditions.
 *
 * Here, we are implementing two stop conditions: StopConditionTime sets a stop condition
 * based on the elapsed time (only relevant if the -maxh flag was set), while
 * StopConditionSignal sets stop conditions via signals received from the operating
 * systems (SIGINT / SIGTERM).
 *
 * The stop conditions are stored as function pointers created by a lambda expression.
 * They bind to required local data, in the case of StopConditionTime and StopConditionSignal
 * these are partially owned by do_md. This requires these function pointers to be deleted
 * at the end of do_md(). This is achieved by having the do_md() specific function pointers
 * owned by StopHandler, which in turn is owned (via unique_ptr) by do_md().
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_STOPHANDLER_H
#define GMX_MDLIB_STOPHANDLER_H

#include <cstdint>
#include <cstdio>

#include <functional>
#include <memory>
#include <vector>

#include "gromacs/compat/pointers.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_walltime_accounting;

namespace gmx
{
/*! \brief Stop signals
 *
 * Signals that stop conditions can send to all ranks. Possible signals include
 *   * nothing to signal
 *   * stop at the next neighbor-searching step
 *   * stop as soon as signal is received
 */
enum class StopSignal : int
{
    noSignal         = 0,
    stopAtNextNSStep = 1,
    stopImmediately  = -1
};

/*! \brief Convert signed char (as used by SimulationSignal) to StopSignal enum
 *
 * * Expected values are
 *   \p sig ==  0 -- no signal
 *   \p sig >=  1 -- stop at next NS
 *   \p sig <= -1 -- stop asap
 */
static inline StopSignal convertToStopSignal(signed char sig)
{
    if (sig <= -1)
    {
        return StopSignal::stopImmediately;
    }
    else if (sig >= 1)
    {
        return StopSignal::stopAtNextNSStep;
    }
    else // sig == 0
    {
        return StopSignal::noSignal;
    }
}

/*! \libinternal
 * \brief Class handling the stop signal
 *
 * Loops over the registered stop conditions and sets a signal if
 * requested (currently only done by main rank).
 * All ranks receive the stop signal and set the respective flag.
 * The functions are implemented within this header file to avoid leaving
 * the translation unit unnecessarily.
 */
class StopHandler final
{
public:
    /*! \brief StopHandler constructor (will be called by StopHandlerBuilder)
     *
     * @param signal Non-null pointer to a signal used for reading and writing of signals
     * @param simulationShareState Whether this signal needs to be shared across multiple simulations
     * @param stopConditions Vector of callback functions setting the signal
     * @param neverUpdateNeighborList Whether simulation keeps same neighbor list forever
     *
     * Note: As the StopHandler does not work without this signal, it keeps a non-const reference
     * to it as a member variable.
     */
    StopHandler(compat::not_null<SimulationSignal*>      signal,
                bool                                     simulationShareState,
                std::vector<std::function<StopSignal()>> stopConditions,
                bool                                     neverUpdateNeighborList);

    /*! \brief Decides whether a stop signal shall be sent
     *
     * Loops over the stopCondition vector passed at build time (consisting of conditions
     * registered with StopHandlerBuilder, and conditions built by StopHandlerBuilder by
     * default), and sets any signal obtained.
     * Returns as soon as a StopSignal::stopImmediately signal was obtained, or after
     * checking all registered stop conditions.
     */
    void setSignal() const
    {
        for (const auto& condition : stopConditions_)
        {
            const StopSignal sig = condition();
            if (sig != StopSignal::noSignal)
            {
                signal_.sig = static_cast<signed char>(sig);
                if (sig == StopSignal::stopImmediately)
                {
                    // We don't want this to be overwritten by a less urgent stop
                    break;
                }
            }
        }
    }

    /*! \brief Decides whether the simulation shall be stopped after the current step
     *
     * The simulation is stopped after the current step if
     *   * the signal for immediate stop was received, or
     *   * the signal for stop at the next neighbor-searching step was received, and
     *     the current step is a neighbor-searching step.
     */
    bool stoppingAfterCurrentStep(bool bNS) const
    {
        return convertToStopSignal(signal_.set) == StopSignal::stopImmediately
               || (convertToStopSignal(signal_.set) == StopSignal::stopAtNextNSStep
                   && (bNS || neverUpdateNeighborlist_));
    }

private:
    SimulationSignal&                              signal_;
    const std::vector<std::function<StopSignal()>> stopConditions_;
    const bool                                     neverUpdateNeighborlist_;
};

/*! \libinternal
 * \brief Class setting the stop signal based on gmx_get_stop_condition()
 *
 * Main rank sets the stop signal if required (generally due to SIGINT).
 */
class StopConditionSignal final
{
public:
    /*! \brief StopConditionSignal constructor
     */
    StopConditionSignal(int nstList, bool makeBinaryReproducibleSimulation, int nstSignalComm);

    /*! \brief Decides whether a stopping signal needs to be set
     *
     * Stop signal is set based on the value of gmx_get_stop_condition(): Set signal for
     * stop at the next neighbor-searching step at first SIGINT / SIGTERM, set signal
     * for stop at the next step at second SIGINT / SIGTERM.
     */
    StopSignal getSignal(FILE* fplog);

private:
    StopCondition handledStopCondition_;
    const bool    makeBinaryReproducibleSimulation_;
    const int     nstSignalComm_;
    const int     nstList_;
};

/*! \libinternal
 * \brief Class setting the stop signal based on maximal run time
 *
 * Main rank sets the stop signal if run time exceeds maximal run time.
 */
class StopConditionTime final
{
public:
    /*! \brief StopConditionTime constructor
     */
    StopConditionTime(int nstList, real maximumHoursToRun, int nstSignalComm);

    /*! \brief Decides whether a stopping signal needs to be set
     *
     * Stop signal is set if run time is greater than 99% of maximal run time. Signal will
     * trigger stopping of the simulation at the next neighbor-searching step.
     */
    StopSignal getSignal(bool bNS, int64_t step, FILE* fplog, gmx_walltime_accounting* walltime_accounting);

private:
    bool signalSent_;

    const real maximumHoursToRun_;
    const int  nstList_;
    const int  nstSignalComm_;
    const bool neverUpdateNeighborlist_;
};

/*! \libinternal
 * \brief Class preparing the creation of a StopHandler
 *
 * An object of this helper class (owned by the runner) allows to register stop conditions
 * outside of the actual simulation run via registerStopCondition(), accepting a std::function
 * object. It then builds a StopHandler object inside do_md, once it can bind to required
 * local data.
 *
 * The registered stop conditions plus the standard MD stop conditions (stop based on
 * received signal from OS [SIGINT / SIGTERM] or based on maximal run time) are then called
 * by the Stophandler every step to determine whether the simulation should be stopped.
 * The registered functions need to be of type `std::function<StopSignal()>`, i.e. not
 * taking any input arguments and returning a `StopSignal` signal which will get propagated
 * to all ranks. If the function needs input arguments, these need to be bound (e.g. via
 * lambda capturing) before being registered with the StopHandlerBuilder.
 */
class StopHandlerBuilder final
{
public:
    /*! \brief Register stop condition
     *
     * This allows code in the scope of the StopHandlerBuilder (runner level) to inject
     * stop conditions in simulations. Stop conditions are defined as argument-less functions
     * which return a StopSignal. The return value of this function is then propagated to all
     * ranks, and allows to stop the simulation at the next global communication step (returned
     * signal StopSignal::stopImmediately), or at the next NS step (returned signal
     * StopSignal::stopAtNextNSStep, allows for exact continuation).
     *
     * Arguments needed by the stop condition function need to be bound / captured. If these
     * arguments are captured by reference or using a pointer, it is the registrant's
     * responsibility to ensure that these arguments do not go out of scope during the lifetime
     * of the StopHandlerBuilder.
     */
    void registerStopCondition(std::function<StopSignal()> stopCondition);

    /*! \brief Create StopHandler
     *
     * Gets called in the scope of the integrator (aka do_md()) to get a pointer to the
     * StopHandler for the current simulations. Adds the standard MD stop conditions
     * (e.g. gmx::StopConditionTime, gmx::StopConditionSignal) to the currently registered
     * stop conditions. Initializes a new StopHandler with this extended vector of
     * stop conditions. It is the caller's responsibility to make sure arguments passed by
     * pointer or reference remain valid for the lifetime of the returned StopHandler.
     */
    std::unique_ptr<StopHandler> getStopHandlerMD(compat::not_null<SimulationSignal*> signal,
                                                  bool            simulationShareState,
                                                  bool            isMain,
                                                  int             nstList,
                                                  bool            makeBinaryReproducibleSimulation,
                                                  int             nstSignalComm,
                                                  real            maximumHoursToRun,
                                                  bool            neverUpdateNeighborList,
                                                  FILE*           fplog,
                                                  const int64_t&  step,
                                                  const gmx_bool& bNS,
                                                  gmx_walltime_accounting* walltime_accounting);

private:
    /*! \brief Initial stopConditions
     *
     * StopConditions registered via registerStopCondition(). getStopHandlerMD will
     * copy this vector and add additional conditions before passing the new vector
     * to the built StopHandler object.
     */
    std::vector<std::function<StopSignal()>> stopConditions_;
};

} // namespace gmx

#endif // GMX_MDLIB_STOPHANDLER_H
