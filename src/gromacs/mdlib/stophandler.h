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
/*! \libinternal \file
 * \brief Declares StopHandler, a helper class and two stop conditions.
 *
 * These classes encapsulate the setting and handling of stop signals.
 *
 * StopHandler lives during the lifetime of do_md. It checks via registered stop
 * conditions whether the simulation should be stopped at the next possible step or
 * at the next possible neighbor-searching step. It communicates this via signal to
 * all ranks and communicates this to do_md via stoppingThisStep().
 *
 * StopHandlerHelper is owned by the runner, and allows to register stop conditions
 * at a higher level, outside of do_md. Within do_md, it is creating a StopHandler
 * object by binding local data and passing a reference to the stop conditions.
 *
 * Here, we are implementing two stop conditions: StopConditionTime sets a stop condition
 * based on the elapsed time (only relevant if the -maxh flag was set), while
 * StopConditionSignal sets stop conditions via signals received from the operating
 * systems (SIGINT / SIGTERM).
 *
 * The stop conditions are stored as function pointers created by a lambda expression.
 * They bind to required local data, the case of StopConditionTime and StopConditionSignal
 * these are partially owned by do_md. This requires the StopHandlerHelper to delete them
 * when do_md() terminates, which is handled via a call-back function that StopHandler
 * calls at destruction time.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_STOPHANDLER_H
#define GMX_MDLIB_STOPHANDLER_H

#include <functional>
#include <memory>
#include <vector>

#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/timing/walltime_accounting.h"

namespace gmx
{
/*! \libinternal
 * \brief Class handling the stop signal
 *
 * Loops over the registered stop conditions and sets a signal if
 * requested (currently only done by master rank).
 * All ranks receive the stop signal and set the respective flag.
 * The functions are implemented within this header file to avoid leaving
 * the translation unit unnecessarily.
 */
class StopHandler
{
    public:
        /*! \brief StopHandler constructor
         */
        StopHandler(
            const std::vector < std::shared_ptr < std::function<signed char()>>> &stopConditions,
            SimulationSignal                                                     *sig,
            const t_inputrec                                                     *ir,
            bool                                                                  needSync,
            std::function<void()>                                                 destructorCallback);

        /*! \brief StopHandler destructor
         *
         * When StopHandler gets out of scope, signal this back to StopHandlerHelper. Currently,
         * StopHandlerHelper builds stop conditions which are binding data local to do_md().
         * The scope of StopHandlerHelper is, however, outside of do_md(), and reusing the same
         * stop conditions outside of do_md(), or on a consecutive do_md() run, would fail
         * badly.
         *
         * TODO: Remove this when new data model is in place - if we have a class representing
         * the current MD state, this can be passed to all stop conditions, allowing for
         * a single function signature and getting rid of the need of binding to data.
         */
        ~StopHandler()
        {
            destructorCallback_();
        }

        /*! \brief StopHandler copy constructor
         *
         * Need to define this because we have a custom destructor. Copying this is dangerous,
         * since we are signalling to StopHandlerHelper to delete certain stop conditions when
         * destructing this.
         */
        StopHandler(const StopHandler&) = delete;

        /*! \brief Decides whether a stop signal shall be sent
         *
         * Loops over the stopConditions stored earlier, sets any signal obtained.
         * Returns as soon as a signal < 0 was obtained, or after checking all signal
         * setters.
         */
        void setSignal(SimulationSignal *signal)
        {
            for (auto &condition : stopConditions_)
            {
                const signed char sig = (*condition)();
                if (sig != 0)
                {
                    signal->sig = sig;
                    if (sig < 0)
                    {
                        // < 0 means immediate stop - we don't want this to be overwritten
                        // by a less urgent stop
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
        bool stoppingAfterCurrentStep(SimulationSignal *signal, bool bNS)
        {
            /* < 0 means stop at next step, > 0 means stop at next NS step */
            return signal->set < 0 || (signal->set > 0 && ( bNS || doNsEveryStep_));
        }

    private:
        const std::vector < std::shared_ptr < std::function<signed char()>>> &stopConditions_;
        const bool                                                            doNsEveryStep_;

        std::function<void()> destructorCallback_;
};

/*! \libinternal
 * \brief Class setting the stop signal based on gmx_get_stop_condition()
 *
 * Master rank sets the stop signal if required (generally due to SIGINT).
 */
class StopConditionSignal
{
    public:
        /*! \brief StopConditionSignal constructor
         */
        StopConditionSignal(
            const t_inputrec          *ir,
            const MdrunOptions        &mdrunOptions,
            int                        nstSignalComm);

        /*! \brief Decides whether a stopping signal needs to be set
         *
         * Stop signal is set based on the value of gmx_get_stop_condition(): Set signal for
         * stop at the next neighbor-searching step at first SIGINT / SIGTERM, set signal
         * for stop at the next step at second SIGINT / SIGTERM.
         */
        signed char getSignal(FILE *fplog);

    private:
        int        handled_stop_condition;
        const bool reproducible;
        const int  nstSignalComm;
        const int  nstlist;

};

/*! \libinternal
 * \brief Class setting the stop signal based on maximal run time
 *
 * Master rank sets the stop signal if run time exceeds maximal run time.
 */
class StopConditionTime
{
    public:
        /*! \brief StopSignalSetterTime constructor
         */
        StopConditionTime(
            const t_inputrec          *ir,
            const MdrunOptions        &mdrunOptions,
            int                        nstSignalComm);

        /*! \brief Decides whether a stopping signal needs to be set
         *
         * Stop signal is set if run time is greater than 99% of maximal run time. Signal will
         * trigger stopping of the simulation at the next neighbor-searching step.
         */
        signed char getSignal(
            bool                      bNS,
            int64_t                   step,
            FILE                     *fplog,
            gmx_walltime_accounting_t walltime_accounting);

    private:
        bool signalSent = false;

        const real                maximumHoursToRun;
        const int                 nstlist;
        const int                 nstSignalComm;
};

/*! \libinternal
 * \brief Class preparing the creation of a StopHandler
 *
 * An object of this helper class (owned by the runner) allows to register stop conditions
 * outside of the actual simulation run via registerStopCondition(). It then builds a
 * StopHandler object inside do_md, once it can bind to required local data.
 *
 * It needs to delete the stop conditions which were bound to local do_md data to avoid
 * failures when calling functions bound to data which got out of scope. This is done by
 * deleteStopHandlerMD().
 */
class StopHandlerHelper
{
    public:
        //! Add stop condition
        void registerStopCondition(std::shared_ptr < std::function < signed char()>> stopCondition);

        //! Get stop conditions
        const std::vector < std::shared_ptr < std::function<signed char()>>> &getStopConditions() const
        {
            return stopConditions_;
        }

        //! Create StopHandler
        std::unique_ptr<StopHandler> getStopHandlerMD(
            gmx::SimulationSignal    *signal,
            bool                      needSync,
            const t_inputrec         *ir,
            const t_commrec          *cr,
            const MdrunOptions       &mdrunOptions,
            int                       nstSignalComm,
            FILE                     *fplog,
            const int64_t            &step,
            const gmx_bool           &bNS,
            gmx_walltime_accounting_t walltime_accounting);

        //! Delete MD-specific stop conditions
        void deleteStopHandlerMD();

    private:
        std::vector < std::shared_ptr < std::function<signed char()>>> stopConditions_;
        std::unique_ptr<StopConditionSignal>                           stopConditionSignal_;
        std::shared_ptr < std::function < signed char()>>              stopConditionSignalFunction_;
        std::unique_ptr<StopConditionTime>                             stopConditionTime_;
        std::shared_ptr < std::function < signed char()>>              stopConditionTimeFunction_;
};

}      // namespace gmx

#endif //GMX_MDLIB_STOPHANDLER_H
