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
 * \brief
 * Declares signal setter and handler classes.
 *
 * These classes encapsulate the setting and handling of signals.
 * Any class setting or handling signals needs to implement ISignalSetter
 * or ISignalHandler, respectively.
 *
 * Current ISignalSetter classes set the checkpoint signal (CheckpointSignalSetter),
 * the stop signal (StopSignalSetterCondition and StopSignalSetterTime), and
 * the reset signal (ResetSignalSetter). Current ISignalHandler classes handle the
 * respective signals (CheckpointSignalHandler, StopSignalHandler, ResetSignalSetter).
 *
 * The approach is ready for a task-based design: The signal setters and handlers are
 * only created if they are needed (e.g. only if counter resetting is needed, only if
 * we're on master node, etc). The setter and handlers bind to data they need to access
 * at run time during construction time via const references. This allows the task to
 * run later without any inputs.
 *
 * Currently, the setter and handler objects are accessed via pointers, which are set
 * to nullptr if no object is needed. This is handled in the (static) functions of the
 * SignalHelper class. In the future, this will be achieved by adding (or not adding)
 * handlers / setters to the task graph.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_SIGNALHANDLER
#define GMX_MDLIB_SIGNALHANDLER

#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/logger.h"

namespace gmx
{

/*! \libinternal
 * \brief Interface for classes that set signals
 *
 * Classes implementing this interface need to redefine the setSignal() function, which
 * takes no arguments. Binding to data needs therefore to happen at build time of the
 * objects.
 */
class ISignalSetter
{
    public:
        //! Sets a signal - needs to be redefined by implementations
        virtual void setSignal() = 0;
        //! Default destructor
        virtual ~ISignalSetter() = default;
};

/*! \libinternal
 * \brief Interface for classes that handle signals
 *
 * Classes implementing this interface need to redefine the handleSignal() function, which
 * takes no arguments. Binding to data needs therefore to happen at build time of the
 * objects.
 */
class ISignalHandler
{
    public:
        //! Handles a received signal - needs to be redefined by implementations
        virtual void handleSignal() = 0;
        //! Default destructor
        virtual ~ISignalHandler()   = default;
};

/*! \libinternal
 * \brief Class setting the checkpoint signal
 *
 * Master rank sets the checkpointing signal periodically
 */
class CheckpointSignalSetter : public ISignalSetter
{
    public:
        /*! \brief CheckpointSignalSetter constructor
         *
         * Needs a pointer to the signal to be set, information on whether multiple simulations
         * need to be synchronized, and (const) references to data it needs to determine the
         * setting of a signal. The latter allows the setSignal() routine to run without
         * additional arguments, making it easier to be ran in a task-based environment.
         */
        CheckpointSignalSetter(gmx::SimulationSignal    *sig,
                               bool                      needSync,
                               const t_commrec          *cr,
                               const MdrunOptions       &mdrunOptions,
                               const gmx_bool           &bGStat,
                               gmx_walltime_accounting_t walltime_accounting
                               );

        /*! \brief Decides whether a checkpointing signal needs to be set
         *
         * Checkpointing signal is set based on the elapsed run time and the checkpointing
         * interval.
         */
        void setSignal() override;

    private:
        SimulationSignal         *signal;
        int                       nchkpt = 1;
        const real                cpt_period;
        const bool                isParallel;
        const gmx_bool           &bGStat;
        gmx_walltime_accounting_t walltime_accounting;
};

/*! \libinternal
 * \brief Class handling the checkpoint signal
 *
 * All ranks receive checkpointing signal and set the respective flag
 */
class CheckpointSignalHandler : public ISignalHandler
{
    public:
        /*! \brief CheckpointSignalHandler constructor
         *
         * Needs a pointer to the signal which is set by CheckpointSignalSetter, and
         * (const) references to data it needs to determine whether a set signal needs
         * to be handled. The latter allows the setSignal() routine to run without
         * additional arguments, making it easier to be ran in a task-based environment.
         */
        CheckpointSignalHandler(gmx::SimulationSignal    *sig,
                                const t_inputrec         *ir,
                                const MdrunOptions       &mdrunOptions,
                                const gmx_bool           &bNS,
                                const gmx_bool           &bLastStep,
                                const int64_t            &step);

        /*! \brief Decides whether a checkpoint shall be written at this step
         *
         * Checkpointing is done if this is not the initial step, and
         *   * a signal has been set and the current step is a neighborlist creation
         *     step, or
         *   * the current step is the last step and a the simulation is writing
         *     configurations.
         */
        void handleSignal() override;

        //! Whether a checkpoint should be written in the current step
        bool doCheckpointThisStep();

    private:
        bool               checkpointThisStep;
        SimulationSignal  *signal;
        const int          nstlist;
        const int64_t      init_step;
        const bool         writeConfout;
        const gmx_bool    &bNS;
        const gmx_bool    &bLastStep;
        const int64_t     &step;
};

/*! \libinternal
 * \brief Class setting the stop signal based on gmx_get_stop_condition()
 *
 * Master rank sets the stop signal if required (generally due to SIGINT).
 */
class StopSignalSetterCondition : public ISignalSetter
{
    public:
        /*! \brief StopSignalSetterCondition constructor
         *
         * Needs a pointer to the signal to be set, information on whether multiple simulations
         * need to be synchronized, and (const) references to data it needs to determine the
         * setting of a signal. The latter allows the setSignal() routine to run without
         * additional arguments, making it easier to be ran in a task-based environment.
         */
        StopSignalSetterCondition(
            gmx::SimulationSignal     *sig,
            bool                       needSync,
            const t_inputrec          *ir,
            const MdrunOptions        &mdrunOptions,
            int                        nstSignalComm,
            FILE                      *fplog);

        /*! \brief Decides whether a stopping signal needs to be set
         *
         * Stop signal is set based on the value of gmx_get_stop_condition(): Set signal for
         * stop at the next neighbor-searching step at first SIGINT / SIGTERM, set signal
         * for stop at the next step at second SIGINT / SIGTERM.
         */
        void setSignal() override;

    private:
        SimulationSignal         *signal;
        int                       handled_stop_condition = gmx_stop_cond_none;
        const bool                reproducible;
        const int                 nstSignalComm;
        const int                 nstlist;
        FILE                     *fplog;

};

/*! \libinternal
 * \brief Class setting the stop signal based on maximal run time
 *
 * Master rank sets the stop signal if run time exceeds maximal run time.
 */
class StopSignalSetterTime : public ISignalSetter
{
    public:
        /*! \brief StopSignalSetterTime constructor
         *
         * Needs a pointer to the signal to be set, information on whether multiple simulations
         * need to be synchronized, and (const) references to data it needs to determine the
         * setting of a signal. The latter allows the setSignal() routine to run without
         * additional arguments, making it easier to be ran in a task-based environment.
         */
        StopSignalSetterTime(
            gmx::SimulationSignal     *sig,
            bool                       needSync,
            const t_inputrec          *ir,
            const MdrunOptions        &mdrunOptions,
            FILE                      *fplog,
            const gmx_bool            &bNS,
            const int64_t             &step,
            gmx_walltime_accounting_t  walltime_accounting);

        /*! \brief Decides whether a stopping signal needs to be set
         *
         * Stop signal is set if run time is greater than 99% of maximal run time. Signal will
         * trigger stopping of the simulation at the next neighbor-searching step.
         */
        void setSignal() override;

    private:
        SimulationSignal         *signal;
        const real                maximumHoursToRun;
        const int                 nstlist;
        const gmx_bool           &bNS;
        const int64_t            &step;
        FILE                     *fplog;
        gmx_walltime_accounting_t walltime_accounting;

};

/*! \libinternal
 * \brief Class handling the stop signal
 *
 * All ranks receive the stop signal and set the respective flag
 */
class StopSignalHandler : public ISignalHandler
{
    public:
        /*! \brief StopSignalHandler constructor
         *
         * Needs a pointer to the signal which is set by StopSignalSetter, and
         * (const) references to data it needs to determine whether a set signal needs
         * to be handled. The latter allows the setSignal() routine to run without
         * additional arguments, making it easier to be ran in a task-based environment.
         */
        StopSignalHandler(gmx::SimulationSignal    *sig,
                          const t_inputrec         *ir,
                          const gmx_bool           &bNS);

        /*! \brief Decides whether the simulation shall be stopped at the next step
         *
         * The simulation is stopped at the next step if
         *   * the signal for immediate stop was received, or
         *   * the signal for stop at the next neighbor-searching step was received, and
         *     the current step is a neighbor-searching step.
         */
        void handleSignal() override;

        //! Whether the simulation should be stopped in the current step
        bool stoppingThisStep();

    private:
        bool              lastStep = false;
        SimulationSignal *signal;
        const int         nstlist;
        const gmx_bool   &bNS;

};

/*! \libinternal
 * \brief Class setting the reset signal
 *
 * Master rank sets the reset signal if half the run time is reached.
 */
class ResetSignalSetter : public ISignalSetter
{
    public:
        /*! \brief ResetSignalSetter constructor
         *
         * Needs a pointer to the signal to be set, information on whether multiple simulations
         * need to be synchronized, and (const) references to data it needs to determine the
         * setting of a signal. The latter allows the setSignal() routine to run without
         * additional arguments, making it easier to be ran in a task-based environment.
         */
        ResetSignalSetter(gmx::SimulationSignal      *sig,
                          bool                        needSync,
                          const MdrunOptions         &mdrunOptions,
                          gmx_walltime_accounting_t   walltime_accounting);

        /*! \brief Decides whether a reset signal needs to be set
         *
         * Reset signal is set if run time is greater than 49.5% of maximal run time.
         */
        void setSignal() override;

    private:
        bool                      active = true;
        SimulationSignal         *signal;
        const real                maximumHoursToRun;
        gmx_walltime_accounting_t walltime_accounting;
};

/*! \libinternal
 * \brief Class handling the reset signal
 *
 * All ranks receive the reset signal and reset their respective counters.
 * This also resets the counters if half the time steps have passed (no communication needed).
 */
class ResetSignalHandler : public ISignalHandler
{
    public:
        /*! \brief ResetSignalHandler constructor
         *
         * Needs a pointer to the signal which is set by ResetSignalSetter, and
         * (const) references to data it needs to determine whether a set signal needs
         * to be handled. The latter allows the setSignal() routine to run without
         * additional arguments, making it easier to be ran in a task-based environment.
         */
        ResetSignalHandler(gmx::SimulationSignal      *sig,
                           const t_commrec            *cr,
                           const int64_t              &step,
                           const int64_t              &step_rel,
                           gmx_walltime_accounting_t   walltime_accounting,
                           gmx_wallcycle_t             wcycle,
                           const gmx::MDLogger        &mdlog,
                           FILE                       *fplog,
                           nonbonded_verlet_t         *nbv,
                           const gmx_pme_t            *pme,
                           t_nrnb                     *nrnb,
                           const pme_load_balancing_t *pme_loadbal);

        /*! \brief Decides whether the counters are reset, and performs the reset if needed
         *
         * The counters are reset if
         *   * the signal for resetting was received, or
         *   * the (local) number of steps reached the defined counter reset step.
         * Note that even if two reset conditions are present (at a specific step and a
         * spedific time), the reset will only take place once, whenever the first condition
         * is met.
         */
        void handleSignal() override;

    private:
        SimulationSignal                 *signal;
        const t_commrec* const            cr;
        const int64_t                    &step;
        const int64_t                    &step_rel;
        gmx_walltime_accounting_t         walltime_accounting;
        gmx_wallcycle_t                   wcycle;
        const gmx::MDLogger              &mdlog;
        FILE* const                       fplog;
        nonbonded_verlet_t* const         nbv;
        const gmx_pme_t* const            pme;
        t_nrnb* const                     nrnb;
        const pme_load_balancing_t* const pme_loadbal;
};

/*! \libinternal
 * \brief Class constructing the signal setter / handler pairs
 *
 * The functions of this class determine (based on the input parameters of the simulation, the current
 * rank, command line options, etc) whether signal setters and / or handlers are needed, and constructs
 * the required ones. The built signal objects are returned via unique pointers. If no signal objects
 * of a specific class is needed, the functions return nullptr instead. By design, the signal objects
 * need to bind to required data at construction time, requiring the makeXXX() functions in this class
 * to have access (via constant references) to this data.
 *
 * \todo In a later stage, the functionality of this class will likely be handled by the owner of
 * the task graph, adding only the required signal objects to the task graph.
 */
class SignalHelper
{
    public:
        //! Create checkpoint signal setter and / or handler, if needed. Return nullptr otherwise.
        static std::pair < std::unique_ptr<CheckpointSignalSetter>, std::unique_ptr < CheckpointSignalHandler>>
        makeCheckpointSignalClasses(gmx::SimulationSignal     *sig,
                                    const t_inputrec          *ir,
                                    const t_commrec           *cr,
                                    const MdrunOptions        &mdrunOptions,
                                    bool                       needSync,
                                    const gmx_bool            &bNS,
                                    const gmx_bool            &bLastStep,
                                    const int64_t             &step,
                                    const gmx_bool            &bGStat,
                                    gmx_walltime_accounting_t  walltime_accounting);

        //! Create stop signal setter(s) and / or handler, if needed. Return nullptr / empty vector otherwise.
        static std::pair < std::vector < std::unique_ptr < ISignalSetter>>, std::unique_ptr < StopSignalHandler>>
        makeStopSignalClasses(gmx::SimulationSignal        *sig,
                              bool                          needSync,
                              const t_inputrec             *ir,
                              const t_commrec              *cr,
                              const MdrunOptions           &mdrunOptions,
                              int                           nstSignalComm,
                              FILE                         *fplog,
                              const gmx_bool               &bNS,
                              const int64_t                &step,
                              gmx_walltime_accounting_t     walltime_accounting);

        //! Create reset signal setter and / or handler, if needed. Return nullptr otherwise.
        static std::pair < std::unique_ptr<ResetSignalSetter>, std::unique_ptr < ResetSignalHandler>>
        makeResetSignalClasses(gmx::SimulationSignal      *sig,
                               bool                        needSync,
                               const t_inputrec           *ir,
                               const t_commrec            *cr,
                               const MdrunOptions         &mdrunOptions,
                               gmx_walltime_accounting_t   walltime_accounting,
                               gmx_wallcycle_t             wcycle,
                               const gmx::MDLogger        &mdlog,
                               const int64_t              &step,
                               const int64_t              &step_rel,
                               FILE                       *fplog,
                               nonbonded_verlet_t         *nbv,
                               const gmx_pme_t            *pme,
                               t_nrnb                     *nrnb,
                               const pme_load_balancing_t *pme_loadbal);

};

}

#endif //GMX_MDLIB_SIGNALHANDLER
