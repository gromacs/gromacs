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
        virtual void setSignal() = 0;
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
        virtual void handleSignal() = 0;
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
        CheckpointSignalSetter(gmx::SimulationSignal    *sig,
                               bool                      needSync,
                               const t_commrec          *cr,
                               const MdrunOptions       &mdrunOptions,
                               const gmx_bool           &bGStat,
                               gmx_walltime_accounting_t walltime_accounting
                               );
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
        CheckpointSignalHandler(gmx::SimulationSignal    *sig,
                                const t_inputrec         *ir,
                                const MdrunOptions       &mdrunOptions,
                                const gmx_bool           &bNS,
                                const gmx_bool           &bLastStep,
                                const int64_t            &step);

        void handleSignal() override;

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
        StopSignalSetterCondition(
            gmx::SimulationSignal     *sig,
            bool                       needSync,
            const t_inputrec          *ir,
            const MdrunOptions        &mdrunOptions,
            int                        nstSignalComm,
            FILE                      *fplog);

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
        StopSignalSetterTime(
            gmx::SimulationSignal     *sig,
            bool                       needSync,
            const t_inputrec          *ir,
            const MdrunOptions        &mdrunOptions,
            FILE                      *fplog,
            const gmx_bool            &bNS,
            const int64_t             &step,
            gmx_walltime_accounting_t  walltime_accounting);

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
        StopSignalHandler(gmx::SimulationSignal    *sig,
                          const t_inputrec         *ir,
                          const gmx_bool           &bNS);

        void handleSignal() override;

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
        ResetSignalSetter(gmx::SimulationSignal      *sig,
                          bool                        needSync,
                          const MdrunOptions         &mdrunOptions,
                          gmx_walltime_accounting_t   walltime_accounting);

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
