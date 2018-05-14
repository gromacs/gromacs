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

#ifndef GMX_MDLIB_SIGNALHANDLER
#define GMX_MDLIB_SIGNALHANDLER


#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/logger.h"
#include "simulationsignal.h"
#include "mdrun.h"
#include "sighandler.h"

namespace gmx
{

class ISignalSetter
{
    public:
        virtual void setSignal() = 0;
        virtual ~ISignalSetter() = default;
};

class ISignalHandler
{
    public:
        virtual void handleSignal() = 0;
        virtual ~ISignalHandler()   = default;
};

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

class CheckpointSignalHandler : public ISignalHandler
{
    public:
        CheckpointSignalHandler(gmx::SimulationSignal    *sig,
                                const t_inputrec         *ir,
                                const MdrunOptions       &mdrunOptions,
                                const gmx_bool           &bNS,
                                const gmx_bool           &bLastStep,
                                const gmx_int64_t        &step);

        void handleSignal() override;

        bool doCheckpointThisStep();

    private:
        bool               checkpointThisStep;
        SimulationSignal  *signal;
        const int          nstlist;
        const gmx_int64_t  init_step;
        const bool         writeConfout;
        const gmx_bool    &bNS;
        const gmx_bool    &bLastStep;
        const gmx_int64_t &step;
};

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
            const gmx_int64_t         &step,
            gmx_walltime_accounting_t  walltime_accounting);

        void setSignal() override;

    private:
        SimulationSignal         *signal;
        const real                maximumHoursToRun;
        const int                 nstlist;
        const gmx_bool           &bNS;
        const gmx_int64_t        &step;
        FILE                     *fplog;
        gmx_walltime_accounting_t walltime_accounting;

};

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

class ResetSignalHandler : public ISignalHandler
{
    public:
        ResetSignalHandler(gmx::SimulationSignal      *sig,
                           const t_commrec            *cr,
                           const gmx_int64_t          &step,
                           const gmx_int64_t          &step_rel,
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
        const gmx_int64_t                &step;
        const gmx_int64_t                &step_rel;
        gmx_walltime_accounting_t         walltime_accounting;
        gmx_wallcycle_t                   wcycle;
        const gmx::MDLogger              &mdlog;
        FILE* const                       fplog;
        nonbonded_verlet_t* const         nbv;
        const gmx_pme_t* const            pme;
        t_nrnb* const                     nrnb;
        const pme_load_balancing_t* const pme_loadbal;
};

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
                                    const gmx_int64_t         &step,
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
                              const gmx_int64_t            &step,
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
                               const gmx_int64_t          &step,
                               const gmx_int64_t          &step_rel,
                               FILE                       *fplog,
                               nonbonded_verlet_t         *nbv,
                               const gmx_pme_t            *pme,
                               t_nrnb                     *nrnb,
                               const pme_load_balancing_t *pme_loadbal);

};

}

#endif //GMX_MDLIB_SIGNALHANDLER
