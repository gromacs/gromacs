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

#ifndef GROMACS_SIGNALHANDLERS_H
#define GROMACS_SIGNALHANDLERS_H


#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/logger.h"
#include "simulationsignal.h"
#include "mdrun.h"
#include "sighandler.h"

namespace gmx
{

class SignalHandler
{
    public:
        explicit SignalHandler(gmx::SimulationSignal* sig) :
            signal(sig) {}

        virtual ~SignalHandler() = default;

        void need_sync()
        {
            // sync = true;
            signal->isLocal = false;
        }

        virtual void set_signal() = 0;

        virtual void handle_signal() = 0;

        void reset_signal()
        {
            signal->set = 0;
        }

        signed char set() const
        {
            return signal->set;
        }

        signed char sig() const
        {
            return signal->sig;
        }

    protected:
        // bool sync = false;
        gmx::SimulationSignal *signal;
};

class CheckpointSignalHandler : public SignalHandler
{
    public:
        CheckpointSignalHandler(gmx::SimulationSignal   * sig,
                                const t_inputrec        * ir,
                                const t_commrec         * cr,
                                const MdrunOptions       &mdrunOptions,
                                const gmx_bool           &bGStat,
                                const gmx_bool           &bNS,
                                const gmx_bool           &bLastStep,
                                const gmx_int64_t        &step,
                                gmx_bool                 &bCPT,
                                gmx_walltime_accounting_t walltime_accounting);

        void set_signal() override;
        void handle_signal() override;

    private:
        int                       nchkpt = 1;
        const int                 nstlist;
        const gmx_int64_t         init_step;
        const real                cpt_period;
        const bool                isMaster;
        const bool                isParallel;
        const bool                writeConfout;
        const bool                isRerun;
        const gmx_bool           &bGStat;
        const gmx_bool           &bNS;
        const gmx_bool           &bLastStep;
        const gmx_int64_t        &step;
        gmx_bool                 &bCPT;
        gmx_walltime_accounting_t walltime_accounting;
};

class StopSignalHandler : public SignalHandler
{
    public:
        StopSignalHandler(gmx::SimulationSignal    *sig,
                          const t_inputrec        * ir,
                          const t_commrec         * cr,
                          const MdrunOptions       &mdrunOptions,
                          int                       nstSignalComm,
                          FILE                    * fplog,
                          const gmx_bool           &bNS,
                          const gmx_int64_t        &step,
                          gmx_bool                 &bLastStep,
                          gmx_walltime_accounting_t walltime_accounting);

        void set_signal() override;
        void handle_signal() override;

    private:
        int                       handled_stop_condition = gmx_stop_cond_none;
        const real                maximumHoursToRun;
        const bool                isMaster;
        const bool                reproducible;
        const int                 nstSignalComm;
        const int                 nstlist;
        const gmx_bool           &bNS;
        const gmx_int64_t        &step;
        gmx_bool                 &bLastStep;
        FILE                     *fplog;
        gmx_walltime_accounting_t walltime_accounting;

};

class ResetSignalHandler : public SignalHandler
{
    public:
        ResetSignalHandler(gmx::SimulationSignal     * sig,
                           const t_inputrec          * ir,
                           const t_commrec           * cr,
                           const MdrunOptions         &mdrunOptions,
                           const gmx_int64_t          &step,
                           const gmx_int64_t          &step_rel,
                           gmx_walltime_accounting_t   walltime_accounting,
                           gmx_wallcycle_t             wcycle,
                           const gmx::MDLogger        &mdlog,
                           FILE                      * fplog,
                           nonbonded_verlet_t        * nbv,
                           const gmx_pme_t           * pme,
                           t_nrnb                    * nrnb,
                           const pme_load_balancing_t* pme_loadbal);

        ~ResetSignalHandler() override
        {
            if (valid_finish())
            {
                walltime_accounting_set_valid_finish(walltime_accounting);
            }
        }

        void set_signal() override;
        void handle_signal() override;

    private:
        bool valid_finish() const;

        bool                              bResetCountersHalfMaxH;
        const real                        maximumHoursToRun;
        const bool                        isMaster;
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

}

#endif //GROMACS_SIGNALHANDLERS_H
