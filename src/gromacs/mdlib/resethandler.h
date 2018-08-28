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
 * Declares the reset handler class.
 *
 * This class resets the various counters based on either the time (master rank sends
 * checkpointing signal after 49.5% or run time), or based on the number of elapsed
 * steps (handled locally by slaves).
 *
 * The approach is ready for a task-based design: The class binds to data it needs to
 * access at run time during construction time via const references. This allows the
 * task to run later without any inputs.
 *
 * The setting and handling is implemented in private functions. They are only called
 * if a respective boolean is true. For the trivial case of no reset needed (or no reset
 * signal setting on any other rank than master), the translation unit of the calling
 * function is therefore never left. The current implementation also allows the handler
 * and setters to be ignored once a reset has been done, as a reset is only allowed to
 * happen once. In the future, many of these cases this will be achieved by adding
 * (or not adding) handlers / setters to the task graph.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_RESETHANDLER_H
#define GMX_MDLIB_RESETHANDLER_H

#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/timing/walltime_accounting.h"

namespace gmx
{
/*! \libinternal
 * \brief Class handling the reset of counters
 *
 * Master rank sets the reset signal if half the run time is reached.
 * All ranks receive the reset signal and reset their respective counters.
 * This also resets the counters if half the time steps have passed (no communication needed).
 */
class ResetHandler
{
    public:
        /*! \brief ResetHandler constructor
         *
         * Needs a pointer to the signal to communicate between ranks, information on whether
         * multiple simulations need to be synchronized, and (const) references to data it
         * needs to determine the when the counters should be reset. The latter allows the
         * setSignal() and handleSignal() routines to run without additional arguments,
         * making it easier to be ran in a task-based environment.
         */
        ResetHandler(
            gmx::SimulationSignal      *sig,
            bool                        needSync,
            const t_inputrec           *ir,
            const t_commrec            *cr,
            const MdrunOptions         &mdrunOptions,
            const gmx::MDLogger        &mdlog,
            FILE                       *fplog,
            const int64_t              &step,
            const int64_t              &step_rel,
            nonbonded_verlet_t         *nbv,
            t_nrnb                     *nrnb,
            const gmx_pme_t            *pme,
            const pme_load_balancing_t *pme_loadbal,
            gmx_walltime_accounting_t   walltime_accounting,
            gmx_wallcycle_t             wcycle);

        /*! \brief Decides whether a reset signal needs to be set
         *
         * Reset signal is set if run time is greater than 49.5% of maximal run time.
         */
        void setSignal()
        {
            if (doSet)
            {
                if (setSignalImpl())
                {
                    // need to set the reset signal only once
                    doSet = false;
                }
            }
        }

        /*! \brief Decides whether the counters are reset, and performs the reset if needed
         *
         * The counters are reset if
         *   * the signal for resetting was received, or
         *   * the (local) number of steps reached the defined counter reset step.
         * Note that even if two reset conditions are present (at a specific step and a
         * spedific time), the reset will only take place once, whenever the first condition
         * is met.
         */
        void handleSignal()
        {
            if (doHandle)
            {
                if (handleSignalImpl())
                {
                    // need to reset the counters only once
                    doHandle = false;
                    doSet    = false;
                }
            }
        }

    private:
        bool setSignalImpl();
        bool handleSignalImpl();

        bool doSet    = false;
        bool doHandle = false;

        const real                        maximumHoursToRun;

        SimulationSignal                 *signal;

        const int64_t                    &step;
        const int64_t                    &step_rel;

        const gmx::MDLogger              &mdlog;
        FILE* const                       fplog;
        const t_commrec* const            cr;
        nonbonded_verlet_t* const         nbv;
        t_nrnb* const                     nrnb;
        const gmx_pme_t* const            pme;
        const pme_load_balancing_t* const pme_loadbal;

        gmx_wallcycle_t                   wcycle;
        gmx_walltime_accounting_t         walltime_accounting;
};
}      // namespace gmx

#endif // GMX_MDLIB_RESETHANDLER_H
