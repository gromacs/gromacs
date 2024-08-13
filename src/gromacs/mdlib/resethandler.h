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
 * \brief
 * Declares the reset handler class.
 *
 * This class resets the various counters based on either the time (main rank sends
 * checkpointing signal after 49.5% or run time), or based on the number of elapsed
 * steps (handled locally by all ranks independently). Resets can happen in different
 * ways:
 *
 *     * at a predetermined step (gmx mdrun -resetstep XXX)
 *     * at half of the number of steps (gmx mdrun -resethway and nsteps set)
 *     * at half of the max wall time (gmx mdrun -resethway -maxh XX), which is
 *       implemented triggered when walltime >= 49.5% of max
 *
 * If two or more of these reset conditions are set, the first condition which is met
 * resets the counters, there is no second reset happening. Note also that
 * -resethway with nsteps set overwrites -resetstep
 * (gmx mdrun -resethway -nsteps 100000 -resetstep 1000 will result in a reset at step
 * 50000, not 1000).
 *
 * The setting and handling is implemented in private functions. They are only called
 * if a respective boolean is true. For the trivial case of no reset needed (or no reset
 * signal setting on any other rank than main), the translation unit of the calling
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

#include <cstdint>
#include <cstdio>

#include "gromacs/compat/pointers.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

struct gmx_pme_t;
struct gmx_wallcycle;
struct gmx_walltime_accounting;
struct t_nrnb;
struct t_commrec;
struct pme_load_balancing_t;

namespace gmx
{
struct nonbonded_verlet_t;

/*! \brief Reset signals
 *
 * Signals set and read by ResetHandler. Possible signals include
 *   * nothing to signal
 *   * reset counters (as soon as signal is received)
 */
enum class ResetSignal
{
    noSignal        = 0,
    doResetCounters = 1
};

/*! \libinternal
 * \brief Class handling the reset of counters
 *
 * Main rank sets the reset signal if half the run time is reached.
 * All ranks receive the reset signal and reset their respective counters.
 * This also resets the counters if half the time steps have passed (no communication needed).
 */
class ResetHandler final
{
public:
    /*! \brief ResetHandler constructor
     *
     * Needs a pointer to the signal to communicate between ranks, information on whether
     * multiple simulations need to be synchronized, and additional data to determine
     * whether counter resetting takes place at all, and whether the current rank can set
     * the resetting signal.
     */
    ResetHandler(compat::not_null<SimulationSignal*> signal,
                 bool                                simulationsShareState,
                 int64_t                             nsteps,
                 bool                                isMain,
                 bool                                resetHalfway,
                 real                                maximumHoursToRun,
                 const MDLogger&                     mdlog,
                 gmx_wallcycle*                      wcycle,
                 gmx_walltime_accounting*            walltime_accounting);

    /*! \brief Decides whether a reset signal needs to be set
     *
     * Reset signal is set if run time is greater than 49.5% of maximal run time.
     */
    void setSignal(gmx_walltime_accounting* walltime_accounting)
    {
        if (rankCanSetSignal_)
        {
            if (setSignalImpl(walltime_accounting))
            {
                // need to set the reset signal only once
                rankCanSetSignal_ = false;
            }
        }
    }

    /*! \brief Decides whether the counters are reset, and performs the reset if needed
     *
     * The counters are reset if
     *
     *     * the signal for resetting was received, or
     *     * the (local) number of steps reached the defined counter reset step.
     *
     * Note that even if two reset conditions are present (at a specific step and a
     * specific time), the reset will only take place once, whenever the first condition
     * is met.
     */
    void resetCounters(int64_t                     step,
                       int64_t                     step_rel,
                       const MDLogger&             mdlog,
                       FILE*                       fplog,
                       const t_commrec*            cr,
                       nonbonded_verlet_t*         nbv,
                       t_nrnb*                     nrnb,
                       const gmx_pme_t*            pme,
                       const pme_load_balancing_t* pme_loadbal,
                       gmx_wallcycle*              wcycle,
                       gmx_walltime_accounting*    walltime_accounting)
    {
        if (simulationNeedsReset_)
        {
            if (resetCountersImpl(step, step_rel, mdlog, fplog, cr, nbv, nrnb, pme, pme_loadbal, wcycle, walltime_accounting))
            {
                // need to reset the counters only once
                simulationNeedsReset_ = false;
                rankCanSetSignal_     = false;
            }
        }
    }

private:
    //! Implementation of the setSignal() function
    bool setSignalImpl(gmx_walltime_accounting* walltime_accounting);

    //! Implementation of the resetCounters() function
    bool resetCountersImpl(int64_t                     step,
                           int64_t                     step_rel,
                           const MDLogger&             mdlog,
                           FILE*                       fplog,
                           const t_commrec*            cr,
                           nonbonded_verlet_t*         nbv,
                           t_nrnb*                     nrnb,
                           const gmx_pme_t*            pme,
                           const pme_load_balancing_t* pme_loadbal,
                           gmx_wallcycle*              wcycle,
                           gmx_walltime_accounting*    walltime_accounting);

    SimulationSignal& signal_;

    bool       rankCanSetSignal_;
    bool       simulationNeedsReset_;
    const real maximumHoursToRun_;
};
} // namespace gmx

#endif // GMX_MDLIB_RESETHANDLER_H
