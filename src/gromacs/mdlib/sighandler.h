/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#ifndef GMX_MDLIB_SIGHANDLER_H
#define GMX_MDLIB_SIGHANDLER_H

#include <csignal>

#include "gromacs/utility/basedefinitions.h"

/* NOTE: the terminology is:
   incoming signals (provided by the operating system, or transmitted from
   other nodes) lead to stop conditions. These stop conditions should be
   checked for and acted on by the outer loop of the simulation */

/* the stop conditions. They are explicitly allowed to be compared against
   each other. */
enum class StopCondition : sig_atomic_t
{
    None = 0,
    NextNS, /* stop a the next neighbour searching step */
    Next,   /* stop a the next step */
    Abort,  /* stop now. (this should never be seen) */
    Count
};

/* Our names for the stop conditions.
   These must match the number given in gmx_stop_cond_t.*/
const char* enumValueToString(StopCondition enumValue);

/* the externally visible functions: */

/* install the signal handlers that can set the stop condition. */
void signal_handler_install();

/* get the current stop condition */
StopCondition gmx_get_stop_condition();

/* set the stop condition upon receiving a remote one */
void gmx_set_stop_condition(StopCondition recvd_stop_cond);

/*!
 * \brief Reinitializes the global stop condition.
 *
 * Resets any stop condition currently stored in global library state as read or
 * written with gmx_get_stop_condition() and gmx_set_stop_condition(). Does not
 * affect the result of gmx_got_usr_signal() gmx_get_signal_name() for
 * previously terminated simulations.
 *
 * The reset is necessary between simulation segments performed in the same
 * process and should be called only while simulation is idle, such as after
 * a gmx::Mdrunner has finished its work and simulation results have been processed.
 */
void gmx_reset_stop_condition();

/* get the signal name that lead to the current stop condition. */
const char* gmx_get_signal_name();

/* check whether we received a USR1 signal.
   The condition is reset once a TRUE value is returned, so this function
   only returns TRUE once for a single signal. */
gmx_bool gmx_got_usr_signal();

#endif
