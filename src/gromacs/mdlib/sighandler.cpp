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
#include "gmxpre.h"

#include "sighandler.h"

#include "config.h"

#include <csignal>
#include <cstdio>
#include <cstdlib>

#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

const char* enumValueToString(StopCondition enumValue)
{
    constexpr gmx::EnumerationArray<StopCondition, const char*> stopConditionNames = {
        "None", "Stop at the next neighbor search step", "Stop at the next step", "Abort"
    };
    return stopConditionNames[enumValue];
}

/* these do not neccesarily match the stop condition, but are
   referred to in the signal handler. */
static const char* const gmx_signal_name[] = {
    "None", "INT",  "TERM", "second INT/TERM", "remote INT/TERM", "remote second INT/TERM",
    "USR1", "Abort"
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static volatile StopCondition stop_condition = StopCondition::None;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static volatile sig_atomic_t last_signal_name = 0;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static volatile sig_atomic_t usr_condition = 0;

void gmx_reset_stop_condition()
{
    stop_condition = StopCondition::None;
    // last_signal_name and usr_condition are left untouched by reset.
}

static void signal_handler(int n)
{
    switch (n)
    {
            /* windows doesn't do SIGINT correctly according to ANSI (yes, signals are in
               ANSI C89, and windows spawns a thread specifically to run the INT signal
               handler), but that doesn't matter for a simple signal handler like this. */
        case SIGTERM:
        case SIGINT:
            /* we explicitly set things up to allow this: */
            switch (stop_condition)
            {
                case StopCondition::None: stop_condition = StopCondition::NextNS; break;
                case StopCondition::NextNS: stop_condition = StopCondition::Next; break;
                case StopCondition::Next: stop_condition = StopCondition::Abort; break;
                default: GMX_THROW(gmx::InternalError("Stop condition increased beyond abort"));
            }
            if (n == SIGINT)
            {
                last_signal_name = 1;
            }
            if (n == SIGTERM)
            {
                last_signal_name = 2;
            }
            if (stop_condition == StopCondition::Next)
            {
                last_signal_name = 3;
            }
            if (stop_condition >= StopCondition::Abort)
            {
                abort();
            }
            break;
#if HAVE_SIGUSR1
        case SIGUSR1: usr_condition = 1; break;
#endif
        default: break;
    }
}

static void gmx_signal(int signum)
{
#if HAVE_SIGACTION
    struct sigaction act;
    act.sa_handler = signal_handler;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_RESTART;
    sigaction(signum, &act, nullptr);
#else
    signal(signum, signal_handler);
#endif
}

void signal_handler_install()
{
    if (getenv("GMX_NO_TERM") == nullptr)
    {
        if (debug)
        {
            fprintf(debug, "Installing signal handler for SIGTERM\n");
        }
        gmx_signal(SIGTERM);
    }
    if (getenv("GMX_NO_INT") == nullptr)
    {
        if (debug)
        {
            fprintf(debug, "Installing signal handler for SIGINT\n");
        }
        gmx_signal(SIGINT);
    }
#if HAVE_SIGUSR1
    if (getenv("GMX_NO_USR1") == nullptr)
    {
        if (debug)
        {
            fprintf(debug, "Installing signal handler for SIGUSR1\n");
        }
        gmx_signal(SIGUSR1);
    }
#endif
}

StopCondition gmx_get_stop_condition()
{
    return stop_condition;
}

void gmx_set_stop_condition(StopCondition recvd_stop_cond)
{
    if (recvd_stop_cond > stop_condition)
    {
        stop_condition = recvd_stop_cond;
        if (stop_condition == StopCondition::NextNS)
        {
            last_signal_name = 4;
        }
        if (stop_condition == StopCondition::Next)
        {
            last_signal_name = 5;
        }
    }
}

const char* gmx_get_signal_name()
{
    return gmx_signal_name[last_signal_name];
}

gmx_bool gmx_got_usr_signal()
{
#if HAVE_SIGUSR1
    gmx_bool ret  = static_cast<gmx_bool>(usr_condition);
    usr_condition = 0;
    return ret;
#else
    return FALSE;
#endif
}
