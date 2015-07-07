/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/legacyheaders/sighandler.h"

#include "config.h"

#include <stdlib.h>

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/utility/fatalerror.h"

const char *gmx_stop_cond_name[] =
{
    "None",
    "Stop at the next neighbor search step",
    "Stop at the next step",
    "Abort"
};

/* these do not neccesarily match the stop condition, but are
   referred to in the signal handler. */
const char *gmx_signal_name[] =
{
    "None",
    "INT",
    "TERM",
    "second INT/TERM",
    "remote INT/TERM",
    "remote second INT/TERM",
    "USR1",
    "Abort"
};

static volatile sig_atomic_t stop_condition   = gmx_stop_cond_none;
static volatile sig_atomic_t last_signal_name = 0;

static volatile sig_atomic_t usr_condition = 0;

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
            stop_condition++;
            if (n == SIGINT)
            {
                last_signal_name = 1;
            }
            if (n == SIGTERM)
            {
                last_signal_name = 2;
            }
            if (stop_condition == gmx_stop_cond_next)
            {
                last_signal_name = 3;
            }
            if (stop_condition >= gmx_stop_cond_abort)
            {
                abort();
            }
            break;
#ifdef HAVE_SIGUSR1
        case SIGUSR1:
            usr_condition = 1;
            break;
#endif
        default:
            break;
    }
}

static void gmx_signal(int signum)
{
#ifdef HAVE_SIGACTION
    struct sigaction act;
    act.sa_handler = signal_handler;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_RESTART;
    sigaction(signum, &act, NULL);
#else
    signal(signum, signal_handler);
#endif
}

void signal_handler_install(void)
{
    if (getenv("GMX_NO_TERM") == NULL)
    {
        if (debug)
        {
            fprintf(debug, "Installing signal handler for SIGTERM\n");
        }
        gmx_signal(SIGTERM);
    }
    if (getenv("GMX_NO_INT") == NULL)
    {
        if (debug)
        {
            fprintf(debug, "Installing signal handler for SIGINT\n");
        }
        gmx_signal(SIGINT);
    }
#ifdef HAVE_SIGUSR1
    if (getenv("GMX_NO_USR1") == NULL)
    {
        if (debug)
        {
            fprintf(debug, "Installing signal handler for SIGUSR1\n");
        }
        gmx_signal(SIGUSR1);
    }
#endif
}

gmx_stop_cond_t gmx_get_stop_condition(void)
{
    return (gmx_stop_cond_t)stop_condition;
}

void gmx_set_stop_condition(gmx_stop_cond_t recvd_stop_cond)
{
    if (recvd_stop_cond > stop_condition)
    {
        stop_condition = recvd_stop_cond;
        if (stop_condition == gmx_stop_cond_next_ns)
        {
            last_signal_name = 4;
        }
        if (stop_condition == gmx_stop_cond_next)
        {
            last_signal_name = 5;
        }
    }
}

const char *gmx_get_signal_name(void)
{
    return gmx_signal_name[last_signal_name];
}

gmx_bool gmx_got_usr_signal(void)
{
#ifdef HAVE_SIGUSR1
    gmx_bool ret = (gmx_bool)usr_condition;
    usr_condition = 0;
    return ret;
#else
    return FALSE;
#endif
}
