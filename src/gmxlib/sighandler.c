/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "gmx_fatal.h"
#include "sighandler.h"


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

static volatile sig_atomic_t stop_condition=gmx_stop_cond_none;
static volatile sig_atomic_t last_signal_name=0;

static volatile sig_atomic_t usr_condition=0;

static RETSIGTYPE signal_handler(int n)
{
    switch (n) {
/* windows doesn't do SIGINT correctly according to ANSI (yes, signals are in 
   ANSI C89, and windows spawns a thread specifically to run the INT signal 
   handler), but that doesn't matter for a simple signal handler like this. */
        case SIGTERM:
        case SIGINT:
            /* we explicitly set things up to allow this: */
            stop_condition++;
            if (n==SIGINT)
                last_signal_name=1;
            if (n==SIGTERM)
                last_signal_name=2;
            if (stop_condition == gmx_stop_cond_next)
                last_signal_name=3;
            if (stop_condition >= gmx_stop_cond_abort)
                abort();
            break;
#ifdef HAVE_SIGUSR1
        case SIGUSR1:
            usr_condition=1;
            break;
#endif
        default:
            break;
    }
}


void signal_handler_install(void)
{
    if (getenv("GMX_NO_TERM") == NULL)
    {
        if (debug)
        {
            fprintf(debug,"Installing signal handler for SIGTERM\n");
        }
        signal(SIGTERM,signal_handler);
    }
    if (getenv("GMX_NO_INT") == NULL)
    {
        if (debug)
        {
            fprintf(debug,"Installing signal handler for SIGINT\n");
        }
        signal(SIGINT,signal_handler);
    }
#ifdef HAVE_SIGUSR1
    if (getenv("GMX_NO_USR1") == NULL)
    {
        if (debug)
        {
            fprintf(debug,"Installing signal handler for SIGUSR1\n");
        }
        signal(SIGUSR1,signal_handler);
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
        stop_condition=recvd_stop_cond;
        if (stop_condition == gmx_stop_cond_next_ns)
            last_signal_name=4;
        if (stop_condition == gmx_stop_cond_next)
            last_signal_name=5;
    }
}

const char *gmx_get_signal_name(void)
{
    return gmx_signal_name[last_signal_name];
}

gmx_bool gmx_got_usr_signal(void)
{
#ifdef HAVE_SIGUSR1
    gmx_bool ret=(gmx_bool)usr_condition;
    usr_condition=0;
    return ret;
#else
    return FALSE;
#endif
}


