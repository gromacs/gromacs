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

/* The following two variables and the signal_handler function
 * are used from md.c and pme.c as well 
 *
 * Do not fear these global variables: they represent inherently process-global
 * information that needs to be shared across threads 
 */


/* we got a signal to stop in the next step: */
volatile sig_atomic_t bGotStopNextStepSignal=FALSE;
/* we got a signal to stop in the next neighbour search step: */
volatile sig_atomic_t bGotStopNextNSStepSignal=FALSE;

/* our names for the handled signals. These must match the number given
   in signal_handler. */
const char *signal_name[] =
{
    "TERM",
    "INT",
    "second INT"
};

/* the last signal received, according to the numbering
   we use in signal_name */
volatile sig_atomic_t last_signal_number_recvd=-1;

RETSIGTYPE signal_handler(int n)
{
    switch (n) {
        case SIGTERM:
            bGotStopNextStepSignal = TRUE;
            last_signal_number_recvd = 0;
            break;
/* windows doesn't do SIGINT correctly according to ANSI (yes, signals are in 
   ANSI C89, and windows spawns a thread specifically to run the INT signal 
   handler), but that doesn't matter for a simple signal handler like this. */
        case SIGINT:
            if (!bGotStopNextNSStepSignal)
            {
                bGotStopNextNSStepSignal = TRUE;
                last_signal_number_recvd = 1;
            }            else if (!bGotStopNextStepSignal)
            {
                bGotStopNextStepSignal = TRUE;
                last_signal_number_recvd = 2;
            }
            else
                abort();
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

}



