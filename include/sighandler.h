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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _sighandler_h
#define _sighandler_h

#include <signal.h>

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* NOTE: the terminology is:
   incoming signals (provided by the operating system, or transmitted from 
   other nodes) lead to stop conditions. These stop conditions should be 
   checked for and acted on by the outer loop of the simulation */

/* the stop conditions. They are explicitly allowed to be compared against
   each other. */
typedef enum
{
    gmx_stop_cond_none=0,
    gmx_stop_cond_next_ns, /* stop a the next neighbour searching step */
    gmx_stop_cond_next, /* stop a the next step */
    gmx_stop_cond_abort  /* stop now. (this should never be seen) */
} gmx_stop_cond_t;

/* Our names for the stop conditions. 
   These must match the number given in gmx_stop_cond_t.*/
extern const char *gmx_stop_cond_name[];

/* the externally visible functions: */

/* install the signal handlers that can set the stop condition. */
void signal_handler_install(void);

/* get the current stop condition */
gmx_stop_cond_t gmx_get_stop_condition(void);

/* set the stop condition upon receiving a remote one */
void gmx_set_stop_condition(gmx_stop_cond_t recvd_stop_cond);

/* get the signal name that lead to the current stop condition. */
const char *gmx_get_signal_name(void);

/* check whether we received a USR1 signal. 
   The condition is reset once a TRUE value is returned, so this function
   only returns TRUE once for a single signal. */
gmx_bool gmx_got_usr_signal(void);


#ifdef __cplusplus
}
#endif


#endif	/* _sighandler_h */
