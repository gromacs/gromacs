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

#ifndef _runtime_h
#define _runtime_h

#include "gromacs/legacyheaders/types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gmx_runtime *gmx_runtime_t;

void runtime_start(gmx_runtime_t runtime);

void runtime_end(gmx_runtime_t runtime);

double runtime_get_elapsed_run_time(gmx_runtime_t runtime);

double runtime_get_start_time_stamp(gmx_runtime_t runtime);

double runtime_get_nsteps_done(gmx_runtime_t runtime);

void runtime_set_nsteps_done(gmx_runtime_t   runtime,
                             gmx_large_int_t _nsteps_done);
// TODO consider whether this should get done in runtime_end

double gmx_gettime();

#ifdef __cplusplus
}
#endif

#endif  /* _runtime_h */
