/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
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

#ifndef _gmx_membed_h
#define _gmx_membed_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* initialisation of membed code */
gmx_membed_t init_membed(FILE *fplog, int nfile, const t_filenm fnm[], gmx_mtop_t *mtop,
                         t_inputrec *inputrec, t_state *state, t_commrec *cr, real *cpt);

/* rescaling the coordinates voor de membed code */
void rescale_membed(int step_rel, gmx_membed_t membed, rvec *x);

#ifdef __cplusplus
}
#endif

#endif
