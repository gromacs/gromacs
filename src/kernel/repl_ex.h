/*
 * $Id$
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _repl_ex_h
#define _repl_ex_h

#include "typedefs.h"

typedef struct {
  int  repl;
  int  nrepl;
  bool bNPT;
  real *temp;
  real *pres;
  int  *ind;
  int  nst;
  int  seed;
  int  nattempt[2];
  real *prob_sum;
  int  *nexchange;
} gmx_repl_ex_t;

extern gmx_repl_ex_t *init_replica_exchange(FILE *fplog,
					    const t_commrec *mcr,
					    const t_state *state,
					    const t_inputrec *ir,
					    int nst,int init_seed);

extern bool replica_exchange(FILE *fplog,
			     const t_commrec *mcr,
			     gmx_repl_ex_t *re,
			     t_state *state,real epot,int step,real time);
/* Attempts replica exchange.
 * Returns TRUE if this state has been exchanged.
 */

extern void print_replica_exchange_statistics(FILE *fplog,gmx_repl_ex_t *re);

#endif	/* _repl_ex_h */
