/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.6.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2011, The GROMACS development team,
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

#ifndef _pme_loadbal_h
#define _pme_loadbal_h

typedef struct pme_loadbal *pme_loadbal_t;

/* Initialze the PME grid tuning data and infrastructure */
void pme_loadbal_init(pme_loadbal_t *pmes_p,
                     const t_inputrec *ir,matrix box,
                     const interaction_const_t *ic,
                     gmx_pme_t pmedata);

/* Adjust the PME grid and Coulomb cut-off.
 * Returns TRUE the tuning continues, FALSE is the tuning is done.
 */
gmx_bool pme_loadbalance(pme_loadbal_t pmes,
                         t_commrec *cr,
                         FILE *fp_err,
                         FILE *fp_log,
                         t_inputrec *ir,
                         t_state *state,
                         double cycles,
                         interaction_const_t *ic,
                         nonbonded_verlet_t *nbv,
                         gmx_pme_t *pmedata,
                         int step);

/* Restart the PME tuning, discarding all timings gathered up till now */
void restart_pme_loadbal(pme_loadbal_t pmes, int n);

/* Finish the PME tuning and print the settings when fplog!=NULL */
void pme_loadbal_done(pme_loadbal_t pmes, FILE *fplog);

#endif /* _pme_loadbal_h */
