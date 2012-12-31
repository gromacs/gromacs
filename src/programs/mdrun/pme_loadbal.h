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

typedef struct pme_load_balancing *pme_load_balancing_t;

/* Initialze the PP-PME load balacing data and infrastructure */
void pme_loadbal_init(pme_load_balancing_t *pme_lb_p,
                      const t_inputrec *ir,matrix box,
                      const interaction_const_t *ic,
                      gmx_pme_t pmedata);

/* Try to adjust the PME grid and Coulomb cut-off.
 * The adjustment is done to generate a different non-bonded PP and PME load.
 * With separate PME nodes (PP and PME on different processes) or with
 * a GPU (PP on GPU, PME on CPU), PP and PME run on different resources
 * and changing the load will affect the load balance and performance.
 * The total time for a set of integration steps is monitored and a range
 * of grid/cut-off setups is scanned. After calling pme_load_balance many
 * times and acquiring enough statistics, the best performing setup is chosen.
 * Here we try to take into account fluctuations and changes due to external
 * factors as well as DD load balancing.
 * Returns TRUE the load balancing continues, FALSE is the balancing is done.
 */
gmx_bool pme_load_balance(pme_load_balancing_t pme_lb,
                          t_commrec *cr,
                          FILE *fp_err,
                          FILE *fp_log,
                          t_inputrec *ir,
                          t_state *state,
                          double cycles,
                          interaction_const_t *ic,
                          nonbonded_verlet_t *nbv,
                          gmx_pme_t *pmedata,
                          gmx_large_int_t step);

/* Restart the PME load balancing discarding all timings gathered up till now */
void restart_pme_loadbal(pme_load_balancing_t pme_lb, int n);

/* Finish the PME load balancing and print the settings when fplog!=NULL */
void pme_loadbal_done(pme_load_balancing_t pme_lb, FILE *fplog);

#endif /* _pme_loadbal_h */
