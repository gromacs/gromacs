/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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

#ifndef _pme_loadbal_h
#define _pme_loadbal_h

typedef struct pme_load_balancing *pme_load_balancing_t;

/* Initialze the PP-PME load balacing data and infrastructure */
void pme_loadbal_init(pme_load_balancing_t *pme_lb_p,
                      const t_inputrec *ir, matrix box,
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
                          t_commrec           *cr,
                          FILE                *fp_err,
                          FILE                *fp_log,
                          t_inputrec          *ir,
                          t_state             *state,
                          double               cycles,
                          interaction_const_t *ic,
                          nonbonded_verlet_t  *nbv,
                          gmx_pme_t           *pmedata,
                          gmx_int64_t          step);

/* Restart the PME load balancing discarding all timings gathered up till now */
void restart_pme_loadbal(pme_load_balancing_t pme_lb, int n);

/* Finish the PME load balancing and print the settings when fplog!=NULL */
void pme_loadbal_done(pme_load_balancing_t pme_lb,
                      t_commrec *cr, FILE *fplog,
                      gmx_bool bNonBondedOnGPU);

#endif /* _pme_loadbal_h */
