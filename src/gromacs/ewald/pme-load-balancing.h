/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_EWALD_PME_LOAD_BALANCING_H
#define GMX_EWALD_PME_LOAD_BALANCING_H

#include "gromacs/legacyheaders/types/commrec_fwd.h"
#include "gromacs/legacyheaders/types/forcerec.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/legacyheaders/types/state.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pme_load_balancing *pme_load_balancing_t;

/* Initialze the PP-PME load balacing data and infrastructure.
 * Returns if the load balancing is activated in bPMELoadBalActive.
 * Active means either already balancing or that we should measure timings
 * and potentially activate it later.
 */
void pme_loadbal_init(pme_load_balancing_t *pme_lb_p,
                      t_commrec *cr, FILE *fp_log,
                      const t_inputrec *ir, matrix box,
                      const interaction_const_t *ic,
                      struct gmx_pme_t *pmedata,
                      gmx_bool bUseGPU,
                      gmx_bool *bPMELoadBalActive);

/* Do PME load balancing, or only check if trigger PME load balancing */
void pme_loadbal_do(pme_load_balancing_t  pme_lb,
                    t_commrec            *cr,
                    FILE                 *fp_err,
                    FILE                 *fp_log,
                    t_inputrec           *ir,
                    t_forcerec           *fr,
                    t_state              *state,
                    double                cycles,
                    gmx_int64_t           step,
                    gmx_int64_t           step_rel,
                    gmx_bool             *bActive);

/* Finish the PME load balancing and print the settings when fplog!=NULL */
void pme_loadbal_done(pme_load_balancing_t pme_lb,
                      t_commrec *cr, FILE *fplog,
                      gmx_bool bNonBondedOnGPU);

#ifdef __cplusplus
}
#endif

#endif
