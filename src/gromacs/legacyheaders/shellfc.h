/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/timing/wallcycle.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_graph;

/* Initialization function, also predicts the initial shell postions.
 * If x!=NULL, the shells are predict for the global coordinates x.
 */
gmx_shellfc_t init_shell_flexcon(FILE *fplog,
                                 gmx_mtop_t *mtop, int nflexcon,
                                 rvec *x);

/* Get the local shell with domain decomposition */
void make_local_shells(t_commrec *cr, t_mdatoms *md,
                       gmx_shellfc_t shfc);

/* Optimize shell positions */
int relax_shell_flexcon(FILE *log, t_commrec *cr, gmx_bool bVerbose,
                        gmx_int64_t mdstep, t_inputrec *inputrec,
                        gmx_bool bDoNS, int force_flags,
                        gmx_localtop_t *top,
                        gmx_constr_t constr,
                        gmx_enerdata_t *enerd, t_fcdata *fcd,
                        t_state *state, rvec f[],
                        tensor force_vir,
                        t_mdatoms *md,
                        t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                        struct t_graph *graph,
                        gmx_groups_t *groups,
                        gmx_shellfc_t shfc,
                        t_forcerec *fr,
                        gmx_bool bBornRadii,
                        double t, rvec mu_tot,
                        gmx_bool *bConverged,
                        gmx_vsite_t *vsite,
                        FILE *fp_field);


#ifdef __cplusplus
}
#endif
