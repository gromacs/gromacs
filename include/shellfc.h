/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include "visibility.h"
#include "typedefs.h"
#include "vsite.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Initialization function, also predicts the initial shell postions.
 * If x!=NULL, the shells are predict for the global coordinates x.
 */
GMX_LIBMD_EXPORT
gmx_shellfc_t init_shell_flexcon(FILE *fplog,
                                 gmx_mtop_t *mtop, int nflexcon,
                                 rvec *x);

/* Get the local shell with domain decomposition */
GMX_LIBMD_EXPORT
void make_local_shells(t_commrec *cr, t_mdatoms *md,
                       gmx_shellfc_t shfc);

/* Optimize shell positions */
GMX_LIBMD_EXPORT
int relax_shell_flexcon(FILE *log, t_commrec *cr, gmx_bool bVerbose,
                        gmx_large_int_t mdstep, t_inputrec *inputrec,
                        gmx_bool bDoNS, int force_flags,
                        gmx_bool bStopCM,
                        gmx_localtop_t *top,
                        gmx_mtop_t *mtop,
                        gmx_constr_t constr,
                        gmx_enerdata_t *enerd, t_fcdata *fcd,
                        t_state *state, rvec f[],
                        tensor force_vir,
                        t_mdatoms *md,
                        t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                        t_graph *graph,
                        gmx_groups_t *groups,
                        gmx_shellfc_t shfc,
                        t_forcerec *fr,
                        gmx_bool bBornRadii,
                        double t, rvec mu_tot,
                        int natoms, gmx_bool *bConverged,
                        gmx_vsite_t *vsite,
                        FILE *fp_field);


#ifdef __cplusplus
}
#endif
