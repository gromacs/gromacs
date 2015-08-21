/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_MINIMIZE_H
#define GMX_MINIMIZE_H

#include "gmxpre.h"
#include "config.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/types/state.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"

namespace gmx
{

typedef struct {
    t_state  s;
    rvec    *f;
    real     epot;
    real     fnorm;
    real     fmax;
    int      a_fmax;
} em_state_t;

em_state_t *init_em_state();

void init_em(FILE *fplog, const char *title,
             t_commrec *cr, t_inputrec *ir,
             t_state *state_global, gmx_mtop_t *top_global,
             em_state_t *ems, gmx_localtop_t **top,
             rvec **f, rvec **f_global,
             t_nrnb *nrnb, rvec mu_tot,
             t_forcerec *fr, gmx_enerdata_t **enerd,
             t_graph **graph, t_mdatoms *mdatoms, gmx_global_stat_t *gstat,
             gmx_vsite_t *vsite, gmx_constr_t constr,
             int nfile, const t_filenm fnm[],
             gmx_mdoutf_t *outf, t_mdebin **mdebin,
             int imdport, unsigned long gmx_unused Flags,
             gmx_wallcycle_t wcycle);

void finish_em(t_commrec *cr, gmx_mdoutf_t outf,
               gmx_walltime_accounting_t walltime_accounting,
               gmx_wallcycle_t wcycle);

void write_em_traj(FILE *fplog, t_commrec *cr,
                   gmx_mdoutf_t outf,
                   gmx_bool bX, gmx_bool bF, const char *confout,
                   gmx_mtop_t *top_global,
                   t_inputrec *ir, gmx_int64_t step,
                   em_state_t *state,
                   t_state *state_global, rvec *f_global);

void do_em_step(t_commrec *cr, t_inputrec *ir, t_mdatoms *md,
                gmx_bool bMolPBC,
                em_state_t *ems1, real a, rvec *f, em_state_t *ems2,
                gmx_constr_t constr, gmx_localtop_t *top,
                t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                gmx_int64_t count);

void evaluate_em_energy(FILE *fplog, t_commrec *cr,
                        gmx_mtop_t *top_global,
                        em_state_t *ems, gmx_localtop_t *top,
                        t_inputrec *inputrec,
                        t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                        gmx_global_stat_t gstat,
                        gmx_vsite_t *vsite, gmx_constr_t constr,
                        t_fcdata *fcd,
                        t_graph *graph, t_mdatoms *mdatoms,
                        t_forcerec *fr, rvec mu_tot,
                        gmx_enerdata_t *enerd, tensor vir, tensor pres,
                        gmx_int64_t count, gmx_bool bFirst);

} /* namespace gmx */

#endif
