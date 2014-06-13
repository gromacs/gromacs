/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares code for calling GROMACS routines from an external program
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_qmmm
 */
#ifndef GMX_MINIMIZE_H
#define GMX_MINIMIZE_H

#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/legacyheaders/network.h"

#ifdef __cplusplus
extern "C" {
#endif

//! Datastructure used for minimizing the energy
typedef struct {
    t_state  s;
    rvec    *f;
    rvec    *A;
    real    *phi;
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
                 gmx_mdoutf_t **outf, t_mdebin **mdebin);

    void evaluate_energy(FILE *fplog, t_commrec *cr,
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
                         gmx_large_int_t count, gmx_bool bFirst);



#ifdef __cplusplus
}
#endif

#endif
