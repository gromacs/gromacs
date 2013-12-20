/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * Copyright (c) 2013, by the GROMACS development team, led by
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

#ifndef GMX_FILEIO_TRAJECTORY_WRITING_H
#define GMX_FILEIO_TRAJECTORY_WRITING_H

#include <stdio.h>
#include "filenm.h"
#include "mdoutf.h"
#include "../legacyheaders/types/simple.h"
#include "../legacyheaders/types/commrec.h"
#include "../legacyheaders/update.h"
#include "../legacyheaders/mdebin.h"

void
do_trajectory_writing(FILE           *fplog,
                      t_commrec      *cr,
                      int             nfile,
                      const t_filenm  fnm[],
                      gmx_int64_t     step,
                      gmx_int64_t     step_rel,
                      double          t,
                      t_inputrec     *ir,
                      t_state        *state,
                      t_state        *state_global,
                      gmx_mtop_t     *top_global,
                      t_forcerec     *fr,
                      gmx_update_t    upd,
                      gmx_mdoutf_t   *outf,
                      t_mdebin       *mdebin,
                      gmx_ekindata_t *ekind,
                      rvec           *f,
                      rvec           *f_global,
                      gmx_wallcycle_t wcycle,
                      gmx_rng_t       mcrng,
                      int            *nchkpt,
                      gmx_bool        bCPT,
                      gmx_bool        bRerunMD,
                      gmx_bool        bLastStep,
                      gmx_bool        bDoConfOut,
                      gmx_bool        bSumEkinhOld
                      );
/* Wrapper routine for trajectory writing */

void write_traj(FILE *fplog, t_commrec *cr,
                gmx_mdoutf_t *of,
                int mdof_flags,
                gmx_mtop_t *top_global,
                gmx_int64_t step, double t,
                t_state *state_local, t_state *state_global,
                rvec *f_local, rvec *f_global,
                int *n_xtc, rvec **x_xtc);
/* Routine that writes frames to trn, xtc and/or checkpoint.
 * What is written is determined by the mdof_flags defined above.
 * Data is collected to the master node only when necessary.
 */

#endif /* GMX_FILEIO_TRAJECTORY_WRITING_H */
