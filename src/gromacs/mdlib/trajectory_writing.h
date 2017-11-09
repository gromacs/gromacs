/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_TRAJECTORY_WRITING_H
#define GMX_MDLIB_TRAJECTORY_WRITING_H

#include <stdio.h>

#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/timing/wallcycle.h"

struct gmx_ekindata_t;
struct gmx_mtop_t;
struct ObservablesHistory;
struct t_commrec;
struct t_filenm;

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Wrapper routine for writing trajectories during mdrun
 *
 * This routine does communication (e.g. collecting distributed coordinates)
 */
void
do_md_trajectory_writing(FILE                     *fplog,
                         struct t_commrec         *cr,
                         int                       nfile,
                         const t_filenm            fnm[],
                         gmx_int64_t               step,
                         gmx_int64_t               step_rel,
                         double                    t,
                         t_inputrec               *ir,
                         t_state                  *state,
                         t_state                  *state_global,
                         ObservablesHistory       *observablesHistory,
                         struct gmx_mtop_t        *top_global,
                         t_forcerec               *fr,
                         gmx_mdoutf_t              outf,
                         t_mdebin                 *mdebin,
                         struct gmx_ekindata_t    *ekind,
                         gmx::ArrayRef<gmx::RVec>  f,
                         int                      *nchkpt,
                         gmx_bool                  bCPT,
                         gmx_bool                  bRerunMD,
                         gmx_bool                  bLastStep,
                         gmx_bool                  bDoConfOut,
                         gmx_bool                  bSumEkinhOld
                         );


#ifdef __cplusplus
}
#endif

#endif
