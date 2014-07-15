/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "sim_util.h"
#include "mdrun.h"
#include "confio.h"
#include "trajectory_writing.h"
#include "mdoutf.h"

#include "gromacs/timing/wallcycle.h"

void
do_md_trajectory_writing(FILE           *fplog,
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
                         gmx_mdoutf_t    outf,
                         t_mdebin       *mdebin,
                         gmx_ekindata_t *ekind,
                         rvec           *f,
                         rvec           *f_global,
                         int            *nchkpt,
                         gmx_bool        bCPT,
                         gmx_bool        bRerunMD,
                         gmx_bool        bLastStep,
                         gmx_bool        bDoConfOut,
                         gmx_bool        bSumEkinhOld
                         )
{
    int   mdof_flags;

    mdof_flags = 0;
    if (do_per_step(step, ir->nstxout))
    {
        mdof_flags |= MDOF_X;
    }
    if (do_per_step(step, ir->nstvout))
    {
        mdof_flags |= MDOF_V;
    }
    if (do_per_step(step, ir->nstfout))
    {
        mdof_flags |= MDOF_F;
    }
    if (do_per_step(step, ir->nstxout_compressed))
    {
        mdof_flags |= MDOF_X_COMPRESSED;
    }
    if (bCPT)
    {
        mdof_flags |= MDOF_CPT;
    }
    ;

#if defined(GMX_FAHCORE) || defined(GMX_WRITELASTSTEP)
    if (bLastStep)
    {
        /* Enforce writing positions and velocities at end of run */
        mdof_flags |= (MDOF_X | MDOF_V);
    }
#endif
#ifdef GMX_FAHCORE
    if (MASTER(cr))
    {
        fcReportProgress( ir->nsteps, step );
    }

#if defined(__native_client__)
    fcCheckin(MASTER(cr));
#endif

    /* sync bCPT and fc record-keeping */
    if (bCPT && MASTER(cr))
    {
        fcRequestCheckPoint();
    }
#endif

    if (mdof_flags != 0)
    {
        wallcycle_start(mdoutf_get_wcycle(outf), ewcTRAJ);
        if (bCPT)
        {
            if (MASTER(cr))
            {
                if (bSumEkinhOld)
                {
                    state_global->ekinstate.bUpToDate = FALSE;
                }
                else
                {
                    update_ekinstate(&state_global->ekinstate, ekind);
                    state_global->ekinstate.bUpToDate = TRUE;
                }
                update_energyhistory(&state_global->enerhist, mdebin);
            }
        }
        mdoutf_write_to_trajectory_files(fplog, cr, outf, mdof_flags, top_global,
                                         step, t, state, state_global, f, f_global);
        if (bCPT)
        {
            (*nchkpt)++;
            bCPT = FALSE;
        }
        debug_gmx();
        if (bLastStep && step_rel == ir->nsteps &&
            bDoConfOut && MASTER(cr) &&
            !bRerunMD)
        {
            /* x and v have been collected in mdoutf_write_to_trajectory_files,
             * because a checkpoint file will always be written
             * at the last step.
             */
            fprintf(stderr, "\nWriting final coordinates.\n");
            if (fr->bMolPBC)
            {
                /* Make molecules whole only for confout writing */
                do_pbc_mtop(fplog, ir->ePBC, state->box, top_global, state_global->x);
            }
            write_sto_conf_mtop(ftp2fn(efSTO, nfile, fnm),
                                *top_global->name, top_global,
                                state_global->x, state_global->v,
                                ir->ePBC, state->box);
            debug_gmx();
        }
        wallcycle_stop(mdoutf_get_wcycle(outf), ewcTRAJ);
    }
}
