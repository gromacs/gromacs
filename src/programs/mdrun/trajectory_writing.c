/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "sim_util.h"
#include "gmx_wallcycle.h"
#include "mdrun.h"
#include "confio.h"

void
do_trajectory_writing(FILE           *fplog,
                      t_commrec      *cr,
                      int             nfile,
                      const t_filenm  fnm[],
                      gmx_large_int_t step,
                      gmx_large_int_t step_rel,
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
                      df_history_t    df_history,
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
                      )
{
    int   mdof_flags;
    int   n_xtc    = -1;
    rvec *x_xtc    = NULL;

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
    if (do_per_step(step, ir->nstxtcout))
    {
        mdof_flags |= MDOF_XTC;
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

    /* sync bCPT and fc record-keeping */
    if (bCPT && MASTER(cr))
    {
        fcRequestCheckPoint();
    }
#endif

    if (mdof_flags != 0)
    {
        wallcycle_start(wcycle, ewcTRAJ);
        if (bCPT)
        {
            if (state->flags & (1<<estLD_RNG))
            {
                get_stochd_state(upd, state);
            }
            if (state->flags  & (1<<estMC_RNG))
            {
                get_mc_state(mcrng, state);
            }
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
                if (ir->efep != efepNO || ir->bSimTemp)
                {
                    state_global->fep_state = state->fep_state;     /* MRS: seems kludgy. The code should be
                                                                       structured so this isn't necessary.
                                                                       Note this reassignment is only necessary
                                                                       for single threads.*/
                    copy_df_history(&state_global->dfhist, &df_history);
                }
            }
        }
        write_traj(fplog, cr, outf, mdof_flags, top_global,
                   step, t, state, state_global, f, f_global, &n_xtc, &x_xtc);
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
            /* x and v have been collected in write_traj,
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
        wallcycle_stop(wcycle, ewcTRAJ);
    }
}
