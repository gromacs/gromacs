/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
#include "gmxpre.h"

#include "trajectory_writing.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"

void
do_md_trajectory_writing(FILE               *fplog,
                         t_commrec          *cr,
                         int                 nfile,
                         const t_filenm      fnm[],
                         gmx_int64_t         step,
                         gmx_int64_t         step_rel,
                         double              t,
                         t_inputrec         *ir,
                         t_state            *state,
                         t_state            *state_global,
                         ObservablesHistory *observablesHistory,
                         gmx_mtop_t         *top_global,
                         t_forcerec         *fr,
                         gmx_mdoutf_t        outf,
                         t_mdebin           *mdebin,
                         gmx_ekindata_t     *ekind,
                         PaddedRVecVector   *f,
                         int                *nchkpt,
                         gmx_bool            bCPT,
                         gmx_bool            bRerunMD,
                         gmx_bool            bLastStep,
                         gmx_bool            bDoConfOut,
                         gmx_bool            bSumEkinhOld
                         )
{
    int   mdof_flags;
    rvec *x_for_confout = nullptr;

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

#if defined(GMX_FAHCORE)
    if (bLastStep)
    {
        /* Enforce writing positions and velocities at end of run */
        mdof_flags |= (MDOF_X | MDOF_V);
    }
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

                update_energyhistory(observablesHistory->energyHistory.get(), mdebin);
            }
        }
        mdoutf_write_to_trajectory_files(fplog, cr, outf, mdof_flags, top_global,
                                         step, t, state, state_global, observablesHistory, f);
        if (bCPT)
        {
            (*nchkpt)++;
        }
        if (bLastStep && step_rel == ir->nsteps &&
            bDoConfOut && MASTER(cr) &&
            !bRerunMD)
        {
            if (fr->bMolPBC && state == state_global)
            {
                /* This (single-rank) run needs to allocate a
                   temporary array of size natoms so that any
                   periodicity removal for mdrun -confout does not
                   perturb the update and thus the final .edr
                   output. This makes .cpt restarts look binary
                   identical, and makes .edr restarts binary
                   identical. */
                snew(x_for_confout, state_global->natoms);
                copy_rvecn(as_rvec_array(state_global->x.data()), x_for_confout, 0, state_global->natoms);
            }
            else
            {
                /* With DD, or no bMolPBC, it doesn't matter if
                   we change as_rvec_array(state_global->x.data()) */
                x_for_confout = as_rvec_array(state_global->x.data());
            }

            /* x and v have been collected in mdoutf_write_to_trajectory_files,
             * because a checkpoint file will always be written
             * at the last step.
             */
            fprintf(stderr, "\nWriting final coordinates.\n");
            if (fr->bMolPBC && !ir->bPeriodicMols)
            {
                /* Make molecules whole only for confout writing */
                do_pbc_mtop(fplog, ir->ePBC, state->box, top_global, x_for_confout);
            }
            write_sto_conf_mtop(ftp2fn(efSTO, nfile, fnm),
                                *top_global->name, top_global,
                                x_for_confout, as_rvec_array(state_global->v.data()),
                                ir->ePBC, state->box);
            if (fr->bMolPBC && state == state_global)
            {
                sfree(x_for_confout);
            }
        }
        wallcycle_stop(mdoutf_get_wcycle(outf), ewcTRAJ);
    }
}
