/*
 * This file is part of the GROMACS molecular simulation package.
 *
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

#include "trajectory_writing.h"
#include "gromacs/legacyheaders/mvdata.h"
#include "gromacs/legacyheaders/domdec.h"
#include "checkpoint.h"
#include "trnio.h"
#include "xtcio.h"
#include "tngio.h"
#include "gromacs/legacyheaders/smalloc.h"

static void moveit(t_commrec *cr, rvec xx[])
{
    if (!xx)
    {
        return;
    }

    move_rvecs(cr, FALSE, FALSE, xx, NULL, (cr->nnodes-cr->npmenodes)-1, NULL);
}

void write_traj(FILE *fplog, t_commrec *cr,
                gmx_mdoutf_t *of,
                int mdof_flags,
                gmx_mtop_t *top_global,
                gmx_large_int_t step, double t,
                t_state *state_local, t_state *state_global,
                rvec *f_local, rvec *f_global,
                int *n_xtc, rvec **x_xtc)
{
    int           i, j;
    gmx_groups_t *groups;
    rvec         *xxtc;
    rvec         *local_v;
    rvec         *global_v;

#define MX(xvf) moveit(cr, xvf)

    /* MRS -- defining these variables is to manage the difference
     * between half step and full step velocities, but there must be a better way . . . */

    local_v  = state_local->v;
    global_v = state_global->v;

    if (DOMAINDECOMP(cr))
    {
        if (mdof_flags & MDOF_CPT)
        {
            dd_collect_state(cr->dd, state_local, state_global);
        }
        else
        {
            if (mdof_flags & (MDOF_X | MDOF_XTC))
            {
                dd_collect_vec(cr->dd, state_local, state_local->x,
                               state_global->x);
            }
            if (mdof_flags & MDOF_V)
            {
                dd_collect_vec(cr->dd, state_local, local_v,
                               global_v);
            }
        }
        if (mdof_flags & MDOF_F)
        {
            dd_collect_vec(cr->dd, state_local, f_local, f_global);
        }
    }
    else
    {
        if (mdof_flags & MDOF_CPT)
        {
            /* All pointers in state_local are equal to state_global,
             * but we need to copy the non-pointer entries.
             */
            state_global->lambda = state_local->lambda;
            state_global->veta   = state_local->veta;
            state_global->vol0   = state_local->vol0;
            copy_mat(state_local->box, state_global->box);
            copy_mat(state_local->boxv, state_global->boxv);
            copy_mat(state_local->svir_prev, state_global->svir_prev);
            copy_mat(state_local->fvir_prev, state_global->fvir_prev);
            copy_mat(state_local->pres_prev, state_global->pres_prev);
        }
        if (cr->nnodes > 1)
        {
            /* Particle decomposition, collect the data on the master node */
            if (mdof_flags & MDOF_CPT)
            {
                if (state_local->flags & (1<<estX))
                {
                    MX(state_global->x);
                }
                if (state_local->flags & (1<<estV))
                {
                    MX(state_global->v);
                }
                if (state_local->flags & (1<<estSDX))
                {
                    MX(state_global->sd_X);
                }
                if (state_global->nrngi > 1)
                {
                    if (state_local->flags & (1<<estLD_RNG))
                    {
#ifdef GMX_MPI
                        MPI_Gather(state_local->ld_rng,
                                   state_local->nrng*sizeof(state_local->ld_rng[0]), MPI_BYTE,
                                   state_global->ld_rng,
                                   state_local->nrng*sizeof(state_local->ld_rng[0]), MPI_BYTE,
                                   MASTERRANK(cr), cr->mpi_comm_mygroup);
#endif
                    }
                    if (state_local->flags & (1<<estLD_RNGI))
                    {
#ifdef GMX_MPI
                        MPI_Gather(state_local->ld_rngi,
                                   sizeof(state_local->ld_rngi[0]), MPI_BYTE,
                                   state_global->ld_rngi,
                                   sizeof(state_local->ld_rngi[0]), MPI_BYTE,
                                   MASTERRANK(cr), cr->mpi_comm_mygroup);
#endif
                    }
                }
            }
            else
            {
                if (mdof_flags & (MDOF_X | MDOF_XTC))
                {
                    MX(state_global->x);
                }
                if (mdof_flags & MDOF_V)
                {
                    MX(global_v);
                }
            }
            if (mdof_flags & MDOF_F)
            {
                MX(f_global);
            }
        }
    }

    if (MASTER(cr))
    {
        if (mdof_flags & MDOF_CPT)
        {
            write_checkpoint(of->fn_cpt, of->bKeepAndNumCPT,
                             fplog, cr, of->eIntegrator, of->simulation_part,
                             of->bExpanded, of->elamstats, step, t, state_global);
        }

        if (mdof_flags & (MDOF_X | MDOF_V | MDOF_F))
        {
            if (of->fp_trn)
            {
                fwrite_trn(of->fp_trn, step, t, state_local->lambda[efptFEP],
                           state_local->box, top_global->natoms,
                           (mdof_flags & MDOF_X) ? state_global->x : NULL,
                           (mdof_flags & MDOF_V) ? global_v : NULL,
                           (mdof_flags & MDOF_F) ? f_global : NULL);
                if (gmx_fio_flush(of->fp_trn) != 0)
                {
                    gmx_file("Cannot write trajectory; maybe you are out of disk space?");
                }
                gmx_fio_check_file_position(of->fp_trn);
            }

            fwrite_tng(of->tng, FALSE, step, t, state_local->lambda[efptFEP],
                       (const rvec *) state_local->box,
                       top_global->natoms,
                       (mdof_flags & MDOF_X) ? (const rvec *) state_global->x : NULL,
                       (mdof_flags & MDOF_V) ? (const rvec *) global_v : NULL,
                       (mdof_flags & MDOF_F) ? (const rvec *) f_global : NULL);
            /* TODO: implement gmx_fio_check_file_position kind of
               functionality? What does checkpointing really need? */
        }
        if (mdof_flags & MDOF_XTC)
        {
            groups = &top_global->groups;
            if (*n_xtc == -1)
            {
                *n_xtc = 0;
                for (i = 0; (i < top_global->natoms); i++)
                {
                    if (ggrpnr(groups, egcXTC, i) == 0)
                    {
                        (*n_xtc)++;
                    }
                }
                if (*n_xtc != top_global->natoms)
                {
                    snew(*x_xtc, *n_xtc);
                }
            }
            if (*n_xtc == top_global->natoms)
            {
                xxtc = state_global->x;
            }
            else
            {
                xxtc = *x_xtc;
                j    = 0;
                for (i = 0; (i < top_global->natoms); i++)
                {
                    if (ggrpnr(groups, egcXTC, i) == 0)
                    {
                        copy_rvec(state_global->x[i], xxtc[j++]);
                    }
                }
            }
            if (write_xtc(of->fp_xtc, *n_xtc, step, t,
                          state_local->box, xxtc, of->xtc_prec) == 0)
            {
                gmx_fatal(FARGS, "XTC error - maybe you are out of disk space?");
            }
            fwrite_tng(of->tng_low_prec,
                       TRUE,
                       step,
                       t,
                       state_local->lambda[efptFEP],
                       (const rvec *) state_local->box,
                       *n_xtc,
                       (const rvec *) xxtc,
                       NULL,
                       NULL);
            gmx_fio_check_file_position(of->fp_xtc);
        }
    }
}
