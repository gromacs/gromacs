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
#include "mdoutf.h"

#include "gromacs/legacyheaders/xvgr.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/mvdata.h"
#include "gromacs/legacyheaders/domdec.h"
#include "trnio.h"
#include "xtcio.h"
#include "tngio.h"
#include "trajectory_writing.h"
#include "checkpoint.h"

struct gmx_mdoutf {
    t_fileio         *fp_trn;
    t_fileio         *fp_xtc;
    tng_trajectory_t  tng;
    tng_trajectory_t  tng_low_prec;
    int               xtc_prec;
    ener_file_t       fp_ene;
    const char       *fn_cpt;
    gmx_bool          bKeepAndNumCPT;
    int               eIntegrator;
    gmx_bool          bExpanded;
    int               elamstats;
    int               simulation_part;
    FILE             *fp_dhdl;
    FILE             *fp_field;
};


gmx_mdoutf_t init_mdoutf(int nfile, const t_filenm fnm[], int mdrun_flags,
                         const t_commrec *cr, const t_inputrec *ir,
                         gmx_mtop_t *mtop,
                         const output_env_t oenv)
{
    gmx_mdoutf_t  of;
    char          filemode[3];
    gmx_bool      bAppendFiles;

    snew(of, 1);

    of->fp_trn       = NULL;
    of->fp_ene       = NULL;
    of->fp_xtc       = NULL;
    of->tng          = NULL;
    of->tng_low_prec = NULL;
    of->fp_dhdl      = NULL;
    of->fp_field     = NULL;

    of->eIntegrator     = ir->eI;
    of->bExpanded       = ir->bExpanded;
    of->elamstats       = ir->expandedvals->elamstats;
    of->simulation_part = ir->simulation_part;

    if (MASTER(cr))
    {
        bAppendFiles = (mdrun_flags & MD_APPENDFILES);

        of->bKeepAndNumCPT = (mdrun_flags & MD_KEEPANDNUMCPT);

        sprintf(filemode, bAppendFiles ? "a+" : "w+");

        if ((EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI))
#ifndef GMX_FAHCORE
            &&
            !(EI_DYNAMICS(ir->eI) &&
              ir->nstxout == 0 &&
              ir->nstvout == 0 &&
              ir->nstfout == 0)
#endif
            )
        {
            const char *filename;
            filename = ftp2fn(efTRN, nfile, fnm);
            switch (fn2ftp(filename))
            {
                case efTRR:
                case efTRN:
                    of->fp_trn = open_trn(filename, filemode);
                    break;
                case efTNG:
                    /* Initialize topology and then set writing frequency. */
                    tng_open(filename, filemode[0], &of->tng);
                    if (filemode[0] == 'w')
                    {
                        tng_prepare_md_writing(of->tng, mtop, ir);
                    }
                    break;
                default:
                    gmx_incons("Invalid full precision file format");
            }
        }
        if (EI_DYNAMICS(ir->eI) &&
            ir->nstxtcout > 0)
        {
            const char *filename;
            filename = ftp2fn(efCOMPRESSED, nfile, fnm);
            switch (fn2ftp(filename))
            {
                case efXTC:
                    of->fp_xtc   = open_xtc(filename, filemode);
                    of->xtc_prec = ir->xtcprec;
                    break;
                case efTNG:
                    /* Initialize topology and then set writing frequency. */
                    tng_open(filename, filemode[0], &of->tng_low_prec);
                    if (filemode[0] == 'w')
                    {
                        tng_prepare_low_prec_writing(of->tng_low_prec, mtop, ir);
                    }
                    break;
                default:
                    gmx_incons("Invalid reduced precision file format");
            }
        }
        if (EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI))
        {
            of->fp_ene = open_enx(ftp2fn(efEDR, nfile, fnm), filemode);
        }
        of->fn_cpt = opt2fn("-cpo", nfile, fnm);

        if ((ir->efep != efepNO || ir->bSimTemp) && ir->fepvals->nstdhdl > 0 &&
            (ir->fepvals->separate_dhdl_file == esepdhdlfileYES ) &&
            EI_DYNAMICS(ir->eI))
        {
            if (bAppendFiles)
            {
                of->fp_dhdl = gmx_fio_fopen(opt2fn("-dhdl", nfile, fnm), filemode);
            }
            else
            {
                of->fp_dhdl = open_dhdl(opt2fn("-dhdl", nfile, fnm), ir, oenv);
            }
        }

        if (opt2bSet("-field", nfile, fnm) &&
            (ir->ex[XX].n || ir->ex[YY].n || ir->ex[ZZ].n))
        {
            if (bAppendFiles)
            {
                of->fp_dhdl = gmx_fio_fopen(opt2fn("-field", nfile, fnm),
                                            filemode);
            }
            else
            {
                of->fp_field = xvgropen(opt2fn("-field", nfile, fnm),
                                        "Applied electric field", "Time (ps)",
                                        "E (V/nm)", oenv);
            }
        }
    }

    return of;
}

FILE *mdoutf_get_fp_field(gmx_mdoutf_t of)
{
    return of->fp_field;
}

ener_file_t mdoutf_get_fp_ene(gmx_mdoutf_t of)
{
    return of->fp_ene;
}

FILE *mdoutf_get_fp_dhdl(gmx_mdoutf_t of)
{
    return of->fp_dhdl;
}

static void moveit(t_commrec *cr, rvec xx[])
{
    if (!xx)
    {
        return;
    }

    move_rvecs(cr, FALSE, FALSE, xx, NULL, (cr->nnodes-cr->npmenodes)-1, NULL);
}

void mdoutf_write_to_trajectory_files(FILE *fplog, t_commrec *cr,
                                      gmx_mdoutf_t of,
                                      int mdof_flags,
                                      gmx_mtop_t *top_global,
                                      gmx_int64_t step, double t,
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
            fflush_tng(of->tng);
            fflush_tng(of->tng_low_prec);
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
            if (of->fp_xtc)
            {
                gmx_fio_check_file_position(of->fp_xtc);
            }
        }
    }
}

void done_mdoutf(gmx_mdoutf_t of)
{
    if (of->fp_ene != NULL)
    {
        close_enx(of->fp_ene);
    }
    if (of->fp_xtc)
    {
        close_xtc(of->fp_xtc);
    }
    if (of->fp_trn)
    {
        close_trn(of->fp_trn);
    }
    if (of->fp_dhdl != NULL)
    {
        gmx_fio_fclose(of->fp_dhdl);
    }
    if (of->fp_field != NULL)
    {
        gmx_fio_fclose(of->fp_field);
    }
    /* We have to check that of->tng and of->tng_low_prec are set because
     * tng_close returns a NULL pointer. */
    if (of->tng)
    {
        tng_close(&of->tng);
    }
    if (of->tng_low_prec)
    {
        tng_close(&of->tng_low_prec);
    }

    sfree(of);
}
