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

#include "mdoutf.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

struct gmx_mdoutf {
    t_fileio               *fp_trn;
    t_fileio               *fp_xtc;
    tng_trajectory_t        tng;
    tng_trajectory_t        tng_low_prec;
    int                     x_compression_precision; /* only used by XTC output */
    ener_file_t             fp_ene;
    const char             *fn_cpt;
    gmx_bool                bKeepAndNumCPT;
    int                     eIntegrator;
    gmx_bool                bExpanded;
    int                     elamstats;
    int                     simulation_part;
    FILE                   *fp_dhdl;
    int                     natoms_global;
    int                     natoms_x_compressed;
    gmx_groups_t           *groups; /* for compressed position writing */
    gmx_wallcycle_t         wcycle;
    rvec                   *f_global;
    gmx::IMDOutputProvider *outputProvider;
};


gmx_mdoutf_t init_mdoutf(FILE *fplog, int nfile, const t_filenm fnm[],
                         const MdrunOptions &mdrunOptions,
                         const t_commrec *cr,
                         gmx::IMDOutputProvider *outputProvider,
                         const t_inputrec *ir, gmx_mtop_t *top_global,
                         const gmx_output_env_t *oenv, gmx_wallcycle_t wcycle)
{
    gmx_mdoutf_t   of;
    const char    *appendMode = "a+", *writeMode = "w+", *filemode;
    gmx_bool       bAppendFiles, bCiteTng = FALSE;
    int            i;

    snew(of, 1);

    of->fp_trn       = nullptr;
    of->fp_ene       = nullptr;
    of->fp_xtc       = nullptr;
    of->tng          = nullptr;
    of->tng_low_prec = nullptr;
    of->fp_dhdl      = nullptr;

    of->eIntegrator             = ir->eI;
    of->bExpanded               = ir->bExpanded;
    of->elamstats               = ir->expandedvals->elamstats;
    of->simulation_part         = ir->simulation_part;
    of->x_compression_precision = static_cast<int>(ir->x_compression_precision);
    of->wcycle                  = wcycle;
    of->f_global                = nullptr;
    of->outputProvider          = outputProvider;

    if (MASTER(cr))
    {
        bAppendFiles = mdrunOptions.continuationOptions.appendFiles;

        of->bKeepAndNumCPT = mdrunOptions.checkpointOptions.keepAndNumberCheckpointFiles;

        filemode = bAppendFiles ? appendMode : writeMode;

        if (EI_DYNAMICS(ir->eI) &&
            ir->nstxout_compressed > 0)
        {
            const char *filename;
            filename = ftp2fn(efCOMPRESSED, nfile, fnm);
            switch (fn2ftp(filename))
            {
                case efXTC:
                    of->fp_xtc                  = open_xtc(filename, filemode);
                    break;
                case efTNG:
                    gmx_tng_open(filename, filemode[0], &of->tng_low_prec);
                    if (filemode[0] == 'w')
                    {
                        gmx_tng_prepare_low_prec_writing(of->tng_low_prec, top_global, ir);
                    }
                    bCiteTng = TRUE;
                    break;
                default:
                    gmx_incons("Invalid reduced precision file format");
            }
        }
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
                    /* If there is no uncompressed coordinate output and
                       there is compressed TNG output write forces
                       and/or velocities to the TNG file instead. */
                    if (ir->nstxout != 0 || ir->nstxout_compressed == 0 ||
                        !of->tng_low_prec)
                    {
                        of->fp_trn = gmx_trr_open(filename, filemode);
                    }
                    break;
                case efTNG:
                    gmx_tng_open(filename, filemode[0], &of->tng);
                    if (filemode[0] == 'w')
                    {
                        gmx_tng_prepare_md_writing(of->tng, top_global, ir);
                    }
                    bCiteTng = TRUE;
                    break;
                default:
                    gmx_incons("Invalid full precision file format");
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

        outputProvider->initOutput(fplog, nfile, fnm, bAppendFiles, oenv);

        /* Set up atom counts so they can be passed to actual
           trajectory-writing routines later. Also, XTC writing needs
           to know what (and how many) atoms might be in the XTC
           groups, and how to look up later which ones they are. */
        of->natoms_global       = top_global->natoms;
        of->groups              = &top_global->groups;
        of->natoms_x_compressed = 0;
        for (i = 0; (i < top_global->natoms); i++)
        {
            if (ggrpnr(of->groups, egcCompressedX, i) == 0)
            {
                of->natoms_x_compressed++;
            }
        }

        if (ir->nstfout && DOMAINDECOMP(cr))
        {
            snew(of->f_global, top_global->natoms);
        }
    }

    if (bCiteTng)
    {
        please_cite(fplog, "Lundborg2014");
    }

    return of;
}

ener_file_t mdoutf_get_fp_ene(gmx_mdoutf_t of)
{
    return of->fp_ene;
}

FILE *mdoutf_get_fp_dhdl(gmx_mdoutf_t of)
{
    return of->fp_dhdl;
}

gmx_wallcycle_t mdoutf_get_wcycle(gmx_mdoutf_t of)
{
    return of->wcycle;
}

void mdoutf_write_to_trajectory_files(FILE *fplog, t_commrec *cr,
                                      gmx_mdoutf_t of,
                                      int mdof_flags,
                                      gmx_mtop_t *top_global,
                                      gmx_int64_t step, double t,
                                      t_state *state_local, t_state *state_global,
                                      ObservablesHistory *observablesHistory,
                                      gmx::ArrayRef<gmx::RVec> f_local)
{
    rvec *f_global;

    if (DOMAINDECOMP(cr))
    {
        if (mdof_flags & MDOF_CPT)
        {
            dd_collect_state(cr->dd, state_local, state_global);
        }
        else
        {
            if (mdof_flags & (MDOF_X | MDOF_X_COMPRESSED))
            {
                gmx::ArrayRef<gmx::RVec> globalXRef = MASTER(cr) ? gmx::makeArrayRef(state_global->x) : gmx::EmptyArrayRef();
                dd_collect_vec(cr->dd, state_local, state_local->x, globalXRef);
            }
            if (mdof_flags & MDOF_V)
            {
                gmx::ArrayRef<gmx::RVec> globalVRef = MASTER(cr) ? gmx::makeArrayRef(state_global->v) : gmx::EmptyArrayRef();
                dd_collect_vec(cr->dd, state_local, state_local->v, globalVRef);
            }
        }
        f_global = of->f_global;
        if (mdof_flags & MDOF_F)
        {
            dd_collect_vec(cr->dd, state_local, f_local, gmx::arrayRefFromArray(reinterpret_cast<gmx::RVec *>(f_global), f_local.size()));
        }
    }
    else
    {
        /* We have the whole state locally: copy the local state pointer */
        state_global = state_local;

        f_global     = as_rvec_array(f_local.data());
    }

    if (MASTER(cr))
    {
        if (mdof_flags & MDOF_CPT)
        {
            fflush_tng(of->tng);
            fflush_tng(of->tng_low_prec);
            ivec one_ivec = { 1, 1, 1 };
            write_checkpoint(of->fn_cpt, of->bKeepAndNumCPT,
                             fplog, cr,
                             DOMAINDECOMP(cr) ? cr->dd->nc : one_ivec,
                             DOMAINDECOMP(cr) ? cr->dd->nnodes : cr->nnodes,
                             of->eIntegrator, of->simulation_part,
                             of->bExpanded, of->elamstats, step, t,
                             state_global, observablesHistory);
        }

        if (mdof_flags & (MDOF_X | MDOF_V | MDOF_F))
        {
            const rvec *x = (mdof_flags & MDOF_X) ? as_rvec_array(state_global->x.data()) : nullptr;
            const rvec *v = (mdof_flags & MDOF_V) ? as_rvec_array(state_global->v.data()) : nullptr;
            const rvec *f = (mdof_flags & MDOF_F) ? f_global : nullptr;

            if (of->fp_trn)
            {
                gmx_trr_write_frame(of->fp_trn, step, t, state_local->lambda[efptFEP],
                                    state_local->box, top_global->natoms,
                                    x, v, f);
                if (gmx_fio_flush(of->fp_trn) != 0)
                {
                    gmx_file("Cannot write trajectory; maybe you are out of disk space?");
                }
            }

            /* If a TNG file is open for uncompressed coordinate output also write
               velocities and forces to it. */
            else if (of->tng)
            {
                gmx_fwrite_tng(of->tng, FALSE, step, t, state_local->lambda[efptFEP],
                               state_local->box,
                               top_global->natoms,
                               x, v, f);
            }
            /* If only a TNG file is open for compressed coordinate output (no uncompressed
               coordinate output) also write forces and velocities to it. */
            else if (of->tng_low_prec)
            {
                gmx_fwrite_tng(of->tng_low_prec, FALSE, step, t, state_local->lambda[efptFEP],
                               state_local->box,
                               top_global->natoms,
                               x, v, f);
            }
        }
        if (mdof_flags & MDOF_X_COMPRESSED)
        {
            rvec *xxtc = nullptr;

            if (of->natoms_x_compressed == of->natoms_global)
            {
                /* We are writing the positions of all of the atoms to
                   the compressed output */
                xxtc = as_rvec_array(state_global->x.data());
            }
            else
            {
                /* We are writing the positions of only a subset of
                   the atoms to the compressed output, so we have to
                   make a copy of the subset of coordinates. */
                int i, j;

                snew(xxtc, of->natoms_x_compressed);
                for (i = 0, j = 0; (i < of->natoms_global); i++)
                {
                    if (ggrpnr(of->groups, egcCompressedX, i) == 0)
                    {
                        copy_rvec(state_global->x[i], xxtc[j++]);
                    }
                }
            }
            if (write_xtc(of->fp_xtc, of->natoms_x_compressed, step, t,
                          state_local->box, xxtc, of->x_compression_precision) == 0)
            {
                gmx_fatal(FARGS, "XTC error - maybe you are out of disk space?");
            }
            gmx_fwrite_tng(of->tng_low_prec,
                           TRUE,
                           step,
                           t,
                           state_local->lambda[efptFEP],
                           state_local->box,
                           of->natoms_x_compressed,
                           xxtc,
                           nullptr,
                           nullptr);
            if (of->natoms_x_compressed != of->natoms_global)
            {
                sfree(xxtc);
            }
        }
    }
}

void mdoutf_tng_close(gmx_mdoutf_t of)
{
    if (of->tng || of->tng_low_prec)
    {
        wallcycle_start(of->wcycle, ewcTRAJ);
        gmx_tng_close(&of->tng);
        gmx_tng_close(&of->tng_low_prec);
        wallcycle_stop(of->wcycle, ewcTRAJ);
    }
}

void done_mdoutf(gmx_mdoutf_t of)
{
    if (of->fp_ene != nullptr)
    {
        close_enx(of->fp_ene);
    }
    if (of->fp_xtc)
    {
        close_xtc(of->fp_xtc);
    }
    if (of->fp_trn)
    {
        gmx_trr_close(of->fp_trn);
    }
    if (of->fp_dhdl != nullptr)
    {
        gmx_fio_fclose(of->fp_dhdl);
    }
    of->outputProvider->finishOutput();
    if (of->f_global != nullptr)
    {
        sfree(of->f_global);
    }

    gmx_tng_close(&of->tng);
    gmx_tng_close(&of->tng_low_prec);

    sfree(of);
}
