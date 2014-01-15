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
#include "mdoutf.h"

#include "gromacs/legacyheaders/xvgr.h"
#include "trnio.h"
#include "xtcio.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/smalloc.h"

gmx_mdoutf_t *init_mdoutf(int nfile, const t_filenm fnm[], int mdrun_flags,
                          const t_commrec *cr, const t_inputrec *ir,
                          gmx_mtop_t *top_global,
                          const output_env_t oenv)
{
    gmx_mdoutf_t *of;
    char          filemode[3];
    gmx_bool      bAppendFiles;
    int           i;

    snew(of, 1);

    of->fp_trn   = NULL;
    of->fp_ene   = NULL;
    of->fp_xtc   = NULL;
    of->fp_dhdl  = NULL;
    of->fp_field = NULL;

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
            of->fp_trn = open_trn(ftp2fn(efTRN, nfile, fnm), filemode);
        }
        if (EI_DYNAMICS(ir->eI) &&
            ir->nstxtcout > 0)
        {
            of->fp_xtc   = open_xtc(ftp2fn(efXTC, nfile, fnm), filemode);
            of->xtc_prec = ir->xtcprec;
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

        /* Set up atom counts so they can be passed to actual
           trajectory-writing routines later. Also, XTC writing needs
           to know what (and how many) atoms might be in the XTC
           groups, and how to look up later which ones they are. */
        of->natoms_global = top_global->natoms;
        of->groups        = &top_global->groups;
        of->natoms_xtc    = 0;
        for (i = 0; (i < top_global->natoms); i++)
        {
            if (ggrpnr(of->groups, egcXTC, i) == 0)
            {
                of->natoms_xtc++;
            }
        }
    }

    return of;
}

void done_mdoutf(gmx_mdoutf_t *of)
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

    sfree(of);
}
