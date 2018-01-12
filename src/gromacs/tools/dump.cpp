/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2013, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "dump.h"

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/mtxio.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/gmxpreprocess/gmxcpp.h"
#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

static void list_tpx(const char *fn,
                     gmx_bool    bShowNumbers,
                     gmx_bool    bShowParameters,
                     const char *mdpfn,
                     gmx_bool    bSysTop,
                     gmx_bool    bOriginalInputrec)
{
    FILE         *gp;
    int           indent, i, j, **gcount, atot;
    t_state       state;
    t_tpxheader   tpx;
    gmx_mtop_t    mtop;
    gmx_groups_t *groups;
    t_topology    top;

    read_tpxheader(fn, &tpx, TRUE);
    t_inputrec     ir;
    read_tpx_state(fn,
                   tpx.bIr ? &ir : nullptr,
                   &state,
                   tpx.bTop ? &mtop : nullptr);
    if (tpx.bIr && !bOriginalInputrec)
    {
        gmx::MDModules().adjustInputrecBasedOnModules(&ir);
    }

    if (mdpfn && tpx.bIr)
    {
        gp = gmx_fio_fopen(mdpfn, "w");
        pr_inputrec(gp, 0, nullptr, &ir, TRUE);
        gmx_fio_fclose(gp);
    }

    if (!mdpfn)
    {
        if (bSysTop)
        {
            top = gmx_mtop_t_to_t_topology(&mtop, false);
        }

        if (available(stdout, &tpx, 0, fn))
        {
            indent = 0;
            pr_title(stdout, indent, fn);
            pr_inputrec(stdout, 0, "inputrec", tpx.bIr ? &ir : nullptr, FALSE);

            pr_tpxheader(stdout, indent, "header", &(tpx));

            if (!bSysTop)
            {
                pr_mtop(stdout, indent, "topology", &(mtop), bShowNumbers, bShowParameters);
            }
            else
            {
                pr_top(stdout, indent, "topology", &(top), bShowNumbers, bShowParameters);
            }

            pr_rvecs(stdout, indent, "box", tpx.bBox ? state.box : nullptr, DIM);
            pr_rvecs(stdout, indent, "box_rel", tpx.bBox ? state.box_rel : nullptr, DIM);
            pr_rvecs(stdout, indent, "boxv", tpx.bBox ? state.boxv : nullptr, DIM);
            pr_rvecs(stdout, indent, "pres_prev", tpx.bBox ? state.pres_prev : nullptr, DIM);
            pr_rvecs(stdout, indent, "svir_prev", tpx.bBox ? state.svir_prev : nullptr, DIM);
            pr_rvecs(stdout, indent, "fvir_prev", tpx.bBox ? state.fvir_prev : nullptr, DIM);
            /* leave nosehoover_xi in for now to match the tpr version */
            pr_doubles(stdout, indent, "nosehoover_xi", state.nosehoover_xi.data(), state.ngtc);
            /*pr_doubles(stdout,indent,"nosehoover_vxi",state.nosehoover_vxi,state.ngtc);*/
            /*pr_doubles(stdout,indent,"therm_integral",state.therm_integral,state.ngtc);*/
            pr_rvecs(stdout, indent, "x", tpx.bX ? as_rvec_array(state.x.data()) : nullptr, state.natoms);
            pr_rvecs(stdout, indent, "v", tpx.bV ? as_rvec_array(state.v.data()) : nullptr, state.natoms);
        }

        groups = &mtop.groups;

        snew(gcount, egcNR);
        for (i = 0; (i < egcNR); i++)
        {
            snew(gcount[i], groups->grps[i].nr);
        }

        for (i = 0; (i < mtop.natoms); i++)
        {
            for (j = 0; (j < egcNR); j++)
            {
                gcount[j][ggrpnr(groups, j, i)]++;
            }
        }
        printf("Group statistics\n");
        for (i = 0; (i < egcNR); i++)
        {
            atot = 0;
            printf("%-12s: ", gtypes[i]);
            for (j = 0; (j < groups->grps[i].nr); j++)
            {
                printf("  %5d", gcount[i][j]);
                atot += gcount[i][j];
            }
            printf("  (total %d atoms)\n", atot);
            sfree(gcount[i]);
        }
        sfree(gcount);
    }
}

static void list_top(const char *fn)
{
    int       status, done;
#define BUFLEN 256
    char      buf[BUFLEN];
    gmx_cpp_t handle;
    char     *cppopts[] = { nullptr };

    status = cpp_open_file(fn, &handle, cppopts);
    if (status != 0)
    {
        gmx_fatal(FARGS, cpp_error(&handle, status));
    }
    do
    {
        status = cpp_read_line(&handle, BUFLEN, buf);
        done   = (status == eCPP_EOF);
        if (!done)
        {
            if (status != eCPP_OK)
            {
                gmx_fatal(FARGS, cpp_error(&handle, status));
            }
            else
            {
                printf("%s\n", buf);
            }
        }
    }
    while (!done);
    status = cpp_close_file(&handle);
    if (status != eCPP_OK)
    {
        gmx_fatal(FARGS, cpp_error(&handle, status));
    }
}

static void list_trr(const char *fn)
{
    t_fileio         *fpread;
    int               nframe, indent;
    char              buf[256];
    rvec             *x, *v, *f;
    matrix            box;
    gmx_trr_header_t  trrheader;
    gmx_bool          bOK;

    fpread  = gmx_trr_open(fn, "r");

    nframe = 0;
    while (gmx_trr_read_frame_header(fpread, &trrheader, &bOK))
    {
        snew(x, trrheader.natoms);
        snew(v, trrheader.natoms);
        snew(f, trrheader.natoms);
        if (gmx_trr_read_frame_data(fpread, &trrheader,
                                    trrheader.box_size ? box : nullptr,
                                    trrheader.x_size   ? x : nullptr,
                                    trrheader.v_size   ? v : nullptr,
                                    trrheader.f_size   ? f : nullptr))
        {
            sprintf(buf, "%s frame %d", fn, nframe);
            indent = 0;
            indent = pr_title(stdout, indent, buf);
            pr_indent(stdout, indent);
            fprintf(stdout, "natoms=%10d  step=%10" GMX_PRId64 "  time=%12.7e  lambda=%10g\n",
                    trrheader.natoms, trrheader.step, trrheader.t, trrheader.lambda);
            if (trrheader.box_size)
            {
                pr_rvecs(stdout, indent, "box", box, DIM);
            }
            if (trrheader.x_size)
            {
                pr_rvecs(stdout, indent, "x", x, trrheader.natoms);
            }
            if (trrheader.v_size)
            {
                pr_rvecs(stdout, indent, "v", v, trrheader.natoms);
            }
            if (trrheader.f_size)
            {
                pr_rvecs(stdout, indent, "f", f, trrheader.natoms);
            }
        }
        else
        {
            fprintf(stderr, "\nWARNING: Incomplete frame: nr %d, t=%g\n",
                    nframe, trrheader.t);
        }

        sfree(x);
        sfree(v);
        sfree(f);
        nframe++;
    }
    if (!bOK)
    {
        fprintf(stderr, "\nWARNING: Incomplete frame header: nr %d, t=%g\n",
                nframe, trrheader.t);
    }
    gmx_trr_close(fpread);
}

static void list_xtc(const char *fn)
{
    t_fileio   *xd;
    int         indent;
    char        buf[256];
    rvec       *x;
    matrix      box;
    int         nframe, natoms;
    gmx_int64_t step;
    real        prec, time;
    gmx_bool    bOK;

    xd = open_xtc(fn, "r");
    read_first_xtc(xd, &natoms, &step, &time, box, &x, &prec, &bOK);

    nframe = 0;
    do
    {
        sprintf(buf, "%s frame %d", fn, nframe);
        indent = 0;
        indent = pr_title(stdout, indent, buf);
        pr_indent(stdout, indent);
        fprintf(stdout, "natoms=%10d  step=%10" GMX_PRId64 "  time=%12.7e  prec=%10g\n",
                natoms, step, time, prec);
        pr_rvecs(stdout, indent, "box", box, DIM);
        pr_rvecs(stdout, indent, "x", x, natoms);
        nframe++;
    }
    while (read_next_xtc(xd, natoms, &step, &time, box, x, &prec, &bOK));
    if (!bOK)
    {
        fprintf(stderr, "\nWARNING: Incomplete frame at time %g\n", time);
    }
    sfree(x);
    close_xtc(xd);
}

/*! \brief Callback used by list_tng_for_gmx_dump. */
static void list_tng_inner(const char *fn,
                           gmx_bool    bFirstFrame,
                           real       *values,
                           gmx_int64_t step,
                           double      frame_time,
                           gmx_int64_t n_values_per_frame,
                           gmx_int64_t n_atoms,
                           real        prec,
                           gmx_int64_t nframe,
                           char       *block_name)
{
    char                 buf[256];
    int                  indent = 0;

    if (bFirstFrame)
    {
        sprintf(buf, "%s frame %" GMX_PRId64, fn, nframe);
        indent = 0;
        indent = pr_title(stdout, indent, buf);
        pr_indent(stdout, indent);
        fprintf(stdout, "natoms=%10" GMX_PRId64 "  step=%10" GMX_PRId64 "  time=%12.7e",
                n_atoms, step, frame_time);
        if (prec > 0)
        {
            fprintf(stdout, "  prec=%10g", prec);
        }
        fprintf(stdout, "\n");
    }
    pr_reals_of_dim(stdout, indent, block_name, values, n_atoms, n_values_per_frame);
}

static void list_tng(const char gmx_unused *fn)
{
#ifdef GMX_USE_TNG
    gmx_tng_trajectory_t tng;
    gmx_int64_t          nframe = 0;
    gmx_int64_t          i, *block_ids = nullptr, step, ndatablocks;
    gmx_bool             bOK;
    real                *values = nullptr;

    gmx_tng_open(fn, 'r', &tng);
    gmx_print_tng_molecule_system(tng, stdout);

    bOK    = gmx_get_tng_data_block_types_of_next_frame(tng, -1,
                                                        0,
                                                        nullptr,
                                                        &step, &ndatablocks,
                                                        &block_ids);
    do
    {
        for (i = 0; i < ndatablocks; i++)
        {
            double               frame_time;
            real                 prec;
            gmx_int64_t          n_values_per_frame, n_atoms;
            char                 block_name[STRLEN];

            gmx_get_tng_data_next_frame_of_block_type(tng, block_ids[i], &values,
                                                      &step, &frame_time,
                                                      &n_values_per_frame, &n_atoms,
                                                      &prec,
                                                      block_name, STRLEN, &bOK);
            if (!bOK)
            {
                /* Can't write any output because we don't know what
                   arrays are valid. */
                fprintf(stderr, "\nWARNING: Incomplete frame at time %g, will not write output\n", frame_time);
            }
            else
            {
                list_tng_inner(fn, (0 == i), values, step, frame_time,
                               n_values_per_frame, n_atoms, prec, nframe, block_name);
            }
        }
        nframe++;
    }
    while (gmx_get_tng_data_block_types_of_next_frame(tng, step,
                                                      0,
                                                      nullptr,
                                                      &step,
                                                      &ndatablocks,
                                                      &block_ids));

    if (block_ids)
    {
        sfree(block_ids);
    }
    sfree(values);
    gmx_tng_close(&tng);
#endif
}

static void list_trx(const char *fn)
{
    switch (fn2ftp(fn))
    {
        case efXTC:
            list_xtc(fn);
            break;
        case efTRR:
            list_trr(fn);
            break;
        case efTNG:
            list_tng(fn);
            break;
        default:
            fprintf(stderr, "File %s is of an unsupported type. Try using the command\n 'less %s'\n",
                    fn, fn);
    }
}

static void list_ene(const char *fn)
{
    ener_file_t    in;
    gmx_bool       bCont;
    gmx_enxnm_t   *enm = nullptr;
    t_enxframe    *fr;
    int            i, j, nre, b;
    char           buf[22];

    printf("gmx dump: %s\n", fn);
    in = open_enx(fn, "r");
    do_enxnms(in, &nre, &enm);
    assert(enm);

    printf("energy components:\n");
    for (i = 0; (i < nre); i++)
    {
        printf("%5d  %-24s (%s)\n", i, enm[i].name, enm[i].unit);
    }

    snew(fr, 1);
    do
    {
        bCont = do_enx(in, fr);

        if (bCont)
        {
            printf("\n%24s  %12.5e  %12s  %12s\n", "time:",
                   fr->t, "step:", gmx_step_str(fr->step, buf));
            printf("%24s  %12s  %12s  %12s\n",
                   "", "", "nsteps:", gmx_step_str(fr->nsteps, buf));
            printf("%24s  %12.5e  %12s  %12s\n",
                   "delta_t:", fr->dt, "sum steps:", gmx_step_str(fr->nsum, buf));
            if (fr->nre == nre)
            {
                printf("%24s  %12s  %12s  %12s\n",
                       "Component", "Energy", "Av. Energy", "Sum Energy");
                if (fr->nsum > 0)
                {
                    for (i = 0; (i < nre); i++)
                    {
                        printf("%24s  %12.5e  %12.5e  %12.5e\n",
                               enm[i].name, fr->ener[i].e, fr->ener[i].eav,
                               fr->ener[i].esum);
                    }
                }
                else
                {
                    for (i = 0; (i < nre); i++)
                    {
                        printf("%24s  %12.5e\n",
                               enm[i].name, fr->ener[i].e);
                    }
                }
            }
            for (b = 0; b < fr->nblock; b++)
            {
                const char *typestr = "";

                t_enxblock *eb = &(fr->block[b]);
                printf("Block data %2d (%3d subblocks, id=%d)\n",
                       b, eb->nsub, eb->id);

                if (eb->id < enxNR)
                {
                    typestr = enx_block_id_name[eb->id];
                }
                printf("  id='%s'\n", typestr);
                for (i = 0; i < eb->nsub; i++)
                {
                    t_enxsubblock *sb = &(eb->sub[i]);
                    printf("  Sub block %3d (%5d elems, type=%s) values:\n",
                           i, sb->nr, xdr_datatype_names[sb->type]);

                    switch (sb->type)
                    {
                        case xdr_datatype_float:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d   %8.4f\n", j, sb->fval[j]);
                            }
                            break;
                        case xdr_datatype_double:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d   %10.6f\n", j, sb->dval[j]);
                            }
                            break;
                        case xdr_datatype_int:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d %10d\n", j, sb->ival[j]);
                            }
                            break;
                        case xdr_datatype_int64:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d %s\n",
                                       j, gmx_step_str(sb->lval[j], buf));
                            }
                            break;
                        case xdr_datatype_char:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d %1c\n", j, sb->cval[j]);
                            }
                            break;
                        case xdr_datatype_string:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d %80s\n", j, sb->sval[j]);
                            }
                            break;
                        default:
                            gmx_incons("Unknown subblock type");
                    }
                }
            }
        }
    }
    while (bCont);

    close_enx(in);

    free_enxframe(fr);
    sfree(fr);
    sfree(enm);
}

static void list_mtx(const char *fn)
{
    int                  nrow, ncol, i, j, k;
    real                *full   = nullptr, value;
    gmx_sparsematrix_t * sparse = nullptr;

    gmx_mtxio_read(fn, &nrow, &ncol, &full, &sparse);

    if (full == nullptr)
    {
        snew(full, nrow*ncol);
        for (i = 0; i < nrow*ncol; i++)
        {
            full[i] = 0;
        }

        for (i = 0; i < sparse->nrow; i++)
        {
            for (j = 0; j < sparse->ndata[i]; j++)
            {
                k              = sparse->data[i][j].col;
                value          = sparse->data[i][j].value;
                full[i*ncol+k] = value;
                full[k*ncol+i] = value;
            }
        }
        gmx_sparsematrix_destroy(sparse);
    }

    printf("%d %d\n", nrow, ncol);
    for (i = 0; i < nrow; i++)
    {
        for (j = 0; j < ncol; j++)
        {
            printf(" %g", full[i*ncol+j]);
        }
        printf("\n");
    }

    sfree(full);
}

int gmx_dump(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] reads a run input file ([REF].tpr[ref]),",
        "a trajectory ([REF].trr[ref]/[REF].xtc[ref]/[TT]/tng[tt]), an energy",
        "file ([REF].edr[ref]) or a checkpoint file ([REF].cpt[ref])",
        "and prints that to standard output in a readable format.",
        "This program is essential for checking your run input file in case of",
        "problems.[PAR]",
        "The program can also preprocess a topology to help finding problems.",
        "Note that currently setting [TT]GMXLIB[tt] is the only way to customize",
        "directories used for searching include files.",
    };
    const char *bugs[] = {
        "Position restraint output from -sys -s is broken"
    };
    t_filenm    fnm[] = {
        { efTPR, "-s", nullptr, ffOPTRD },
        { efTRX, "-f", nullptr, ffOPTRD },
        { efEDR, "-e", nullptr, ffOPTRD },
        { efCPT, nullptr, nullptr, ffOPTRD },
        { efTOP, "-p", nullptr, ffOPTRD },
        { efMTX, "-mtx", "hessian", ffOPTRD },
        { efMDP, "-om", nullptr, ffOPTWR }
    };
#define NFILE asize(fnm)

    gmx_output_env_t *oenv;
    /* Command line options */
    gmx_bool          bShowNumbers      = TRUE;
    gmx_bool          bShowParams       = FALSE;
    gmx_bool          bSysTop           = FALSE;
    gmx_bool          bOriginalInputrec = FALSE;
    t_pargs           pa[]              = {
        { "-nr", FALSE, etBOOL, {&bShowNumbers}, "Show index numbers in output (leaving them out makes comparison easier, but creates a useless topology)" },
        { "-param", FALSE, etBOOL, {&bShowParams}, "Show parameters for each bonded interaction (for comparing dumps, it is useful to combine this with -nonr)" },
        { "-sys", FALSE, etBOOL, {&bSysTop}, "List the atoms and bonded interactions for the whole system instead of for each molecule type" },
        { "-orgir", FALSE, etBOOL, {&bOriginalInputrec}, "Show input parameters from tpr as they were written by the version that produced the file, instead of how the current version reads them" }
    };

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }


    if (ftp2bSet(efTPR, NFILE, fnm))
    {
        list_tpx(ftp2fn(efTPR, NFILE, fnm), bShowNumbers, bShowParams,
                 ftp2fn_null(efMDP, NFILE, fnm), bSysTop, bOriginalInputrec);
    }
    else if (ftp2bSet(efTRX, NFILE, fnm))
    {
        list_trx(ftp2fn(efTRX, NFILE, fnm));
    }
    else if (ftp2bSet(efEDR, NFILE, fnm))
    {
        list_ene(ftp2fn(efEDR, NFILE, fnm));
    }
    else if (ftp2bSet(efCPT, NFILE, fnm))
    {
        list_checkpoint(ftp2fn(efCPT, NFILE, fnm), stdout);
    }
    else if (ftp2bSet(efTOP, NFILE, fnm))
    {
        list_top(ftp2fn(efTOP, NFILE, fnm));
    }
    else if (ftp2bSet(efMTX, NFILE, fnm))
    {
        list_mtx(ftp2fn(efMTX, NFILE, fnm));
    }

    return 0;
}
