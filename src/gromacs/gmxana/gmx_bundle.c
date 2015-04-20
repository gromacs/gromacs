/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include <math.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#define MAX_ENDS 3

typedef struct {
    int   n;
    int   nend;
    rvec *end[MAX_ENDS];
    rvec *mid;
    rvec *dir;
    real *len;
} t_bundle;

static void rotate_ends(t_bundle *bun, rvec axis, int c0, int c1)
{
    int  end, i;
    rvec ax, tmp;

    unitv(axis, ax);
    for (end = 0; end < bun->nend; end++)
    {
        for (i = 0; i < bun->n; i++)
        {
            copy_rvec(bun->end[end][i], tmp);
            bun->end[end][i][c0] = ax[c1]*tmp[c0] - ax[c0]*tmp[c1];
            bun->end[end][i][c1] = ax[c0]*tmp[c0] + ax[c1]*tmp[c1];
        }
    }
    copy_rvec(axis, tmp);
    axis[c0] = ax[c1]*tmp[c0] - ax[c0]*tmp[c1];
    axis[c1] = ax[c0]*tmp[c0] + ax[c1]*tmp[c1];
}

static void calc_axes(rvec x[], t_atom atom[], int gnx[], atom_id *index[],
                      gmx_bool bRot, t_bundle *bun)
{
    int   end, i, div, d;
    real *mtot, m;
    rvec  axis[MAX_ENDS], cent;

    snew(mtot, bun->n);

    for (end = 0; end < bun->nend; end++)
    {
        for (i = 0; i < bun->n; i++)
        {
            clear_rvec(bun->end[end][i]);
            mtot[i] = 0;
        }
        div = gnx[end]/bun->n;
        for (i = 0; i < gnx[end]; i++)
        {
            m = atom[index[end][i]].m;
            for (d = 0; d < DIM; d++)
            {
                bun->end[end][i/div][d] += m*x[index[end][i]][d];
            }
            mtot[i/div] += m;
        }
        clear_rvec(axis[end]);
        for (i = 0; i < bun->n; i++)
        {
            svmul(1.0/mtot[i], bun->end[end][i], bun->end[end][i]);
            rvec_inc(axis[end], bun->end[end][i]);
        }
        svmul(1.0/bun->n, axis[end], axis[end]);
    }
    sfree(mtot);

    rvec_add(axis[0], axis[1], cent);
    svmul(0.5, cent, cent);
    /* center the bundle on the origin */
    for (end = 0; end < bun->nend; end++)
    {
        rvec_dec(axis[end], cent);
        for (i = 0; i < bun->n; i++)
        {
            rvec_dec(bun->end[end][i], cent);
        }
    }
    if (bRot)
    {
        /* rotate the axis parallel to the z-axis */
        rotate_ends(bun, axis[0], YY, ZZ);
        rotate_ends(bun, axis[0], XX, ZZ);
    }
    for (i = 0; i < bun->n; i++)
    {
        rvec_add(bun->end[0][i], bun->end[1][i], bun->mid[i]);
        svmul(0.5, bun->mid[i], bun->mid[i]);
        rvec_sub(bun->end[0][i], bun->end[1][i], bun->dir[i]);
        bun->len[i] = norm(bun->dir[i]);
        unitv(bun->dir[i], bun->dir[i]);
    }
}

static void dump_axes(t_trxstatus *status, t_trxframe *fr, t_atoms *outat,
                      t_bundle *bun)
{
    t_trxframe   frout;
    static rvec *xout = NULL;
    int          i, end;

    if (xout == NULL)
    {
        snew(xout, outat->nr);
    }

    for (i = 0; i < bun->n; i++)
    {
        copy_rvec(bun->end[0][i], xout[3*i]);
        if (bun->nend >= 3)
        {
            copy_rvec(bun->end[2][i], xout[3*i+1]);
        }
        else
        {
            copy_rvec(bun->mid[i], xout[3*i+1]);
        }
        copy_rvec(bun->end[1][i], xout[3*i+2]);
    }
    frout        = *fr;
    frout.bV     = FALSE;
    frout.bF     = FALSE;
    frout.bBox   = FALSE;
    frout.bAtoms = TRUE;
    frout.natoms = outat->nr;
    frout.atoms  = outat;
    frout.x      = xout;
    write_trxframe(status, &frout, NULL);
}

int gmx_bundle(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] analyzes bundles of axes. The axes can be for instance",
        "helix axes. The program reads two index groups and divides both",
        "of them in [TT]-na[tt] parts. The centers of mass of these parts",
        "define the tops and bottoms of the axes.",
        "Several quantities are written to file:",
        "the axis length, the distance and the z-shift of the axis mid-points",
        "with respect to the average center of all axes, the total tilt,",
        "the radial tilt and the lateral tilt with respect to the average axis.",
        "[PAR]",
        "With options [TT]-ok[tt], [TT]-okr[tt] and [TT]-okl[tt] the total,",
        "radial and lateral kinks of the axes are plotted. An extra index",
        "group of kink atoms is required, which is also divided into [TT]-na[tt]",
        "parts. The kink angle is defined as the angle between the kink-top and",
        "the bottom-kink vectors.",
        "[PAR]",
        "With option [TT]-oa[tt] the top, mid (or kink when [TT]-ok[tt] is set)",
        "and bottom points of each axis",
        "are written to a [REF].pdb[ref] file each frame. The residue numbers correspond",
        "to the axis numbers. When viewing this file with Rasmol, use the",
        "command line option [TT]-nmrpdb[tt], and type [TT]set axis true[tt] to",
        "display the reference axis."
    };
    static int      n    = 0;
    static gmx_bool bZ   = FALSE;
    t_pargs         pa[] = {
        { "-na", FALSE, etINT, {&n},
          "Number of axes" },
        { "-z", FALSE, etBOOL, {&bZ},
          "Use the [IT]z[it]-axis as reference instead of the average axis" }
    };
    FILE           *out, *flen, *fdist, *fz, *ftilt, *ftiltr, *ftiltl;
    FILE           *fkink = NULL, *fkinkr = NULL, *fkinkl = NULL;
    t_trxstatus    *status;
    t_trxstatus    *fpdb;
    t_topology      top;
    int             ePBC;
    rvec           *xtop;
    matrix          box;
    t_trxframe      fr;
    t_atoms         outatoms;
    real            t, comp;
    int             natoms;
    char           *grpname[MAX_ENDS], title[256];
    /* FIXME: The constness should not be cast away */
    char           *anm = (char *)"CA", *rnm = (char *)"GLY";
    int             i, j, gnx[MAX_ENDS];
    atom_id        *index[MAX_ENDS];
    t_bundle        bun;
    gmx_bool        bKink;
    rvec            va, vb, vc, vr, vl;
    output_env_t    oenv;
    gmx_rmpbc_t     gpbc = NULL;

#define NLEG asize(leg)
    t_filenm fnm[] = {
        { efTRX, "-f", NULL, ffREAD },
        { efTPS, NULL, NULL, ffREAD },
        { efNDX, NULL, NULL, ffOPTRD },
        { efXVG, "-ol", "bun_len", ffWRITE },
        { efXVG, "-od", "bun_dist", ffWRITE },
        { efXVG, "-oz", "bun_z", ffWRITE },
        { efXVG, "-ot", "bun_tilt", ffWRITE },
        { efXVG, "-otr", "bun_tiltr", ffWRITE },
        { efXVG, "-otl", "bun_tiltl", ffWRITE },
        { efXVG, "-ok", "bun_kink", ffOPTWR },
        { efXVG, "-okr", "bun_kinkr", ffOPTWR },
        { efXVG, "-okl", "bun_kinkl", ffOPTWR },
        { efPDB, "-oa", "axes", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &xtop, NULL, box, TRUE);

    bKink = opt2bSet("-ok", NFILE, fnm) || opt2bSet("-okr", NFILE, fnm)
        || opt2bSet("-okl", NFILE, fnm);
    if (bKink)
    {
        bun.nend = 3;
    }
    else
    {
        bun.nend = 2;
    }

    fprintf(stderr, "Select a group of top and a group of bottom ");
    if (bKink)
    {
        fprintf(stderr, "and a group of kink ");
    }
    fprintf(stderr, "atoms\n");
    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), bun.nend,
              gnx, index, grpname);

    if (n <= 0 || gnx[0] % n || gnx[1] % n || (bKink && gnx[2] % n))
    {
        gmx_fatal(FARGS,
                  "The size of one of your index groups is not a multiple of n");
    }
    bun.n = n;
    snew(bun.end[0], n);
    snew(bun.end[1], n);
    if (bKink)
    {
        snew(bun.end[2], n);
    }
    snew(bun.mid, n);
    snew(bun.dir, n);
    snew(bun.len, n);

    flen   = xvgropen(opt2fn("-ol", NFILE, fnm), "Axis lengths",
                      output_env_get_xvgr_tlabel(oenv), "(nm)", oenv);
    fdist  = xvgropen(opt2fn("-od", NFILE, fnm), "Distance of axis centers",
                      output_env_get_xvgr_tlabel(oenv), "(nm)", oenv);
    fz     = xvgropen(opt2fn("-oz", NFILE, fnm), "Z-shift of axis centers",
                      output_env_get_xvgr_tlabel(oenv), "(nm)", oenv);
    ftilt  = xvgropen(opt2fn("-ot", NFILE, fnm), "Axis tilts",
                      output_env_get_xvgr_tlabel(oenv), "(degrees)", oenv);
    ftiltr = xvgropen(opt2fn("-otr", NFILE, fnm), "Radial axis tilts",
                      output_env_get_xvgr_tlabel(oenv), "(degrees)", oenv);
    ftiltl = xvgropen(opt2fn("-otl", NFILE, fnm), "Lateral axis tilts",
                      output_env_get_xvgr_tlabel(oenv), "(degrees)", oenv);

    if (bKink)
    {
        fkink  = xvgropen(opt2fn("-ok", NFILE, fnm), "Kink angles",
                          output_env_get_xvgr_tlabel(oenv), "(degrees)", oenv);
        fkinkr = xvgropen(opt2fn("-okr", NFILE, fnm), "Radial kink angles",
                          output_env_get_xvgr_tlabel(oenv), "(degrees)", oenv);
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(fkinkr, "@ subtitle \"+ = ) (   - = ( )\"\n");
        }
        fkinkl = xvgropen(opt2fn("-okl", NFILE, fnm), "Lateral kink angles",
                          output_env_get_xvgr_tlabel(oenv), "(degrees)", oenv);
    }

    if (opt2bSet("-oa", NFILE, fnm))
    {
        init_t_atoms(&outatoms, 3*n, FALSE);
        outatoms.nr = 3*n;
        for (i = 0; i < 3*n; i++)
        {
            outatoms.atomname[i]       = &anm;
            outatoms.atom[i].resind    = i/3;
            outatoms.resinfo[i/3].name = &rnm;
            outatoms.resinfo[i/3].nr   = i/3 + 1;
            outatoms.resinfo[i/3].ic   = ' ';
        }
        fpdb = open_trx(opt2fn("-oa", NFILE, fnm), "w");
    }
    else
    {
        fpdb = NULL;
    }

    read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, TRX_NEED_X);
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, fr.natoms);

    do
    {
        gmx_rmpbc_trxfr(gpbc, &fr);
        calc_axes(fr.x, top.atoms.atom, gnx, index, !bZ, &bun);
        t = output_env_conv_time(oenv, fr.time);
        fprintf(flen, " %10g", t);
        fprintf(fdist, " %10g", t);
        fprintf(fz, " %10g", t);
        fprintf(ftilt, " %10g", t);
        fprintf(ftiltr, " %10g", t);
        fprintf(ftiltl, " %10g", t);
        if (bKink)
        {
            fprintf(fkink, " %10g", t);
            fprintf(fkinkr, " %10g", t);
            fprintf(fkinkl, " %10g", t);
        }

        for (i = 0; i < bun.n; i++)
        {
            fprintf(flen, " %6g", bun.len[i]);
            fprintf(fdist, " %6g", norm(bun.mid[i]));
            fprintf(fz, " %6g", bun.mid[i][ZZ]);
            fprintf(ftilt, " %6g", RAD2DEG*acos(bun.dir[i][ZZ]));
            comp = bun.mid[i][XX]*bun.dir[i][XX]+bun.mid[i][YY]*bun.dir[i][YY];
            fprintf(ftiltr, " %6g", RAD2DEG*
                    asin(comp/sqrt(sqr(comp)+sqr(bun.dir[i][ZZ]))));
            comp = bun.mid[i][YY]*bun.dir[i][XX]-bun.mid[i][XX]*bun.dir[i][YY];
            fprintf(ftiltl, " %6g", RAD2DEG*
                    asin(comp/sqrt(sqr(comp)+sqr(bun.dir[i][ZZ]))));
            if (bKink)
            {
                rvec_sub(bun.end[0][i], bun.end[2][i], va);
                rvec_sub(bun.end[2][i], bun.end[1][i], vb);
                unitv_no_table(va, va);
                unitv_no_table(vb, vb);
                fprintf(fkink, " %6g", RAD2DEG*acos(iprod(va, vb)));
                cprod(va, vb, vc);
                copy_rvec(bun.mid[i], vr);
                vr[ZZ] = 0;
                unitv_no_table(vr, vr);
                fprintf(fkinkr, " %6g", RAD2DEG*asin(iprod(vc, vr)));
                vl[XX] = vr[YY];
                vl[YY] = -vr[XX];
                vl[ZZ] = 0;
                fprintf(fkinkl, " %6g", RAD2DEG*asin(iprod(vc, vl)));
            }
        }
        fprintf(flen, "\n");
        fprintf(fdist, "\n");
        fprintf(fz, "\n");
        fprintf(ftilt, "\n");
        fprintf(ftiltr, "\n");
        fprintf(ftiltl, "\n");
        if (bKink)
        {
            fprintf(fkink, "\n");
            fprintf(fkinkr, "\n");
            fprintf(fkinkl, "\n");
        }
        if (fpdb)
        {
            dump_axes(fpdb, &fr, &outatoms, &bun);
        }
    }
    while (read_next_frame(oenv, status, &fr));
    gmx_rmpbc_done(gpbc);

    close_trx(status);

    if (fpdb)
    {
        close_trx(fpdb);
    }
    xvgrclose(flen);
    xvgrclose(fdist);
    xvgrclose(fz);
    xvgrclose(ftilt);
    xvgrclose(ftiltr);
    xvgrclose(ftiltl);
    if (bKink)
    {
        xvgrclose(fkink);
        xvgrclose(fkinkr);
        xvgrclose(fkinkl);
    }

    return 0;
}
