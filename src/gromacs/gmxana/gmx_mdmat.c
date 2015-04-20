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
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


#define FARAWAY 10000

int *res_ndx(t_atoms *atoms)
{
    int *rndx;
    int  i, r0;

    if (atoms->nr <= 0)
    {
        return NULL;
    }
    snew(rndx, atoms->nr);
    r0 = atoms->atom[0].resind;
    for (i = 0; (i < atoms->nr); i++)
    {
        rndx[i] = atoms->atom[i].resind-r0;
    }

    return rndx;
}

int *res_natm(t_atoms *atoms)
{
    int *natm;
    int  i, j, r0;

    if (atoms->nr <= 0)
    {
        return NULL;
    }
    snew(natm, atoms->nres);
    r0 = atoms->atom[0].resind;
    j  = 0;
    for (i = 0; (i < atoms->nres); i++)
    {
        while ((atoms->atom[j].resind)-r0 == i)
        {
            natm[i]++;
            j++;
        }
    }

    return natm;
}

static void calc_mat(int nres, int natoms, int rndx[],
                     rvec x[], atom_id *index,
                     real trunc, real **mdmat, int **nmat, int ePBC, matrix box)
{
    int   i, j, resi, resj;
    real  trunc2, r, r2;
    t_pbc pbc;
    rvec  ddx;

    set_pbc(&pbc, ePBC, box);
    trunc2 = sqr(trunc);
    for (resi = 0; (resi < nres); resi++)
    {
        for (resj = 0; (resj < nres); resj++)
        {
            mdmat[resi][resj] = FARAWAY;
        }
    }
    for (i = 0; (i < natoms); i++)
    {
        resi = rndx[i];
        for (j = i+1; (j < natoms); j++)
        {
            resj = rndx[j];
            pbc_dx(&pbc, x[index[i]], x[index[j]], ddx);
            r2 = norm2(ddx);
            if (r2 < trunc2)
            {
                nmat[resi][j]++;
                nmat[resj][i]++;
            }
            mdmat[resi][resj] = min(r2, mdmat[resi][resj]);
        }
    }

    for (resi = 0; (resi < nres); resi++)
    {
        mdmat[resi][resi] = 0;
        for (resj = resi+1; (resj < nres); resj++)
        {
            r                 = sqrt(mdmat[resi][resj]);
            mdmat[resi][resj] = r;
            mdmat[resj][resi] = r;
        }
    }
}

static void tot_nmat(int nres, int natoms, int nframes, int **nmat,
                     int *tot_n, real *mean_n)
{
    int i, j;

    for (i = 0; (i < nres); i++)
    {
        for (j = 0; (j < natoms); j++)
        {
            if (nmat[i][j] != 0)
            {
                tot_n[i]++;
                mean_n[i] += nmat[i][j];
            }
        }
        mean_n[i] /= nframes;
    }
}

int gmx_mdmat(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] makes distance matrices consisting of the smallest distance",
        "between residue pairs. With [TT]-frames[tt], these distance matrices can be",
        "stored in order to see differences in tertiary structure as a",
        "function of time. If you choose your options unwisely, this may generate",
        "a large output file. By default, only an averaged matrix over the whole",
        "trajectory is output.",
        "Also a count of the number of different atomic contacts between",
        "residues over the whole trajectory can be made.",
        "The output can be processed with [gmx-xpm2ps] to make a PostScript (tm) plot."
    };
    static real     truncate = 1.5;
    static gmx_bool bAtom    = FALSE;
    static int      nlevels  = 40;
    t_pargs         pa[]     = {
        { "-t",   FALSE, etREAL, {&truncate},
          "trunc distance" },
        { "-nlevels",   FALSE, etINT,  {&nlevels},
          "Discretize distance in this number of levels" }
    };
    t_filenm        fnm[] = {
        { efTRX, "-f",  NULL, ffREAD },
        { efTPS, NULL,  NULL, ffREAD },
        { efNDX, NULL,  NULL, ffOPTRD },
        { efXPM, "-mean", "dm", ffWRITE },
        { efXPM, "-frames", "dmf", ffOPTWR },
        { efXVG, "-no", "num", ffOPTWR },
    };
#define NFILE asize(fnm)

    FILE          *out = NULL, *fp;
    t_topology     top;
    int            ePBC;
    t_atoms        useatoms;
    int            isize;
    atom_id       *index;
    char          *grpname;
    int           *rndx, *natm, prevres, newres;

    int            i, j, nres, natoms, nframes, it, trxnat;
    t_trxstatus   *status;
    int            nr0;
    gmx_bool       bCalcN, bFrames;
    real           t, ratio;
    char           title[256], label[234];
    t_rgb          rlo, rhi;
    rvec          *x;
    real         **mdmat, *resnr, **totmdmat;
    int          **nmat, **totnmat;
    real          *mean_n;
    int           *tot_n;
    matrix         box = {{0}};
    output_env_t   oenv;
    gmx_rmpbc_t    gpbc = NULL;

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME, NFILE, fnm,
                           asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    fprintf(stderr, "Will truncate at %f nm\n", truncate);
    bCalcN  = opt2bSet("-no", NFILE, fnm);
    bFrames = opt2bSet("-frames", NFILE, fnm);
    if (bCalcN)
    {
        fprintf(stderr, "Will calculate number of different contacts\n");
    }

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &x, NULL, box, FALSE);

    fprintf(stderr, "Select group for analysis\n");
    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);

    natoms = isize;
    snew(useatoms.atom, natoms);
    snew(useatoms.atomname, natoms);

    useatoms.nres = 0;
    snew(useatoms.resinfo, natoms);

    prevres = top.atoms.atom[index[0]].resind;
    newres  = 0;
    for (i = 0; (i < isize); i++)
    {
        int ii = index[i];
        useatoms.atomname[i] = top.atoms.atomname[ii];
        if (top.atoms.atom[ii].resind != prevres)
        {
            prevres = top.atoms.atom[ii].resind;
            newres++;
            useatoms.resinfo[i] = top.atoms.resinfo[prevres];
            if (debug)
            {
                fprintf(debug, "New residue: atom %5s %5s %6d, index entry %5d, newres %5d\n",
                        *(top.atoms.resinfo[top.atoms.atom[ii].resind].name),
                        *(top.atoms.atomname[ii]),
                        ii, i, newres);
            }
        }
        useatoms.atom[i].resind = newres;
    }
    useatoms.nres = newres+1;
    useatoms.nr   = isize;

    rndx = res_ndx(&(useatoms));
    natm = res_natm(&(useatoms));
    nres = useatoms.nres;
    fprintf(stderr, "There are %d residues with %d atoms\n", nres, natoms);

    snew(resnr, nres);
    snew(mdmat, nres);
    snew(nmat, nres);
    snew(totnmat, nres);
    snew(mean_n, nres);
    snew(tot_n, nres);
    for (i = 0; (i < nres); i++)
    {
        snew(mdmat[i], nres);
        snew(nmat[i], natoms);
        snew(totnmat[i], natoms);
        resnr[i] = i+1;
    }
    snew(totmdmat, nres);
    for (i = 0; (i < nres); i++)
    {
        snew(totmdmat[i], nres);
    }

    trxnat = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    nframes = 0;

    rlo.r = 1.0, rlo.g = 1.0, rlo.b = 1.0;
    rhi.r = 0.0, rhi.g = 0.0, rhi.b = 0.0;

    gpbc = gmx_rmpbc_init(&top.idef, ePBC, trxnat);

    if (bFrames)
    {
        out = opt2FILE("-frames", NFILE, fnm, "w");
    }
    do
    {
        gmx_rmpbc(gpbc, trxnat, box, x);
        nframes++;
        calc_mat(nres, natoms, rndx, x, index, truncate, mdmat, nmat, ePBC, box);
        for (i = 0; (i < nres); i++)
        {
            for (j = 0; (j < natoms); j++)
            {
                if (nmat[i][j])
                {
                    totnmat[i][j]++;
                }
            }
        }
        for (i = 0; (i < nres); i++)
        {
            for (j = 0; (j < nres); j++)
            {
                totmdmat[i][j] += mdmat[i][j];
            }
        }
        if (bFrames)
        {
            sprintf(label, "t=%.0f ps", t);
            write_xpm(out, 0, label, "Distance (nm)", "Residue Index", "Residue Index",
                      nres, nres, resnr, resnr, mdmat, 0, truncate, rlo, rhi, &nlevels);
        }
    }
    while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "\n");
    close_trj(status);
    gmx_rmpbc_done(gpbc);
    if (bFrames)
    {
        gmx_ffclose(out);
    }

    fprintf(stderr, "Processed %d frames\n", nframes);

    for (i = 0; (i < nres); i++)
    {
        for (j = 0; (j < nres); j++)
        {
            totmdmat[i][j] /= nframes;
        }
    }
    write_xpm(opt2FILE("-mean", NFILE, fnm, "w"), 0, "Mean smallest distance",
              "Distance (nm)", "Residue Index", "Residue Index",
              nres, nres, resnr, resnr, totmdmat, 0, truncate, rlo, rhi, &nlevels);

    if (bCalcN)
    {
        char **legend;

        snew(legend, 5);
        for (i = 0; i < 5; i++)
        {
            snew(legend[i], STRLEN);
        }
        tot_nmat(nres, natoms, nframes, totnmat, tot_n, mean_n);
        fp = xvgropen(ftp2fn(efXVG, NFILE, fnm),
                      "Increase in number of contacts", "Residue", "Ratio", oenv);
        sprintf(legend[0], "Total/mean");
        sprintf(legend[1], "Total");
        sprintf(legend[2], "Mean");
        sprintf(legend[3], "# atoms");
        sprintf(legend[4], "Mean/# atoms");
        xvgr_legend(fp, 5, (const char**)legend, oenv);
        for (i = 0; (i < nres); i++)
        {
            if (mean_n[i] == 0)
            {
                ratio = 1;
            }
            else
            {
                ratio = tot_n[i]/mean_n[i];
            }
            fprintf(fp, "%3d  %8.3f  %3d  %8.3f  %3d  %8.3f\n",
                    i+1, ratio, tot_n[i], mean_n[i], natm[i], mean_n[i]/natm[i]);
        }
        xvgrclose(fp);
    }

    return 0;
}
