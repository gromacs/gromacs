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
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    char *label;
    int   cg;
    real  q;
} t_charge;

static t_charge *mk_charge(t_atoms *atoms, t_block *cgs, int *nncg)
{
    t_charge *cg = NULL;
    char      buf[32];
    int       i, j, ncg, resnr, anr;
    real      qq;

    /* Find the charged groups */
    ncg = 0;
    for (i = 0; (i < cgs->nr); i++)
    {
        qq = 0.0;
        for (j = cgs->index[i]; (j < cgs->index[i+1]); j++)
        {
            qq += atoms->atom[j].q;
        }
        if (fabs(qq) > 1.0e-5)
        {
            srenew(cg, ncg+1);
            cg[ncg].q  = qq;
            cg[ncg].cg = i;
            anr        = cgs->index[i];
            resnr      = atoms->atom[anr].resind;
            sprintf(buf, "%s%d-%d",
                    *(atoms->resinfo[resnr].name),
                    atoms->resinfo[resnr].nr,
                    anr+1);
            cg[ncg].label = gmx_strdup(buf);
            ncg++;
        }
    }
    *nncg = ncg;

    for (i = 0; (i < ncg); i++)
    {
        printf("CG: %10s Q: %6g  Atoms:",
               cg[i].label, cg[i].q);
        for (j = cgs->index[cg[i].cg]; (j < cgs->index[cg[i].cg+1]); j++)
        {
            printf(" %4d", j);
        }
        printf("\n");
    }

    return cg;
}

static real calc_dist(t_pbc *pbc, rvec x[], t_block *cgs, int icg, int jcg)
{
    int  i, j;
    rvec dx;
    real d2, mindist2 = 1000;

    for (i = cgs->index[icg]; (i < cgs->index[icg+1]); i++)
    {
        for (j = cgs->index[jcg]; (j < cgs->index[jcg+1]); j++)
        {
            pbc_dx(pbc, x[i], x[j], dx);
            d2 = norm2(dx);
            if (d2 < mindist2)
            {
                mindist2 = d2;
            }
        }
    }
    return sqrt(mindist2);
}

int gmx_saltbr(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] plots the distance between all combination of charged groups",
        "as a function of time. The groups are combined in different ways.",
        "A minimum distance can be given (i.e. a cut-off), such that groups",
        "that are never closer than that distance will not be plotted.[PAR]",
        "Output will be in a number of fixed filenames, [TT]min-min.xvg[tt], [TT]plus-min.xvg[tt]",
        "and [TT]plus-plus.xvg[tt], or files for every individual ion pair if the [TT]-sep[tt]",
        "option is selected. In this case, files are named as [TT]sb-(Resname)(Resnr)-(Atomnr)[tt].",
        "There may be [BB]many[bb] such files."
    };
    static gmx_bool bSep     = FALSE;
    static real     truncate = 1000.0;
    t_pargs         pa[]     = {
        { "-t",   FALSE, etREAL, {&truncate},
          "Groups that are never closer than this distance are not plotted" },
        { "-sep", FALSE, etBOOL, {&bSep},
          "Use separate files for each interaction (may be MANY)" }
    };
    t_filenm        fnm[] = {
        { efTRX, "-f",  NULL, ffREAD },
        { efTPR, NULL,  NULL, ffREAD },
    };
#define NFILE asize(fnm)

    FILE              *out[3], *fp;
    static const char *title[3] = {
        "Distance between positively charged groups",
        "Distance between negatively charged groups",
        "Distance between oppositely charged groups"
    };
    static const char *fn[3] = {
        "plus-plus.xvg",
        "min-min.xvg",
        "plus-min.xvg"
    };
    int                nset[3] = {0, 0, 0};

    t_topology        *top;
    int                ePBC;
    char              *buf;
    t_trxstatus       *status;
    int                i, j, k, m, nnn, teller, ncg, n1, n2, n3, natoms;
    real               t, *time, qi, qj;
    t_charge          *cg;
    real            ***cgdist;
    int              **nWithin;

    double             t0, dt;
    char               label[234];
    t_pbc              pbc;
    rvec              *x;
    matrix             box;
    output_env_t       oenv;

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC);
    cg  = mk_charge(&top->atoms, &(top->cgs), &ncg);
    snew(cgdist, ncg);
    snew(nWithin, ncg);
    for (i = 0; (i < ncg); i++)
    {
        snew(cgdist[i], ncg);
        snew(nWithin[i], ncg);
    }

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    teller = 0;
    time   = NULL;
    do
    {
        srenew(time, teller+1);
        time[teller] = t;

        set_pbc(&pbc, ePBC, box);

        for (i = 0; (i < ncg); i++)
        {
            for (j = i+1; (j < ncg); j++)
            {
                srenew(cgdist[i][j], teller+1);
                cgdist[i][j][teller] =
                    calc_dist(&pbc, x, &(top->cgs), cg[i].cg, cg[j].cg);
                if (cgdist[i][j][teller] < truncate)
                {
                    nWithin[i][j] = 1;
                }
            }
        }

        teller++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "\n");
    close_trj(status);

    if (bSep)
    {
        snew(buf, 256);
        for (i = 0; (i < ncg); i++)
        {
            for (j = i+1; (j < ncg); j++)
            {
                if (nWithin[i][j])
                {
                    sprintf(buf, "sb-%s:%s.xvg", cg[i].label, cg[j].label);
                    fp = xvgropen(buf, buf, "Time (ps)", "Distance (nm)", oenv);
                    for (k = 0; (k < teller); k++)
                    {
                        fprintf(fp, "%10g  %10g\n", time[k], cgdist[i][j][k]);
                    }
                    xvgrclose(fp);
                }
            }
        }
        sfree(buf);
    }
    else
    {

        for (m = 0; (m < 3); m++)
        {
            out[m] = xvgropen(fn[m], title[m], "Time (ps)", "Distance (nm)", oenv);
        }

        snew(buf, 256);
        for (i = 0; (i < ncg); i++)
        {
            qi = cg[i].q;
            for (j = i+1; (j < ncg); j++)
            {
                qj = cg[j].q;
                if (nWithin[i][j])
                {
                    sprintf(buf, "%s:%s", cg[i].label, cg[j].label);
                    if (qi*qj < 0)
                    {
                        nnn = 2;
                    }
                    else if (qi+qj > 0)
                    {
                        nnn = 0;
                    }
                    else
                    {
                        nnn = 1;
                    }

                    if (nset[nnn] == 0)
                    {
                        xvgr_legend(out[nnn], 1, (const char**)&buf, oenv);
                    }
                    else
                    {
                        if (output_env_get_xvg_format(oenv) == exvgXMGR)
                        {
                            fprintf(out[nnn], "@ legend string %d \"%s\"\n", nset[nnn], buf);
                        }
                        else if (output_env_get_xvg_format(oenv) == exvgXMGRACE)
                        {
                            fprintf(out[nnn], "@ s%d legend \"%s\"\n", nset[nnn], buf);
                        }
                    }
                    nset[nnn]++;
                    nWithin[i][j] = nnn+1;
                }
            }
        }
        for (k = 0; (k < teller); k++)
        {
            for (m = 0; (m < 3); m++)
            {
                fprintf(out[m], "%10g", time[k]);
            }

            for (i = 0; (i < ncg); i++)
            {
                for (j = i+1; (j < ncg); j++)
                {
                    nnn = nWithin[i][j];
                    if (nnn > 0)
                    {
                        fprintf(out[nnn-1], "  %10g", cgdist[i][j][k]);
                    }
                }
            }
            for (m = 0; (m < 3); m++)
            {
                fprintf(out[m], "\n");
            }
        }
        for (m = 0; (m < 3); m++)
        {
            xvgrclose(out[m]);
            if (nset[m] == 0)
            {
                remove(fn[m]);
            }
        }
    }

    return 0;
}
