/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

enum class PbcType : int;
struct gmx_output_env_t;

typedef struct
{
    char* label;
    int   cg;
    real  q;
} t_charge;

static t_charge* mk_charge(const t_atoms* atoms, int* nncg)
{
    t_charge* cg = nullptr;
    char      buf[32];
    int       i, ncg, resnr, anr;
    real      qq;

    /* Find the charged groups */
    ncg = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        qq = atoms->atom[i].q;
        if (std::abs(qq) > 1.0e-5)
        {
            srenew(cg, ncg + 1);
            cg[ncg].q  = qq;
            cg[ncg].cg = i;
            anr        = i;
            resnr      = atoms->atom[anr].resind;
            sprintf(buf, "%s%d-%d", *(atoms->resinfo[resnr].name), atoms->resinfo[resnr].nr, anr + 1);
            cg[ncg].label = gmx_strdup(buf);
            ncg++;
        }
    }
    *nncg = ncg;

    for (i = 0; (i < ncg); i++)
    {
        printf("CG: %10s Q: %6g  Atoms:", cg[i].label, cg[i].q);
        printf(" %4d", cg[i].cg);
        printf("\n");
    }

    return cg;
}

int gmx_saltbr(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] plots the distance between all combination of charged groups",
        "as a function of time. The groups are combined in different ways.",
        "A minimum distance can be given (i.e. a cut-off), such that groups",
        "that are never closer than that distance will not be plotted.[PAR]",
        "Output will be in a number of fixed filenames, [TT]min-min.xvg[tt], [TT]plus-min.xvg[tt]",
        "and [TT]plus-plus.xvg[tt], or files for every individual ion pair if the [TT]-sep[tt]",
        "option is selected. In this case, files are named as ",
        "[TT]sb-(Resname)(Resnr)-(Atomnr)[tt].",
        "There may be [BB]many[bb] such files."
    };
    static gmx_bool bSep     = FALSE;
    static real     truncate = 1000.0;
    t_pargs         pa[]     = { { "-t",
                       FALSE,
                       etREAL,
                       { &truncate },
                       "Groups that are never closer than this distance are not plotted" },
                     { "-sep",
                       FALSE,
                       etBOOL,
                       { &bSep },
                       "Use separate files for each interaction (may be MANY)" } };
    t_filenm        fnm[]    = {
        { efTRX, "-f", nullptr, ffREAD },
        { efTPR, nullptr, nullptr, ffREAD },
    };
#define NFILE asize(fnm)

    FILE *             out[3], *fp;
    static const char* title[3] = { "Distance between positively charged groups",
                                    "Distance between negatively charged groups",
                                    "Distance between oppositely charged groups" };
    static const char* fn[3]    = { "plus-plus.xvg", "min-min.xvg", "plus-min.xvg" };
    int                nset[3]  = { 0, 0, 0 };

    t_topology*  top;
    PbcType      pbcType;
    t_trxstatus* status;
    int          i, j, k, m, nnn, teller, ncg;
    real         t, *time, qi, qj;
    t_charge*    cg;
    real***      cgdist;
    int**        nWithin;

    t_pbc             pbc;
    rvec*             x;
    matrix            box;
    gmx_output_env_t* oenv;

    if (!parse_common_args(
                &argc, argv, PCA_CAN_TIME, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &pbcType);
    cg  = mk_charge(&top->atoms, &ncg);
    snew(cgdist, ncg);
    snew(nWithin, ncg);
    for (i = 0; (i < ncg); i++)
    {
        snew(cgdist[i], ncg);
        snew(nWithin[i], ncg);
    }

    read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    teller = 0;
    time   = nullptr;
    do
    {
        srenew(time, teller + 1);
        time[teller] = t;

        set_pbc(&pbc, pbcType, box);

        for (i = 0; (i < ncg); i++)
        {
            for (j = i + 1; (j < ncg); j++)
            {
                srenew(cgdist[i][j], teller + 1);
                rvec dx;
                pbc_dx(&pbc, x[cg[i].cg], x[cg[j].cg], dx);
                cgdist[i][j][teller] = norm(dx);
                if (cgdist[i][j][teller] < truncate)
                {
                    nWithin[i][j] = 1;
                }
            }
        }

        teller++;
    } while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "\n");
    close_trx(status);

    if (bSep)
    {
        for (i = 0; (i < ncg); i++)
        {
            for (j = i + 1; (j < ncg); j++)
            {
                if (nWithin[i][j])
                {
                    std::string buf = gmx::formatString("sb-%s:%s.xvg", cg[i].label, cg[j].label);
                    fp = xvgropen(buf.c_str(), buf.c_str(), "Time (ps)", "Distance (nm)", oenv);
                    for (k = 0; (k < teller); k++)
                    {
                        fprintf(fp, "%10g  %10g\n", time[k], cgdist[i][j][k]);
                    }
                    xvgrclose(fp);
                }
            }
        }
    }
    else
    {

        for (m = 0; (m < 3); m++)
        {
            out[m] = xvgropen(fn[m], title[m], "Time (ps)", "Distance (nm)", oenv);
        }

        for (i = 0; (i < ncg); i++)
        {
            qi = cg[i].q;
            for (j = i + 1; (j < ncg); j++)
            {
                qj = cg[j].q;
                if (nWithin[i][j])
                {
                    auto buf = gmx::formatString("%s:%s", cg[i].label, cg[j].label);
                    if (qi * qj < 0)
                    {
                        nnn = 2;
                    }
                    else if (qi + qj > 0)
                    {
                        nnn = 0;
                    }
                    else
                    {
                        nnn = 1;
                    }

                    if (nset[nnn] == 0)
                    {
                        xvgrLegend(out[nnn], gmx::arrayRefFromArray(&buf, 1), oenv);
                    }
                    else
                    {
                        if (output_env_get_xvg_format(oenv) == XvgFormat::Xmgr)
                        {
                            fprintf(out[nnn], "@ legend string %d \"%s\"\n", nset[nnn], buf.c_str());
                        }
                        else if (output_env_get_xvg_format(oenv) == XvgFormat::Xmgrace)
                        {
                            fprintf(out[nnn], "@ s%d legend \"%s\"\n", nset[nnn], buf.c_str());
                        }
                    }
                    nset[nnn]++;
                    nWithin[i][j] = nnn + 1;
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
                for (j = i + 1; (j < ncg); j++)
                {
                    nnn = nWithin[i][j];
                    if (nnn > 0)
                    {
                        fprintf(out[nnn - 1], "  %10g", cgdist[i][j][k]);
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
