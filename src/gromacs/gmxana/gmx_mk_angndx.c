/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static int calc_ntype(int nft, int *ft, t_idef *idef)
{
    int  i, f, nf = 0;

    for (i = 0; (i < idef->ntypes); i++)
    {
        for (f = 0; f < nft; f++)
        {
            if (idef->functype[i] == ft[f])
            {
                nf++;
            }
        }
    }

    return nf;
}

static void fill_ft_ind(int nft, int *ft, t_idef *idef,
                        int ft_ind[], char *grpnames[])
{
    char buf[125];
    int  i, f, ftype, ind = 0;

    /* Loop over all the function types in the topology */
    for (i = 0; (i < idef->ntypes); i++)
    {
        ft_ind[i] = -1;
        /* Check all the selected function types */
        for (f = 0; f < nft; f++)
        {
            ftype = ft[f];
            if (idef->functype[i] == ftype)
            {
                ft_ind[i] = ind;
                switch (ftype)
                {
                    case F_ANGLES:
                        sprintf(buf, "Theta=%.1f_%.2f", idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case F_G96ANGLES:
                        sprintf(buf, "Cos_th=%.1f_%.2f", idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case F_UREY_BRADLEY:
                        sprintf(buf, "UB_th=%.1f_%.2f2f", idef->iparams[i].u_b.thetaA,
                                idef->iparams[i].u_b.kthetaA);
                        break;
                    case F_QUARTIC_ANGLES:
                        sprintf(buf, "Q_th=%.1f_%.2f_%.2f", idef->iparams[i].qangle.theta,
                                idef->iparams[i].qangle.c[0], idef->iparams[i].qangle.c[1]);
                        break;
                    case F_TABANGLES:
                        sprintf(buf, "Table=%d_%.2f", idef->iparams[i].tab.table,
                                idef->iparams[i].tab.kA);
                        break;
                    case F_PDIHS:
                        sprintf(buf, "Phi=%.1f_%d_%.2f", idef->iparams[i].pdihs.phiA,
                                idef->iparams[i].pdihs.mult, idef->iparams[i].pdihs.cpA);
                        break;
                    case F_IDIHS:
                        sprintf(buf, "Xi=%.1f_%.2f", idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case F_RBDIHS:
                        sprintf(buf, "RB-A1=%.2f", idef->iparams[i].rbdihs.rbcA[1]);
                        break;
                    case  F_RESTRANGLES:
                        sprintf(buf, "Theta=%.1f_%.2f", idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case  F_RESTRDIHS:
                        sprintf(buf, "Theta=%.1f_%.2f", idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case  F_CBTDIHS:
                        sprintf(buf, "CBT-A1=%.2f", idef->iparams[i].cbtdihs.cbtcA[1]);
                        break;

                    default:
                        gmx_fatal(FARGS, "Unsupported function type '%s' selected",
                                  interaction_function[ftype].longname);
                }
                grpnames[ind] = gmx_strdup(buf);
                ind++;
            }
        }
    }
}

static void fill_ang(int nft, int *ft, int fac,
                     int nr[], int *index[], int ft_ind[], t_topology *top,
                     gmx_bool bNoH, real hq)
{
    int         f, ftype, i, j, indg, nr_fac;
    gmx_bool    bUse;
    t_idef     *idef;
    t_atom     *atom;
    t_iatom    *ia;


    idef = &top->idef;
    atom = top->atoms.atom;

    for (f = 0; f < nft; f++)
    {
        ftype = ft[f];
        ia    = idef->il[ftype].iatoms;
        for (i = 0; (i < idef->il[ftype].nr); )
        {
            indg = ft_ind[ia[0]];
            if (indg == -1)
            {
                gmx_incons("Routine fill_ang");
            }
            bUse = TRUE;
            if (bNoH)
            {
                for (j = 0; j < fac; j++)
                {
                    if (atom[ia[1+j]].m < 1.5)
                    {
                        bUse = FALSE;
                    }
                }
            }
            if (hq)
            {
                for (j = 0; j < fac; j++)
                {
                    if (atom[ia[1+j]].m < 1.5 && fabs(atom[ia[1+j]].q) < hq)
                    {
                        bUse = FALSE;
                    }
                }
            }
            if (bUse)
            {
                if (nr[indg] % 1000 == 0)
                {
                    srenew(index[indg], fac*(nr[indg]+1000));
                }
                nr_fac = fac*nr[indg];
                for (j = 0; (j < fac); j++)
                {
                    index[indg][nr_fac+j] = ia[j+1];
                }
                nr[indg]++;
            }
            ia += interaction_function[ftype].nratoms+1;
            i  += interaction_function[ftype].nratoms+1;
        }
    }
}

static int *select_ftype(const char *opt, int *nft, int *mult)
{
    int *ft = NULL, ftype;

    if (opt[0] == 'a')
    {
        *mult = 3;
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if ((interaction_function[ftype].flags & IF_ATYPE) ||
                ftype == F_TABANGLES)
            {
                (*nft)++;
                srenew(ft, *nft);
                ft[*nft-1] = ftype;
            }
        }
    }
    else
    {
        *mult = 4;
        *nft  = 1;
        snew(ft, *nft);
        switch (opt[0])
        {
            case 'd':
                ft[0] = F_PDIHS;
                break;
            case 'i':
                ft[0] = F_IDIHS;
                break;
            case 'r':
                ft[0] = F_RBDIHS;
                break;
            default:
                break;
        }
    }

    return ft;
}

int gmx_mk_angndx(int argc, char *argv[])
{
    static const char *desc[] = {
        "[THISMODULE] makes an index file for calculation of",
        "angle distributions etc. It uses a run input file ([REF].tpx[ref]) for the",
        "definitions of the angles, dihedrals etc."
    };
    static const char *opt[] = { NULL, "angle", "dihedral", "improper", "ryckaert-bellemans", NULL };
    static gmx_bool    bH    = TRUE;
    static real        hq    = -1;
    t_pargs            pa[]  = {
        { "-type", FALSE, etENUM, {opt},
          "Type of angle" },
        { "-hyd", FALSE, etBOOL, {&bH},
          "Include angles with atoms with mass < 1.5" },
        { "-hq", FALSE, etREAL, {&hq},
          "Ignore angles with atoms with mass < 1.5 and magnitude of their charge less than this value" }
    };

    output_env_t       oenv;
    FILE              *out;
    t_topology        *top;
    int                i, j, ntype;
    int                nft = 0, *ft, mult = 0;
    int              **index;
    int               *ft_ind;
    int               *nr;
    char             **grpnames;
    t_filenm           fnm[] = {
        { efTPR, NULL, NULL, ffREAD  },
        { efNDX, NULL, "angle", ffWRITE }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }


    ft = select_ftype(opt[0], &nft, &mult);

    top = read_top(ftp2fn(efTPR, NFILE, fnm), NULL);

    ntype = calc_ntype(nft, ft, &(top->idef));
    snew(grpnames, ntype);
    snew(ft_ind, top->idef.ntypes);
    fill_ft_ind(nft, ft, &top->idef, ft_ind, grpnames);

    snew(nr, ntype);
    snew(index, ntype);
    fill_ang(nft, ft, mult, nr, index, ft_ind, top, !bH, hq);

    out = ftp2FILE(efNDX, NFILE, fnm, "w");
    for (i = 0; (i < ntype); i++)
    {
        if (nr[i] > 0)
        {
            fprintf(out, "[ %s ]\n", grpnames[i]);
            for (j = 0; (j < nr[i]*mult); j++)
            {
                fprintf(out, " %5d", index[i][j]+1);
                if ((j % 12) == 11)
                {
                    fprintf(out, "\n");
                }
            }
            fprintf(out, "\n");
        }
    }
    gmx_ffclose(out);

    return 0;
}
