/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include "toputil.h"

#include <climits>
#include <cmath>
#include <cstring>

#include <algorithm>

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/* UTILITIES */

void set_p_string(t_param *p, const char *s)
{
    if (s)
    {
        if (strlen(s) < sizeof(p->s)-1)
        {
            strncpy(p->s, s, sizeof(p->s));
        }
        else
        {
            gmx_fatal(FARGS, "Increase MAXSLEN in the grompp code to at least %zu,"
                      " or shorten your definition of bonds like %s to at most %d",
                      strlen(s)+1, s, MAXSLEN-1);
        }
    }
    else
    {
        strcpy(p->s, "");
    }
}

void cp_param(t_param *dest, t_param *src)
{
    int j;

    for (j = 0; (j < MAXATOMLIST); j++)
    {
        dest->a[j] = src->a[j];
    }
    for (j = 0; (j < MAXFORCEPARAM); j++)
    {
        dest->c[j] = src->c[j];
    }
    strncpy(dest->s, src->s, sizeof(dest->s));
}

void add_param_to_list(SystemParameters *list, t_param *b)
{
    /* allocate one position extra */
    list->param.emplace_back();
    initializeTparam(&list->param.back());

    /* fill the arrays */
    for (int j = 0; (j < MAXFORCEPARAM); j++)
    {
        list->param.back().c[j]   = b->c[j];
    }
    for (int j = 0; (j < MAXATOMLIST); j++)
    {
        list->param.back().a[j]   = b->a[j];
    }
}

void initializeTparam(t_param *param)
{
    for (int j = 0; (j < MAXATOMLIST); j++)
    {
        param->a[j] = 0;
    }
    for (int j = 0; (j < MAXFORCEPARAM); j++)
    {
        param->c[j] = 0;
    }
    set_p_string(param, "");
}

/* FREEING MEMORY */

void done_mi(MoleculeInformation *mi)
{
    done_atom (&(mi->atoms));
    done_block(&(mi->cgs));
    done_block(&(mi->mols));
}

/* PRINTING STRUCTURES */

static void print_bt(FILE *out, Directive d, gpp_atomtype *at,
                     int ftype, int fsubtype, gmx::ArrayRef<const SystemParameters> plist,
                     bool bFullDih)
{
    /* This dihp is a DIRTY patch because the dih-types do not use
     * all four atoms to determine the type.
     */
    const int               dihp[2][2] = { { 1, 2 }, { 0, 3 } };
    int                     f, nral, nrfp;
    bool                    bDih = false, bSwapParity;

    const SystemParameters *bt = &(plist[ftype]);

    if (bt->nr() == 0)
    {
        return;
    }

    f = 0;
    switch (ftype)
    {
        case F_G96ANGLES:
            f = 1;
            break;
        case F_G96BONDS:
            f = 1;
            break;
        case F_MORSE:
            f = 2;
            break;
        case F_CUBICBONDS:
            f = 3;
            break;
        case F_CONNBONDS:
            f = 4;
            break;
        case F_HARMONIC:
            f = 5;
            break;
        case F_CROSS_BOND_ANGLES:
            f = 2;
            break;
        case F_CROSS_BOND_BONDS:
            f = 3;
            break;
        case F_UREY_BRADLEY:
            f = 4;
            break;
        case F_PDIHS:
        case F_RBDIHS:
        case F_FOURDIHS:
            bDih = TRUE;
            break;
        case F_IDIHS:
            f    = 1;
            bDih = TRUE;
            break;
        case F_CONSTRNC:
            f = 1;
            break;
        case F_VSITE3FD:
            f = 1;
            break;
        case F_VSITE3FAD:
            f = 2;
            break;
        case F_VSITE3OUT:
            f = 3;
            break;
        case F_VSITE4FDN:
            f = 1;
            break;
        case F_CMAP:
            f = 1;
            break;

        default:
            bDih = FALSE;
    }
    if (bFullDih)
    {
        bDih = FALSE;
    }
    if (fsubtype)
    {
        f = fsubtype-1;
    }

    nral = NRAL(ftype);
    nrfp = NRFP(ftype);

    /* header */
    fprintf(out, "[ %s ]\n", dir2str(d));
    fprintf(out, "; ");
    if (!bDih)
    {
        fprintf (out, "%3s  %4s", "ai", "aj");
        for (int j = 2; (j < nral); j++)
        {
            fprintf (out, "  %3c%c", 'a', 'i'+j);
        }
    }
    else
    {
        for (int j = 0; (j < 2); j++)
        {
            fprintf (out, "%3c%c", 'a', 'i'+dihp[f][j]);
        }
    }

    fprintf (out, " funct");
    for (int j = 0; (j < nrfp); j++)
    {
        fprintf (out, " %12c%1d", 'c', j);
    }
    fprintf (out, "\n");

    /* print bondtypes */
    for (const auto &par : bt->param)
    {
        bSwapParity = (par.c0() == NOTSET) && (par.c1() == -1);
        if (!bDih)
        {
            for (int j = 0; (j < nral); j++)
            {
                fprintf (out, "%5s ", get_atomtype_name(par.a[j], at));
            }
        }
        else
        {
            for (int j = 0; (j < 2); j++)
            {
                fprintf (out, "%5s ", get_atomtype_name(par.a[dihp[f][j]], at));
            }
        }
        fprintf (out, "%5d ", bSwapParity ? -f-1 : f+1);

        if (par.s[0])
        {
            fprintf(out, "   %s", par.s);
        }
        else
        {
            for (int j = 0; (j < nrfp && (par.c[j] != NOTSET)); j++)
            {
                fprintf (out, "%13.6e ", par.c[j]);
            }
        }

        fprintf (out, "\n");
    }
    fprintf (out, "\n");
    fflush (out);
}

void print_blocka(FILE *out, const char *szName,
                  const char *szIndex, const char *szA,
                  t_blocka *block)
{
    int i, j;

    fprintf (out, "; %s\n", szName);
    fprintf (out, "; %4s    %s\n", szIndex, szA);
    for (i = 0; (i < block->nr); i++)
    {
        for (i = 0; (i < block->nr); i++)
        {
            fprintf (out, "%6d", i+1);
            for (j = block->index[i]; (j < (block->index[i+1])); j++)
            {
                fprintf (out, "%5d", block->a[j]+1);
            }
            fprintf (out, "\n");
        }
        fprintf (out, "\n");
    }
}

void print_excl(FILE *out, int natoms, t_excls excls[])
{
    int         i;
    bool        have_excl;
    int         j;

    have_excl = FALSE;
    for (i = 0; i < natoms && !have_excl; i++)
    {
        have_excl = (excls[i].nr > 0);
    }

    if (have_excl)
    {
        fprintf (out, "[ %s ]\n", dir2str(Directive::d_exclusions));
        fprintf (out, "; %4s    %s\n", "i", "excluded from i");
        for (i = 0; i < natoms; i++)
        {
            if (excls[i].nr > 0)
            {
                fprintf (out, "%6d ", i+1);
                for (j = 0; j < excls[i].nr; j++)
                {
                    fprintf (out, " %5d", excls[i].e[j]+1);
                }
                fprintf (out, "\n");
            }
        }
        fprintf (out, "\n");
        fflush(out);
    }
}

static double get_residue_charge(const t_atoms *atoms, int at)
{
    int    ri;
    double q;

    ri = atoms->atom[at].resind;
    q  = 0;
    while (at < atoms->nr && atoms->atom[at].resind == ri)
    {
        q += atoms->atom[at].q;
        at++;
    }

    return q;
}

void print_atoms(FILE *out, gpp_atomtype *atype, t_atoms *at, int *cgnr,
                 bool bRTPresname)
{
    int         i, ri;
    int         tpA, tpB;
    const char *as;
    char       *tpnmA, *tpnmB;
    double      qres, qtot;

    as = dir2str(Directive::d_atoms);
    fprintf(out, "[ %s ]\n", as);
    fprintf(out, "; %4s %10s %6s %7s%6s %6s %10s %10s %6s %10s %10s\n",
            "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass", "typeB", "chargeB", "massB");

    qtot  = 0;

    if (at->nres)
    {
        /* if the information is present... */
        for (i = 0; (i < at->nr); i++)
        {
            ri = at->atom[i].resind;
            if ((i == 0 || ri != at->atom[i-1].resind) &&
                at->resinfo[ri].rtp != nullptr)
            {
                qres = get_residue_charge(at, i);
                fprintf(out, "; residue %3d %-3s rtp %-4s q ",
                        at->resinfo[ri].nr,
                        *at->resinfo[ri].name,
                        *at->resinfo[ri].rtp);
                if (fabs(qres) < 0.001)
                {
                    fprintf(out, " %s", "0.0");
                }
                else
                {
                    fprintf(out, "%+3.1f", qres);
                }
                fprintf(out, "\n");
            }
            tpA = at->atom[i].type;
            if ((tpnmA = get_atomtype_name(tpA, atype)) == nullptr)
            {
                gmx_fatal(FARGS, "tpA = %d, i= %d in print_atoms", tpA, i);
            }

            /* This is true by construction, but static analysers don't know */
            GMX_ASSERT(!bRTPresname || at->resinfo[at->atom[i].resind].rtp, "-rtpres did not have residue name available");
            fprintf(out, "%6d %10s %6d%c %5s %6s %6d %10g %10g",
                    i+1, tpnmA,
                    at->resinfo[ri].nr,
                    at->resinfo[ri].ic,
                    bRTPresname ?
                    *(at->resinfo[at->atom[i].resind].rtp) :
                    *(at->resinfo[at->atom[i].resind].name),
                    *(at->atomname[i]), cgnr[i],
                    at->atom[i].q, at->atom[i].m);
            if (PERTURBED(at->atom[i]))
            {
                tpB = at->atom[i].typeB;
                if ((tpnmB = get_atomtype_name(tpB, atype)) == nullptr)
                {
                    gmx_fatal(FARGS, "tpB = %d, i= %d in print_atoms", tpB, i);
                }
                fprintf(out, " %6s %10g %10g",
                        tpnmB, at->atom[i].qB, at->atom[i].mB);
            }
            // Accumulate the total charge to help troubleshoot issues.
            qtot += static_cast<double>(at->atom[i].q);
            // Round it to zero if it is close to zero, because
            // printing -9.34e-5 confuses users.
            if (fabs(qtot) < 0.0001)
            {
                qtot = 0;
            }
            // Write the total charge for the last atom of the system
            // and/or residue, because generally that's where it is
            // expected to be an integer.
            if (i == at->nr-1 || ri != at->atom[i+1].resind)
            {
                fprintf(out, "   ; qtot %.4g\n", qtot);
            }
            else
            {
                fputs("\n", out);
            }
        }
    }
    fprintf(out, "\n");
    fflush(out);
}

void print_bondeds(FILE *out, int natoms, Directive d,
                   int ftype, int fsubtype, gmx::ArrayRef<const SystemParameters> plist)
{
    t_symtab       stab;
    gpp_atomtype  *atype;
    t_param       *param;
    t_atom        *a;
    int            i;

    atype = init_atomtype();
    snew(a, 1);
    snew(param, 1);
    open_symtab(&stab);
    for (i = 0; (i < natoms); i++)
    {
        char buf[12];
        sprintf(buf, "%4d", (i+1));
        add_atomtype(atype, &stab, a, buf, param, 0, 0);
    }
    print_bt(out, d, atype, ftype, fsubtype, plist, TRUE);

    done_symtab(&stab);
    sfree(a);
    sfree(param);
    done_atomtype(atype);
}
