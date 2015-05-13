/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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

#include <assert.h>
#include <math.h>
#include <string.h>

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/fatalerror.h"
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
            gmx_fatal(FARGS, "Increase MAXSLEN in the grompp code to at least %d,"
                      " or shorten your definition of bonds like %s to at most %d",
                      strlen(s)+1, s, MAXSLEN-1);
        }
    }
    else
    {
        strcpy(p->s, "");
    }
}

void pr_alloc (int extra, t_params *pr)
{
    int i, j;

    /* get new space for arrays */
    if (extra < 0)
    {
        gmx_fatal(FARGS, "Trying to make array smaller.\n");
    }
    if (extra == 0)
    {
        return;
    }
    assert(!((pr->nr == 0) && (pr->param != NULL)));
    if (pr->nr+extra > pr->maxnr)
    {
        pr->maxnr = max(1.2*pr->maxnr, pr->maxnr + extra);
        srenew(pr->param, pr->maxnr);
        for (i = pr->nr; (i < pr->maxnr); i++)
        {
            for (j = 0; (j < MAXATOMLIST); j++)
            {
                pr->param[i].a[j] = 0;
            }
            for (j = 0; (j < MAXFORCEPARAM); j++)
            {
                pr->param[i].c[j] = 0;
            }
            set_p_string(&(pr->param[i]), "");
        }
    }
}

void init_plist(t_params plist[])
{
    int i;

    for (i = 0; (i < F_NRE); i++)
    {
        plist[i].nr    = 0;
        plist[i].maxnr = 0;
        plist[i].param = NULL;

        /* CMAP */
        plist[i].ncmap        = 0;
        plist[i].cmap         = NULL;
        plist[i].grid_spacing = 0;
        plist[i].nc           = 0;
        plist[i].nct          = 0;
        plist[i].cmap_types   = NULL;
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

void add_param_to_list(t_params *list, t_param *b)
{
    int j;

    /* allocate one position extra */
    pr_alloc (1, list);

    /* fill the arrays */
    for (j = 0; (j < MAXFORCEPARAM); j++)
    {
        list->param[list->nr].c[j]   = b->c[j];
    }
    for (j = 0; (j < MAXATOMLIST); j++)
    {
        list->param[list->nr].a[j]   = b->a[j];
    }
    memset(list->param[list->nr].s, 0, sizeof(list->param[list->nr].s));

    list->nr++;
}


void init_molinfo(t_molinfo *mol)
{
    mol->nrexcl     = 0;
    mol->excl_set   = FALSE;
    mol->bProcessed = FALSE;
    init_plist(mol->plist);
    init_block(&mol->cgs);
    init_block(&mol->mols);
    init_blocka(&mol->excls);
    init_atom(&mol->atoms);
}

/* FREEING MEMORY */

void done_bt (t_params *pl)
{
    sfree(pl->param);
}

void done_mi(t_molinfo *mi)
{
    int i;

    done_atom (&(mi->atoms));
    done_block(&(mi->cgs));
    done_block(&(mi->mols));
    for (i = 0; (i < F_NRE); i++)
    {
        done_bt(&(mi->plist[i]));
    }
}

/* PRINTING STRUCTURES */

void print_bt(FILE *out, directive d, gpp_atomtype_t at,
              int ftype, int fsubtype, t_params plist[],
              gmx_bool bFullDih)
{
    /* This dihp is a DIRTY patch because the dih-types do not use
     * all four atoms to determine the type.
     */
    const int    dihp[2][2] = { { 1, 2 }, { 0, 3 } };
    t_params    *bt;
    int          i, j, f, nral, nrfp;
    gmx_bool     bDih = FALSE, bSwapParity;

    bt = &(plist[ftype]);

    if (!bt->nr)
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
        for (j = 2; (j < nral); j++)
        {
            fprintf (out, "  %3c%c", 'a', 'i'+j);
        }
    }
    else
    {
        for (j = 0; (j < 2); j++)
        {
            fprintf (out, "%3c%c", 'a', 'i'+dihp[f][j]);
        }
    }

    fprintf (out, " funct");
    for (j = 0; (j < nrfp); j++)
    {
        fprintf (out, " %12c%1d", 'c', j);
    }
    fprintf (out, "\n");

    /* print bondtypes */
    for (i = 0; (i < bt->nr); i++)
    {
        bSwapParity = (bt->param[i].C0 == NOTSET) && (bt->param[i].C1 == -1);
        if (!bDih)
        {
            for (j = 0; (j < nral); j++)
            {
                fprintf (out, "%5s ", get_atomtype_name(bt->param[i].a[j], at));
            }
        }
        else
        {
            for (j = 0; (j < 2); j++)
            {
                fprintf (out, "%5s ", get_atomtype_name(bt->param[i].a[dihp[f][j]], at));
            }
        }
        fprintf (out, "%5d ", bSwapParity ? -f-1 : f+1);

        if (bt->param[i].s[0])
        {
            fprintf(out, "   %s", bt->param[i].s);
        }
        else
        {
            for (j = 0; (j < nrfp && (bt->param[i].c[j] != NOTSET)); j++)
            {
                fprintf (out, "%13.6e ", bt->param[i].c[j]);
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
            for (j = block->index[i]; (j < ((int)block->index[i+1])); j++)
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
    atom_id     i;
    gmx_bool    have_excl;
    int         j;

    have_excl = FALSE;
    for (i = 0; i < natoms && !have_excl; i++)
    {
        have_excl = (excls[i].nr > 0);
    }

    if (have_excl)
    {
        fprintf (out, "[ %s ]\n", dir2str(d_exclusions));
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

void print_atoms(FILE *out, gpp_atomtype_t atype, t_atoms *at, int *cgnr,
                 gmx_bool bRTPresname)
{
    int         i, ri;
    int         tpA, tpB;
    const char *as;
    char       *tpnmA, *tpnmB;
    double      qres, qtot;

    as = dir2str(d_atoms);
    fprintf(out, "[ %s ]\n", as);
    fprintf(out, "; %4s %10s %6s %7s%6s %6s %10s %10s %6s %10s %10s\n",
            "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass", "typeB", "chargeB", "massB");

    qtot  = 0;

    if (debug)
    {
        fprintf(debug, "This molecule has %d atoms and %d residues\n",
                at->nr, at->nres);
    }

    if (at->nres)
    {
        /* if the information is present... */
        for (i = 0; (i < at->nr); i++)
        {
            ri = at->atom[i].resind;
            if ((i == 0 || ri != at->atom[i-1].resind) &&
                at->resinfo[ri].rtp != NULL)
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
            if ((tpnmA = get_atomtype_name(tpA, atype)) == NULL)
            {
                gmx_fatal(FARGS, "tpA = %d, i= %d in print_atoms", tpA, i);
            }

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
                if ((tpnmB = get_atomtype_name(tpB, atype)) == NULL)
                {
                    gmx_fatal(FARGS, "tpB = %d, i= %d in print_atoms", tpB, i);
                }
                fprintf(out, " %6s %10g %10g",
                        tpnmB, at->atom[i].qB, at->atom[i].mB);
            }
            qtot += (double)at->atom[i].q;
            if (fabs(qtot) < 4*GMX_REAL_EPS)
            {
                qtot = 0;
            }
            fprintf(out, "   ; qtot %.4g\n", qtot);
        }
    }
    fprintf(out, "\n");
    fflush(out);
}

void print_bondeds(FILE *out, int natoms, directive d,
                   int ftype, int fsubtype, t_params plist[])
{
    t_symtab       stab;
    gpp_atomtype_t atype;
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
        add_atomtype(atype, &stab, a, buf, param, 0, 0, 0, 0, 0, 0, 0);
    }
    print_bt(out, d, atype, ftype, fsubtype, plist, TRUE);

    done_symtab(&stab);
    sfree(a);
    sfree(param);
    done_atomtype(atype);
}
