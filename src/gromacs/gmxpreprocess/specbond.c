/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
#include "gmxpre.h"

#include "specbond.h"

#include <ctype.h>
#include <math.h>
#include <string.h>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

gmx_bool yesno(void)
{
    char c;

    do
    {
        c = toupper(fgetc(stdin));
    }
    while ((c != 'Y') && (c != 'N'));

    return (c == 'Y');
}

t_specbond *get_specbonds(int *nspecbond)
{
    const char  *sbfile = "specbond.dat";

    t_specbond  *sb = NULL;
    char         r1buf[32], r2buf[32], a1buf[32], a2buf[32], nr1buf[32], nr2buf[32];
    double       length;
    int          nb1, nb2;
    char       **lines;
    int          nlines, i, n;

    nlines = get_lines(sbfile, &lines);
    if (nlines > 0)
    {
        snew(sb, nlines);
    }

    n = 0;
    for (i = 0; (i < nlines); i++)
    {
        if (sscanf(lines[i], "%s%s%d%s%s%d%lf%s%s",
                   r1buf, a1buf, &nb1, r2buf, a2buf, &nb2, &length, nr1buf, nr2buf) != 9)
        {
            fprintf(stderr, "Invalid line '%s' in %s\n", lines[i], sbfile);
        }
        else
        {
            sb[n].res1    = gmx_strdup(r1buf);
            sb[n].res2    = gmx_strdup(r2buf);
            sb[n].newres1 = gmx_strdup(nr1buf);
            sb[n].newres2 = gmx_strdup(nr2buf);
            sb[n].atom1   = gmx_strdup(a1buf);
            sb[n].atom2   = gmx_strdup(a2buf);
            sb[n].nbond1  = nb1;
            sb[n].nbond2  = nb2;
            sb[n].length  = length;
            n++;
        }
        sfree(lines[i]);
    }
    if (nlines > 0)
    {
        sfree(lines);
    }
    fprintf(stderr, "%d out of %d lines of %s converted successfully\n",
            n, nlines, sbfile);

    *nspecbond = n;

    return sb;
}

void done_specbonds(int nsb, t_specbond sb[])
{
    int i;

    for (i = 0; (i < nsb); i++)
    {
        sfree(sb[i].res1);
        sfree(sb[i].res2);
        sfree(sb[i].atom1);
        sfree(sb[i].atom2);
        sfree(sb[i].newres1);
        sfree(sb[i].newres2);
    }
}

static gmx_bool is_special(int nsb, t_specbond sb[], char *res, char *atom)
{
    int i;

    for (i = 0; (i < nsb); i++)
    {
        if (((strncmp(sb[i].res1, res, 3) == 0) &&
             (gmx_strcasecmp(sb[i].atom1, atom) == 0)) ||
            ((strncmp(sb[i].res2, res, 3) == 0) &&
             (gmx_strcasecmp(sb[i].atom2, atom) == 0)))
        {
            return TRUE;
        }
    }
    return FALSE;
}

static gmx_bool is_bond(int nsb, t_specbond sb[], t_atoms *pdba, int a1, int a2,
                        real d, int *index_sb, gmx_bool *bSwap)
{
    int   i;
    char *at1, *at2, *res1, *res2;

    at1  = *pdba->atomname[a1];
    at2  = *pdba->atomname[a2];
    res1 = *pdba->resinfo[pdba->atom[a1].resind].name;
    res2 = *pdba->resinfo[pdba->atom[a2].resind].name;

    if (debug)
    {
        fprintf(stderr, "Checking %s-%d %s-%d and %s-%d %s-%d: %g ",
                res1, pdba->resinfo[pdba->atom[a1].resind].nr, at1, a1+1,
                res2, pdba->resinfo[pdba->atom[a2].resind].nr, at2, a2+1, d);
    }

    for (i = 0; (i < nsb); i++)
    {
        *index_sb = i;
        if (((strncmp(sb[i].res1, res1, 3) == 0)  &&
             (gmx_strcasecmp(sb[i].atom1, at1) == 0) &&
             (strncmp(sb[i].res2, res2, 3) == 0)  &&
             (gmx_strcasecmp(sb[i].atom2, at2) == 0)))
        {
            *bSwap = FALSE;
            if ((0.9*sb[i].length < d) && (1.1*sb[i].length > d))
            {
                if (debug)
                {
                    fprintf(stderr, "%g\n", sb[i].length);
                }
                return TRUE;
            }
        }
        if (((strncmp(sb[i].res1, res2, 3) == 0)  &&
             (gmx_strcasecmp(sb[i].atom1, at2) == 0) &&
             (strncmp(sb[i].res2, res1, 3) == 0)  &&
             (gmx_strcasecmp(sb[i].atom2, at1) == 0)))
        {
            *bSwap = TRUE;
            if ((0.9*sb[i].length < d) && (1.1*sb[i].length > d))
            {
                if (debug)
                {
                    fprintf(stderr, "%g\n", sb[i].length);
                }
                return TRUE;
            }
        }
    }
    if (debug)
    {
        fprintf(stderr, "\n");
    }
    return FALSE;
}

static void rename_1res(t_atoms *pdba, int resind, char *newres, gmx_bool bVerbose)
{
    if (bVerbose)
    {
        printf("Using rtp entry %s for %s %d\n",
               newres,
               *pdba->resinfo[resind].name,
               pdba->resinfo[resind].nr);
    }
    /* this used to free *resname, which messes up the symtab! */
    snew(pdba->resinfo[resind].rtp, 1);
    *pdba->resinfo[resind].rtp = gmx_strdup(newres);
}

int mk_specbonds(t_atoms *pdba, rvec x[], gmx_bool bInteractive,
                 t_ssbond **specbonds, gmx_bool bVerbose)
{
    t_specbond *sb    = NULL;
    t_ssbond   *bonds = NULL;
    int         nsb;
    int         nspec, nbonds;
    int        *specp, *sgp;
    gmx_bool    bDoit, bSwap;
    int         i, j, b, e, e2;
    int         ai, aj, index_sb;
    real      **d;
    char        buf[10];

    nbonds = 0;
    sb     = get_specbonds(&nsb);

    if (nsb > 0)
    {
        snew(specp, pdba->nr);
        snew(sgp, pdba->nr);

        nspec = 0;
        for (i = 0; (i < pdba->nr); i++)
        {
            /* Check if this atom is special and if it is not a double atom
             * in the input that still needs to be removed.
             */
            if (is_special(nsb, sb, *pdba->resinfo[pdba->atom[i].resind].name,
                           *pdba->atomname[i]) &&
                !(nspec > 0 &&
                  pdba->atom[sgp[nspec-1]].resind == pdba->atom[i].resind &&
                  gmx_strcasecmp(*pdba->atomname[sgp[nspec-1]],
                                 *pdba->atomname[i]) == 0))
            {
                specp[nspec] = pdba->atom[i].resind;
                sgp[nspec]   = i;
                nspec++;
            }
        }
        /* distance matrix d[nspec][nspec] */
        snew(d, nspec);
        for (i = 0; (i < nspec); i++)
        {
            snew(d[i], nspec);
        }

        for (i = 0; (i < nspec); i++)
        {
            ai = sgp[i];
            for (j = 0; (j < nspec); j++)
            {
                aj      = sgp[j];
                d[i][j] = sqrt(distance2(x[ai], x[aj]));
            }
        }
        if (nspec > 1)
        {
#define MAXCOL 7
            fprintf(stderr, "Special Atom Distance matrix:\n");
            for (b = 0; (b < nspec); b += MAXCOL)
            {
                /* print resname/number column headings */
                fprintf(stderr, "%8s%8s", "", "");
                e = min(b+MAXCOL, nspec-1);
                for (i = b; (i < e); i++)
                {
                    sprintf(buf, "%s%d", *pdba->resinfo[pdba->atom[sgp[i]].resind].name,
                            pdba->resinfo[specp[i]].nr);
                    fprintf(stderr, "%8s", buf);
                }
                fprintf(stderr, "\n");
                /* print atomname/number column headings */
                fprintf(stderr, "%8s%8s", "", "");
                e = min(b+MAXCOL, nspec-1);
                for (i = b; (i < e); i++)
                {
                    sprintf(buf, "%s%d", *pdba->atomname[sgp[i]], sgp[i]+1);
                    fprintf(stderr, "%8s", buf);
                }
                fprintf(stderr, "\n");
                /* print matrix */
                e = min(b+MAXCOL, nspec);
                for (i = b+1; (i < nspec); i++)
                {
                    sprintf(buf, "%s%d", *pdba->resinfo[pdba->atom[sgp[i]].resind].name,
                            pdba->resinfo[specp[i]].nr);
                    fprintf(stderr, "%8s", buf);
                    sprintf(buf, "%s%d", *pdba->atomname[sgp[i]], sgp[i]+1);
                    fprintf(stderr, "%8s", buf);
                    e2 = min(i, e);
                    for (j = b; (j < e2); j++)
                    {
                        fprintf(stderr, " %7.3f", d[i][j]);
                    }
                    fprintf(stderr, "\n");
                }
            }
        }

        snew(bonds, nspec);

        for (i = 0; (i < nspec); i++)
        {
            ai = sgp[i];
            for (j = i+1; (j < nspec); j++)
            {
                aj = sgp[j];
                if (is_bond(nsb, sb, pdba, ai, aj, d[i][j], &index_sb, &bSwap))
                {
                    fprintf(stderr, "%s %s-%d %s-%d and %s-%d %s-%d%s",
                            bInteractive ? "Link" : "Linking",
                            *pdba->resinfo[pdba->atom[ai].resind].name,
                            pdba->resinfo[specp[i]].nr,
                            *pdba->atomname[ai], ai+1,
                            *pdba->resinfo[pdba->atom[aj].resind].name,
                            pdba->resinfo[specp[j]].nr,
                            *pdba->atomname[aj], aj+1,
                            bInteractive ? " (y/n) ?" : "...\n");
                    bDoit = bInteractive ? yesno() : TRUE;

                    if (bDoit)
                    {
                        /* Store the residue numbers in the bonds array */
                        bonds[nbonds].res1 = specp[i];
                        bonds[nbonds].res2 = specp[j];
                        bonds[nbonds].a1   = gmx_strdup(*pdba->atomname[ai]);
                        bonds[nbonds].a2   = gmx_strdup(*pdba->atomname[aj]);
                        /* rename residues */
                        if (bSwap)
                        {
                            rename_1res(pdba, specp[i], sb[index_sb].newres2, bVerbose);
                            rename_1res(pdba, specp[j], sb[index_sb].newres1, bVerbose);
                        }
                        else
                        {
                            rename_1res(pdba, specp[i], sb[index_sb].newres1, bVerbose);
                            rename_1res(pdba, specp[j], sb[index_sb].newres2, bVerbose);
                        }
                        nbonds++;
                    }
                }
            }
        }

        for (i = 0; (i < nspec); i++)
        {
            sfree(d[i]);
        }
        sfree(d);
        sfree(sgp);
        sfree(specp);

        done_specbonds(nsb, sb);
        sfree(sb);
    }

    *specbonds = bonds;

    return nbonds;
}
