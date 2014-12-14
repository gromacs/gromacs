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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "tomorse.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    char *ai, *aj;
    real  e_diss;
} t_2morse;

static t_2morse *read_dissociation_energies(int *n2morse)
{
    FILE       *fp;
    char        ai[32], aj[32];
    double      e_diss;
    const char *fn     = "edissoc.dat";
    t_2morse   *t2m    = NULL;
    int         maxn2m = 0, n2m = 0;
    int         nread;

    /* Open the file with dissociation energies */
    fp = libopen(fn);
    do
    {
        /* Try and read two atom names and an energy term from it */
        nread = fscanf(fp, "%s%s%lf", ai, aj, &e_diss);
        if (nread == 3)
        {
            /* If we got three terms, it probably was OK, no further checking */
            if (n2m >= maxn2m)
            {
                /* Increase memory for 16 at once, some mallocs are stupid */
                maxn2m += 16;
                srenew(t2m, maxn2m);
            }
            /* Copy the values */
            t2m[n2m].ai     = gmx_strdup(ai);
            t2m[n2m].aj     = gmx_strdup(aj);
            t2m[n2m].e_diss = e_diss;
            /* Increment counter */
            n2m++;
        }
        /* If we did not read three items, quit reading */
    }
    while (nread == 3);
    gmx_ffclose(fp);

    /* Set the return values */
    *n2morse = n2m;

    return t2m;
}

static int nequal(char *a1, char *a2)
{
    int i;

    /* Count the number of (case insensitive) characters that are equal in
     * two strings. If they are equally long their respective null characters are
     * counted also.
     */
    for (i = 0; (a1[i] != '\0') && (a2[i] != '\0'); i++)
    {
        if (toupper(a1[i]) != toupper(a2[i]))
        {
            break;
        }
    }
    if ((a1[i] == '\0') && (a2[i] == '\0'))
    {
        i++;
    }

    return i;
}

static real search_e_diss(int n2m, t_2morse t2m[], char *ai, char *aj)
{
    int  i;
    int  ibest = -1;
    int  nii, njj, nbstii = 0, nbstjj = 0;
    real ediss = 400;

    /* Do a best match search for dissociation energies */
    for (i = 0; (i < n2m); i++)
    {
        /* Check for a perfect match */
        if (((gmx_strcasecmp(t2m[i].ai, ai) == 0) && (gmx_strcasecmp(t2m[i].aj, aj) == 0)) ||
            ((gmx_strcasecmp(t2m[i].aj, ai) == 0) && (gmx_strcasecmp(t2m[i].ai, aj) == 0)))
        {
            ibest = i;
            break;
        }
        else
        {
            /* Otherwise count the number of equal characters in the strings ai and aj
             * and the ones from the file
             */
            nii = nequal(t2m[i].ai, ai);
            njj = nequal(t2m[i].aj, aj);
            if (((nii >  nbstii) && (njj >= nbstjj)) ||
                ((nii >= nbstii) && (njj >  nbstjj)))
            {
                if ((nii > 0) && (njj > 0))
                {
                    ibest  = i;
                    nbstii = nii;
                    nbstjj = njj;
                }
            }
            else
            {
                /* Swap ai and aj (at least in counting the number of equal chars) */
                nii = nequal(t2m[i].ai, aj);
                njj = nequal(t2m[i].aj, ai);
                if (((nii >  nbstii) && (njj >= nbstjj)) ||
                    ((nii >= nbstii) && (njj >  nbstjj)))
                {
                    if ((nii > 0) && (njj > 0))
                    {
                        ibest  = i;
                        nbstii = nii;
                        nbstjj = njj;
                    }
                }
            }
        }
    }
    /* Return the dissocation energy corresponding to the best match, if we have
     * found one. Do some debug output anyway.
     */
    if (ibest == -1)
    {
        if (debug)
        {
            fprintf(debug, "MORSE: Couldn't find E_diss for bond %s - %s, using default %g\n", ai, aj, ediss);
        }
        return ediss;
    }
    else
    {
        if (debug)
        {
            fprintf(debug, "MORSE: Dissoc. E (%10.3f) for bond %4s-%4s taken from bond %4s-%4s\n",
                    t2m[ibest].e_diss, ai, aj, t2m[ibest].ai, t2m[ibest].aj);
        }
        return t2m[ibest].e_diss;
    }
}

void convert_harmonics(int nrmols, t_molinfo mols[], gpp_atomtype_t atype)
{
    int       n2m;
    t_2morse *t2m;

    int       i, j, k, last, ni, nj;
    int       nrharm, nrmorse, bb;
    real      edis, kb, b0, beta;
    gmx_bool *bRemoveHarm;

    /* First get the data */
    t2m = read_dissociation_energies(&n2m);
    if (debug)
    {
        fprintf(debug, "MORSE: read %d dissoc energies\n", n2m);
    }
    if (n2m <= 0)
    {
        fprintf(stderr, "No dissocation energies read\n");
        return;
    }

    /* For all the molecule types */
    for (i = 0; (i < nrmols); i++)
    {
        /* Check how many morse and harmonic BONDSs there are, increase size of
         * morse with the number of harmonics
         */
        nrmorse = mols[i].plist[F_MORSE].nr;

        for (bb = 0; (bb < F_NRE); bb++)
        {
            if ((interaction_function[bb].flags & IF_BTYPE) && (bb != F_MORSE))
            {
                nrharm  = mols[i].plist[bb].nr;
                pr_alloc(nrharm, &(mols[i].plist[F_MORSE]));
                snew(bRemoveHarm, nrharm);

                /* Now loop over the harmonics, trying to convert them */
                for (j = 0; (j < nrharm); j++)
                {
                    ni   = mols[i].plist[bb].param[j].AI;
                    nj   = mols[i].plist[bb].param[j].AJ;
                    edis =
                        search_e_diss(n2m, t2m,
                                      get_atomtype_name(mols[i].atoms.atom[ni].type, atype),
                                      get_atomtype_name(mols[i].atoms.atom[nj].type, atype));
                    if (edis != 0)
                    {
                        bRemoveHarm[j] = TRUE;
                        b0             = mols[i].plist[bb].param[j].c[0];
                        kb             = mols[i].plist[bb].param[j].c[1];
                        beta           = sqrt(kb/(2*edis));
                        mols[i].plist[F_MORSE].param[nrmorse].a[0] = ni;
                        mols[i].plist[F_MORSE].param[nrmorse].a[1] = nj;
                        mols[i].plist[F_MORSE].param[nrmorse].c[0] = b0;
                        mols[i].plist[F_MORSE].param[nrmorse].c[1] = edis;
                        mols[i].plist[F_MORSE].param[nrmorse].c[2] = beta;
                        nrmorse++;
                    }
                }
                mols[i].plist[F_MORSE].nr = nrmorse;

                /* Now remove the harmonics */
                for (j = last = 0; (j < nrharm); j++)
                {
                    if (!bRemoveHarm[j])
                    {
                        /* Copy it to the last position */
                        for (k = 0; (k < MAXATOMLIST); k++)
                        {
                            mols[i].plist[bb].param[last].a[k] =
                                mols[i].plist[bb].param[j].a[k];
                        }
                        for (k = 0; (k < MAXFORCEPARAM); k++)
                        {
                            mols[i].plist[bb].param[last].c[k] =
                                mols[i].plist[bb].param[j].c[k];
                        }
                        last++;
                    }
                }
                sfree(bRemoveHarm);
                fprintf(stderr, "Converted %d out of %d %s to morse bonds for mol %d\n",
                        nrharm-last, nrharm, interaction_function[bb].name, i);
                mols[i].plist[bb].nr = last;
            }
        }
    }
    sfree(t2m);
}
