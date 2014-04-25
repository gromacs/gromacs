/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gromacs/energyanalysis/select.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/fatalerror.h"

int *select_it(int nre, char *nm[], int *nset)
{
    gmx_bool *bE;
    int       n, k, j, i;
    int      *set;
    gmx_bool  bVerbose = TRUE;

    if ((getenv("VERBOSE")) != NULL)
    {
        bVerbose = FALSE;
    }

    fprintf(stderr, "Select the terms you want from the following list\n");
    fprintf(stderr, "End your selection with 0\n");

    if (bVerbose)
    {
        for (k = 0; (k < nre); )
        {
            for (j = 0; (j < 4) && (k < nre); j++, k++)
            {
                fprintf(stderr, " %3d=%14s", k+1, nm[k]);
            }
            fprintf(stderr, "\n");
        }
    }

    snew(bE, nre);
    do
    {
        if (1 != scanf("%d", &n))
        {
            gmx_fatal(FARGS, "Error reading user input");
        }
        if ((n > 0) && (n <= nre))
        {
            bE[n-1] = TRUE;
        }
    }
    while (n != 0);

    snew(set, nre);
    for (i = (*nset) = 0; (i < nre); i++)
    {
        if (bE[i])
        {
            set[(*nset)++] = i;
        }
    }

    sfree(bE);

    return set;
}

static void chomp(char *buf)
{
    int len = strlen(buf);

    while ((len > 0) && (buf[len-1] == '\n'))
    {
        buf[len-1] = '\0';
        len--;
    }
}

int *select_by_name(int nre, gmx_enxnm_t *nm, int *nset)
{
    gmx_bool   *bE;
    int         n, k, kk, j, i, nmatch, nind, nss;
    int        *set;
    gmx_bool    bEOF, bVerbose = TRUE, bLong = FALSE;
    char       *ptr, buf[STRLEN];
    const char *fm4   = "%3d  %-14s";
    const char *fm2   = "%3d  %-34s";
    char      **newnm = NULL;

    if ((getenv("VERBOSE")) != NULL)
    {
        bVerbose = FALSE;
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Select the terms you want from the following list by\n");
    fprintf(stderr, "selecting either (part of) the name or the number or a combination.\n");
    fprintf(stderr, "End your selection with an empty line or a zero.\n");
    fprintf(stderr, "-------------------------------------------------------------------\n");

    snew(newnm, nre);
    j = 0;
    for (k = 0; k < nre; k++)
    {
        newnm[k] = strdup(nm[k].name);
        /* Insert dashes in all the names */
        while ((ptr = strchr(newnm[k], ' ')) != NULL)
        {
            *ptr = '-';
        }
        if (bVerbose)
        {
            if (j == 0)
            {
                if (k > 0)
                {
                    fprintf(stderr, "\n");
                }
                bLong = FALSE;
                for (kk = k; kk < k+4; kk++)
                {
                    if (kk < nre && strlen(nm[kk].name) > 14)
                    {
                        bLong = TRUE;
                    }
                }
            }
            else
            {
                fprintf(stderr, " ");
            }
            if (!bLong)
            {
                fprintf(stderr, fm4, k+1, newnm[k]);
                j++;
                if (j == 4)
                {
                    j = 0;
                }
            }
            else
            {
                fprintf(stderr, fm2, k+1, newnm[k]);
                j++;
                if (j == 2)
                {
                    j = 0;
                }
            }
        }
    }
    if (bVerbose)
    {
        fprintf(stderr, "\n\n");
    }

    snew(bE, nre);

    bEOF = FALSE;
    while (!bEOF && (fgets2(buf, STRLEN-1, stdin)))
    {
        /* Remove newlines */
        chomp(buf);

        /* Remove spaces */
        trim(buf);

        /* Empty line means end of input */
        bEOF = (strlen(buf) == 0);
        if (!bEOF)
        {
            ptr = buf;
            do
            {
                if (!bEOF)
                {
                    /* First try to read an integer */
                    nss   = sscanf(ptr, "%d", &nind);
                    if (nss == 1)
                    {
                        /* Zero means end of input */
                        if (nind == 0)
                        {
                            bEOF = TRUE;
                        }
                        else if ((1 <= nind) && (nind <= nre))
                        {
                            bE[nind-1] = TRUE;
                        }
                        else
                        {
                            fprintf(stderr, "number %d is out of range\n", nind);
                        }
                    }
                    else
                    {
                        /* Now try to read a string */
                        i      = strlen(ptr);
                        nmatch = 0;
                        for (nind = 0; nind < nre; nind++)
                        {
                            if (gmx_strcasecmp(newnm[nind], ptr) == 0)
                            {
                                bE[nind] = TRUE;
                                nmatch++;
                            }
                        }
                        if (nmatch == 0)
                        {
                            i      = strlen(ptr);
                            nmatch = 0;
                            for (nind = 0; nind < nre; nind++)
                            {
                                if (gmx_strncasecmp(newnm[nind], ptr, i) == 0)
                                {
                                    bE[nind] = TRUE;
                                    nmatch++;
                                }
                            }
                            if (nmatch == 0)
                            {
                                fprintf(stderr, "String '%s' does not match anything\n", ptr);
                            }
                        }
                    }
                }
                /* Look for the first space, and remove spaces from there */
                if ((ptr = strchr(ptr, ' ')) != NULL)
                {
                    trim(ptr);
                }
            }
            while (!bEOF && (ptr && (strlen(ptr) > 0)));
        }
    }

    snew(set, nre);
    for (i = (*nset) = 0; (i < nre); i++)
    {
        if (bE[i])
        {
            set[(*nset)++] = i;
        }
    }

    sfree(bE);

    if (*nset == 0)
    {
        gmx_fatal(FARGS, "No energy terms selected");
    }

    for (i = 0; (i < nre); i++)
    {
        sfree(newnm[i]);
    }
    sfree(newnm);

    return set;
}
