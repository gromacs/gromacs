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

#include "strdb.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

gmx_bool get_a_line(FILE *fp, char line[], int n)
{
    char *line0;
    char *dum;

    snew(line0, n+1);

    do
    {
        if (!fgets(line0, n+1, fp))
        {
            sfree(line0);
            return FALSE;
        }
        dum = strchr(line0, '\n');
        if (dum)
        {
            dum[0] = '\0';
        }
        else if (strlen(line0) == n)
        {
            fprintf(stderr, "Warning: line length exceeds buffer length (%d), data might be corrupted\n", n);
            line0[n-1] = '\0';
        }
        else
        {
            fprintf(stderr, "Warning: file does not end with a newline, last line:\n%s\n",
                    line0);
        }
        dum = strchr(line0, ';');
        if (dum)
        {
            dum[0] = '\0';
        }
        strncpy(line, line0, n);
        dum = line0;
        ltrim(dum);
    }
    while (dum[0] == '\0');

    sfree(line0);
    return TRUE;
}

gmx_bool get_header(char line[], char *header)
{
    char temp[STRLEN], *dum;

    strcpy(temp, line);
    dum = strchr(temp, '[');
    if (dum == NULL)
    {
        return FALSE;
    }
    dum[0] = ' ';
    dum    = strchr(temp, ']');
    if (dum == NULL)
    {
        gmx_fatal(FARGS, "header is not terminated on line:\n'%s'\n", line);
        return FALSE;
    }
    dum[0] = '\0';
    if (sscanf(temp, "%s%*s", header) != 1)
    {
        return FALSE;
    }

    return TRUE;
}

int get_strings(const char *db, char ***strings)
{
    FILE  *in;
    char **ptr;
    char   buf[STRLEN];
    int    i, nstr;

    in = libopen(db);

    if (fscanf(in, "%d", &nstr) != 1)
    {
        gmx_warning("File %s is empty", db);
        gmx_ffclose(in);
        return 0;
    }
    snew(ptr, nstr);
    for (i = 0; (i < nstr); i++)
    {
        if (NULL == fgets2(buf, STRLEN, in))
        {
            gmx_fatal(FARGS, "Cannot read string from buffer");
        }
#ifdef DEBUG
        fprintf(stderr, "Have read: %s\n", buf);
#endif
        ptr[i] = gmx_strdup(buf);
    }
    gmx_ffclose(in);

    *strings = ptr;

    return nstr;
}

int search_str(int nstr, char **str, char *key)
{
    int i;

    /* Linear search */
    for (i = 0; (i < nstr); i++)
    {
        if (gmx_strcasecmp(str[i], key) == 0)
        {
            return i;
        }
    }

    return -1;
}

int fget_lines(FILE *in, char ***strings)
{
    char **ptr;
    char   buf[STRLEN];
    int    i, nstr;
    char  *pret;

    pret = fgets(buf, STRLEN, in);
    if (pret == NULL  || sscanf(buf, "%d", &nstr) != 1)
    {
        gmx_warning("File is empty");
        gmx_ffclose(in);

        return 0;
    }
    snew(ptr, nstr);
    for (i = 0; (i < nstr); i++)
    {
        fgets2(buf, STRLEN, in);
        ptr[i] = gmx_strdup(buf);
    }

    (*strings) = ptr;

    return nstr;
}

int get_lines(const char *db, char ***strings)
{
    FILE *in;
    int   nstr;

    in   = libopen(db);
    nstr = fget_lines(in, strings);
    gmx_ffclose(in);

    return nstr;
}

int get_file(const char *db, char ***strings)
{
    FILE  *in;
    char **ptr = NULL;
    char   buf[STRLEN];
    int    i, nstr, maxi;

    in = libopen(db);

    i = maxi = 0;
    while (fgets2(buf, STRLEN-1, in))
    {
        if (i >= maxi)
        {
            maxi += 50;
            srenew(ptr, maxi);
        }
        ptr[i] = gmx_strdup(buf);
        i++;
    }
    nstr = i;
    gmx_ffclose(in);
    srenew(ptr, nstr);
    *strings = ptr;

    return nstr;
}
