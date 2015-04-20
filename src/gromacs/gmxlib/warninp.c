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

#include "gromacs/legacyheaders/warninp.h"

#include <string.h>

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

typedef struct warninp {
    gmx_bool bAllowWarnings;
    int      nwarn_note;
    int      nwarn_warn;
    int      nwarn_error;
    int      maxwarn;
    int      lineno;
    char     filenm[256];
} t_warninp;

warninp_t init_warning(gmx_bool bAllowWarnings, int maxwarning)
{
    warninp_t wi;

    snew(wi, 1);

    wi->bAllowWarnings = bAllowWarnings;
    wi->maxwarn        = maxwarning;
    wi->nwarn_note     = 0;
    wi->nwarn_warn     = 0;
    wi->nwarn_error    = 0;
    strcpy(wi->filenm, "unknown");
    wi->lineno         = 0;

    return wi;
}

void set_warning_line(warninp_t wi, const char *s, int line)
{
    if (s != NULL)
    {
        strcpy(wi->filenm, s);
    }
    wi->lineno = line;
}

int get_warning_line(warninp_t wi)
{
    return wi->lineno;
}

const char *get_warning_file(warninp_t wi)
{
    return wi->filenm;
}

static void low_warning(warninp_t wi, const char *wtype, int n, const char *s)
{
#define indent 2
    char *temp, *temp2;
    int   i;

    if (s == NULL)
    {
        s = "Empty error message.";
    }
    snew(temp, strlen(s)+indent+1);
    for (i = 0; i < indent; i++)
    {
        temp[i] = ' ';
    }
    temp[indent] = '\0';
    strcat(temp, s);
    temp2 = wrap_lines(temp, 78-indent, indent, FALSE);
    if (strlen(wi->filenm) > 0)
    {
        if (wi->lineno != -1)
        {
            fprintf(stderr, "\n%s %d [file %s, line %d]:\n%s\n\n",
                    wtype, n, wi->filenm, wi->lineno, temp2);
        }
        else
        {
            fprintf(stderr, "\n%s %d [file %s]:\n%s\n\n",
                    wtype, n, wi->filenm, temp2);
        }
    }
    else
    {
        fprintf(stderr, "\n%s %d:\n%s\n\n", wtype, n, temp2);
    }
    sfree(temp);
    sfree(temp2);
}

void warning(warninp_t wi, const char *s)
{
    if (wi->bAllowWarnings)
    {
        wi->nwarn_warn++;
        low_warning(wi, "WARNING", wi->nwarn_warn, s);
    }
    else
    {
        warning_error(wi, s);
    }
}

void warning_note(warninp_t wi, const char *s)
{
    wi->nwarn_note++;
    low_warning(wi, "NOTE", wi->nwarn_note, s);
}

void warning_error(warninp_t wi, const char *s)
{
    wi->nwarn_error++;
    low_warning(wi, "ERROR", wi->nwarn_error, s);
}

static void print_warn_count(const char *type, int n)
{
    if (n > 0)
    {
        fprintf(stderr, "\nThere %s %d %s%s\n",
                (n == 1) ? "was" : "were", n, type, (n == 1) ? "" : "s");
    }
}

void check_warning_error(warninp_t wi, int f_errno, const char *file, int line)
{
    if (wi->nwarn_error > 0)
    {
        print_warn_count("note", wi->nwarn_note);
        print_warn_count("warning", wi->nwarn_warn);

        gmx_fatal(f_errno, file, line, "There %s %d error%s in input file(s)",
                  (wi->nwarn_error == 1) ? "was" : "were", wi->nwarn_error,
                  (wi->nwarn_error == 1) ? ""    : "s");
    }
}

void done_warning(warninp_t wi, int f_errno, const char *file, int line)
{
    print_warn_count("note", wi->nwarn_note);
    print_warn_count("warning", wi->nwarn_warn);

    check_warning_error(wi, f_errno, file, line);

    if (wi->maxwarn >= 0 && wi->nwarn_warn > wi->maxwarn)
    {
        gmx_fatal(f_errno, file, line,
                  "Too many warnings (%d), %s terminated.\n"
                  "If you are sure all warnings are harmless, use the -maxwarn option.",
                  wi->nwarn_warn, ShortProgram());
    }

    sfree(wi);
}

void _too_few(warninp_t wi, const char *fn, int line)
{
    char buf[STRLEN];

    sprintf(buf,
            "Too few parameters on line (source file %s, line %d)",
            fn, line);
    warning(wi, buf);
}

void _incorrect_n_param(warninp_t wi, const char *fn, int line)
{
    char buf[STRLEN];

    sprintf(buf,
            "Incorrect number of parameters on line (source file %s, line %d)",
            fn, line);
    warning(wi, buf);
}
