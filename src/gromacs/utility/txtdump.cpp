/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

/* This file is completely threadsafe - please keep it that way! */

#include "txtdump.h"

#include <cstdio>
#include <cstdlib>

#include "gromacs/utility/cstringutil.h"

int pr_indent(FILE* fp, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        fprintf(fp, " ");
    }
    return n;
}

bool available(FILE* fp, const void* p, int indent, const char* title)
{
    if (!p)
    {
        if (indent > 0)
        {
            pr_indent(fp, indent);
        }
        fprintf(fp, "%s: not available\n", title);
    }
    return (p != nullptr);
}

int pr_title(FILE* fp, int indent, const char* title)
{
    pr_indent(fp, indent);
    fprintf(fp, "%s:\n", title);
    return (indent + INDENT);
}

int pr_title_n(FILE* fp, int indent, const char* title, int n)
{
    pr_indent(fp, indent);
    fprintf(fp, "%s (%d):\n", title, n);
    return (indent + INDENT);
}

int pr_title_nxn(FILE* fp, int indent, const char* title, int n1, int n2)
{
    pr_indent(fp, indent);
    fprintf(fp, "%s (%dx%d):\n", title, n1, n2);
    return (indent + INDENT);
}

void pr_reals(FILE* fp, int indent, const char* title, const real* vec, int n)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        pr_indent(fp, indent);
        fprintf(fp, "%s:\t", title);
        for (i = 0; i < n; i++)
        {
            fprintf(fp, "  %10g", vec[i]);
        }
        fprintf(fp, "\n");
    }
}

void pr_doubles(FILE* fp, int indent, const char* title, const double* vec, int n)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        pr_indent(fp, indent);
        fprintf(fp, "%s:\t", title);
        for (i = 0; i < n; i++)
        {
            fprintf(fp, "  %10g", vec[i]);
        }
        fprintf(fp, "\n");
    }
}

void pr_reals_of_dim(FILE* fp, int indent, const char* title, const real* vec, int n, int dim)
{
    int         i, j;
    const char* fshort = "%12.5e";
    const char* flong  = "%15.8e";
    const char* format;

    if (getenv("GMX_PRINT_LONGFORMAT") != nullptr)
    {
        format = flong;
    }
    else
    {
        format = fshort;
    }

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, dim);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < dim; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, format, vec[i * dim + j]);
            }
            fprintf(fp, "}\n");
        }
    }
}

void pr_int(FILE* fp, int indent, const char* title, int i)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %d\n", title, i);
}

void pr_int64(FILE* fp, int indent, const char* title, int64_t i)
{
    char buf[STEPSTRSIZE];

    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %s\n", title, gmx_step_str(i, buf));
}

void pr_real(FILE* fp, int indent, const char* title, real r)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %g\n", title, r);
}

void pr_double(FILE* fp, int indent, const char* title, double d)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %g\n", title, d);
}

void pr_str(FILE* fp, int indent, const char* title, const char* s)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %s\n", title, s);
}

void pr_strings(FILE* fp, int indent, const char* title, char*** nm, int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, nm, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\"}\n", title, bShowNumbers ? i : -1, *(nm[i]));
        }
    }
}
