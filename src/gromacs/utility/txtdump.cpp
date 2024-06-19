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

#include "gromacs/utility/txtdump.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/* This file is completely threadsafe - please keep it that way! */

#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "gromacs/utility/cstringutil.h"

int pr_indent(FILE* fp, int n)
{
    for (int i = 0; i < n; i++)
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
    if (available(fp, vec, indent, title))
    {
        pr_indent(fp, indent);
        fprintf(fp, "%s:\t", title);
        for (int i = 0; i < n; i++)
        {
            fprintf(fp, "  %10g", vec[i]);
        }
        fprintf(fp, "\n");
    }
}

void pr_doubles(FILE* fp, int indent, const char* title, const double* vec, int n)
{
    if (available(fp, vec, indent, title))
    {
        pr_indent(fp, indent);
        fprintf(fp, "%s:\t", title);
        for (int i = 0; i < n; i++)
        {
            fprintf(fp, "  %10g", vec[i]);
        }
        fprintf(fp, "\n");
    }
}

void pr_reals_of_dim(FILE* fp, int indent, const char* title, const real* vec, int n, int dim)
{
    const char* fshort = "%12.5e";
    const char* flong  = "%15.8e";
    const char* format = (getenv("GMX_PRINT_LONGFORMAT") != nullptr) ? flong : fshort;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, dim);
        for (int i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (int j = 0; j < dim; j++)
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

void pr_strings(FILE* fp, int indent, const char* title, const char* const* const* nm, int n, gmx_bool bShowNumbers)
{
    if (available(fp, nm, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (int i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\"}\n", title, bShowNumbers ? i : -1, *(nm[i]));
        }
    }
}
