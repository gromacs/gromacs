/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "vecdump.h"

#include <cstdio>
#include <cstdlib>

#include "gromacs/math/vec.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/txtdump.h"

void pr_ivec(FILE *fp, int indent, const char *title, const int vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%d\n", title, bShowNumbers ? i : -1, vec[i]);
        }
    }
}

void pr_ivec_block(FILE *fp, int indent, const char *title, const int vec[], int n, gmx_bool bShowNumbers)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        i      = 0;
        while (i < n)
        {
            j = i+1;
            while (j < n && vec[j] == vec[j-1]+1)
            {
                j++;
            }
            /* Print consecutive groups of 3 or more as blocks */
            if (j - i < 3)
            {
                while (i < j)
                {
                    pr_indent(fp, indent);
                    fprintf(fp, "%s[%d]=%d\n",
                            title, bShowNumbers ? i : -1, vec[i]);
                    i++;
                }
            }
            else
            {
                pr_indent(fp, indent);
                fprintf(fp, "%s[%d,...,%d] = {%d,...,%d}\n",
                        title,
                        bShowNumbers ? i : -1,
                        bShowNumbers ? j-1 : -1,
                        vec[i], vec[j-1]);
                i = j;
            }
        }
    }
}

void pr_bvec(FILE *fp, int indent, const char *title, const gmx_bool vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%s\n", title, bShowNumbers ? i : -1,
                    gmx::boolToString(vec[i]));
        }
    }
}

void pr_ivecs(FILE *fp, int indent, const char *title, const ivec vec[], int n, gmx_bool bShowNumbers)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={", title, bShowNumbers ? i : -1);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, "%d", vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}

template <typename T>
static void printRealVector(FILE *fp, int indent, const char *title, const T vec[], int n, gmx_bool bShowNumbers)
{
    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (int i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%12.5e\n", title, bShowNumbers ? i : -1, vec[i]);
        }
    }
}

void pr_rvec(FILE *fp, int indent, const char *title, const real vec[], int n, gmx_bool bShowNumbers)
{
    printRealVector<real>(fp, indent, title, vec, n, bShowNumbers);
}

void pr_fvec(FILE *fp, int indent, const char *title, const float vec[], int n, gmx_bool bShowNumbers)
{
    printRealVector<float>(fp, indent, title, vec, n, bShowNumbers);
}

void pr_dvec(FILE *fp, int indent, const char *title, const double vec[], int n, gmx_bool bShowNumbers)
{
    printRealVector<double>(fp, indent, title, vec, n, bShowNumbers);
}


void pr_rvecs_len(FILE *fp, int indent, const char *title, const rvec vec[], int n)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, "%12.5e", vec[i][j]);
            }
            fprintf(fp, "} len=%12.5e\n", norm(vec[i]));
        }
    }
}

void pr_rvecs(FILE *fp, int indent, const char *title, const rvec vec[], int n)
{
    const char *fshort = "%12.5e";
    const char *flong  = "%15.8e";
    const char *format;
    int         i, j;

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
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, format, vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}


void pr_rvecs_of_dim(FILE *fp, int indent, const char *title, const rvec vec[], int n, int dim)
{
    const char *fshort = "%12.5e";
    const char *flong  = "%15.8e";
    const char *format;
    int         i, j;

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
                fprintf(fp, format, vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}
