/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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

#include "forcefieldparameters.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/txtdump.h"

static void pr_cmap(FILE* fp, int indent, const char* title, const gmx_cmap_t* cmap_grid, gmx_bool bShowNumbers)
{
    int  j, nelem;
    real dx, idx;

    if (cmap_grid->grid_spacing != 0)
    {
        dx = 360.0 / cmap_grid->grid_spacing;
    }
    else
    {
        dx = 0;
    }
    nelem = cmap_grid->grid_spacing * cmap_grid->grid_spacing;

    if (available(fp, cmap_grid, indent, title))
    {
        fprintf(fp, "%s\n", title);

        for (gmx::index i = 0; i < gmx::ssize(cmap_grid->cmapdata); i++)
        {
            idx = -180.0;
            fprintf(fp, "%8s %8s %8s %8s\n", "V", "dVdx", "dVdy", "d2dV");

            fprintf(fp, "grid[%3zd]={\n", bShowNumbers ? i : -1);

            for (j = 0; j < nelem; j++)
            {
                if ((j % cmap_grid->grid_spacing) == 0)
                {
                    fprintf(fp, "%8.1f\n", idx);
                    idx += dx;
                }

                fprintf(fp, "%8.3f ", cmap_grid->cmapdata[i].cmap[j * 4]);
                fprintf(fp, "%8.3f ", cmap_grid->cmapdata[i].cmap[j * 4 + 1]);
                fprintf(fp, "%8.3f ", cmap_grid->cmapdata[i].cmap[j * 4 + 2]);
                fprintf(fp, "%8.3f\n", cmap_grid->cmapdata[i].cmap[j * 4 + 3]);
            }
            fprintf(fp, "\n");
        }
    }
}

void pr_ffparams(FILE* fp, int indent, const char* title, const gmx_ffparams_t* ffparams, gmx_bool bShowNumbers)
{
    int i;

    indent = pr_title(fp, indent, title);
    pr_indent(fp, indent);
    fprintf(fp, "atnr=%d\n", ffparams->atnr);
    pr_indent(fp, indent);
    fprintf(fp, "ntypes=%d\n", ffparams->numTypes());
    for (i = 0; i < ffparams->numTypes(); i++)
    {
        pr_indent(fp, indent + INDENT);
        fprintf(fp, "functype[%d]=%s, ", bShowNumbers ? i : -1,
                interaction_function[ffparams->functype[i]].name);
        pr_iparams(fp, ffparams->functype[i], ffparams->iparams[i]);
    }
    pr_double(fp, indent, "reppow", ffparams->reppow);
    pr_real(fp, indent, "fudgeQQ", ffparams->fudgeQQ);
    pr_cmap(fp, indent, "cmap", &ffparams->cmap_grid, bShowNumbers);
}
