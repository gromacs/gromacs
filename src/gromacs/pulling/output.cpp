/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "output.h"

#include <cstdio>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "pull_internal.h"

static std::string append_before_extension(const std::string &pathname,
                                           const std::string &to_append)
{
    /* Appends to_append before last '.' in pathname */
    size_t extPos = pathname.find_last_of('.');
    if (extPos == std::string::npos)
    {
        return pathname + to_append;
    }
    else
    {
        return pathname.substr(0, extPos) + to_append +
               pathname.substr(extPos, std::string::npos);
    }
}

static void pull_print_group_x(FILE *out, const ivec dim,
                               const pull_group_work_t *pgrp)
{
    int m;

    for (m = 0; m < DIM; m++)
    {
        if (dim[m])
        {
            fprintf(out, "\t%g", pgrp->x[m]);
        }
    }
}

static void pull_print_coord_dr_components(FILE *out, const ivec dim, const dvec dr)
{
    for (int m = 0; m < DIM; m++)
    {
        if (dim[m])
        {
            fprintf(out, "\t%g", dr[m]);
        }
    }
}

static void pull_print_coord_dr(FILE *out, const pull_coord_work_t *pcrd,
                                gmx_bool bPrintRefValue,
                                gmx_bool bPrintComponents)
{
    double unit_factor = pull_conversion_factor_internal2userinput(&pcrd->params);

    fprintf(out, "\t%g", pcrd->spatialData.value*unit_factor);

    if (bPrintRefValue)
    {
        fprintf(out, "\t%g", pcrd->value_ref*unit_factor);
    }

    if (bPrintComponents)
    {
        pull_print_coord_dr_components(out, pcrd->params.dim, pcrd->spatialData.dr01);
        if (pcrd->params.ngroup >= 4)
        {
            pull_print_coord_dr_components(out, pcrd->params.dim, pcrd->spatialData.dr23);
        }
        if (pcrd->params.ngroup >= 6)
        {
            pull_print_coord_dr_components(out, pcrd->params.dim, pcrd->spatialData.dr45);
        }
    }
}

static void pull_print_x(FILE *out, struct pull_t *pull, double t)
{
    fprintf(out, "%.4f", t);

    for (size_t c = 0; c < pull->coord.size(); c++)
    {
        pull_coord_work_t *pcrd;

        pcrd = &pull->coord[c];

        pull_print_coord_dr(out, pcrd,
                            pull->params.bPrintRefValue && pcrd->params.eType != epullEXTERNAL,
                            pull->params.bPrintComp);

        if (pull->params.bPrintCOM)
        {
            if (pcrd->params.eGeom == epullgCYL)
            {
                pull_print_group_x(out, pcrd->params.dim, &pull->dyna[c]);
            }
            else
            {
                pull_print_group_x(out, pcrd->params.dim,
                                   &pull->group[pcrd->params.group[0]]);
            }
            for (int g = 1; g < pcrd->params.ngroup; g++)
            {
                pull_print_group_x(out, pcrd->params.dim, &pull->group[pcrd->params.group[g]]);
            }
        }
    }
    fprintf(out, "\n");
}

static void pull_print_f(FILE *out, const pull_t *pull, double t)
{
    fprintf(out, "%.4f", t);

    for (const pull_coord_work_t &coord : pull->coord)
    {
        fprintf(out, "\t%g", coord.scalarForce);
    }
    fprintf(out, "\n");
}

void pull_print_output(struct pull_t *pull, int64_t step, double time)
{
    GMX_ASSERT(pull->numExternalPotentialsStillToBeAppliedThisStep == 0, "pull_print_output called before all external pull potentials have been applied");

    if ((pull->params.nstxout != 0) && (step % pull->params.nstxout == 0))
    {
        pull_print_x(pull->out_x, pull, time);
    }

    if ((pull->params.nstfout != 0) && (step % pull->params.nstfout == 0))
    {
        pull_print_f(pull->out_f, pull, time);
    }
}

static void set_legend_for_coord_components(const pull_coord_work_t *pcrd, int coord_index, char **setname, int *nsets_ptr)
{
    /*  Loop over the distance vectors and print their components. Each vector is made up of two consecutive groups. */
    for (int g = 0; g < pcrd->params.ngroup; g += 2)
    {
        /* Loop over the components */
        for (int m  = 0; m < DIM; m++)
        {
            if (pcrd->params.dim[m])
            {
                char legend[STRLEN];

                if (g == 0 && pcrd->params.ngroup <= 2)
                {
                    /*  For the simplest case we print a simplified legend without group indices, just the cooordinate index
                        and which dimensional component it is. */
                    sprintf(legend, "%d d%c", coord_index + 1, 'X' + m);
                }
                else
                {
                    /* Otherwise, print also the group indices (relative to the coordinate, not the global ones). */
                    sprintf(legend, "%d g %d-%d d%c", coord_index + 1, g + 1, g + 2, 'X' + m);
                }

                setname[*nsets_ptr] = gmx_strdup(legend);
                (*nsets_ptr)++;
            }
        }
    }
}

static FILE *open_pull_out(const char *fn, struct pull_t *pull,
                           const gmx_output_env_t *oenv,
                           gmx_bool bCoord,
                           const ContinuationOptions &continuationOptions)
{
    FILE  *fp;
    int    nsets, m;
    char **setname, buf[50];

    if (continuationOptions.appendFiles)
    {
        fp = gmx_fio_fopen(fn, "a+");
    }
    else
    {
        fp = gmx_fio_fopen(fn, "w+");
        if (bCoord)
        {
            sprintf(buf, "Position (nm%s)", pull->bAngle ? ", deg" : "");
            xvgr_header(fp, "Pull COM",  "Time (ps)", buf,
                        exvggtXNY, oenv);
        }
        else
        {
            sprintf(buf, "Force (kJ/mol/nm%s)", pull->bAngle ? ", kJ/mol/rad" : "");
            xvgr_header(fp, "Pull force", "Time (ps)", buf,
                        exvggtXNY, oenv);
        }

        /* With default mdp options only the actual coordinate value is printed (1),
         * but optionally the reference value (+ 1),
         * the group COMs for all the groups (+ ngroups_max*DIM)
         * and the components of the distance vectors can be printed (+ (ngroups_max/2)*DIM).
         */
        snew(setname, pull->coord.size()*(1 + 1 + c_pullCoordNgroupMax*DIM + c_pullCoordNgroupMax/2*DIM));

        nsets = 0;
        for (size_t c = 0; c < pull->coord.size(); c++)
        {
            if (bCoord)
            {
                /* The order of this legend should match the order of printing
                 * the data in print_pull_x above.
                 */

                /* The pull coord distance */
                sprintf(buf, "%zu", c+1);
                setname[nsets] = gmx_strdup(buf);
                nsets++;
                if (pull->params.bPrintRefValue &&
                    pull->coord[c].params.eType != epullEXTERNAL)
                {
                    sprintf(buf, "%zu ref", c+1);
                    setname[nsets] = gmx_strdup(buf);
                    nsets++;
                }
                if (pull->params.bPrintComp)
                {
                    set_legend_for_coord_components(&pull->coord[c], c, setname, &nsets);
                }

                if (pull->params.bPrintCOM)
                {
                    for (int g = 0; g < pull->coord[c].params.ngroup; g++)
                    {
                        /* Legend for reference group position */
                        for (m = 0; m < DIM; m++)
                        {
                            if (pull->coord[c].params.dim[m])
                            {
                                sprintf(buf, "%zu g %d %c", c+1, g + 1, 'X'+m);
                                setname[nsets] = gmx_strdup(buf);
                                nsets++;
                            }
                        }
                    }
                }
            }
            else
            {
                /* For the pull force we always only use one scalar */
                sprintf(buf, "%zu", c+1);
                setname[nsets] = gmx_strdup(buf);
                nsets++;
            }
        }
        if (nsets > 1)
        {
            xvgr_legend(fp, nsets, setname, oenv);
        }
        for (int c = 0; c < nsets; c++)
        {
            sfree(setname[c]);
        }
        sfree(setname);
    }

    return fp;
}

void init_pull_output_files(pull_t                    *pull,
                            int                        nfile,
                            const t_filenm             fnm[],
                            const gmx_output_env_t    *oenv,
                            const ContinuationOptions &continuationOptions)
{
    /* Check for px and pf filename collision, if we are writing
       both files */
    std::string px_filename, pf_filename;
    std::string px_appended, pf_appended;
    try
    {
        px_filename  = std::string(opt2fn("-px", nfile, fnm));
        pf_filename  = std::string(opt2fn("-pf", nfile, fnm));
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    if ((pull->params.nstxout != 0) &&
        (pull->params.nstfout != 0) &&
        (px_filename == pf_filename))
    {
        if (!opt2bSet("-px", nfile, fnm) && !opt2bSet("-pf", nfile, fnm))
        {
            /* We are writing both pull files but neither set directly. */
            try
            {
                px_appended   = append_before_extension(px_filename, "_pullx");
                pf_appended   = append_before_extension(pf_filename, "_pullf");
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            pull->out_x = open_pull_out(px_appended.c_str(), pull, oenv,
                                        TRUE, continuationOptions);
            pull->out_f = open_pull_out(pf_appended.c_str(), pull, oenv,
                                        FALSE, continuationOptions);
            return;
        }
        else
        {
            /* If at least one of -px and -pf is set but the filenames are identical: */
            gmx_fatal(FARGS, "Identical pull_x and pull_f output filenames %s",
                      px_filename.c_str());
        }
    }
    if (pull->params.nstxout != 0)
    {
        pull->out_x = open_pull_out(opt2fn("-px", nfile, fnm), pull, oenv,
                                    TRUE, continuationOptions);
    }
    if (pull->params.nstfout != 0)
    {
        pull->out_f = open_pull_out(opt2fn("-pf", nfile, fnm), pull, oenv,
                                    FALSE, continuationOptions);
    }
}
