/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file declares functions for mdrun to call to manage the
 * details of doing a restart (ie. reading checkpoints, appending
 * output files).
 *
 * \todo Clean up the error-prone logic here. Add doxygen.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrunutility
 */

#include "gmxpre.h"

#include "handlerestart.h"

#include <string.h>

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/checkpoint.h"
#include "gromacs/legacyheaders/main.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/smalloc.h"

/*! \brief Search for \p fnm_cp in fnm and return true iff found
 *
 * \todo This could be implemented sanely with a for loop. */
static gmx_bool exist_output_file(const char *fnm_cp, int nfile, const t_filenm fnm[])
{
    int i;

    /* Check if the output file name stored in the checkpoint file
     * is one of the output file names of mdrun.
     */
    i = 0;
    while (i < nfile &&
           !(is_output(&fnm[i]) && strcmp(fnm_cp, fnm[i].fns[0]) == 0))
    {
        i++;
    }

    return (i < nfile && gmx_fexist(fnm_cp));
}

/*! \brief Support handling restarts
 *
 * \todo Clean this up (next patch)
 *
 * Read just the simulation 'generation' and with bTryToAppendFiles check files.
 * This is is needed at the beginning of mdrun,
 * to be able to rename the logfile correctly.
 * When file appending is requested, checks which output files are present,
 * and returns TRUE/FALSE in bDoAppendFiles if all or none are present.
 * If only some output files are present, give a fatal error.
 * When bDoAppendFiles is TRUE upon return, bAddPart will tell whether the simulation part
 * needs to be added to the output file name.
 *
 * This routine cannot print tons of data, since it is called before
 * the log file is opened. */
static void
read_checkpoint_data(const char *filename, int *simulation_part,
                     t_commrec *cr,
                     gmx_bool bTryToAppendFiles,
                     int nfile, const t_filenm fnm[],
                     const char *part_suffix,
                     gmx_bool *bAddPart,
                     gmx_bool *bDoAppendFiles)
{
    t_fileio            *fp;
    int                  nfiles;
    gmx_file_position_t *outputfiles;
    int                  nexist, f;
    char                *fn, suf_up[STRLEN];

    *bDoAppendFiles = FALSE;

    if (SIMMASTER(cr))
    {
        if (!gmx_fexist(filename) || (!(fp = gmx_fio_open(filename, "r")) ))
        {
            *simulation_part = 0;
        }
        else
        {
            read_checkpoint_simulation_part_and_filenames(fp,
                                                          simulation_part,
                                                          &nfiles,
                                                          &outputfiles);

            if (bTryToAppendFiles)
            {
                nexist = 0;
                for (f = 0; f < nfiles; f++)
                {
                    if (exist_output_file(outputfiles[f].filename, nfile, fnm))
                    {
                        nexist++;
                    }
                }
                if (nexist == nfiles)
                {
                    *bDoAppendFiles = bTryToAppendFiles;
                }
                else if (nexist > 0)
                {
                    fprintf(stderr,
                            "Output file appending has been requested,\n"
                            "but some output files listed in the checkpoint file %s\n"
                            "are not present or are named differently by the current program:\n",
                            filename);
                    fprintf(stderr, "output files present:");
                    for (f = 0; f < nfiles; f++)
                    {
                        if (exist_output_file(outputfiles[f].filename,
                                              nfile, fnm))
                        {
                            fprintf(stderr, " %s", outputfiles[f].filename);
                        }
                    }
                    fprintf(stderr, "\n");
                    fprintf(stderr, "output files not present or named differently:");
                    for (f = 0; f < nfiles; f++)
                    {
                        if (!exist_output_file(outputfiles[f].filename,
                                               nfile, fnm))
                        {
                            fprintf(stderr, " %s", outputfiles[f].filename);
                        }
                    }
                    fprintf(stderr, "\n");

                    gmx_fatal(FARGS, "File appending requested, but %d of the %d output files are not present or are named differently", nfiles-nexist, nfiles);
                }
            }

            if (*bDoAppendFiles)
            {
                if (nfiles == 0)
                {
                    gmx_fatal(FARGS, "File appending requested, but no output file information is stored in the checkpoint file");
                }
                fn = outputfiles[0].filename;
                if (strlen(fn) < 4 ||
                    gmx_strcasecmp(fn+strlen(fn)-4, ftp2ext(efLOG)) == 0)
                {
                    gmx_fatal(FARGS, "File appending requested, but the log file is not the first file listed in the checkpoint file");
                }
                /* Set bAddPart to whether the suffix string '.part' is present
                 * in the log file name.
                 */
                strcpy(suf_up, part_suffix);
                upstring(suf_up);
                *bAddPart = (strstr(fn, part_suffix) != NULL ||
                             strstr(fn, suf_up) != NULL);
            }

            sfree(outputfiles);
        }
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(*simulation_part), simulation_part, cr);

        if (*simulation_part > 0 && bTryToAppendFiles)
        {
            gmx_bcast(sizeof(*bDoAppendFiles), bDoAppendFiles, cr);
            gmx_bcast(sizeof(*bAddPart), bAddPart, cr);
        }
    }
}

/* This routine cannot print tons of data, since it is called before the log file is opened. */
void
handleRestart(t_commrec *cr,
              gmx_bool   bTryToAppendFiles,
              const int  NFILE,
              t_filenm   fnm[],
              gmx_bool  *bDoAppendFiles,
              gmx_bool  *bStartFromCpt)
{
    gmx_bool        bAddPart;
    int             sim_part, sim_part_fn;
    const char     *part_suffix = ".part";
    FILE           *fpmulti;

    bAddPart = !bTryToAppendFiles;

    /* Check if there is ANY checkpoint file available */
    sim_part    = 1;
    sim_part_fn = sim_part;
    if (opt2bSet("-cpi", NFILE, fnm))
    {
        read_checkpoint_data(opt2fn_master("-cpi", NFILE, fnm, cr),
                             &sim_part_fn, cr,
                             bTryToAppendFiles, NFILE, fnm,
                             part_suffix, &bAddPart,
                             bDoAppendFiles);
        if (sim_part_fn == 0 && MULTIMASTER(cr))
        {
            fprintf(stdout, "No previous checkpoint file present, assuming this is a new run.\n");
        }
        else
        {
            sim_part = sim_part_fn + 1;
        }

        if (MULTISIM(cr) && MASTER(cr))
        {
            if (MULTIMASTER(cr))
            {
                /* Log file is not yet available, so if there's a
                 * problem we can only write to stderr. */
                fpmulti = stderr;
            }
            else
            {
                fpmulti = NULL;
            }
            check_multi_int(fpmulti, cr->ms, sim_part, "simulation part", TRUE);
        }
    }
    else
    {
        *bDoAppendFiles = FALSE;
    }
    *bStartFromCpt = sim_part > 1;

    if (!*bDoAppendFiles)
    {
        sim_part_fn = sim_part;
    }

    if (bAddPart)
    {
        char suffix[STRLEN];

        /* Rename all output files (except checkpoint files) */
        /* create new part name first (zero-filled) */
        sprintf(suffix, "%s%04d", part_suffix, sim_part_fn);

        add_suffix_to_output_names(fnm, NFILE, suffix);
        if (MULTIMASTER(cr))
        {
            fprintf(stdout, "Checkpoint file is from part %d, new output files will be suffixed '%s'.\n", sim_part-1, suffix);
        }
    }
}
