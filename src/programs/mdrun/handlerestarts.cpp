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
 * details of doing restarts (ie. reading checkpoints, appending
 * output files).
 *
 * \todo Clean up the error-prone logic here. Add doxygen.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrun
 */

#include "gmxpre.h"

#include "handlerestarts.h"

#include <string.h>

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/checkpoint.h"
#include "gromacs/legacyheaders/main.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/smalloc.h"

/*! \brief Handle reading of mdrun -cpi to get basic info
 *
 * If a checkpoint file was supplied to mdrun and this is a master
 * rank, read that file to get the old simulation part number, and the
 * names of files used in that run.
 *
 * Does not alter the output parameters if there's no checkpoint file,
 * or it cannot be read, or it should not be read on this non-master
 * rank. This anticipates global broadcast of simulation part to all
 * non-master ranks in a simulation.
 *
 * \param[in]    bIsSimMaster           Flag to permit only one rank to write to read the checkpoint file; should normally be set to SIMMASTER(cr)
 * \param[in]    bIsMultiMaster         Flag to permit only one rank to write to stdout; should normally be set to MULTIMASTER(cr)
 * \param[in]    NFILE                  Size of fnm struct
 * \param[in]    fnm                    Filename parameters to mdrun
 * \param[out]   partForThisSimulation  The number of the simulation part that will be run
 * \param[out]   nfiles                 Number of output files from the previous run
 * \param[out]   outputfiles            Pointer to array of output file names from the previous run.
 *
 * The outputfiles pointer is allocated in this function if the
 * checkpoint is actually read. The caller has responsibility for
 * calling sfree on \p outputfiles when it is non-NULL.
 */
static void
getDataFromCheckpoint(const gmx_bool        bIsSimMaster,
                      const gmx_bool        bIsMultiMaster,
                      const int             NFILE,
                      const t_filenm        fnm[],
                      int                  *partForThisSimulation,
                      int                  *nfiles,
                      gmx_file_position_t **outputfiles)
{
    t_fileio *fp;

    if (!opt2bSet("-cpi", NFILE, fnm) || !bIsSimMaster)
    {
        return;
    }

    const char *filename = opt2fn("-cpi", NFILE, fnm);

    if (!gmx_fexist(filename) || (!(fp = gmx_fio_open(filename, "r")) ))
    {
        if (bIsMultiMaster)
        {
            // TODO This might be a user input error if a .log file also exists, but currently we don't check for that.
            fprintf(stdout, "No previous checkpoint file present, assuming this is a new run.\n");
        }
        return;
    }

    int partFromCpt;
    // A checkpoint file exists and is open in fp
    read_checkpoint_simulation_part_and_filenames(fp, &partFromCpt,
                                                  nfiles, outputfiles);

    // The new part number is the one from the checkpoint file, plus one
    *partForThisSimulation = partFromCpt + 1;
}

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

/*! \brief Check whether all the files for appending exist
 *
 * Returns TRUE/FALSE if they all/none exist. Gives fatal error
 * if only some of them exist. */
static gmx_bool checkIfFilesToAppendExist(const char                *checkpointFilename,
                                          const int                  nfile,
                                          const t_filenm             fnm[],
                                          const int                  nfiles,
                                          const gmx_file_position_t *outputfiles)
{
    int nexist = 0;

    for (int f = 0; f < nfiles; f++)
    {
        if (exist_output_file(outputfiles[f].filename, nfile, fnm))
        {
            nexist++;
        }
    }
    if (nexist == nfiles)
    {
        return TRUE;
    }
    if (nexist > 0)
    {
        fprintf(stderr,
                "Output file appending has been requested,\n"
                "but some output files listed in the checkpoint file %s\n"
                "are not present or are named differently by the current program:\n",
                checkpointFilename);
        fprintf(stderr, "output files present:");
        for (int f = 0; f < nfiles; f++)
        {
            if (exist_output_file(outputfiles[f].filename,
                                  nfile, fnm))
            {
                fprintf(stderr, " %s", outputfiles[f].filename);
            }
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "output files not present or named differently:");
        for (int f = 0; f < nfiles; f++)
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
    return FALSE;
}

/*! \brief Return whether the conditions on this rank permit appending
 *
 * If using mdrun -append and the old simulation part number is > 1,
 * appending might occur. Check that the files on disk form a
 * complete/empty set, and prepare to return true/false
 * respectively. If only some output files are present, give a fatal
 * error. However, if the .log filename in the checkpoint file
 * contains \c part_suffix, then do not append any files.
 *
 * Gives fatal errors if \p nfiles and \p outputfiles are malformed.
 *
 * \param[in]    bTryToAppendFiles      Whether mdrun -append was used
 * \param[in]    NFILE                  Size of fnm struct
 * \param[in]    fnm                    Filename parameters to mdrun
 * \param[in]    part_suffix            Suffix string used to produce filenames like "md.part0003.log"
 * \param[in]    partForThisSimulation  The number of the simulation part that will be run
 * \param[in]    nfiles                 Number of output files from the previous run
 * \param[in]    outputfiles            Pointer to array of output file names from the previous run
 *
 * \return On master rank, whether the conditions are correct for the
 * run to actually do appending; otherwise false. It is expected that
 * non-master ranks will get the information via broadcast.
 */
static gmx_bool
chooseWhetherToAppend(const gmx_bool             bTryToAppendFiles,
                      const int                  NFILE,
                      const t_filenm             fnm[],
                      const char                *part_suffix,
                      const int                  partForThisSimulation,
                      const int                  nfiles,
                      const gmx_file_position_t *outputfiles)
{
    gmx_bool bDoAppendFiles;

    if (!bTryToAppendFiles || partForThisSimulation == 1)
    {
        return FALSE;
    }

    if (nfiles == 0)
    {
        gmx_fatal(FARGS, "File appending requested, but no output file information is stored in the checkpoint file");
    }
    bDoAppendFiles = checkIfFilesToAppendExist(opt2fn("-cpi", NFILE, fnm),
                                               NFILE, fnm,
                                               nfiles, outputfiles);

    if (bDoAppendFiles)
    {
        /* Check that the .log file is the first one in the checkpoint
           file, and see whether it contains the part suffix, in which
           case we'll have to continue doing that. */
        char        suf_up[STRLEN];
        const char *fn = outputfiles[0].filename;
        gmx_bool    bNeedToAddPart;

        if (strlen(fn) < 4 ||
            gmx_strcasecmp(fn+strlen(fn)-4, ftp2ext(efLOG)) == 0)
        {
            gmx_fatal(FARGS, "File appending requested, but the name of the previous log file in the checkpoint file is either invalid not listed first");
        }
        /* Set bNeedToAddPart to whether the suffix string '.part' is
         * present in the log file name. */
        strcpy(suf_up, part_suffix);
        upstring(suf_up);
        bNeedToAddPart = (strstr(fn, part_suffix) != NULL ||
                          strstr(fn, suf_up) != NULL);
        /* We can't start appending to files that previously had
           explicit .part names, even if somehow all the other
           conditions are satisfied. */
        bDoAppendFiles = !bNeedToAddPart;
    }

    return bDoAppendFiles;
}

/* This routine cannot print tons of data, since it is called before the log file is opened. */
void
handleRestarts(t_commrec *cr,
               gmx_bool   bTryToAppendFiles,
               const int  NFILE,
               t_filenm   fnm[],
               gmx_bool  *bDoAppendFiles,
               gmx_bool  *bStartFromCpt)
{
    const char *part_suffix = ".part";
    // Set safe default
    int         partForThisSimulation = 1;

    {
        int                  nfiles;
        gmx_file_position_t *outputfiles = NULL;

        getDataFromCheckpoint(SIMMASTER(cr),
                              MULTIMASTER(cr),
                              NFILE,
                              fnm,
                              &partForThisSimulation,
                              &nfiles,
                              &outputfiles);

        *bDoAppendFiles = chooseWhetherToAppend(bTryToAppendFiles,
                                                NFILE, fnm, part_suffix,
                                                partForThisSimulation,
                                                nfiles,
                                                outputfiles);
        if (outputfiles)
        {
            sfree(outputfiles);
        }
    }

    /* Communicate results of checks */
    if (PAR(cr))
    {
        gmx_bcast(sizeof(partForThisSimulation), &partForThisSimulation, cr);

        if (partForThisSimulation > 1 && bTryToAppendFiles)
        {
            /* In this case, we need to coordinate whether we'll do
               appending. Otherwise, the default for bDoAppendFiles
               will do. */
            gmx_bcast(sizeof(*bDoAppendFiles), bDoAppendFiles, cr);
        }
    }
    *bStartFromCpt = partForThisSimulation > 1;

    if (*bDoAppendFiles && MULTISIM(cr) && MASTER(cr))
    {
        /* Log file is not yet available, so if there's a
         * problem we can only write to stderr. */
        check_multi_int(stderr, cr->ms, partForThisSimulation, "simulation part", TRUE);
    }

    if (*bStartFromCpt && !*bDoAppendFiles)
    {
        /* Now modify filenames on all ranks (which happens only when
           we're not appending) */
        char suffix[STRLEN];

        /* Rename all output files (except checkpoint files) using the
           suffix containing the new part number. */
        sprintf(suffix, "%s%04d", part_suffix, partForThisSimulation);

        add_suffix_to_output_names(fnm, NFILE, suffix);
        if (MULTIMASTER(cr))
        {
            fprintf(stdout, "Checkpoint file is from part %d, new output files will be suffixed '%s'.\n", partForThisSimulation-1, suffix);
        }
    }
}
