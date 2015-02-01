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
 * \brief This file defines implementation details for
 * gmx::handleRestart()
 *
 * These handle the details of doing restarts (ie. reading
 * checkpoints, appending output files).
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrunutility
 */

#include "gmxpre.h"

#include "handlerestart-impl.h"

#include <string>
#include <vector>

#include "gromacs/fileio/fileiohandler.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/mpihandler.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/scoped_cptr.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace HandleRestart
{

Impl::Impl(FileIOHandlerInterface   *fileIOHandler,
           MpiHandlerInterface      *mpiHandler,
           ConstCommandLineFilenames mdrunFilenames) :
    fileIOHandler_(fileIOHandler),
    mpiHandler_(mpiHandler),
    mdrunFilenames_(mdrunFilenames)
{
}

Impl::~Impl()
{
}

DataFromCheckpoint Impl::getDataFromCheckpoint(const bool bIsSimMaster,
                                               const bool bIsMultiMaster) const
{
    DataFromCheckpoint data;
    bool               bMdrunIsUsingCheckpointFile =
        opt2bSet("-cpi", mdrunFilenames_.size(), mdrunFilenames_.data());

    if (!bMdrunIsUsingCheckpointFile || !bIsSimMaster)
    {
        return data;
    }

    const char *filename = opt2fn("-cpi", mdrunFilenames_.size(), mdrunFilenames_.data());

    if (!fileIOHandler_->gmx_fexist(filename))
    {
        const char *logFilename = ftp2fn(efLOG, mdrunFilenames_.size(), mdrunFilenames_.data());
        if (fileIOHandler_->gmx_fexist(logFilename))
        {
            GMX_THROW
                (InconsistentInputError
                    (formatString
                        ("Checkpoint file '%s' given with -cpi is not present, but log file '%s' is present. "
                        "Stopping to avoid overwriting output files of a previous run.",
                        filename, logFilename)));
        }
        if (bIsMultiMaster)
        {
            fprintf(stdout, "No previous checkpoint file present, assuming this is a new run.\n");
        }
        return data;
    }

    int                  numOutputFilenames;
    gmx_file_position_t *outputFilenames;

    // A checkpoint file exists and is open in fp
    fileIOHandler_->read_checkpoint_simulation_part_and_filenames(filename, &data.partNumber_,
                                                                  &numOutputFilenames,
                                                                  &outputFilenames);
    // Use RAII for exception safety; fp is already closed
    scoped_cptr<gmx_file_position_t> outputfilesPtr(outputFilenames);
    if (numOutputFilenames == 0)
    {
        GMX_THROW(InternalError("File appending requested, but no output file information is stored in the "
                                "checkpoint file. This should never happen with a valid checkpoint file."));
    }
    // Store the checkpoint filenames more conveniently
    for (int i = 0; i < numOutputFilenames; ++i)
    {
        data.expectedOutputFilenames_.push_back(outputFilenames[i].filename);
    }

    return data;
}

bool
Impl::outputFileExists(const std::string &filename) const
{
    for (size_t i = 0; i < mdrunFilenames_.size(); ++i)
    {
        if (is_output(&mdrunFilenames_[i]) &&
            strcmp(filename.c_str(), mdrunFilenames_[i].fns[0]) == 0)
        {
            return fileIOHandler_->gmx_fexist(filename.c_str());
        }
    }
    return false;
}


std::string
Impl::makeMessageWhenSetsOfFilesDontMatch(const std::vector<std::string> &expectedOutputFilenames,
                                          size_t                          numFilesThatExist) const
{
    std::string message;
    const char *filename = opt2fn("-cpi", mdrunFilenames_.size(), mdrunFilenames_.data());

    message = formatString("Output file appending has been requested,\n"
                           "but some output files listed in the checkpoint file %s\n"
                           "are not present or are named differently by the current program:\n",
                           filename);
    message += "output files present:";
    for (size_t f = 0; f < expectedOutputFilenames.size(); ++f)
    {
        if (outputFileExists(expectedOutputFilenames[f]))
        {
            message += " ";
            message += expectedOutputFilenames[f];
        }
    }
    message += "\noutput files not present or named differently:";
    for (size_t f = 0; f < expectedOutputFilenames.size(); ++f)
    {
        if (!outputFileExists(expectedOutputFilenames[f]))
        {
            message += " ";
            message += expectedOutputFilenames[f];
        }
    }
    message += formatString("\nFile appending requested, but %d of the %d output "
                            "files are not present or are named differently",
                            expectedOutputFilenames.size()-numFilesThatExist,
                            expectedOutputFilenames.size());

    return message;
}

void
Impl::checkAllFilesForAppendingExist(const std::vector<std::string> &expectedOutputFilenames) const
{
    size_t numFilesThatExist = 0;
    GMX_ASSERT(expectedOutputFilenames.size(), "expectedOutputFilenames should contain filenames");

    for (size_t f = 0; f < expectedOutputFilenames.size(); ++f)
    {
        if (outputFileExists(expectedOutputFilenames[f]))
        {
            numFilesThatExist++;
        }
    }
    if (numFilesThatExist == 0)
    {
        /* Zero of the files that are in the checkpoint exist, so we can't
         * append. We don't know what the user actually wants us to do, so
         * give an error. */
        GMX_THROW(InconsistentInputError("Cannot restart from checkpoint file with appending because none of the files used for the old simulation can be found. Either use the old files, don't use mdrun -cpi, or use mdrun -noappend"));
    }
    else if (numFilesThatExist < expectedOutputFilenames.size())
    {
        /* Only some of the files in the checkpoint exist, so we
         * don't know what the user wants to do. */
        GMX_THROW(InconsistentInputError
                      (makeMessageWhenSetsOfFilesDontMatch
                          (expectedOutputFilenames,
                          numFilesThatExist)));
    }
    GMX_ASSERT(numFilesThatExist == expectedOutputFilenames.size(), "code error: more files were found on disk than were used in the run that produced the checkpoint file");
}

bool
Impl::canAppend(const DataFromCheckpoint &dataFromCheckpoint,
                const bool                bTryToAppendFiles) const
{
    if (!bTryToAppendFiles || dataFromCheckpoint.partNumber_ == 0)
    {
        return false;
    }

    checkAllFilesForAppendingExist(dataFromCheckpoint.expectedOutputFilenames_);

    /* The set of mdrun output files matches the checkpoint and exist
       on disk. Check that the .log file is the first one in the
       checkpoint file, and see whether it contains the part suffix,
       in which case we'll have to continue doing that. */
    std::string logFilename = dataFromCheckpoint.expectedOutputFilenames_[0];
    bool        bNeedToAddPart;

    if (Path::getFilename(logFilename).find(ftp2ext(efLOG)) ==
        std::string::npos)
    {
        GMX_THROW(InternalError("File appending requested, but the name of the previous log "
                                "file in the checkpoint file is either invalid not listed first"));
    }
    /* Does the log filename have the part suffix? */
    bNeedToAddPart = (Path::getFilename(logFilename).find(".part") !=
                      std::string::npos);

    /* We can't start appending to files that previously had
       explicit .part names, even if somehow all the other
       conditions are satisfied. */
    return !bNeedToAddPart;
}

} // namespace

} // namespace
