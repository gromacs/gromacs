/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2012,2013,2014,2015,2017,2018, by the GROMACS development team, led by
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

#include "fflibutil.h"

#include <cstring>

#include <string>
#include <vector>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/datafilefinder.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/directoryenumerator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

const char *fflib_forcefield_dir_ext()
{
    return ".ff";
}

const char *fflib_forcefield_itp()
{
    return "forcefield.itp";
}

const char *fflib_forcefield_doc()
{
    return "forcefield.doc";
}

void fflib_filename_base(const char *filename, char *filebase, int maxlen)
{
    const char *cptr;
    char       *ptr;

    cptr = strrchr(filename, DIR_SEPARATOR);
    if (cptr != nullptr)
    {
        /* Skip the separator */
        cptr += 1;
    }
    else
    {
        cptr = filename;
    }
    if (strlen(filename) >= static_cast<size_t>(maxlen))
    {
        gmx_fatal(FARGS, "filename is longer (%zu) than maxlen (%d)",
                  strlen(filename), maxlen);
    }
    strcpy(filebase, cptr);
    /* Remove the extension */
    ptr = strrchr(filebase, '.');
    if (ptr != nullptr)
    {
        ptr[0] = '\0';
    }
}

std::vector<std::string> fflib_search_file_end(const char *ffdir,
                                               const char *file_end,
                                               bool        bFatalError)
{
    try
    {
        std::string              ffdirFull(gmx::getLibraryFileFinder().findFile(ffdir));
        std::vector<std::string> result
            = gmx::DirectoryEnumerator::enumerateFilesWithExtension(
                        ffdirFull.c_str(), file_end, true);
        if (result.empty() && bFatalError)
        {
            std::string message
                = gmx::formatString("Could not find any files ending on '%s' "
                                    "in the force field directory '%s'",
                                    file_end, ffdir);
            GMX_THROW(gmx::InvalidInputError(message));
        }
        for (std::string &filename : result)
        {
            filename = gmx::Path::join(ffdir, filename);
        }
        return result;
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

std::vector<gmx::DataFileInfo> fflib_enumerate_forcefields()
{
    const char *const              dirend   = fflib_forcefield_dir_ext();
    const char *const              filename = fflib_forcefield_itp();
    std::vector<gmx::DataFileInfo> candidates
        = gmx::getLibraryFileFinder().enumerateFiles(
                    gmx::DataFileOptions(dirend)
                        .throwIfNotFound(false));

    std::vector<gmx::DataFileInfo> result;
    for (size_t i = 0; i < candidates.size(); ++i)
    {
        std::string testPath(gmx::Path::join(
                                     candidates[i].dir, candidates[i].name, filename));
        // TODO: Consider also checking that the directory can be listed.
        if (gmx::File::exists(testPath, gmx::File::returnFalseOnError))
        {
            result.push_back(candidates[i]);
        }
    }

    // TODO: Consider merging this into enumerateFiles(), such that the error
    // could also list the directories searched.
    if (result.empty())
    {
        std::string message
            = gmx::formatString("No force fields found (files with name '%s' "
                                "in subdirectories ending on '%s')",
                                filename, dirend);
        GMX_THROW(gmx::InvalidInputError(message));
    }

    return result;
}

bool fflib_fexist(const std::string &file)
{
    return !gmx::findLibraryFile(file, true, false).empty();
}


FILE *fflib_open(const std::string &file)
{
    std::string fileFullPath = gmx::findLibraryFile(file);
    fprintf(stderr, "Opening force field file %s\n", fileFullPath.c_str());
    return gmx_ffopen(fileFullPath, "r");
}
